#include "packing/improver.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <deque>
#include <mutex>
#include <optional>
#include <random>
#include <string>
#include <thread>

#include "packing/annealer.h"
#include "packing/crossover.h"
#include "packing/nfp.h"
#include "packing/rotation.h"
#include "packing/squeezer.h"
#include "packing/starts.h"

namespace packer::packing {
namespace {
using Clock = std::chrono::steady_clock;

enum class Op {
    Visit,
    Squeeze,
    Shake,
    Refine,
    Remove,
    Add,
    Random,
    Border,
    Compress,
    Crossover,
    Crystal,
    Dense,
    Symmetric,
    Ruin
};
constexpr std::size_t kOps = static_cast<std::size_t>(Op::Ruin) + 1;

// Bands of n, by their largest n, for the operator stats: what works differs with n.
constexpr std::array<int, 6> kBands = {12, 30, 60, 100, 150, 200};

std::size_t bandOf(int n) {
    for (std::size_t b = 0; b < kBands.size(); ++b) {
        if (n <= kBands[b]) return b;
    }
    return kBands.size() - 1;
}

const char *opName(Op op) {
    switch (op) {
    case Op::Visit:
        return "visit";
    case Op::Squeeze:
        return "squeeze";
    case Op::Shake:
        return "shake";
    case Op::Refine:
        return "refine";
    case Op::Remove:
        return "from n+1 minus a tree";
    case Op::Add:
        return "from n-1 plus a tree";
    case Op::Random:
        return "random start";
    case Op::Border:
        return "border squeeze";
    case Op::Compress:
        return "compress";
    case Op::Crossover:
        return "crossover";
    case Op::Crystal:
        return "crystal cut";
    case Op::Dense:
        return "dense start";
    case Op::Symmetric:
        return "symmetric start";
    case Op::Ruin:
        return "ruin and recreate";
    }

    return "?";
}

// A visit runs every operator on n; the other ops are single neighbour tries set off by an
// improvement next door. `visit` counts the visits to n, so each one gets fresh randomness.
struct Task {
    int n = 0;
    Op op = Op::Visit;
    int visit = 0;
};

// Mixes the run seed with where we are, so every step gets its own random stream. Runs still
// can't be repeated exactly: which n each thread works on, and on which best layout, depends
// on timing.
std::uint64_t mix(std::uint64_t seed, std::uint64_t a, std::uint64_t b, std::uint64_t c) {
    std::uint64_t z = seed;
    for (const std::uint64_t v : {a, b, c}) {
        z += 0x9e3779b97f4a7c15ULL + v; // the splitmix64 finaliser
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        z ^= z >> 31;
    }

    return z;
}

// Shared state of a run. Everything below the mutex is guarded by it.
class Run {
  public:
    Run(std::vector<Layout> &best, std::vector<std::string> &origins, const ImproveConfig &cfg,
        std::atomic<bool> &stop, std::uint64_t seed)
        : best_(best), origins_(origins), cfg_(cfg), stop_(stop), seed_(seed), pick_(seed) {
        const size_t count = best.size();
        busy_.assign(count, false);
        shrink_.assign(count, cfg.squeezeStartShrink);
        elites_.assign(count, {});
        seconds_.assign(count, 0);
        gain_.assign(count, 0);
        visits_.assign(count, 0);
        cursor_ = cfg.maxN;
    }

    void work() {
        while (!stop_.load()) {
            const std::optional<Task> task = next();
            if (!task) {
                // Every n is taken (more threads than n): wait for one to come free.
                std::this_thread::sleep_for(std::chrono::milliseconds(20));
                continue;
            }

            const Clock::time_point began = Clock::now();
            if (task->op == Op::Visit) {
                visit(*task);
            } else {
                single(*task);
            }

            std::lock_guard lock(mutex);
            busy_[task->n] = false;
            seconds_[task->n] += std::chrono::duration<double>(Clock::now() - began).count();
        }
    }

    double elapsed() const {
        return std::chrono::duration<double>(Clock::now() - started_).count();
    }

    // Per operator this run: tries, thread-seconds, improvements and score gained, most time
    // first; then the gain per hour of each in each band of n. Printed, and written to
    // cfg.statsPath if set. The caller holds the mutex.
    void printStats() const {
        printStats(stdout);
        if (cfg_.statsPath.empty()) return;
        if (std::FILE *f = std::fopen(cfg_.statsPath.c_str(), "w")) {
            printStats(f);
            std::fclose(f);
        }
    }

    void printStats(std::FILE *out) const {
        const std::array<OpStats, kOps> all = totals();
        double total = 0;
        std::array<std::size_t, kOps> order{};
        for (std::size_t k = 0; k < kOps; ++k) {
            order[k] = k;
            total += all[k].seconds;
        }
        std::sort(order.begin(), order.end(),
                  [&](std::size_t a, std::size_t b) { return all[a].seconds > all[b].seconds; });

        std::fprintf(out, "%-22s %8s %10s %6s %6s %10s %12s\n", "operator", "tries", "seconds",
                     "time", "wins", "gain", "gain/hour");
        for (const std::size_t k : order) {
            const OpStats &s = all[k];
            if (s.tries == 0) continue;
            std::fprintf(out, "%-22s %8ld %10.0f %5.1f%% %6d %10.6f %12.6f\n",
                         opName(static_cast<Op>(k)), s.tries, s.seconds,
                         total > 0 ? 100 * s.seconds / total : 0.0, s.wins, s.gain,
                         s.seconds > 0 ? 3600 * s.gain / s.seconds : 0.0);
        }

        std::fprintf(out, "\n%-22s", "gain/hour by n");
        for (std::size_t b = 0; b < kBands.size(); ++b) {
            const int lo = b == 0 ? 1 : kBands[b - 1] + 1;
            std::fprintf(out, " %9s",
                         (std::to_string(lo) + "-" + std::to_string(kBands[b])).c_str());
        }
        std::fprintf(out, "\n");
        for (const std::size_t k : order) {
            if (all[k].tries == 0) continue;
            std::fprintf(out, "%-22s", opName(static_cast<Op>(k)));
            for (std::size_t b = 0; b < kBands.size(); ++b) {
                const OpStats &s = stats_[b][k];
                if (s.seconds > 0) {
                    std::fprintf(out, " %9.6f", 3600 * s.gain / s.seconds);
                } else {
                    std::fprintf(out, " %9s", "-");
                }
            }
            std::fprintf(out, "\n");
        }
        std::fflush(out);
    }

    std::mutex mutex;

  private:
    std::vector<Layout> &best_;
    std::vector<std::string> &origins_;
    const ImproveConfig &cfg_;
    std::atomic<bool> &stop_;
    const std::uint64_t seed_;
    const Clock::time_point started_ = Clock::now();

    std::vector<bool> busy_;     // n some thread is working on
    std::vector<double> shrink_; // how much the next squeeze of each n tries to take off
    std::vector<std::vector<Layout>> elites_;   // per n, the best few distinct layouts, smallest first
    std::vector<double> seconds_; // per n, time spent on it this run
    std::vector<double> gain_;    // per n, score gained this run
    std::vector<int> visits_;     // per n, tasks started this run
    std::deque<Task> urgent_;     // neighbour tries set off by an improvement
    std::mt19937_64 pick_;        // for choosing between the sweep and the rates
    int cursor_ = 0;              // the next n the sweep visits
    int pass_ = 1;

    struct OpStats {
        long tries = 0;
        double seconds = 0; // thread-seconds
        int wins = 0;
        double gain = 0; // score gained
    };
    std::array<std::array<OpStats, kOps>, kBands.size()> stats_{}; // per band of n, per op

    // The stats of every band added up.
    std::array<OpStats, kOps> totals() const {
        std::array<OpStats, kOps> all{};
        for (const auto &band : stats_) {
            for (std::size_t k = 0; k < kOps; ++k) {
                all[k].tries += band[k].tries;
                all[k].seconds += band[k].seconds;
                all[k].wins += band[k].wins;
                all[k].gain += band[k].gain;
            }
        }
        return all;
    }

    // Charges the time from its creation to its end to one try of `op` on n.
    class Timed {
      public:
        Timed(Run &run, Op op, int n) : run_(run), op_(op), band_(bandOf(n)), began_(Clock::now()) {
        }
        ~Timed() {
            const double s = std::chrono::duration<double>(Clock::now() - began_).count();
            std::lock_guard lock(run_.mutex);
            OpStats &st = run_.stats_[band_][static_cast<std::size_t>(op_)];
            ++st.tries;
            st.seconds += s;
        }
        Timed(const Timed &) = delete;
        Timed &operator=(const Timed &) = delete;

      private:
        Run &run_;
        const Op op_;
        const std::size_t band_;
        const Clock::time_point began_;
    };

    // The smallest n worth working on: n = 1 is already exact.
    int firstN() const {
        return std::max(cfg_.minN, 2);
    }

    Task start(int n, Op op) {
        busy_[n] = true;
        return {n, op, ++visits_[n]};
    }

    // The next task: urgent neighbour tries first; otherwise the n with the best gain rate or
    // the next n of the sweep, which goes down so that every n sees n+1 freshly visited. An n
    // another thread is on is skipped. Nothing if every n is taken.
    std::optional<Task> next() {
        std::lock_guard lock(mutex);

        // One round through the urgent queue: tries for a busy n go back for later.
        for (std::size_t k = urgent_.size(); k-- > 0;) {
            const Task t = urgent_.front();
            urgent_.pop_front();
            if (!busy_[t.n]) return start(t.n, t.op);
            urgent_.push_back(t);
        }

        const bool byRate = std::uniform_real_distribution<double>(0, 1)(pick_) < cfg_.banditShare;
        if (byRate) {
            const int n = bestRate();
            if (n > 0) return start(n, Op::Visit);
        }

        for (int tried = 0; tried <= cfg_.maxN - firstN(); ++tried) {
            const int n = cursor_;
            advance();
            if (!busy_[n]) return start(n, Op::Visit);
        }

        return std::nullopt;
    }

    // The free n with the most score gained per second, each counted as if it had already run
    // banditPriorSeconds at the average rate (so an n that hasn't had its turn yet still
    // gets one). Ties are broken at random. 0 if every n is busy.
    int bestRate() {
        double gained = 0, spent = 0;
        for (int n = firstN(); n <= cfg_.maxN; ++n) {
            gained += gain_[n];
            spent += seconds_[n];
        }
        const double average = spent > 0 ? gained / spent : 0;
        const double prior = cfg_.banditPriorSeconds;

        std::uniform_real_distribution<double> tiny(0, 1e-15);
        int best = 0;
        double bestRate = -1;
        for (int n = firstN(); n <= cfg_.maxN; ++n) {
            if (busy_[n]) continue;
            const double rate = (gain_[n] + average * prior) / (seconds_[n] + prior) + tiny(pick_);
            if (rate > bestRate) {
                bestRate = rate;
                best = n;
            }
        }
        return best;
    }

    void advance() {
        if (--cursor_ >= firstN()) return;

        std::printf("[%5.0fs] sweep %d done, score %.6f\n", elapsed(), pass_, totalScore(best_));
        std::fflush(stdout);
        cursor_ = cfg_.maxN;
        ++pass_;
    }

    Layout bestCopy(int n) {
        std::lock_guard lock(mutex);
        return best_[n];
    }

    std::mt19937_64 rngFor(const Task &task, int step) const {
        return std::mt19937_64(mix(seed_, static_cast<std::uint64_t>(task.visit),
                                   static_cast<std::uint64_t>(task.n),
                                   static_cast<std::uint64_t>(step)));
    }

    void anneal(Layout &layout, long itersPerTree, double temp, double step, std::mt19937_64 &rng) {
        AnnealConfig ac;
        ac.iterations = itersPerTree * static_cast<long>(layout.size());
        ac.startTemp = temp;
        ac.startStep = step;
        ac.stop = &stop_;

        packing::anneal(layout, ac, rng);
        rotateWhole(layout);
    }

    SqueezeConfig squeezeConfig(const Layout &layout, long itersPerTree = 0) const {
        SqueezeConfig sc;
        if (itersPerTree == 0) itersPerTree = cfg_.squeezeItersPerTree;
        sc.iterations = itersPerTree * static_cast<long>(layout.size());
        sc.relaxIters = cfg_.squeezeRelaxIters;
        sc.glsRounds = cfg_.squeezeGlsRounds;
        sc.lineShare = cfg_.squeezeLineShare;
        sc.stop = &stop_;
        return sc;
    }

    CompressConfig compressConfig() const {
        CompressConfig cc;
        cc.iterations = cfg_.compressIters;
        cc.rounds = cfg_.compressRounds;
        cc.stop = &stop_;
        return cc;
    }

    std::optional<Layout> squeezeWith(const Layout &layout, double target, const SqueezeConfig &sc,
                                      std::mt19937_64 &rng) {
        std::optional<Layout> l = squeeze(layout, target, sc, rng);
        if (l) rotateWhole(*l);
        return l;
    }

    std::optional<Layout> squeezeInto(const Layout &layout, double target, std::mt19937_64 &rng,
                                      long itersPerTree = 0) {
        return squeezeWith(layout, target, squeezeConfig(layout, itersPerTree), rng);
    }

    double bestSide(int n) {
        std::lock_guard lock(mutex);
        return best_[n].side();
    }

    // The square the next squeeze of n aims for.
    double squeezeTarget(int n) {
        std::lock_guard lock(mutex);
        return best_[n].side() * (1 - shrink_[n]);
    }

    void adaptShrink(int n, bool worked) {
        std::lock_guard lock(mutex);
        const double s = worked ? shrink_[n] * 2 : shrink_[n] / 2;
        shrink_[n] = std::clamp(s, cfg_.squeezeStartShrink / 64, cfg_.squeezeMaxShrink);
    }

    // Keeps a valid layout among n's elites if it is close enough to the best and new.
    void addEliteLocked(int n, const Layout &l) {
        if (cfg_.eliteSize <= 0) return;
        const double side = l.side();
        if (side > best_[n].side() * (1 + cfg_.eliteSlack)) return;

        std::vector<Layout> &pool = elites_[n];
        // Different layouts practically never share a side to 1e-10: same side, same layout.
        auto same = [&](const Layout &e) { return std::abs(e.side() - side) < 1e-10; };
        if (std::any_of(pool.begin(), pool.end(), same)) return;

        pool.push_back(l);
        std::sort(pool.begin(), pool.end(),
                  [](const Layout &a, const Layout &b) { return a.side() < b.side(); });
        if (static_cast<int>(pool.size()) > cfg_.eliteSize) pool.pop_back();
    }

    void addElite(int n, const Layout &l) {
        std::lock_guard lock(mutex);
        addEliteLocked(n, l);
    }

    // A random elite of n other than the best itself, if there is one.
    std::optional<Layout> otherElite(int n, std::mt19937_64 &rng) {
        std::lock_guard lock(mutex);
        std::vector<const Layout *> others;
        for (const Layout &e : elites_[n]) {
            if (std::abs(e.side() - best_[n].side()) >= 1e-10) others.push_back(&e);
        }
        if (others.empty()) return std::nullopt;
        return *others[std::uniform_int_distribution<std::size_t>(0, others.size() - 1)(rng)];
    }

    // Whether to go ahead with a try of `op` on n (see ImproveConfig::opBandit). Rates are
    // those of n's band; the totals stand in until the band has gained anything.
    bool worth(Op op, int n, std::mt19937_64 &rng) {
        if (!cfg_.opBandit) return true;
        double chance = 1;
        {
            std::lock_guard lock(mutex);
            auto rateOf = [&](const std::array<OpStats, kOps> &st, double &average) {
                double gained = 0, spent = 0;
                for (const OpStats &s : st) {
                    gained += s.gain;
                    spent += s.seconds;
                }
                if (gained <= 0 || spent <= 0) return false;
                average = gained / spent;
                const double prior = cfg_.opPriorSeconds;
                const OpStats &s = st[static_cast<std::size_t>(op)];
                chance = (s.gain + average * prior) / (s.seconds + prior) / average;
                return true;
            };
            double average = 0;
            const bool known = rateOf(stats_[bandOf(n)], average) || rateOf(totals(), average);
            if (!known) return true;
            chance = std::clamp(chance, cfg_.opMinChance, 1.0);
        }
        return std::uniform_real_distribution<double>(0.0, 1.0)(rng) < chance;
    }

    // Squeezes into ever smaller squares while that works, offering every success. Returns
    // the smallest layout reached.
    Layout descend(int n, Op op, Layout l, int tries, std::mt19937_64 &rng,
                   bool symmetric = false) {
        double shrink = 4e-3;
        for (int k = 0; k < tries && !stop_.load(); ++k) {
            SqueezeConfig sc = squeezeConfig(l);
            sc.symmetric = symmetric;
            std::optional<Layout> r = squeezeWith(l, l.side() * (1 - shrink), sc, rng);
            if (r) {
                l = std::move(*r);
                offer(n, op, l);
                addElite(n, l);
                shrink = std::min(shrink * 1.5, 1e-2);
            } else {
                shrink /= 3;
                if (shrink < 1e-5) break;
            }
        }
        return l;
    }

    // Every operator on n once, in order. `step` numbers the random streams they draw, so
    // each try gets its own.
    void visit(const Task &task) {
        int step = 0;
        tryRemove(task, step);
        tryAdd(task, step);
        squeezeLoop(task, step);
        tryBorder(task, step);
        compressLoop(task, step);
        shakeWalk(task, step);
        tryCrossover(task, step);
        ruinWalk(task, step);
        tryCrystal(task, step);
        denseStarts(task, step);
        refineWalk(task, step);
        tryRandom(task, step);
    }

    // n+1 minus each of its best few trees to remove, squeezed into the best square so far.
    void tryRemove(const Task &task, int &step) {
        const int n = task.n;
        const bool hasBiggerN = n + 1 < static_cast<int>(best_.size());
        if (!hasBiggerN) return;

        const Layout from = bestCopy(n + 1);
        const std::vector<std::size_t> order = from.removalOrder();
        const int tries = std::min<int>(cfg_.removeCandidates, static_cast<int>(order.size()));
        for (int k = 0; k < tries && !stop_.load(); ++k) {
            std::mt19937_64 rng = rngFor(task, step++);
            if (!worth(Op::Remove, n, rng)) continue;
            const Timed timed(*this, Op::Remove, n);
            Layout l = from;
            l.remove(order[k]);
            offerIf(n, Op::Remove, squeezeInto(l, squeezeTarget(n), rng));
        }
    }

    // n-1 plus a tree, squeezed into the best square so far.
    void tryAdd(const Task &task, int &step) {
        const int n = task.n;
        std::mt19937_64 rng = rngFor(task, step++);
        if (n - 1 < 1 || stop_.load() || !worth(Op::Add, n, rng)) return;

        const Timed timed(*this, Op::Add, n);
        offerIf(n, Op::Add, addAndSqueeze(n, rng));
    }

    // Squeeze the best layout until it fails often enough.
    void squeezeLoop(const Task &task, int &step) {
        const int n = task.n;
        for (int failures = 0; failures < cfg_.squeezeFailures && !stop_.load();) {
            const Timed timed(*this, Op::Squeeze, n);
            std::mt19937_64 rng = rngFor(task, step++);
            const bool worked =
                offerIf(n, Op::Squeeze, squeezeInto(bestCopy(n), squeezeTarget(n), rng));
            adaptShrink(n, worked);
            failures = worked ? 0 : failures + 1;
        }
    }

    // Only the trees near the walls move; the inside stays as it is.
    void tryBorder(const Task &task, int &step) {
        const int n = task.n;
        std::mt19937_64 rng = rngFor(task, step++);
        if (stop_.load() || !worth(Op::Border, n, rng)) return;

        const Timed timed(*this, Op::Border, n);
        const Layout b = bestCopy(n);
        SqueezeConfig sc = squeezeConfig(b);
        sc.movable = borderTrees(b, cfg_.borderWidth);
        const bool hasInside = std::count(sc.movable.begin(), sc.movable.end(), false) > 0;
        if (hasInside) offerIf(n, Op::Border, squeezeWith(b, squeezeTarget(n), sc, rng));
    }

    // The side as a variable: the walls close in as far as the trees let them. Each win is the
    // next start, until it fails often enough; every other try works twice as long.
    void compressLoop(const Task &task, int &step) {
        const int n = task.n;
        for (int k = 0, failures = 0;
             k < cfg_.compressTries && failures < cfg_.compressFailures && !stop_.load(); ++k) {
            std::mt19937_64 rng = rngFor(task, step++);
            if (!worth(Op::Compress, n, rng)) break;
            const Timed timed(*this, Op::Compress, n);
            Layout c = bestCopy(n);
            CompressConfig cc = compressConfig();
            if (k % 2 == 1) {
                cc.iterations *= 2;
                cc.rounds *= 2;
            }
            const double side = compress(c, cc);
            const bool worked = side < bestSide(n) &&
                                offerIf(n, Op::Compress, squeezeInto(c, side * (1 + 1e-7), rng));
            failures = worked ? 0 : failures + 1;
        }
    }

    // Shake and settle: a walk through arrangements as good as the best, each one a fresh
    // start for a squeeze into a smaller square.
    void shakeWalk(const Task &task, int &step) {
        walk(task, step, Op::Shake, cfg_.shakeTries,
             [&](const Layout &current, int k, std::mt19937_64 &rng) {
                 return shake(current, 1 + k % cfg_.shakeKicks, rng);
             });
    }

    // Crossover: half of the best with half of another good layout, settled.
    void tryCrossover(const Task &task, int &step) {
        const int n = task.n;
        for (int k = 0; k < cfg_.crossoverTries && !stop_.load(); ++k) {
            std::mt19937_64 rng = rngFor(task, step++);
            if (!worth(Op::Crossover, n, rng)) continue;
            const Timed timed(*this, Op::Crossover, n);
            const Layout a = bestCopy(n);
            const Layout b = partnerFor(n, rng);
            bool worked = false;
            settleAndOffer(n, Op::Crossover, crossover(a, b, n, rng), rng, worked);
        }
    }

    // Ruin and recreate: a patch of trees taken out and put back at the best spots, settled
    // into the best square; a walk like shaking, each step squeezed once.
    void ruinWalk(const Task &task, int &step) {
        const int n = task.n;
        if (n < cfg_.ruinMinN) return;

        const std::size_t most =
            std::max<std::size_t>(3, static_cast<std::size_t>(cfg_.ruinShare * n));
        walk(task, step, Op::Ruin, cfg_.ruinTries,
             [&](const Layout &current, int, std::mt19937_64 &rng) {
                 const std::size_t count = std::uniform_int_distribution<std::size_t>(3, most)(rng);
                 const Box square =
                     geometry::squareAround(geometry::centre(current.bounds()), bestSide(n));
                 return ruinAndRecreate(current, square, count, rng);
             });
    }

    // A crystal cut: compress its border first (the inside is already a perfect pattern), then
    // all of it, and squeeze down from there if it came close.
    void tryCrystal(const Task &task, int &step) {
        const int n = task.n;
        const bool cutCrystal = !cfg_.crystals.empty() && n >= cfg_.crystalMinN;
        std::mt19937_64 rng = rngFor(task, step++);
        if (!cutCrystal || stop_.load() || !worth(Op::Crystal, n, rng)) return;

        const Timed timed(*this, Op::Crystal, n);
        const std::size_t which = static_cast<std::size_t>(task.visit) % cfg_.crystals.size();
        const Crystal &c = cfg_.crystals[which];
        Layout cut = bestCrystalCut(c, n, cfg_.crystalAngles, rng).first;
        offer(n, Op::Crystal, cut);

        CompressConfig cc = compressConfig();
        cc.movable = borderTrees(cut, cfg_.borderWidth);
        compress(cut, cc);
        cc.movable.clear();
        const double side = compress(cut, cc);
        std::optional<Layout> r = squeezeInto(cut, side * (1 + 1e-7), rng, cfg_.shakeItersPerTree);
        if (!r || stop_.load()) return;

        offer(n, Op::Crystal, *r);
        addElite(n, *r); // a different arrangement: a good crossover partner
        if (r->side() < bestSide(n) * (1 + cfg_.crystalSlack)) {
            descend(n, Op::Crystal, *r, cfg_.crystalDescent, rng);
        }
    }

    // Dense starts: random overlapping trees in a roomy square, compressed and squeezed valid,
    // then down if they came close. Half are point-symmetric, kept so on the way.
    void denseStarts(const Task &task, int &step) {
        const int n = task.n;
        if (n < 2 || n > cfg_.denseMax) return;

        std::uniform_real_distribution<double> unit(0.0, 1.0);
        const int starts = cfg_.denseStarts * (n <= cfg_.denseSmallMax ? cfg_.denseSmallBoost : 1);
        for (int k = 0; k < starts && !stop_.load(); ++k) {
            std::mt19937_64 rng = rngFor(task, step++);
            const bool symmetric = n >= 4 && unit(rng) < 0.5;
            const Op op = symmetric ? Op::Symmetric : Op::Dense;
            if (!worth(op, n, rng)) continue;
            const Timed timed(*this, op, n);
            const double roomy = 1.3 * bestSide(n);
            Layout l = symmetric ? symmetricStart(n, roomy, rng) : denseStart(n, roomy, rng);

            CompressConfig cc = compressConfig();
            cc.symmetric = symmetric;
            const double side = compress(l, cc);
            SqueezeConfig sc = squeezeConfig(l);
            sc.symmetric = symmetric;
            sc.glsRounds = cfg_.denseGlsRounds;
            std::optional<Layout> r = squeezeWith(l, side * (1 + 1e-6), sc, rng);
            if (!r) continue;
            offer(n, op, *r);
            if (r->side() < bestSide(n) * (1 + cfg_.denseSlack)) {
                descend(n, op, *r, cfg_.denseDescent, rng, symmetric);
            }
        }
    }

    // Refine: a walk of re-anneals, hotter at each step of the ladder. The anneal can drop a
    // tree into a hole, which a squeeze can't; the squeeze then takes up the room.
    void refineWalk(const Task &task, int &step) {
        const int n = task.n;
        Layout walk = bestCopy(n);
        for (int k = 0; k < cfg_.refineTries && !stop_.load(); ++k) {
            std::mt19937_64 rng = rngFor(task, step++);
            if (!worth(Op::Refine, n, rng)) continue;
            const Timed timed(*this, Op::Refine, n);
            const double temp =
                cfg_.refineTemps[static_cast<std::size_t>(k) % cfg_.refineTemps.size()];
            Layout l = walk;
            anneal(l, cfg_.refineItersPerTree, temp, cfg_.refineStep, rng);
            const bool worked = offer(n, Op::Refine, l);
            std::optional<Layout> r = squeezeInto(l, squeezeTarget(n), rng);
            if (offerIf(n, Op::Refine, r) || worked) {
                walk = bestCopy(n);
            } else if (l.isValid() && l.side() <= bestSide(n) * (1 + cfg_.eliteSlack)) {
                addElite(n, l);
                walk = std::move(l); // as good as the best, near enough: walk on from it
            }
        }
    }

    // Small n: a fresh random start every other visit, annealed hot and long.
    void tryRandom(const Task &task, int &step) {
        const int n = task.n;
        const bool randomTurn = n <= cfg_.smallMax && task.visit % 2 == 1;
        std::mt19937_64 rng = rngFor(task, step++);
        if (!randomTurn || stop_.load() || !worth(Op::Random, n, rng)) return;

        const Timed timed(*this, Op::Random, n);
        Layout l = randomStart(n, rng);
        anneal(l, cfg_.randomItersPerTree, cfg_.randomTemp, cfg_.randomStep, rng);
        offer(n, Op::Random, std::move(l));
    }

    // A walk through arrangements as good as the best: each try perturbs the current layout
    // (perturb(current, k, rng)), settles it (see settleAndOffer()) and walks on from there,
    // or from the new best after a win.
    template <class Perturb>
    void walk(const Task &task, int &step, Op op, int tries, Perturb perturb) {
        const int n = task.n;
        Layout current = bestCopy(n);
        for (int k = 0; k < tries && !stop_.load(); ++k) {
            std::mt19937_64 rng = rngFor(task, step++);
            if (!worth(op, n, rng)) continue;
            const Timed timed(*this, op, n);
            bool worked = false;
            std::optional<Layout> settled =
                settleAndOffer(n, op, perturb(current, k, rng), rng, worked);
            if (!settled) continue;

            current = worked ? bestCopy(n) : std::move(*settled);
        }
    }

    // Settles a perturbed layout into the best square with a short squeeze, keeps it among the
    // elites, and offers it, or else a squeeze of it into a smaller square. Returns the settled
    // layout, or nothing if it couldn't be settled; `worked` says whether best[n] improved.
    std::optional<Layout> settleAndOffer(int n, Op op, const Layout &perturbed,
                                         std::mt19937_64 &rng, bool &worked) {
        worked = false;
        std::optional<Layout> settled =
            squeezeInto(perturbed, bestSide(n), rng, cfg_.shakeItersPerTree);
        if (!settled) return std::nullopt;

        addElite(n, *settled);
        worked =
            offer(n, op, *settled) || offerIf(n, op, squeezeInto(*settled, squeezeTarget(n), rng));
        return settled;
    }

    // The other parent for a crossover with best[n]: another elite, the best itself (turned
    // by crossover()), or a neighbour n±1 (crossover() fixes the count).
    Layout partnerFor(int n, std::mt19937_64 &rng) {
        const double r = std::uniform_real_distribution<double>(0.0, 1.0)(rng);
        if (r < 0.5) {
            if (std::optional<Layout> e = otherElite(n, rng)) return *e;
        }
        if (r < 0.8) return bestCopy(n);

        const bool up = n + 1 < static_cast<int>(best_.size()) && (r < 0.9 || n - 1 < 1);
        return bestCopy(up ? n + 1 : n - 1);
    }

    // One neighbour try.
    void single(const Task &task) {
        const int n = task.n;
        const Timed timed(*this, task.op, n);
        std::mt19937_64 rng = rngFor(task, 1000 + static_cast<int>(task.op));

        if (task.op == Op::Remove) {
            const bool hasBiggerN = n + 1 < static_cast<int>(best_.size());
            if (!hasBiggerN) return;

            const Layout l = withoutTree(bestCopy(n + 1));
            offerIf(n, Op::Remove, squeezeInto(l, squeezeTarget(n), rng));
        } else if (task.op == Op::Add) {
            offerIf(n, Op::Add, addAndSqueeze(n, rng));
        }
    }

    // best[n-1] plus a tree, overlap and all, squeezed into a square a bit smaller than
    // best[n]'s. Half the time the tree goes to the best spot (see bestInsertion()), half the
    // time anywhere.
    std::optional<Layout> addAndSqueeze(int n, std::mt19937_64 &rng) {
        Layout l = bestCopy(n - 1);
        const double target = squeezeTarget(n);
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        const Box b = l.bounds();
        if (unit(rng) < 0.5) {
            const Box square =
                geometry::squareAround(geometry::centre(b), std::max(target, l.side()));
            l.add(bestInsertion(l, square, 8, rng));
        } else {
            l.add({b.minX + unit(rng) * (b.maxX - b.minX), b.minY + unit(rng) * (b.maxY - b.minY),
                   360.0 * unit(rng)});
        }

        return squeezeInto(l, target, rng);
    }

    bool offerIf(int n, Op op, std::optional<Layout> l) {
        return l && offer(n, op, std::move(*l));
    }

    // Replaces best[n] if the layout is valid and smaller. Returns whether it did.
    bool offer(int n, Op op, Layout result) {
        // Check validity outside the lock; it is the slowest part of accepting.
        const bool usable = static_cast<int>(result.size()) == n && result.isValid();
        if (!usable) return false;

        std::lock_guard lock(mutex);
        const double before = best_[n].side();
        const double after = result.side();
        const bool improved = after < before - 1e-12;
        if (!improved) return false;

        best_[n] = std::move(result);
        origins_[n] = opName(op);
        const double gain = (before * before - after * after) / n;
        gain_[n] += gain;
        OpStats &st = stats_[bandOf(n)][static_cast<std::size_t>(op)];
        ++st.wins;
        st.gain += gain;
        addEliteLocked(n, best_[n]);

        std::printf("[%5.0fs] n=%-3d %.6f -> %.6f  (%s)  score %.6f\n", elapsed(), n,
                    before * before / n, after * after / n, opName(op), totalScore(best_));
        std::fflush(stdout);

        // A smaller n can now come from this one minus a tree, and a bigger one from it plus a
        // tree. The sweep is going down, so it reaches n-1 soon anyway: only n+1 needs a nudge,
        // and only after a real gain (compress wins tiny amounts often).
        if (n + 1 <= cfg_.maxN && gain > cfg_.urgentMinGain) urgent_.push_back({n + 1, Op::Add, 0});

        return true;
    }
};
} // namespace

double totalScore(const std::vector<Layout> &layouts) {
    double total = 0;

    for (std::size_t n = 1; n < layouts.size(); ++n) {
        const double s = layouts[n].side();
        total += s * s / static_cast<double>(n);
    }

    return total;
}

void improve(std::vector<Layout> &best, std::vector<std::string> &origins, const ImproveConfig &cfg,
             std::atomic<bool> &stop, const SaveFn &save) {
    const unsigned halfTheCores = std::max(1u, std::thread::hardware_concurrency() / 2);
    const int threads = cfg.threads > 0 ? cfg.threads : halfTheCores;

    const auto clockSeed = static_cast<unsigned>(Clock::now().time_since_epoch().count());
    const unsigned seed = cfg.seed != 0 ? cfg.seed : clockSeed;

    std::printf("Improving n=%d..%d for %.1f min on %d threads (seed %u), from score %.6f\n",
                cfg.minN, cfg.maxN, cfg.minutes, threads, seed, totalScore(best));
    std::fflush(stdout);

    Run run(best, origins, cfg, stop, seed);
    std::vector<std::thread> workers;
    for (int t = 0; t < threads; ++t) {
        workers.emplace_back([&run] { run.work(); });
    }

    // The main thread keeps time and saves now and then.
    double lastSave = 0;
    const double budget = cfg.minutes * 60;
    while (!stop.load()) {
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        const double now = run.elapsed();

        if (now >= budget) stop = true;

        if (now - lastSave >= cfg.saveSeconds) {
            std::lock_guard lock(run.mutex);
            save(best, origins);
            run.printStats();
            lastSave = now;
        }
    }

    for (std::thread &w : workers) {
        w.join();
    }

    save(best, origins);
    run.printStats();
}
} // namespace packer::packing
