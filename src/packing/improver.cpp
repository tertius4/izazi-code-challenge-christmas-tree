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
    for (std::size_t band = 0; band < kBands.size(); ++band) {
        if (n <= kBands[band]) return band;
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
std::uint64_t mix(std::uint64_t seed, std::uint64_t visit, std::uint64_t n, std::uint64_t step) {
    std::uint64_t hash = seed;
    for (const std::uint64_t value : {visit, n, step}) {
        hash += 0x9e3779b97f4a7c15ULL + value; // the splitmix64 finaliser
        hash = (hash ^ (hash >> 30)) * 0xbf58476d1ce4e5b9ULL;
        hash = (hash ^ (hash >> 27)) * 0x94d049bb133111ebULL;
        hash ^= hash >> 31;
    }

    return hash;
}

// Shared state of a run. Everything below the mutex is guarded by it.
class Run {
  public:
    Run(std::vector<Layout> &best, std::vector<std::string> &origins, const ImproveConfig &cfg,
        std::atomic<bool> &stop, std::uint64_t seed)
        : best_(best), origins_(origins), cfg_(cfg), stop_(stop), seed_(seed), pick_(seed) {
        const std::size_t count = best.size();
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
                neighbourTry(*task);
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
        if (std::FILE *file = std::fopen(cfg_.statsPath.c_str(), "w")) {
            printStats(file);
            std::fclose(file);
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
            const OpStats &stats = all[k];
            if (stats.tries == 0) continue;
            std::fprintf(out, "%-22s %8ld %10.0f %5.1f%% %6d %10.6f %12.6f\n",
                         opName(static_cast<Op>(k)), stats.tries, stats.seconds,
                         total > 0 ? 100 * stats.seconds / total : 0.0, stats.wins, stats.gain,
                         stats.seconds > 0 ? 3600 * stats.gain / stats.seconds : 0.0);
        }

        std::fprintf(out, "\n%-22s", "gain/hour by n");
        for (std::size_t band = 0; band < kBands.size(); ++band) {
            const int first = band == 0 ? 1 : kBands[band - 1] + 1;
            std::fprintf(out, " %9s",
                         (std::to_string(first) + "-" + std::to_string(kBands[band])).c_str());
        }
        std::fprintf(out, "\n");
        for (const std::size_t k : order) {
            if (all[k].tries == 0) continue;
            std::fprintf(out, "%-22s", opName(static_cast<Op>(k)));
            for (std::size_t band = 0; band < kBands.size(); ++band) {
                const OpStats &stats = stats_[band][k];
                if (stats.seconds > 0) {
                    std::fprintf(out, " %9.6f", 3600 * stats.gain / stats.seconds);
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

    // Per n, the best few distinct layouts, smallest first.
    std::vector<std::vector<Layout>> elites_;
    std::vector<bool> busy_;      // n some thread is working on
    std::vector<double> shrink_;  // how much the next squeeze of each n tries to take off
    std::vector<double> seconds_; // per n, time spent on it this run
    std::vector<double> gain_;    // per n, score gained this run
    std::vector<int> visits_;     // per n, tasks started this run
    std::deque<Task> urgent_;     // neighbour tries set off by an improvement
    std::mt19937_64 pick_;        // for choosing between the sweep and the rates
    int cursor_ = 0;              // the next n the sweep visits
    int pass_ = 1;                // the sweep's number, counted from 1

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
            const double seconds = std::chrono::duration<double>(Clock::now() - began_).count();
            std::lock_guard lock(run_.mutex);
            OpStats &stats = run_.stats_[band_][static_cast<std::size_t>(op_)];
            ++stats.tries;
            stats.seconds += seconds;
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
            const Task queued = urgent_.front();
            urgent_.pop_front();
            if (!busy_[queued.n]) return start(queued.n, queued.op);
            urgent_.push_back(queued);
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
        int bestN = 0;
        double highest = -1;
        for (int n = firstN(); n <= cfg_.maxN; ++n) {
            if (busy_[n]) continue;
            const double rate = (gain_[n] + average * prior) / (seconds_[n] + prior) + tiny(pick_);
            if (rate > highest) {
                highest = rate;
                bestN = n;
            }
        }
        return bestN;
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
        AnnealConfig annealCfg;
        annealCfg.iterations = itersPerTree * static_cast<long>(layout.size());
        annealCfg.startTemp = temp;
        annealCfg.startStep = step;
        annealCfg.stop = &stop_;

        packing::anneal(layout, annealCfg, rng);
        rotateWhole(layout);
    }

    SqueezeConfig squeezeConfig(const Layout &layout, long itersPerTree = 0) const {
        SqueezeConfig squeezeCfg;
        if (itersPerTree == 0) itersPerTree = cfg_.squeezeItersPerTree;
        squeezeCfg.iterations = itersPerTree * static_cast<long>(layout.size());
        squeezeCfg.relaxIters = cfg_.squeezeRelaxIters;
        squeezeCfg.glsRounds = cfg_.squeezeGlsRounds;
        squeezeCfg.lineShare = cfg_.squeezeLineShare;
        squeezeCfg.stop = &stop_;
        return squeezeCfg;
    }

    CompressConfig compressConfig() const {
        CompressConfig compressCfg;
        compressCfg.iterations = cfg_.compressIters;
        compressCfg.rounds = cfg_.compressRounds;
        compressCfg.stop = &stop_;
        return compressCfg;
    }

    std::optional<Layout> squeezeWith(const Layout &layout, double target,
                                      const SqueezeConfig &squeezeCfg, std::mt19937_64 &rng) {
        std::optional<Layout> squeezed = squeeze(layout, target, squeezeCfg, rng);
        if (squeezed) rotateWhole(*squeezed);
        return squeezed;
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
        const double next = worked ? shrink_[n] * 2 : shrink_[n] / 2;
        shrink_[n] = std::clamp(next, cfg_.squeezeStartShrink / 64, cfg_.squeezeMaxShrink);
    }

    // Keeps a valid layout among n's elites if it is close enough to the best and new.
    void addEliteLocked(int n, const Layout &layout) {
        if (cfg_.eliteSize <= 0) return;
        const double side = layout.side();
        if (side > best_[n].side() * (1 + cfg_.eliteSlack)) return;

        std::vector<Layout> &pool = elites_[n];
        // Different layouts practically never share a side to 1e-10: same side, same layout.
        auto same = [&](const Layout &elite) { return std::abs(elite.side() - side) < 1e-10; };
        if (std::any_of(pool.begin(), pool.end(), same)) return;

        pool.push_back(layout);
        std::sort(pool.begin(), pool.end(),
                  [](const Layout &a, const Layout &b) { return a.side() < b.side(); });
        if (static_cast<int>(pool.size()) > cfg_.eliteSize) pool.pop_back();
    }

    void addElite(int n, const Layout &layout) {
        std::lock_guard lock(mutex);
        addEliteLocked(n, layout);
    }

    // A random elite of n other than the best itself, if there is one.
    std::optional<Layout> otherElite(int n, std::mt19937_64 &rng) {
        std::lock_guard lock(mutex);
        std::vector<const Layout *> others;
        for (const Layout &elite : elites_[n]) {
            if (std::abs(elite.side() - best_[n].side()) >= 1e-10) others.push_back(&elite);
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
            auto rateOf = [&](const std::array<OpStats, kOps> &allOps, double &average) {
                double gained = 0, spent = 0;
                for (const OpStats &stats : allOps) {
                    gained += stats.gain;
                    spent += stats.seconds;
                }
                if (gained <= 0 || spent <= 0) return false;
                average = gained / spent;
                const double prior = cfg_.opPriorSeconds;
                const OpStats &stats = allOps[static_cast<std::size_t>(op)];
                chance = (stats.gain + average * prior) / (stats.seconds + prior) / average;
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
    Layout descend(int n, Op op, Layout layout, int tries, std::mt19937_64 &rng,
                   bool symmetric = false) {
        double shrink = 4e-3;
        for (int k = 0; k < tries && !stop_.load(); ++k) {
            SqueezeConfig squeezeCfg = squeezeConfig(layout);
            squeezeCfg.symmetric = symmetric;
            std::optional<Layout> squeezed =
                squeezeWith(layout, layout.side() * (1 - shrink), squeezeCfg, rng);
            if (squeezed) {
                layout = std::move(*squeezed);
                offer(n, op, layout);
                addElite(n, layout);
                shrink = std::min(shrink * 1.5, 1e-2);
            } else {
                shrink /= 3;
                if (shrink < 1e-5) break;
            }
        }
        return layout;
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
            Layout layout = from;
            layout.remove(order[k]);
            offerIf(n, Op::Remove, squeezeInto(layout, squeezeTarget(n), rng));
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
        const Layout best = bestCopy(n);
        SqueezeConfig squeezeCfg = squeezeConfig(best);
        squeezeCfg.movable = borderTrees(best, cfg_.borderWidth);
        const auto &movable = squeezeCfg.movable;
        const bool hasInside = std::count(movable.begin(), movable.end(), false) > 0;
        if (hasInside) offerIf(n, Op::Border, squeezeWith(best, squeezeTarget(n), squeezeCfg, rng));
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
            Layout compressed = bestCopy(n);
            CompressConfig compressCfg = compressConfig();
            if (k % 2 == 1) {
                compressCfg.iterations *= 2;
                compressCfg.rounds *= 2;
            }
            const double side = compress(compressed, compressCfg);
            const bool worked =
                side < bestSide(n) &&
                offerIf(n, Op::Compress, squeezeInto(compressed, side * (1 + 1e-7), rng));
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
            const Layout best = bestCopy(n);
            const Layout partner = partnerFor(n, rng);
            bool worked = false;
            settleAndOffer(n, Op::Crossover, crossover(best, partner, n, rng), rng, worked);
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
        const Crystal &crystal = cfg_.crystals[which];
        Layout cut = bestCrystalCut(crystal, n, cfg_.crystalAngles, rng).first;
        offer(n, Op::Crystal, cut);

        CompressConfig compressCfg = compressConfig();
        compressCfg.movable = borderTrees(cut, cfg_.borderWidth);
        compress(cut, compressCfg);
        compressCfg.movable.clear();
        const double side = compress(cut, compressCfg);
        std::optional<Layout> squeezed =
            squeezeInto(cut, side * (1 + 1e-7), rng, cfg_.shakeItersPerTree);
        if (!squeezed || stop_.load()) return;

        offer(n, Op::Crystal, *squeezed);
        addElite(n, *squeezed); // a different arrangement: a good crossover partner
        if (squeezed->side() < bestSide(n) * (1 + cfg_.crystalSlack)) {
            descend(n, Op::Crystal, *squeezed, cfg_.crystalDescent, rng);
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
            Layout start = symmetric ? symmetricStart(n, roomy, rng) : denseStart(n, roomy, rng);

            CompressConfig compressCfg = compressConfig();
            compressCfg.symmetric = symmetric;
            const double side = compress(start, compressCfg);
            SqueezeConfig squeezeCfg = squeezeConfig(start);
            squeezeCfg.symmetric = symmetric;
            squeezeCfg.glsRounds = cfg_.denseGlsRounds;
            std::optional<Layout> squeezed = squeezeWith(start, side * (1 + 1e-6), squeezeCfg, rng);
            if (!squeezed) continue;
            offer(n, op, *squeezed);
            if (squeezed->side() < bestSide(n) * (1 + cfg_.denseSlack)) {
                descend(n, op, *squeezed, cfg_.denseDescent, rng, symmetric);
            }
        }
    }

    // Refine: a walk of re-anneals, hotter at each step of the ladder. The anneal can drop a
    // tree into a hole, which a squeeze can't; the squeeze then takes up the room.
    void refineWalk(const Task &task, int &step) {
        const int n = task.n;
        Layout current = bestCopy(n);
        for (int k = 0; k < cfg_.refineTries && !stop_.load(); ++k) {
            std::mt19937_64 rng = rngFor(task, step++);
            if (!worth(Op::Refine, n, rng)) continue;
            const Timed timed(*this, Op::Refine, n);
            const double temp =
                cfg_.refineTemps[static_cast<std::size_t>(k) % cfg_.refineTemps.size()];
            Layout annealed = current;
            anneal(annealed, cfg_.refineItersPerTree, temp, cfg_.refineStep, rng);
            const bool worked = offer(n, Op::Refine, annealed);
            std::optional<Layout> squeezed = squeezeInto(annealed, squeezeTarget(n), rng);
            if (offerIf(n, Op::Refine, squeezed) || worked) {
                current = bestCopy(n);
            } else if (annealed.isValid() &&
                       annealed.side() <= bestSide(n) * (1 + cfg_.eliteSlack)) {
                addElite(n, annealed);
                current = std::move(annealed); // as good as the best, near enough: walk on from it
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
        Layout layout = randomStart(n, rng);
        anneal(layout, cfg_.randomItersPerTree, cfg_.randomTemp, cfg_.randomStep, rng);
        offer(n, Op::Random, std::move(layout));
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
        const double roll = std::uniform_real_distribution<double>(0.0, 1.0)(rng);
        if (roll < 0.5) {
            if (std::optional<Layout> elite = otherElite(n, rng)) return *elite;
        }
        if (roll < 0.8) return bestCopy(n);

        const bool up = n + 1 < static_cast<int>(best_.size()) && (roll < 0.9 || n - 1 < 1);
        return bestCopy(up ? n + 1 : n - 1);
    }

    // One neighbour try: n from n+1 minus a tree, or from n-1 plus one.
    void neighbourTry(const Task &task) {
        const int n = task.n;
        const Timed timed(*this, task.op, n);
        std::mt19937_64 rng = rngFor(task, 1000 + static_cast<int>(task.op));

        if (task.op == Op::Remove) {
            const bool hasBiggerN = n + 1 < static_cast<int>(best_.size());
            if (!hasBiggerN) return;

            const Layout layout = withoutTree(bestCopy(n + 1));
            offerIf(n, Op::Remove, squeezeInto(layout, squeezeTarget(n), rng));
        } else if (task.op == Op::Add) {
            offerIf(n, Op::Add, addAndSqueeze(n, rng));
        }
    }

    // best[n-1] plus a tree, overlap and all, squeezed into a square a bit smaller than
    // best[n]'s. Half the time the tree goes to the best spot (see bestInsertion()), half the
    // time anywhere.
    std::optional<Layout> addAndSqueeze(int n, std::mt19937_64 &rng) {
        Layout layout = bestCopy(n - 1);
        const double target = squeezeTarget(n);
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        const Box bounds = layout.bounds();
        if (unit(rng) < 0.5) {
            const Box square =
                geometry::squareAround(geometry::centre(bounds), std::max(target, layout.side()));
            layout.add(bestInsertion(layout, square, 8, rng));
        } else {
            layout.add({bounds.minX + unit(rng) * (bounds.maxX - bounds.minX),
                        bounds.minY + unit(rng) * (bounds.maxY - bounds.minY), 360.0 * unit(rng)});
        }

        return squeezeInto(layout, target, rng);
    }

    bool offerIf(int n, Op op, std::optional<Layout> layout) {
        return layout && offer(n, op, std::move(*layout));
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
        OpStats &stats = stats_[bandOf(n)][static_cast<std::size_t>(op)];
        ++stats.wins;
        stats.gain += gain;
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
