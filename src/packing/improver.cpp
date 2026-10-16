#include "packing/improver.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <deque>
#include <mutex>
#include <optional>
#include <random>
#include <thread>

#include "packing/annealer.h"
#include "packing/rotation.h"
#include "packing/squeezer.h"
#include "packing/starts.h"

namespace packer::packing {
namespace {
using Clock = std::chrono::steady_clock;

enum class Op { Visit, Squeeze, Shake, Refine, Remove, Add, Random };

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
    }

    return "?";
}

// A visit runs every operator on n; the other ops are single neighbour tries set off by an
// improvement next door.
struct Task {
    int n = 0;
    Op op = Op::Visit;
    int pass = 0;
};

// Mixes the run seed with where we are, so every step gets its own random stream and a run can
// be repeated with the same seed.
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
        : best_(best), origins_(origins), cfg_(cfg), stop_(stop), seed_(seed) {
        const size_t count = best.size();
        busy_.assign(count, false);
        shrink_.assign(count, cfg.squeezeStartShrink);
        cursor_ = cfg.maxN;
    }

    void work() {
        while (!stop_.load()) {
            const Task task = next();
            if (task.op == Op::Visit) {
                visit(task);
            } else {
                single(task);
            }

            std::lock_guard lock(mutex);
            busy_[task.n] = false;
        }
    }

    double elapsed() const {
        return std::chrono::duration<double>(Clock::now() - started_).count();
    }

    std::mutex mutex;

  private:
    std::vector<Layout> &best_;
    std::vector<std::string> &origins_;
    const ImproveConfig &cfg_;
    std::atomic<bool> &stop_;
    const std::uint64_t seed_;
    const Clock::time_point started_ = Clock::now();

    std::vector<bool> busy_;      // n some thread is working on
    std::vector<double> shrink_;  // how much the next squeeze of each n tries to take off
    std::deque<Task> urgent_;     // neighbour tries set off by an improvement
    int cursor_ = 0;              // the next n this pass visits
    int pass_ = 1;

    // The smallest n worth working on: n = 1 is already exact.
    int firstN() const {
        return std::max(cfg_.minN, 2);
    }

    // The next task: urgent neighbour tries first, otherwise the next n of the pass, going
    // down so that every n sees n+1 freshly visited. An n another thread is still on is left
    // for the next pass.
    Task next() {
        std::lock_guard lock(mutex);

        while (!urgent_.empty()) {
            const Task t = urgent_.front();
            urgent_.pop_front();
            if (!busy_[t.n]) {
                busy_[t.n] = true;
                return t;
            }
        }

        for (int tried = 0; tried <= cfg_.maxN - firstN(); ++tried) {
            const int n = cursor_;
            const int pass = pass_;
            advance();
            if (busy_[n]) continue;

            busy_[n] = true;
            return {n, Op::Visit, pass};
        }

        // More threads than n: share one.
        busy_[firstN()] = true;
        return {firstN(), Op::Visit, pass_};
    }

    void advance() {
        if (--cursor_ >= firstN()) return;

        std::printf("[%5.0fs] pass %d done, score %.6f\n", elapsed(), pass_, totalScore(best_));
        std::fflush(stdout);
        cursor_ = cfg_.maxN;
        ++pass_;
    }

    Layout bestCopy(int n) {
        std::lock_guard lock(mutex);
        return best_[n];
    }

    std::mt19937_64 rngFor(const Task &task, int step) const {
        return std::mt19937_64(mix(seed_, static_cast<std::uint64_t>(task.pass),
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

    std::optional<Layout> squeezeInto(const Layout &layout, double target, std::mt19937_64 &rng,
                                      long itersPerTree = 0) {
        SqueezeConfig sc;
        if (itersPerTree == 0) itersPerTree = cfg_.squeezeItersPerTree;
        sc.iterations = itersPerTree * static_cast<long>(layout.size());
        sc.relaxIters = cfg_.squeezeRelaxIters;
        sc.stop = &stop_;

        std::optional<Layout> l = squeeze(layout, target, sc, rng);
        if (l) rotateWhole(*l);
        return l;
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

    // Every operator on n once, in order.
    void visit(const Task &task) {
        const int n = task.n;
        int step = 0;

        // n+1 minus each of its best few trees to remove, squeezed into the best square so far.
        const bool hasBiggerN = n + 1 < static_cast<int>(best_.size());
        if (hasBiggerN) {
            const Layout from = bestCopy(n + 1);
            const std::vector<std::size_t> order = from.removalOrder();
            const int tries = std::min<int>(cfg_.removeCandidates, static_cast<int>(order.size()));
            for (int k = 0; k < tries && !stop_.load(); ++k) {
                std::mt19937_64 rng = rngFor(task, step++);
                Layout l = from;
                l.remove(order[k]);
                offerIf(n, Op::Remove, squeezeInto(l, squeezeTarget(n), rng));
            }
        }

        // n-1 plus a tree, squeezed into the best square so far.
        if (n - 1 >= 1 && !stop_.load()) {
            std::mt19937_64 rng = rngFor(task, step++);
            offerIf(n, Op::Add, addAndSqueeze(n, rng));
        }

        // Squeeze the best layout until it fails often enough.
        for (int failures = 0; failures < cfg_.squeezeFailures && !stop_.load();) {
            std::mt19937_64 rng = rngFor(task, step++);
            const bool worked =
                offerIf(n, Op::Squeeze, squeezeInto(bestCopy(n), squeezeTarget(n), rng));
            adaptShrink(n, worked);
            failures = worked ? 0 : failures + 1;
        }

        // Shake and settle: a walk through arrangements as good as the best, each one a fresh
        // start for a squeeze into a smaller square.
        Layout current = bestCopy(n);
        for (int k = 0; k < cfg_.shakeTries && !stop_.load(); ++k) {
            std::mt19937_64 rng = rngFor(task, step++);
            const Layout kicked = shake(current, 1 + k % cfg_.shakeKicks, rng);
            std::optional<Layout> settled =
                squeezeInto(kicked, bestSide(n), rng, cfg_.shakeItersPerTree);
            if (!settled) continue;

            current = std::move(*settled);
            const bool worked = offer(n, Op::Shake, current) ||
                                offerIf(n, Op::Shake, squeezeInto(current, squeezeTarget(n), rng));
            if (worked) current = bestCopy(n);
        }

        if (!stop_.load()) {
            std::mt19937_64 rng = rngFor(task, step++);
            Layout l = bestCopy(n);
            anneal(l, cfg_.refineItersPerTree, cfg_.refineTemp, cfg_.refineStep, rng);
            offer(n, Op::Refine, std::move(l));
        }

        const bool randomPass = n <= cfg_.smallMax && task.pass % 2 == 1;
        if (randomPass && !stop_.load()) {
            std::mt19937_64 rng = rngFor(task, step++);
            Layout l = randomStart(n, rng);
            anneal(l, cfg_.randomItersPerTree, cfg_.randomTemp, cfg_.randomStep, rng);
            offer(n, Op::Random, std::move(l));
        }
    }

    // One neighbour try.
    void single(const Task &task) {
        const int n = task.n;
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

    // best[n-1] plus a tree at a random spot in its square, overlap and all, squeezed into a
    // square a bit smaller than best[n]'s.
    std::optional<Layout> addAndSqueeze(int n, std::mt19937_64 &rng) {
        Layout l = bestCopy(n - 1);
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        const Box b = l.bounds();
        l.add({b.minX + unit(rng) * (b.maxX - b.minX), b.minY + unit(rng) * (b.maxY - b.minY),
               360.0 * unit(rng)});

        return squeezeInto(l, squeezeTarget(n), rng);
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

        std::printf("[%5.0fs] n=%-3d %.6f -> %.6f  (%s)  score %.6f\n", elapsed(), n,
                    before * before / n, after * after / n, opName(op), totalScore(best_));
        std::fflush(stdout);

        // A smaller n can now come from this one minus a tree, and a bigger one from it plus a
        // tree. The pass is going down, so it reaches n-1 soon anyway: only n+1 needs a nudge.
        if (n + 1 <= cfg_.maxN) urgent_.push_back({n + 1, Op::Add, pass_});

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
            lastSave = now;
        }
    }

    for (std::thread &w : workers) {
        w.join();
    }

    save(best, origins);
}
} // namespace packer::packing
