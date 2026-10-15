#include "packing/improver.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <deque>
#include <mutex>
#include <random>
#include <thread>

#include "packing/annealer.h"
#include "packing/rotation.h"
#include "packing/starts.h"

namespace packer::packing {
namespace {
using Clock = std::chrono::steady_clock;

enum class Op { Random, Refine, Remove, Add };

const char *opName(Op op) {
    switch (op) {
    case Op::Random: return "random start";
    case Op::Refine: return "refine";
    case Op::Remove: return "from n+1 minus a tree";
    case Op::Add: return "from n-1 plus a tree";
    }
    return "?";
}

struct Task {
    int n = 0;
    Op op = Op::Refine;
};

// Shared state of a run. Everything below the mutex is guarded by it.
class Run {
  public:
    Run(std::vector<Layout> &best, std::vector<std::string> &origins, const ImproveConfig &cfg,
        std::atomic<bool> &stop)
        : best_(best), origins_(origins), cfg_(cfg), stop_(stop), started_(Clock::now()),
          spent_(best.size(), 0.0), wins_(best.size(), 0), busy_(best.size(), false) {
    }

    void work(unsigned seed) {
        std::mt19937_64 rng(seed);
        while (!stop_.load()) {
            const Task task = next(rng);
            const auto t0 = Clock::now();
            std::optional<Layout> result = attempt(task, rng);
            finish(task, std::chrono::duration<double>(Clock::now() - t0).count(),
                   std::move(result));
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
    const Clock::time_point started_;

    std::vector<double> spent_; // seconds spent on each n this run
    std::vector<int> wins_;     // improvements of each n this run
    std::vector<bool> busy_;    // n some thread is working on
    std::deque<Task> urgent_;   // neighbour tasks set off by an improvement

    // The next task: urgent neighbour tasks first; otherwise the n with the least time spent
    // for its weight (n that improved this run get more), with a band-dependent operator.
    Task next(std::mt19937_64 &rng) {
        std::lock_guard lock(mutex);
        while (!urgent_.empty()) {
            const Task t = urgent_.front();
            urgent_.pop_front();
            if (!busy_[t.n]) {
                busy_[t.n] = true;
                return t;
            }
        }
        int pick = -1;
        double pickCost = INFINITY;
        for (int n = std::max(cfg_.minN, 2); n <= cfg_.maxN; ++n) {
            const double cost = spent_[n] / (1.0 + 3.0 * wins_[n]);
            if (!busy_[n] && cost < pickCost) {
                pickCost = cost;
                pick = n;
            }
        }
        if (pick < 0) pick = std::max(cfg_.minN, 2); // every n is busy: share one
        busy_[pick] = true;

        std::uniform_real_distribution<double> unit(0.0, 1.0);
        const double r = unit(rng);
        Op op = Op::Refine;
        if (pick <= cfg_.smallMax && r < 0.5) {
            op = Op::Random;
        } else if (r > 0.8) {
            op = unit(rng) < 0.5 ? Op::Remove : Op::Add;
        }
        return {pick, op};
    }

    Layout bestCopy(int n) {
        std::lock_guard lock(mutex);
        return best_[n];
    }

    void anneal(Layout &layout, long itersPerTree, double temp, double step,
                std::mt19937_64 &rng) {
        AnnealConfig ac;
        ac.iterations = itersPerTree * static_cast<long>(layout.size());
        ac.startTemp = temp;
        ac.startStep = step;
        ac.stop = &stop_;
        packing::anneal(layout, ac, rng);
        rotateWhole(layout);
    }

    std::optional<Layout> attempt(const Task &task, std::mt19937_64 &rng) {
        const int n = task.n;
        switch (task.op) {
        case Op::Random: {
            Layout l = randomStart(n, rng);
            anneal(l, cfg_.randomItersPerTree, cfg_.randomTemp, cfg_.randomStep, rng);
            return l;
        }
        case Op::Refine: {
            Layout l = bestCopy(n);
            anneal(l, cfg_.refineItersPerTree, cfg_.refineTemp, cfg_.refineStep, rng);
            return l;
        }
        case Op::Remove: {
            if (n + 1 >= static_cast<int>(best_.size())) return std::nullopt;
            Layout l = withoutTree(bestCopy(n + 1));
            anneal(l, cfg_.neighbourItersPerTree, cfg_.refineTemp, cfg_.refineStep, rng);
            return l;
        }
        case Op::Add: {
            auto l = withExtraTree(bestCopy(n - 1), rng);
            if (!l) return std::nullopt;
            anneal(*l, cfg_.neighbourItersPerTree, cfg_.refineTemp, cfg_.refineStep, rng);
            return l;
        }
        }
        return std::nullopt;
    }

    void finish(const Task &task, double seconds, std::optional<Layout> result) {
        const int n = task.n;
        // Check validity outside the lock; it is the slowest part of accepting.
        const bool usable = result && static_cast<int>(result->size()) == n && result->isValid();

        std::lock_guard lock(mutex);
        spent_[n] += seconds;
        busy_[n] = false;
        if (!usable) return;
        const double before = best_[n].side();
        const double after = result->side();
        if (!(after < before - 1e-12)) return;

        best_[n] = std::move(*result);
        origins_[n] = opName(task.op);
        ++wins_[n];
        std::printf("[%5.0fs] n=%-3d %.6f -> %.6f  (%s)  score %.6f\n", elapsed(), n,
                    before * before / n, after * after / n, opName(task.op), totalScore(best_));
        std::fflush(stdout);
        // A smaller n can now come from this one minus a tree, and a bigger one from it plus a
        // tree: try those soon.
        if (n - 1 >= std::max(cfg_.minN, 2)) urgent_.push_back({n - 1, Op::Remove});
        if (n + 1 <= cfg_.maxN) urgent_.push_back({n + 1, Op::Add});
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

void improve(std::vector<Layout> &best, std::vector<std::string> &origins,
             const ImproveConfig &cfg, std::atomic<bool> &stop, const SaveFn &save) {
    const int threads = cfg.threads > 0
                            ? cfg.threads
                            : std::max(1u, std::thread::hardware_concurrency() / 2);
    const unsigned seed = cfg.seed != 0 ? cfg.seed
                                        : static_cast<unsigned>(
                                              Clock::now().time_since_epoch().count());
    std::printf("Improving n=%d..%d for %.1f min on %d threads (seed %u), from score %.6f\n",
                cfg.minN, cfg.maxN, cfg.minutes, threads, seed, totalScore(best));
    std::fflush(stdout);

    Run run(best, origins, cfg, stop);
    std::vector<std::thread> workers;
    for (int t = 0; t < threads; ++t) {
        workers.emplace_back([&run, seed, t] { run.work(seed + 7919u * static_cast<unsigned>(t)); });
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
    for (std::thread &w : workers) w.join();
    save(best, origins);
}
} // namespace packer::packing
