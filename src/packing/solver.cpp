#include "packing/solver.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <mutex>
#include <optional>
#include <random>
#include <sstream>
#include <utility>

#include "packing/lattice.h"
#include "packing/layout.h"
#include "packing/rotation.h"
#include "packing/starts.h"
#include "util/parallel.h"

namespace packer::packing {
namespace {
using geometry::kInf;

// A short label for the report, such as "lattice 12.3deg".
std::string latticeName(const LatticeParams &lattice) {
    std::ostringstream out;
    out << "lattice " << std::fixed << std::setprecision(1) << geometry::normalizeDeg(lattice.deg)
        << "deg";

    return out.str();
}

// True if the two lattices are the same pattern in the same orientation, up to a mirror
// image or a half turn. Such copies cut into exactly the same squares.
bool sameCut(const LatticeParams &a, const LatticeParams &b) {
    auto fold = [](double deg) {
        const double halfTurn = std::fmod(geometry::normalizeDeg(deg), 180.0);
        return std::min(halfTurn, 180.0 - halfTurn);
    };

    const bool sameAngle = std::abs(fold(a.deg) - fold(b.deg)) < 0.5;
    const bool sameA = std::abs(a.a - b.a) < 1e-3;
    const bool sameB = std::abs(a.b - b.b) < 1e-3;

    return sameAngle && sameA && sameB;
}

// Searches for dense lattices from random cells and keeps the `latticeKeep` densest distinct
// ones. The searches end in different local optima, and the densest pattern is often not the
// one that cuts best into a square, so keeping several pays off.
std::vector<LatticeParams> findLattices(const SolverConfig &cfg) {
    std::vector<std::optional<LatticeParams>> found(cfg.latticeStarts);
    util::parallelFor(found.size(), cfg.threads, [&](std::size_t k) {
        std::mt19937_64 rng(cfg.seed * 1000003ULL + k);
        found[k] = findLattice(cfg.latticeIters, rng);
    });

    std::vector<LatticeParams> sorted;
    for (const std::optional<LatticeParams> &lattice : found) {
        if (lattice) sorted.push_back(*lattice);
    }

    std::sort(sorted.begin(), sorted.end(), [](const LatticeParams &a, const LatticeParams &b) {
        return areaPerTree(a) < areaPerTree(b);
    });

    std::vector<LatticeParams> lattices;
    for (const LatticeParams &lattice : sorted) {
        const bool haveEnough = static_cast<int>(lattices.size()) == cfg.latticeKeep;
        if (haveEnough) break;

        auto sameAsThis = [&](const LatticeParams &kept) { return sameCut(lattice, kept); };
        const bool isNew = std::none_of(lattices.begin(), lattices.end(), sameAsThis);
        if (isNew) {
            lattices.push_back(lattice);
        }
    }

    for (const LatticeParams &lattice : lattices) {
        std::printf("%s: area per tree %g\n", latticeName(lattice).c_str(), areaPerTree(lattice));
    }

    return lattices;
}

struct Start {
    Layout layout;
    std::string origin;
    AnnealConfig anneal; // how hard to anneal this start
};

// The strategy bands: which starting layouts get annealed for this n, and how.
std::vector<Start> startsFor(int n, const std::vector<LatticeParams> &lattices,
                             const SolverConfig &cfg, std::mt19937_64 &rng) {
    std::vector<Start> starts;
    if (n == 1) {
        starts.push_back({bestSingle(), "exact rotation", {}});
        starts.back().anneal.iterations = 0; // already the best single tree
        return starts;
    }

    // Random starts are loose, so they need a hot, long anneal to pack at all.
    AnnealConfig hotAnneal = cfg.anneal;
    hotAnneal.iterations = cfg.randomItersPerTree * n;
    hotAnneal.startTemp = cfg.randomTemp;
    hotAnneal.startStep = cfg.randomStep;

    // A lattice cut is already dense: refine it, gently (cold) or with the usual heat.
    AnnealConfig cutAnneal = cfg.anneal;
    cutAnneal.iterations = (n > cfg.mediumMax  ? cfg.largeItersPerTree
                            : n > cfg.smallMax ? cfg.mediumItersPerTree
                                               : cfg.itersPerTree) *
                           n;

    AnnealConfig coldAnneal = cutAnneal;
    coldAnneal.startTemp = cfg.latticeTemp;
    coldAnneal.startStep = cfg.latticeStep;

    std::vector<Start> cuts;
    for (const LatticeParams &lattice : lattices) {
        cuts.push_back({cutLattice(lattice, n), latticeName(lattice), cutAnneal});
    }

    std::sort(cuts.begin(), cuts.end(),
              [](const Start &a, const Start &b) { return a.layout.side() < b.layout.side(); });

    auto addRandom = [&](int count) {
        for (int k = 0; k < count; ++k) {
            starts.push_back({randomStart(n, rng), "random", hotAnneal});
        }
    };

    if (n <= cfg.smallMax) {
        addRandom(cfg.restarts);

        if (!cuts.empty()) {
            starts.push_back({cuts.front().layout, cuts.front().origin + " (cold)", coldAnneal});
            starts.push_back(std::move(cuts.front()));
        }
    } else if (n <= cfg.mediumMax) {
        if (!cuts.empty()) {
            starts.push_back({cuts.front().layout, cuts.front().origin + " (cold)", coldAnneal});
        }

        const int count = std::min(static_cast<int>(cuts.size()), cfg.mediumCuts);
        for (int k = 0; k < count; ++k) {
            starts.push_back(std::move(cuts[k]));
        }

        addRandom(cfg.mediumRandom);
    } else {
        // Large n: the cut is already close to the best we know, so refine it gently.
        const int count = std::min(static_cast<int>(cuts.size()), cfg.largeCuts);
        for (int k = 0; k < count; ++k) {
            cuts[k].anneal = coldAnneal;
            starts.push_back(std::move(cuts[k]));
        }
    }

    return starts;
}

// A short, cold anneal of a layout that is already good, then the best whole-layout turn.
void refine(Layout &layout, const SolverConfig &cfg, std::mt19937_64 &rng) {
    AnnealConfig annealCfg = cfg.anneal;
    annealCfg.iterations = cfg.polishItersPerTree * static_cast<long>(layout.size());
    annealCfg.startTemp = cfg.polishTemp;
    annealCfg.startStep = cfg.polishStep;

    anneal(layout, annealCfg, rng);
    rotateWhole(layout);
}

} // namespace

Solution solveAll(const SolverConfig &cfg) {
    const auto started = std::chrono::steady_clock::now();
    auto elapsed = [started] { // whole seconds since the start, for the progress log
        const auto spent = std::chrono::steady_clock::now() - started;
        return static_cast<long>(std::chrono::duration_cast<std::chrono::seconds>(spent).count());
    };

    // Lattice search is the same for every n, so do it once up front.
    const std::vector<LatticeParams> lattices = findLattices(cfg);
    std::printf("[%lds] lattices found\n", elapsed());
    std::fflush(stdout);

    Solution solution;
    solution.layouts.assign(cfg.maxN + 1, Layout());
    solution.origins.assign(cfg.maxN + 1, "");

    // Each n is independent: spread them over all cores, slowest first so no core is left
    // with a long job at the end. The small n, with their many random starts, are slowest.
    std::vector<int> order;
    for (int n = std::min(cfg.smallMax, cfg.maxN); n >= 1; --n) {
        order.push_back(n);
    }

    for (int n = cfg.maxN; n > cfg.smallMax; --n) {
        order.push_back(n);
    }

    util::parallelFor(order.size(), cfg.threads, [&](std::size_t k) {
        const int n = order[k];
        std::mt19937_64 rng(cfg.seed + n);

        double bestSide = kInf;
        for (Start &start : startsFor(n, lattices, cfg, rng)) {
            // Never anneal from a bad start.
            if (!start.layout.isValid()) continue;

            anneal(start.layout, start.anneal, rng);
            const bool smaller = start.layout.side() < bestSide;
            const bool isBest = smaller && start.layout.isValid();
            if (isBest) {
                bestSide = start.layout.side();
                solution.layouts[n] = std::move(start.layout);
                solution.origins[n] = start.origin;
            }
        }

        const bool hasAllTrees = static_cast<int>(solution.layouts[n].size()) == n;
        if (n > 1 && hasAllTrees) {
            rotateWhole(solution.layouts[n]);
        }

        if (!hasAllTrees) { // fallback, should not happen
            Layout layout = randomStart(n, rng);
            AnnealConfig annealCfg = cfg.anneal;
            annealCfg.iterations = cfg.randomItersPerTree * n;
            anneal(layout, annealCfg, rng);
            solution.layouts[n] = std::move(layout);
            solution.origins[n] = "random (fallback)";
        }

        // One printf call writes its line whole, so lines from different threads don't mix.
        std::printf("[%lds] n=%d side=%g (%s)\n", elapsed(), n, solution.layouts[n].side(),
                    solution.origins[n].c_str());
        std::fflush(stdout);
    });

    // Neighbouring n help each other. Backward: n trees minus the worst tree can beat n-1.
    // Forward: n-1 trees plus a tree in a free spot can beat n. Each candidate is polished
    // before it replaces anything. Every pass reads a snapshot, so the n run in parallel.
    constexpr double kTol = 1e-12;
    for (int pass = 0; pass < cfg.passes; ++pass) {
        const std::vector<Layout> snapshot = solution.layouts;
        std::vector<std::optional<Layout>> better(cfg.maxN + 1);
        std::vector<std::string> betterOrigin(cfg.maxN + 1);

        // Keeps a candidate layout of `trees` trees if it beats the one kept so far.
        std::mutex betterMutex;
        auto offer = [&](int trees, Layout candidate, const std::string &origin) {
            const std::lock_guard lock(betterMutex);
            std::optional<Layout> &kept = better[trees];
            const bool beatsKept = !kept || candidate.side() < kept->side();
            if (beatsKept) {
                kept = std::move(candidate);
                betterOrigin[trees] = origin;
            }
        };

        // n from maxN down to 2.
        const std::size_t count = cfg.maxN >= 2 ? static_cast<std::size_t>(cfg.maxN - 1) : 0;
        util::parallelFor(count, cfg.threads, [&](std::size_t k) {
            const int n = cfg.maxN - static_cast<int>(k);
            std::mt19937_64 rng(cfg.seed * 7919ULL + n * 131ULL + pass);
            const double current = snapshot[n].side();

            // Backward: build the (n-1)-tree layout, so it is judged against n-1.
            {
                Layout candidate = withoutTree(snapshot[n]);
                const double target = snapshot[n - 1].side();
                // Worth a polish if it is already close: polishing usually gains ~1%.
                const bool closeEnough = candidate.side() < target * 1.02;
                if (closeEnough) {
                    refine(candidate, cfg, rng);

                    const bool smaller = candidate.side() < target - kTol;
                    const bool improves = smaller && candidate.isValid();
                    if (improves) {
                        const std::string origin = "removed 1 tree from n=" + std::to_string(n);
                        offer(n - 1, std::move(candidate), origin);
                    }
                }
            }

            // Forward: add a tree to the (n-1)-tree layout and compare with n.
            if (auto grown = withExtraTree(snapshot[n - 1], rng)) {
                refine(*grown, cfg, rng);

                const bool smaller = grown->side() < current - kTol;
                const bool improves = smaller && grown->isValid();
                if (improves) {
                    const std::string origin = "added 1 tree to n=" + std::to_string(n - 1);
                    offer(n, std::move(*grown), origin);
                }
            }
        });

        bool changed = false;
        for (int n = 1; n <= cfg.maxN; ++n) {
            std::optional<Layout> &found = better[n];
            if (!found) continue;

            const bool smaller = found->side() < solution.layouts[n].side() - kTol;
            if (smaller) {
                solution.layouts[n] = std::move(*found);
                solution.origins[n] = betterOrigin[n];
                changed = true;
            }
        }

        std::printf("[%lds] pass %d: score %g%s\n", elapsed(), pass, totalScore(solution.layouts),
                    changed ? "" : " (no change)");
        std::fflush(stdout);

        if (!changed) break;
    }

    return solution;
}
} // namespace packer::packing
