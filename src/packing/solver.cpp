#include "packing/solver.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <random>
#include <sstream>
#include <utility>

#include "packing/improver.h"
#include "packing/lattice.h"
#include "packing/rotation.h"
#include "packing/starts.h"

namespace packer::packing {
namespace {
constexpr double kInf = std::numeric_limits<double>::infinity();

// A short label for the report, such as "lattice 12.3deg".
std::string latticeName(const LatticeParams &p) {
    std::ostringstream out;
    out << "lattice " << std::fixed << std::setprecision(1) << geometry::normalizeDeg(p.deg)
        << "deg";

    return out.str();
}

// True if the two lattices are the same pattern in the same orientation, up to a mirror
// image or a half turn. Such copies cut into exactly the same squares.
bool sameCut(const LatticeParams &p, const LatticeParams &q) {
    auto fold = [](double deg) {
        const double d = std::fmod(geometry::normalizeDeg(deg), 180.0);
        return std::min(d, 180.0 - d);
    };

    const bool sameAngle = std::abs(fold(p.deg) - fold(q.deg)) < 0.5;
    const bool sameA = std::abs(p.a - q.a) < 1e-3;
    const bool sameB = std::abs(p.b - q.b) < 1e-3;

    return sameAngle && sameA && sameB;
}

// Searches for dense lattices from random cells and keeps the `latticeKeep` densest distinct
// ones. The searches end in different local optima, and the densest pattern is often not the
// one that cuts best into a square, so keeping several pays off.
std::vector<LatticeParams> findLattices(const SolverConfig &cfg) {
    std::vector<std::optional<LatticeParams>> found(cfg.latticeStarts);
#pragma omp parallel for schedule(dynamic)
    for (int k = 0; k < cfg.latticeStarts; ++k) {
        std::mt19937_64 rng(cfg.seed * 1000003ULL + k);
        found[k] = findLattice(cfg.latticeIters, rng);
    }

    std::vector<LatticeParams> sorted;
    for (const auto &p : found) {
        if (p) sorted.push_back(*p);
    }

    std::sort(sorted.begin(), sorted.end(), [](const LatticeParams &a, const LatticeParams &b) {
        return areaPerTree(a) < areaPerTree(b);
    });

    std::vector<LatticeParams> lattices;
    for (const LatticeParams &p : sorted) {
        const bool haveEnough = static_cast<int>(lattices.size()) == cfg.latticeKeep;
        if (haveEnough) break;

        auto sameAsP = [&](const LatticeParams &q) { return sameCut(p, q); };
        const bool isNew = std::none_of(lattices.begin(), lattices.end(), sameAsP);
        if (isNew) {
            lattices.push_back(p);
        }
    }

    for (const LatticeParams &p : lattices) {
        std::cout << latticeName(p) << ": area per tree " << areaPerTree(p) << "\n";
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
    AnnealConfig random = cfg.anneal;
    random.iterations = cfg.randomItersPerTree * n;
    random.startTemp = cfg.randomTemp;
    random.startStep = cfg.randomStep;

    // A lattice cut is already dense: refine it, gently (cold) or with the usual heat.
    AnnealConfig lattice = cfg.anneal;
    lattice.iterations = (n > cfg.mediumMax  ? cfg.largeItersPerTree
                          : n > cfg.smallMax ? cfg.mediumItersPerTree
                                             : cfg.itersPerTree) *
                         n;

    AnnealConfig cold = lattice;
    cold.startTemp = cfg.latticeTemp;
    cold.startStep = cfg.latticeStep;

    std::vector<Start> grids;
    for (const LatticeParams &p : lattices) {
        grids.push_back({cutLattice(p, n), latticeName(p), lattice});
    }

    std::sort(grids.begin(), grids.end(),
              [](const Start &a, const Start &b) { return a.layout.side() < b.layout.side(); });

    auto addRandom = [&](int count) {
        for (int r = 0; r < count; ++r) {
            starts.push_back({randomStart(n, rng), "random", random});
        }
    };

    if (n <= cfg.smallMax) {
        addRandom(cfg.restarts);

        if (!grids.empty()) {
            starts.push_back({grids.front().layout, grids.front().origin + " (cold)", cold});
            starts.push_back(std::move(grids.front()));
        }
    } else if (n <= cfg.mediumMax) {
        if (!grids.empty()) {
            starts.push_back({grids.front().layout, grids.front().origin + " (cold)", cold});
        }

        const int cuts = std::min(static_cast<int>(grids.size()), cfg.mediumCuts);
        for (int i = 0; i < cuts; ++i) {
            starts.push_back(std::move(grids[i]));
        }

        addRandom(cfg.mediumRandom);
    } else {
        // Large n: the cut is already close to the best we know, so refine it gently.
        const int cuts = std::min(static_cast<int>(grids.size()), cfg.largeCuts);
        for (int i = 0; i < cuts; ++i) {
            grids[i].anneal = cold;
            starts.push_back(std::move(grids[i]));
        }
    }

    return starts;
}

// A short, cold anneal of a layout that is already good, then the best whole-layout turn.
void refine(Layout &layout, const SolverConfig &cfg, std::mt19937_64 &rng) {
    AnnealConfig ac = cfg.anneal;
    ac.iterations = cfg.polishItersPerTree * static_cast<long>(layout.size());
    ac.startTemp = cfg.polishTemp;
    ac.startStep = cfg.polishStep;

    anneal(layout, ac, rng);
    rotateWhole(layout);
}

} // namespace

Solution solveAll(const SolverConfig &cfg) {
    const auto started = std::chrono::steady_clock::now();
    auto elapsed = [started] { // whole seconds since the start, for the progress log
        return std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() -
                                                                started)
            .count();
    };

    // Lattice search is the same for every n, so do it once up front.
    const std::vector<LatticeParams> lattices = findLattices(cfg);
    std::cout << "[" << elapsed() << "s] lattices found" << std::endl;

    Solution sol;
    sol.layouts.assign(cfg.maxN + 1, Layout());
    sol.origins.assign(cfg.maxN + 1, "");

    // Each n is independent: spread them over all cores, slowest first so no core is left
    // with a long job at the end. The small n, with their many random starts, are slowest.
    std::vector<int> order;
    for (int n = std::min(cfg.smallMax, cfg.maxN); n >= 1; --n) {
        order.push_back(n);
    }

    for (int n = cfg.maxN; n > cfg.smallMax; --n) {
        order.push_back(n);
    }

#pragma omp parallel for schedule(dynamic)
    for (std::size_t k = 0; k < order.size(); ++k) {
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
                sol.layouts[n] = std::move(start.layout);
                sol.origins[n] = start.origin;
            }
        }

        const bool hasAllTrees = static_cast<int>(sol.layouts[n].size()) == n;
        if (n > 1 && hasAllTrees) {
            rotateWhole(sol.layouts[n]);
        }

        if (!hasAllTrees) { // fallback, should not happen
            Layout layout = randomStart(n, rng);
            AnnealConfig ac = cfg.anneal;
            ac.iterations = cfg.randomItersPerTree * n;
            anneal(layout, ac, rng);
            sol.layouts[n] = std::move(layout);
            sol.origins[n] = "random (fallback)";
        }

#pragma omp critical
        std::cout << "[" << elapsed() << "s] n=" << n << " side=" << sol.layouts[n].side() << " ("
                  << sol.origins[n] << ")" << std::endl;
    }

    // Neighbouring n help each other. Backward: n trees minus the worst tree can beat n-1.
    // Forward: n-1 trees plus a tree in a free spot can beat n. Each candidate is polished
    // before it replaces anything. Every pass reads a snapshot, so the n run in parallel.
    constexpr double kTol = 1e-12;
    for (int pass = 0; pass < cfg.passes; ++pass) {
        const std::vector<Layout> snapshot = sol.layouts;
        std::vector<std::optional<Layout>> better(cfg.maxN + 1);
        std::vector<std::string> how(cfg.maxN + 1);

        // Keeps a candidate layout for m trees if it beats the one kept so far.
        auto offer = [&](int m, Layout candidate, const std::string &origin) {
#pragma omp critical(better)
            {
                std::optional<Layout> &kept = better[m];
                const bool beatsKept = !kept || candidate.side() < kept->side();
                if (beatsKept) {
                    kept = std::move(candidate);
                    how[m] = origin;
                }
            }
        };

#pragma omp parallel for schedule(dynamic)
        for (int n = cfg.maxN; n >= 2; --n) {
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
        }

        bool changed = false;
        for (int n = 1; n <= cfg.maxN; ++n) {
            std::optional<Layout> &found = better[n];
            if (!found) continue;

            const bool smaller = found->side() < sol.layouts[n].side() - kTol;
            if (smaller) {
                sol.layouts[n] = std::move(*found);
                sol.origins[n] = how[n];
                changed = true;
            }
        }

        std::cout << "[" << elapsed() << "s] pass " << pass << ": score " << totalScore(sol.layouts)
                  << (changed ? "" : " (no change)") << std::endl;

        if (!changed) break;
    }

    return sol;
}
} // namespace packer::packing
