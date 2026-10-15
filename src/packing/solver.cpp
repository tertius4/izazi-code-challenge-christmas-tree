#include "packing/solver.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <optional>
#include <random>
#include <sstream>
#include <utility>

#include "packing/rotation.h"
#include "packing/lattice.h"
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
    return std::abs(fold(p.deg) - fold(q.deg)) < 0.5 && std::abs(p.a - q.a) < 1e-3 &&
           std::abs(p.b - q.b) < 1e-3;
}
} // namespace

// The searches end in different local optima, and the densest pattern is often not the one
// that cuts best into a square, so keeping several pays off.
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
        if (static_cast<int>(lattices.size()) == cfg.latticeKeep) break;
        if (std::none_of(lattices.begin(), lattices.end(),
                         [&](const LatticeParams &q) { return sameCut(p, q); })) {
            lattices.push_back(p);
        }
    }

    for (const LatticeParams &p : lattices) {
        std::cout << latticeName(p) << ": area per tree " << areaPerTree(p) << "\n";
    }
    return lattices;
}

namespace {
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
    lattice.iterations = (n > cfg.mediumMax   ? cfg.largeItersPerTree
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
        for (std::size_t i = 0; i < grids.size() && static_cast<int>(i) < cfg.mediumCuts; ++i) {
            starts.push_back(std::move(grids[i]));
        }
        addRandom(cfg.mediumRandom);
    } else {
        // Large n: the cut is already close to the best we know, so refine it gently.
        for (std::size_t i = 0; i < grids.size() && static_cast<int>(i) < cfg.largeCuts; ++i) {
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
    for (int n = std::min(cfg.smallMax, cfg.maxN); n >= 1; --n) order.push_back(n);
    for (int n = cfg.maxN; n > cfg.smallMax; --n) order.push_back(n);

#pragma omp parallel for schedule(dynamic)
    for (std::size_t k = 0; k < order.size(); ++k) {
        const int n = order[k];
        std::mt19937_64 rng(cfg.seed + n);

        double bestSide = kInf;
        for (Start &start : startsFor(n, lattices, cfg, rng)) {
            if (!start.layout.isValid()) continue; // never anneal from a bad start

            anneal(start.layout, start.anneal, rng);
            if (start.layout.isValid() && start.layout.side() < bestSide) {
                bestSide = start.layout.side();
                sol.layouts[n] = std::move(start.layout);
                sol.origins[n] = start.origin;
            }
        }

        if (n > 1 && static_cast<int>(sol.layouts[n].size()) == n) {
            rotateWhole(sol.layouts[n]);
        }

        if (static_cast<int>(sol.layouts[n].size()) != n) { // fallback, should not happen
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

#pragma omp parallel for schedule(dynamic)
        for (int n = cfg.maxN; n >= 2; --n) {
            std::mt19937_64 rng(cfg.seed * 7919ULL + n * 131ULL + pass);
            const double current = snapshot[n].side();

            // Backward: build the (n-1)-tree layout, so it is judged against n-1.
            {
                Layout candidate = withoutTree(snapshot[n]);
                const double target = snapshot[n - 1].side();
                // Worth a polish if it is already close: polishing usually gains ~1%.
                if (candidate.side() < target * 1.02) {
                    refine(candidate, cfg, rng);
                    if (candidate.isValid() && candidate.side() < target - kTol) {
#pragma omp critical(better)
                        if (!better[n - 1] || candidate.side() < better[n - 1]->side()) {
                            better[n - 1] = std::move(candidate);
                            how[n - 1] = "removed 1 tree from n=" + std::to_string(n);
                        }
                    }
                }
            }

            // Forward: add a tree to the (n-1)-tree layout and compare with n.
            if (auto grown = withExtraTree(snapshot[n - 1], rng)) {
                refine(*grown, cfg, rng);
                if (grown->isValid() && grown->side() < current - kTol) {
#pragma omp critical(better)
                    if (!better[n] || grown->side() < better[n]->side()) {
                        better[n] = std::move(*grown);
                        how[n] = "added 1 tree to n=" + std::to_string(n - 1);
                    }
                }
            }
        }

        bool changed = false;
        for (int n = 1; n <= cfg.maxN; ++n) {
            if (better[n] && better[n]->side() < sol.layouts[n].side() - kTol) {
                sol.layouts[n] = std::move(*better[n]);
                sol.origins[n] = how[n];
                changed = true;
            }
        }
        std::cout << "[" << elapsed() << "s] pass " << pass << ": score " << score(sol)
                  << (changed ? "" : " (no change)") << std::endl;
        if (!changed) break;
    }
    return sol;
}

double score(const Solution &solution) {
    double total = 0;
    for (std::size_t n = 1; n < solution.layouts.size(); ++n) {
        const double s = solution.layouts[n].side();
        total += s * s / static_cast<double>(n);
    }
    return total;
}
} // namespace packer::packing
