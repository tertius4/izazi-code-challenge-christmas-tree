#include "packing/solver.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <random>
#include <utility>

#include "packing/lattice.h"
#include "packing/starts.h"

namespace packer::packing {
namespace {
constexpr double kInf = std::numeric_limits<double>::infinity();

struct Start {
    Layout layout;
    std::string origin;
};

// The strategy bands: which starting layouts get annealed for this n.
std::vector<Start> startsFor(int n, const Polygon &shape,
                             const std::vector<LatticeParams> &lattices, const SolverConfig &cfg,
                             std::mt19937_64 &rng) {
    std::vector<Start> starts;
    if (n == 1) {
        starts.push_back({bestSingle(shape), "exact rotation"});
        return starts;
    }

    std::vector<Start> grids;
    for (const LatticeParams &p : lattices) {
        grids.push_back({buildLattice(shape, p, n),
                         "lattice " + std::to_string(static_cast<int>(p.baseDeg)) + "deg"});
    }
    std::sort(grids.begin(), grids.end(),
              [](const Start &a, const Start &b) { return a.layout.side() < b.layout.side(); });

    if (n <= cfg.smallMax) {
        for (int r = 0; r < cfg.restarts; ++r) {
            starts.push_back({randomStart(shape, n, rng), "random"});
        }

        starts.push_back(std::move(grids.front()));
    } else if (n <= cfg.mediumMax) {
        for (Start &g : grids) {
            starts.push_back(std::move(g));
        }

        for (int r = 0; r < 2; ++r) {
            starts.push_back({randomStart(shape, n, rng), "random"});
        }
    } else {
        for (std::size_t i = 0; i < grids.size() && i < 2; ++i) {
            starts.push_back(std::move(grids[i]));
        }
    }
    return starts;
}

// The best (n-1)-tree layout we can get by deleting one tree.
Layout bestRemoval(const Layout &layout) {
    std::size_t bestIndex = 0;
    double bestSide = kInf;
    for (std::size_t skip = 0; skip < layout.size(); ++skip) {
        Box b{kInf, kInf, -kInf, -kInf};
        for (std::size_t j = 0; j < layout.size(); ++j) {
            if (j == skip) continue;
            const Box &t = layout.boxes[j];
            b.minX = std::min(b.minX, t.minX);
            b.minY = std::min(b.minY, t.minY);
            b.maxX = std::max(b.maxX, t.maxX);
            b.maxY = std::max(b.maxY, t.maxY);
        }
        const double side = std::max(b.maxX - b.minX, b.maxY - b.minY);
        if (side < bestSide) {
            bestSide = side;
            bestIndex = skip;
        }
    }
    Layout result = layout;
    result.remove(bestIndex);
    return result;
}
} // namespace

Solution solveAll(const Polygon &shape, const SolverConfig &cfg) {
    // Lattice search is the same for every n, so do it once per angle up front.
    std::vector<LatticeParams> lattices;
    for (double deg : cfg.latticeAngles) {
        lattices.push_back(findLattice(shape, deg));
    }

    Solution sol;
    sol.layouts.assign(cfg.maxN + 1, Layout(shape));
    sol.origins.assign(cfg.maxN + 1, "");

    // Each n is independent: spread them over all cores, biggest (slowest) first.
#pragma omp parallel for schedule(dynamic)
    for (int n = cfg.maxN; n >= 1; --n) {
        std::mt19937_64 rng(cfg.seed + n);
        AnnealConfig ac = cfg.anneal;
        ac.iterations = cfg.itersPerTree * n;

        double bestSide = kInf;
        for (Start &start : startsFor(n, shape, lattices, cfg, rng)) {
            if (!start.layout.isValid()) continue; // never anneal from a bad start

            anneal(start.layout, ac, rng);
            if (start.layout.isValid() && start.layout.side() < bestSide) {
                bestSide = start.layout.side();
                sol.layouts[n] = std::move(start.layout);
                sol.origins[n] = start.origin;
            }
        }

        if (static_cast<int>(sol.layouts[n].size()) != n) { // fallback, should not happen
            Layout layout = randomStart(shape, n, rng);
            anneal(layout, ac, rng);
            sol.layouts[n] = std::move(layout);
            sol.origins[n] = "random (fallback)";
        }

#pragma omp critical
        std::cout << "n=" << n << " side=" << sol.layouts[n].side() << " (" << sol.origins[n]
                  << ")\n";
    }

    // Backward pass: a good layout for n trees minus its worst tree can beat our n-1.
    for (int n = cfg.maxN; n >= 2; --n) {
        Layout candidate = bestRemoval(sol.layouts[n]);
        if (candidate.side() < sol.layouts[n - 1].side()) {
            sol.layouts[n - 1] = std::move(candidate);
            sol.origins[n - 1] = "removed 1 tree from n=" + std::to_string(n);
        }
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
