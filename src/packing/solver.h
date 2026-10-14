#pragma once

#include <string>
#include <vector>

#include "packing/annealer.h"
#include "packing/layout.h"

namespace packer::packing {
struct SolverConfig {
    int maxN = 200;
    int smallMax = 12;        // 2 → smallMax: many random starts
    int mediumMax = 60;       // smallMax+1 → mediumMax: all lattice angles + some random starts
                              // > mediumMax: the two best lattices only
    int restarts = 8;         // random starts per small N
    long itersPerTree = 2000; // annealing effort: iterations = itersPerTree * N per start
    std::vector<double> latticeAngles = {0, 15, 30, 45, 60, 75};
    AnnealConfig anneal;
    unsigned seed = 42;
};

struct Solution {
    std::vector<Layout> layouts;      // layouts[n] holds n trees; index 0 is unused
    std::vector<std::string> origins; // which strategy produced layouts[n]
};

Solution solveAll(const Polygon &shape, const SolverConfig &cfg);

// Kaggle score: sum over n of side^2 / n. Lower is better.
double score(const Solution &solution);
} // namespace packer::packing
