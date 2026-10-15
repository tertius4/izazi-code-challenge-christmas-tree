#pragma once

#include <string>
#include <vector>

#include "packing/annealer.h"
#include "packing/layout.h"

namespace packer::packing {
struct SolverConfig {
    int maxN = 200;
    int smallMax = 12;  // 2 → smallMax: many random starts + the best lattice cut
    int mediumMax = 30; // smallMax+1 → mediumMax: the best lattice cuts + some random starts
                        // > mediumMax: the best lattice cuts only
    // Random starts are loose, so they are annealed hot and long.
    int restarts = 100;              // random starts per small N
    int mediumRandom = 16;           // random starts per medium N
    long randomItersPerTree = 20000; // annealing moves per random start, per tree
    double randomTemp = 0.05;
    double randomStep = 0.2;
    // Lattice cuts are already dense, so they are refined.
    long itersPerTree = 2500;       // annealing moves per cut, per tree, for small and medium N
    long largeItersPerTree = 40000; // the same for large N: the main use of the time budget
    int mediumCuts = 4;             // how many of the smallest cuts are annealed for medium N...
    int largeCuts = 2;              // ...and for large N
    // Annealing from a lattice cut can begin cold, so the dense pattern is refined, not melted.
    double latticeTemp = 1e-3;
    double latticeStep = 0.01;
    // Lattice search, done once: each search anneals a random cell.
    int latticeStarts = 48;
    long latticeIters = 100000; // annealing moves per search
    int latticeKeep = 24;       // how many of the densest distinct lattices are cut for each N
    AnnealConfig anneal;
    // Whole-layout rotation runs on every n's best layout, followed by shrink-and-repair if
    // repairIters > 0 (repair moves per tree, per shrink attempt). It is off by default:
    // after annealing it found nothing, and it was the slowest part of the passes.
    long repairIters = 0;
    // Backward and forward passes: build n from n+1 minus a tree, or n-1 plus a tree, then
    // polish it with a short cold anneal. Repeats until nothing improves, up to `passes`.
    int passes = 10;
    long polishItersPerTree = 400;
    double polishTemp = 5e-4;
    double polishStep = 0.01;
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
