#pragma once

#include <string>
#include <vector>

#include "packing/annealer.h"
#include "packing/lattice.h"
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
    long itersPerTree = 2500;        // annealing moves per cut, per tree, for small N...
    long mediumItersPerTree = 40000; // ...medium N...
    long largeItersPerTree = 40000;  // ...and large N: the main use of the time budget
    int mediumCuts = 4;              // how many of the smallest cuts are annealed for medium N...
    int largeCuts = 2;               // ...and for large N
    // Annealing from a lattice cut can begin cold, so the dense pattern is refined, not melted.
    double latticeTemp = 1e-3;
    double latticeStep = 0.01;
    // Lattice search, done once: each search anneals a random cell.
    int latticeStarts = 48;
    long latticeIters = 100000; // annealing moves per search
    int latticeKeep = 24;       // how many of the densest distinct lattices are cut for each N
    AnnealConfig anneal;
    // Backward and forward passes: build n from n+1 minus a tree, or n-1 plus a tree, then
    // polish it with a short cold anneal. Repeats until nothing improves, up to `passes`.
    int passes = 10;
    long polishItersPerTree = 400;
    double polishTemp = 5e-4;
    double polishStep = 0.01;
    unsigned seed = 42;
    int threads = 0; // 0: every core
};

struct Solution {
    std::vector<Layout> layouts;      // layouts[n] holds n trees; index 0 is unused
    std::vector<std::string> origins; // which strategy produced layouts[n]
};

Solution solveAll(const SolverConfig &cfg);
} // namespace packer::packing
