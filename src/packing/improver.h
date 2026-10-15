#pragma once

#include <atomic>
#include <functional>
#include <string>
#include <vector>

#include "packing/layout.h"

namespace packer::packing {
// Settings for improving a set of best-known layouts over a time budget.
struct ImproveConfig {
    double minutes = 10; // time budget
    int threads = 0;     // 0: half the cores, so the computer stays usable
    int minN = 1;        // only work on n in [minN, maxN]
    int maxN = 200;
    double saveSeconds = 120; // save this often, so an interrupted run loses little
    unsigned seed = 0;        // 0: from the clock, so every run tries something new
    // Small n are worth fresh random starts, annealed hot and long.
    int smallMax = 30;
    long randomItersPerTree = 20000;
    double randomTemp = 0.05;
    double randomStep = 0.2;
    // Every n: re-anneal the best layout cold, with a new seed.
    long refineItersPerTree = 20000;
    double refineTemp = 1e-3;
    double refineStep = 0.01;
    // Neighbours: n+1 minus a tree or n-1 plus a tree, then a short cold anneal.
    long neighbourItersPerTree = 2000;
};

// best[n] holds the best known n-tree layout (index 0 unused); origins[n] says where it came
// from. Both are updated in place: a layout is only ever replaced by a valid, smaller one.
// `save` is called every cfg.saveSeconds and at the end. Stops at the time budget or when
// `stop` turns true (for example on Ctrl-C).
using SaveFn =
    std::function<void(const std::vector<Layout> &best, const std::vector<std::string> &origins)>;
void improve(std::vector<Layout> &best, std::vector<std::string> &origins, const ImproveConfig &cfg,
             std::atomic<bool> &stop, const SaveFn &save);

// Kaggle score of a set of layouts: sum over n of side^2 / n.
double totalScore(const std::vector<Layout> &layouts);
} // namespace packer::packing
