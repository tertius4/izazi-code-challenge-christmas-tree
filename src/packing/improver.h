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
    // Squeeze: fit the best layout into a slightly smaller square, letting trees overlap and
    // then clearing the overlap. The shrink adapts per n: doubled on success, halved on failure.
    long squeezeItersPerTree = 2000; // the single-tree fallback; relax does most of the work
    double squeezeStartShrink = 2e-4;
    double squeezeMaxShrink = 1e-2;
    int squeezeFailures = 10; // a visit stops squeezing after this many failures in a row
    long squeezeRelaxIters = 2000; // all-trees relax steps per squeeze (see SqueezeConfig)
    // Shake: kick the current layout (swap, flip or move a few trees, see shake()), settle it
    // back into the best square with a short squeeze, and walk on from there; each settled
    // layout then gets one squeeze into a smaller square.
    int shakeTries = 20;          // per visit
    int shakeKicks = 3;           // up to this many trees changed per try
    long shakeItersPerTree = 200; // single-tree fallback when settling; relax does the rest
    // Neighbours: n+1 minus one of its best few trees to remove, and n-1 plus a tree anywhere
    // in its square, each squeezed into the best square so far.
    int removeCandidates = 3;
    // Re-anneal the best layout cold, with a new seed.
    long refineItersPerTree = 4000;
    double refineTemp = 1e-3;
    double refineStep = 0.01;
    // Small n are also worth fresh random starts (every other pass), annealed hot and long.
    int smallMax = 30;
    long randomItersPerTree = 20000;
    double randomTemp = 0.05;
    double randomStep = 0.2;
};

// best[n] holds the best known n-tree layout (index 0 unused); origins[n] says where it came
// from. Both are updated in place: a layout is only ever replaced by a valid, smaller one.
// Passes go through n = maxN down to minN, each visit running every operator once.
// `save` is called every cfg.saveSeconds and at the end. Stops at the time budget or when
// `stop` turns true (for example on Ctrl-C).
using SaveFn =
    std::function<void(const std::vector<Layout> &best, const std::vector<std::string> &origins)>;
void improve(std::vector<Layout> &best, std::vector<std::string> &origins, const ImproveConfig &cfg,
             std::atomic<bool> &stop, const SaveFn &save);

// Kaggle score of a set of layouts: sum over n of side^2 / n.
double totalScore(const std::vector<Layout> &layouts);
} // namespace packer::packing
