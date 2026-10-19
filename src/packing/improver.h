#pragma once

#include <atomic>
#include <functional>
#include <string>
#include <vector>

#include "packing/crystal.h"
#include "packing/layout.h"

namespace packer::packing {
// Settings for improving a set of best-known layouts over a time budget.
struct ImproveConfig {
    double minutes = 10; // time budget
    int threads = 0;     // 0: half the cores, so the computer stays usable
    int minN = 1;        // only work on n in [minN, maxN]
    int maxN = 200;
    double saveSeconds = 120; // save this often, so an interrupted run loses little
    std::string statsPath;    // if set, the operator stats are also written here on every save
    unsigned seed = 0;        // 0: from the clock, so every run tries something new
    // Squeeze: fit the best layout into a slightly smaller square, letting trees overlap and
    // then clearing the overlap. The shrink adapts per n: doubled on success, halved on failure.
    long squeezeItersPerTree = 2000; // the single-tree fallback; relax does most of the work
    double squeezeStartShrink = 2e-4;
    double squeezeMaxShrink = 1e-2;
    int squeezeFailures = 10; // a visit stops squeezing after this many failures in a row
    long squeezeRelaxIters = 2000; // all-trees relax steps per squeeze (see SqueezeConfig)
    int squeezeGlsRounds = 2;      // penalty raises when relax is stuck (see SqueezeConfig)
    double squeezeLineShare = 0.1; // share of best-of-several single-tree moves
    // Border: squeeze with only the trees within this distance of the walls free to move.
    double borderWidth = 1.0;
    // Compress: make the side a variable (see compress()), then squeeze into where it got.
    // Repeats from each win until it fails compressFailures times in a row, at most
    // compressTries times a visit (it wins tiny amounts nearly every time).
    long compressIters = 1500;
    int compressRounds = 4;
    int compressFailures = 2;
    int compressTries = 4;
    // Shake: kick the current layout (swap, flip or move a few trees, see shake()), settle it
    // back into the best square with a short squeeze, and walk on from there; each settled
    // layout then gets one squeeze into a smaller square.
    int shakeTries = 8;           // per visit
    int shakeKicks = 3;           // up to this many trees changed per try
    long shakeItersPerTree = 200; // single-tree fallback when settling; relax does the rest
    // Neighbours: n+1 minus one of its best few trees to remove, and n-1 plus a tree (at the
    // best spot or a random one), each squeezed into the best square so far.
    int removeCandidates = 1;
    // A win on n queues "n plus a tree" for n+1 right away, if it gained more than this.
    double urgentMinGain = 1e-5;
    // Elites: per n, the few best distinct valid layouts seen, within this fraction of the
    // best side. Crossover mixes the best with an elite (or a turned copy of itself, or a
    // neighbour) and settles the child.
    int eliteSize = 6;
    double eliteSlack = 0.003;
    int crossoverTries = 2; // per visit
    // Crystal cuts, for n >= crystalMinN: cut the densest patterns into squares at many
    // angles, compress, and squeeze down from there.
    std::vector<Crystal> crystals;
    int crystalMinN = 20;
    int crystalAngles = 12;
    double crystalSlack = 0.005; // squeeze down only from cuts within this of the best side
    int crystalDescent = 8;     // squeezes when going down from a cut
    // Dense starts, for n <= denseMax: random overlapping trees (half the time point-symmetric)
    // in a roomy square, compressed and squeezed valid, then squeezed down if within
    // denseSlack of the best side. Small n get more of them.
    int denseMax = 40;
    int denseStarts = 1;     // per visit...
    int denseSmallMax = 10;  // ...times denseSmallBoost for n up to this
    int denseSmallBoost = 6;
    double denseSlack = 0.05;
    int denseDescent = 12;
    int denseGlsRounds = 4;
    // Ruin and recreate (see ruinAndRecreate()): take out a patch of trees and put them back
    // at the best spots, settle into the best square and walk on from there, like shaking but
    // a whole patch at once. For n >= ruinMinN; the patch is 3 trees up to ruinShare of n.
    int ruinTries = 2; // per visit
    int ruinMinN = 8;
    double ruinShare = 0.12;
    // Refine: a walk of re-anneals (no overlap allowed, but trees can relocate into holes),
    // each at the next temperature of refineTemps and then squeezed into a smaller square.
    int refineTries = 4; // per visit
    long refineItersPerTree = 4000;
    std::vector<double> refineTemps = {1e-3, 3e-3, 1e-2};
    double refineStep = 0.01;
    // Small n are also worth fresh random starts (every other pass), annealed hot and long.
    int smallMax = 30;
    long randomItersPerTree = 20000;
    double randomTemp = 0.05;
    double randomStep = 0.2;
    // Which n to visit next: this share of visits goes to the n with the best recent score
    // gain per second (each n starts as if it had run banditPriorSeconds at the average rate),
    // the rest to the plain sweep from maxN down to minN.
    double banditShare = 0.5;
    double banditPriorSeconds = 60;
    // Which operators run: each try of an operator (other than the squeeze, which always runs)
    // goes ahead with probability its score gain per second this run over the average of all
    // of them, at least opMinChance. Each starts as if it had run opPriorSeconds at the average.
    // Rates are kept per band of n (see kBands in improver.cpp), so an operator that pays off
    // only for small n is not judged by its record for large n.
    bool opBandit = true;
    double opMinChance = 0.03;
    double opPriorSeconds = 120;
};

// best[n] holds the best known n-tree layout (index 0 unused); origins[n] says where it came
// from. Both are updated in place: a layout is only ever replaced by a valid, smaller one.
// Visits run every operator on one n; see ImproveConfig for the operators and how n are
// picked. `save` is called every cfg.saveSeconds and at the end. Stops at the time budget or
// when `stop` turns true (for example on Ctrl-C).
using SaveFn =
    std::function<void(const std::vector<Layout> &best, const std::vector<std::string> &origins)>;
void improve(std::vector<Layout> &best, std::vector<std::string> &origins, const ImproveConfig &cfg,
             std::atomic<bool> &stop, const SaveFn &save);

// Kaggle score of a set of layouts: sum over n of side^2 / n.
double totalScore(const std::vector<Layout> &layouts);
} // namespace packer::packing
