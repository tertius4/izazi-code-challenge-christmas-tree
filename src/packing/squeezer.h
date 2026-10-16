#pragma once

#include <atomic>
#include <optional>
#include <random>
#include <vector>

#include "packing/layout.h"

namespace packer::packing {
struct SqueezeConfig {
    long iterations = 1000000;
    // Temperatures are relative to the starting penalty per tree, so they suit any n and any
    // amount of squeezing.
    double startTemp = 0.3;
    double endTemp = 1e-4;
    double startStep = 0.005; // largest shift per move at the start (rotation: 100x, in degrees)
    double endStep = 1e-4;
    double badShare = 0.7; // share of moves spent on trees that still overlap or stick out
    // Share of single-tree moves that try several distances along one direction (x, y or
    // turning) and take the best of them, instead of one random step.
    double lineShare = 0.1;
    // Relax: move every tree at once down the slope of the overlap (L-BFGS), before the
    // single-tree moves and again after them if overlap is left. 0 turns it off.
    long relaxIters = 0;
    double relaxMargin = 1e-8; // how far apart relax aims to leave trees and walls
    // Guided local search: when relax gets stuck with overlap left, raise the penalty of the
    // pairs that still overlap (most for the deepest) and carry on, up to this many times.
    int glsRounds = 0;
    // If not empty, only trees marked true move; the others stay exactly where they are (and
    // must already fit the square).
    std::vector<bool> movable;
    // Relax keeps the layout point-symmetric about the square's centre: tree i + n/2 stays the
    // mirror image of tree i (turned 180° and reflected through the centre). If n is odd, the
    // last tree moves freely. The input must already be built that way (symmetricStart()).
    bool symmetric = false;
    // If set, the squeeze gives up early once this turns true.
    const std::atomic<bool> *stop = nullptr;
};

// Tries to fit the trees into a square of side `target`. The trees are first pulled toward the
// centre to fit, overlapping where they must; then single-tree moves shrink the overlap (and
// how far trees stick out of the square) until none is left. The input may overlap. Returns a
// valid layout with side() <= target, or nothing if the overlap could not be cleared.
std::optional<Layout> squeeze(const Layout &layout, double target, const SqueezeConfig &cfg,
                              std::mt19937_64 &rng);

struct CompressConfig {
    long iterations = 3000; // L-BFGS steps per round
    int rounds = 5;
    // Weight of the overlap against the side: it rises from the first to the last value over
    // the rounds, so the trees first crowd in, then push each other apart.
    double startWeight = 1e3;
    double endWeight = 1e8;
    double margin = 1e-8;
    std::vector<bool> movable; // as in SqueezeConfig
    bool symmetric = false;    // as in SqueezeConfig
    const std::atomic<bool> *stop = nullptr;
};

// Treats the square's side as one more variable: minimises side + weight·overlap² over the
// trees' poses and the side together, so the walls close in as far as the trees let them.
// The layout is left at the end point, which may still overlap a little; returns the side it
// reached. Squeeze the layout into that side (or a hair more) to make it valid.
double compress(Layout &layout, const CompressConfig &cfg);

// Changes the arrangement a little, so that a jammed layout can leave its basin: `kicks` times,
// one of: swap where two nearby, differently turned trees stand (each keeps its angle), turn a
// tree by 180°, or move a tree to a random spot in the box at a random angle. The result
// usually overlaps: squeeze it to clear that.
Layout shake(const Layout &layout, int kicks, std::mt19937_64 &rng);

// Trees whose box comes within `width` of the edge of the square around the layout: marked
// true, the rest false. For SqueezeConfig::movable.
std::vector<bool> borderTrees(const Layout &layout, double width);
} // namespace packer::packing
