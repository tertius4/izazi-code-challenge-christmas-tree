#pragma once

#include <atomic>
#include <optional>
#include <random>

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
    // Relax: move every tree at once down the slope of the overlap (L-BFGS), before the
    // single-tree moves and again after them if overlap is left. 0 turns it off.
    long relaxIters = 0;
    double relaxMargin = 1e-8; // how far apart relax aims to leave trees and walls
    // If set, the squeeze gives up early once this turns true.
    const std::atomic<bool> *stop = nullptr;
};

// Tries to fit the trees into a square of side `target`. The trees are first pulled toward the
// centre to fit, overlapping where they must; then single-tree moves shrink the overlap (and
// how far trees stick out of the square) until none is left. The input may overlap. Returns a
// valid layout with side() <= target, or nothing if the overlap could not be cleared.
std::optional<Layout> squeeze(const Layout &layout, double target, const SqueezeConfig &cfg,
                              std::mt19937_64 &rng);

// Changes the arrangement a little, so that a jammed layout can leave its basin: `kicks` times,
// one of: swap where two nearby, differently turned trees stand (each keeps its angle), turn a
// tree by 180°, or move a tree to a random spot in the box at a random angle. The result
// usually overlaps: squeeze it to clear that.
Layout shake(const Layout &layout, int kicks, std::mt19937_64 &rng);
} // namespace packer::packing
