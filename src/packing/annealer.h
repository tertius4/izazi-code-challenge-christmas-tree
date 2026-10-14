#pragma once

#include <random>

#include "packing/layout.h"

namespace packer::packing {
struct AnnealConfig {
    long iterations = 100000;
    double startTemp = 0.02; // how willing to accept worse moves at the start...
    double endTemp = 1e-5;   // ...and at the end
    double startStep = 0.1;  // largest shift per move at the start (rotation: 100x, in degrees)
    double endStep = 0.001;
    double pullWeight = 0.1; // strength of the pull toward the centre
    // Compaction moves, as shares of all moves. The rest are plain shifts, turns and flips.
    double edgeShare = 0.05;     // push a tree on the edge of the square toward the centre
    double relocateShare = 0.05; // drop a tree on the edge somewhere inside the square
    double squeezeShare = 0.01;  // pull every tree toward the centre at once
};

// Improves the layout in place. Never accepts a move that makes trees collide, and
// returns the best layout seen, so the square never gets bigger.
void anneal(Layout &layout, const AnnealConfig &cfg, std::mt19937_64 &rng);
} // namespace packer::packing
