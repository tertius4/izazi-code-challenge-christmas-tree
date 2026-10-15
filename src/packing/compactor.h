#pragma once

#include <random>

#include "packing/layout.h"

namespace packer::packing {
// Turns the whole layout about its centre to the angle that gives the smallest axis-aligned
// square. Trees keep their relative positions, so the layout stays valid; if rounding makes
// the turned layout collide, it is left as it was. Returns true if the square got smaller.
bool rotateWhole(Layout &layout);

// Shrink and repair: pulls every tree toward the centre so the square gets a little smaller,
// then tries to move the overlapping trees apart without leaving the new square. Keeps the
// shrink if every overlap is resolved, otherwise undoes it and tries a smaller shrink.
// `repairIters` is the repair budget per attempt, per tree. The layout is always valid after.
void compact(Layout &layout, std::mt19937_64 &rng, long repairIters);
} // namespace packer::packing
