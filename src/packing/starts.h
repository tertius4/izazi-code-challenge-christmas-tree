#pragma once

#include <random>

#include "packing/layout.h"

namespace packer::packing {
// One tree at the rotation (0.01° steps) that gives the smallest square.
Layout bestSingle(const Polygon &shape);

// n trees at random positions and angles, without collisions, in a roomy square.
Layout randomStart(const Polygon &shape, int n, std::mt19937_64 &rng);
} // namespace packer::packing
