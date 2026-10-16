#pragma once

#include <random>

#include "packing/layout.h"

namespace packer::packing {
// A child of two layouts: b is turned by a random multiple of 90° (and maybe mirrored) and
// centred on a; then a random line near the centre splits the plane, and the child takes a's
// trees on one side and b's on the other. Trees are then dropped (the most overlapped first)
// or added (at the best spots, see bestInsertion()) until there are exactly n, all within a
// square of a's side. The child usually overlaps where the halves meet: squeeze it.
Layout crossover(const Layout &a, const Layout &b, int n, std::mt19937_64 &rng);
} // namespace packer::packing
