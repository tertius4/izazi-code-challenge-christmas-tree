#pragma once

#include <optional>
#include <random>

#include "packing/layout.h"

namespace packer::packing {
// One tree at the rotation (0.01° steps) that gives the smallest square.
Layout bestSingle();

// n trees at random positions and angles, without collisions, in a roomy square.
Layout randomStart(int n, std::mt19937_64 &rng);

// The layout minus the one tree whose removal leaves the smallest square.
Layout withoutTree(const Layout &layout);

// The layout plus one more tree, dropped into a free spot inside its current square, or
// nothing if no spot was found.
std::optional<Layout> withExtraTree(const Layout &layout, std::mt19937_64 &rng);
} // namespace packer::packing
