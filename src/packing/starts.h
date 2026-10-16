#pragma once

#include <optional>
#include <random>

#include "packing/layout.h"

namespace packer::packing {
// One tree at the rotation (0.01° steps) that gives the smallest square.
Layout bestSingle();

// n trees at random positions and angles, without collisions, in a roomy square.
Layout randomStart(int n, std::mt19937_64 &rng);

// n trees at random poses in a square of this side around the origin, overlapping freely.
// A start for squeeze(), which clears the overlap.
Layout denseStart(int n, double side, std::mt19937_64 &rng);

// Like denseStart, but point-symmetric about the origin: tree i + n/2 is tree i turned 180°
// and reflected through the origin (if n is odd, the last tree is free). The order that
// SqueezeConfig::symmetric expects.
Layout symmetricStart(int n, double side, std::mt19937_64 &rng);

// The layout minus the one tree whose removal leaves the smallest square.
Layout withoutTree(const Layout &layout);

// The layout plus one more tree, dropped into a free spot inside its current square, or
// nothing if no spot was found.
std::optional<Layout> withExtraTree(const Layout &layout, std::mt19937_64 &rng);
} // namespace packer::packing
