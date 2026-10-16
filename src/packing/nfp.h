#pragma once

#include <random>
#include <vector>

#include "packing/layout.h"

namespace packer::packing {
// The corners of the no-fit polygons of a tree turned by `deg` against the placed tree
// `fixed`. The tree is a union of convex parts, so its no-fit polygon is the union of one
// convex polygon per pair of parts: the places for the moving tree's origin at which those
// two parts overlap (the Minkowski sum of the fixed part and the mirrored moving part). Each
// corner is pushed out from its polygon by `gap`, so the moving tree there just clears those
// two parts. Corners are where a new tree can nestle against an old one.
std::vector<geometry::Point> nfpCorners(const tree::Placed &fixed, double deg, double gap);

// Where to add one more tree to the layout: tries the no-fit polygon corners against every
// tree and a grid over the square, each at `angles` angles (evenly spaced from a random
// start), and returns the pose with the least penalty: overlap² with the other trees plus how
// far it sticks out of the square². Zero penalty means a free spot. With `region`, the origin
// only goes inside it (and the square), and only trees near it are tried against: much
// quicker in a large layout.
Pose bestInsertion(const Layout &layout, const Box &square, int angles, std::mt19937_64 &rng,
                   const Box *region = nullptr);

// Ruin and recreate: removes the `count` trees nearest to one tree (a tree on the edge of the
// layout half the time) and puts them back one by one at the best spots in the hole they
// left, inside `square`. The result usually overlaps a little: squeeze it to clear that.
Layout ruinAndRecreate(const Layout &layout, const Box &square, std::size_t count,
                       std::mt19937_64 &rng);
} // namespace packer::packing
