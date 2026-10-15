#pragma once

#include <array>
#include <cstddef>

#include "geometry/polygon.h"

// The one shape this program packs: the competition's Christmas tree.
namespace packer::tree {
using geometry::Box;
using geometry::Point;
using geometry::Polygon;
using geometry::Pose;

// The outline, clockwise from the tip: 15 vertices, 1.0 tall, 0.7 wide, with the trunk's top
// at the origin.
const Polygon &outline();

// The tree is exactly the union of 4 convex parts: the three tiers, each extended up to a
// triangle, and the trunk. Collisions between trees reduce to collisions between parts.
inline constexpr std::size_t kParts = 4;

// Every vertex lies within kRadius of the point (0, kCentreY) of the unturned tree. So two
// trees can only collide if those points (moved with the trees) are within 2 * kRadius.
// (The tip and the trunk's bottom corners are exactly 0.5028125 away; the rest are closer.)
inline constexpr double kCentreY = 0.2971875;
inline constexpr double kRadius = 0.5028126;

// A convex part placed in the plane: its vertices counter-clockwise, the unit outward normal
// of each edge (edge k runs from v[k] to v[k+1]), and its box.
struct Part {
    std::array<Point, 4> v{};
    std::array<Point, 4> normal{};
    std::size_t count = 0;
    Box box{};
};

// A tree placed at a pose, ready for collision checks.
struct Placed {
    std::array<Part, kParts> parts{};
    Box box{};          // around the whole tree
    Point centre{};     // the turned and moved (0, kCentreY)
};

void place(const Pose &pose, Placed &out);
Placed place(const Pose &pose);

// True if the trees overlap or touch. Touching counts, so every accepted layout keeps a real
// gap that survives rounding in the output (the same rule as geometry::polygonsCollide).
bool collide(const Placed &a, const Placed &b);
} // namespace packer::tree
