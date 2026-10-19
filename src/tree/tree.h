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
// of each edge (edge k runs from vertex[k] to vertex[k+1]), and its box.
struct Part {
    std::array<Point, 4> vertex{};
    std::array<Point, 4> normal{};
    std::size_t count = 0;
    Box box{};
};

// A tree placed at a pose, ready for collision checks.
struct Placed {
    std::array<Part, kParts> parts{};
    Box box{};      // around the whole tree
    Point centre{}; // the turned and moved (0, kCentreY)
};

void place(const Pose &pose, Placed &out);
Placed place(const Pose &pose);

// True if the trees overlap or touch. Touching counts, so every accepted layout keeps a real
// gap that survives rounding in the output (the same rule as geometry::polygonsCollide).
bool collide(const Placed &a, const Placed &b);

// How deep the trees overlap: for each pair of colliding parts, the shortest distance one must
// move to clear the other, summed over the pairs. Positive exactly when collide() is true, so
// a layout with no penetration anywhere is valid.
double penetration(const Placed &a, const Placed &b);

// penetration() with a margin, and its slope. `depth` sums, over the part pairs closer than
// `margin` to touching, how deep they overlap plus the margin; the other fields say how fast
// that sum grows as tree a moves in x, in y, or turns (per degree) about `origin`, the x and y
// of its pose. Tree b stays put. Zero exactly when every part pair is more than `margin` apart.
struct Overlap {
    double depth = 0;
    double dx = 0;
    double dy = 0;
    double ddeg = 0;
};
Overlap overlap(const Placed &a, Point origin, const Placed &b, double margin);
} // namespace packer::tree
