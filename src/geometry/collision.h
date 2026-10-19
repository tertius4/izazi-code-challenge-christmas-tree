#pragma once

#include "geometry/polygon.h"

namespace packer::geometry {
// Anything closer than about this counts as touching. Touching is treated as a collision,
// so every accepted layout has a real gap and survives rounding in the output.
inline constexpr double kEps = 1e-13;

// Only used by the tests, as an independent check of tree::collide.
// True if the point is strictly inside poly (even-odd ray casting).
bool pointInPolygon(const Point &point, const Polygon &poly);

// Only used by the tests, as an independent check of tree::collide.
// True if the polygons overlap or touch. Pass each polygon's bounding box so the
// expensive test is skipped when the boxes are apart.
bool polygonsCollide(const Polygon &a, const Box &boxA, const Polygon &b, const Box &boxB);
} // namespace packer::geometry
