#pragma once

#include <vector>

namespace packer::geometry {
struct Point {
    double x = 0;
    double y = 0;
};

// Vertices in order around the outline. The last vertex connects back to the first.
using Polygon = std::vector<Point>;

// Where a tree goes: rotate by `deg` (counter-clockwise, about (0,0)), then shift by (x, y).
struct Pose {
    double x = 0;
    double y = 0;
    double deg = 0;
};

// Axis-aligned bounding box.
struct Box {
    double minX, minY, maxX, maxY;
};

// Only used by the tests, as an independent check of tree::collide.
Polygon transform(const Polygon &shape, const Pose &pose);
// Same as transform, but writes into `out` and reuses its storage.
void transformInto(const Polygon &shape, const Pose &pose, Polygon &out);
// Only used by the tests, as an independent check of tree::collide.
Box boundingBox(const Polygon &poly);
// True if the boxes overlap or touch.
inline bool boxesOverlap(const Box &a, const Box &b) {
    return a.minX <= b.maxX && b.minX <= a.maxX && a.minY <= b.maxY && b.minY <= a.maxY;
}
double polygonArea(const Polygon &poly);
// Maps any angle into [0, 360).
double normalizeDeg(double deg);
} // namespace packer::geometry
