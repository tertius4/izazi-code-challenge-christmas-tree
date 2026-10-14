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

Polygon transform(const Polygon &shape, const Pose &pose);
Box boundingBox(const Polygon &poly);
// True if the boxes overlap or touch.
bool boxesOverlap(const Box &a, const Box &b);
double polygonArea(const Polygon &poly);
// Maps any angle into [0, 360).
double normalizeDeg(double deg);
} // namespace packer::geometry