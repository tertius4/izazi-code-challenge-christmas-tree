#include "geometry/polygon.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>

namespace packer::geometry {
// Rotate (around 0,0) and shift (on the x and y axis) the polygon
void transformInto(const Polygon &shape, const Pose &pose, Polygon &out) {
    const double rad = pose.deg * std::numbers::pi / 180.0;
    const double c = std::cos(rad);
    const double s = std::sin(rad);
    out.resize(shape.size());
    for (std::size_t i = 0; i < shape.size(); ++i) {
        const Point &p = shape[i];
        // Standard 2D rotation matrix, then translation.
        out[i] = {p.x * c - p.y * s + pose.x, p.x * s + p.y * c + pose.y};
    }
}

Polygon transform(const Polygon &shape, const Pose &pose) {
    Polygon out;
    transformInto(shape, pose, out);
    return out;
}

// Get the min and max of x and y axis of the polygon.
Box boundingBox(const Polygon &poly) {
    Box b{poly[0].x, poly[0].y, poly[0].x, poly[0].y};
    for (const Point &p : poly) {
        b.minX = std::min(b.minX, p.x);
        b.minY = std::min(b.minY, p.y);
        b.maxX = std::max(b.maxX, p.x);
        b.maxY = std::max(b.maxY, p.y);
    }
    return b;
}

double polygonArea(const Polygon &poly) {
    // Shoelace formula.
    double sum = 0;
    for (std::size_t i = 0; i < poly.size(); ++i) {
        const Point &a = poly[i];
        const Point &b = poly[(i + 1) % poly.size()];
        sum += a.x * b.y - b.x * a.y;
    }
    return std::abs(sum) / 2.0;
}

double normalizeDeg(double deg) {
    deg = std::fmod(deg, 360.0);
    return deg < 0 ? deg + 360.0 : deg;
}
} // namespace packer::geometry