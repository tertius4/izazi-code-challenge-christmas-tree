#include "geometry/polygon.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>

namespace packer::geometry {
// Rotate (around 0,0) and shift (on the x and y axis) the polygon
void transformInto(const Polygon &shape, const Pose &pose, Polygon &out) {
    const double rad = pose.deg * std::numbers::pi / 180.0;
    const double cosine = std::cos(rad);
    const double sine = std::sin(rad);

    out.resize(shape.size());
    for (std::size_t i = 0; i < shape.size(); ++i) {
        const Point &point = shape[i];
        // Standard 2D rotation matrix, then translation.
        out[i] = {point.x * cosine - point.y * sine + pose.x,
                  point.x * sine + point.y * cosine + pose.y};
    }
}

Polygon transform(const Polygon &shape, const Pose &pose) {
    Polygon out;
    transformInto(shape, pose, out);

    return out;
}

// Get the min and max of x and y axis of the polygon.
Box boundingBox(const Polygon &poly) {
    Box box{poly[0].x, poly[0].y, poly[0].x, poly[0].y};
    for (const Point &point : poly) {
        box.minX = std::min(box.minX, point.x);
        box.minY = std::min(box.minY, point.y);
        box.maxX = std::max(box.maxX, point.x);
        box.maxY = std::max(box.maxY, point.y);
    }

    return box;
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
