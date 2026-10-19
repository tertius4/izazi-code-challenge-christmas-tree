#include "geometry/collision.h"

#include <algorithm>
#include <cstddef>

namespace packer::geometry {
namespace {
// > 0 if c is left of the line a->b, < 0 if right, about 0 if on the line.
double cross(const Point &a, const Point &b, const Point &c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

int sign(double value) {
    return value > kEps ? 1 : (value < -kEps ? -1 : 0);
}

// For p known to lie on the line a-b: is it between a and b?
bool withinSegment(const Point &a, const Point &b, const Point &p) {
    return std::min(a.x, b.x) - kEps <= p.x && p.x <= std::max(a.x, b.x) + kEps &&
           std::min(a.y, b.y) - kEps <= p.y && p.y <= std::max(a.y, b.y) + kEps;
}

// True if segments p1-p2 and q1-q2 share any point, including just touching.
bool segmentsIntersect(const Point &p1, const Point &p2, const Point &q1, const Point &q2) {
    // Which side of the other segment's line each endpoint is on.
    const int p1Side = sign(cross(q1, q2, p1));
    const int p2Side = sign(cross(q1, q2, p2));
    const int q1Side = sign(cross(p1, p2, q1));
    const int q2Side = sign(cross(p1, p2, q2));

    const bool crossProperly = p1Side * p2Side < 0 && q1Side * q2Side < 0;
    if (crossProperly) return true;

    // Otherwise they can only meet if an endpoint lies on the other segment.
    const bool p1OnQ = p1Side == 0 && withinSegment(q1, q2, p1);
    const bool p2OnQ = p2Side == 0 && withinSegment(q1, q2, p2);
    const bool q1OnP = q1Side == 0 && withinSegment(p1, p2, q1);
    const bool q2OnP = q2Side == 0 && withinSegment(p1, p2, q2);

    return p1OnQ || p2OnQ || q1OnP || q2OnP;
}
} // namespace

bool pointInPolygon(const Point &point, const Polygon &poly) {
    // Shoot a ray to the right; an odd number of edge crossings means inside.
    bool inside = false;
    for (std::size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++) {
        const Point &a = poly[i];
        const Point &b = poly[j];

        const bool edgeSpansRay = (a.y > point.y) != (b.y > point.y);
        if (edgeSpansRay) {
            const double xCross = a.x + (point.y - a.y) * (b.x - a.x) / (b.y - a.y);
            if (point.x < xCross) inside = !inside;
        }
    }

    return inside;
}

bool polygonsCollide(const Polygon &a, const Box &boxA, const Polygon &b, const Box &boxB) {
    if (!boxesOverlap(boxA, boxB)) return false;

    // An edge lying wholly outside the other polygon's box can't meet any of its edges.
    // The margin keeps this a pure speedup: it only skips edges that are clearly apart.
    constexpr double kMargin = 1e-9;
    auto edgeOutside = [&](const Point &start, const Point &end, const Box &box) {
        const bool leftOf = std::max(start.x, end.x) < box.minX - kMargin;
        const bool rightOf = std::min(start.x, end.x) > box.maxX + kMargin;
        const bool below = std::max(start.y, end.y) < box.minY - kMargin;
        const bool above = std::min(start.y, end.y) > box.maxY + kMargin;

        return leftOf || rightOf || below || above;
    };

    for (std::size_t i = 0; i < a.size(); ++i) {
        const Point &a1 = a[i];
        const Point &a2 = a[i + 1 == a.size() ? 0 : i + 1];

        if (edgeOutside(a1, a2, boxB)) continue;
        for (std::size_t j = 0; j < b.size(); ++j) {
            const Point &b1 = b[j];
            const Point &b2 = b[j + 1 == b.size() ? 0 : j + 1];

            if (edgeOutside(b1, b2, boxA)) continue;
            if (segmentsIntersect(a1, a2, b1, b2)) return true;
        }
    }

    // No edges meet: the polygons are apart, or one is completely inside the other.
    return pointInPolygon(a[0], b) || pointInPolygon(b[0], a);
}
} // namespace packer::geometry
