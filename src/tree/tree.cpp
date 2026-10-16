#include "tree/tree.h"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace packer::tree {
namespace {
// The tiers' slanted edges, extended up to the middle: the middle tier's edge runs from
// (0.0625, 0.5) to (0.2, 0.25), so it meets x = 0 at this height.
constexpr double kMiddleApex = 0.5 + 0.0625 * 0.25 / 0.1375;

// The parts of the unturned tree, counter-clockwise, with their outward unit normals.
struct BasePart {
    std::array<Point, 4> v{};
    std::array<Point, 4> normal{};
    std::size_t count = 0;
};

BasePart makePart(std::initializer_list<Point> pts) {
    BasePart p;
    for (const Point &q : pts) {
        p.v[p.count++] = q;
    }

    for (std::size_t k = 0; k < p.count; ++k) {
        const Point &a = p.v[k];
        const Point &b = p.v[(k + 1) % p.count];
        const double len = std::hypot(b.x - a.x, b.y - a.y);
        p.normal[k] = {(b.y - a.y) / len, -(b.x - a.x) / len}; // right of a->b: outward for CCW
    }

    return p;
}

const std::array<BasePart, kParts> &baseParts() {
    static const std::array<BasePart, kParts> parts = {
        makePart({{0, 0.8}, {-0.125, 0.5}, {0.125, 0.5}}),                  // top tier
        makePart({{0, kMiddleApex}, {-0.2, 0.25}, {0.2, 0.25}}),            // middle tier
        makePart({{0, 0.35}, {-0.35, 0}, {0.35, 0}}),                       // bottom tier
        makePart({{-0.075, -0.2}, {0.075, -0.2}, {0.075, 0}, {-0.075, 0}}), // trunk
    };

    return parts;
}

// Separating-axis test on one part's edge normals: true if one of them separates a from b
// by more than a hair. Touching (a gap of about 0) does not separate.
bool separatedAlong(const Part &a, const Part &b) {
    constexpr double kGap = 1e-12;
    for (std::size_t k = 0; k < a.count; ++k) {
        const Point &n = a.normal[k];
        // a lies on the inner side of its own edge k: its maximum along n is that edge.
        const double aMax = n.x * a.v[k].x + n.y * a.v[k].y;
        double bMin = INFINITY;
        for (std::size_t j = 0; j < b.count; ++j) {
            bMin = std::min(bMin, n.x * b.v[j].x + n.y * b.v[j].y);
        }

        if (bMin > aMax + kGap) return true;
    }

    return false;
}

// How far b reaches past a's edges: the smallest push along one of a's edge normals that
// would clear b. Below -kGap means that normal separates them.
double reachPast(const Part &a, const Part &b) {
    double least = INFINITY;
    for (std::size_t k = 0; k < a.count; ++k) {
        const Point &n = a.normal[k];
        const double aMax = n.x * a.v[k].x + n.y * a.v[k].y;
        double bMin = INFINITY;
        for (std::size_t j = 0; j < b.count; ++j) {
            bMin = std::min(bMin, n.x * b.v[j].x + n.y * b.v[j].y);
        }

        least = std::min(least, aMax - bMin);
    }

    return least;
}

bool partsCollide(const Part &a, const Part &b) {
    if (!geometry::boxesOverlap(a.box, b.box)) return false;

    return !separatedAlong(a, b) && !separatedAlong(b, a);
}
} // namespace

const Polygon &outline() {
    static const Polygon tree = {
        {0, 0.8},   {0.125, 0.5}, {0.0625, 0.5}, {0.2, 0.25},    {0.1, 0.25},
        {0.35, 0},  {0.075, 0},   {0.075, -0.2}, {-0.075, -0.2}, {-0.075, 0},
        {-0.35, 0}, {-0.1, 0.25}, {-0.2, 0.25},  {-0.0625, 0.5}, {-0.125, 0.5},
    };

    return tree;
}

void place(const Pose &pose, Placed &out) {
    const double rad = pose.deg * std::numbers::pi / 180.0;
    const double c = std::cos(rad);
    const double s = std::sin(rad);

    out.box = {INFINITY, INFINITY, -INFINITY, -INFINITY};
    const auto &base = baseParts();
    for (std::size_t p = 0; p < kParts; ++p) {
        const BasePart &from = base[p];
        Part &to = out.parts[p];
        to.count = from.count;
        to.box = {INFINITY, INFINITY, -INFINITY, -INFINITY};

        for (std::size_t k = 0; k < from.count; ++k) {
            const Point &v = from.v[k];
            const Point &n = from.normal[k];
            to.v[k] = {v.x * c - v.y * s + pose.x, v.x * s + v.y * c + pose.y};
            to.normal[k] = {n.x * c - n.y * s, n.x * s + n.y * c};
            to.box.minX = std::min(to.box.minX, to.v[k].x);
            to.box.minY = std::min(to.box.minY, to.v[k].y);
            to.box.maxX = std::max(to.box.maxX, to.v[k].x);
            to.box.maxY = std::max(to.box.maxY, to.v[k].y);
        }

        out.box.minX = std::min(out.box.minX, to.box.minX);
        out.box.minY = std::min(out.box.minY, to.box.minY);
        out.box.maxX = std::max(out.box.maxX, to.box.maxX);
        out.box.maxY = std::max(out.box.maxY, to.box.maxY);
    }

    out.centre = {-kCentreY * s + pose.x, kCentreY * c + pose.y};
}

Placed place(const Pose &pose) {
    Placed out;
    place(pose, out);

    return out;
}

bool collide(const Placed &a, const Placed &b) {
    const double dx = a.centre.x - b.centre.x;
    const double dy = a.centre.y - b.centre.y;
    const bool tooFar = dx * dx + dy * dy > 4 * kRadius * kRadius;
    if (tooFar) return false;
    if (!geometry::boxesOverlap(a.box, b.box)) return false;

    for (const Part &pa : a.parts) {
        if (!geometry::boxesOverlap(pa.box, b.box)) continue;
        for (const Part &pb : b.parts) {
            if (partsCollide(pa, pb)) return true;
        }
    }

    return false;
}

double penetration(const Placed &a, const Placed &b) {
    constexpr double kGap = 1e-12; // as in separatedAlong: touching counts as overlap
    const double dx = a.centre.x - b.centre.x;
    const double dy = a.centre.y - b.centre.y;
    const bool tooFar = dx * dx + dy * dy > 4 * kRadius * kRadius;
    if (tooFar) return 0;
    if (!geometry::boxesOverlap(a.box, b.box)) return 0;

    double total = 0;
    for (const Part &pa : a.parts) {
        if (!geometry::boxesOverlap(pa.box, b.box)) continue;
        for (const Part &pb : b.parts) {
            if (!geometry::boxesOverlap(pa.box, pb.box)) continue;

            const double depth = std::min(reachPast(pa, pb), reachPast(pb, pa));
            // Exactly the parts that collide: no normal separates them by more than kGap.
            if (depth >= -kGap) total += depth + 2 * kGap;
        }
    }

    return total;
}

Overlap overlap(const Placed &a, Point origin, const Placed &b, double margin) {
    Overlap out;
    const double dx = a.centre.x - b.centre.x;
    const double dy = a.centre.y - b.centre.y;
    const double reach = 2 * kRadius + margin;
    if (dx * dx + dy * dy > reach * reach) return out;

    auto near = [margin](const Box &p, const Box &q) {
        return p.minX <= q.maxX + margin && q.minX <= p.maxX + margin &&
               p.minY <= q.maxY + margin && q.minY <= p.maxY + margin;
    };
    if (!near(a.box, b.box)) return out;

    constexpr double kRadPerDeg = std::numbers::pi / 180.0;
    auto perp = [](const Point &v) { return Point{-v.y, v.x}; }; // turned 90° counter-clockwise
    auto dot = [](const Point &u, const Point &v) { return u.x * v.x + u.y * v.y; };
    auto from = [&](const Point &v) { return Point{v.x - origin.x, v.y - origin.y}; };

    for (const Part &pa : a.parts) {
        if (!near(pa.box, b.box)) continue;
        for (const Part &pb : b.parts) {
            if (!near(pa.box, pb.box)) continue;

            // The shallowest axis decides the depth, as in reachPast. Along one of a's normals
            // n, a's edge and n turn with a, b's deepest vertex w stays: the depth moves with
            // n and turns by perp(n)·(origin - w). Along one of b's normals m, a's deepest
            // vertex u moves away from b's edge: by -m, and turns by -m·perp(u - origin).
            double depth = INFINITY;
            Point slope{};
            double turn = 0;
            for (std::size_t k = 0; k < pa.count; ++k) {
                const Point &n = pa.normal[k];
                std::size_t w = 0;
                for (std::size_t j = 1; j < pb.count; ++j) {
                    if (dot(n, pb.v[j]) < dot(n, pb.v[w])) w = j;
                }
                const double d = dot(n, pa.v[k]) - dot(n, pb.v[w]);
                if (d < depth) {
                    depth = d;
                    slope = n;
                    turn = dot(perp(n), Point{origin.x - pb.v[w].x, origin.y - pb.v[w].y});
                }
            }
            for (std::size_t k = 0; k < pb.count; ++k) {
                const Point &m = pb.normal[k];
                std::size_t u = 0;
                for (std::size_t j = 1; j < pa.count; ++j) {
                    if (dot(m, pa.v[j]) < dot(m, pa.v[u])) u = j;
                }
                const double d = dot(m, pb.v[k]) - dot(m, pa.v[u]);
                if (d < depth) {
                    depth = d;
                    slope = {-m.x, -m.y};
                    turn = -dot(m, perp(from(pa.v[u])));
                }
            }

            if (depth <= -margin) continue;
            out.depth += depth + margin;
            out.dx += slope.x;
            out.dy += slope.y;
            out.ddeg += turn * kRadPerDeg;
        }
    }

    return out;
}
} // namespace packer::tree
