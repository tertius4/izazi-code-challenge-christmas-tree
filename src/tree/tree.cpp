#include "tree/tree.h"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace packer::tree {
namespace {
using geometry::kInf;

// The tiers' slanted edges, extended up to the middle: the middle tier's edge runs from
// (0.0625, 0.5) to (0.2, 0.25), so it meets x = 0 at this height.
constexpr double kMiddleApex = 0.5 + 0.0625 * 0.25 / 0.1375;

// The parts of the unturned tree, counter-clockwise, with their outward unit normals.
struct BasePart {
    std::array<Point, 4> vertex{};
    std::array<Point, 4> normal{};
    std::size_t count = 0;
};

BasePart makePart(std::initializer_list<Point> vertices) {
    BasePart part;
    for (const Point &vertex : vertices) {
        part.vertex[part.count++] = vertex;
    }

    for (std::size_t k = 0; k < part.count; ++k) {
        const Point &start = part.vertex[k];
        const Point &end = part.vertex[(k + 1) % part.count];
        const double len = std::hypot(end.x - start.x, end.y - start.y);
        // Right of start->end: outward for counter-clockwise vertices.
        part.normal[k] = {(end.y - start.y) / len, -(end.x - start.x) / len};
    }

    return part;
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
        const Point &normal = a.normal[k];
        // a lies on the inner side of its own edge k: its maximum along the normal is that edge.
        const double aMax = normal.x * a.vertex[k].x + normal.y * a.vertex[k].y;
        double bMin = kInf;
        for (std::size_t j = 0; j < b.count; ++j) {
            bMin = std::min(bMin, normal.x * b.vertex[j].x + normal.y * b.vertex[j].y);
        }

        if (bMin > aMax + kGap) return true;
    }

    return false;
}

// How far b reaches past a's edges: the smallest push along one of a's edge normals that
// would clear b. Below -kGap means that normal separates them.
double reachPast(const Part &a, const Part &b) {
    double least = kInf;
    for (std::size_t k = 0; k < a.count; ++k) {
        const Point &normal = a.normal[k];
        const double aMax = normal.x * a.vertex[k].x + normal.y * a.vertex[k].y;
        double bMin = kInf;
        for (std::size_t j = 0; j < b.count; ++j) {
            bMin = std::min(bMin, normal.x * b.vertex[j].x + normal.y * b.vertex[j].y);
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
    const double cosine = std::cos(rad);
    const double sine = std::sin(rad);

    out.box = {kInf, kInf, -kInf, -kInf};
    const auto &base = baseParts();
    for (std::size_t index = 0; index < kParts; ++index) {
        const BasePart &from = base[index];
        Part &to = out.parts[index];
        to.count = from.count;
        to.box = {kInf, kInf, -kInf, -kInf};

        for (std::size_t k = 0; k < from.count; ++k) {
            const Point &vertex = from.vertex[k];
            const Point &normal = from.normal[k];
            to.vertex[k] = {vertex.x * cosine - vertex.y * sine + pose.x,
                            vertex.x * sine + vertex.y * cosine + pose.y};
            to.normal[k] = {normal.x * cosine - normal.y * sine,
                            normal.x * sine + normal.y * cosine};
            to.box.minX = std::min(to.box.minX, to.vertex[k].x);
            to.box.minY = std::min(to.box.minY, to.vertex[k].y);
            to.box.maxX = std::max(to.box.maxX, to.vertex[k].x);
            to.box.maxY = std::max(to.box.maxY, to.vertex[k].y);
        }

        out.box.minX = std::min(out.box.minX, to.box.minX);
        out.box.minY = std::min(out.box.minY, to.box.minY);
        out.box.maxX = std::max(out.box.maxX, to.box.maxX);
        out.box.maxY = std::max(out.box.maxY, to.box.maxY);
    }

    out.centre = {-kCentreY * sine + pose.x, kCentreY * cosine + pose.y};
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

    for (const Part &partA : a.parts) {
        if (!geometry::boxesOverlap(partA.box, b.box)) continue;
        for (const Part &partB : b.parts) {
            if (partsCollide(partA, partB)) return true;
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
    for (const Part &partA : a.parts) {
        if (!geometry::boxesOverlap(partA.box, b.box)) continue;
        for (const Part &partB : b.parts) {
            if (!geometry::boxesOverlap(partA.box, partB.box)) continue;

            const double depth = std::min(reachPast(partA, partB), reachPast(partB, partA));
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

    auto near = [margin](const Box &boxA, const Box &boxB) {
        return boxA.minX <= boxB.maxX + margin && boxB.minX <= boxA.maxX + margin &&
               boxA.minY <= boxB.maxY + margin && boxB.minY <= boxA.maxY + margin;
    };
    if (!near(a.box, b.box)) return out;

    constexpr double kRadPerDeg = std::numbers::pi / 180.0;
    auto perp = [](const Point &v) { return Point{-v.y, v.x}; }; // turned 90° counter-clockwise
    auto dot = [](const Point &u, const Point &v) { return u.x * v.x + u.y * v.y; };
    auto fromOrigin = [&](const Point &v) { return Point{v.x - origin.x, v.y - origin.y}; };

    for (const Part &partA : a.parts) {
        if (!near(partA.box, b.box)) continue;
        for (const Part &partB : b.parts) {
            if (!near(partA.box, partB.box)) continue;

            // The shallowest axis decides the depth, as in reachPast. Along one of a's normals
            // n, a's edge and n turn with a, b's deepest vertex w stays: the depth moves with
            // n and turns by perp(n)·(origin - w). Along one of b's normals m, a's deepest
            // vertex u moves away from b's edge: by -m, and turns by -m·perp(u - origin).
            double depth = kInf;
            Point slope{};
            double turn = 0;
            for (std::size_t k = 0; k < partA.count; ++k) {
                const Point &normal = partA.normal[k]; // n
                std::size_t deepest = 0;               // w
                for (std::size_t j = 1; j < partB.count; ++j) {
                    if (dot(normal, partB.vertex[j]) < dot(normal, partB.vertex[deepest]))
                        deepest = j;
                }
                const Point &w = partB.vertex[deepest];
                const double along = dot(normal, partA.vertex[k]) - dot(normal, w);
                if (along < depth) {
                    depth = along;
                    slope = normal;
                    turn = dot(perp(normal), Point{origin.x - w.x, origin.y - w.y});
                }
            }
            for (std::size_t k = 0; k < partB.count; ++k) {
                const Point &normal = partB.normal[k]; // m
                std::size_t deepest = 0;               // u
                for (std::size_t j = 1; j < partA.count; ++j) {
                    if (dot(normal, partA.vertex[j]) < dot(normal, partA.vertex[deepest]))
                        deepest = j;
                }
                const double along =
                    dot(normal, partB.vertex[k]) - dot(normal, partA.vertex[deepest]);
                if (along < depth) {
                    depth = along;
                    slope = {-normal.x, -normal.y};
                    turn = -dot(normal, perp(fromOrigin(partA.vertex[deepest])));
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
