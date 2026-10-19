#include "packing/nfp.h"

#include <algorithm>
#include <cmath>

namespace packer::packing {
namespace {
using geometry::kInf;
using geometry::outside;
using geometry::Point;

// > 0 if b is left of the line o->a, < 0 if right, 0 if on it.
double cross(const Point &o, const Point &a, const Point &b) {
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

// Convex hull, counter-clockwise (Andrew's monotone chain).
std::vector<Point> hull(std::vector<Point> points) {
    std::sort(points.begin(), points.end(), [](const Point &a, const Point &b) {
        return a.x < b.x || (a.x == b.x && a.y < b.y);
    });
    if (points.size() < 3) return points;

    std::vector<Point> corners(2 * points.size());
    std::size_t count = 0;
    for (const Point &point : points) {
        while (count >= 2 && cross(corners[count - 2], corners[count - 1], point) <= 0)
            --count;
        corners[count++] = point;
    }
    const std::size_t lower = count + 1;
    for (std::size_t i = points.size() - 1; i-- > 0;) {
        while (count >= lower && cross(corners[count - 2], corners[count - 1], points[i]) <= 0) {
            --count;
        }
        corners[count++] = points[i];
    }
    corners.resize(count - 1);
    return corners;
}
} // namespace

std::vector<Point> nfpCorners(const tree::Placed &fixed, double deg, double gap) {
    const tree::Placed moving = tree::place({0, 0, deg});
    std::vector<Point> corners;
    std::vector<Point> sums;

    for (const tree::Part &a : fixed.parts) {
        for (const tree::Part &b : moving.parts) {
            // The moving origin at t puts b's vertex w at w + t: a and b + t overlap exactly
            // when t lies in a ⊕ (−b).
            sums.clear();
            for (std::size_t i = 0; i < a.count; ++i) {
                for (std::size_t j = 0; j < b.count; ++j) {
                    sums.push_back({a.vertex[i].x - b.vertex[j].x, a.vertex[i].y - b.vertex[j].y});
                }
            }
            const std::vector<Point> polygon = hull(sums);

            Point centroid{};
            for (const Point &corner : polygon) {
                centroid.x += corner.x / polygon.size();
                centroid.y += corner.y / polygon.size();
            }
            for (const Point &corner : polygon) {
                const double len = std::hypot(corner.x - centroid.x, corner.y - centroid.y);
                if (len == 0) continue;
                corners.push_back({corner.x + gap * (corner.x - centroid.x) / len,
                                   corner.y + gap * (corner.y - centroid.y) / len});
            }
        }
    }

    return corners;
}

Pose bestInsertion(const Layout &layout, const Box &square, int angles, std::mt19937_64 &rng,
                   const Box *region) {
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    // Where the new tree's origin may go: the square, or the part of it inside the region.
    const Box area =
        region == nullptr
            ? square
            : Box{std::max(square.minX, region->minX), std::max(square.minY, region->minY),
                  std::min(square.maxX, region->maxX), std::min(square.maxY, region->maxY)};
    const double side = std::max(area.maxX - area.minX, area.maxY - area.minY);
    const auto [cx, cy] = geometry::centre(area);

    Layout trial = layout;
    trial.add({cx, cy, 0});
    const std::size_t last = trial.size() - 1;

    auto penalty = [&](const Pose &pose) {
        trial.setPose(last, pose);
        const double stickOut = outside(trial.boxes[last], square);
        double sum = stickOut * stickOut;
        trial.anyNeighbour(last, [&](std::size_t j) {
            const double depth = tree::penetration(trial.placed[last], trial.placed[j]);
            sum += depth * depth;
            return false;
        });
        return sum;
    };

    // Only spots with the origin near the square can be any good.
    auto nearSquare = [&](const Point &point) {
        return point.x > area.minX - 0.5 && point.x < area.maxX + 0.5 &&
               point.y > area.minY - 0.5 && point.y < area.maxY + 0.5;
    };
    // Only trees that can touch a tree with its origin there have corners worth trying.
    constexpr double kTouch = 2 * 0.8; // twice the farthest a vertex is from the origin
    auto canTouch = [&](const Box &box) {
        return box.maxX > area.minX - kTouch && box.minX < area.maxX + kTouch &&
               box.maxY > area.minY - kTouch && box.minY < area.maxY + kTouch;
    };

    Pose best{cx, cy, 0};
    double least = kInf;
    const int count = std::max(1, angles);
    const double offset = 360.0 * unit(rng);
    const int grid = std::max(4, static_cast<int>(std::ceil(side / 0.2)));

    for (int k = 0; k < count && least > 0; ++k) {
        const double deg = geometry::normalizeDeg(offset + k * 360.0 / count);
        auto tryAt = [&](const Point &point) {
            if (!nearSquare(point)) return;
            const double cost = penalty({point.x, point.y, deg});
            if (cost < least) {
                least = cost;
                best = {point.x, point.y, deg};
            }
        };

        for (std::size_t j = 0; j < layout.size() && least > 0; ++j) {
            if (!canTouch(layout.boxes[j])) continue;
            for (const Point &corner : nfpCorners(layout.placed[j], deg, 1e-7)) {
                tryAt(corner);
            }
        }
        for (int gx = 0; gx <= grid && least > 0; ++gx) {
            for (int gy = 0; gy <= grid && least > 0; ++gy) {
                tryAt({area.minX + side * gx / grid, area.minY + side * gy / grid});
            }
        }
    }

    return best;
}

Layout ruinAndRecreate(const Layout &layout, const Box &square, std::size_t count,
                       std::mt19937_64 &rng) {
    Layout result = layout;
    const std::size_t n = result.size();
    count = std::min(count, n - 1);
    if (count == 0) return result;
    std::uniform_real_distribution<double> unit(0.0, 1.0);

    // The patch: around a tree on the edge of the square (where the side is decided) half the
    // time, around any tree otherwise.
    const std::vector<std::size_t> edge = result.edgeTrees();
    const bool atEdge = !edge.empty() && unit(rng) < 0.5;
    const std::size_t seedTree =
        atEdge ? edge[std::uniform_int_distribution<std::size_t>(0, edge.size() - 1)(rng)]
               : std::uniform_int_distribution<std::size_t>(0, n - 1)(rng);
    const Point centre = result.placed[seedTree].centre;

    std::vector<std::size_t> order(n);
    for (std::size_t i = 0; i < n; ++i)
        order[i] = i;
    auto distSq = [&](std::size_t i) {
        const Point &point = result.placed[i].centre;
        return (point.x - centre.x) * (point.x - centre.x) +
               (point.y - centre.y) * (point.y - centre.y);
    };
    const auto cut = order.begin() + static_cast<std::ptrdiff_t>(count);
    std::partial_sort(order.begin(), cut, order.end(),
                      [&](std::size_t a, std::size_t b) { return distSq(a) < distSq(b); });

    // Remove from the highest index down, so the others keep theirs.
    std::vector<std::size_t> gone(order.begin(), cut);
    double reach = 0;
    for (const std::size_t i : gone)
        reach = std::max(reach, std::sqrt(distSq(i)));
    std::sort(gone.rbegin(), gone.rend());
    for (const std::size_t i : gone)
        result.remove(i);

    // Put them back one by one at the best spot within the patch (plus a little room).
    const double radius = reach + 0.4;
    const Box patch{centre.x - radius, centre.y - radius, centre.x + radius, centre.y + radius};
    for (std::size_t k = 0; k < count; ++k)
        result.add(bestInsertion(result, square, 6, rng, &patch));

    return result;
}
} // namespace packer::packing
