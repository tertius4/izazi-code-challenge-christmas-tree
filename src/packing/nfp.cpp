#include "packing/nfp.h"

#include <algorithm>
#include <cmath>

namespace packer::packing {
namespace {
using geometry::outside;
using geometry::Point;

double cross(const Point &o, const Point &a, const Point &b) {
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

// Convex hull, counter-clockwise (Andrew's monotone chain).
std::vector<Point> hull(std::vector<Point> pts) {
    std::sort(pts.begin(), pts.end(),
              [](const Point &a, const Point &b) { return a.x < b.x || (a.x == b.x && a.y < b.y); });
    if (pts.size() < 3) return pts;

    std::vector<Point> h(2 * pts.size());
    std::size_t k = 0;
    for (const Point &p : pts) {
        while (k >= 2 && cross(h[k - 2], h[k - 1], p) <= 0) --k;
        h[k++] = p;
    }
    const std::size_t lower = k + 1;
    for (std::size_t i = pts.size() - 1; i-- > 0;) {
        while (k >= lower && cross(h[k - 2], h[k - 1], pts[i]) <= 0) --k;
        h[k++] = pts[i];
    }
    h.resize(k - 1);
    return h;
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
                    sums.push_back({a.v[i].x - b.v[j].x, a.v[i].y - b.v[j].y});
                }
            }
            const std::vector<Point> h = hull(sums);

            Point c{};
            for (const Point &p : h) {
                c.x += p.x / h.size();
                c.y += p.y / h.size();
            }
            for (const Point &p : h) {
                const double len = std::hypot(p.x - c.x, p.y - c.y);
                if (len == 0) continue;
                corners.push_back({p.x + gap * (p.x - c.x) / len, p.y + gap * (p.y - c.y) / len});
            }
        }
    }

    return corners;
}

Pose bestInsertion(const Layout &layout, const Box &square, int angles, std::mt19937_64 &rng,
                   const Box *region) {
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    // Where the new tree's origin may go: the square, or the part of it inside the region.
    const Box area = region == nullptr
                         ? square
                         : Box{std::max(square.minX, region->minX), std::max(square.minY, region->minY),
                               std::min(square.maxX, region->maxX), std::min(square.maxY, region->maxY)};
    const double side = std::max(area.maxX - area.minX, area.maxY - area.minY);
    const auto [cx, cy] = geometry::centre(area);

    Layout l = layout;
    l.add({cx, cy, 0});
    const std::size_t last = l.size() - 1;

    auto penalty = [&](const Pose &p) {
        l.setPose(last, p);
        const double out = outside(l.boxes[last], square);
        double sum = out * out;
        l.anyNeighbour(last, [&](std::size_t j) {
            const double d = tree::penetration(l.placed[last], l.placed[j]);
            sum += d * d;
            return false;
        });
        return sum;
    };

    // Only spots with the origin near the square can be any good.
    auto nearSquare = [&](const Point &p) {
        return p.x > area.minX - 0.5 && p.x < area.maxX + 0.5 && p.y > area.minY - 0.5 &&
               p.y < area.maxY + 0.5;
    };
    // Only trees that can touch a tree with its origin there have corners worth trying.
    constexpr double kTouch = 2 * 0.8; // twice the farthest a vertex is from the origin
    auto canTouch = [&](const Box &b) {
        return b.maxX > area.minX - kTouch && b.minX < area.maxX + kTouch &&
               b.maxY > area.minY - kTouch && b.minY < area.maxY + kTouch;
    };

    Pose best{cx, cy, 0};
    double least = INFINITY;
    const int count = std::max(1, angles);
    const double offset = 360.0 * unit(rng);
    const int grid = std::max(4, static_cast<int>(std::ceil(side / 0.2)));

    for (int a = 0; a < count && least > 0; ++a) {
        const double deg = geometry::normalizeDeg(offset + a * 360.0 / count);
        auto tryAt = [&](const Point &p) {
            if (!nearSquare(p)) return;
            const double cost = penalty({p.x, p.y, deg});
            if (cost < least) {
                least = cost;
                best = {p.x, p.y, deg};
            }
        };

        for (std::size_t j = 0; j < layout.size() && least > 0; ++j) {
            if (!canTouch(layout.boxes[j])) continue;
            for (const Point &p : nfpCorners(layout.placed[j], deg, 1e-7)) {
                tryAt(p);
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
    Layout l = layout;
    const std::size_t n = l.size();
    count = std::min(count, n - 1);
    if (count == 0) return l;
    std::uniform_real_distribution<double> unit(0.0, 1.0);

    // The patch: around a tree on the edge of the square (where the side is decided) half the
    // time, around any tree otherwise.
    const std::vector<std::size_t> edge = l.edgeTrees();
    const bool atEdge = !edge.empty() && unit(rng) < 0.5;
    const std::size_t centre =
        atEdge ? edge[std::uniform_int_distribution<std::size_t>(0, edge.size() - 1)(rng)]
               : std::uniform_int_distribution<std::size_t>(0, n - 1)(rng);
    const Point c = l.placed[centre].centre;

    std::vector<std::size_t> order(n);
    for (std::size_t i = 0; i < n; ++i) order[i] = i;
    auto dist = [&](std::size_t i) {
        const Point &p = l.placed[i].centre;
        return (p.x - c.x) * (p.x - c.x) + (p.y - c.y) * (p.y - c.y);
    };
    std::partial_sort(order.begin(), order.begin() + static_cast<std::ptrdiff_t>(count), order.end(),
                      [&](std::size_t a, std::size_t b) { return dist(a) < dist(b); });

    // Remove from the highest index down, so the others keep theirs.
    std::vector<std::size_t> gone(order.begin(), order.begin() + static_cast<std::ptrdiff_t>(count));
    double reach = 0;
    for (const std::size_t i : gone) reach = std::max(reach, std::sqrt(dist(i)));
    std::sort(gone.rbegin(), gone.rend());
    for (const std::size_t i : gone) l.remove(i);

    // Put them back one by one at the best spot within the patch (plus a little room).
    const double r = reach + 0.4;
    const Box patch{c.x - r, c.y - r, c.x + r, c.y + r};
    for (std::size_t k = 0; k < count; ++k) l.add(bestInsertion(l, square, 6, rng, &patch));

    return l;
}
} // namespace packer::packing
