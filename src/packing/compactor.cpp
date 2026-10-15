#include "packing/compactor.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

namespace packer::packing {
namespace {
// Side of the square around the layout turned by `deg` about (cx, cy).
double sideTurnedBy(const Layout &layout, double cx, double cy, double deg, Polygon &scratch) {
    constexpr double kInf = 1e300;
    const double rad = deg * std::numbers::pi / 180.0;
    const double c = std::cos(rad);
    const double s = std::sin(rad);
    double minX = kInf, minY = kInf, maxX = -kInf, maxY = -kInf;
    for (const Pose &p : layout.poses) {
        const double dx = p.x - cx;
        const double dy = p.y - cy;
        geometry::transformInto(layout.shape,
                                {cx + dx * c - dy * s, cy + dx * s + dy * c, p.deg + deg}, scratch);
        for (const geometry::Point &q : scratch) {
            minX = std::min(minX, q.x);
            maxX = std::max(maxX, q.x);
            minY = std::min(minY, q.y);
            maxY = std::max(maxY, q.y);
        }
    }
    return std::max(maxX - minX, maxY - minY);
}
} // namespace

bool rotateWhole(Layout &layout) {
    if (layout.size() == 0) return false;
    const Box b = layout.bounds();
    const double cx = (b.minX + b.maxX) / 2;
    const double cy = (b.minY + b.maxY) / 2;
    const double before = layout.side();
    Polygon scratch;

    // A square looks the same after a 90° turn, so 0..90° covers every case.
    double bestDeg = 0;
    double bestSide = before;
    for (double deg = 0.25; deg < 90.0; deg += 0.25) {
        const double side = sideTurnedBy(layout, cx, cy, deg, scratch);
        if (side < bestSide) {
            bestSide = side;
            bestDeg = deg;
        }
    }
    // Refine around the best angle.
    for (double step = 0.125; step > 1e-5; step /= 2) {
        for (const double deg : {bestDeg - step, bestDeg + step}) {
            const double side = sideTurnedBy(layout, cx, cy, deg, scratch);
            if (side < bestSide) {
                bestSide = side;
                bestDeg = deg;
            }
        }
    }
    if (bestSide >= before - 1e-12) return false;

    const std::vector<Pose> saved = layout.poses;
    const double rad = bestDeg * std::numbers::pi / 180.0;
    const double c = std::cos(rad);
    const double s = std::sin(rad);
    for (std::size_t i = 0; i < layout.size(); ++i) {
        const Pose &p = saved[i];
        const double dx = p.x - cx;
        const double dy = p.y - cy;
        layout.setPose(i, {cx + dx * c - dy * s, cy + dx * s + dy * c,
                           geometry::normalizeDeg(p.deg + bestDeg)});
    }
    if (!layout.isValid() || layout.side() >= before) {
        for (std::size_t i = 0; i < layout.size(); ++i) layout.setPose(i, saved[i]);
        return false;
    }
    return true;
}

void compact(Layout &layout, std::mt19937_64 &rng, long repairIters) {
    const std::size_t n = layout.size();
    if (n < 2) return;

    std::uniform_real_distribution<double> unit(0.0, 1.0);
    auto jiggle = [&](double scale) {
        return scale * (2 * unit(rng) - 1);
    };

    double shrink = 0.01; // fraction of the side removed per attempt
    int failures = 0;
    while (shrink > 2e-4 && failures < 6) {
        const std::vector<Pose> saved = layout.poses;
        const Box b = layout.bounds();
        const double side = std::max(b.maxX - b.minX, b.maxY - b.minY);
        const double cx = (b.minX + b.maxX) / 2;
        const double cy = (b.minY + b.maxY) / 2;
        const double f = 1 - shrink;
        // The new square: trees may not leave it during repair.
        const Box target{cx - side * f / 2, cy - side * f / 2, cx + side * f / 2,
                         cy + side * f / 2};

        for (std::size_t i = 0; i < n; ++i) {
            const Pose &p = saved[i];
            layout.setPose(i, {cx + (p.x - cx) * f, cy + (p.y - cy) * f, p.deg});
        }
        // Scaling toward the centre can leave a tree poking out at the edge: that is the
        // part of the layout that does not fit, so treat it as a failure unless repaired.
        auto inside = [&](const Box &t) {
            return t.minX >= target.minX && t.maxX <= target.maxX && t.minY >= target.minY &&
                   t.maxY <= target.maxY;
        };

        std::vector<int> hits(n, 0);
        for (std::size_t i = 0; i < n; ++i) {
            hits[i] = static_cast<int>(layout.collidingWith(i).size());
        }
        std::vector<std::size_t> bad;
        auto rebuildBad = [&] {
            bad.clear();
            for (std::size_t i = 0; i < n; ++i) {
                if (hits[i] > 0 || !inside(layout.boxes[i])) bad.push_back(i);
            }
        };
        rebuildBad();

        const long budget = repairIters * static_cast<long>(n);
        for (long it = 0; it < budget && !bad.empty(); ++it) {
            const double progress = static_cast<double>(it) / budget;
            const double step = side * 0.02 * (1 - 0.9 * progress);
            const std::size_t i = bad[std::uniform_int_distribution<std::size_t>(0, bad.size() - 1)(rng)];

            Pose moved = layout.poses[i];
            const Pose old = moved;
            moved.x += jiggle(step);
            moved.y += jiggle(step);
            if (unit(rng) < 0.5) moved.deg += jiggle(300 * step / side);
            moved.deg = geometry::normalizeDeg(moved.deg);

            const std::vector<std::size_t> oldHits = layout.collidingWith(i);
            layout.setPose(i, moved);
            const std::vector<std::size_t> newHits = layout.collidingWith(i);
            // Sideways moves are accepted too, so overlaps can drift apart. A tree must end up
            // inside the new square and never in a worse spot than before.
            if (inside(layout.boxes[i]) && newHits.size() <= oldHits.size()) {
                for (std::size_t j : oldHits) --hits[j];
                for (std::size_t j : newHits) ++hits[j];
                hits[i] = static_cast<int>(newHits.size());
                rebuildBad();
            } else {
                layout.setPose(i, old);
            }
        }

        if (bad.empty() && layout.isValid() && layout.side() < side) {
            shrink = std::min(shrink * 1.5, 0.02);
            failures = 0;
        } else {
            for (std::size_t i = 0; i < n; ++i) layout.setPose(i, saved[i]);
            shrink /= 2;
            ++failures;
        }
    }
}
} // namespace packer::packing
