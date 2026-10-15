#include "packing/rotation.h"

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
        geometry::transformInto(tree::outline(),
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

    const bool smaller = layout.side() < before;
    const bool keepTurn = smaller && layout.isValid();
    if (!keepTurn) {
        for (std::size_t i = 0; i < layout.size(); ++i) {
            layout.setPose(i, saved[i]);
        }
        return false;
    }

    return true;
}
} // namespace packer::packing
