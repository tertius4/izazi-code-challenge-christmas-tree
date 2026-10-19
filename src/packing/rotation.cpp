#include "packing/rotation.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

namespace packer::packing {
namespace {
using geometry::kInf;

// Side of the square around the layout turned by `deg` about (cx, cy).
double sideTurnedBy(const Layout &layout, double cx, double cy, double deg, Polygon &scratch) {
    const double rad = deg * std::numbers::pi / 180.0;
    const double cosine = std::cos(rad);
    const double sine = std::sin(rad);
    double minX = kInf, minY = kInf, maxX = -kInf, maxY = -kInf;

    for (const Pose &pose : layout.poses) {
        const double dx = pose.x - cx;
        const double dy = pose.y - cy;
        const Pose turned{cx + dx * cosine - dy * sine, cy + dx * sine + dy * cosine,
                          pose.deg + deg};
        geometry::transformInto(tree::outline(), turned, scratch);

        for (const geometry::Point &vertex : scratch) {
            minX = std::min(minX, vertex.x);
            maxX = std::max(maxX, vertex.x);
            minY = std::min(minY, vertex.y);
            maxY = std::max(maxY, vertex.y);
        }
    }

    return std::max(maxX - minX, maxY - minY);
}
} // namespace

bool rotateWhole(Layout &layout) {
    if (layout.size() == 0) return false;

    const auto [cx, cy] = geometry::centre(layout.bounds());
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
    const double cosine = std::cos(rad);
    const double sine = std::sin(rad);

    for (std::size_t i = 0; i < layout.size(); ++i) {
        const Pose &pose = saved[i];
        const double dx = pose.x - cx;
        const double dy = pose.y - cy;
        layout.setPose(i, {cx + dx * cosine - dy * sine, cy + dx * sine + dy * cosine,
                           geometry::normalizeDeg(pose.deg + bestDeg)});
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
