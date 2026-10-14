#include "packing/annealer.h"

#include <algorithm>
#include <cmath>

namespace packer::packing {
namespace {
// Square side, plus a small reward for trees sitting near the centre. The pull opens up
// room inside so trees on the edge can move in later.
double energy(const Layout &layout, double pullWeight) {
    const Box b = layout.bounds();
    const double cx = (b.minX + b.maxX) / 2;
    const double cy = (b.minY + b.maxY) / 2;
    double pull = 0;
    for (const Box &t : layout.boxes) {
        pull += std::max({t.maxX - cx, cx - t.minX, t.maxY - cy, cy - t.minY});
    }
    const double side = std::max(b.maxX - b.minX, b.maxY - b.minY);
    return side + pullWeight * pull / layout.size();
}
} // namespace

void anneal(Layout &layout, const AnnealConfig &cfg, std::mt19937_64 &rng) {
    if (layout.size() == 0 || cfg.iterations <= 0) return;
    std::uniform_int_distribution<std::size_t> pickTree(0, layout.size() - 1);
    std::uniform_real_distribution<double> unit(0.0, 1.0);

    double current = energy(layout, cfg.pullWeight);
    Layout best = layout;
    double bestSide = layout.side();

    for (long it = 0; it < cfg.iterations; ++it) {
        // Progress 0 -> 1. Temperature and step size shrink geometrically.
        const double t = static_cast<double>(it) / cfg.iterations;
        const double temp = cfg.startTemp * std::pow(cfg.endTemp / cfg.startTemp, t);
        const double step = cfg.startStep * std::pow(cfg.endStep / cfg.startStep, t);

        const std::size_t i = pickTree(rng);
        const Pose old = layout.poses[i];
        Pose moved = old;
        const double r = unit(rng);
        if (r < 0.70) {
            moved.x += step * (2 * unit(rng) - 1);
            moved.y += step * (2 * unit(rng) - 1);
        } else if (r < 0.95) {
            moved.deg += 100 * step * (2 * unit(rng) - 1);
        } else {
            moved.deg += 180;
        }
        moved.deg = geometry::normalizeDeg(moved.deg);

        layout.setPose(i, moved);
        if (layout.collidesWithOthers(i)) { // illegal: undo right away
            layout.setPose(i, old);
            continue;
        }

        const double next = energy(layout, cfg.pullWeight);
        if (next <= current || unit(rng) < std::exp((current - next) / temp)) {
            current = next;
            if (layout.side() < bestSide) {
                bestSide = layout.side();
                best = layout;
            }
        } else {
            layout.setPose(i, old);
        }
    }
    layout = std::move(best);
}
} // namespace packer::packing
