#include "packing/annealer.h"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

namespace packer::packing {
namespace {
// Square side, plus a small reward for trees sitting near the centre. The pull opens up
// room inside so trees on the edge can move in later.
double energy(const Layout &layout, double pullWeight, double &side) {
    const Box b = layout.bounds();
    const double cx = (b.minX + b.maxX) / 2;
    const double cy = (b.minY + b.maxY) / 2;
    double pull = 0;
    for (const Box &t : layout.boxes) {
        pull += std::max({t.maxX - cx, cx - t.minX, t.maxY - cy, cy - t.minY});
    }
    side = std::max(b.maxX - b.minX, b.maxY - b.minY);
    return side + pullWeight * pull / layout.size();
}
} // namespace

void anneal(Layout &layout, const AnnealConfig &cfg, std::mt19937_64 &rng) {
    if (layout.size() == 0 || cfg.iterations <= 0) return;

    std::uniform_int_distribution<std::size_t> pickTree(0, layout.size() - 1);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    auto jiggle = [&](double scale) { return scale * (2 * unit(rng) - 1); };

    // What the current move changed, so it can be undone.
    std::vector<std::pair<std::size_t, Pose>> undo;
    auto moveTree = [&](std::size_t i, Pose pose) {
        pose.deg = geometry::normalizeDeg(pose.deg);
        undo.push_back({i, layout.poses[i]});
        layout.setPose(i, pose);
    };
    auto revert = [&] {
        for (auto it = undo.rbegin(); it != undo.rend(); ++it) {
            layout.setPose(it->first, it->second);
        }
    };

    const double edgeEnd = cfg.squeezeShare + cfg.edgeShare;
    const double relocateEnd = edgeEnd + cfg.relocateShare;

    double currentSide = 0;
    double current = energy(layout, cfg.pullWeight, currentSide);
    // Only the poses are saved; the cached polygons are rebuilt once at the end.
    std::vector<Pose> bestPoses = layout.poses;
    double bestSide = currentSide;

    for (long it = 0; it < cfg.iterations; ++it) {
        // Progress 0 -> 1. Temperature and step size shrink geometrically.
        const double t = static_cast<double>(it) / cfg.iterations;
        const double temp = cfg.startTemp * std::pow(cfg.endTemp / cfg.startTemp, t);
        const double step = cfg.startStep * std::pow(cfg.endStep / cfg.startStep, t);

        undo.clear();
        const Box b = layout.bounds();
        const double cx = (b.minX + b.maxX) / 2;
        const double cy = (b.minY + b.maxY) / 2;
        const double r = unit(rng);
        bool legal = false;

        if (r < cfg.squeezeShare) {
            // Squeeze: pull every tree toward the centre at once, keeping the angles.
            const double side = std::max(b.maxX - b.minX, b.maxY - b.minY);
            const double f = 1 - step * unit(rng) / side;
            for (std::size_t i = 0; i < layout.size(); ++i) {
                const Pose p = layout.poses[i];
                moveTree(i, {cx + (p.x - cx) * f, cy + (p.y - cy) * f, p.deg});
            }
            legal = layout.isValid();
        } else {
            std::size_t i = 0;
            Pose moved;
            if (r < relocateEnd) {
                // Only trees on the edge decide the side, so work on those.
                const std::vector<std::size_t> edge = layout.edgeTrees();
                i = edge[std::uniform_int_distribution<std::size_t>(0, edge.size() - 1)(rng)];
                moved = layout.poses[i];
                if (r < edgeEnd) {
                    // Edge push: a small step toward the centre, with a small turn.
                    const double dx = cx - moved.x;
                    const double dy = cy - moved.y;
                    const double len = std::hypot(dx, dy);
                    if (len > 0) {
                        const double d = step * unit(rng);
                        moved.x += d * dx / len;
                        moved.y += d * dy / len;
                    }
                    moved.deg += jiggle(30 * step);
                } else {
                    // Relocate: drop it at a random spot and angle inside the square.
                    moved = {b.minX + unit(rng) * (b.maxX - b.minX),
                             b.minY + unit(rng) * (b.maxY - b.minY), 360 * unit(rng)};
                }
            } else {
                i = pickTree(rng);
                moved = layout.poses[i];
                const double kind = unit(rng);
                if (kind < 0.70) {
                    moved.x += jiggle(step);
                    moved.y += jiggle(step);
                } else if (kind < 0.95) {
                    moved.deg += jiggle(100 * step);
                } else {
                    moved.deg += 180;
                }
            }
            moveTree(i, moved);
            legal = !layout.collidesWithOthers(i);
        }

        if (!legal) { // a collision: undo right away
            revert();
            continue;
        }

        double nextSide = 0;
        const double next = energy(layout, cfg.pullWeight, nextSide);
        if (next <= current || unit(rng) < std::exp((current - next) / temp)) {
            current = next;
            if (nextSide < bestSide) {
                bestSide = nextSide;
                bestPoses = layout.poses;
            }
        } else {
            revert();
        }
    }
    for (std::size_t i = 0; i < layout.size(); ++i) {
        layout.setPose(i, bestPoses[i]);
    }
}
} // namespace packer::packing
