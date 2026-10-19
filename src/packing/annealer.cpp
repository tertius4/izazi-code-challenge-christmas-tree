#include "packing/annealer.h"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

namespace packer::packing {
namespace {
// How far a tree's box reaches from the centre (cx, cy), along x or y.
double reach(const Box &t, double cx, double cy) {
    return std::max({t.maxX - cx, cx - t.minX, t.maxY - cy, cy - t.minY});
}

// The square side, plus a small reward for trees sitting near the centre. The pull opens up
// room inside so trees on the edge can move in later.
struct Energy {
    Box box{};       // around all trees
    double pull = 0; // sum of reach() over all trees
    double side = 0;
    double value = 0; // what the annealer minimises
};

Energy measure(const Layout &layout, double pullWeight) {
    Energy e;
    e.box = layout.bounds();
    const auto [cx, cy] = geometry::centre(e.box);
    for (const Box &t : layout.boxes) {
        e.pull += reach(t, cx, cy);
    }

    e.side = std::max(e.box.maxX - e.box.minX, e.box.maxY - e.box.minY);
    e.value = e.side + pullWeight * e.pull / static_cast<double>(layout.size());

    return e;
}

bool strictlyInside(const Box &t, const Box &b) {
    return t.minX > b.minX && t.maxX < b.maxX && t.minY > b.minY && t.maxY < b.maxY;
}

// The energy after one tree's box changed from `before` to `after`. Most moves stay clear of
// the edge, so the box and its centre stay put and only that tree's pull changes. Anything
// else is measured from scratch.
Energy afterMove(const Energy &e, const Box &before, const Box &after, const Layout &layout,
                 double pullWeight) {
    const bool wasInside = strictlyInside(before, e.box);
    const bool isInside = strictlyInside(after, e.box);
    if (!wasInside || !isInside) {
        return measure(layout, pullWeight);
    }

    const auto [cx, cy] = geometry::centre(e.box);
    Energy next = e;
    next.pull += reach(after, cx, cy) - reach(before, cx, cy);
    next.value = next.side + pullWeight * next.pull / static_cast<double>(layout.size());

    return next;
}
} // namespace

void anneal(Layout &layout, const AnnealConfig &cfg, std::mt19937_64 &rng) {
    const bool nothingToDo = layout.size() == 0 || cfg.iterations <= 0;
    if (nothingToDo) return;

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

    Energy current = measure(layout, cfg.pullWeight);
    // Only the poses are saved; the cached parts are rebuilt once at the end.
    std::vector<Pose> bestPoses = layout.poses;
    double bestSide = current.side;

    for (long it = 0; it < cfg.iterations; ++it) {
        // Look for a stop request only now and then: the atomic load is not free.
        const bool checkStop = cfg.stop && it % 4096 == 0;
        if (checkStop && cfg.stop->load(std::memory_order_relaxed)) break;

        // Progress 0 -> 1. Temperature and step size shrink geometrically.
        const double t = static_cast<double>(it) / cfg.iterations;
        const double temp = cfg.startTemp * std::pow(cfg.endTemp / cfg.startTemp, t);
        const double step = cfg.startStep * std::pow(cfg.endStep / cfg.startStep, t);

        undo.clear();
        const Box b = current.box;
        const auto [cx, cy] = geometry::centre(b);
        const double r = unit(rng);
        bool legal = false;
        bool single = false; // a one-tree move: tree `i`, whose box was `before`
        std::size_t i = 0;
        Box before{};

        if (r < cfg.squeezeShare) {
            // Squeeze: pull every tree toward the centre at once, keeping the angles.
            const double side = std::max(b.maxX - b.minX, b.maxY - b.minY);
            const double f = 1 - step * unit(rng) / side;

            for (std::size_t k = 0; k < layout.size(); ++k) {
                const Pose p = layout.poses[k];
                moveTree(k, {cx + (p.x - cx) * f, cy + (p.y - cy) * f, p.deg});
            }

            legal = !layout.anyCollision();
        } else {
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

            before = layout.boxes[i];
            moveTree(i, moved);
            legal = !layout.collidesWithOthers(i);
            single = true;
        }

        if (!legal) { // a collision: undo right away
            revert();
            continue;
        }

        const Energy next =
            single ? afterMove(current, before, layout.boxes[i], layout, cfg.pullWeight)
                   : measure(layout, cfg.pullWeight);

        // Always accept a move that is no worse; accept a worse one with a chance that shrinks
        // as it gets worse and as the temperature drops.
        const double worseBy = next.value - current.value;
        const bool noWorse = worseBy <= 0;
        const bool accept = noWorse || unit(rng) < std::exp(-worseBy / temp);
        if (accept) {
            current = next;
            if (next.side < bestSide) {
                bestSide = next.side;
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
