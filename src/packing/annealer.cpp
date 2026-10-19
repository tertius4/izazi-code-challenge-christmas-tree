#include "packing/annealer.h"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

namespace packer::packing {
namespace {
// How far a tree's box reaches from the centre (cx, cy), along x or y.
double reach(const Box &box, double cx, double cy) {
    return std::max({box.maxX - cx, cx - box.minX, box.maxY - cy, cy - box.minY});
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
    Energy energy;
    energy.box = layout.bounds();
    const auto [cx, cy] = geometry::centre(energy.box);
    for (const Box &box : layout.boxes) {
        energy.pull += reach(box, cx, cy);
    }

    energy.side = std::max(energy.box.maxX - energy.box.minX, energy.box.maxY - energy.box.minY);
    energy.value = energy.side + pullWeight * energy.pull / static_cast<double>(layout.size());

    return energy;
}

bool strictlyInside(const Box &inner, const Box &outer) {
    return inner.minX > outer.minX && inner.maxX < outer.maxX && inner.minY > outer.minY &&
           inner.maxY < outer.maxY;
}

// The energy after one tree's box changed from `before` to `after`. Most moves stay clear of
// the edge, so the box and its centre stay put and only that tree's pull changes. Anything
// else is measured from scratch.
Energy afterMove(const Energy &energy, const Box &before, const Box &after, const Layout &layout,
                 double pullWeight) {
    const bool wasInside = strictlyInside(before, energy.box);
    const bool isInside = strictlyInside(after, energy.box);
    if (!wasInside || !isInside) {
        return measure(layout, pullWeight);
    }

    const auto [cx, cy] = geometry::centre(energy.box);
    Energy next = energy;
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
        const double progress = static_cast<double>(it) / cfg.iterations;
        const double temp = cfg.startTemp * std::pow(cfg.endTemp / cfg.startTemp, progress);
        const double step = cfg.startStep * std::pow(cfg.endStep / cfg.startStep, progress);

        undo.clear();
        const Box bounds = current.box;
        const auto [cx, cy] = geometry::centre(bounds);
        const double roll = unit(rng); // picks the kind of move
        bool legal = false;
        bool single = false; // a one-tree move: tree `i`, whose box was `before`
        std::size_t i = 0;
        Box before{};

        if (roll < cfg.squeezeShare) {
            // Squeeze: pull every tree toward the centre at once, keeping the angles.
            const double side = std::max(bounds.maxX - bounds.minX, bounds.maxY - bounds.minY);
            const double scale = 1 - step * unit(rng) / side;

            for (std::size_t k = 0; k < layout.size(); ++k) {
                const Pose pose = layout.poses[k];
                moveTree(k, {cx + (pose.x - cx) * scale, cy + (pose.y - cy) * scale, pose.deg});
            }

            legal = !layout.anyCollision();
        } else {
            Pose moved;
            if (roll < relocateEnd) {
                // Only trees on the edge decide the side, so work on those.
                const std::vector<std::size_t> edge = layout.edgeTrees();
                i = edge[std::uniform_int_distribution<std::size_t>(0, edge.size() - 1)(rng)];
                moved = layout.poses[i];

                if (roll < edgeEnd) {
                    // Edge push: a small step toward the centre, with a small turn.
                    const double dx = cx - moved.x;
                    const double dy = cy - moved.y;
                    const double len = std::hypot(dx, dy);

                    if (len > 0) {
                        const double dist = step * unit(rng);
                        moved.x += dist * dx / len;
                        moved.y += dist * dy / len;
                    }

                    moved.deg += jiggle(30 * step);
                } else {
                    // Relocate: drop it at a random spot and angle inside the square.
                    moved = {bounds.minX + unit(rng) * (bounds.maxX - bounds.minX),
                             bounds.minY + unit(rng) * (bounds.maxY - bounds.minY),
                             360 * unit(rng)};
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
