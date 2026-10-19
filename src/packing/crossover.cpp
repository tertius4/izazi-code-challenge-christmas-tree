#include "packing/crossover.h"

#include <cmath>
#include <numbers>

#include "packing/nfp.h"

namespace packer::packing {
namespace {
using geometry::Point;

// How much tree i overlaps the others: its summed penetration.
double overlapOf(const Layout &layout, std::size_t i) {
    double sum = 0;
    layout.anyNeighbour(i, [&](std::size_t j) {
        sum += tree::penetration(layout.placed[i], layout.placed[j]);
        return false;
    });
    return sum;
}
} // namespace

Layout crossover(const Layout &a, const Layout &b, int n, std::mt19937_64 &rng) {
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const Point centreA = geometry::centre(a.bounds());
    const Point centreB = geometry::centre(b.bounds());
    const double side = a.side();

    // b turned by a quarter turns (the square looks the same), maybe mirrored left to right
    // (a mirrored tree is the tree turned the other way), and centred on a.
    const int quarters = static_cast<int>(4 * unit(rng));
    const bool mirror = unit(rng) < 0.5;
    const double rad = quarters * std::numbers::pi / 2;
    const double cosine = std::cos(rad), sine = std::sin(rad);
    auto fromB = [&](const Pose &pose) {
        double dx = pose.x - centreB.x, dy = pose.y - centreB.y, deg = pose.deg;
        if (mirror) {
            dx = -dx;
            deg = -deg;
        }
        return Pose{centreA.x + dx * cosine - dy * sine, centreA.y + dx * sine + dy * cosine,
                    deg + quarters * 90.0};
    };

    // A line through a point near the centre, at a random angle.
    const double angle = 2 * std::numbers::pi * unit(rng);
    const Point normal{std::cos(angle), std::sin(angle)};
    const double offset = (unit(rng) - 0.5) * 0.5 * side;
    auto onASide = [&](const Pose &pose) {
        return (pose.x - centreA.x) * normal.x + (pose.y - centreA.y) * normal.y < offset;
    };

    Layout child;
    for (const Pose &pose : a.poses) {
        if (onASide(pose)) child.add(pose);
    }
    for (const Pose &pose : b.poses) {
        const Pose moved = fromB(pose);
        if (!onASide(moved)) child.add(moved);
    }

    // Too many: drop the most overlapped tree, or if none overlaps, the one that shrinks the
    // square most.
    while (static_cast<int>(child.size()) > n) {
        std::size_t worst = 0;
        double most = 0;
        for (std::size_t i = 0; i < child.size(); ++i) {
            const double overlap = overlapOf(child, i);
            if (overlap > most) {
                most = overlap;
                worst = i;
            }
        }
        child.remove(most > 0 ? worst : child.bestTreeToRemove());
    }

    // Too few: add trees at the best spots in a's square.
    const Box square{centreA.x - side / 2, centreA.y - side / 2, centreA.x + side / 2,
                     centreA.y + side / 2};
    while (static_cast<int>(child.size()) < n) {
        if (child.size() == 0) {
            child.add({centreA.x, centreA.y, 360 * unit(rng)});
            continue;
        }
        child.add(bestInsertion(child, square, 4, rng));
    }

    return child;
}
} // namespace packer::packing
