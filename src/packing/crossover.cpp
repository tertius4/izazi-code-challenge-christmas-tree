#include "packing/crossover.h"

#include <cmath>
#include <numbers>

#include "packing/nfp.h"

namespace packer::packing {
namespace {
using geometry::Point;

Point centreOf(const Layout &l) {
    const Box b = l.bounds();
    return {(b.minX + b.maxX) / 2, (b.minY + b.maxY) / 2};
}

// How much tree i overlaps the others: its summed penetration.
double overlapOf(const Layout &l, std::size_t i) {
    double sum = 0;
    l.anyNeighbour(i, [&](std::size_t j) {
        sum += tree::penetration(l.placed[i], l.placed[j]);
        return false;
    });
    return sum;
}
} // namespace

Layout crossover(const Layout &a, const Layout &b, int n, std::mt19937_64 &rng) {
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const Point ca = centreOf(a);
    const Point cb = centreOf(b);
    const double side = a.side();

    // b turned by a quarter turns (the square looks the same), maybe mirrored left to right
    // (a mirrored tree is the tree turned the other way), and centred on a.
    const int quarters = static_cast<int>(4 * unit(rng));
    const bool mirror = unit(rng) < 0.5;
    const double rad = quarters * std::numbers::pi / 2;
    const double cs = std::cos(rad), sn = std::sin(rad);
    auto fromB = [&](const Pose &p) {
        double dx = p.x - cb.x, dy = p.y - cb.y, deg = p.deg;
        if (mirror) {
            dx = -dx;
            deg = -deg;
        }
        return Pose{ca.x + dx * cs - dy * sn, ca.y + dx * sn + dy * cs, deg + quarters * 90.0};
    };

    // A line through a point near the centre, at a random angle.
    const double angle = 2 * std::numbers::pi * unit(rng);
    const Point normal{std::cos(angle), std::sin(angle)};
    const double offset = (unit(rng) - 0.5) * 0.5 * side;
    auto onASide = [&](const Pose &p) {
        return (p.x - ca.x) * normal.x + (p.y - ca.y) * normal.y < offset;
    };

    Layout child;
    for (const Pose &p : a.poses) {
        if (onASide(p)) child.add(p);
    }
    for (const Pose &p : b.poses) {
        const Pose q = fromB(p);
        if (!onASide(q)) child.add(q);
    }

    // Too many: drop the most overlapped tree, or if none overlaps, the one that shrinks the
    // square most.
    while (static_cast<int>(child.size()) > n) {
        std::size_t worst = 0;
        double most = 0;
        for (std::size_t i = 0; i < child.size(); ++i) {
            const double o = overlapOf(child, i);
            if (o > most) {
                most = o;
                worst = i;
            }
        }
        child.remove(most > 0 ? worst : child.bestTreeToRemove());
    }

    // Too few: add trees at the best spots in a's square.
    const Box square{ca.x - side / 2, ca.y - side / 2, ca.x + side / 2, ca.y + side / 2};
    while (static_cast<int>(child.size()) < n) {
        if (child.size() == 0) {
            child.add({ca.x, ca.y, 360 * unit(rng)});
            continue;
        }
        child.add(bestInsertion(child, square, 4, rng));
    }

    return child;
}
} // namespace packer::packing
