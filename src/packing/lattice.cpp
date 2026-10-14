#include "packing/lattice.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "geometry/collision.h"

namespace packer::packing {
namespace {
using geometry::Point;

constexpr double kMinStep = 0.05; // far below any valid step; keeps the neighbour loops bounded
constexpr double kTol = 1e-9;     // slack when testing whether a box lies inside a square

// Where B sits relative to A in the same cell: s·u + t·v.
Point offsetB(const LatticeParams &p) {
    return {p.s * p.a + p.t * p.h * p.a, p.t * p.b};
}

// Largest distance from the tree's origin to a vertex. The same at any rotation.
double reachOf(const Polygon &shape) {
    double reach = 0;
    for (const Point &v : shape) {
        reach = std::max(reach, std::hypot(v.x, v.y));
    }
    return reach;
}

Box shifted(const Box &b, double dx, double dy) {
    return {b.minX + dx, b.minY + dy, b.maxX + dx, b.maxY + dy};
}

// True if `fixed` collides with `other` moved by (dx, dy). Builds the moved copy only when the
// boxes overlap, which is rare.
bool collidesShifted(const Polygon &fixed, const Box &fixedBox, const Polygon &other,
                     const Box &otherBox, double dx, double dy) {
    const Box movedBox = shifted(otherBox, dx, dy);
    if (!geometry::boxesOverlap(fixedBox, movedBox)) return false;
    Polygon moved = other;
    for (Point &v : moved) {
        v.x += dx;
        v.y += dy;
    }
    return geometry::polygonsCollide(fixed, fixedBox, moved, movedBox);
}

// The smallest value in [lo, hi] for which isValid holds. `hi` must be valid. Interlocking
// trees can collide at some size and fit again at a smaller one, so a plain bisection would
// stop at the first collision. Instead: scan downward on a geometric grid for the smallest
// valid point, then bisect between it and the next (invalid) grid point below.
template <class IsValid> double smallestValid(double lo, double hi, IsValid isValid) {
    constexpr int kScan = 120;
    const double ratio = std::pow(lo / hi, 1.0 / kScan);
    double best = hi;
    double below = hi * ratio;
    double v = hi;
    for (int k = 1; k <= kScan; ++k) {
        v *= ratio;
        if (isValid(v)) {
            best = v;
            below = v * ratio;
        }
    }
    for (int i = 0; i < 40; ++i) {
        const double mid = (below + best) / 2;
        if (isValid(mid)) {
            best = mid;
        } else {
            below = mid;
        }
    }
    return best;
}

double wrap01(double x) {
    return x - std::floor(x);
}

double wrapHalf(double x) {
    return x - std::floor(x + 0.5);
}

// Makes the cell as small as possible for p's deg, h, s and t: shrink a, then b, and repeat
// once, because a smaller b can free up room for a smaller a. Returns false if even a roomy
// cell collides (B sits on top of A).
bool shrink(const Polygon &shape, LatticeParams &p) {
    const double roomy = 4.0 * reachOf(shape);
    p.a = roomy;
    p.b = roomy;
    if (!latticeIsValid(shape, p)) return false;
    for (int round = 0; round < 2; ++round) {
        p.a = smallestValid(kMinStep, p.a, [&](double v) {
            LatticeParams q = p;
            q.a = v;
            return latticeIsValid(shape, q);
        });
        p.b = smallestValid(kMinStep, p.b, [&](double v) {
            LatticeParams q = p;
            q.b = v;
            return latticeIsValid(shape, q);
        });
    }
    // Back off by a hair, so rounding far from the origin can never turn a touch into an overlap.
    LatticeParams safer = p;
    safer.a *= 1 + 1e-9;
    safer.b *= 1 + 1e-9;
    if (latticeIsValid(shape, safer)) p = safer;
    return true;
}
} // namespace

double areaPerTree(const LatticeParams &p) {
    return p.a * p.b / 2;
}

Pose latticePose(const LatticeParams &p, int i, int j, bool second) {
    const double x = i * p.a + j * p.h * p.a;
    const double y = j * p.b;
    if (!second) return {x, y, p.deg};
    const Point d = offsetB(p);
    return {x + d.x, y + d.y, p.deg + 180};
}

bool latticeIsValid(const Polygon &shape, const LatticeParams &p) {
    // The pattern repeats, so checking A(0,0) against its neighbours covers every pair:
    // A against A(i,j) covers all A-A pairs, A against B(i,j) all A-B pairs, and B-B pairs
    // are A-A pairs turned 180°.
    const Polygon treeA = geometry::transform(shape, {0, 0, p.deg});
    const Polygon treeB = geometry::transform(shape, {0, 0, p.deg + 180});
    const Box boxA = geometry::boundingBox(treeA);
    const Box boxB = geometry::boundingBox(treeB);
    const Point d = offsetB(p);

    // Two trees can only touch if their origins are within 2 reaches; B is a further |d| away.
    const double range = 2 * reachOf(shape) + std::hypot(d.x, d.y);
    const int maxJ = static_cast<int>(std::ceil(range / p.b)) + 1;
    for (int j = -maxJ; j <= maxJ; ++j) {
        const double rowX = j * p.h * p.a;
        const int minI = static_cast<int>(std::floor((-range - rowX) / p.a)) - 1;
        const int maxI = static_cast<int>(std::ceil((range - rowX) / p.a)) + 1;
        for (int i = minI; i <= maxI; ++i) {
            const Pose other = latticePose(p, i, j, false);
            if ((i != 0 || j != 0) &&
                collidesShifted(treeA, boxA, treeA, boxA, other.x, other.y)) {
                return false;
            }
            const Pose otherB = latticePose(p, i, j, true);
            if (collidesShifted(treeA, boxA, treeB, boxB, otherB.x, otherB.y)) return false;
        }
    }
    return true;
}

std::optional<LatticeParams> findLattice(const Polygon &shape, double startDeg, int starts,
                                         int tweaks, std::mt19937_64 &rng) {
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    auto jiggle = [&](double step) {
        return step * (2 * unit(rng) - 1);
    };

    std::optional<LatticeParams> best;
    for (int k = 0; k < starts; ++k) {
        LatticeParams p;
        p.deg = startDeg;
        p.h = unit(rng) - 0.5;
        p.s = unit(rng);
        p.t = unit(rng);
        if (!shrink(shape, p)) continue;

        // Local search: nudge the shape of the cell, re-shrink, keep it if the area went down.
        for (int it = 0; it < tweaks; ++it) {
            // Big nudges first, fine ones at the end: 0.1 down to 0.001.
            const double step = 0.1 * std::pow(0.01, static_cast<double>(it) / tweaks);
            LatticeParams q = p;
            q.deg += 30 * jiggle(step); // ±3° down to ±0.03°
            q.h = wrapHalf(q.h + jiggle(step));
            q.s = wrap01(q.s + jiggle(step));
            q.t = wrap01(q.t + jiggle(step));
            if (shrink(shape, q) && areaPerTree(q) < areaPerTree(p)) p = q;
        }

        if (!best || areaPerTree(p) < areaPerTree(*best)) best = p;
    }
    return best;
}

Layout cutLattice(const Polygon &shape, const LatticeParams &p, int n) {
    const Box boxA = geometry::boundingBox(geometry::transform(shape, {0, 0, p.deg}));
    const Box boxB = geometry::boundingBox(geometry::transform(shape, {0, 0, p.deg + 180}));
    const Point d = offsetB(p);
    const double margin = 2 * reachOf(shape) + std::hypot(d.x, d.y);

    struct Tree {
        Pose pose;
        Box box;
    };
    std::vector<Tree> trees;

    // Every tree that could lie inside a square of this side whose bottom edge is near row 0
    // and whose left edge is in [0, a).
    auto generate = [&](double side) {
        trees.clear();
        const int minJ = static_cast<int>(std::floor(-margin / p.b)) - 1;
        const int maxJ = static_cast<int>(std::ceil((side + margin) / p.b)) + 1;
        for (int j = minJ; j <= maxJ; ++j) {
            const double rowX = j * p.h * p.a;
            const int minI = static_cast<int>(std::floor((-margin - rowX) / p.a)) - 1;
            const int maxI = static_cast<int>(std::ceil((p.a + side + margin - rowX) / p.a)) + 1;
            for (int i = minI; i <= maxI; ++i) {
                for (bool second : {false, true}) {
                    const Pose pose = latticePose(p, i, j, second);
                    trees.push_back({pose, shifted(second ? boxB : boxA, pose.x, pose.y)});
                }
            }
        }
    };

    auto inside = [](const Box &b, double x0, double y0, double side) {
        return b.minX >= x0 - kTol && b.maxX <= x0 + side + kTol && b.minY >= y0 - kTol &&
               b.maxY <= y0 + side + kTol;
    };

    struct Window {
        double x0 = 0;
        double y0 = 0;
        int count = -1;
    };

    // Where to put a square of this side so it holds the most whole trees. The best square can
    // always slide down until its bottom edge touches a tree, and left until its left edge does.
    // Moving by whole cells changes nothing, so the bottom edge only needs trying at the lowest
    // points of A(0,0) and B(0,0), and the left edge at the leftmost points of the trees in the
    // band, reduced into [0, a).
    auto bestWindow = [&](double side) {
        Window best;
        for (double y0 : {boxA.minY, boxB.minY + d.y}) {
            std::vector<double> xs;
            for (const Tree &t : trees) {
                if (t.box.minY >= y0 - kTol && t.box.maxY <= y0 + side + kTol) {
                    xs.push_back(t.box.minX - p.a * std::floor(t.box.minX / p.a));
                }
            }
            std::sort(xs.begin(), xs.end());
            xs.erase(std::unique(xs.begin(), xs.end(),
                                 [](double l, double r) { return r - l < kTol; }),
                     xs.end());
            for (double x0 : xs) {
                int count = 0;
                for (const Tree &t : trees) {
                    if (inside(t.box, x0, y0, side)) ++count;
                }
                if (count > best.count) best = {x0, y0, count};
            }
        }
        return best;
    };

    // A side that surely holds n trees, then bisect down to the smallest one that does.
    double hi = std::sqrt(n * areaPerTree(p)) * 1.5 + 2.0;
    generate(hi);
    while (bestWindow(hi).count < n) {
        hi *= 1.5;
        generate(hi);
    }
    double lo = 0;
    for (int it = 0; it < 45; ++it) {
        const double mid = (lo + hi) / 2;
        if (bestWindow(mid).count >= n) {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    const Window w = bestWindow(hi);
    Layout layout(shape);
    for (const Tree &t : trees) {
        if (inside(t.box, w.x0, w.y0, hi)) layout.add(t.pose);
    }
    // The square may hold a few more than n: drop the ones that shrink it most.
    while (static_cast<int>(layout.size()) > n) {
        layout.remove(layout.bestTreeToRemove());
    }
    return layout;
}
} // namespace packer::packing
