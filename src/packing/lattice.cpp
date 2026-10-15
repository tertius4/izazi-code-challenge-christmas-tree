#include "packing/lattice.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "tree/tree.h"

namespace packer::packing {
namespace {
using geometry::Point;

constexpr double kTol = 1e-9; // slack when testing whether a box lies inside a square

// Where B sits relative to A in the same cell: s·u + t·v.
Point offsetB(const LatticeParams &p) {
    return {p.s * p.a + p.t * p.h * p.a, p.t * p.b};
}

// Largest distance from the tree's origin to a vertex (the tip), the same at any rotation.
constexpr double kReach = 0.8;

Box shifted(const Box &b, double dx, double dy) {
    return {b.minX + dx, b.minY + dy, b.maxX + dx, b.maxY + dy};
}

double wrap01(double x) {
    return x - std::floor(x);
}

double wrapHalf(double x) {
    return x - std::floor(x + 0.5);
}

// Scales the cell (a and b together, so B keeps its place in the cell) down to about the
// smallest size that is still valid. A colliding cell is first grown until it is valid.
// Returns false if no valid size was found.
bool tighten(LatticeParams &p) {
    auto scaled = [&](double f) {
        LatticeParams q = p;
        q.a *= f;
        q.b *= f;
        return q;
    };

    double hi = 1.0;
    while (!latticeIsValid(scaled(hi))) {
        hi *= 1.1;
        if (hi > 20) return false; // B sits on top of A: no size works
    }

    double lo = hi * 0.7;
    for (int i = 0; i < 32; ++i) {
        const double mid = (lo + hi) / 2;
        if (latticeIsValid(scaled(mid))) {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    // Back off by a hair, so rounding far from the origin can never turn a touch into an overlap.
    const LatticeParams safer = scaled(hi * (1 + 1e-9));
    p = latticeIsValid(safer) ? safer : scaled(hi);

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

bool latticeIsValid(const LatticeParams &p) {
    // The pattern repeats, so checking A(0,0) against its neighbours covers every pair:
    // A against A(i,j) covers all A-A pairs, A against B(i,j) all A-B pairs, and B-B pairs
    // are A-A pairs turned 180°.
    const tree::Placed treeA = tree::place({0, 0, p.deg});
    const Point d = offsetB(p);
    tree::Placed other;

    // Two trees can only touch if their origins are within 2 reaches; B is a further |d| away.
    const double range = 2 * kReach + std::hypot(d.x, d.y);
    const int maxJ = static_cast<int>(std::ceil(range / p.b)) + 1;

    for (int j = -maxJ; j <= maxJ; ++j) {
        const double rowX = j * p.h * p.a;
        const int minI = static_cast<int>(std::floor((-range - rowX) / p.a)) - 1;
        const int maxI = static_cast<int>(std::ceil((range - rowX) / p.a)) + 1;

        for (int i = minI; i <= maxI; ++i) {
            for (bool second : {false, true}) {
                const bool isTreeA = !second && i == 0 && j == 0;
                if (isTreeA) continue;

                const Pose pose = latticePose(p, i, j, second);
                // Most neighbours are too far away to touch: skip placing them.
                const double dx = pose.x;
                const double dy = pose.y;
                const bool tooFar = dx * dx + dy * dy > 4 * kReach * kReach;
                if (tooFar) continue;

                tree::place(pose, other);
                if (tree::collide(treeA, other)) return false;
            }
        }
    }

    return true;
}

std::optional<LatticeParams> findLattice(long iterations, std::mt19937_64 &rng) {
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    auto jiggle = [&](double step) { return step * (2 * unit(rng) - 1); };

    // A random cell, roomy enough that tighten() only has to shrink it.
    LatticeParams p;
    for (int attempt = 0;; ++attempt) {
        if (attempt == 100) return std::nullopt;

        p.deg = 360 * unit(rng);
        p.a = 2 * kReach * (1 + unit(rng));
        p.b = 2 * kReach * (1 + unit(rng));
        p.h = unit(rng) - 0.5;
        p.s = unit(rng);
        p.t = unit(rng);

        if (tighten(p)) break;
    }

    // Simulated annealing on the shape of the cell. Every candidate is tightened, so the
    // search only ever compares cells that are as small as they can be.
    LatticeParams best = p;
    double current = areaPerTree(p);

    for (long it = 0; it < iterations; ++it) {
        const double progress = static_cast<double>(it) / iterations;
        const double step = 0.05 * std::pow(0.002, progress); // 0.05 down to 0.0001
        const double temp = 0.01 * std::pow(1e-4, progress);

        LatticeParams q = p;
        const int kind = static_cast<int>(6 * unit(rng));
        const bool all = kind == 5;
        const bool turn = kind == 0 || all;
        const bool stretch = kind == 1 || all; // one way, squashing the other
        const bool shiftRows = kind == 2 || all;
        const bool moveBAlongU = kind == 3 || all;
        const bool moveBAlongV = kind == 4 || all;

        if (turn) q.deg += 30 * jiggle(step);
        if (stretch) {
            const double r = std::exp(jiggle(step));
            q.a *= r;
            q.b /= r;
        }
        if (shiftRows) q.h = wrapHalf(q.h + jiggle(step));
        if (moveBAlongU) q.s = wrap01(q.s + jiggle(step));
        if (moveBAlongV) q.t = wrap01(q.t + jiggle(step));
        if (!tighten(q)) continue;

        const double next = areaPerTree(q);
        const double worseBy = next - current;
        const bool smaller = worseBy < 0;
        const bool accept = smaller || unit(rng) < std::exp(-worseBy / temp);
        if (accept) {
            p = q;
            current = next;
            if (next < areaPerTree(best)) best = q;
        }
    }

    best.deg = geometry::normalizeDeg(best.deg);

    return best;
}

Layout cutLattice(const LatticeParams &p, int n) {
    const Box boxA = tree::place({0, 0, p.deg}).box;
    const Box boxB = tree::place({0, 0, p.deg + 180}).box;
    const Point d = offsetB(p);
    const double margin = 2 * kReach + std::hypot(d.x, d.y);

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
        const bool fitsX = b.minX >= x0 - kTol && b.maxX <= x0 + side + kTol;
        const bool fitsY = b.minY >= y0 - kTol && b.maxY <= y0 + side + kTol;

        return fitsX && fitsY;
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
                const bool fitsY = t.box.minY >= y0 - kTol && t.box.maxY <= y0 + side + kTol;
                if (fitsY) {
                    xs.push_back(t.box.minX - p.a * std::floor(t.box.minX / p.a));
                }
            }

            std::sort(xs.begin(), xs.end());
            auto nearlyEqual = [](double l, double r) { return r - l < kTol; };
            xs.erase(std::unique(xs.begin(), xs.end(), nearlyEqual), xs.end());

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
    Layout layout;
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
