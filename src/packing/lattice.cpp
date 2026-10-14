#include "packing/lattice.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace packer::packing {
namespace {
Pose latticePose(const LatticeParams &p, int col, int row) {
    const bool flipped = col % 2 == 1;
    return {col * p.colStep, row * p.rowStep + (flipped ? p.shift : 0.0),
            p.baseDeg + (flipped ? 180.0 : 0.0)};
}

// A 6x4 block covers every kind of neighbour, so if it is clean, the whole lattice is.
bool patchIsValid(const Polygon &shape, const LatticeParams &p) {
    Layout patch(shape);
    for (int col = 0; col < 6; ++col) {
        for (int row = 0; row < 4; ++row) {
            patch.add(latticePose(p, col, row));
        }
    }
    return patch.isValid();
}

// Shrinks `hi` toward `lo` while isValid(value) holds. `hi` must be valid. Returns a valid value.
template <class IsValid> double bisect(double lo, double hi, IsValid isValid) {
    for (int i = 0; i < 50; ++i) {
        const double mid = (lo + hi) / 2;
        if (isValid(mid)) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    return hi;
}
} // namespace

LatticeParams findLattice(const Polygon &shape, double baseDeg) {
    // Steps of 2x the tree's reach from its origin can never collide: a safe upper bound.
    double reach = 0;
    for (const auto &v : shape)
        reach = std::max(reach, std::hypot(v.x, v.y));
    const double safe = 2.0 * reach * 1.01;

    LatticeParams best;
    double bestArea = std::numeric_limits<double>::infinity();
    const int kShifts = 40;
    for (int k = 0; k < kShifts; ++k) {
        LatticeParams p;
        p.baseDeg = baseDeg;
        p.shift = safe * k / kShifts;
        p.rowStep = safe;
        p.colStep = bisect(0.0, safe, [&](double v) {
            LatticeParams q = p;
            q.colStep = v;
            return patchIsValid(shape, q);
        });
        p.rowStep = bisect(0.0, safe, [&](double v) {
            LatticeParams q = p;
            q.rowStep = v;
            return patchIsValid(shape, q);
        });
        const double area = p.colStep * p.rowStep; // area per 2 trees (one of each orientation)
        if (area < bestArea) {
            bestArea = area;
            best = p;
        }
    }
    return best;
}

Layout buildLattice(const Polygon &shape, const LatticeParams &params, int n) {
    Layout best(shape);
    double bestSide = std::numeric_limits<double>::infinity();
    for (int cols = 1; cols <= n; ++cols) {
        const int rows = (n + cols - 1) / cols; // ceil(n / cols)
        if ((cols - 1) * rows >= n) continue;   // fewer columns would already hold n trees
        Layout layout(shape);
        for (int k = 0; k < n; ++k)
            layout.add(latticePose(params, k / rows, k % rows));
        if (layout.side() < bestSide) {
            bestSide = layout.side();
            best = std::move(layout);
        }
    }
    return best;
}
} // namespace packer::packing
