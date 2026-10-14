#include "packing/starts.h"
#include "geometry/polygon.h"
#include <cmath>

namespace packer::packing {
Layout bestSingle(const Polygon &shape) {
    Layout layout(shape);
    layout.add({0, 0, 0});
    double bestDeg = 0;
    double bestSide = layout.side();
    // A square looks the same after a 90° turn, so 0..90° covers every case.
    for (int k = 1; k < 90 * 100; ++k) {
        const double deg = k / 100.0;
        layout.setPose(0, {0, 0, deg});
        if (layout.side() < bestSide) {
            bestSide = layout.side();
            bestDeg = deg;
        }
    }
    layout.setPose(0, {0, 0, bestDeg});
    return layout;
}

Layout randomStart(const Polygon &shape, int n, std::mt19937_64 &rng) {
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    // Square with about 4x the trees' total area: ~25% full, so trees fit easily.
    double side = 2.0 * std::sqrt(n * geometry::polygonArea(shape));
    Layout layout(shape);
    while (static_cast<int>(layout.size()) < n) {
        bool placed = false;
        for (int attempt = 0; attempt < 1000 && !placed; ++attempt) {
            layout.add({side * unit(rng), side * unit(rng), 360.0 * unit(rng)});
            if (layout.collidesWithOthers(layout.size() - 1)) {
                layout.remove(layout.size() - 1);
            } else {
                placed = true;
            }
        }
        if (!placed) side *= 1.2; // too crowded: enlarge the area and keep going
    }
    return layout;
}
} // namespace packer::packing