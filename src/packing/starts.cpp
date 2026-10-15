#include "packing/starts.h"
#include "geometry/polygon.h"
#include "tree/tree.h"
#include <cmath>

namespace packer::packing {
Layout bestSingle() {
    Layout layout;
    layout.add({0, 0, 0});

    double bestDeg = 0;
    double bestSide = layout.side();
    // A square looks the same after a 90° turn, so 0° → 90° covers every case.
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

Layout randomStart(int n, std::mt19937_64 &rng) {
    std::uniform_real_distribution<double> unit(0.0, 1.0);

    // Square with about 4x the trees' total area: ~25% full, so trees fit easily.
    double side = 2.0 * std::sqrt(n * geometry::polygonArea(tree::outline()));
    Layout layout;
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

        // Too crowded: enlarge the area and keep going.
        const bool crowded = !placed;
        if (crowded) side *= 1.2;
    }

    return layout;
}

Layout withoutTree(const Layout &layout) {
    Layout result = layout;
    result.remove(layout.bestTreeToRemove());
    return result;
}

std::optional<Layout> withExtraTree(const Layout &layout, std::mt19937_64 &rng) {
    if (layout.size() == 0) return std::nullopt;

    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const Box b = layout.bounds();
    Layout out = layout;
    out.add({0, 0, 0});
    const std::size_t last = out.size() - 1;

    for (int attempt = 0; attempt < 4000; ++attempt) {
        out.setPose(last, {b.minX + unit(rng) * (b.maxX - b.minX),
                           b.minY + unit(rng) * (b.maxY - b.minY), 360.0 * unit(rng)});
        const Box &t = out.boxes[last];
        const bool fitsX = t.minX >= b.minX && t.maxX <= b.maxX;
        const bool fitsY = t.minY >= b.minY && t.maxY <= b.maxY;
        if (!fitsX || !fitsY) continue;

        if (!out.collidesWithOthers(last)) return out;
    }

    return std::nullopt;
}
} // namespace packer::packing
