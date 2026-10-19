#include "packing/starts.h"

#include <cmath>

#include "geometry/polygon.h"
#include "tree/tree.h"

namespace packer::packing {
Layout bestSingle() {
    Layout layout;
    layout.add({0, 0, 0});

    double bestDeg = 0;
    double bestSide = layout.side();
    // A square looks the same after a 90° turn, so 0° → 90° covers every case.
    for (int hundredths = 1; hundredths < 90 * 100; ++hundredths) {
        const double deg = hundredths / 100.0;
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

Layout denseStart(int n, double side, std::mt19937_64 &rng) {
    std::uniform_real_distribution<double> unit(-0.5, 0.5);
    Layout layout;
    for (int i = 0; i < n; ++i) {
        layout.add({side * unit(rng), side * unit(rng), 360.0 * (unit(rng) + 0.5)});
    }

    return layout;
}

Layout symmetricStart(int n, double side, std::mt19937_64 &rng) {
    std::uniform_real_distribution<double> unit(-0.5, 0.5);
    const int pairs = n / 2;
    std::vector<Pose> half;
    for (int i = 0; i < pairs; ++i) {
        half.push_back({side * unit(rng), side * unit(rng), 360.0 * (unit(rng) + 0.5)});
    }

    Layout layout;
    for (const Pose &pose : half) {
        layout.add(pose);
    }
    for (const Pose &pose : half) {
        layout.add({-pose.x, -pose.y, geometry::normalizeDeg(pose.deg + 180)});
    }
    if (n % 2 == 1) {
        layout.add({0.1 * side * unit(rng), 0.1 * side * unit(rng), 360.0 * (unit(rng) + 0.5)});
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
    const Box bounds = layout.bounds();
    Layout grown = layout;
    grown.add({0, 0, 0});
    const std::size_t last = grown.size() - 1;

    for (int attempt = 0; attempt < 4000; ++attempt) {
        grown.setPose(last,
                      {bounds.minX + unit(rng) * (bounds.maxX - bounds.minX),
                       bounds.minY + unit(rng) * (bounds.maxY - bounds.minY), 360.0 * unit(rng)});
        const Box &box = grown.boxes[last];
        const bool fitsX = box.minX >= bounds.minX && box.maxX <= bounds.maxX;
        const bool fitsY = box.minY >= bounds.minY && box.maxY <= bounds.maxY;
        if (!fitsX || !fitsY) continue;

        if (!grown.collidesWithOthers(last)) return grown;
    }

    return std::nullopt;
}
} // namespace packer::packing
