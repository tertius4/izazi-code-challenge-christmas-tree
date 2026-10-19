#include "packing/crystal.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numbers>

#include "tree/tree.h"
#include "util/parallel.h"

namespace packer::packing {
namespace {
using geometry::kInf;
using geometry::Point;

constexpr double kTol = 1e-9;  // slack when testing whether a box lies inside a square
constexpr double kReach = 0.8; // largest distance from a tree's origin to its outline (the tip)

double cross(const Point &a, const Point &b) {
    return a.x * b.y - a.y * b.x;
}

// The point at these fractions of u and v.
Point at(const Crystal &crystal, double fracU, double fracV) {
    return {fracU * crystal.u.x + fracV * crystal.v.x, fracU * crystal.u.y + fracV * crystal.v.y};
}

// The point's coordinates as fractions of u and v.
Point fractions(const Crystal &crystal, const Point &point) {
    const double cellArea = cross(crystal.u, crystal.v); // signed
    return {cross(point, crystal.v) / cellArea, cross(crystal.u, point) / cellArea};
}

Box shifted(const Box &box, double dx, double dy) {
    return {box.minX + dx, box.minY + dy, box.maxX + dx, box.maxY + dy};
}

double wrap01(double x) {
    return x - std::floor(x);
}

Crystal scaled(const Crystal &crystal, double factor) {
    Crystal out = crystal;
    out.u = {crystal.u.x * factor, crystal.u.y * factor};
    out.v = {crystal.v.x * factor, crystal.v.y * factor};
    return out;
}

// Scales the cell (the trees keep their places in it) to about the smallest size that is still
// valid, growing it first if it collides. Returns false if no valid size was found.
bool tighten(Crystal &crystal) {
    double hi = 1.0;
    double lo = 0;
    if (crystalIsValid(crystal)) {
        // Shrink in steps until it collides, so the bisection below starts with a tight bracket.
        lo = hi;
        for (int i = 0; i < 60; ++i) {
            lo *= 0.97;
            if (!crystalIsValid(scaled(crystal, lo))) break;
            hi = lo;
        }
    } else {
        for (int i = 0;; ++i) {
            if (i == 100) return false; // two trees sit on top of each other: no size works
            lo = hi;
            hi *= 1.03;
            if (crystalIsValid(scaled(crystal, hi))) break;
        }
    }

    for (int i = 0; i < 28; ++i) {
        const double mid = (lo + hi) / 2;
        if (crystalIsValid(scaled(crystal, mid))) {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    // Back off by a hair, so rounding far from the origin can never turn a touch into an overlap.
    const Crystal safer = scaled(crystal, hi * (1 + 1e-9));
    crystal = crystalIsValid(safer) ? safer : scaled(crystal, hi);

    return true;
}

// The cell of the search: u = (a, 0) and v = (h·a, b). Turning the whole pattern changes
// nothing about its density, so u stays level.
struct Shape {
    double a = 0;
    double b = 0;
    double h = 0;
};

Crystal build(const Shape &shape, const std::vector<Pose> &cell) {
    return {{shape.a, 0}, {shape.h * shape.a, shape.b}, cell};
}

Shape shapeOf(const Crystal &crystal) {
    return {crystal.u.x, crystal.v.y, crystal.v.x / crystal.u.x};
}
} // namespace

double areaPerTree(const Crystal &crystal) {
    return std::abs(cross(crystal.u, crystal.v)) / static_cast<double>(crystal.cell.size());
}

Pose crystalPose(const Crystal &crystal, int i, int j, std::size_t index) {
    const Pose &tree = crystal.cell[index];
    const Point point = at(crystal, i + tree.x, j + tree.y);
    return {point.x, point.y, tree.deg};
}

bool crystalIsValid(const Crystal &crystal) {
    const double cellArea = cross(crystal.u, crystal.v); // signed
    if (std::abs(cellArea) < 1e-12) return false;

    // A vector w has fractions (w×v, u×w) / cellArea, so |fraction along u| <= |w|·|v| / |area|
    // and |fraction along v| <= |w|·|u| / |area|. Trees only touch if their origins are within
    // 2 reaches: that bounds the cells to check.
    const double range = 2 * kReach;
    const double spanI = range * std::hypot(crystal.v.x, crystal.v.y) / std::abs(cellArea);
    const double spanJ = range * std::hypot(crystal.u.x, crystal.u.y) / std::abs(cellArea);

    const std::size_t perCell = crystal.cell.size();
    std::vector<tree::Placed> base(perCell);
    for (std::size_t treeA = 0; treeA < perCell; ++treeA) {
        base[treeA] = tree::place(crystalPose(crystal, 0, 0, treeA));
    }

    tree::Placed other;
    for (std::size_t treeA = 0; treeA < perCell; ++treeA) {
        const Pose poseA = crystalPose(crystal, 0, 0, treeA);
        // Pairs (A, B) and (B, A) are the same pairs of trees, moved by a whole cell.
        for (std::size_t treeB = treeA; treeB < perCell; ++treeB) {
            const double fracU = crystal.cell[treeB].x - crystal.cell[treeA].x;
            const double fracV = crystal.cell[treeB].y - crystal.cell[treeA].y;
            const int minI = static_cast<int>(std::ceil(-fracU - spanI));
            const int maxI = static_cast<int>(std::floor(-fracU + spanI));
            const int minJ = static_cast<int>(std::ceil(-fracV - spanJ));
            const int maxJ = static_cast<int>(std::floor(-fracV + spanJ));

            for (int i = minI; i <= maxI; ++i) {
                for (int j = minJ; j <= maxJ; ++j) {
                    if (treeA == treeB && i == 0 && j == 0) continue;

                    const Pose poseB = crystalPose(crystal, i, j, treeB);
                    const double dx = poseB.x - poseA.x;
                    const double dy = poseB.y - poseA.y;
                    if (dx * dx + dy * dy > range * range) continue;

                    tree::place(poseB, other);
                    if (tree::collide(base[treeA], other)) return false;
                }
            }
        }
    }

    return true;
}

Crystal turned(const Crystal &crystal, double deg) {
    const double rad = deg * std::numbers::pi / 180.0;
    const double cosine = std::cos(rad);
    const double sine = std::sin(rad);
    auto turn = [&](const Point &point) {
        return Point{point.x * cosine - point.y * sine, point.x * sine + point.y * cosine};
    };

    Crystal out = crystal;
    out.u = turn(crystal.u);
    out.v = turn(crystal.v);
    for (Pose &tree : out.cell) {
        tree.deg = geometry::normalizeDeg(tree.deg + deg); // fractions turn with u and v
    }

    return out;
}

Crystal doubled(const Crystal &crystal) {
    Crystal out = crystal;
    out.u = {2 * crystal.u.x, 2 * crystal.u.y};
    out.cell.clear();
    for (const Pose &tree : crystal.cell) {
        out.cell.push_back({tree.x / 2, tree.y, tree.deg});
    }
    for (const Pose &tree : crystal.cell) {
        out.cell.push_back({tree.x / 2 + 0.5, tree.y, tree.deg});
    }
    return out;
}

Crystal toCrystal(const LatticeParams &lattice) {
    Crystal crystal;
    crystal.u = {lattice.a, 0};
    crystal.v = {lattice.h * lattice.a, lattice.b};
    crystal.cell = {{0, 0, lattice.deg}, {lattice.s, lattice.t, lattice.deg + 180}};
    return crystal;
}

std::optional<Crystal> findCrystal(int perCell, bool paired, long iterations, std::mt19937_64 &rng,
                                   const std::atomic<bool> *stop, const Crystal *start) {
    if (start) perCell = static_cast<int>(start->cell.size());
    if (perCell < 1) return std::nullopt;

    std::uniform_real_distribution<double> unit(0.0, 1.0);
    auto jiggle = [&](double step) { return step * (2 * unit(rng) - 1); };
    const auto count = static_cast<std::size_t>(perCell);
    auto partnerOf = [&](std::size_t index) -> std::size_t { // the tree paired with this one
        if (!paired) return index;
        const std::size_t partner = index ^ 1u;
        return partner < count ? partner : index;
    };

    // A random roomy cell, so that tighten() only has to shrink it.
    Crystal current;
    for (int attempt = 0; start == nullptr; ++attempt) {
        if (attempt == 100) return std::nullopt;

        const double side = std::sqrt(static_cast<double>(perCell)) * (1 + unit(rng));
        const Shape shape{side, side * (0.7 + 0.6 * unit(rng)), unit(rng) - 0.5};
        std::vector<Pose> cell(count);
        for (std::size_t k = 0; k < count; ++k) {
            cell[k] = {unit(rng), unit(rng), 360 * unit(rng)};
            if (paired && k % 2 == 1) cell[k].deg = cell[k - 1].deg + 180;
        }
        cell[0].x = cell[0].y = 0; // moving every tree together changes nothing

        current = build(shape, cell);
        if (tighten(current)) break;
    }
    if (start) {
        // Level u, as the moves below expect: turning the pattern changes nothing.
        current = turned(*start, -std::atan2(start->u.y, start->u.x) * 180 / std::numbers::pi);
        if (!tighten(current)) return std::nullopt;
    }

    Crystal best = current;
    double currentArea = areaPerTree(current);

    for (long it = 0; it < iterations; ++it) {
        if (stop && it % 256 == 0 && stop->load(std::memory_order_relaxed)) break;

        const double progress = static_cast<double>(it) / iterations;
        // A given start is already good: search near it, cooler, rather than melt it.
        const double heat = start ? 0.1 : 1.0;
        const double step = heat * 0.05 * std::pow(0.002, progress); // 0.05 down to 0.0001
        const double temp = heat * heat * 0.01 * std::pow(1e-4, progress);

        Shape shape = shapeOf(current);
        std::vector<Pose> cell = current.cell;
        const std::size_t picked = std::uniform_int_distribution<std::size_t>(0, count - 1)(rng);
        const double kind = unit(rng);

        if (kind < 0.2) { // stretch one way, squash the other
            const double ratio = std::exp(jiggle(step));
            shape.a *= ratio;
            shape.b /= ratio;
        } else if (kind < 0.35) {
            shape.h += jiggle(step);
        } else if (kind < 0.65 && count > 1) {
            const std::size_t moved = 1 + picked % (count - 1); // tree 0 stays at the corner
            cell[moved].x = wrap01(cell[moved].x + jiggle(step));
            cell[moved].y = wrap01(cell[moved].y + jiggle(step));
        } else if (kind < 0.9) {
            cell[picked].deg += 30 * jiggle(step);
            const std::size_t partner = partnerOf(picked);
            if (partner != picked) cell[partner].deg = cell[picked].deg + 180;
        } else if (kind < 0.93) {
            cell[picked].deg += 180;
            const std::size_t partner = partnerOf(picked);
            if (partner != picked) cell[partner].deg += 180;
        } else { // turn every tree together, against the cell
            const double turn = 30 * jiggle(step);
            for (Pose &tree : cell)
                tree.deg += turn;
        }

        Crystal candidate = build(shape, cell);
        if (!tighten(candidate)) continue;

        const double area = areaPerTree(candidate);
        const double worseBy = area - currentArea;
        const bool accept = worseBy < 0 || unit(rng) < std::exp(-worseBy / temp);
        if (accept) {
            current = std::move(candidate);
            currentArea = area;
            if (area < areaPerTree(best)) best = current;
        }
    }

    for (Pose &tree : best.cell) {
        tree.deg = geometry::normalizeDeg(tree.deg);
    }

    return best;
}

std::vector<Crystal> findCrystals(int searches, long iterations, int keep, int threads,
                                  unsigned seed, const std::atomic<bool> *stop) {
    // The kinds of cell to try, in turn: the double lattice, then freer and bigger cells.
    struct Kind {
        int perCell;
        bool paired;
    };
    const Kind kinds[] = {{2, true}, {4, true}, {2, false}, {3, false}, {4, false}};

    std::vector<std::optional<Crystal>> found(static_cast<std::size_t>(2 * searches));

    // Runs search(index) for every index in [first, last) on the threads.
    auto runAll = [&](int first, int last, const auto &search) {
        const auto count = static_cast<std::size_t>(last - first);
        util::parallelFor(count, std::max(1, threads), [&](std::size_t k) {
            if (stop && stop->load()) return;
            const int index = first + static_cast<int>(k);
            std::mt19937_64 rng(seed * 1000003ULL + static_cast<unsigned>(index));
            found[static_cast<std::size_t>(index)] = search(index, rng);
        });
    };

    // First from random cells. Bigger cells have more to arrange: give them longer.
    runAll(0, searches, [&](int index, std::mt19937_64 &rng) {
        const Kind &kind = kinds[index % std::size(kinds)];
        return findCrystal(kind.perCell, kind.paired, iterations * kind.perCell / 2, rng, stop);
    });

    // Then from the densest found so far, doubled: a cell of twice the trees is at least as
    // dense, and its copies are free to drift apart into something denser.
    std::vector<Crystal> seeds;
    for (int index = 0; index < searches; ++index) {
        const std::optional<Crystal> &crystal = found[static_cast<std::size_t>(index)];
        if (crystal) seeds.push_back(*crystal);
    }
    std::sort(seeds.begin(), seeds.end(),
              [](const Crystal &a, const Crystal &b) { return areaPerTree(a) < areaPerTree(b); });
    if (!seeds.empty()) {
        const std::size_t densest = std::min<std::size_t>(seeds.size(), 3);
        runAll(searches, 2 * searches, [&](int index, std::mt19937_64 &rng) {
            const Crystal start = doubled(seeds[static_cast<std::size_t>(index) % densest]);
            return findCrystal(0, false, iterations, rng, stop, &start);
        });
    }

    std::vector<Crystal> sorted;
    for (const std::optional<Crystal> &crystal : found) {
        if (crystal) sorted.push_back(*crystal);
    }
    std::sort(sorted.begin(), sorted.end(),
              [](const Crystal &a, const Crystal &b) { return areaPerTree(a) < areaPerTree(b); });

    // Searches that end in about the same pattern end at about the same density: keep one of
    // each, so the kept ones cut differently.
    std::vector<Crystal> kept;
    for (const Crystal &crystal : sorted) {
        if (static_cast<int>(kept.size()) == keep) break;
        const bool isNew =
            kept.empty() || areaPerTree(crystal) > areaPerTree(kept.back()) * (1 + 1e-3);
        if (isNew) kept.push_back(crystal);
    }

    for (const Crystal &crystal : kept) {
        std::printf("crystal with %zu trees per cell: area per tree %.6f\n", crystal.cell.size(),
                    areaPerTree(crystal));
    }
    std::fflush(stdout);

    return kept;
}

Layout cutCrystal(const Crystal &crystal, int n) {
    const std::size_t perCell = crystal.cell.size();

    // A tree of the pattern that a square might hold.
    struct Candidate {
        Pose pose;
        Box box;
    };

    std::vector<Box> baseBox(perCell);
    std::vector<Candidate> anchors(perCell); // the trees of cell (0, 0)
    for (std::size_t k = 0; k < perCell; ++k) {
        baseBox[k] = tree::place({0, 0, crystal.cell[k].deg}).box;
        const Pose pose = crystalPose(crystal, 0, 0, k);
        anchors[k] = {pose, shifted(baseBox[k], pose.x, pose.y)};
    }

    Box anchorSpan = anchors[0].box;
    for (const Candidate &anchor : anchors) {
        anchorSpan.minX = std::min(anchorSpan.minX, anchor.pose.x);
        anchorSpan.minY = std::min(anchorSpan.minY, anchor.pose.y);
        anchorSpan.maxX = std::max(anchorSpan.maxX, anchor.pose.x);
        anchorSpan.maxY = std::max(anchorSpan.maxY, anchor.pose.y);
    }

    std::vector<Candidate> candidates;

    // A square that fits whole trees can slide down until its bottom edge touches a tree, then
    // move by whole cells until that tree is one of cell (0, 0). So the square only needs
    // trying with its bottom on an anchor, holding that anchor: every tree such a square of
    // this side could hold lies within this region.
    auto generate = [&](double side) {
        candidates.clear();
        const Box region{anchorSpan.minX - side - 2 * kReach, anchorSpan.minY - 2 * kReach,
                         anchorSpan.maxX + side + 2 * kReach, anchorSpan.maxY + side + 2 * kReach};

        double minI = kInf, maxI = -kInf, minJ = kInf, maxJ = -kInf;
        const Point corners[] = {{region.minX, region.minY},
                                 {region.maxX, region.minY},
                                 {region.minX, region.maxY},
                                 {region.maxX, region.maxY}};
        for (const Point &corner : corners) {
            const Point frac = fractions(crystal, corner);
            minI = std::min(minI, frac.x);
            maxI = std::max(maxI, frac.x);
            minJ = std::min(minJ, frac.y);
            maxJ = std::max(maxJ, frac.y);
        }

        const int firstI = static_cast<int>(std::floor(minI)) - 1;
        const int lastI = static_cast<int>(std::ceil(maxI)) + 1;
        const int firstJ = static_cast<int>(std::floor(minJ)) - 1;
        const int lastJ = static_cast<int>(std::ceil(maxJ)) + 1;
        for (int i = firstI; i <= lastI; ++i) {
            for (int j = firstJ; j <= lastJ; ++j) {
                for (std::size_t k = 0; k < perCell; ++k) {
                    const Pose pose = crystalPose(crystal, i, j, k);
                    const bool inRegion = pose.x >= region.minX && pose.x <= region.maxX &&
                                          pose.y >= region.minY && pose.y <= region.maxY;
                    if (inRegion) candidates.push_back({pose, shifted(baseBox[k], pose.x, pose.y)});
                }
            }
        }
    };

    struct Window {
        double x0 = 0;
        double y0 = 0;
        int count = -1;
    };

    // The square of this side, its bottom on an anchor and holding that anchor, that holds the
    // most whole trees. Its left edge can slide right until it touches a tree.
    std::vector<std::pair<double, double>> band; // (minX, maxX) of the trees within the rows
    auto bestWindow = [&](double side) {
        Window best;
        for (const Candidate &anchor : anchors) {
            const double y0 = anchor.box.minY;
            band.clear();
            for (const Candidate &candidate : candidates) {
                const Box &box = candidate.box;
                const bool fitsY = box.minY >= y0 - kTol && box.maxY <= y0 + side + kTol;
                if (fitsY) band.push_back({box.minX, box.maxX});
            }
            std::sort(band.begin(), band.end());

            const double firstX = anchor.box.maxX - side - kTol;
            const double lastX = anchor.box.minX + kTol;
            const auto from = std::lower_bound(band.begin(), band.end(), std::pair{firstX, -kInf});
            for (auto left = from; left != band.end() && left->first <= lastX; ++left) {
                const double x0 = left->first;
                int count = 0;
                const auto firstIn =
                    std::lower_bound(band.begin(), band.end(), std::pair{x0 - kTol, -kInf});
                for (auto tree = firstIn; tree != band.end() && tree->first <= x0 + side + kTol;
                     ++tree) {
                    if (tree->second <= x0 + side + kTol) ++count;
                }
                if (count > best.count) best = {x0, y0, count};
            }
        }
        return best;
    };

    // A side that surely holds n trees, then bisect down to the smallest one that does.
    double hi = std::sqrt(n * areaPerTree(crystal)) * 1.3 + 1.5;
    generate(hi);
    while (bestWindow(hi).count < n) {
        hi *= 1.3;
        generate(hi);
    }

    double lo = 0;
    for (int it = 0; it < 40; ++it) {
        const double mid = (lo + hi) / 2;
        if (bestWindow(mid).count >= n) {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    const Window window = bestWindow(hi);
    Layout layout;
    for (const Candidate &candidate : candidates) {
        const Box &box = candidate.box;
        const bool fitsX = box.minX >= window.x0 - kTol && box.maxX <= window.x0 + hi + kTol;
        const bool fitsY = box.minY >= window.y0 - kTol && box.maxY <= window.y0 + hi + kTol;
        if (fitsX && fitsY) layout.add(candidate.pose);
    }

    // The square may hold a few more than n: drop the ones that shrink it most.
    while (static_cast<int>(layout.size()) > n) {
        layout.remove(layout.bestTreeToRemove());
    }

    return layout;
}

std::pair<Layout, double> bestCrystalCut(const Crystal &crystal, int n, int angles,
                                         std::mt19937_64 &rng) {
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const double spacing = 90.0 / std::max(1, angles);
    const double offset = spacing * unit(rng);

    // Evenly spaced angles, and the ones that lay a short lattice vector level: those cut
    // rows cleanly.
    std::vector<double> degs;
    for (int k = 0; k < std::max(1, angles); ++k) {
        degs.push_back(offset + k * spacing);
    }
    const Point &u = crystal.u, &v = crystal.v;
    for (const Point &edge : {u, v, Point{u.x + v.x, u.y + v.y}, Point{u.x - v.x, u.y - v.y}}) {
        degs.push_back(-std::atan2(edge.y, edge.x) * 180 / std::numbers::pi);
    }

    Layout best;
    double bestDeg = 0;
    for (const double deg : degs) {
        Layout cut = cutCrystal(turned(crystal, deg), n);
        if (best.size() == 0 || cut.side() < best.side()) {
            best = std::move(cut);
            bestDeg = deg;
        }
    }

    return {std::move(best), bestDeg};
}
} // namespace packer::packing
