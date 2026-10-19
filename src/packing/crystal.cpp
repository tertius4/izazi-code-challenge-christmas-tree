#include "packing/crystal.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numbers>

#include "tree/tree.h"
#include "util/parallel.h"

namespace packer::packing {
namespace {
using geometry::Point;

constexpr double kTol = 1e-9; // slack when testing whether a box lies inside a square
constexpr double kReach = 0.8; // largest distance from a tree's origin to its outline (the tip)

double cross(const Point &a, const Point &b) {
    return a.x * b.y - a.y * b.x;
}

Point at(const Crystal &c, double fx, double fy) {
    return {fx * c.u.x + fy * c.v.x, fx * c.u.y + fy * c.v.y};
}

// The point's coordinates as fractions of u and v.
Point fractions(const Crystal &c, const Point &p) {
    const double det = cross(c.u, c.v);
    return {cross(p, c.v) / det, cross(c.u, p) / det};
}

Box shifted(const Box &b, double dx, double dy) {
    return {b.minX + dx, b.minY + dy, b.maxX + dx, b.maxY + dy};
}

double wrap01(double x) {
    return x - std::floor(x);
}

Crystal scaled(const Crystal &c, double f) {
    Crystal out = c;
    out.u = {c.u.x * f, c.u.y * f};
    out.v = {c.v.x * f, c.v.y * f};
    return out;
}

// Scales the cell (the trees keep their places in it) to about the smallest size that is still
// valid, growing it first if it collides. Returns false if no valid size was found.
bool tighten(Crystal &c) {
    double hi = 1.0;
    double lo = 0;
    if (crystalIsValid(c)) {
        // Shrink in steps until it collides, so the bisection below starts with a tight bracket.
        lo = hi;
        for (int i = 0; i < 60; ++i) {
            lo *= 0.97;
            if (!crystalIsValid(scaled(c, lo))) break;
            hi = lo;
        }
    } else {
        for (int i = 0;; ++i) {
            if (i == 100) return false; // two trees sit on top of each other: no size works
            lo = hi;
            hi *= 1.03;
            if (crystalIsValid(scaled(c, hi))) break;
        }
    }

    for (int i = 0; i < 28; ++i) {
        const double mid = (lo + hi) / 2;
        if (crystalIsValid(scaled(c, mid))) {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    // Back off by a hair, so rounding far from the origin can never turn a touch into an overlap.
    const Crystal safer = scaled(c, hi * (1 + 1e-9));
    c = crystalIsValid(safer) ? safer : scaled(c, hi);

    return true;
}

// The cell of the search: u = (a, 0) and v = (h·a, b). Turning the whole pattern changes
// nothing about its density, so u stays level.
struct Shape {
    double a = 0;
    double b = 0;
    double h = 0;
};

Crystal build(const Shape &s, const std::vector<Pose> &cell) {
    return {{s.a, 0}, {s.h * s.a, s.b}, cell};
}

Shape shapeOf(const Crystal &c) {
    return {c.u.x, c.v.y, c.v.x / c.u.x};
}
} // namespace

double areaPerTree(const Crystal &c) {
    return std::abs(cross(c.u, c.v)) / static_cast<double>(c.cell.size());
}

Pose crystalPose(const Crystal &c, int i, int j, std::size_t k) {
    const Pose &t = c.cell[k];
    const Point p = at(c, i + t.x, j + t.y);
    return {p.x, p.y, t.deg};
}

bool crystalIsValid(const Crystal &c) {
    const double det = cross(c.u, c.v);
    if (std::abs(det) < 1e-12) return false;

    // A vector w has fractions (w×v, u×w) / det, so |fraction along u| <= |w|·|v| / |det| and
    // |fraction along v| <= |w|·|u| / |det|. Trees only touch if their origins are within
    // 2 reaches: that bounds the cells to check.
    const double range = 2 * kReach;
    const double spanI = range * std::hypot(c.v.x, c.v.y) / std::abs(det);
    const double spanJ = range * std::hypot(c.u.x, c.u.y) / std::abs(det);

    const std::size_t k = c.cell.size();
    std::vector<tree::Placed> base(k);
    for (std::size_t a = 0; a < k; ++a) {
        base[a] = tree::place(crystalPose(c, 0, 0, a));
    }

    tree::Placed other;
    for (std::size_t a = 0; a < k; ++a) {
        const Pose pa = crystalPose(c, 0, 0, a);
        // Pairs (a, b) and (b, a) are the same pairs of trees, moved by a whole cell.
        for (std::size_t b = a; b < k; ++b) {
            const double fx = c.cell[b].x - c.cell[a].x;
            const double fy = c.cell[b].y - c.cell[a].y;
            const int minI = static_cast<int>(std::ceil(-fx - spanI));
            const int maxI = static_cast<int>(std::floor(-fx + spanI));
            const int minJ = static_cast<int>(std::ceil(-fy - spanJ));
            const int maxJ = static_cast<int>(std::floor(-fy + spanJ));

            for (int i = minI; i <= maxI; ++i) {
                for (int j = minJ; j <= maxJ; ++j) {
                    if (a == b && i == 0 && j == 0) continue;

                    const Pose pb = crystalPose(c, i, j, b);
                    const double dx = pb.x - pa.x;
                    const double dy = pb.y - pa.y;
                    if (dx * dx + dy * dy > range * range) continue;

                    tree::place(pb, other);
                    if (tree::collide(base[a], other)) return false;
                }
            }
        }
    }

    return true;
}

Crystal turned(const Crystal &c, double deg) {
    const double rad = deg * std::numbers::pi / 180.0;
    const double cs = std::cos(rad);
    const double sn = std::sin(rad);
    auto turn = [&](const Point &p) { return Point{p.x * cs - p.y * sn, p.x * sn + p.y * cs}; };

    Crystal out = c;
    out.u = turn(c.u);
    out.v = turn(c.v);
    for (Pose &t : out.cell) {
        t.deg = geometry::normalizeDeg(t.deg + deg); // fractions turn with u and v
    }

    return out;
}

Crystal doubled(const Crystal &c) {
    Crystal out = c;
    out.u = {2 * c.u.x, 2 * c.u.y};
    out.cell.clear();
    for (const Pose &t : c.cell) {
        out.cell.push_back({t.x / 2, t.y, t.deg});
    }
    for (const Pose &t : c.cell) {
        out.cell.push_back({t.x / 2 + 0.5, t.y, t.deg});
    }
    return out;
}

Crystal toCrystal(const LatticeParams &p) {
    Crystal c;
    c.u = {p.a, 0};
    c.v = {p.h * p.a, p.b};
    c.cell = {{0, 0, p.deg}, {p.s, p.t, p.deg + 180}};
    return c;
}

std::optional<Crystal> findCrystal(int k, bool paired, long iterations, std::mt19937_64 &rng,
                                   const std::atomic<bool> *stop, const Crystal *start) {
    if (start) k = static_cast<int>(start->cell.size());
    if (k < 1) return std::nullopt;

    std::uniform_real_distribution<double> unit(0.0, 1.0);
    auto jiggle = [&](double step) { return step * (2 * unit(rng) - 1); };
    const auto count = static_cast<std::size_t>(k);
    auto partnerOf = [&](std::size_t m) -> std::size_t {
        if (!paired) return m;
        const std::size_t p = m ^ 1u;
        return p < count ? p : m;
    };

    // A random roomy cell, so that tighten() only has to shrink it.
    Shape shape;
    std::vector<Pose> cell(count);
    Crystal c;
    for (int attempt = 0; start == nullptr; ++attempt) {
        if (attempt == 100) return std::nullopt;

        const double side = std::sqrt(static_cast<double>(k)) * (1 + unit(rng));
        shape = {side, side * (0.7 + 0.6 * unit(rng)), unit(rng) - 0.5};
        for (std::size_t m = 0; m < count; ++m) {
            cell[m] = {unit(rng), unit(rng), 360 * unit(rng)};
            if (paired && m % 2 == 1) cell[m].deg = cell[m - 1].deg + 180;
        }
        cell[0].x = cell[0].y = 0; // moving every tree together changes nothing

        c = build(shape, cell);
        if (tighten(c)) break;
    }
    if (start) {
        // Level u, as the moves below expect: turning the pattern changes nothing.
        c = turned(*start, -std::atan2(start->u.y, start->u.x) * 180 / std::numbers::pi);
        if (!tighten(c)) return std::nullopt;
    }

    Crystal best = c;
    double current = areaPerTree(c);

    for (long it = 0; it < iterations; ++it) {
        if (stop && it % 256 == 0 && stop->load(std::memory_order_relaxed)) break;

        const double progress = static_cast<double>(it) / iterations;
        // A given start is already good: search near it, cooler, rather than melt it.
        const double heat = start ? 0.1 : 1.0;
        const double step = heat * 0.05 * std::pow(0.002, progress); // 0.05 down to 0.0001
        const double temp = heat * heat * 0.01 * std::pow(1e-4, progress);

        Shape s = shapeOf(c);
        std::vector<Pose> q = c.cell;
        const std::size_t m = std::uniform_int_distribution<std::size_t>(0, count - 1)(rng);
        const double kind = unit(rng);

        if (kind < 0.2) { // stretch one way, squash the other
            const double r = std::exp(jiggle(step));
            s.a *= r;
            s.b /= r;
        } else if (kind < 0.35) {
            s.h += jiggle(step);
        } else if (kind < 0.65 && count > 1) {
            const std::size_t t = 1 + m % (count - 1); // tree 0 stays at the corner
            q[t].x = wrap01(q[t].x + jiggle(step));
            q[t].y = wrap01(q[t].y + jiggle(step));
        } else if (kind < 0.9) {
            q[m].deg += 30 * jiggle(step);
            if (partnerOf(m) != m) q[partnerOf(m)].deg = q[m].deg + 180;
        } else if (kind < 0.93) {
            q[m].deg += 180;
            if (partnerOf(m) != m) q[partnerOf(m)].deg += 180;
        } else { // turn every tree together, against the cell
            const double turn = 30 * jiggle(step);
            for (Pose &t : q) t.deg += turn;
        }

        Crystal next = build(s, q);
        if (!tighten(next)) continue;

        const double area = areaPerTree(next);
        const double worseBy = area - current;
        const bool accept = worseBy < 0 || unit(rng) < std::exp(-worseBy / temp);
        if (accept) {
            c = std::move(next);
            current = area;
            if (area < areaPerTree(best)) best = c;
        }
    }

    for (Pose &t : best.cell) {
        t.deg = geometry::normalizeDeg(t.deg);
    }

    return best;
}

std::vector<Crystal> findCrystals(int searches, long iterations, int keep, int threads,
                                  unsigned seed, const std::atomic<bool> *stop) {
    // The kinds of cell to try, in turn: the double lattice, then freer and bigger cells.
    struct Kind {
        int k;
        bool paired;
    };
    const Kind kinds[] = {{2, true}, {4, true}, {2, false}, {3, false}, {4, false}};

    std::vector<std::optional<Crystal>> found(static_cast<std::size_t>(2 * searches));

    // Runs search(s) for s in [first, last) on the threads.
    auto runAll = [&](int first, int last, const auto &search) {
        const auto count = static_cast<std::size_t>(last - first);
        util::parallelFor(count, std::max(1, threads), [&](std::size_t k) {
            if (stop && stop->load()) return;
            const int s = first + static_cast<int>(k);
            std::mt19937_64 rng(seed * 1000003ULL + static_cast<unsigned>(s));
            found[static_cast<std::size_t>(s)] = search(s, rng);
        });
    };

    // First from random cells. Bigger cells have more to arrange: give them longer.
    runAll(0, searches, [&](int s, std::mt19937_64 &rng) {
        const Kind &kind = kinds[s % std::size(kinds)];
        return findCrystal(kind.k, kind.paired, iterations * kind.k / 2, rng, stop);
    });

    // Then from the densest found so far, doubled: a cell of twice the trees is at least as
    // dense, and its copies are free to drift apart into something denser.
    std::vector<Crystal> seeds;
    for (int s = 0; s < searches; ++s) {
        if (found[static_cast<std::size_t>(s)]) seeds.push_back(*found[static_cast<std::size_t>(s)]);
    }
    std::sort(seeds.begin(), seeds.end(), [](const Crystal &a, const Crystal &b) {
        return areaPerTree(a) < areaPerTree(b);
    });
    if (!seeds.empty()) {
        runAll(searches, 2 * searches, [&](int s, std::mt19937_64 &rng) {
            const Crystal start = doubled(seeds[static_cast<std::size_t>(s) % std::min<std::size_t>(seeds.size(), 3)]);
            return findCrystal(0, false, iterations, rng, stop, &start);
        });
    }

    std::vector<Crystal> sorted;
    for (const auto &c : found) {
        if (c) sorted.push_back(*c);
    }
    std::sort(sorted.begin(), sorted.end(), [](const Crystal &a, const Crystal &b) {
        return areaPerTree(a) < areaPerTree(b);
    });

    // Searches that end in about the same pattern end at about the same density: keep one of
    // each, so the kept ones cut differently.
    std::vector<Crystal> kept;
    for (const Crystal &c : sorted) {
        if (static_cast<int>(kept.size()) == keep) break;
        const bool isNew = kept.empty() || areaPerTree(c) > areaPerTree(kept.back()) * (1 + 1e-3);
        if (isNew) kept.push_back(c);
    }

    for (const Crystal &c : kept) {
        std::printf("crystal with %zu trees per cell: area per tree %.6f\n", c.cell.size(),
                    areaPerTree(c));
    }
    std::fflush(stdout);

    return kept;
}

Layout cutCrystal(const Crystal &c, int n) {
    const std::size_t k = c.cell.size();

    struct Tree {
        Pose pose;
        Box box;
    };

    std::vector<Box> baseBox(k);
    std::vector<Tree> anchors(k); // the trees of cell (0, 0)
    for (std::size_t m = 0; m < k; ++m) {
        baseBox[m] = tree::place({0, 0, c.cell[m].deg}).box;
        const Pose p = crystalPose(c, 0, 0, m);
        anchors[m] = {p, shifted(baseBox[m], p.x, p.y)};
    }

    Box anchorSpan = anchors[0].box;
    for (const Tree &t : anchors) {
        anchorSpan.minX = std::min(anchorSpan.minX, t.pose.x);
        anchorSpan.minY = std::min(anchorSpan.minY, t.pose.y);
        anchorSpan.maxX = std::max(anchorSpan.maxX, t.pose.x);
        anchorSpan.maxY = std::max(anchorSpan.maxY, t.pose.y);
    }

    std::vector<Tree> trees;

    // A square that fits whole trees can slide down until its bottom edge touches a tree, then
    // move by whole cells until that tree is one of cell (0, 0). So the square only needs
    // trying with its bottom on an anchor, holding that anchor: every tree such a square of
    // this side could hold lies within this region.
    auto generate = [&](double side) {
        trees.clear();
        const Box region{anchorSpan.minX - side - 2 * kReach, anchorSpan.minY - 2 * kReach,
                         anchorSpan.maxX + side + 2 * kReach, anchorSpan.maxY + side + 2 * kReach};

        double minI = INFINITY, maxI = -INFINITY, minJ = INFINITY, maxJ = -INFINITY;
        for (const Point &corner : {Point{region.minX, region.minY}, Point{region.maxX, region.minY},
                                    Point{region.minX, region.maxY}, Point{region.maxX, region.maxY}}) {
            const Point f = fractions(c, corner);
            minI = std::min(minI, f.x);
            maxI = std::max(maxI, f.x);
            minJ = std::min(minJ, f.y);
            maxJ = std::max(maxJ, f.y);
        }

        for (int i = static_cast<int>(std::floor(minI)) - 1; i <= static_cast<int>(std::ceil(maxI)) + 1; ++i) {
            for (int j = static_cast<int>(std::floor(minJ)) - 1; j <= static_cast<int>(std::ceil(maxJ)) + 1; ++j) {
                for (std::size_t m = 0; m < k; ++m) {
                    const Pose p = crystalPose(c, i, j, m);
                    const bool inRegion = p.x >= region.minX && p.x <= region.maxX &&
                                          p.y >= region.minY && p.y <= region.maxY;
                    if (inRegion) trees.push_back({p, shifted(baseBox[m], p.x, p.y)});
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
        for (const Tree &anchor : anchors) {
            const double y0 = anchor.box.minY;
            band.clear();
            for (const Tree &t : trees) {
                const bool fitsY = t.box.minY >= y0 - kTol && t.box.maxY <= y0 + side + kTol;
                if (fitsY) band.push_back({t.box.minX, t.box.maxX});
            }
            std::sort(band.begin(), band.end());

            const double firstX = anchor.box.maxX - side - kTol;
            const double lastX = anchor.box.minX + kTol;
            const auto from = std::lower_bound(band.begin(), band.end(), std::pair{firstX, -HUGE_VAL});
            for (auto it = from; it != band.end() && it->first <= lastX; ++it) {
                const double x0 = it->first;
                int count = 0;
                for (auto t = std::lower_bound(band.begin(), band.end(), std::pair{x0 - kTol, -HUGE_VAL});
                     t != band.end() && t->first <= x0 + side + kTol; ++t) {
                    if (t->second <= x0 + side + kTol) ++count;
                }
                if (count > best.count) best = {x0, y0, count};
            }
        }
        return best;
    };

    // A side that surely holds n trees, then bisect down to the smallest one that does.
    double hi = std::sqrt(n * areaPerTree(c)) * 1.3 + 1.5;
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

    const Window w = bestWindow(hi);
    Layout layout;
    for (const Tree &t : trees) {
        const bool fitsX = t.box.minX >= w.x0 - kTol && t.box.maxX <= w.x0 + hi + kTol;
        const bool fitsY = t.box.minY >= w.y0 - kTol && t.box.maxY <= w.y0 + hi + kTol;
        if (fitsX && fitsY) layout.add(t.pose);
    }

    // The square may hold a few more than n: drop the ones that shrink it most.
    while (static_cast<int>(layout.size()) > n) {
        layout.remove(layout.bestTreeToRemove());
    }

    return layout;
}

std::pair<Layout, double> bestCrystalCut(const Crystal &c, int n, int angles,
                                         std::mt19937_64 &rng) {
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const double spacing = 90.0 / std::max(1, angles);
    const double offset = spacing * unit(rng);

    // Evenly spaced angles, and the ones that lay a short lattice vector level: those cut
    // rows cleanly.
    std::vector<double> degs;
    for (int a = 0; a < std::max(1, angles); ++a) {
        degs.push_back(offset + a * spacing);
    }
    for (const Point &w : {c.u, c.v, Point{c.u.x + c.v.x, c.u.y + c.v.y},
                           Point{c.u.x - c.v.x, c.u.y - c.v.y}}) {
        degs.push_back(-std::atan2(w.y, w.x) * 180 / std::numbers::pi);
    }

    Layout best;
    double bestDeg = 0;
    for (const double deg : degs) {
        Layout cut = cutCrystal(turned(c, deg), n);
        if (best.size() == 0 || cut.side() < best.side()) {
            best = std::move(cut);
            bestDeg = deg;
        }
    }

    return {std::move(best), bestDeg};
}
} // namespace packer::packing
