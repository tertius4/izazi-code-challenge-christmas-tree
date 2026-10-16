#include "packing/squeezer.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <numbers>
#include <utility>
#include <vector>

namespace packer::packing {
namespace {
using geometry::Point;

// How far a box sticks out of the square, summed over its four sides.
double outside(const Box &t, const Box &square) {
    return std::max(0.0, square.minX - t.minX) + std::max(0.0, t.maxX - square.maxX) +
           std::max(0.0, square.minY - t.minY) + std::max(0.0, t.maxY - square.maxY);
}

// The trees that still overlap or stick out, kept so that any of them can be picked, added or
// dropped in constant time.
class BadSet {
  public:
    explicit BadSet(std::size_t n) : at_(n, kNone) {}

    void set(std::size_t i, bool bad) {
        const bool listed = at_[i] != kNone;
        if (bad && !listed) {
            at_[i] = list_.size();
            list_.push_back(i);
        } else if (!bad && listed) {
            const std::size_t last = list_.back();
            list_[at_[i]] = last;
            at_[last] = at_[i];
            list_.pop_back();
            at_[i] = kNone;
        }
    }

    bool empty() const {
        return list_.empty();
    }
    std::size_t size() const {
        return list_.size();
    }
    std::size_t operator[](std::size_t k) const {
        return list_[k];
    }

  private:
    static constexpr std::size_t kNone = static_cast<std::size_t>(-1);
    std::vector<std::size_t> list_;
    std::vector<std::size_t> at_;
};

// One unit of relax's angle variable turns the tree by 2 radians, which moves its tip (about
// 0.5 from where it turns) about as far as one unit of x or y moves it.
constexpr double kDegPerUnit = 2 * 180 / std::numbers::pi;

// The overlap energy of the squeeze as a function of every tree's pose at once:
// sum over pairs of overlap², plus sum over trees of how far they stick out², both counted
// with relaxMargin so that zero energy leaves a real gap.
class Relaxer {
  public:
    Relaxer(Layout &l, const Box &square, double margin)
        : l_(l), square_(square), margin_(margin) {}

    std::vector<double> get() const {
        std::vector<double> z(3 * l_.size());
        for (std::size_t i = 0; i < l_.size(); ++i) {
            z[3 * i] = l_.poses[i].x;
            z[3 * i + 1] = l_.poses[i].y;
            z[3 * i + 2] = l_.poses[i].deg / kDegPerUnit;
        }
        return z;
    }

    void set(const std::vector<double> &z) {
        for (std::size_t i = 0; i < l_.size(); ++i) {
            l_.setPose(i, {z[3 * i], z[3 * i + 1], z[3 * i + 2] * kDegPerUnit});
        }
    }

    // Energy and slope at the current poses; `clear` says whether the layout is already valid
    // and inside the square (by collide(), not the margin).
    double eval(std::vector<double> &g, bool &clear) const {
        const std::size_t n = l_.size();
        g.assign(3 * n, 0.0);
        clear = true;
        double energy = 0;

        for (std::size_t i = 0; i < n; ++i) {
            const Point o{l_.poses[i].x, l_.poses[i].y};
            double *gi = &g[3 * i];

            if (outside(l_.boxes[i], square_) > 0) clear = false;
            const Wall w = wall(i, o);
            if (w.depth > 0) {
                energy += w.depth * w.depth;
                gi[0] += 2 * w.depth * w.dx;
                gi[1] += 2 * w.depth * w.dy;
                gi[2] += 2 * w.depth * w.ddeg * kDegPerUnit;
            }

            l_.anyNeighbour(i, [&](std::size_t j) {
                const tree::Overlap ov = tree::overlap(l_.placed[i], o, l_.placed[j], margin_);
                if (ov.depth <= 0) return false;
                if (j > i) energy += ov.depth * ov.depth;
                // d(depth²)/d(pose of i): each pair is visited once from each side.
                gi[0] += 2 * ov.depth * ov.dx;
                gi[1] += 2 * ov.depth * ov.dy;
                gi[2] += 2 * ov.depth * ov.ddeg * kDegPerUnit;
                if (clear && tree::collide(l_.placed[i], l_.placed[j])) clear = false;
                return false;
            });
        }

        return energy;
    }

  private:
    Layout &l_;
    const Box square_;
    const double margin_;

    // How far tree i sticks out of the square shrunk by the margin, summed over the sides,
    // and its slope; the same shape as tree::Overlap.
    using Wall = tree::Overlap;
    Wall wall(std::size_t i, const Point &o) const {
        Wall w;
        const Box &b = l_.boxes[i];
        const double minX = square_.minX + margin_, maxX = square_.maxX - margin_;
        const double minY = square_.minY + margin_, maxY = square_.maxY - margin_;
        if (b.minX >= minX && b.maxX <= maxX && b.minY >= minY && b.maxY <= maxY) return w;

        // The vertex that reaches furthest each way; turning moves it by perp(v - o).
        constexpr double kRadPerDeg = std::numbers::pi / 180.0;
        Point lo{INFINITY, INFINITY}, hi{-INFINITY, -INFINITY};
        Point loX{}, hiX{}, loY{}, hiY{};
        for (const tree::Part &part : l_.placed[i].parts) {
            for (std::size_t k = 0; k < part.count; ++k) {
                const Point &v = part.v[k];
                if (v.x < lo.x) lo.x = v.x, loX = v;
                if (v.x > hi.x) hi.x = v.x, hiX = v;
                if (v.y < lo.y) lo.y = v.y, loY = v;
                if (v.y > hi.y) hi.y = v.y, hiY = v;
            }
        }

        if (b.minX < minX) {
            w.depth += minX - b.minX;
            w.dx -= 1;
            w.ddeg += (loX.y - o.y) * kRadPerDeg; // d(v.x)/dθ = -(v.y - o.y)
        }
        if (b.maxX > maxX) {
            w.depth += b.maxX - maxX;
            w.dx += 1;
            w.ddeg -= (hiX.y - o.y) * kRadPerDeg;
        }
        if (b.minY < minY) {
            w.depth += minY - b.minY;
            w.dy -= 1;
            w.ddeg -= (loY.x - o.x) * kRadPerDeg; // d(v.y)/dθ = v.x - o.x
        }
        if (b.maxY > maxY) {
            w.depth += b.maxY - maxY;
            w.dy += 1;
            w.ddeg += (hiY.x - o.x) * kRadPerDeg;
        }

        return w;
    }
};

double dot(const std::vector<double> &u, const std::vector<double> &v) {
    double sum = 0;
    for (std::size_t k = 0; k < u.size(); ++k) {
        sum += u[k] * v[k];
    }
    return sum;
}

// Moves every tree at once down the slope of the overlap energy (L-BFGS with a backtracking
// line search) until the layout is valid inside the square, or it gets stuck. Returns whether
// the layout came out valid; either way `l` holds the lowest-energy poses found.
bool relax(Layout &l, const Box &square, const SqueezeConfig &cfg) {
    constexpr std::size_t kMemory = 8;
    Relaxer r(l, square, cfg.relaxMargin);

    std::vector<double> z = r.get(), g, zNext, gNext, d;
    bool clear = false;
    double energy = r.eval(g, clear);
    if (clear) return true;

    struct Pair {
        std::vector<double> s, y;
        double rho;
    };
    std::deque<Pair> memory;
    std::vector<double> alpha(kMemory);

    for (long it = 0; it < cfg.relaxIters; ++it) {
        if (cfg.stop && cfg.stop->load(std::memory_order_relaxed)) break;

        // Direction: -H·g by the two-loop recursion.
        d = g;
        for (std::size_t k = memory.size(); k-- > 0;) {
            alpha[k] = memory[k].rho * dot(memory[k].s, d);
            for (std::size_t q = 0; q < d.size(); ++q) d[q] -= alpha[k] * memory[k].y[q];
        }
        double scale = 0;
        if (memory.empty()) {
            // No curvature known yet: a first step of at most 1e-3 along any coordinate.
            double most = 0;
            for (const double v : g) most = std::max(most, std::abs(v));
            scale = most > 0 ? 1e-3 / most : 0;
        } else {
            const Pair &last = memory.back();
            scale = dot(last.s, last.y) / dot(last.y, last.y);
        }
        for (double &v : d) v *= scale;
        for (std::size_t k = 0; k < memory.size(); ++k) {
            const double beta = memory[k].rho * dot(memory[k].y, d);
            for (std::size_t q = 0; q < d.size(); ++q) d[q] += (alpha[k] - beta) * memory[k].s[q];
        }
        for (double &v : d) v = -v;

        // The energy is nearly flat in places, which can make the curvature estimate huge:
        // never move anything further than a small step at once, or trees fly off (and the
        // layout's grid grows to cover them).
        constexpr double kMaxStep = 0.02;
        double longest = 0;
        for (const double v : d) longest = std::max(longest, std::abs(v));
        if (longest > kMaxStep) {
            for (double &v : d) v *= kMaxStep / longest;
        }

        double slope = dot(g, d);
        if (!(slope < 0)) {
            if (memory.empty()) break; // flat: nothing left to follow
            memory.clear();
            continue;
        }

        // Backtrack until the energy drops enough.
        double t = 1;
        bool moved = false;
        double eNext = 0;
        bool clearNext = false;
        zNext.resize(z.size());
        for (int tries = 0; tries < 40; ++tries, t /= 2) {
            for (std::size_t q = 0; q < z.size(); ++q) zNext[q] = z[q] + t * d[q];
            r.set(zNext);
            eNext = r.eval(gNext, clearNext);
            if (eNext <= energy + 1e-4 * t * slope || clearNext) {
                moved = true;
                break;
            }
        }

        if (!moved) {
            r.set(z);
            if (memory.empty()) return false; // stuck at a local minimum with overlap left
            memory.clear();
            continue;
        }

        Pair p{std::vector<double>(z.size()), std::vector<double>(z.size()), 0};
        for (std::size_t q = 0; q < z.size(); ++q) {
            p.s[q] = zNext[q] - z[q];
            p.y[q] = gNext[q] - g[q];
        }
        const double sy = dot(p.s, p.y);
        if (sy > 1e-30) {
            p.rho = 1 / sy;
            memory.push_back(std::move(p));
            if (memory.size() > kMemory) memory.pop_front();
        }

        std::swap(z, zNext);
        std::swap(g, gNext);
        energy = eNext;
        if (clearNext) return true;
    }

    r.set(z);
    return false;
}
} // namespace

std::optional<Layout> squeeze(const Layout &layout, double target, const SqueezeConfig &cfg,
                              std::mt19937_64 &rng) {
    const std::size_t n = layout.size();
    if (n == 0 || cfg.iterations <= 0) return std::nullopt;

    // Pull the trees toward the centre so their box fits the square; angles stay. Each axis
    // shrinks only as much as it needs, so a narrow layout gets no overlap across its width.
    Layout l = layout;
    const Box b = l.bounds();
    const double cx = (b.minX + b.maxX) / 2;
    const double cy = (b.minY + b.maxY) / 2;
    const double fx = std::min(1.0, target / (b.maxX - b.minX));
    const double fy = std::min(1.0, target / (b.maxY - b.minY));
    for (std::size_t i = 0; i < n; ++i) {
        const Pose p = l.poses[i];
        l.setPose(i, {cx + (p.x - cx) * fx, cy + (p.y - cy) * fy, p.deg});
    }

    const Box square{cx - target / 2, cy - target / 2, cx + target / 2, cy + target / 2};

    // A valid layout inside the square, checked independently of the penalty.
    auto fits = [&] { return l.side() <= target && l.isValid() ? std::optional(l) : std::nullopt; };

    if (cfg.relaxIters > 0 && relax(l, square, cfg)) return fits();

    // Penalty: the sum of every overlapping pair's penetration squared, plus how far each tree
    // sticks out squared. Each tree counts its own violations, so a tree is clear at count 0.
    std::vector<int> violations(n, 0);
    BadSet bad(n);
    double energy = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const double out = outside(l.boxes[i], square);
        energy += out * out;
        violations[i] += out > 0;
        l.anyNeighbour(i, [&](std::size_t j) {
            const double d = tree::penetration(l.placed[i], l.placed[j]);
            if (j > i) energy += d * d;
            violations[i] += d > 0;
            return false;
        });
        bad.set(i, violations[i] > 0);
    }

    if (bad.empty()) return fits();

    // Temperatures scale with the starting penalty per bad tree.
    const double scale = energy / static_cast<double>(bad.size());

    std::uniform_int_distribution<std::size_t> pickTree(0, n - 1);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    auto jiggle = [&](double s) { return s * (2 * unit(rng) - 1); };

    // Each move's penetration with its neighbours, before and after: (j, depth).
    std::vector<std::pair<std::size_t, double>> before, after;
    auto contacts = [&](std::size_t i, std::vector<std::pair<std::size_t, double>> &into) {
        into.clear();
        double sum = 0;
        l.anyNeighbour(i, [&](std::size_t j) {
            const double d = tree::penetration(l.placed[i], l.placed[j]);
            if (d > 0) {
                into.push_back({j, d});
                sum += d * d;
            }
            return false;
        });
        return sum;
    };

    for (long it = 0; it < cfg.iterations; ++it) {
        const bool checkStop = cfg.stop && it % 4096 == 0;
        if (checkStop && cfg.stop->load(std::memory_order_relaxed)) break;

        const double t = static_cast<double>(it) / cfg.iterations;
        const double temp = scale * cfg.startTemp * std::pow(cfg.endTemp / cfg.startTemp, t);
        const double step = cfg.startStep * std::pow(cfg.endStep / cfg.startStep, t);

        const bool fixBad = unit(rng) < cfg.badShare;
        const std::size_t i =
            fixBad ? bad[std::uniform_int_distribution<std::size_t>(0, bad.size() - 1)(rng)]
                   : pickTree(rng);
        const Pose old = l.poses[i];
        Pose moved = old;
        const double kind = unit(rng);
        if (kind < 0.6) {
            moved.x += jiggle(step);
            moved.y += jiggle(step);
        } else if (kind < 0.98) {
            moved.deg += jiggle(100 * step);
        } else {
            moved.deg += 180;
        }
        moved.deg = geometry::normalizeDeg(moved.deg);

        const double outBefore = outside(l.boxes[i], square);
        const double pairsBefore = contacts(i, before);
        l.setPose(i, moved);
        const double outAfter = outside(l.boxes[i], square);
        const double pairsAfter = contacts(i, after);

        const double worseBy =
            pairsAfter - pairsBefore + outAfter * outAfter - outBefore * outBefore;
        const bool accept = worseBy <= 0 || unit(rng) < std::exp(-worseBy / temp);
        if (!accept) {
            l.setPose(i, old);
            continue;
        }

        energy += worseBy;
        // Recount the violations of tree i and of every neighbour it touched before or after.
        for (const auto &[j, d] : before) {
            violations[j] -= 1;
        }
        for (const auto &[j, d] : after) {
            violations[j] += 1;
        }
        for (const auto &[j, d] : before) {
            bad.set(j, violations[j] > 0);
        }
        for (const auto &[j, d] : after) {
            bad.set(j, violations[j] > 0);
        }
        violations[i] = static_cast<int>(after.size()) + (outAfter > 0);
        bad.set(i, violations[i] > 0);

        if (bad.empty()) break;
    }

    // Single-tree moves left overlap: one more relax from where they got to.
    if (!bad.empty() && !(cfg.relaxIters > 0 && relax(l, square, cfg))) return std::nullopt;

    return fits();
}

Layout shake(const Layout &layout, int kicks, std::mt19937_64 &rng) {
    Layout l = layout;
    const std::size_t n = l.size();
    if (n < 2) return l;

    std::uniform_int_distribution<std::size_t> pickTree(0, n - 1);
    std::uniform_real_distribution<double> unit(0.0, 1.0);

    for (int k = 0; k < kicks; ++k) {
        const std::size_t i = pickTree(rng);
        const Pose pi = l.poses[i];
        const double kind = unit(rng);

        if (kind < 0.5) {
            // Swapping two trees turned the same way changes nothing: only differently turned
            // neighbours are worth it.
            std::vector<std::size_t> others;
            l.anyNeighbour(i, [&](std::size_t j) {
                const double turn = geometry::normalizeDeg(l.poses[j].deg - pi.deg);
                if (std::min(turn, 360 - turn) > 20) others.push_back(j);
                return false;
            });

            if (!others.empty()) {
                const std::size_t j =
                    others[std::uniform_int_distribution<std::size_t>(0, others.size() - 1)(rng)];
                const Pose pj = l.poses[j];
                l.setPose(i, {pj.x, pj.y, pi.deg});
                l.setPose(j, {pi.x, pi.y, pj.deg});
                continue;
            }
        }

        if (kind < 0.8) {
            l.setPose(i, {pi.x, pi.y, geometry::normalizeDeg(pi.deg + 180)});
        } else {
            const Box b = l.bounds();
            l.setPose(i, {b.minX + unit(rng) * (b.maxX - b.minX),
                          b.minY + unit(rng) * (b.maxY - b.minY), 360 * unit(rng)});
        }
    }

    return l;
}
} // namespace packer::packing
