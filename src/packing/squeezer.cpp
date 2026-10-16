#include "packing/squeezer.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <numbers>
#include <unordered_map>
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

// One unit of the angle variable turns the tree by 2 radians, which moves its tip (about
// 0.5 from where it turns) about as far as one unit of x or y moves it.
constexpr double kDegPerUnit = 2 * 180 / std::numbers::pi;

// The overlap energy as a function of every free tree's pose at once (and, when compressing,
// of the square's side): sum over pairs of weight·overlap², plus sum over trees of how far
// they stick out², all counted with a margin so that zero energy leaves a real gap.
//
// The variables come in blocks of three (x, y, angle), one per free tree. A tree either has
// its own block, copies a block mirrored through the square's centre (symmetric mode), or has
// none and stays put.
class Energy {
  public:
    struct Options {
        double margin = 1e-8;
        const std::vector<bool> *movable = nullptr; // null or empty: every tree moves
        bool symmetric = false;
        bool sideIsVariable = false; // compress: the energy is side + weight·penalty
    };

    Energy(Layout &l, const Box &square, const Options &o)
        : l_(l), margin_(o.margin), sideIsVariable_(o.sideIsVariable),
          cx_((square.minX + square.maxX) / 2), cy_((square.minY + square.maxY) / 2),
          half_((square.maxX - square.minX) / 2) {
        const std::size_t n = l.size();
        block_.assign(n, -1);
        mirrored_.assign(n, false);
        const bool all = o.movable == nullptr || o.movable->empty();
        auto moves = [&](std::size_t i) { return all || (*o.movable)[i]; };

        const std::size_t pairs = o.symmetric ? n / 2 : 0;
        for (std::size_t i = 0; i < pairs; ++i) {
            block_[i] = block_[i + pairs] = static_cast<int>(primary_.size());
            mirrored_[i + pairs] = true;
            primary_.push_back(i);
        }
        for (std::size_t i = 2 * pairs; i < n; ++i) {
            if (!moves(i)) continue;
            block_[i] = static_cast<int>(primary_.size());
            primary_.push_back(i);
        }
    }

    std::size_t size() const {
        return 3 * primary_.size() + (sideIsVariable_ ? 1 : 0);
    }

    double side() const {
        return 2 * half_;
    }

    Box square() const {
        return {cx_ - half_, cy_ - half_, cx_ + half_, cy_ + half_};
    }

    void setWeight(double w) {
        weight_ = w;
    }

    std::vector<double> get() const {
        std::vector<double> z(size());
        for (std::size_t b = 0; b < primary_.size(); ++b) {
            const Pose &p = l_.poses[primary_[b]];
            z[3 * b] = p.x;
            z[3 * b + 1] = p.y;
            z[3 * b + 2] = p.deg / kDegPerUnit;
        }
        if (sideIsVariable_) z.back() = 2 * half_;
        return z;
    }

    void set(const std::vector<double> &z) {
        if (sideIsVariable_) half_ = z.back() / 2;
        for (std::size_t i = 0; i < l_.size(); ++i) {
            const int b = block_[i];
            if (b < 0) continue;
            const double x = z[3 * b], y = z[3 * b + 1], deg = z[3 * b + 2] * kDegPerUnit;
            l_.setPose(i, mirrored_[i] ? Pose{2 * cx_ - x, 2 * cy_ - y, deg + 180}
                                       : Pose{x, y, deg});
        }
    }

    // Energy and slope at the current poses; `clear` says whether the layout is already valid
    // and inside the square (by collide(), not the margin).
    double eval(std::vector<double> &g, bool &clear) const {
        const std::size_t n = l_.size();
        const Box sq = square();
        perTree_.assign(3 * n, 0.0);
        clear = true;
        double penalty = 0;
        double sideSlope = 0; // d(penalty)/d(side)

        for (std::size_t i = 0; i < n; ++i) {
            const Point o{l_.poses[i].x, l_.poses[i].y};
            double *gi = &perTree_[3 * i];

            if (outside(l_.boxes[i], sq) > 0) clear = false;
            const Wall w = wall(i, o, sq);
            if (w.depth > 0) {
                penalty += w.depth * w.depth;
                gi[0] += 2 * w.depth * w.dx;
                gi[1] += 2 * w.depth * w.dy;
                gi[2] += 2 * w.depth * w.ddeg * kDegPerUnit;
                sideSlope += 2 * w.depth * (-0.5 * w.sides); // each wall moves out by side/2
            }

            l_.anyNeighbour(i, [&](std::size_t j) {
                if (block_[i] < 0 && block_[j] < 0) return false; // neither moves
                const tree::Overlap ov = tree::overlap(l_.placed[i], o, l_.placed[j], margin_);
                if (ov.depth <= 0) return false;
                const double w = pairWeight(i, j);
                if (j > i) penalty += w * ov.depth * ov.depth;
                // d(depth²)/d(pose of i): each pair is visited once from each side.
                gi[0] += 2 * w * ov.depth * ov.dx;
                gi[1] += 2 * w * ov.depth * ov.dy;
                gi[2] += 2 * w * ov.depth * ov.ddeg * kDegPerUnit;
                if (clear && tree::collide(l_.placed[i], l_.placed[j])) clear = false;
                return false;
            });
        }

        // Fold each tree's slope into its block: a mirrored tree moves against its block's x
        // and y, and turns with it.
        const double scale = sideIsVariable_ ? weight_ : 1.0;
        g.assign(size(), 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            const int b = block_[i];
            if (b < 0) continue;
            const double sign = mirrored_[i] ? -1.0 : 1.0;
            g[3 * b] += scale * sign * perTree_[3 * i];
            g[3 * b + 1] += scale * sign * perTree_[3 * i + 1];
            g[3 * b + 2] += scale * perTree_[3 * i + 2];
        }

        if (!sideIsVariable_) return penalty;
        g.back() = 1 + weight_ * sideSlope;
        return 2 * half_ + weight_ * penalty;
    }

    // Guided local search: raise the weight of every pair that still overlaps, the deepest by
    // 1 and the others in proportion. Returns false if nothing overlaps.
    bool raiseWeights() {
        std::vector<std::pair<std::uint64_t, double>> hits;
        double deepest = 0;
        for (std::size_t i = 0; i < l_.size(); ++i) {
            const Point o{l_.poses[i].x, l_.poses[i].y};
            l_.anyNeighbour(i, [&](std::size_t j) {
                if (j < i) return false;
                const double d = tree::overlap(l_.placed[i], o, l_.placed[j], margin_).depth;
                if (d > 0) {
                    hits.push_back({key(i, j), d});
                    deepest = std::max(deepest, d);
                }
                return false;
            });
        }

        for (const auto &[k, d] : hits) {
            auto it = weights_.try_emplace(k, 1.0).first;
            it->second += d / deepest;
        }

        return !hits.empty();
    }

  private:
    Layout &l_;
    const double margin_;
    const bool sideIsVariable_;
    const double cx_, cy_;
    double half_;
    double weight_ = 1;
    std::vector<int> block_;
    std::vector<bool> mirrored_;
    std::vector<std::size_t> primary_; // the tree each block places as it is
    std::unordered_map<std::uint64_t, double> weights_; // raised pair weights; the rest are 1
    mutable std::vector<double> perTree_;

    std::uint64_t key(std::size_t i, std::size_t j) const {
        if (j < i) std::swap(i, j);
        return static_cast<std::uint64_t>(i) * l_.size() + j;
    }

    double pairWeight(std::size_t i, std::size_t j) const {
        if (weights_.empty()) return 1;
        const auto it = weights_.find(key(i, j));
        return it == weights_.end() ? 1 : it->second;
    }

    // How far tree i sticks out of the square shrunk by the margin, summed over the sides,
    // its slope, and how many sides it sticks out of.
    struct Wall {
        double depth = 0;
        double dx = 0;
        double dy = 0;
        double ddeg = 0;
        int sides = 0;
    };
    Wall wall(std::size_t i, const Point &o, const Box &sq) const {
        Wall w;
        const Box &b = l_.boxes[i];
        const double minX = sq.minX + margin_, maxX = sq.maxX - margin_;
        const double minY = sq.minY + margin_, maxY = sq.maxY - margin_;
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
            ++w.sides;
        }
        if (b.maxX > maxX) {
            w.depth += b.maxX - maxX;
            w.dx += 1;
            w.ddeg -= (hiX.y - o.y) * kRadPerDeg;
            ++w.sides;
        }
        if (b.minY < minY) {
            w.depth += minY - b.minY;
            w.dy -= 1;
            w.ddeg -= (loY.x - o.x) * kRadPerDeg; // d(v.y)/dθ = v.x - o.x
            ++w.sides;
        }
        if (b.maxY > maxY) {
            w.depth += b.maxY - maxY;
            w.dy += 1;
            w.ddeg += (hiY.x - o.x) * kRadPerDeg;
            ++w.sides;
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

// L-BFGS with a backtracking line search on the energy. With `untilClear`, stops as soon as
// the layout is valid inside the square and returns true; otherwise runs until stuck or out
// of iterations. When stuck with `glsRounds` left, raises the weights of the overlapping
// pairs and carries on. Either way the layout holds the lowest-energy poses found.
bool minimise(Energy &e, long iterations, int glsRounds, bool untilClear,
              const std::atomic<bool> *stop) {
    constexpr std::size_t kMemory = 8;
    constexpr long kCheckEvery = 100; // how often to look for a stall worth a GLS round

    std::vector<double> z = e.get(), g, zNext, gNext, d;
    bool clear = false;
    double energy = e.eval(g, clear);
    if (untilClear && clear) return true;

    struct Pair {
        std::vector<double> s, y;
        double rho;
    };
    std::deque<Pair> memory;
    std::vector<double> alpha(kMemory);
    double checkpoint = energy;

    // Stuck: either raise the weights and go on, or give up.
    auto guide = [&]() {
        if (glsRounds <= 0 || !e.raiseWeights()) return false;
        --glsRounds;
        energy = e.eval(g, clear);
        checkpoint = energy;
        memory.clear();
        return true;
    };

    for (long it = 0; it < iterations; ++it) {
        if (stop && stop->load(std::memory_order_relaxed)) break;

        if (it > 0 && it % kCheckEvery == 0) {
            const bool stalled = energy > 0.99 * checkpoint;
            checkpoint = energy;
            if (stalled && untilClear) guide();
        }

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

        const double slope = dot(g, d);
        if (!(slope < 0)) {
            if (!memory.empty()) {
                memory.clear();
                continue;
            }
            if (guide()) continue;
            break; // flat: nothing left to follow
        }

        // Backtrack until the energy drops enough.
        double t = 1;
        bool moved = false;
        double eNext = 0;
        bool clearNext = false;
        zNext.resize(z.size());
        for (int tries = 0; tries < 40; ++tries, t /= 2) {
            for (std::size_t q = 0; q < z.size(); ++q) zNext[q] = z[q] + t * d[q];
            e.set(zNext);
            eNext = e.eval(gNext, clearNext);
            if (eNext <= energy + 1e-4 * t * slope || (untilClear && clearNext)) {
                moved = true;
                break;
            }
        }

        if (!moved) {
            e.set(z);
            if (!memory.empty()) {
                memory.clear();
                continue;
            }
            if (guide()) continue;
            return false; // stuck at a local minimum
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
        if (untilClear && clearNext) return true;
    }

    e.set(z);
    return false;
}

// True if tree i + n/2 is tree i turned 180° and reflected through (cx, cy), for every
// i < n/2: the layout SqueezeConfig::symmetric works on.
bool pointSymmetric(const Layout &l, double cx, double cy) {
    constexpr double kTol = 1e-6;
    const std::size_t pairs = l.size() / 2;
    if (pairs == 0) return false;
    for (std::size_t i = 0; i < pairs; ++i) {
        const Pose &p = l.poses[i];
        const Pose &q = l.poses[i + pairs];
        const double turn = geometry::normalizeDeg(q.deg - p.deg - 180);
        const bool matches = std::abs(p.x + q.x - 2 * cx) < kTol &&
                             std::abs(p.y + q.y - 2 * cy) < kTol &&
                             std::min(turn, 360 - turn) < kTol;
        if (!matches) return false;
    }
    return true;
}

// Moves every free tree at once down the slope of the overlap until the layout is valid
// inside the square, or it gets stuck. Returns whether the layout came out valid; either way
// `l` holds the lowest-energy poses found.
bool relax(Layout &l, const Box &square, const SqueezeConfig &cfg, bool symmetric) {
    Energy::Options o;
    o.margin = cfg.relaxMargin;
    o.movable = &cfg.movable;
    o.symmetric = symmetric;
    Energy e(l, square, o);
    return minimise(e, cfg.relaxIters, cfg.glsRounds, true, cfg.stop);
}
} // namespace

std::optional<Layout> squeeze(const Layout &layout, double target, const SqueezeConfig &cfg,
                              std::mt19937_64 &rng) {
    const std::size_t n = layout.size();
    if (n == 0 || cfg.iterations <= 0) return std::nullopt;

    const bool all = cfg.movable.empty();
    auto moves = [&](std::size_t i) { return all || cfg.movable[i]; };
    std::vector<std::size_t> free;
    for (std::size_t i = 0; i < n; ++i) {
        if (moves(i)) free.push_back(i);
    }
    if (free.empty()) return std::nullopt;

    // Pull the free trees toward the centre so their box fits the square; angles stay. Each
    // axis shrinks only as much as it needs, so a narrow layout gets no overlap across its
    // width. Trees that stay put must already fit.
    Layout l = layout;
    const Box b = l.bounds();
    const double cx = (b.minX + b.maxX) / 2;
    const double cy = (b.minY + b.maxY) / 2;
    const double fx = std::min(1.0, target / (b.maxX - b.minX));
    const double fy = std::min(1.0, target / (b.maxY - b.minY));
    for (const std::size_t i : free) {
        const Pose p = l.poses[i];
        l.setPose(i, {cx + (p.x - cx) * fx, cy + (p.y - cy) * fy, p.deg});
    }

    const Box square{cx - target / 2, cy - target / 2, cx + target / 2, cy + target / 2};

    // A valid layout inside the square, checked independently of the penalty.
    auto fits = [&] { return l.side() <= target && l.isValid() ? std::optional(l) : std::nullopt; };

    // Symmetric relax only on a layout that is symmetric about the square's centre.
    const bool symmetric = cfg.symmetric && pointSymmetric(l, cx, cy);
    if (cfg.relaxIters > 0 && relax(l, square, cfg, symmetric)) return fits();

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

    std::uniform_int_distribution<std::size_t> pickFree(0, free.size() - 1);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    auto jiggle = [&](double s) { return s * (2 * unit(rng) - 1); };

    // Each move's penetration with its neighbours, before and after: (j, depth).
    std::vector<std::pair<std::size_t, double>> before, after, scratch;
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
    auto costAt = [&](std::size_t i, const Pose &p) {
        l.setPose(i, p);
        const double out = outside(l.boxes[i], square);
        return contacts(i, scratch) + out * out;
    };

    constexpr int kLineTries = 8;
    for (long it = 0; it < cfg.iterations; ++it) {
        const bool checkStop = cfg.stop && it % 4096 == 0;
        if (checkStop && cfg.stop->load(std::memory_order_relaxed)) break;

        const double t = static_cast<double>(it) / cfg.iterations;
        const double temp = scale * cfg.startTemp * std::pow(cfg.endTemp / cfg.startTemp, t);
        const double step = cfg.startStep * std::pow(cfg.endStep / cfg.startStep, t);

        std::size_t i = unit(rng) < cfg.badShare
                            ? bad[std::uniform_int_distribution<std::size_t>(0, bad.size() - 1)(rng)]
                            : free[pickFree(rng)];
        if (!moves(i)) i = free[pickFree(rng)];

        const Pose old = l.poses[i];
        const double outBefore = outside(l.boxes[i], square);
        const double pairsBefore = contacts(i, before);

        Pose moved = old;
        if (unit(rng) < cfg.lineShare) {
            // Several distances along one direction; take the cheapest.
            const int axis = static_cast<int>(3 * unit(rng));
            const double span = axis == 2 ? 400 * step : 4 * step;
            double cheapest = INFINITY;
            for (int k = 0; k < kLineTries; ++k) {
                Pose p = old;
                const double shift = jiggle(span);
                (axis == 0 ? p.x : axis == 1 ? p.y : p.deg) += shift;
                const double cost = costAt(i, p);
                if (cost < cheapest) {
                    cheapest = cost;
                    moved = p;
                }
            }
        } else {
            const double kind = unit(rng);
            if (kind < 0.6) {
                moved.x += jiggle(step);
                moved.y += jiggle(step);
            } else if (kind < 0.98) {
                moved.deg += jiggle(100 * step);
            } else {
                moved.deg += 180;
            }
        }
        moved.deg = geometry::normalizeDeg(moved.deg);

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

    // Single-tree moves left overlap: one more relax from where they got to. The single-tree
    // moves broke any symmetry, so this one runs free.
    if (!bad.empty() && !(cfg.relaxIters > 0 && relax(l, square, cfg, false))) return std::nullopt;

    return fits();
}

double compress(Layout &layout, const CompressConfig &cfg) {
    if (layout.size() == 0) return 0;

    const Box b = layout.bounds();
    const double cx = (b.minX + b.maxX) / 2;
    const double cy = (b.minY + b.maxY) / 2;
    const double side = layout.side();

    Energy::Options o;
    o.margin = cfg.margin;
    o.movable = &cfg.movable;
    o.symmetric = cfg.symmetric && pointSymmetric(layout, cx, cy);
    o.sideIsVariable = true;
    Energy e(layout, {cx - side / 2, cy - side / 2, cx + side / 2, cy + side / 2}, o);

    for (int r = 0; r < cfg.rounds; ++r) {
        if (cfg.stop && cfg.stop->load()) break;
        const double f = cfg.rounds > 1 ? static_cast<double>(r) / (cfg.rounds - 1) : 1.0;
        e.setWeight(cfg.startWeight * std::pow(cfg.endWeight / cfg.startWeight, f));
        minimise(e, cfg.iterations, 0, false, cfg.stop);
    }

    return e.side();
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

std::vector<bool> borderTrees(const Layout &layout, double width) {
    std::vector<bool> near(layout.size(), false);
    if (layout.size() == 0) return near;

    const Box b = layout.bounds();
    const double half = layout.side() / 2;
    const double cx = (b.minX + b.maxX) / 2;
    const double cy = (b.minY + b.maxY) / 2;
    for (std::size_t i = 0; i < layout.size(); ++i) {
        const Box &t = layout.boxes[i];
        near[i] = t.minX < cx - half + width || t.maxX > cx + half - width ||
                  t.minY < cy - half + width || t.maxY > cy + half - width;
    }

    return near;
}
} // namespace packer::packing
