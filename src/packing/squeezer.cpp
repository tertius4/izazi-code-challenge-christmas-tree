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
using geometry::centre;
using geometry::kInf;
using geometry::outside;
using geometry::Point;
using geometry::squareAround;

// The trees that still overlap or stick out, kept so that any of them can be picked, added or
// dropped in constant time.
class BadSet {
  public:
    explicit BadSet(std::size_t n) : position_(n, kNone) {
    }

    void set(std::size_t i, bool bad) {
        const bool listed = position_[i] != kNone;
        if (bad && !listed) {
            position_[i] = trees_.size();
            trees_.push_back(i);
        } else if (!bad && listed) {
            const std::size_t last = trees_.back();
            trees_[position_[i]] = last;
            position_[last] = position_[i];
            trees_.pop_back();
            position_[i] = kNone;
        }
    }

    bool empty() const {
        return trees_.empty();
    }
    std::size_t size() const {
        return trees_.size();
    }
    std::size_t operator[](std::size_t k) const {
        return trees_[k];
    }

  private:
    static constexpr std::size_t kNone = static_cast<std::size_t>(-1);
    std::vector<std::size_t> trees_;    // the bad trees, in no order
    std::vector<std::size_t> position_; // position_[i]: where tree i is in trees_, or kNone
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

    Energy(Layout &layout, const Box &square, const Options &options)
        : layout_(layout), margin_(options.margin), sideIsVariable_(options.sideIsVariable),
          cx_(centre(square).x), cy_(centre(square).y), halfSide_((square.maxX - square.minX) / 2) {
        const std::size_t n = layout.size();
        block_.assign(n, -1);
        mirrored_.assign(n, false);
        const bool all = options.movable == nullptr || options.movable->empty();
        auto moves = [&](std::size_t i) { return all || (*options.movable)[i]; };

        const std::size_t pairs = options.symmetric ? n / 2 : 0;
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
        return 2 * halfSide_;
    }

    Box square() const {
        return squareAround({cx_, cy_}, 2 * halfSide_);
    }

    void setWeight(double weight) {
        weight_ = weight;
    }

    std::vector<double> get() const {
        std::vector<double> vars(size());
        for (std::size_t block = 0; block < primary_.size(); ++block) {
            const Pose &pose = layout_.poses[primary_[block]];
            vars[3 * block] = pose.x;
            vars[3 * block + 1] = pose.y;
            vars[3 * block + 2] = pose.deg / kDegPerUnit;
        }
        if (sideIsVariable_) vars.back() = 2 * halfSide_;
        return vars;
    }

    void set(const std::vector<double> &vars) {
        if (sideIsVariable_) halfSide_ = vars.back() / 2;
        for (std::size_t i = 0; i < layout_.size(); ++i) {
            const int block = block_[i];
            if (block < 0) continue;
            const double x = vars[3 * block], y = vars[3 * block + 1];
            const double deg = vars[3 * block + 2] * kDegPerUnit;
            layout_.setPose(i, mirrored_[i] ? Pose{2 * cx_ - x, 2 * cy_ - y, deg + 180}
                                            : Pose{x, y, deg});
        }
    }

    // Energy and its gradient at the current poses; `clear` says whether the layout is already
    // valid and inside the square (by collide(), not the margin).
    double eval(std::vector<double> &grad, bool &clear) const {
        const std::size_t n = layout_.size();
        const Box walls = square();
        treeGrad_.assign(3 * n, 0.0);
        clear = true;
        double penalty = 0;
        double sideSlope = 0; // d(penalty)/d(side)

        for (std::size_t i = 0; i < n; ++i) {
            const Point origin{layout_.poses[i].x, layout_.poses[i].y};
            double *gradI = &treeGrad_[3 * i];

            if (outside(layout_.boxes[i], walls) > 0) clear = false;
            const Wall wall = wallOverlap(i, origin, walls);
            if (wall.depth > 0) {
                penalty += wall.depth * wall.depth;
                gradI[0] += 2 * wall.depth * wall.dx;
                gradI[1] += 2 * wall.depth * wall.dy;
                gradI[2] += 2 * wall.depth * wall.ddeg * kDegPerUnit;
                sideSlope += 2 * wall.depth * (-0.5 * wall.sides); // each wall moves out by side/2
            }

            layout_.anyNeighbour(i, [&](std::size_t j) {
                if (block_[i] < 0 && block_[j] < 0) return false; // neither moves
                const tree::Overlap pair =
                    tree::overlap(layout_.placed[i], origin, layout_.placed[j], margin_);
                if (pair.depth <= 0) return false;
                const double weight = pairWeight(i, j);
                if (j > i) penalty += weight * pair.depth * pair.depth;
                // d(depth²)/d(pose of i): each pair is visited once from each side.
                gradI[0] += 2 * weight * pair.depth * pair.dx;
                gradI[1] += 2 * weight * pair.depth * pair.dy;
                gradI[2] += 2 * weight * pair.depth * pair.ddeg * kDegPerUnit;
                if (clear && tree::collide(layout_.placed[i], layout_.placed[j])) clear = false;
                return false;
            });
        }

        // Fold each tree's gradient into its block: a mirrored tree moves against its block's x
        // and y, and turns with it.
        const double scale = sideIsVariable_ ? weight_ : 1.0;
        grad.assign(size(), 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            const int block = block_[i];
            if (block < 0) continue;
            const double sign = mirrored_[i] ? -1.0 : 1.0;
            grad[3 * block] += scale * sign * treeGrad_[3 * i];
            grad[3 * block + 1] += scale * sign * treeGrad_[3 * i + 1];
            grad[3 * block + 2] += scale * treeGrad_[3 * i + 2];
        }

        if (!sideIsVariable_) return penalty;
        grad.back() = 1 + weight_ * sideSlope;
        return 2 * halfSide_ + weight_ * penalty;
    }

    // Guided local search: raise the weight of every pair that still overlaps, the deepest by
    // 1 and the others in proportion. Returns false if nothing overlaps.
    bool raiseWeights() {
        std::vector<std::pair<std::uint64_t, double>> hits; // (pair key, depth)
        double deepest = 0;
        for (std::size_t i = 0; i < layout_.size(); ++i) {
            const Point origin{layout_.poses[i].x, layout_.poses[i].y};
            layout_.anyNeighbour(i, [&](std::size_t j) {
                if (j < i) return false;
                const double depth =
                    tree::overlap(layout_.placed[i], origin, layout_.placed[j], margin_).depth;
                if (depth > 0) {
                    hits.push_back({pairKey(i, j), depth});
                    deepest = std::max(deepest, depth);
                }
                return false;
            });
        }

        for (const auto &[key, depth] : hits) {
            auto it = pairWeights_.try_emplace(key, 1.0).first;
            it->second += depth / deepest;
        }

        return !hits.empty();
    }

  private:
    Layout &layout_;
    const double margin_;
    const bool sideIsVariable_;
    const double cx_, cy_; // the square's centre
    double halfSide_;
    double weight_ = 1;
    std::vector<int> block_;           // block_[i]: tree i's block of variables, or -1
    std::vector<bool> mirrored_;       // tree i mirrors its block instead of placing it as is
    std::vector<std::size_t> primary_; // the tree each block places as it is
    std::unordered_map<std::uint64_t, double> pairWeights_; // raised pair weights; the rest are 1
    mutable std::vector<double> treeGrad_; // per tree, the gradient by its own x, y and angle

    std::uint64_t pairKey(std::size_t i, std::size_t j) const {
        if (j < i) std::swap(i, j);
        return static_cast<std::uint64_t>(i) * layout_.size() + j;
    }

    double pairWeight(std::size_t i, std::size_t j) const {
        if (pairWeights_.empty()) return 1;
        const auto it = pairWeights_.find(pairKey(i, j));
        return it == pairWeights_.end() ? 1 : it->second;
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
    Wall wallOverlap(std::size_t i, const Point &origin, const Box &walls) const {
        Wall wall;
        const Box &box = layout_.boxes[i];
        const double minX = walls.minX + margin_, maxX = walls.maxX - margin_;
        const double minY = walls.minY + margin_, maxY = walls.maxY - margin_;
        if (box.minX >= minX && box.maxX <= maxX && box.minY >= minY && box.maxY <= maxY) {
            return wall;
        }

        // The vertex that reaches furthest each way; turning moves it by perp(vertex - origin).
        constexpr double kRadPerDeg = std::numbers::pi / 180.0;
        Point least{kInf, kInf}, most{-kInf, -kInf};
        Point leftmost{}, rightmost{}, lowest{}, highest{};
        for (const tree::Part &part : layout_.placed[i].parts) {
            for (std::size_t k = 0; k < part.count; ++k) {
                const Point &vertex = part.vertex[k];
                if (vertex.x < least.x) least.x = vertex.x, leftmost = vertex;
                if (vertex.x > most.x) most.x = vertex.x, rightmost = vertex;
                if (vertex.y < least.y) least.y = vertex.y, lowest = vertex;
                if (vertex.y > most.y) most.y = vertex.y, highest = vertex;
            }
        }

        if (box.minX < minX) {
            wall.depth += minX - box.minX;
            wall.dx -= 1;
            wall.ddeg += (leftmost.y - origin.y) * kRadPerDeg; // d(vertex.x)/dθ = -(y - origin.y)
            ++wall.sides;
        }
        if (box.maxX > maxX) {
            wall.depth += box.maxX - maxX;
            wall.dx += 1;
            wall.ddeg -= (rightmost.y - origin.y) * kRadPerDeg;
            ++wall.sides;
        }
        if (box.minY < minY) {
            wall.depth += minY - box.minY;
            wall.dy -= 1;
            wall.ddeg -= (lowest.x - origin.x) * kRadPerDeg; // d(vertex.y)/dθ = x - origin.x
            ++wall.sides;
        }
        if (box.maxY > maxY) {
            wall.depth += box.maxY - maxY;
            wall.dy += 1;
            wall.ddeg += (highest.x - origin.x) * kRadPerDeg;
            ++wall.sides;
        }

        return wall;
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
bool minimise(Energy &energy, long iterations, int glsRounds, bool untilClear,
              const std::atomic<bool> *stop) {
    constexpr std::size_t kMemory = 8;
    constexpr long kCheckEvery = 100; // how often to look for a stall worth a GLS round

    std::vector<double> vars = energy.get(), grad, nextVars, nextGrad, dir;
    bool clear = false;
    double value = energy.eval(grad, clear);
    if (untilClear && clear) return true;

    // One remembered step: how the variables and the gradient changed.
    struct Update {
        std::vector<double> dVars, dGrad;
        double invDot; // 1 / (dVars · dGrad)
    };
    std::deque<Update> memory;
    std::vector<double> coeff(kMemory);
    double checkpoint = value;

    // Stuck: either raise the weights and go on, or give up.
    auto guide = [&]() {
        if (glsRounds <= 0 || !energy.raiseWeights()) return false;
        --glsRounds;
        value = energy.eval(grad, clear);
        checkpoint = value;
        memory.clear();
        return true;
    };

    for (long it = 0; it < iterations; ++it) {
        if (stop && stop->load(std::memory_order_relaxed)) break;

        if (it > 0 && it % kCheckEvery == 0) {
            const bool stalled = value > 0.99 * checkpoint;
            checkpoint = value;
            if (stalled && untilClear) guide();
        }

        // Direction: -H·grad by the two-loop recursion.
        dir = grad;
        for (std::size_t k = memory.size(); k-- > 0;) {
            coeff[k] = memory[k].invDot * dot(memory[k].dVars, dir);
            for (std::size_t i = 0; i < dir.size(); ++i)
                dir[i] -= coeff[k] * memory[k].dGrad[i];
        }
        double scale = 0;
        if (memory.empty()) {
            // No curvature known yet: a first step of at most 1e-3 along any coordinate.
            double most = 0;
            for (const double component : grad)
                most = std::max(most, std::abs(component));
            scale = most > 0 ? 1e-3 / most : 0;
        } else {
            const Update &last = memory.back();
            scale = dot(last.dVars, last.dGrad) / dot(last.dGrad, last.dGrad);
        }
        for (double &component : dir)
            component *= scale;
        for (std::size_t k = 0; k < memory.size(); ++k) {
            const double correction = memory[k].invDot * dot(memory[k].dGrad, dir);
            for (std::size_t i = 0; i < dir.size(); ++i) {
                dir[i] += (coeff[k] - correction) * memory[k].dVars[i];
            }
        }
        for (double &component : dir)
            component = -component;

        // The energy is nearly flat in places, which can make the curvature estimate huge:
        // never move anything further than a small step at once, or trees fly off (and the
        // layout's grid grows to cover them).
        constexpr double kMaxStep = 0.02;
        double longest = 0;
        for (const double component : dir)
            longest = std::max(longest, std::abs(component));
        if (longest > kMaxStep) {
            for (double &component : dir)
                component *= kMaxStep / longest;
        }

        const double slope = dot(grad, dir);
        if (!(slope < 0)) {
            if (!memory.empty()) {
                memory.clear();
                continue;
            }
            if (guide()) continue;
            break; // flat: nothing left to follow
        }

        // Backtrack until the energy drops enough.
        double stepSize = 1;
        bool moved = false;
        double nextValue = 0;
        bool nextClear = false;
        nextVars.resize(vars.size());
        for (int tries = 0; tries < 40; ++tries, stepSize /= 2) {
            for (std::size_t i = 0; i < vars.size(); ++i)
                nextVars[i] = vars[i] + stepSize * dir[i];
            energy.set(nextVars);
            nextValue = energy.eval(nextGrad, nextClear);
            if (nextValue <= value + 1e-4 * stepSize * slope || (untilClear && nextClear)) {
                moved = true;
                break;
            }
        }

        if (!moved) {
            energy.set(vars);
            if (!memory.empty()) {
                memory.clear();
                continue;
            }
            if (guide()) continue;
            return false; // stuck at a local minimum
        }

        Update update{std::vector<double>(vars.size()), std::vector<double>(vars.size()), 0};
        for (std::size_t i = 0; i < vars.size(); ++i) {
            update.dVars[i] = nextVars[i] - vars[i];
            update.dGrad[i] = nextGrad[i] - grad[i];
        }
        const double curvature = dot(update.dVars, update.dGrad);
        if (curvature > 1e-30) {
            update.invDot = 1 / curvature;
            memory.push_back(std::move(update));
            if (memory.size() > kMemory) memory.pop_front();
        }

        std::swap(vars, nextVars);
        std::swap(grad, nextGrad);
        value = nextValue;
        if (untilClear && nextClear) return true;
    }

    energy.set(vars);
    return false;
}

// True if tree i + n/2 is tree i turned 180° and reflected through (cx, cy), for every
// i < n/2: the layout SqueezeConfig::symmetric works on.
bool pointSymmetric(const Layout &layout, double cx, double cy) {
    constexpr double kTol = 1e-6;
    const std::size_t pairs = layout.size() / 2;
    if (pairs == 0) return false;
    for (std::size_t i = 0; i < pairs; ++i) {
        const Pose &pose = layout.poses[i];
        const Pose &mirror = layout.poses[i + pairs];
        const double turn = geometry::normalizeDeg(mirror.deg - pose.deg - 180);
        const bool matches = std::abs(pose.x + mirror.x - 2 * cx) < kTol &&
                             std::abs(pose.y + mirror.y - 2 * cy) < kTol &&
                             std::min(turn, 360 - turn) < kTol;
        if (!matches) return false;
    }
    return true;
}

// Moves every free tree at once down the slope of the overlap until the layout is valid
// inside the square, or it gets stuck. Returns whether the layout came out valid; either way
// `layout` holds the lowest-energy poses found.
bool relax(Layout &layout, const Box &square, const SqueezeConfig &cfg, bool symmetric) {
    Energy::Options options;
    options.margin = cfg.relaxMargin;
    options.movable = &cfg.movable;
    options.symmetric = symmetric;
    Energy energy(layout, square, options);
    return minimise(energy, cfg.relaxIters, cfg.glsRounds, true, cfg.stop);
}
} // namespace

std::optional<Layout> squeeze(const Layout &layout, double target, const SqueezeConfig &cfg,
                              std::mt19937_64 &rng) {
    const std::size_t n = layout.size();
    if (n == 0 || cfg.iterations <= 0) return std::nullopt;

    const bool all = cfg.movable.empty();
    auto moves = [&](std::size_t i) { return all || cfg.movable[i]; };
    std::vector<std::size_t> freeTrees;
    for (std::size_t i = 0; i < n; ++i) {
        if (moves(i)) freeTrees.push_back(i);
    }
    if (freeTrees.empty()) return std::nullopt;

    // Pull the free trees toward the centre so their box fits the square; angles stay. Each
    // axis shrinks only as much as it needs, so a narrow layout gets no overlap across its
    // width. Trees that stay put must already fit.
    Layout current = layout;
    const Box bounds = current.bounds();
    const auto [cx, cy] = centre(bounds);
    const double scaleX = std::min(1.0, target / (bounds.maxX - bounds.minX));
    const double scaleY = std::min(1.0, target / (bounds.maxY - bounds.minY));
    for (const std::size_t i : freeTrees) {
        const Pose pose = current.poses[i];
        current.setPose(i, {cx + (pose.x - cx) * scaleX, cy + (pose.y - cy) * scaleY, pose.deg});
    }

    const Box square = squareAround({cx, cy}, target);

    // A valid layout inside the square, checked independently of the penalty.
    auto fits = [&] {
        return current.side() <= target && current.isValid() ? std::optional(current)
                                                             : std::nullopt;
    };

    // Symmetric relax only on a layout that is symmetric about the square's centre.
    const bool symmetric = cfg.symmetric && pointSymmetric(current, cx, cy);
    if (cfg.relaxIters > 0 && relax(current, square, cfg, symmetric)) return fits();

    // Penalty: the sum of every overlapping pair's penetration squared, plus how far each tree
    // sticks out squared. Each tree counts its own violations, so a tree is clear at count 0.
    std::vector<int> violations(n, 0);
    BadSet bad(n);
    double energy = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const double stickOut = outside(current.boxes[i], square);
        energy += stickOut * stickOut;
        violations[i] += stickOut > 0;
        current.anyNeighbour(i, [&](std::size_t j) {
            const double depth = tree::penetration(current.placed[i], current.placed[j]);
            if (j > i) energy += depth * depth;
            violations[i] += depth > 0;
            return false;
        });
        bad.set(i, violations[i] > 0);
    }

    if (bad.empty()) return fits();

    // Temperatures scale with the starting penalty per bad tree.
    const double tempScale = energy / static_cast<double>(bad.size());

    std::uniform_int_distribution<std::size_t> pickFree(0, freeTrees.size() - 1);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    auto jiggle = [&](double size) { return size * (2 * unit(rng) - 1); };

    // Each move's penetration with its neighbours, before and after: (j, depth).
    std::vector<std::pair<std::size_t, double>> before, after, scratch;
    auto contacts = [&](std::size_t i, std::vector<std::pair<std::size_t, double>> &into) {
        into.clear();
        double sum = 0;
        current.anyNeighbour(i, [&](std::size_t j) {
            const double depth = tree::penetration(current.placed[i], current.placed[j]);
            if (depth > 0) {
                into.push_back({j, depth});
                sum += depth * depth;
            }
            return false;
        });
        return sum;
    };
    auto costAt = [&](std::size_t i, const Pose &pose) {
        current.setPose(i, pose);
        const double stickOut = outside(current.boxes[i], square);
        return contacts(i, scratch) + stickOut * stickOut;
    };

    constexpr int kLineTries = 8;
    for (long it = 0; it < cfg.iterations; ++it) {
        const bool checkStop = cfg.stop && it % 4096 == 0;
        if (checkStop && cfg.stop->load(std::memory_order_relaxed)) break;

        const double progress = static_cast<double>(it) / cfg.iterations;
        const double temp =
            tempScale * cfg.startTemp * std::pow(cfg.endTemp / cfg.startTemp, progress);
        const double step = cfg.startStep * std::pow(cfg.endStep / cfg.startStep, progress);

        std::size_t i =
            unit(rng) < cfg.badShare
                ? bad[std::uniform_int_distribution<std::size_t>(0, bad.size() - 1)(rng)]
                : freeTrees[pickFree(rng)];
        if (!moves(i)) i = freeTrees[pickFree(rng)];

        const Pose oldPose = current.poses[i];
        const double outBefore = outside(current.boxes[i], square);
        const double pairsBefore = contacts(i, before);

        Pose moved = oldPose;
        if (unit(rng) < cfg.lineShare) {
            // Several distances along one direction; take the cheapest.
            const int axis = static_cast<int>(3 * unit(rng));
            const double span = axis == 2 ? 400 * step : 4 * step;
            double cheapest = kInf;
            for (int k = 0; k < kLineTries; ++k) {
                Pose trial = oldPose;
                const double shift = jiggle(span);
                (axis == 0 ? trial.x : axis == 1 ? trial.y : trial.deg) += shift;
                const double cost = costAt(i, trial);
                if (cost < cheapest) {
                    cheapest = cost;
                    moved = trial;
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

        current.setPose(i, moved);
        const double outAfter = outside(current.boxes[i], square);
        const double pairsAfter = contacts(i, after);

        const double worseBy =
            pairsAfter - pairsBefore + outAfter * outAfter - outBefore * outBefore;
        const bool accept = worseBy <= 0 || unit(rng) < std::exp(-worseBy / temp);
        if (!accept) {
            current.setPose(i, oldPose);
            continue;
        }

        energy += worseBy;
        // Recount the violations of tree i and of every neighbour it touched before or after.
        for (const auto &[j, depth] : before) {
            violations[j] -= 1;
        }
        for (const auto &[j, depth] : after) {
            violations[j] += 1;
        }
        for (const auto &[j, depth] : before) {
            bad.set(j, violations[j] > 0);
        }
        for (const auto &[j, depth] : after) {
            bad.set(j, violations[j] > 0);
        }
        violations[i] = static_cast<int>(after.size()) + (outAfter > 0);
        bad.set(i, violations[i] > 0);

        if (bad.empty()) break;
    }

    // Single-tree moves left overlap: one more relax from where they got to. The single-tree
    // moves broke any symmetry, so this one runs free.
    if (!bad.empty() && !(cfg.relaxIters > 0 && relax(current, square, cfg, false))) {
        return std::nullopt;
    }

    return fits();
}

double compress(Layout &layout, const CompressConfig &cfg) {
    if (layout.size() == 0) return 0;

    const auto [cx, cy] = centre(layout.bounds());
    const double side = layout.side();

    Energy::Options options;
    options.margin = cfg.margin;
    options.movable = &cfg.movable;
    options.symmetric = cfg.symmetric && pointSymmetric(layout, cx, cy);
    options.sideIsVariable = true;
    Energy energy(layout, squareAround({cx, cy}, side), options);

    for (int round = 0; round < cfg.rounds; ++round) {
        if (cfg.stop && cfg.stop->load()) break;
        const double progress =
            cfg.rounds > 1 ? static_cast<double>(round) / (cfg.rounds - 1) : 1.0;
        energy.setWeight(cfg.startWeight * std::pow(cfg.endWeight / cfg.startWeight, progress));
        minimise(energy, cfg.iterations, 0, false, cfg.stop);
    }

    return energy.side();
}

Layout shake(const Layout &layout, int kicks, std::mt19937_64 &rng) {
    Layout shaken = layout;
    const std::size_t n = shaken.size();
    if (n < 2) return shaken;

    std::uniform_int_distribution<std::size_t> pickTree(0, n - 1);
    std::uniform_real_distribution<double> unit(0.0, 1.0);

    for (int k = 0; k < kicks; ++k) {
        const std::size_t i = pickTree(rng);
        const Pose poseI = shaken.poses[i];
        const double kind = unit(rng);

        if (kind < 0.5) {
            // Swapping two trees turned the same way changes nothing: only differently turned
            // neighbours are worth it.
            std::vector<std::size_t> others;
            shaken.anyNeighbour(i, [&](std::size_t j) {
                const double turn = geometry::normalizeDeg(shaken.poses[j].deg - poseI.deg);
                if (std::min(turn, 360 - turn) > 20) others.push_back(j);
                return false;
            });

            if (!others.empty()) {
                const std::size_t j =
                    others[std::uniform_int_distribution<std::size_t>(0, others.size() - 1)(rng)];
                const Pose poseJ = shaken.poses[j];
                shaken.setPose(i, {poseJ.x, poseJ.y, poseI.deg});
                shaken.setPose(j, {poseI.x, poseI.y, poseJ.deg});
                continue;
            }
        }

        if (kind < 0.8) {
            shaken.setPose(i, {poseI.x, poseI.y, geometry::normalizeDeg(poseI.deg + 180)});
        } else {
            const Box bounds = shaken.bounds();
            shaken.setPose(i, {bounds.minX + unit(rng) * (bounds.maxX - bounds.minX),
                               bounds.minY + unit(rng) * (bounds.maxY - bounds.minY),
                               360 * unit(rng)});
        }
    }

    return shaken;
}

std::vector<bool> borderTrees(const Layout &layout, double width) {
    std::vector<bool> near(layout.size(), false);
    if (layout.size() == 0) return near;

    const Box square = squareAround(centre(layout.bounds()), layout.side());
    for (std::size_t i = 0; i < layout.size(); ++i) {
        const Box &box = layout.boxes[i];
        near[i] = box.minX < square.minX + width || box.maxX > square.maxX - width ||
                  box.minY < square.minY + width || box.maxY > square.maxY - width;
    }

    return near;
}
} // namespace packer::packing
