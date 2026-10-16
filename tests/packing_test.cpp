#include "packing/annealer.h"
#include "packing/crossover.h"
#include "packing/crystal.h"
#include "packing/improver.h"
#include "packing/lattice.h"
#include "packing/nfp.h"
#include "packing/rotation.h"
#include "packing/solver.h"
#include "packing/squeezer.h"
#include "packing/starts.h"

#include <gtest/gtest.h>

using namespace packer::packing;

TEST(Layout, DetectsOverlapAndSide) {
    Layout l;
    l.add({0, 0, 0});
    EXPECT_NEAR(l.side(), 1.0, 1e-12); // upright tree: 0.7 wide, 1.0 tall
    l.add({0.1, 0, 0});
    EXPECT_FALSE(l.isValid());
    l.setPose(1, {2, 0, 0});
    EXPECT_TRUE(l.isValid());
}

TEST(Layout, GridAgreesWithBruteForce) {
    // Random moves, some far outside the grid, then compare each tree's collision check
    // through the grid with a plain check against every other tree.
    std::mt19937_64 rng(4);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    Layout l;
    for (int i = 0; i < 60; ++i)
        l.add({4 * u(rng), 4 * u(rng), 360 * u(rng)});
    for (int step = 0; step < 3000; ++step) {
        const std::size_t i = static_cast<std::size_t>(u(rng) * l.size());
        const double range = step % 100 == 0 ? 12 : 4;
        l.setPose(
            i, {range * u(rng) - (range - 4) / 2, range * u(rng) - (range - 4) / 2, 360 * u(rng)});
        if (step % 500 == 0) l.remove(static_cast<std::size_t>(u(rng) * l.size()));
        for (std::size_t a = 0; a < l.size(); ++a) {
            bool expected = false;
            for (std::size_t b = 0; b < l.size(); ++b) {
                expected = expected || (a != b && packer::tree::collide(l.placed[a], l.placed[b]));
            }
            ASSERT_EQ(l.collidesWithOthers(a), expected) << "step " << step << " tree " << a;
        }
    }
}

TEST(Layout, RemovalOrderStartsWithBestTreeToRemove) {
    std::mt19937_64 rng(7);
    const Layout l = randomStart(15, rng);
    const std::vector<std::size_t> order = l.removalOrder();
    ASSERT_EQ(order.size(), l.size());
    EXPECT_EQ(order.front(), l.bestTreeToRemove());
    double last = 0;
    for (const std::size_t i : order) {
        Layout without = l;
        without.remove(i);
        EXPECT_GE(without.side(), last - 1e-12);
        last = without.side();
    }
}

TEST(Squeezer, FitsALooseLayoutIntoASmallerSquare) {
    std::mt19937_64 rng(8);
    const Layout loose = randomStart(10, rng);
    const double target = 0.9 * loose.side();
    SqueezeConfig cfg;
    cfg.iterations = 200000;
    const std::optional<Layout> l = squeeze(loose, target, cfg, rng);
    ASSERT_TRUE(l.has_value());
    EXPECT_EQ(l->size(), 10u);
    EXPECT_TRUE(l->isValid());
    EXPECT_LE(l->side(), target);
}

TEST(Squeezer, GivesUpOnAnImpossibleSquare) {
    std::mt19937_64 rng(9);
    SqueezeConfig cfg;
    cfg.iterations = 20000;
    Layout two;
    two.add({0, 0, 0});
    two.add({2, 0, 0});
    // Two trees cannot share a square that only just fits one.
    EXPECT_FALSE(squeeze(two, bestSingle().side(), cfg, rng).has_value());
}

TEST(Squeezer, RelaxFitsALooseLayoutIntoASmallerSquare) {
    std::mt19937_64 rng(10);
    const Layout loose = randomStart(30, rng);
    const double target = 0.95 * loose.side();
    SqueezeConfig cfg;
    cfg.iterations = 1000;
    cfg.relaxIters = 2000;
    const std::optional<Layout> l = squeeze(loose, target, cfg, rng);
    ASSERT_TRUE(l.has_value());
    EXPECT_EQ(l->size(), 30u);
    EXPECT_TRUE(l->isValid());
    EXPECT_LE(l->side(), target);
}

TEST(Squeezer, ShakenLayoutSettlesBackIntoItsSquare) {
    std::mt19937_64 rng(11);
    const Layout l = randomStart(30, rng);
    SqueezeConfig cfg;
    cfg.iterations = 30000;
    cfg.relaxIters = 2000;
    int settled = 0;
    for (int k = 0; k < 10; ++k) {
        const Layout kicked = shake(l, 3, rng);
        ASSERT_EQ(kicked.size(), l.size());
        const std::optional<Layout> back = squeeze(kicked, l.side(), cfg, rng);
        if (!back) continue;
        ++settled;
        EXPECT_TRUE(back->isValid());
        EXPECT_LE(back->side(), l.side());
    }
    EXPECT_GT(settled, 5); // a loose layout has room for the kicked trees
}

TEST(Starts, RandomStartIsValid) {
    std::mt19937_64 rng(1);
    const Layout l = randomStart(20, rng);
    EXPECT_EQ(l.size(), 20u);
    EXPECT_TRUE(l.isValid());
}

TEST(Starts, BestSingleBeatsUpright) {
    EXPECT_LE(bestSingle().side(), 1.0);
}

// One short lattice search shared by the tests below.
const LatticeParams &testLattice() {
    static const LatticeParams lattice = [] {
        std::mt19937_64 rng(1);
        return findLattice(10000, rng).value();
    }();
    return lattice;
}

TEST(Lattice, PeriodicCheckAgreesWithBruteForce) {
    const LatticeParams &p = testLattice();
    EXPECT_TRUE(latticeIsValid(p));

    // Build a 10x10-cell patch and check every pair the slow way.
    Layout patch;
    for (int i = 0; i < 10; ++i) {
        for (int j = 0; j < 10; ++j) {
            patch.add(latticePose(p, i, j, false));
            patch.add(latticePose(p, i, j, true));
        }
    }
    EXPECT_TRUE(patch.isValid());

    // Squashing the cell must make it collide, in both checks.
    LatticeParams squashed = p;
    squashed.a *= 0.9;
    squashed.b *= 0.9;
    EXPECT_FALSE(latticeIsValid(squashed));
}

TEST(Lattice, DenserThanOldColumnPattern) {
    // The old 0° column lattice used about 0.41 per tree once cut into a square for N=200.
    EXPECT_LT(areaPerTree(testLattice()), 0.36);
}

TEST(Lattice, CutHoldsExactlyNValidTrees) {
    for (int n : {1, 2, 17, 50, 200}) {
        const Layout l = cutLattice(testLattice(), n);
        EXPECT_EQ(static_cast<int>(l.size()), n);
        EXPECT_TRUE(l.isValid()) << "n=" << n;
    }
}

TEST(Lattice, CutBeatsOldResultForLargeN) {
    // The old column lattice gave side 9.1 for N=200.
    EXPECT_LT(cutLattice(testLattice(), 200).side(), 9.1);
}

TEST(Crystal, DoubleLatticeIsTheSamePattern) {
    const Crystal c = toCrystal(testLattice());
    EXPECT_TRUE(crystalIsValid(c));
    EXPECT_NEAR(areaPerTree(c), areaPerTree(testLattice()), 1e-12);
    for (int i = -2; i <= 2; ++i) {
        for (int j = -2; j <= 2; ++j) {
            for (std::size_t k = 0; k < 2; ++k) {
                const Pose a = crystalPose(c, i, j, k);
                const Pose b = latticePose(testLattice(), i, j, k == 1);
                EXPECT_NEAR(a.x, b.x, 1e-9);
                EXPECT_NEAR(a.y, b.y, 1e-9);
            }
        }
    }
}

TEST(Crystal, PeriodicCheckAgreesWithBruteForce) {
    std::mt19937_64 rng(3);
    const Crystal c = findCrystal(3, false, 3000, rng).value();
    ASSERT_TRUE(crystalIsValid(c));
    for (const Crystal &t : {c, turned(c, 37), doubled(c)}) {
        EXPECT_TRUE(crystalIsValid(t));
        Layout patch;
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                for (std::size_t k = 0; k < t.cell.size(); ++k) patch.add(crystalPose(t, i, j, k));
            }
        }
        EXPECT_TRUE(patch.isValid());
    }
    EXPECT_NEAR(areaPerTree(doubled(c)), areaPerTree(c), 1e-12);
    EXPECT_FALSE(crystalIsValid({{c.u.x * 0.9, c.u.y * 0.9}, {c.v.x * 0.9, c.v.y * 0.9}, c.cell}));
}

TEST(Crystal, CutsHoldExactlyNValidTreesAtAnyAngle) {
    const Crystal c = toCrystal(testLattice());
    std::mt19937_64 rng(2);
    for (int n : {1, 7, 40, 150}) {
        for (double deg : {0.0, 13.0, 61.5}) {
            const Layout l = cutCrystal(turned(c, deg), n);
            EXPECT_EQ(static_cast<int>(l.size()), n);
            EXPECT_TRUE(l.isValid()) << "n=" << n << " deg=" << deg;
        }
        const auto [best, deg] = bestCrystalCut(c, n, 8, rng);
        EXPECT_TRUE(best.isValid());
        EXPECT_LE(best.side(), cutCrystal(turned(c, 13.0), n).side() + 0.5);
    }
}

TEST(Squeezer, SymmetricSqueezeStaysSymmetric) {
    std::mt19937_64 rng(12);
    const Layout start = symmetricStart(16, 4.0, rng);
    SqueezeConfig cfg;
    cfg.iterations = 1; // relax only: the single-tree moves would break the symmetry
    cfg.relaxIters = 5000;
    cfg.glsRounds = 5;
    cfg.symmetric = true;
    const std::optional<Layout> l = squeeze(start, 3.2, cfg, rng);
    ASSERT_TRUE(l.has_value());
    EXPECT_TRUE(l->isValid());
    EXPECT_LE(l->side(), 3.2);
    const Box b = l->bounds();
    const double cx = (b.minX + b.maxX) / 2, cy = (b.minY + b.maxY) / 2;
    for (std::size_t i = 0; i < 8; ++i) {
        EXPECT_NEAR(l->poses[i].x + l->poses[i + 8].x, 2 * cx, 1e-9);
        EXPECT_NEAR(l->poses[i].y + l->poses[i + 8].y, 2 * cy, 1e-9);
    }
}

TEST(Squeezer, OnlyMovableTreesMove) {
    std::mt19937_64 rng(13);
    const Layout cut = cutLattice(testLattice(), 60);
    SqueezeConfig cfg;
    cfg.iterations = 100000;
    cfg.relaxIters = 2000;
    cfg.movable = borderTrees(cut, 0.8);
    const int fixed = static_cast<int>(std::count(cfg.movable.begin(), cfg.movable.end(), false));
    ASSERT_GT(fixed, 0);
    ASSERT_LT(fixed, 60);
    const std::optional<Layout> l = squeeze(cut, cut.side() * 1.02, cfg, rng);
    ASSERT_TRUE(l.has_value());
    EXPECT_TRUE(l->isValid());
    for (std::size_t i = 0; i < cut.size(); ++i) {
        if (cfg.movable[i]) continue;
        EXPECT_EQ(l->poses[i].x, cut.poses[i].x);
        EXPECT_EQ(l->poses[i].y, cut.poses[i].y);
        EXPECT_EQ(l->poses[i].deg, cut.poses[i].deg);
    }
}

TEST(Squeezer, GuidedSearchUntanglesADenseStart) {
    std::mt19937_64 rng(14);
    const Layout start = denseStart(12, 2.0, rng);
    SqueezeConfig cfg;
    cfg.iterations = 50000;
    cfg.relaxIters = 3000;
    cfg.glsRounds = 10;
    cfg.lineShare = 0.3;
    const std::optional<Layout> l = squeeze(start, 2.6, cfg, rng);
    ASSERT_TRUE(l.has_value());
    EXPECT_TRUE(l->isValid());
    EXPECT_LE(l->side(), 2.6);
}

TEST(Squeezer, CompressClosesTheWallsOnALooseLayout) {
    std::mt19937_64 rng(15);
    Layout l = randomStart(20, rng);
    const double before = l.side();
    const double side = compress(l, CompressConfig{});
    EXPECT_LT(side, 0.9 * before);
    SqueezeConfig cfg;
    cfg.iterations = 20000;
    cfg.relaxIters = 2000;
    const std::optional<Layout> fitted = squeeze(l, side * 1.001, cfg, rng);
    ASSERT_TRUE(fitted.has_value());
    EXPECT_TRUE(fitted->isValid());
}

TEST(Nfp, CornersTouchWithoutOverlap) {
    const packer::tree::Placed fixed = packer::tree::place({0, 0, 30});
    int clear = 0, total = 0;
    for (const auto &p : nfpCorners(fixed, 200, 1e-7)) {
        ++total;
        if (!packer::tree::collide(fixed, packer::tree::place({p.x, p.y, 200}))) ++clear;
    }
    // A corner of one pair's polygon may lie inside another pair's, but some are clear.
    EXPECT_GT(total, 50);
    EXPECT_GT(clear, 10);
    // Pulled into its polygon, every corner overlaps: the polygons are where parts meet.
    for (const auto &p : nfpCorners(fixed, 200, -1e-4)) {
        EXPECT_TRUE(packer::tree::collide(fixed, packer::tree::place({p.x, p.y, 200})));
    }
}

TEST(Nfp, BestInsertionFindsAFreeSpot) {
    std::mt19937_64 rng(16);
    Layout l = randomStart(10, rng);
    const Box b = l.bounds();
    l.add(bestInsertion(l, b, 4, rng));
    EXPECT_TRUE(l.isValid());
    const Box t = l.boxes.back();
    EXPECT_GE(t.minX, b.minX);
    EXPECT_LE(t.maxX, b.maxX);
    EXPECT_GE(t.minY, b.minY);
    EXPECT_LE(t.maxY, b.maxY);
}

TEST(Nfp, RuinAndRecreateKeepsTheTreesAndSettles) {
    std::mt19937_64 rng(18);
    const Layout l = randomStart(30, rng);
    const Box b = l.bounds();
    const double half = l.side() / 2;
    const double cx = (b.minX + b.maxX) / 2, cy = (b.minY + b.maxY) / 2;
    const Box square{cx - half, cy - half, cx + half, cy + half};
    for (int k = 0; k < 3; ++k) {
        const Layout r = ruinAndRecreate(l, square, 6, rng);
        ASSERT_EQ(r.size(), l.size());
        SqueezeConfig sc;
        sc.iterations = 200000;
        sc.relaxIters = 2000;
        const std::optional<Layout> settled = squeeze(r, l.side() * 1.05, sc, rng);
        ASSERT_TRUE(settled.has_value());
        EXPECT_TRUE(settled->isValid());
    }
}

TEST(Crossover, ChildHasNTrees) {
    std::mt19937_64 rng(17);
    const Layout a = randomStart(15, rng);
    const Layout b = randomStart(15, rng);
    const Layout c = randomStart(18, rng);
    for (int k = 0; k < 5; ++k) {
        EXPECT_EQ(crossover(a, b, 15, rng).size(), 15u);
        EXPECT_EQ(crossover(a, c, 15, rng).size(), 15u);
        EXPECT_EQ(crossover(c, a, 18, rng).size(), 18u);
    }
}

TEST(Annealer, StaysValidAndNeverWorse) {
    Layout l = cutLattice(testLattice(), 10);
    const double before = l.side();
    AnnealConfig cfg;
    cfg.iterations = 20000;
    std::mt19937_64 rng(1);
    anneal(l, cfg, rng);
    EXPECT_TRUE(l.isValid());
    EXPECT_LE(l.side(), before);
}

TEST(Rotation, RotateWholeKeepsValidAndNeverWorse) {
    Layout l = cutLattice(testLattice(), 12);
    const double before = l.side();
    rotateWhole(l);
    EXPECT_TRUE(l.isValid());
    EXPECT_LE(l.side(), before);
}

TEST(Solver, SmallRunIsValid) {
    SolverConfig cfg;
    cfg.maxN = 15;
    cfg.itersPerTree = 200;
    cfg.restarts = 2;
    cfg.randomItersPerTree = 200;
    cfg.mediumRandom = 1;
    cfg.latticeStarts = 2;
    cfg.latticeKeep = 2;
    cfg.latticeIters = 200;
    const Solution sol = solveAll(cfg);
    for (int n = 1; n <= cfg.maxN; ++n) {
        EXPECT_EQ(static_cast<int>(sol.layouts[n].size()), n);
        EXPECT_TRUE(sol.layouts[n].isValid()) << "n=" << n;
    }
    EXPECT_GT(totalScore(sol.layouts), 0.0);
}

TEST(Improver, KeepsLayoutsValidAndNeverWorse) {
    std::mt19937_64 rng(5);
    std::vector<Layout> best(9);
    best[1] = bestSingle();
    for (int n = 2; n <= 8; ++n)
        best[n] = randomStart(n, rng);
    std::vector<double> before(best.size());
    for (std::size_t n = 1; n < best.size(); ++n)
        before[n] = best[n].side();
    std::vector<std::string> origins(best.size(), "start");

    ImproveConfig cfg;
    cfg.minutes = 0.03;
    cfg.threads = 2;
    cfg.maxN = 8;
    cfg.seed = 1;
    cfg.randomItersPerTree = 2000;
    cfg.refineItersPerTree = 2000;
    cfg.squeezeItersPerTree = 2000;
    cfg.crystals = {toCrystal(testLattice())};
    cfg.crystalMinN = 5;
    std::atomic<bool> stop{false};
    int saves = 0;
    improve(best, origins, cfg, stop,
            [&](const std::vector<Layout> &, const std::vector<std::string> &) { ++saves; });

    EXPECT_GE(saves, 1);
    EXPECT_LT(
        totalScore(best), 0.9 * [&] {
            double total = 0;
            for (std::size_t n = 1; n < best.size(); ++n)
                total += before[n] * before[n] / n;
            return total;
        }()); // random starts are loose, so a short run must gain a lot
    for (std::size_t n = 1; n < best.size(); ++n) {
        EXPECT_EQ(best[n].size(), n);
        EXPECT_TRUE(best[n].isValid()) << "n=" << n;
        EXPECT_LE(best[n].side(), before[n]) << "n=" << n;
    }
}
