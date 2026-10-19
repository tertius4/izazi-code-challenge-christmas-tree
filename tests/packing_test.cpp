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
    Layout layout;
    layout.add({0, 0, 0});
    EXPECT_NEAR(layout.side(), 1.0, 1e-12); // upright tree: 0.7 wide, 1.0 tall
    layout.add({0.1, 0, 0});
    EXPECT_FALSE(layout.isValid());
    layout.setPose(1, {2, 0, 0});
    EXPECT_TRUE(layout.isValid());
}

TEST(Layout, GridAgreesWithBruteForce) {
    // Random moves, some far outside the grid, then compare each tree's collision check
    // through the grid with a plain check against every other tree.
    std::mt19937_64 rng(4);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    Layout layout;
    for (int i = 0; i < 60; ++i) {
        layout.add({4 * unit(rng), 4 * unit(rng), 360 * unit(rng)});
    }
    for (int step = 0; step < 3000; ++step) {
        const std::size_t i = static_cast<std::size_t>(unit(rng) * layout.size());
        const double range = step % 100 == 0 ? 12 : 4;
        layout.setPose(i, {range * unit(rng) - (range - 4) / 2, range * unit(rng) - (range - 4) / 2,
                           360 * unit(rng)});
        if (step % 500 == 0) layout.remove(static_cast<std::size_t>(unit(rng) * layout.size()));
        for (std::size_t a = 0; a < layout.size(); ++a) {
            bool expected = false;
            for (std::size_t b = 0; b < layout.size(); ++b) {
                expected = expected ||
                           (a != b && packer::tree::collide(layout.placed[a], layout.placed[b]));
            }
            ASSERT_EQ(layout.collidesWithOthers(a), expected) << "step " << step << " tree " << a;
        }
    }
}

TEST(Layout, RemovalOrderStartsWithBestTreeToRemove) {
    std::mt19937_64 rng(7);
    const Layout layout = randomStart(15, rng);
    const std::vector<std::size_t> order = layout.removalOrder();
    ASSERT_EQ(order.size(), layout.size());
    EXPECT_EQ(order.front(), layout.bestTreeToRemove());
    double last = 0;
    for (const std::size_t i : order) {
        Layout without = layout;
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
    const std::optional<Layout> squeezed = squeeze(loose, target, cfg, rng);
    ASSERT_TRUE(squeezed.has_value());
    EXPECT_EQ(squeezed->size(), 10u);
    EXPECT_TRUE(squeezed->isValid());
    EXPECT_LE(squeezed->side(), target);
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
    const std::optional<Layout> squeezed = squeeze(loose, target, cfg, rng);
    ASSERT_TRUE(squeezed.has_value());
    EXPECT_EQ(squeezed->size(), 30u);
    EXPECT_TRUE(squeezed->isValid());
    EXPECT_LE(squeezed->side(), target);
}

TEST(Squeezer, ShakenLayoutSettlesBackIntoItsSquare) {
    std::mt19937_64 rng(11);
    const Layout layout = randomStart(30, rng);
    SqueezeConfig cfg;
    cfg.iterations = 30000;
    cfg.relaxIters = 2000;
    int settled = 0;
    for (int k = 0; k < 10; ++k) {
        const Layout kicked = shake(layout, 3, rng);
        ASSERT_EQ(kicked.size(), layout.size());
        const std::optional<Layout> back = squeeze(kicked, layout.side(), cfg, rng);
        if (!back) continue;
        ++settled;
        EXPECT_TRUE(back->isValid());
        EXPECT_LE(back->side(), layout.side());
    }
    EXPECT_GT(settled, 5); // a loose layout has room for the kicked trees
}

TEST(Starts, RandomStartIsValid) {
    std::mt19937_64 rng(1);
    const Layout layout = randomStart(20, rng);
    EXPECT_EQ(layout.size(), 20u);
    EXPECT_TRUE(layout.isValid());
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
    const LatticeParams &lattice = testLattice();
    EXPECT_TRUE(latticeIsValid(lattice));

    // Build a 10x10-cell patch and check every pair the slow way.
    Layout patch;
    for (int i = 0; i < 10; ++i) {
        for (int j = 0; j < 10; ++j) {
            patch.add(latticePose(lattice, i, j, false));
            patch.add(latticePose(lattice, i, j, true));
        }
    }
    EXPECT_TRUE(patch.isValid());

    // Squashing the cell must make it collide, in both checks.
    LatticeParams squashed = lattice;
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
        const Layout cut = cutLattice(testLattice(), n);
        EXPECT_EQ(static_cast<int>(cut.size()), n);
        EXPECT_TRUE(cut.isValid()) << "n=" << n;
    }
}

TEST(Lattice, CutBeatsOldResultForLargeN) {
    // The old column lattice gave side 9.1 for N=200.
    EXPECT_LT(cutLattice(testLattice(), 200).side(), 9.1);
}

TEST(Crystal, DoubleLatticeIsTheSamePattern) {
    const Crystal crystal = toCrystal(testLattice());
    EXPECT_TRUE(crystalIsValid(crystal));
    EXPECT_NEAR(areaPerTree(crystal), areaPerTree(testLattice()), 1e-12);
    for (int i = -2; i <= 2; ++i) {
        for (int j = -2; j <= 2; ++j) {
            for (std::size_t k = 0; k < 2; ++k) {
                const Pose fromCrystal = crystalPose(crystal, i, j, k);
                const Pose fromLattice = latticePose(testLattice(), i, j, k == 1);
                EXPECT_NEAR(fromCrystal.x, fromLattice.x, 1e-9);
                EXPECT_NEAR(fromCrystal.y, fromLattice.y, 1e-9);
            }
        }
    }
}

TEST(Crystal, PeriodicCheckAgreesWithBruteForce) {
    std::mt19937_64 rng(3);
    const Crystal crystal = findCrystal(3, false, 3000, rng).value();
    ASSERT_TRUE(crystalIsValid(crystal));
    for (const Crystal &variant : {crystal, turned(crystal, 37), doubled(crystal)}) {
        EXPECT_TRUE(crystalIsValid(variant));
        Layout patch;
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                for (std::size_t k = 0; k < variant.cell.size(); ++k) {
                    patch.add(crystalPose(variant, i, j, k));
                }
            }
        }
        EXPECT_TRUE(patch.isValid());
    }
    EXPECT_NEAR(areaPerTree(doubled(crystal)), areaPerTree(crystal), 1e-12);
    const Crystal squashed{{crystal.u.x * 0.9, crystal.u.y * 0.9},
                           {crystal.v.x * 0.9, crystal.v.y * 0.9},
                           crystal.cell};
    EXPECT_FALSE(crystalIsValid(squashed));
}

TEST(Crystal, CutsHoldExactlyNValidTreesAtAnyAngle) {
    const Crystal crystal = toCrystal(testLattice());
    std::mt19937_64 rng(2);
    for (int n : {1, 7, 40, 150}) {
        for (double deg : {0.0, 13.0, 61.5}) {
            const Layout cut = cutCrystal(turned(crystal, deg), n);
            EXPECT_EQ(static_cast<int>(cut.size()), n);
            EXPECT_TRUE(cut.isValid()) << "n=" << n << " deg=" << deg;
        }
        const auto [best, deg] = bestCrystalCut(crystal, n, 8, rng);
        EXPECT_TRUE(best.isValid());
        EXPECT_LE(best.side(), cutCrystal(turned(crystal, 13.0), n).side() + 0.5);
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
    const std::optional<Layout> squeezed = squeeze(start, 3.2, cfg, rng);
    ASSERT_TRUE(squeezed.has_value());
    EXPECT_TRUE(squeezed->isValid());
    EXPECT_LE(squeezed->side(), 3.2);
    const Box bounds = squeezed->bounds();
    const double cx = (bounds.minX + bounds.maxX) / 2, cy = (bounds.minY + bounds.maxY) / 2;
    for (std::size_t i = 0; i < 8; ++i) {
        EXPECT_NEAR(squeezed->poses[i].x + squeezed->poses[i + 8].x, 2 * cx, 1e-9);
        EXPECT_NEAR(squeezed->poses[i].y + squeezed->poses[i + 8].y, 2 * cy, 1e-9);
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
    const std::optional<Layout> squeezed = squeeze(cut, cut.side() * 1.02, cfg, rng);
    ASSERT_TRUE(squeezed.has_value());
    EXPECT_TRUE(squeezed->isValid());
    for (std::size_t i = 0; i < cut.size(); ++i) {
        if (cfg.movable[i]) continue;
        EXPECT_EQ(squeezed->poses[i].x, cut.poses[i].x);
        EXPECT_EQ(squeezed->poses[i].y, cut.poses[i].y);
        EXPECT_EQ(squeezed->poses[i].deg, cut.poses[i].deg);
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
    const std::optional<Layout> squeezed = squeeze(start, 2.6, cfg, rng);
    ASSERT_TRUE(squeezed.has_value());
    EXPECT_TRUE(squeezed->isValid());
    EXPECT_LE(squeezed->side(), 2.6);
}

TEST(Squeezer, CompressClosesTheWallsOnALooseLayout) {
    std::mt19937_64 rng(15);
    Layout layout = randomStart(20, rng);
    const double before = layout.side();
    const double side = compress(layout, CompressConfig{});
    EXPECT_LT(side, 0.9 * before);
    SqueezeConfig cfg;
    cfg.iterations = 20000;
    cfg.relaxIters = 2000;
    const std::optional<Layout> fitted = squeeze(layout, side * 1.001, cfg, rng);
    ASSERT_TRUE(fitted.has_value());
    EXPECT_TRUE(fitted->isValid());
}

TEST(Nfp, CornersTouchWithoutOverlap) {
    const packer::tree::Placed fixed = packer::tree::place({0, 0, 30});
    int clear = 0, total = 0;
    for (const auto &corner : nfpCorners(fixed, 200, 1e-7)) {
        ++total;
        if (!packer::tree::collide(fixed, packer::tree::place({corner.x, corner.y, 200}))) ++clear;
    }
    // A corner of one pair's polygon may lie inside another pair's, but some are clear.
    EXPECT_GT(total, 50);
    EXPECT_GT(clear, 10);
    // Pulled into its polygon, every corner overlaps: the polygons are where parts meet.
    for (const auto &corner : nfpCorners(fixed, 200, -1e-4)) {
        EXPECT_TRUE(packer::tree::collide(fixed, packer::tree::place({corner.x, corner.y, 200})));
    }
}

TEST(Nfp, BestInsertionFindsAFreeSpot) {
    std::mt19937_64 rng(16);
    Layout layout = randomStart(10, rng);
    const Box bounds = layout.bounds();
    layout.add(bestInsertion(layout, bounds, 4, rng));
    EXPECT_TRUE(layout.isValid());
    const Box added = layout.boxes.back();
    EXPECT_GE(added.minX, bounds.minX);
    EXPECT_LE(added.maxX, bounds.maxX);
    EXPECT_GE(added.minY, bounds.minY);
    EXPECT_LE(added.maxY, bounds.maxY);
}

TEST(Nfp, RuinAndRecreateKeepsTheTreesAndSettles) {
    std::mt19937_64 rng(18);
    const Layout layout = randomStart(30, rng);
    const Box bounds = layout.bounds();
    const double half = layout.side() / 2;
    const double cx = (bounds.minX + bounds.maxX) / 2, cy = (bounds.minY + bounds.maxY) / 2;
    const Box square{cx - half, cy - half, cx + half, cy + half};
    for (int k = 0; k < 3; ++k) {
        const Layout rebuilt = ruinAndRecreate(layout, square, 6, rng);
        ASSERT_EQ(rebuilt.size(), layout.size());
        SqueezeConfig squeezeCfg;
        squeezeCfg.iterations = 200000;
        squeezeCfg.relaxIters = 2000;
        const std::optional<Layout> settled =
            squeeze(rebuilt, layout.side() * 1.05, squeezeCfg, rng);
        ASSERT_TRUE(settled.has_value());
        EXPECT_TRUE(settled->isValid());
    }
}

TEST(Crossover, ChildHasNTrees) {
    std::mt19937_64 rng(17);
    const Layout first = randomStart(15, rng);
    const Layout second = randomStart(15, rng);
    const Layout bigger = randomStart(18, rng);
    for (int k = 0; k < 5; ++k) {
        EXPECT_EQ(crossover(first, second, 15, rng).size(), 15u);
        EXPECT_EQ(crossover(first, bigger, 15, rng).size(), 15u);
        EXPECT_EQ(crossover(bigger, first, 18, rng).size(), 18u);
    }
}

TEST(Annealer, StaysValidAndNeverWorse) {
    Layout layout = cutLattice(testLattice(), 10);
    const double before = layout.side();
    AnnealConfig cfg;
    cfg.iterations = 20000;
    std::mt19937_64 rng(1);
    anneal(layout, cfg, rng);
    EXPECT_TRUE(layout.isValid());
    EXPECT_LE(layout.side(), before);
}

TEST(Rotation, RotateWholeKeepsValidAndNeverWorse) {
    Layout layout = cutLattice(testLattice(), 12);
    const double before = layout.side();
    rotateWhole(layout);
    EXPECT_TRUE(layout.isValid());
    EXPECT_LE(layout.side(), before);
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
    const Solution solution = solveAll(cfg);
    for (int n = 1; n <= cfg.maxN; ++n) {
        EXPECT_EQ(static_cast<int>(solution.layouts[n].size()), n);
        EXPECT_TRUE(solution.layouts[n].isValid()) << "n=" << n;
    }
    EXPECT_GT(totalScore(solution.layouts), 0.0);
}

TEST(Improver, KeepsLayoutsValidAndNeverWorse) {
    std::mt19937_64 rng(5);
    std::vector<Layout> best(9);
    best[1] = bestSingle();
    for (int n = 2; n <= 8; ++n) {
        best[n] = randomStart(n, rng);
    }
    std::vector<double> before(best.size());
    for (std::size_t n = 1; n < best.size(); ++n) {
        before[n] = best[n].side();
    }
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
    double scoreBefore = 0;
    for (std::size_t n = 1; n < best.size(); ++n) {
        scoreBefore += before[n] * before[n] / n;
    }
    // Random starts are loose, so a short run must gain a lot.
    EXPECT_LT(totalScore(best), 0.9 * scoreBefore);
    for (std::size_t n = 1; n < best.size(); ++n) {
        EXPECT_EQ(best[n].size(), n);
        EXPECT_TRUE(best[n].isValid()) << "n=" << n;
        EXPECT_LE(best[n].side(), before[n]) << "n=" << n;
    }
}
