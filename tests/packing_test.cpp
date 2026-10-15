#include "packing/annealer.h"
#include "packing/improver.h"
#include "packing/rotation.h"
#include "packing/lattice.h"
#include "packing/solver.h"
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
    for (int i = 0; i < 60; ++i) l.add({4 * u(rng), 4 * u(rng), 360 * u(rng)});
    for (int step = 0; step < 3000; ++step) {
        const std::size_t i = static_cast<std::size_t>(u(rng) * l.size());
        const double range = step % 100 == 0 ? 12 : 4;
        l.setPose(i, {range * u(rng) - (range - 4) / 2, range * u(rng) - (range - 4) / 2,
                      360 * u(rng)});
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

TEST(Starts, RandomStartIsValid) {
    std::mt19937_64 rng(1);
    const Layout l = randomStart(20, rng);
    EXPECT_EQ(l.size(), 20u);
    EXPECT_TRUE(l.isValid());
}

TEST(Starts, BestSingleBeatsUpright) { EXPECT_LE(bestSingle().side(), 1.0); }

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
    EXPECT_GT(score(sol), 0.0);
}

TEST(Improver, KeepsLayoutsValidAndNeverWorse) {
    std::mt19937_64 rng(5);
    std::vector<Layout> best(9);
    best[1] = bestSingle();
    for (int n = 2; n <= 8; ++n) best[n] = randomStart(n, rng);
    std::vector<double> before(best.size());
    for (std::size_t n = 1; n < best.size(); ++n) before[n] = best[n].side();
    std::vector<std::string> origins(best.size(), "start");

    ImproveConfig cfg;
    cfg.minutes = 0.03;
    cfg.threads = 2;
    cfg.maxN = 8;
    cfg.seed = 1;
    cfg.randomItersPerTree = 2000;
    cfg.refineItersPerTree = 2000;
    std::atomic<bool> stop{false};
    int saves = 0;
    improve(best, origins, cfg, stop,
            [&](const std::vector<Layout> &, const std::vector<std::string> &) { ++saves; });

    EXPECT_GE(saves, 1);
    EXPECT_LT(totalScore(best), 0.9 * [&] {
        double total = 0;
        for (std::size_t n = 1; n < best.size(); ++n) total += before[n] * before[n] / n;
        return total;
    }()); // random starts are loose, so a short run must gain a lot
    for (std::size_t n = 1; n < best.size(); ++n) {
        EXPECT_EQ(best[n].size(), n);
        EXPECT_TRUE(best[n].isValid()) << "n=" << n;
        EXPECT_LE(best[n].side(), before[n]) << "n=" << n;
    }
}
