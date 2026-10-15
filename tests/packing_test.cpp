#include "io/csv_reader.h"
#include "io/tree_loader.h"
#include "packing/annealer.h"
#include "packing/compactor.h"
#include "packing/lattice.h"
#include "packing/solver.h"
#include "packing/starts.h"

#include <gtest/gtest.h>

using namespace packer::packing;

namespace {
Polygon tree() {
    return packer::io::loadPolygon(packer::io::readCsv(TREE_CSV_PATH));
}
} // namespace

TEST(Layout, DetectsOverlapAndSide) {
    Layout l(tree());
    l.add({0, 0, 0});
    EXPECT_NEAR(l.side(), 1.0, 1e-12); // upright tree: 0.7 wide, 1.0 tall
    l.add({0.1, 0, 0});
    EXPECT_FALSE(l.isValid());
    l.setPose(1, {2, 0, 0});
    EXPECT_TRUE(l.isValid());
}

TEST(Starts, RandomStartIsValid) {
    std::mt19937_64 rng(1);
    const Layout l = randomStart(tree(), 20, rng);
    EXPECT_EQ(l.size(), 20u);
    EXPECT_TRUE(l.isValid());
}

TEST(Starts, BestSingleBeatsUpright) { EXPECT_LE(bestSingle(tree()).side(), 1.0); }

// One short lattice search shared by the tests below.
const LatticeParams &testLattice() {
    static const LatticeParams lattice = [] {
        std::mt19937_64 rng(1);
        return findLattice(tree(), 10000, rng).value();
    }();
    return lattice;
}

TEST(Lattice, PeriodicCheckAgreesWithBruteForce) {
    const Polygon shape = tree();
    const LatticeParams &p = testLattice();
    EXPECT_TRUE(latticeIsValid(shape, p));

    // Build a 10x10-cell patch and check every pair the slow way.
    Layout patch(shape);
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
    EXPECT_FALSE(latticeIsValid(shape, squashed));
}

TEST(Lattice, DenserThanOldColumnPattern) {
    // The old 0° column lattice used about 0.41 per tree once cut into a square for N=200.
    EXPECT_LT(areaPerTree(testLattice()), 0.36);
}

TEST(Lattice, CutHoldsExactlyNValidTrees) {
    const Polygon shape = tree();
    for (int n : {1, 2, 17, 50, 200}) {
        const Layout l = cutLattice(shape, testLattice(), n);
        EXPECT_EQ(static_cast<int>(l.size()), n);
        EXPECT_TRUE(l.isValid()) << "n=" << n;
    }
}

TEST(Lattice, CutBeatsOldResultForLargeN) {
    // The old column lattice gave side 9.1 for N=200.
    EXPECT_LT(cutLattice(tree(), testLattice(), 200).side(), 9.1);
}

TEST(Annealer, StaysValidAndNeverWorse) {
    const Polygon shape = tree();
    Layout l = cutLattice(shape, testLattice(), 10);
    const double before = l.side();
    AnnealConfig cfg;
    cfg.iterations = 20000;
    std::mt19937_64 rng(1);
    anneal(l, cfg, rng);
    EXPECT_TRUE(l.isValid());
    EXPECT_LE(l.side(), before);
}

TEST(Layout, CollidingWithAgreesWithCollidesWithOthers) {
    Layout l(tree());
    l.add({0, 0, 0});
    l.add({0.1, 0.1, 30});
    l.add({5, 5, 0});
    EXPECT_EQ(l.collidingWith(0), std::vector<std::size_t>{1});
    EXPECT_TRUE(l.collidingWith(2).empty());
    EXPECT_EQ(l.collidesWithOthers(1), !l.collidingWith(1).empty());
}

TEST(Compactor, RotateWholeKeepsValidAndNeverWorse) {
    Layout l = cutLattice(tree(), testLattice(), 12);
    const double before = l.side();
    rotateWhole(l);
    EXPECT_TRUE(l.isValid());
    EXPECT_LE(l.side(), before);
}

TEST(Compactor, CompactKeepsValidAndNeverWorse) {
    std::mt19937_64 rng(3);
    Layout l = randomStart(tree(), 6, rng);
    const double before = l.side();
    compact(l, rng, 100);
    EXPECT_EQ(l.size(), 6u);
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
    const Solution sol = solveAll(tree(), cfg);
    for (int n = 1; n <= cfg.maxN; ++n) {
        EXPECT_EQ(static_cast<int>(sol.layouts[n].size()), n);
        EXPECT_TRUE(sol.layouts[n].isValid()) << "n=" << n;
    }
    EXPECT_GT(score(sol), 0.0);
}
