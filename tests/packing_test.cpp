#include "io/csv_reader.h"
#include "io/tree_loader.h"
#include "packing/annealer.h"
#include "packing/lattice.h"
#include "packing/solver.h"
#include "packing/starts.h"

#include <gtest/gtest.h>

using namespace packer::packing;

namespace {
Polygon tree() { return packer::io::loadPolygon(packer::io::readCsv(TREE_CSV_PATH)); }
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

TEST(Lattice, ValidForManySizes) {
    const Polygon shape = tree();
    const LatticeParams params = findLattice(shape, 0);
    for (int n : {1, 2, 5, 17, 50}) {
        const Layout l = buildLattice(shape, params, n);
        EXPECT_EQ(static_cast<int>(l.size()), n);
        EXPECT_TRUE(l.isValid()) << "n=" << n;
    }
}

TEST(Annealer, StaysValidAndNeverWorse) {
    const Polygon shape = tree();
    Layout l = buildLattice(shape, findLattice(shape, 0), 10);
    const double before = l.side();
    AnnealConfig cfg;
    cfg.iterations = 20000;
    std::mt19937_64 rng(1);
    anneal(l, cfg, rng);
    EXPECT_TRUE(l.isValid());
    EXPECT_LE(l.side(), before);
}

TEST(Solver, SmallRunIsValid) {
    SolverConfig cfg;
    cfg.maxN = 15;
    cfg.itersPerTree = 200;
    cfg.restarts = 2;
    cfg.latticeAngles = {0, 45};
    const Solution sol = solveAll(tree(), cfg);
    for (int n = 1; n <= cfg.maxN; ++n) {
        EXPECT_EQ(static_cast<int>(sol.layouts[n].size()), n);
        EXPECT_TRUE(sol.layouts[n].isValid()) << "n=" << n;
    }
    EXPECT_GT(score(sol), 0.0);
}
