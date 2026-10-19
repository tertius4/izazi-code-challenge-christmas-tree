#include "io/crystals.h"
#include "io/submission.h"

#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

using namespace packer;

TEST(Submission, Format) {
    std::vector<packing::Layout> layouts(3);
    layouts[1].add({0, 0, 0});
    layouts[2].add({0, 0, 0});
    layouts[2].add({1, 0, -90});
    const auto rows = io::submissionTable(layouts);
    ASSERT_EQ(rows.size(), 4u); // header + 1 + 2
    EXPECT_EQ(rows[0], (std::vector<std::string>{"id", "x", "y", "deg"}));
    EXPECT_EQ(rows[1][0], "001_0");
    EXPECT_EQ(rows[3][0], "002_1");
    EXPECT_EQ(rows[3][3], "s270.000000000000000");
    EXPECT_EQ(rows[1][1].front(), 's');
}

TEST(Submission, LoadsWhatItWrites) {
    std::vector<packing::Layout> layouts(3);
    layouts[1].add({0, 0, 10});
    layouts[2].add({-1, 0, 0});
    layouts[2].add({1, 0.5, 90});
    const auto loaded = io::loadSubmission(io::submissionTable(layouts));
    ASSERT_EQ(loaded.size(), 3u);
    ASSERT_EQ(loaded[2].size(), 2u);
    // Written layouts are centred, so compare relative positions and angles.
    EXPECT_NEAR(loaded[2].poses[1].x - loaded[2].poses[0].x, 2.0, 1e-12);
    EXPECT_NEAR(loaded[2].poses[1].y - loaded[2].poses[0].y, 0.5, 1e-12);
    EXPECT_NEAR(loaded[2].poses[1].deg, 90.0, 1e-12);
    EXPECT_NEAR(loaded[1].side(), layouts[1].side(), 1e-12);
}

TEST(Crystals, LoadsWhatItWrites) {
    packing::Crystal twoTrees{{1.25, 0}, {0.3, 0.9}, {{0, 0, 10}, {0.5, 0.25, 190}}};
    packing::Crystal oneTree{{0.7, 0.1}, {-0.2, 1.1}, {{0, 0, 1.0 / 3}}};
    const auto loaded = io::loadCrystals(io::crystalTable({twoTrees, oneTree}));
    ASSERT_EQ(loaded.size(), 2u);
    ASSERT_EQ(loaded[0].cell.size(), 2u);
    ASSERT_EQ(loaded[1].cell.size(), 1u);
    EXPECT_EQ(loaded[0].v.x, 0.3);
    EXPECT_EQ(loaded[0].cell[1].deg, 190);
    EXPECT_EQ(loaded[1].cell[0].deg, 1.0 / 3); // full precision
}
