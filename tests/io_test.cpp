#include "io/csv_reader.h"
#include "io/submission.h"
#include "io/tree_loader.h"

#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

using namespace packer;

TEST(TreeLoader, LoadsTreeCsv) {
    const auto tree = io::loadPolygon(io::readCsv(TREE_CSV_PATH));
    EXPECT_EQ(tree.size(), 15u);
    // 3 tiers + trunk: 0.0375 + 0.065625 + 0.1125 + 0.03
    EXPECT_NEAR(geometry::polygonArea(tree), 0.245625, 1e-12);
}

TEST(TreeLoader, RejectsTooFewVertices) {
    EXPECT_THROW(io::loadPolygon({{"vertex", "x", "y"}, {"1", "0", "0"}}), std::runtime_error);
}

TEST(Submission, KaggleFormat) {
    const auto tree = io::loadPolygon(io::readCsv(TREE_CSV_PATH));
    std::vector<packing::Layout> layouts(3, packing::Layout(tree));
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
