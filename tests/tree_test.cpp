#include "geometry/collision.h"
#include "io/csv_reader.h"
#include "tree/tree.h"

#include <cmath>
#include <numbers>
#include <random>
#include <string>

#include <gtest/gtest.h>

using namespace packer;
using geometry::Point;
using geometry::Pose;

namespace {
bool inPart(const Point &point, const tree::Part &part) {
    for (std::size_t k = 0; k < part.count; ++k) {
        const Point &normal = part.normal[k];
        const Point &vertex = part.vertex[k];
        if (normal.x * (point.x - vertex.x) + normal.y * (point.y - vertex.y) > 0) return false;
    }
    return true;
}

// The general check the program used before: the reference for tree::collide.
bool generalCollide(const Pose &a, const Pose &b) {
    const auto outlineA = geometry::transform(tree::outline(), a);
    const auto outlineB = geometry::transform(tree::outline(), b);
    return geometry::polygonsCollide(outlineA, geometry::boundingBox(outlineA), outlineB,
                                     geometry::boundingBox(outlineB));
}

Pose randomPose(std::mt19937_64 &rng, double range) {
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    return {range * (2 * unit(rng) - 1), range * (2 * unit(rng) - 1), 360 * unit(rng)};
}
} // namespace

TEST(Tree, OutlineMatchesDataFile) {
    // data/input/tree.csv: vertex,x,y,description, one row per vertex after the header.
    const auto rows = io::readCsv(TREE_CSV_PATH);
    ASSERT_EQ(rows.size(), tree::outline().size() + 1);
    for (std::size_t i = 0; i < tree::outline().size(); ++i) {
        EXPECT_DOUBLE_EQ(std::stod(rows[i + 1][1]), tree::outline()[i].x);
        EXPECT_DOUBLE_EQ(std::stod(rows[i + 1][2]), tree::outline()[i].y);
    }
    // 3 tiers + trunk: 0.0375 + 0.065625 + 0.1125 + 0.03
    EXPECT_NEAR(geometry::polygonArea(tree::outline()), 0.245625, 1e-12);
}

TEST(Tree, PartsCoverExactlyTheOutline) {
    const tree::Placed placed = tree::place({0, 0, 0});
    std::mt19937_64 rng(1);
    std::uniform_real_distribution<double> randomX(-0.4, 0.4);
    std::uniform_real_distribution<double> randomY(-0.25, 0.85);
    int inside = 0;
    for (int i = 0; i < 200000; ++i) {
        const Point point{randomX(rng), randomY(rng)};
        bool inAny = false;
        for (const auto &part : placed.parts) {
            inAny = inAny || inPart(point, part);
        }
        const bool inOutline = geometry::pointInPolygon(point, tree::outline());
        EXPECT_EQ(inAny, inOutline) << point.x << ", " << point.y;
        inside += inOutline;
    }
    EXPECT_GT(inside, 1000);
}

TEST(Tree, EveryVertexWithinRadius) {
    for (const Point &vertex : tree::outline()) {
        EXPECT_LE(std::hypot(vertex.x, vertex.y - tree::kCentreY), tree::kRadius);
    }
}

TEST(Tree, BoxMatchesOutline) {
    const Pose pose{1.5, -2, 37};
    const tree::Placed placed = tree::place(pose);
    const auto expected = geometry::boundingBox(geometry::transform(tree::outline(), pose));
    EXPECT_NEAR(placed.box.minX, expected.minX, 1e-12);
    EXPECT_NEAR(placed.box.minY, expected.minY, 1e-12);
    EXPECT_NEAR(placed.box.maxX, expected.maxX, 1e-12);
    EXPECT_NEAR(placed.box.maxY, expected.maxY, 1e-12);
}

TEST(Tree, CollideAgreesWithGeneralCheckOnRandomPairs) {
    std::mt19937_64 rng(2);
    int hits = 0;
    for (int i = 0; i < 1000000; ++i) {
        const Pose a = randomPose(rng, 0.6);
        const Pose b = randomPose(rng, 0.6);
        const bool expected = generalCollide(a, b);
        ASSERT_EQ(tree::collide(tree::place(a), tree::place(b)), expected) << "pair " << i;
        hits += expected;
    }
    EXPECT_GT(hits, 100000); // the range makes both outcomes common
}

TEST(Tree, CollideAgreesNearContact) {
    // Slide b toward a until they just touch, then test a hair either side of contact.
    std::mt19937_64 rng(3);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    for (int i = 0; i < 20000; ++i) {
        const Pose a{0, 0, 360 * unit(rng)};
        const double direction = 2 * std::numbers::pi * unit(rng);
        const double deg = 360 * unit(rng);
        auto at = [&](double dist) {
            return Pose{dist * std::cos(direction), dist * std::sin(direction), deg};
        };
        double lo = 0;
        double hi = 2.1;
        for (int it = 0; it < 60; ++it) {
            const double mid = (lo + hi) / 2;
            (generalCollide(a, at(mid)) ? lo : hi) = mid;
        }
        // Trees can collide again closer in, but just outside the last contact they are apart.
        for (const double dist : {hi + 1e-7, hi + 1e-9}) {
            ASSERT_FALSE(tree::collide(tree::place(a), tree::place(at(dist))))
                << i << " dist=" << dist;
        }
        for (const double dist : {lo - 1e-9, lo - 1e-7}) {
            ASSERT_TRUE(tree::collide(tree::place(a), tree::place(at(dist))))
                << i << " dist=" << dist;
        }
    }
}

TEST(Tree, PenetrationPositiveExactlyOnCollision) {
    std::mt19937_64 rng(6);
    for (int i = 0; i < 200000; ++i) {
        const auto a = tree::place(randomPose(rng, 0.6));
        const auto b = tree::place(randomPose(rng, 0.6));
        const double depth = tree::penetration(a, b);
        ASSERT_GE(depth, 0) << "pair " << i;
        ASSERT_EQ(depth > 0, tree::collide(a, b)) << "pair " << i;
        ASSERT_DOUBLE_EQ(depth, tree::penetration(b, a)) << "pair " << i;
    }
}

TEST(Tree, PenetrationGrowsAsTreesPushTogether) {
    const auto a = tree::place({0, 0, 0});
    double last = 0;
    for (double x = 0.69; x > 0.4; x -= 0.05) { // bottom tiers side by side, sliding in
        const double depth = tree::penetration(a, tree::place({x, 0, 0}));
        EXPECT_GT(depth, last) << "x=" << x;
        last = depth;
    }
}

TEST(Tree, OverlapSlopeMatchesFiniteDifferences) {
    std::mt19937_64 rng(8);
    constexpr double kMargin = 1e-3;
    constexpr double kStep = 1e-7;
    int checked = 0;
    for (int i = 0; i < 20000; ++i) {
        const Pose poseA = randomPose(rng, 0.6);
        const auto b = tree::place(randomPose(rng, 0.6));
        auto depthAt = [&](double dx, double dy, double ddeg) {
            const Pose pose{poseA.x + dx, poseA.y + dy, poseA.deg + ddeg};
            return tree::overlap(tree::place(pose), {pose.x, pose.y}, b, kMargin).depth;
        };
        const tree::Overlap analytic =
            tree::overlap(tree::place(poseA), {poseA.x, poseA.y}, b, kMargin);
        if (analytic.depth <= 0) continue;

        // The depth has kinks where the shallowest axis changes: skip pairs right at one.
        const double numericDx = (depthAt(kStep, 0, 0) - depthAt(-kStep, 0, 0)) / (2 * kStep);
        const double numericDy = (depthAt(0, kStep, 0) - depthAt(0, -kStep, 0)) / (2 * kStep);
        const double numericDdeg = (depthAt(0, 0, kStep) - depthAt(0, 0, -kStep)) / (2 * kStep);
        auto straight = [&](double ahead, double behind) {
            return std::abs(ahead + behind - 2 * analytic.depth) < 1e-9;
        };
        const bool smooth = straight(depthAt(kStep, 0, 0), depthAt(-kStep, 0, 0)) &&
                            straight(depthAt(0, kStep, 0), depthAt(0, -kStep, 0)) &&
                            straight(depthAt(0, 0, kStep), depthAt(0, 0, -kStep));
        if (!smooth) continue;

        ++checked;
        ASSERT_NEAR(analytic.dx, numericDx, 1e-5) << "pair " << i;
        ASSERT_NEAR(analytic.dy, numericDy, 1e-5) << "pair " << i;
        ASSERT_NEAR(analytic.ddeg, numericDdeg, 1e-5) << "pair " << i;
    }
    EXPECT_GT(checked, 1000);
}

TEST(Tree, OverlapWithoutMarginAgreesWithPenetration) {
    std::mt19937_64 rng(9);
    for (int i = 0; i < 100000; ++i) {
        const Pose poseA = randomPose(rng, 0.6);
        const auto a = tree::place(poseA);
        const auto b = tree::place(randomPose(rng, 0.6));
        const double depth = tree::penetration(a, b);
        const double overlapDepth = tree::overlap(a, {poseA.x, poseA.y}, b, 0).depth;
        ASSERT_NEAR(overlapDepth, depth, 1e-9) << "pair " << i;
    }
}
