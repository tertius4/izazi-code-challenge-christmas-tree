#include "geometry/collision.h"
#include "io/csv_reader.h"
#include "tree/tree.h"

#include <cmath>
#include <random>
#include <string>

#include <gtest/gtest.h>

using namespace packer;
using geometry::Point;
using geometry::Pose;

namespace {
bool inPart(const Point &p, const tree::Part &part) {
    for (std::size_t k = 0; k < part.count; ++k) {
        const Point &n = part.normal[k];
        if (n.x * (p.x - part.v[k].x) + n.y * (p.y - part.v[k].y) > 0) return false;
    }
    return true;
}

// The general check the program used before: the reference for tree::collide.
bool generalCollide(const Pose &a, const Pose &b) {
    const auto pa = geometry::transform(tree::outline(), a);
    const auto pb = geometry::transform(tree::outline(), b);
    return geometry::polygonsCollide(pa, geometry::boundingBox(pa), pb, geometry::boundingBox(pb));
}

Pose randomPose(std::mt19937_64 &rng, double range) {
    std::uniform_real_distribution<double> u(0.0, 1.0);
    return {range * (2 * u(rng) - 1), range * (2 * u(rng) - 1), 360 * u(rng)};
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
    const tree::Placed t = tree::place({0, 0, 0});
    std::mt19937_64 rng(1);
    std::uniform_real_distribution<double> x(-0.4, 0.4);
    std::uniform_real_distribution<double> y(-0.25, 0.85);
    int inside = 0;
    for (int i = 0; i < 200000; ++i) {
        const Point p{x(rng), y(rng)};
        bool inAny = false;
        for (const auto &part : t.parts)
            inAny = inAny || inPart(p, part);
        const bool inOutline = geometry::pointInPolygon(p, tree::outline());
        EXPECT_EQ(inAny, inOutline) << p.x << ", " << p.y;
        inside += inOutline;
    }
    EXPECT_GT(inside, 1000);
}

TEST(Tree, EveryVertexWithinRadius) {
    for (const Point &v : tree::outline()) {
        EXPECT_LE(std::hypot(v.x, v.y - tree::kCentreY), tree::kRadius);
    }
}

TEST(Tree, BoxMatchesOutline) {
    const Pose pose{1.5, -2, 37};
    const tree::Placed t = tree::place(pose);
    const auto b = geometry::boundingBox(geometry::transform(tree::outline(), pose));
    EXPECT_NEAR(t.box.minX, b.minX, 1e-12);
    EXPECT_NEAR(t.box.minY, b.minY, 1e-12);
    EXPECT_NEAR(t.box.maxX, b.maxX, 1e-12);
    EXPECT_NEAR(t.box.maxY, b.maxY, 1e-12);
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
    std::uniform_real_distribution<double> u(0.0, 1.0);
    for (int i = 0; i < 20000; ++i) {
        const Pose a{0, 0, 360 * u(rng)};
        const double dir = 2 * M_PI * u(rng);
        const double deg = 360 * u(rng);
        auto at = [&](double d) { return Pose{d * std::cos(dir), d * std::sin(dir), deg}; };
        double lo = 0;
        double hi = 2.1;
        for (int it = 0; it < 60; ++it) {
            const double mid = (lo + hi) / 2;
            (generalCollide(a, at(mid)) ? lo : hi) = mid;
        }
        // Trees can collide again closer in, but just outside the last contact they are apart.
        for (double d : {hi + 1e-7, hi + 1e-9}) {
            ASSERT_FALSE(tree::collide(tree::place(a), tree::place(at(d)))) << i << " d=" << d;
        }
        for (double d : {lo - 1e-9, lo - 1e-7}) {
            ASSERT_TRUE(tree::collide(tree::place(a), tree::place(at(d)))) << i << " d=" << d;
        }
    }
}
