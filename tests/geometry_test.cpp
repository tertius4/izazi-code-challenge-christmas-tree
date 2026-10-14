#include "geometry/collision.h"
#include "geometry/polygon.h"

#include <gtest/gtest.h>

using namespace packer::geometry;

namespace {
Polygon square(double x, double y, double size) {
    return {{x, y}, {x + size, y}, {x + size, y + size}, {x, y + size}};
}
bool collide(const Polygon &a, const Polygon &b) {
    return polygonsCollide(a, boundingBox(a), b, boundingBox(b));
}
} // namespace

TEST(Geometry, RotatesCounterClockwiseThenTranslates) {
    const Polygon p = transform({{1, 0}}, {2, 3, 90});
    EXPECT_NEAR(p[0].x, 2, 1e-12);
    EXPECT_NEAR(p[0].y, 4, 1e-12);
}

TEST(Geometry, SquareArea) { EXPECT_DOUBLE_EQ(polygonArea(square(0, 0, 2)), 4.0); }

TEST(Geometry, NormalizesAngles) {
    EXPECT_DOUBLE_EQ(normalizeDeg(370), 10);
    EXPECT_DOUBLE_EQ(normalizeDeg(-90), 270);
}

TEST(Collision, OverlappingSquares) { EXPECT_TRUE(collide(square(0, 0, 1), square(0.5, 0.5, 1))); }
TEST(Collision, SeparatedSquares) { EXPECT_FALSE(collide(square(0, 0, 1), square(1.5, 0, 1))); }
TEST(Collision, TouchingCountsAsCollision) { EXPECT_TRUE(collide(square(0, 0, 1), square(1, 0, 1))); }

TEST(Collision, ContainedEitherWay) {
    EXPECT_TRUE(collide(square(0, 0, 4), square(1, 1, 1)));
    EXPECT_TRUE(collide(square(1, 1, 1), square(0, 0, 4)));
}

TEST(Collision, CrossingWithNoVertexInside) {
    const Polygon wide = {{0, 1}, {3, 1}, {3, 2}, {0, 2}};
    const Polygon tall = {{1, 0}, {2, 0}, {2, 3}, {1, 3}};
    EXPECT_TRUE(collide(wide, tall));
}
