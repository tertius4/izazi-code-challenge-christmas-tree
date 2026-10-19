#include "util/parallel.h"

#include <atomic>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

using packer::util::parallelFor;

TEST(ParallelForTest, CallsEveryIndexExactlyOnce) {
    std::vector<std::atomic<int>> calls(1000);
    parallelFor(calls.size(), 4, [&](std::size_t i) { ++calls[i]; });

    for (const auto &c : calls) {
        EXPECT_EQ(c.load(), 1);
    }
}

TEST(ParallelForTest, WorksOnOneThreadAndOnEveryCore) {
    for (const int threads : {1, 0}) {
        std::atomic<int> sum{0};
        parallelFor(100, threads, [&](std::size_t i) { sum += static_cast<int>(i); });
        EXPECT_EQ(sum.load(), 4950);
    }
}

TEST(ParallelForTest, DoesNothingForNoIterations) {
    bool called = false;
    parallelFor(0, 4, [&](std::size_t) { called = true; });
    EXPECT_FALSE(called);
}

TEST(ParallelForTest, RethrowsAnExceptionFromABody) {
    auto failAt7 = [](std::size_t i) {
        if (i == 7) throw std::runtime_error("boom");
    };
    EXPECT_THROW(parallelFor(20, 4, failAt7), std::runtime_error);
}
