#include "cli/args.h"

#include <gtest/gtest.h>

#include <sstream>
#include <vector>

using packer::cli::parseArgs;
using packer::cli::printUsage;

namespace {
std::optional<packer::cli::Args> parse(std::vector<const char *> argv) {
    argv.insert(argv.begin(), "main");
    return parseArgs(static_cast<int>(argv.size()), argv.data());
}
} // namespace

TEST(ArgsTest, Defaults) {
    const auto args = parse({});
    ASSERT_TRUE(args.has_value());
    EXPECT_EQ(args->minutes, 10);
    EXPECT_EQ(args->threads, 0);
    EXPECT_EQ(args->minN, 1);
    EXPECT_EQ(args->maxN, 200);
    EXPECT_EQ(args->bestDir, "data/current-best");
    EXPECT_FALSE(args->help);
}

TEST(ArgsTest, AllFlags) {
    const auto args = parse({"--minutes", "2.5", "--threads", "4", "--n", "13-30", "--best", "x"});
    ASSERT_TRUE(args.has_value());
    EXPECT_EQ(args->minutes, 2.5);
    EXPECT_EQ(args->threads, 4);
    EXPECT_EQ(args->minN, 13);
    EXPECT_EQ(args->maxN, 30);
    EXPECT_EQ(args->bestDir, "x");
}

TEST(ArgsTest, SingleN) {
    const auto args = parse({"--n", "50"});
    ASSERT_TRUE(args.has_value());
    EXPECT_EQ(args->minN, 50);
    EXPECT_EQ(args->maxN, 50);
}

TEST(ArgsTest, HelpFlags) {
    for (const char *flag : {"--help", "-h"}) {
        const auto args = parse({flag});
        ASSERT_TRUE(args.has_value());
        EXPECT_TRUE(args->help);
    }
}

TEST(ArgsTest, InvalidUsage) {
    EXPECT_FALSE(parse({"--minutes"}).has_value());        // missing value
    EXPECT_FALSE(parse({"--minutes", "ten"}).has_value()); // not a number
    EXPECT_FALSE(parse({"--minutes", "10x"}).has_value()); // trailing junk
    EXPECT_FALSE(parse({"--minutes", "0"}).has_value());
    EXPECT_FALSE(parse({"--n", "30-13"}).has_value());
    EXPECT_FALSE(parse({"--n", "13-3x"}).has_value());
    EXPECT_FALSE(parse({"--n", "0-5"}).has_value());
    EXPECT_FALSE(parse({"--n", "1-201"}).has_value());
    EXPECT_FALSE(parse({"--colour", "red"}).has_value());
    EXPECT_FALSE(parse({"input.csv"}).has_value());
}

TEST(ArgsTest, UsageMentionsTheFlags) {
    std::ostringstream out;
    printUsage(out);
    for (const char *word : {"Usage:", "--minutes", "--threads", "--n", "--best"}) {
        EXPECT_NE(out.str().find(word), std::string::npos) << word;
    }
}
