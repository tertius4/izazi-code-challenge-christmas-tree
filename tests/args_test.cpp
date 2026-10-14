#include "cli/args.h"

#include <gtest/gtest.h>

#include <sstream>

using packer::cli::parseArgs;
using packer::cli::printUsage;

TEST(ArgsTest, InputOnlyUsesDefaultOutput) {
    const char *argv[] = {"main", "in.csv"};
    const auto args = parseArgs(2, argv);
    ASSERT_TRUE(args.has_value());
    EXPECT_EQ(args->inputPath, "in.csv");
    EXPECT_EQ(args->outputPath, "data/output/results.csv");
    EXPECT_FALSE(args->help);
}

TEST(ArgsTest, InputAndOutput) {
    const char *argv[] = {"main", "in.csv", "out.csv"};
    const auto args = parseArgs(3, argv);
    ASSERT_TRUE(args.has_value());
    EXPECT_EQ(args->inputPath, "in.csv");
    EXPECT_EQ(args->outputPath, "out.csv");
}

TEST(ArgsTest, HelpFlags) {
    for (const char *flag : {"--help", "-h"}) {
        const char *argv[] = {"main", flag};
        const auto args = parseArgs(2, argv);
        ASSERT_TRUE(args.has_value());
        EXPECT_TRUE(args->help);
    }
}

TEST(ArgsTest, TooFewArgsIsInvalid) {
    const char *argv[] = {"main"};
    EXPECT_FALSE(parseArgs(1, argv).has_value());
}

TEST(ArgsTest, TooManyArgsIsInvalid) {
    const char *argv[] = {"main", "a", "b", "c"};
    EXPECT_FALSE(parseArgs(4, argv).has_value());
}

TEST(ArgsTest, UsageMentionsDefaultOutput) {
    std::ostringstream out;
    printUsage(out);
    EXPECT_NE(out.str().find("Usage:"), std::string::npos);
    EXPECT_NE(out.str().find("data/output/results.csv"), std::string::npos);
}
