#pragma once

#include <iosfwd>
#include <optional>
#include <string>

namespace packer::cli {
struct Args {
    std::string inputPath;
    std::string outputPath;
    bool help = false;
};

// Parses command-line arguments. Returns nullopt on invalid usage.
std::optional<Args> parseArgs(int argc, const char *const argv[]);

void printUsage(std::ostream &out);
} // namespace packer::cli
