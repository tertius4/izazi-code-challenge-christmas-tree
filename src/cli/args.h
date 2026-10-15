#pragma once

#include <iosfwd>
#include <optional>
#include <string>

namespace packer::cli {
struct Args {
    double minutes = 10;                       // time budget for improving
    int threads = 0;                           // 0: half the cores
    int minN = 1;                              // only work on n in [minN, maxN]
    int maxN = 200;
    std::string bestDir = "data/current-best"; // best-known layouts: read, improved, saved
    bool help = false;
};

// Parses command-line arguments. Returns nullopt on invalid usage.
std::optional<Args> parseArgs(int argc, const char *const argv[]);

void printUsage(std::ostream &out);
} // namespace packer::cli
