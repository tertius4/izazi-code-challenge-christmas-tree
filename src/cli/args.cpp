#include "cli/args.h"

#include <ostream>

namespace packer::cli {
namespace {
const char *const kDefaultOutput = "data/output/results.csv";
} // namespace

std::optional<Args> parseArgs(int argc, const char *const argv[]) {
    if (argc == 2) {
        const std::string arg = argv[1];
        if (arg == "--help" || arg == "-h") {
            Args args;
            args.help = true;
            return args;
        }
    }
    if (argc < 2 || argc > 3) {
        return std::nullopt;
    }
    Args args;
    args.inputPath = argv[1];
    args.outputPath = argc == 3 ? argv[2] : kDefaultOutput;
    return args;
}

void printUsage(std::ostream &out) {
    out << "Usage: ./build/main <input.csv> [output.csv]\n"
        << "  output.csv defaults to " << kDefaultOutput << "\n";
}
} // namespace packer::cli
