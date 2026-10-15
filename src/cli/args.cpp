#include "cli/args.h"

#include <ostream>
#include <stdexcept>

namespace packer::cli {
std::optional<Args> parseArgs(int argc, const char *const argv[]) {
    Args args;
    try {
        for (int i = 1; i < argc; ++i) {
            const std::string flag = argv[i];
            if (flag == "--help" || flag == "-h") {
                args.help = true;
                continue;
            }
            if (i + 1 >= argc) return std::nullopt; // every other flag takes a value
            const std::string value = argv[++i];
            std::size_t used = 0;
            if (flag == "--minutes") {
                args.minutes = std::stod(value, &used);
                if (args.minutes <= 0) return std::nullopt;
            } else if (flag == "--threads") {
                args.threads = std::stoi(value, &used);
                if (args.threads < 0) return std::nullopt;
            } else if (flag == "--n") { // "50" or "13-30"
                const std::size_t dash = value.find('-');
                const std::string first = value.substr(0, dash);
                args.minN = std::stoi(first, &used);
                if (used != first.size()) return std::nullopt;
                args.maxN = args.minN;
                if (dash != std::string::npos) {
                    const std::string second = value.substr(dash + 1);
                    args.maxN = std::stoi(second, &used);
                    if (used != second.size()) return std::nullopt;
                }
                used = value.size();
                if (args.minN < 1 || args.maxN > 200 || args.minN > args.maxN) return std::nullopt;
            } else if (flag == "--best") {
                args.bestDir = value;
                used = value.size();
            } else {
                return std::nullopt;
            }
            if (used != value.size()) return std::nullopt; // trailing junk, like "10x"
        }
    } catch (const std::logic_error &) { // stod/stoi on a non-number
        return std::nullopt;
    }
    return args;
}

void printUsage(std::ostream &out) {
    out << "Usage: ./build/main [--minutes M] [--threads T] [--n A-B] [--best DIR]\n"
        << "  Improves the best-known layouts in DIR (default data/current-best) for M minutes\n"
        << "  (default 10) on T threads (default: half the cores), only for n in A..B\n"
        << "  (default 1-200). Saves only improvements; Ctrl-C stops early and saves.\n"
        << "  If DIR has no results.csv, the layouts are first built from scratch.\n";
}
} // namespace packer::cli
