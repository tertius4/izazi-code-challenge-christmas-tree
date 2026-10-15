#include <atomic>
#include <csignal>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "cli/args.h"
#include "io/csv_reader.h"
#include "io/csv_writer.h"
#include "io/report.h"
#include "io/submission.h"
#include "packing/improver.h"
#include "packing/solver.h"

namespace {
namespace fs = std::filesystem;
using packer::packing::Layout;

std::atomic<bool> stopRequested{false};

void onInterrupt(int) {
    stopRequested = true; // the run saves and exits as soon as its threads notice
}

// Writes the file next to its final place, then renames it, so an interrupted write never
// leaves a broken store behind.
void writeAtomically(const fs::path &path, const std::vector<std::vector<std::string>> &rows) {
    const fs::path tmp = path.string() + ".tmp";
    if (!packer::io::writeCsv(tmp.string(), rows)) {
        throw std::runtime_error("Failed to write " + tmp.string());
    }

    fs::rename(tmp, path);
}

// The best-known layouts in `dir`, checked: layouts[n] must hold n trees that don't collide.
// The origins come from report.csv when it is there.
packer::packing::Solution load(const fs::path &dir) {
    packer::packing::Solution sol;
    sol.layouts = packer::io::loadSubmission(packer::io::readCsv((dir / "results.csv").string()));
    if (sol.layouts.size() != 201) {
        throw std::runtime_error("results.csv must hold layouts for n = 1..200");
    }

    for (std::size_t n = 1; n < sol.layouts.size(); ++n) {
        const Layout &layout = sol.layouts[n];
        const bool validLayout = layout.size() == n && layout.isValid();
        if (!validLayout) {
            throw std::runtime_error("results.csv: invalid layout for n=" + std::to_string(n));
        }
    }

    sol.origins.assign(sol.layouts.size(), "loaded");
    if (fs::exists(dir / "report.csv")) {
        for (const auto &row : packer::io::readCsv((dir / "report.csv").string())) {
            if (row.size() < 4) continue;

            const bool isHeaderOrTotal = row[0] == "n" || row[0] == "total";
            if (isHeaderOrTotal) continue;

            const std::size_t n = std::stoul(row[0]);
            if (n < sol.origins.size()) sol.origins[n] = row[3];
        }
    }

    return sol;
}

void save(const fs::path &dir, const std::vector<Layout> &layouts,
          const std::vector<std::string> &origins) {
    // reportTable checks every layout again before anything is written.
    const auto report = packer::io::reportTable({layouts, origins});

    fs::create_directories(dir);
    writeAtomically(dir / "results.csv", packer::io::submissionTable(layouts));
    writeAtomically(dir / "report.csv", report);
}
} // namespace

int main(int argc, char *argv[]) {
    const auto args = packer::cli::parseArgs(argc, argv);
    if (!args) {
        packer::cli::printUsage(std::cerr);
        return 1;
    }

    if (args->help) {
        packer::cli::printUsage(std::cout);
        return 0;
    }

    try {
        const fs::path dir = args->bestDir;
        packer::packing::Solution sol;
        if (fs::exists(dir / "results.csv")) {
            sol = load(dir);
        } else {
            std::cout << "No " << (dir / "results.csv").string()
                      << " yet: building every layout from scratch first (this takes a while).\n";
            sol = packer::packing::solveAll(packer::packing::SolverConfig{});
            save(dir, sol.layouts, sol.origins);
        }

        const double before = packer::packing::totalScore(sol.layouts);

        packer::packing::ImproveConfig cfg;
        cfg.minutes = args->minutes;
        cfg.threads = args->threads;
        cfg.minN = args->minN;
        cfg.maxN = args->maxN;

        std::signal(SIGINT, onInterrupt);
        packer::packing::improve(
            sol.layouts, sol.origins, cfg, stopRequested,
            [&](const std::vector<Layout> &layouts, const std::vector<std::string> &origins) {
                save(dir, layouts, origins);
            });

        const double after = packer::packing::totalScore(sol.layouts);
        std::cout << "Score: " << before << " -> " << after << "\n"
                  << "Saved in " << dir.string() << "\n";
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
