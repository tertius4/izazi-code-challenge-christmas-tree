#include <algorithm>
#include <atomic>
#include <csignal>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "cli/args.h"
#include "io/crystals.h"
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
    const fs::path tempPath = path.string() + ".tmp";
    if (!packer::io::writeCsv(tempPath.string(), rows)) {
        throw std::runtime_error("Failed to write " + tempPath.string());
    }

    fs::rename(tempPath, path);
}

// The best-known layouts in `dir`, checked: layouts[n] must hold n trees that don't collide.
// The origins come from report.csv when it is there.
packer::packing::Solution load(const fs::path &dir) {
    packer::packing::Solution solution;
    solution.layouts =
        packer::io::loadSubmission(packer::io::readCsv((dir / "results.csv").string()));
    if (solution.layouts.size() != 201) {
        throw std::runtime_error("results.csv must hold layouts for n = 1..200");
    }

    for (std::size_t n = 1; n < solution.layouts.size(); ++n) {
        const Layout &layout = solution.layouts[n];
        const bool validLayout = layout.size() == n && layout.isValid();
        if (!validLayout) {
            throw std::runtime_error("results.csv: invalid layout for n=" + std::to_string(n));
        }
    }

    solution.origins.assign(solution.layouts.size(), "loaded");
    if (fs::exists(dir / "report.csv")) {
        for (const auto &row : packer::io::readCsv((dir / "report.csv").string())) {
            if (row.size() < 4) continue;

            const bool isHeaderOrTotal = row[0] == "n" || row[0] == "total";
            if (isHeaderOrTotal) continue;

            const std::size_t n = std::stoul(row[0]);
            if (n < solution.origins.size()) solution.origins[n] = row[3];
        }
    }

    return solution;
}

// The densest crystals, from crystals.csv in `dir`; searched for (which takes a few minutes)
// and saved there the first time.
std::vector<packer::packing::Crystal> crystals(const fs::path &dir, int threads) {
    const fs::path path = dir / "crystals.csv";
    if (fs::exists(path)) return packer::io::loadCrystals(packer::io::readCsv(path.string()));

    std::printf("No %s yet: searching for dense crystals first.\n", path.string().c_str());
    std::fflush(stdout);
    std::vector<packer::packing::Crystal> found =
        packer::packing::findCrystals(16, 100000, 4, threads, 1, &stopRequested);
    if (!stopRequested && !found.empty()) {
        fs::create_directories(dir);
        writeAtomically(path, packer::io::crystalTable(found));
    }
    return found;
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
        packer::packing::Solution solution;
        if (fs::exists(dir / "results.csv")) {
            solution = load(dir);
        } else {
            std::printf("No %s yet: building every layout from scratch first (this takes a "
                        "while).\n",
                        (dir / "results.csv").string().c_str());
            std::fflush(stdout);
            solution = packer::packing::solveAll(packer::packing::SolverConfig{});
            save(dir, solution.layouts, solution.origins);
        }

        const double before = packer::packing::totalScore(solution.layouts);

        packer::packing::ImproveConfig cfg;
        cfg.minutes = args->minutes;
        cfg.threads = args->threads;
        cfg.minN = args->minN;
        cfg.maxN = args->maxN;
        cfg.statsPath = (dir / "stats.txt").string(); // this run's operator stats

        std::signal(SIGINT, onInterrupt);
        if (cfg.maxN >= cfg.crystalMinN) {
            const unsigned halfTheCores = std::max(1u, std::thread::hardware_concurrency() / 2);
            const int threads = cfg.threads > 0 ? cfg.threads : static_cast<int>(halfTheCores);
            cfg.crystals = crystals(dir, threads);
        }
        packer::packing::improve(
            solution.layouts, solution.origins, cfg, stopRequested,
            [&](const std::vector<Layout> &layouts, const std::vector<std::string> &origins) {
                save(dir, layouts, origins);
            });

        const double after = packer::packing::totalScore(solution.layouts);
        std::printf("Score: %g -> %g\nSaved in %s\n", before, after, dir.string().c_str());
    } catch (const std::exception &error) {
        std::fprintf(stderr, "Error: %s\n", error.what());
        return 1;
    }

    return 0;
}
