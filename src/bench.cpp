// Small, quick benchmarks for trying out a change without a full run.
//
//   bench patterns                      lattice search, then the best cut for the large batch
//   bench improve <stage> [batch|n...]  apply one stage (so far: anneal) to results.csv
//   bench band [maxN]                   a full solve up to maxN (default 31)
//
// Batches: "medium" (n 13-60) and "large" (n 80-200). Run from the repository root.
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "io/csv_reader.h"
#include "io/submission.h"
#include "packing/annealer.h"
#include "packing/rotation.h"
#include "packing/lattice.h"
#include "packing/solver.h"

namespace {
using namespace packer;
using packing::Layout;

const char *const kResults = "data/current-best/results.csv";
const std::vector<int> kMedium = {13, 17, 20, 24, 30, 40, 50, 60};
const std::vector<int> kLarge = {80, 95, 119, 158, 200};

double seconds(std::chrono::steady_clock::time_point since) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - since).count();
}

double scorePart(const Layout &layout) {
    const double s = layout.side();
    return s * s / static_cast<double>(layout.size());
}

// The n values named by args: batch names or numbers. Defaults to the large batch.
std::vector<int> pickN(const std::vector<std::string> &args) {
    std::vector<int> ns;
    for (const std::string &a : args) {
        if (a == "medium") {
            ns.insert(ns.end(), kMedium.begin(), kMedium.end());
        } else if (a == "large") {
            ns.insert(ns.end(), kLarge.begin(), kLarge.end());
        } else {
            ns.push_back(std::stoi(a));
        }
    }
    return ns.empty() ? kLarge : ns;
}

int patterns() {
    const auto started = std::chrono::steady_clock::now();
    const packing::SolverConfig cfg;
    const auto lattices = packing::findLattices(cfg);
    std::printf("[%.0fs] %zu lattices, densest %.6f per tree\n", seconds(started),
                lattices.size(), lattices.empty() ? 0.0 : packing::areaPerTree(lattices[0]));
    double total = 0;
    for (int n : kLarge) {
        double best = INFINITY;
        for (const auto &p : lattices) {
            best = std::min(best, scorePart(packing::cutLattice(p, n)));
        }
        std::printf("n=%-4d best raw cut %.6f\n", n, best);
        total += best;
    }
    std::printf("mean raw cut %.6f per n\n", total / kLarge.size());
    return 0;
}

int improve(const std::string &stage, const std::vector<int> &ns) {
    if (stage != "anneal") {
        std::cerr << "unknown stage: " << stage << " (anneal)\n";
        return 1;
    }
    const auto layouts = io::loadSubmission(io::readCsv(kResults));
    const packing::SolverConfig cfg;
    std::vector<std::string> lines(ns.size());
    std::vector<double> gains(ns.size(), 0.0);
    const auto started = std::chrono::steady_clock::now();

#pragma omp parallel for schedule(dynamic)
    for (std::size_t k = 0; k < ns.size(); ++k) {
        const int n = ns[k];
        if (n < 1 || n >= static_cast<int>(layouts.size())) continue;
        Layout layout = layouts[n];
        std::mt19937_64 rng(cfg.seed + n);
        const double before = scorePart(layout);
        const auto t0 = std::chrono::steady_clock::now();

        packing::AnnealConfig ac = cfg.anneal; // a long, cold anneal
        ac.iterations = cfg.largeItersPerTree * n;
        ac.startTemp = cfg.latticeTemp;
        ac.startStep = cfg.latticeStep;
        packing::anneal(layout, ac, rng);
        packing::rotateWhole(layout);

        const double after = scorePart(layout);
        gains[k] = before - after;
        char line[160];
        std::snprintf(line, sizeof line, "n=%-4d %.6f -> %.6f  gain %.6f  %5.0fs%s", n, before,
                      after, gains[k], seconds(t0), layout.isValid() ? "" : "  INVALID");
        lines[k] = line;
    }
    double total = 0;
    for (std::size_t k = 0; k < ns.size(); ++k) {
        std::cout << lines[k] << "\n";
        total += gains[k];
    }
    std::printf("%s: mean gain %.6f per n, %.0fs in total\n", stage.c_str(), total / ns.size(),
                seconds(started));
    return 0;
}

int band(int maxN) {
    packing::SolverConfig cfg;
    cfg.maxN = maxN;
    const auto started = std::chrono::steady_clock::now();
    const auto sol = packing::solveAll(cfg);
    double small = 0;
    double medium = 0;
    for (int n = 1; n <= std::min(30, maxN); ++n) {
        (n <= 12 ? small : medium) += scorePart(sol.layouts[n]);
    }
    std::printf("1-12 = %.5f  13-%d = %.5f  (%.0fs)\n", small, std::min(30, maxN), medium,
                seconds(started));
    return 0;
}
} // namespace

int main(int argc, char *argv[]) {
    const std::vector<std::string> args(argv + 1, argv + argc);
    try {
        if (!args.empty() && args[0] == "patterns") return patterns();
        if (args.size() >= 2 && args[0] == "improve") {
            return improve(args[1], pickN({args.begin() + 2, args.end()}));
        }
        if (!args.empty() && args[0] == "band") {
            return band(args.size() >= 2 ? std::stoi(args[1]) : 31);
        }
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    std::cerr << "Usage: bench patterns\n"
              << "       bench improve anneal [medium|large|n...]\n"
              << "       bench band [maxN]\n";
    return 1;
}
