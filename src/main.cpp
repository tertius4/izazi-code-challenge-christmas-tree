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
#include "io/tree_loader.h"
#include "packing/solver.h"

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

    const std::string &outputPath = args->outputPath;

    try {
        const auto inputTable = packer::io::readCsv(args->inputPath);
        const auto polygon = packer::io::loadPolygon(inputTable);

        packer::packing::SolverConfig config; // tweak effort and bands here
        const auto solution = packer::packing::solveAll(polygon, config);

        // Final safety net: throws if any layout lacks n trees or has collisions.
        const auto report = packer::io::reportTable(solution);
        const auto submissionTable = packer::io::submissionTable(solution.layouts);
        
        // Get the needed paths.
        const auto parentPath = std::filesystem::path(outputPath).parent_path();
        if (!parentPath.empty()) std::filesystem::create_directories(parentPath);
        const auto reportPath = (parentPath / "report.csv").string();

        // Write the output
        const bool isOutputWritten = packer::io::writeCsv(outputPath, submissionTable);
        if (!isOutputWritten) {
            std::cerr << "Failed to write " << outputPath << "\n";
            return 1;
        }

        // Write the report
        const bool isReportWritten = packer::io::writeCsv(reportPath, report);
        if (!isReportWritten) {
            std::cerr << "Failed to write " << reportPath << "\n";
            return 1;
        }

        // Print the score at the end.
        std::cout << "Score: " << packer::packing::score(solution) << "\n"
                  << "Wrote " << outputPath << " and " << reportPath << "\n";
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
