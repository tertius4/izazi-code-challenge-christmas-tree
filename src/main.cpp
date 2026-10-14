#include <filesystem>
#include <iostream>
#include <string>

#include "cli/args.h"
#include "io/csv_reader.h"
#include "io/csv_writer.h"

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
        auto table = packer::io::readCsv(args->inputPath);

        // ... packing work goes here ...

        const auto parent = std::filesystem::path(outputPath).parent_path();
        if (!parent.empty()) {
            std::filesystem::create_directories(parent);
        }
        if (!packer::io::writeCsv(outputPath, table)) {
            std::cerr << "Failed to write " << outputPath << "\n";
            return 1;
        }
        std::cout << "Wrote " << outputPath << "\n";
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
