#include <filesystem>
#include <iostream>

#include "io/csv_reader.h"
#include "io/csv_writer.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: ./build/main <input.csv>\n";
        return 1;
    }
    const std::string inputPath = argv[1];
    const std::string outputPath = "data/output/result.csv";

    try {
        auto table = packer::io::readCsv(inputPath);

        // ... packing work goes here ...

        std::filesystem::create_directories("data/output");
        if (!packer::io::writeCsv(outputPath, table)) {
            std::cerr << "Failed to write " << outputPath << "\n";
            return 1;
        }
        std::cout << "Wrote " << outputPath << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}