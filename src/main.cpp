#include <filesystem>
#include <iostream>
#include <string>

#include "io/csv_reader.h"
#include "io/csv_writer.h"

namespace
{
    const char *const kDefaultOutput = "data/output/results.csv";

    void printUsage(std::ostream &out)
    {
        out << "Usage: ./build/main <input.csv> [output.csv]\n"
            << "  output.csv defaults to " << kDefaultOutput << "\n";
    }
} // namespace

int main(int argc, char *argv[])
{
    if (argc == 2 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h"))
    {
        printUsage(std::cout);
        return 0;
    }
    if (argc < 2 || argc > 3)
    {
        printUsage(std::cerr);
        return 1;
    }
    const std::string inputPath = argv[1];
    const std::string outputPath = argc == 3 ? argv[2] : kDefaultOutput;

    try
    {
        auto table = packer::io::readCsv(inputPath);

        // ... packing work goes here ...

        const auto parent = std::filesystem::path(outputPath).parent_path();
        if (!parent.empty())
        {
            std::filesystem::create_directories(parent);
        }
        if (!packer::io::writeCsv(outputPath, table))
        {
            std::cerr << "Failed to write " << outputPath << "\n";
            return 1;
        }
        std::cout << "Wrote " << outputPath << "\n";
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
