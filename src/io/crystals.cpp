#include "io/crystals.h"

#include <cstdio>
#include <stdexcept>

namespace packer::io {
namespace {
std::string exact(double value) {
    char text[32];
    std::snprintf(text, sizeof text, "%.17g", value);
    return text;
}
} // namespace

std::vector<std::vector<std::string>> crystalTable(const std::vector<packing::Crystal> &crystals) {
    std::vector<std::vector<std::string>> rows{
        {"crystal", "ux", "uy", "vx", "vy", "x", "y", "deg"}};
    for (std::size_t id = 0; id < crystals.size(); ++id) {
        const packing::Crystal &crystal = crystals[id];
        for (const geometry::Pose &tree : crystal.cell) {
            rows.push_back({std::to_string(id), exact(crystal.u.x), exact(crystal.u.y),
                            exact(crystal.v.x), exact(crystal.v.y), exact(tree.x), exact(tree.y),
                            exact(tree.deg)});
        }
    }
    return rows;
}

std::vector<packing::Crystal> loadCrystals(const std::vector<std::vector<std::string>> &table) {
    std::vector<packing::Crystal> crystals;
    for (std::size_t rowIndex = 1; rowIndex < table.size(); ++rowIndex) { // row 0 is the header
        const auto &row = table[rowIndex];
        if (row.size() < 8) {
            throw std::runtime_error("bad crystal row " + std::to_string(rowIndex));
        }

        const std::size_t id = std::stoul(row[0]);
        if (id == crystals.size()) {
            crystals.push_back({{std::stod(row[1]), std::stod(row[2])},
                                {std::stod(row[3]), std::stod(row[4])},
                                {}});
        } else if (id + 1 != crystals.size()) {
            throw std::runtime_error("crystal rows out of order at row " +
                                     std::to_string(rowIndex));
        }
        crystals.back().cell.push_back({std::stod(row[5]), std::stod(row[6]), std::stod(row[7])});
    }
    return crystals;
}
} // namespace packer::io
