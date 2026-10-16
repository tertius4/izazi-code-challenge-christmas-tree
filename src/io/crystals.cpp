#include "io/crystals.h"

#include <cstdio>
#include <stdexcept>

namespace packer::io {
namespace {
std::string exact(double v) {
    char buf[32];
    std::snprintf(buf, sizeof buf, "%.17g", v);
    return buf;
}
} // namespace

std::vector<std::vector<std::string>> crystalTable(const std::vector<packing::Crystal> &crystals) {
    std::vector<std::vector<std::string>> rows{{"crystal", "ux", "uy", "vx", "vy", "x", "y", "deg"}};
    for (std::size_t c = 0; c < crystals.size(); ++c) {
        const packing::Crystal &k = crystals[c];
        for (const geometry::Pose &t : k.cell) {
            rows.push_back({std::to_string(c), exact(k.u.x), exact(k.u.y), exact(k.v.x),
                            exact(k.v.y), exact(t.x), exact(t.y), exact(t.deg)});
        }
    }
    return rows;
}

std::vector<packing::Crystal> loadCrystals(const std::vector<std::vector<std::string>> &table) {
    std::vector<packing::Crystal> crystals;
    for (std::size_t r = 1; r < table.size(); ++r) { // row 0 is the header
        const auto &row = table[r];
        if (row.size() < 8) throw std::runtime_error("bad crystal row " + std::to_string(r));

        const std::size_t id = std::stoul(row[0]);
        if (id == crystals.size()) {
            crystals.push_back({{std::stod(row[1]), std::stod(row[2])},
                                {std::stod(row[3]), std::stod(row[4])},
                                {}});
        } else if (id + 1 != crystals.size()) {
            throw std::runtime_error("crystal rows out of order at row " + std::to_string(r));
        }
        crystals.back().cell.push_back({std::stod(row[5]), std::stod(row[6]), std::stod(row[7])});
    }
    return crystals;
}
} // namespace packer::io
