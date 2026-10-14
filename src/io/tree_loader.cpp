#include "io/tree_loader.h"

#include <cstddef>
#include <stdexcept>

namespace packer::io {
geometry::Polygon loadPolygon(const std::vector<std::vector<std::string>> &table) {
    geometry::Polygon poly;
    for (std::size_t r = 1; r < table.size(); ++r) {
        const auto &row = table[r];
        if (row.size() < 3) {
            throw std::runtime_error("Row " + std::to_string(r + 1) + " needs vertex,x,y");
        }
        poly.push_back({std::stod(row[1]), std::stod(row[2])});
    }
    if (poly.size() < 3) throw std::runtime_error("A polygon needs at least 3 vertices");
    return poly;
}
} // namespace packer::io
