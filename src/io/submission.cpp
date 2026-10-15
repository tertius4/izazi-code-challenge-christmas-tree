#include "io/submission.h"

#include <cstddef>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace packer::io {
namespace {
std::string formatId(std::size_t n, std::size_t i) {
    std::ostringstream out;
    out << std::setw(3) << std::setfill('0') << n << '_' << i;
    return out.str();
}

std::string formatValue(double v) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(15) << v;
    return "s" + out.str();
}

double parseValue(const std::string &field) {
    if (field.empty() || field.front() != 's') throw std::runtime_error("bad value: " + field);
    return std::stod(field.substr(1));
}
} // namespace

std::vector<std::vector<std::string>> submissionTable(const std::vector<packing::Layout> &layouts) {
    std::vector<std::vector<std::string>> rows{{"id", "x", "y", "deg"}};
    for (std::size_t n = 1; n < layouts.size(); ++n) {
        const packing::Layout &layout = layouts[n];
        // Moving every tree by the same amount keeps the packing valid and stays within ±100.
        const geometry::Box b = layout.bounds();
        const double cx = (b.minX + b.maxX) / 2;
        const double cy = (b.minY + b.maxY) / 2;
        for (std::size_t i = 0; i < layout.size(); ++i) {
            const geometry::Pose &p = layout.poses[i];
            rows.push_back({formatId(n, i), formatValue(p.x - cx), formatValue(p.y - cy),
                            formatValue(geometry::normalizeDeg(p.deg))});
        }
    }
    return rows;
}

std::vector<packing::Layout> loadSubmission(const std::vector<std::vector<std::string>> &table) {
    std::vector<packing::Layout> layouts(1);
    for (std::size_t r = 1; r < table.size(); ++r) { // row 0 is the header
        const auto &row = table[r];
        if (row.size() < 4 || row[0].find('_') == std::string::npos) {
            throw std::runtime_error("bad submission row " + std::to_string(r));
        }
        const std::size_t n = std::stoul(row[0].substr(0, row[0].find('_')));
        while (layouts.size() <= n) layouts.emplace_back();
        layouts[n].add({parseValue(row[1]), parseValue(row[2]), parseValue(row[3])});
    }
    return layouts;
}
} // namespace packer::io
