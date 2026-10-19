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

std::string formatValue(double value) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(15) << value;
    return "s" + out.str();
}

double parseValue(const std::string &field) {
    const bool hasPrefix = !field.empty() && field.front() == 's';
    if (!hasPrefix) throw std::runtime_error("bad value: " + field);

    return std::stod(field.substr(1));
}
} // namespace

std::vector<std::vector<std::string>> submissionTable(const std::vector<packing::Layout> &layouts) {
    std::vector<std::vector<std::string>> rows{{"id", "x", "y", "deg"}};

    for (std::size_t n = 1; n < layouts.size(); ++n) {
        const packing::Layout &layout = layouts[n];
        // Moving every tree by the same amount keeps the packing valid and stays within ±100.
        const auto [cx, cy] = geometry::centre(layout.bounds());

        for (std::size_t i = 0; i < layout.size(); ++i) {
            const geometry::Pose &pose = layout.poses[i];
            rows.push_back({formatId(n, i), formatValue(pose.x - cx), formatValue(pose.y - cy),
                            formatValue(geometry::normalizeDeg(pose.deg))});
        }
    }

    return rows;
}

std::vector<packing::Layout> loadSubmission(const std::vector<std::vector<std::string>> &table) {
    std::vector<packing::Layout> layouts(1);

    for (std::size_t rowIndex = 1; rowIndex < table.size(); ++rowIndex) { // row 0 is the header
        const auto &row = table[rowIndex];
        const std::string id = row.empty() ? "" : row[0]; // like "050_3"
        const std::size_t underscore = id.find('_');
        const bool validRow = row.size() >= 4 && underscore != std::string::npos;
        if (!validRow) {
            throw std::runtime_error("bad submission row " + std::to_string(rowIndex));
        }

        const std::size_t n = std::stoul(id.substr(0, underscore));
        while (layouts.size() <= n) {
            layouts.emplace_back();
        }

        layouts[n].add({parseValue(row[1]), parseValue(row[2]), parseValue(row[3])});
    }

    return layouts;
}
} // namespace packer::io
