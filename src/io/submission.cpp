#include "io/submission.h"

#include <cstddef>
#include <iomanip>
#include <sstream>

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
} // namespace packer::io
