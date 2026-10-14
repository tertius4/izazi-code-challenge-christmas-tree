#include "io/report.h"

#include <stdexcept>

namespace packer::io {
std::vector<std::vector<std::string>> reportTable(const packing::Solution &solution) {
    std::vector<std::vector<std::string>> report{{"n", "side", "score_part", "origin"}};
    double total = 0;
    for (std::size_t n = 1; n < solution.layouts.size(); ++n) {
        const auto &layout = solution.layouts[n];
        if (layout.size() != n || !layout.isValid()) {
            throw std::runtime_error("Invalid layout for n=" + std::to_string(n));
        }
        const double side = layout.side();
        const double part = side * side / static_cast<double>(n);
        total += part;
        report.push_back({std::to_string(n), std::to_string(side), std::to_string(part),
                          solution.origins[n]});
    }
    report.push_back({"total", "", std::to_string(total), ""});
    return report;
}
} // namespace packer::io
