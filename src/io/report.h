#pragma once

#include <string>
#include <vector>

#include "packing/solver.h"

namespace packer::io {
// Per-N report: header n,side,score_part,origin, one row per layout, and a final "total" row
// holding the sum of all score parts. Throws std::runtime_error if any layout does not have
// n trees or contains collisions.
std::vector<std::vector<std::string>> reportTable(const packing::Solution &solution);
} // namespace packer::io
