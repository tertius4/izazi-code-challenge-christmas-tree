#pragma once

#include <string>
#include <vector>

#include "packing/layout.h"

namespace packer::io {
// Kaggle submission rows: header id,x,y,deg, then one row per tree with ids like "050_3"
// and values like "s0.123456789012345". Each layout is centred on (0,0).
std::vector<std::vector<std::string>> submissionTable(const std::vector<packing::Layout> &layouts);
} // namespace packer::io
