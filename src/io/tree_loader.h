#pragma once

#include <string>
#include <vector>

#include "geometry/polygon.h"

namespace packer::io {
// Builds the polygon from a table with columns vertex,x,y,... (row 0 is the header).
// Throws std::runtime_error if a row is too short or there are fewer than 3 vertices.
geometry::Polygon loadPolygon(const std::vector<std::vector<std::string>> &table);
} // namespace packer::io
