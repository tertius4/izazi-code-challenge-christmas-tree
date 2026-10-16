#pragma once

#include <string>
#include <vector>

#include "packing/crystal.h"

namespace packer::io {
// Crystal rows: header crystal,ux,uy,vx,vy,x,y,deg, then one row per tree of each crystal's
// cell, numbered from 0. Values are written in full precision, so a loaded crystal is exactly
// the one saved.
std::vector<std::vector<std::string>> crystalTable(const std::vector<packing::Crystal> &crystals);

// The reverse. Throws std::runtime_error on a malformed row.
std::vector<packing::Crystal> loadCrystals(const std::vector<std::vector<std::string>> &table);
} // namespace packer::io
