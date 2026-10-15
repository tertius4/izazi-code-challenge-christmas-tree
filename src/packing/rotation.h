#pragma once

#include "packing/layout.h"

namespace packer::packing {
// Turns the whole layout about its centre to the angle that gives the smallest axis-aligned
// square. Trees keep their relative positions, so the layout stays valid; if rounding makes
// the turned layout collide, it is left as it was. Returns true if the square got smaller.
bool rotateWhole(Layout &layout);
} // namespace packer::packing
