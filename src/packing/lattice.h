#pragma once

#include "packing/layout.h"

namespace packer::packing {
// Even columns hold trees at baseDeg; odd columns hold trees turned a further 180°,
// shifted up by `shift`, so neighbouring columns interlock.
struct LatticeParams {
    double colStep = 0;
    double rowStep = 0;
    double shift = 0;
    double baseDeg = 0;
};

// Searches for the tightest collision-free lattice at this base angle. Slow-ish; run once per angle.
LatticeParams findLattice(const Polygon &shape, double baseDeg);

// Takes n trees from the lattice, filled column by column, choosing the column count
// that gives the smallest square. The result may need checking with isValid().
Layout buildLattice(const Polygon &shape, const LatticeParams &params, int n);
} // namespace packer::packing
