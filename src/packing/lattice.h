#pragma once

#include <optional>
#include <random>

#include "packing/layout.h"

namespace packer::packing {
// A double lattice: a repeating cell holding two trees, A and B.
//   A(i, j) sits at i·u + j·v, with u = (a, 0) and v = (h·a, b), turned by deg.
//   B(i, j) sits at A(i, j) + s·u + t·v, turned a further 180° so the two interlock.
// Rows stay horizontal, so a square can be cut out of the pattern cleanly.
struct LatticeParams {
    double deg = 0; // rotation of tree A
    double a = 0;   // step along a row
    double b = 0;   // step from one row to the next (upwards)
    double h = 0;   // sideways shift of the next row, as a fraction of a, in [-0.5, 0.5)
    double s = 0;   // where B sits inside the cell, as fractions of u and v, in [0, 1)
    double t = 0;
};

// Area per tree (a cell holds two trees). Lower means a denser pattern.
double areaPerTree(const LatticeParams &lattice);

// Pose of tree A (second = false) or tree B (second = true) of cell (i, j).
Pose latticePose(const LatticeParams &lattice, int i, int j, bool second);

// True if no two trees of the infinite lattice collide.
bool latticeIsValid(const LatticeParams &lattice);

// Searches for the tightest valid lattice: a random cell improved by simulated annealing
// over its angle, shape and the place of B. Returns nothing if no valid cell was found.
std::optional<LatticeParams> findLattice(long iterations, std::mt19937_64 &rng);

// The n trees of the lattice that fit in the smallest square.
Layout cutLattice(const LatticeParams &lattice, int n);
} // namespace packer::packing
