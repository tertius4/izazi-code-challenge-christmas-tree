#pragma once

#include <atomic>
#include <optional>
#include <random>
#include <vector>

#include "packing/lattice.h"
#include "packing/layout.h"

namespace packer::packing {
// A periodic pattern of trees: a cell spanned by the vectors u and v, repeated forever. Unlike
// LatticeParams, the cell may be turned any way, the trees in it may be turned any way, and it
// may hold any number of them.
struct Crystal {
    geometry::Point u;
    geometry::Point v;
    // The trees of cell (0, 0): x and y are fractions of u and v (a tree stands at
    // x·u + y·v), deg is its turn.
    std::vector<Pose> cell;
};

// Area per tree (the cell's area over the trees it holds). Lower means denser.
double areaPerTree(const Crystal &crystal);

// Pose of tree `index` of cell (i, j).
Pose crystalPose(const Crystal &crystal, int i, int j, std::size_t index);

// True if no two trees of the infinite pattern collide.
bool crystalIsValid(const Crystal &crystal);

// The same pattern, turned by `deg` about the origin.
Crystal turned(const Crystal &crystal, double deg);

// The same pattern with a cell twice as long along u, holding twice the trees.
Crystal doubled(const Crystal &crystal);

// The two-tree double lattice as a crystal.
Crystal toCrystal(const LatticeParams &lattice);

// Searches for the densest valid crystal with `perCell` trees per cell: a random cell improved
// by simulated annealing over its shape and the trees' places and turns, each candidate scaled
// down to as small as it can be. With `paired`, trees come in pairs turned 180° from each
// other (the double lattice is perCell = 2, paired). Returns nothing if no valid cell was
// found. With `start`, the search starts from that crystal instead (and perCell is its size).
std::optional<Crystal> findCrystal(int perCell, bool paired, long iterations, std::mt19937_64 &rng,
                                   const std::atomic<bool> *stop = nullptr,
                                   const Crystal *start = nullptr);

// Runs `searches` crystal searches from random cells, then as many again from the densest
// ones doubled, on `threads` threads, and keeps the `keep` densest distinct
// ones, densest first. Prints each kept crystal's area per tree.
std::vector<Crystal> findCrystals(int searches, long iterations, int keep, int threads,
                                  unsigned seed, const std::atomic<bool> *stop = nullptr);

// The n whole trees of the pattern, as it is turned, that fit in the smallest square.
Layout cutCrystal(const Crystal &crystal, int n);

// The best of cutCrystal over the pattern turned by each of `angles` evenly spaced angles in
// [0°, 90°), starting at a random offset. Returns the cut and the angle it used.
std::pair<Layout, double> bestCrystalCut(const Crystal &crystal, int n, int angles,
                                         std::mt19937_64 &rng);
} // namespace packer::packing
