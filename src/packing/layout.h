#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "geometry/polygon.h"

namespace packer::packing {
using geometry::Box;
using geometry::Polygon;
using geometry::Pose;

// N copies of one shape, each with its own pose. The transformed polygons and their
// bounding boxes are cached so collision checks don't recompute them.
struct Layout {
    Polygon shape;               // the untransformed tree
    std::vector<Pose> poses;     // poses[i] places tree i
    std::vector<Polygon> placed; // shape transformed by poses[i]
    std::vector<Box> boxes;      // bounding box of placed[i]

    explicit Layout(Polygon s = {}) : shape(std::move(s)) {
    }

    std::size_t size() const {
        return poses.size();
    }
    void add(const Pose &pose);
    void setPose(std::size_t i, const Pose &pose);
    void remove(std::size_t i);

    // True if tree i overlaps or touches any other tree.
    bool collidesWithOthers(std::size_t i) const;
    // True if no pair of trees collides. Checks every pair, so use it for verification.
    bool isValid() const;

    // Box around all trees. Requires at least one tree.
    Box bounds() const;
    // Side of the smallest axis-aligned square containing all trees.
    double side() const;

    // Index of the tree whose removal leaves the smallest square. Requires at least one tree.
    std::size_t bestTreeToRemove() const;
    // Trees whose box touches the edge of the box around all trees. Only these decide the side.
    std::vector<std::size_t> edgeTrees() const;
};
} // namespace packer::packing