#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "geometry/polygon.h"
#include "tree/tree.h"

namespace packer::packing {
using geometry::Box;
using geometry::Polygon;
using geometry::Pose;

// N trees, each with its own pose. Each tree's placed parts and box are cached, and a grid
// over the trees' centres means a collision check only looks at nearby trees.
struct Layout {
    std::vector<Pose> poses;          // poses[i] places tree i
    std::vector<tree::Placed> placed; // tree i at poses[i]
    std::vector<Box> boxes;           // placed[i].box, kept separately for quick scans

    std::size_t size() const {
        return poses.size();
    }
    void add(const Pose &pose);
    void setPose(std::size_t i, const Pose &pose);
    void remove(std::size_t i);

    // True if tree i overlaps or touches any other tree.
    bool collidesWithOthers(std::size_t i) const;
    // True if any two trees collide. Uses the grid, so it is quick.
    bool anyCollision() const;
    // True if no pair of trees collides. Checks every pair without the grid, so it is an
    // independent verification.
    bool isValid() const;

    // Box around all trees. Requires at least one tree.
    Box bounds() const;
    // Side of the smallest axis-aligned square containing all trees.
    double side() const;

    // Index of the tree whose removal leaves the smallest square. Requires at least one tree.
    std::size_t bestTreeToRemove() const;
    // Trees whose box touches the edge of the box around all trees. Only these decide the side.
    std::vector<std::size_t> edgeTrees() const;

  private:
    // Grid cells are 2 * tree::kRadius wide, so trees that can collide have centres in the
    // same or neighbouring cells. Centres outside the grid are clamped to its border cells,
    // which keeps that true; the grid is rebuilt when a centre leaves it.
    double originX = 0;
    double originY = 0;
    int cols = 0;
    int rows = 0;
    std::vector<std::vector<std::uint32_t>> cells;
    std::vector<int> cellOf; // cellOf[i]: the cell tree i is listed in

    void rebuildGrid();
    // The cell of point (x, y); `inside` says whether it lies within the grid unclamped.
    int cellAt(double x, double y, bool &inside) const;
    void unlist(std::size_t i);
};
} // namespace packer::packing
