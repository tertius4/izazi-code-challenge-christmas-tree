#include "packing/layout.h"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace packer::packing {
namespace {
using geometry::kInf;

constexpr double kCell = 2 * tree::kRadius;
constexpr int kMargin = 2; // spare cells around the trees when the grid is (re)built
} // namespace

void Layout::add(const Pose &pose) {
    poses.push_back(pose);
    placed.push_back(tree::place(pose));
    boxes.push_back(placed.back().box);

    rebuildGrid();
}

void Layout::setPose(std::size_t i, const Pose &pose) {
    poses[i] = pose;
    tree::place(pose, placed[i]);
    boxes[i] = placed[i].box;

    bool inside = false;
    const int cell = cellAt(placed[i].centre.x, placed[i].centre.y, inside);
    if (!inside) {
        rebuildGrid();
    } else if (cell != cellOf[i]) {
        unlist(i);
        cells[cell].push_back(static_cast<std::uint32_t>(i));
        cellOf[i] = cell;
    }
}

void Layout::remove(std::size_t i) {
    poses.erase(poses.begin() + static_cast<std::ptrdiff_t>(i));
    placed.erase(placed.begin() + static_cast<std::ptrdiff_t>(i));
    boxes.erase(boxes.begin() + static_cast<std::ptrdiff_t>(i));
    rebuildGrid(); // indices after i shifted
}

void Layout::rebuildGrid() {
    cells.clear();
    cellOf.assign(size(), 0);

    if (poses.empty()) {
        cols = rows = 0;
        return;
    }

    double minX = kInf, minY = kInf, maxX = -kInf, maxY = -kInf;
    for (const tree::Placed &tree : placed) {
        minX = std::min(minX, tree.centre.x);
        minY = std::min(minY, tree.centre.y);
        maxX = std::max(maxX, tree.centre.x);
        maxY = std::max(maxY, tree.centre.y);
    }

    originX = minX - kMargin * kCell;
    originY = minY - kMargin * kCell;
    cols = static_cast<int>((maxX - minX) / kCell) + 1 + 2 * kMargin;
    rows = static_cast<int>((maxY - minY) / kCell) + 1 + 2 * kMargin;

    cells.assign(static_cast<std::size_t>(cols * rows), {});
    for (std::size_t i = 0; i < size(); ++i) {
        bool inside = false;
        const int cell = cellAt(placed[i].centre.x, placed[i].centre.y, inside);
        cells[cell].push_back(static_cast<std::uint32_t>(i));
        cellOf[i] = cell;
    }
}

int Layout::cellAt(double x, double y, bool &inside) const {
    int cx = static_cast<int>(std::floor((x - originX) / kCell));
    int cy = static_cast<int>(std::floor((y - originY) / kCell));
    inside = cx >= 0 && cx < cols && cy >= 0 && cy < rows;

    cx = std::clamp(cx, 0, cols - 1);
    cy = std::clamp(cy, 0, rows - 1);

    return cy * cols + cx;
}

void Layout::unlist(std::size_t i) {
    auto &list = cells[cellOf[i]];
    const auto it = std::find(list.begin(), list.end(), static_cast<std::uint32_t>(i));
    *it = list.back();
    list.pop_back();
}

bool Layout::collidesWithOthers(std::size_t i) const {
    const tree::Placed &self = placed[i];
    return anyNeighbour(i, [&](std::size_t j) { return tree::collide(self, placed[j]); });
}

bool Layout::anyCollision() const {
    for (std::size_t i = 0; i < size(); ++i) {
        if (collidesWithOthers(i)) return true;
    }

    return false;
}

bool Layout::isValid() const {
    for (std::size_t i = 0; i < size(); ++i) {
        for (std::size_t j = i + 1; j < size(); ++j) {
            if (tree::collide(placed[i], placed[j])) return false;
        }
    }

    return true;
}

Box Layout::bounds() const {
    Box all = boxes[0];
    for (const Box &box : boxes) {
        all.minX = std::min(all.minX, box.minX);
        all.minY = std::min(all.minY, box.minY);
        all.maxX = std::max(all.maxX, box.maxX);
        all.maxY = std::max(all.maxY, box.maxY);
    }

    return all;
}

double Layout::side() const {
    if (poses.empty()) return 0.0;

    const Box all = bounds();
    return std::max(all.maxX - all.minX, all.maxY - all.minY);
}

std::size_t Layout::bestTreeToRemove() const {
    return removalOrder().front();
}

std::vector<std::size_t> Layout::removalOrder() const {
    std::vector<double> sideWithout(size());

    for (std::size_t skip = 0; skip < size(); ++skip) {
        Box rest{kInf, kInf, -kInf, -kInf};
        for (std::size_t j = 0; j < size(); ++j) {
            if (j == skip) continue;
            rest.minX = std::min(rest.minX, boxes[j].minX);
            rest.minY = std::min(rest.minY, boxes[j].minY);
            rest.maxX = std::max(rest.maxX, boxes[j].maxX);
            rest.maxY = std::max(rest.maxY, boxes[j].maxY);
        }

        sideWithout[skip] = std::max(rest.maxX - rest.minX, rest.maxY - rest.minY);
    }

    std::vector<std::size_t> order(size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t a, std::size_t b) { return sideWithout[a] < sideWithout[b]; });

    return order;
}

std::vector<std::size_t> Layout::edgeTrees() const {
    constexpr double kTol = 1e-9;
    const Box all = bounds();
    std::vector<std::size_t> edge;

    for (std::size_t i = 0; i < size(); ++i) {
        const Box &box = boxes[i];
        const bool onLeft = box.minX <= all.minX + kTol;
        const bool onBottom = box.minY <= all.minY + kTol;
        const bool onRight = box.maxX >= all.maxX - kTol;
        const bool onTop = box.maxY >= all.maxY - kTol;

        if (onLeft || onBottom || onRight || onTop) {
            edge.push_back(i);
        }
    }

    return edge;
}

double totalScore(const std::vector<Layout> &layouts) {
    double total = 0;

    for (std::size_t n = 1; n < layouts.size(); ++n) {
        const double side = layouts[n].side();
        total += side * side / static_cast<double>(n);
    }

    return total;
}
} // namespace packer::packing
