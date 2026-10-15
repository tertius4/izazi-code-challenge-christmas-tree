#include "packing/layout.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace packer::packing {
namespace {
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

    double minX = INFINITY, minY = INFINITY, maxX = -INFINITY, maxY = -INFINITY;
    for (const tree::Placed &t : placed) {
        minX = std::min(minX, t.centre.x);
        minY = std::min(minY, t.centre.y);
        maxX = std::max(maxX, t.centre.x);
        maxY = std::max(maxY, t.centre.y);
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
    const tree::Placed &t = placed[i];
    const int cx = cellOf[i] % cols;
    const int cy = cellOf[i] / cols;

    // Trees close enough to touch tree i sit in its cell or one of the eight around it.
    const int firstRow = std::max(cy - 1, 0);
    const int lastRow = std::min(cy + 1, rows - 1);
    const int firstCol = std::max(cx - 1, 0);
    const int lastCol = std::min(cx + 1, cols - 1);

    for (int y = firstRow; y <= lastRow; ++y) {
        for (int x = firstCol; x <= lastCol; ++x) {
            for (const std::uint32_t j : cells[y * cols + x]) {
                if (j == i) continue;
                if (tree::collide(t, placed[j])) return true;
            }
        }
    }

    return false;
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
    Box b = boxes[0];
    for (const Box &t : boxes) {
        b.minX = std::min(b.minX, t.minX);
        b.minY = std::min(b.minY, t.minY);
        b.maxX = std::max(b.maxX, t.maxX);
        b.maxY = std::max(b.maxY, t.maxY);
    }

    return b;
}

double Layout::side() const {
    if (poses.empty()) return 0.0;

    const Box b = bounds();
    return std::max(b.maxX - b.minX, b.maxY - b.minY);
}

std::size_t Layout::bestTreeToRemove() const {
    constexpr double kInf = std::numeric_limits<double>::infinity();
    std::size_t bestIndex = 0;
    double bestSide = kInf;

    for (std::size_t skip = 0; skip < size(); ++skip) {
        Box b{kInf, kInf, -kInf, -kInf};
        for (std::size_t j = 0; j < size(); ++j) {
            if (j == skip) continue;
            b.minX = std::min(b.minX, boxes[j].minX);
            b.minY = std::min(b.minY, boxes[j].minY);
            b.maxX = std::max(b.maxX, boxes[j].maxX);
            b.maxY = std::max(b.maxY, boxes[j].maxY);
        }

        const double side = std::max(b.maxX - b.minX, b.maxY - b.minY);
        if (side < bestSide) {
            bestSide = side;
            bestIndex = skip;
        }
    }

    return bestIndex;
}

std::vector<std::size_t> Layout::edgeTrees() const {
    constexpr double kTol = 1e-9;
    const Box b = bounds();
    std::vector<std::size_t> edge;

    for (std::size_t i = 0; i < size(); ++i) {
        const Box &t = boxes[i];
        const bool onLeft = t.minX <= b.minX + kTol;
        const bool onBottom = t.minY <= b.minY + kTol;
        const bool onRight = t.maxX >= b.maxX - kTol;
        const bool onTop = t.maxY >= b.maxY - kTol;

        if (onLeft || onBottom || onRight || onTop) {
            edge.push_back(i);
        }
    }

    return edge;
}
} // namespace packer::packing
