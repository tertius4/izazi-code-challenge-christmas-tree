#include "packing/layout.h"
#include "geometry/collision.h"
#include <algorithm>
#include <limits>

namespace packer::packing {
void Layout::add(const Pose &pose) {
    poses.push_back(pose);
    placed.push_back(geometry::transform(shape, pose));
    boxes.push_back(geometry::boundingBox(placed.back()));
}

void Layout::setPose(std::size_t i, const Pose &pose) {
    poses[i] = pose;
    geometry::transformInto(shape, pose, placed[i]);
    boxes[i] = geometry::boundingBox(placed[i]);
}

void Layout::remove(std::size_t i) {
    poses.erase(poses.begin() + i);
    placed.erase(placed.begin() + i);
    boxes.erase(boxes.begin() + i);
}

bool Layout::collidesWithOthers(std::size_t i) const {
    const Box &bi = boxes[i];
    for (std::size_t j = 0; j < size(); ++j) {
        // Cheap box test inline; most trees are far away.
        if (j == i || !geometry::boxesOverlap(bi, boxes[j])) continue;
        if (geometry::polygonsCollide(placed[i], bi, placed[j], boxes[j])) return true;
    }
    return false;
}

std::vector<std::size_t> Layout::collidingWith(std::size_t i) const {
    std::vector<std::size_t> hits;
    const Box &bi = boxes[i];
    for (std::size_t j = 0; j < size(); ++j) {
        if (j == i || !geometry::boxesOverlap(bi, boxes[j])) continue;
        if (geometry::polygonsCollide(placed[i], bi, placed[j], boxes[j])) hits.push_back(j);
    }
    return hits;
}

bool Layout::isValid() const {
    for (std::size_t i = 0; i < size(); ++i) {
        for (std::size_t j = i + 1; j < size(); ++j) {
            if (!geometry::boxesOverlap(boxes[i], boxes[j])) continue;
            if (geometry::polygonsCollide(placed[i], boxes[i], placed[j], boxes[j])) return false;
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
        if (t.minX <= b.minX + kTol || t.minY <= b.minY + kTol || t.maxX >= b.maxX - kTol ||
            t.maxY >= b.maxY - kTol) {
            edge.push_back(i);
        }
    }
    return edge;
}
} // namespace packer::packing