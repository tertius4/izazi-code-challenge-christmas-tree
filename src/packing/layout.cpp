#include "packing/layout.h"
#include "geometry/collision.h"
#include <algorithm>

namespace packer::packing {
void Layout::add(const Pose &pose) {
    poses.push_back(pose);
    placed.push_back(geometry::transform(shape, pose));
    boxes.push_back(geometry::boundingBox(placed.back()));
}

void Layout::setPose(std::size_t i, const Pose &pose) {
    poses[i] = pose;
    placed[i] = geometry::transform(shape, pose);
    boxes[i] = geometry::boundingBox(placed[i]);
}

void Layout::remove(std::size_t i) {
    poses.erase(poses.begin() + i);
    placed.erase(placed.begin() + i);
    boxes.erase(boxes.begin() + i);
}

bool Layout::collidesWithOthers(std::size_t i) const {
    for (std::size_t j = 0; j < size(); ++j) {
        if (j != i && geometry::polygonsCollide(placed[i], boxes[i], placed[j], boxes[j])) {
            return true;
        }
    }
    return false;
}

bool Layout::isValid() const {
    for (std::size_t i = 0; i < size(); ++i) {
        for (std::size_t j = i + 1; j < size(); ++j) {
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
} // namespace packer::packing