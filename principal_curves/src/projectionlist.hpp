#ifndef PROJECTIONLIST_HPP
#define PROJECTIONLIST_HPP

#include <vector>

#include "projection2d.hpp"
#include "curve2d.hpp"

struct ProjectionEntry {
    ProjectionMatrix2D projectionMatrix { 0 };
    std::vector<AugmentedPoint> curve;
};

using ProjectionList = std::vector<ProjectionEntry>;

#endif
