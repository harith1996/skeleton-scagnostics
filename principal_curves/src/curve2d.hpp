#ifndef CURVE2D_HPP
#define CURVE2D_HPP

#include <vector>
#include <cmath>
#include <algorithm>

#include "image2d.hpp"
#include "multival.hpp"

template<class T>
T couplingDistanceSquared(const std::vector<MultiVal<T, 2>> curve1, const std::vector<MultiVal<T, 2>> curve2) {
    Image2D<T> table { curve1.size(), curve2.size() };

    for(int x = 0; x < curve1.size(); x++) {
        for(int y = 0; y < curve2.size(); y++) {
            auto diff = curve1[x] - curve2[y];
            auto distSquared = vectorLengthSquared(diff);

            if(x == 0 && y == 0) {
                table(x, y) = distSquared;
            } else if(x == 0) {
                table(x, y) = std::max(distSquared, table(x, y-1));
            } else if(y == 0) {
                table(x, y) = std::max(distSquared, table(x-1, y));
            } else {
                table(x, y) = std::max(distSquared, std::min(table(x-1, y), std::min(table(x, y-1), table(x-1, y-1))) );
            }
        }
    }
    return table(curve1.size() - 1, curve2.size() - 1);
}

template<class T>
float couplingDistance(const std::vector<MultiVal<T, 2>> curve1, const std::vector<MultiVal<T, 2>> curve2) {
    return sqrt(couplingDistanceSquared(curve1, curve2));
}

void tryLongestPathFrom(const Image2D<char>& mask, Int2 start, float& currentMaxLength, std::vector<Int2>& currentMaxPath, int pixPerSegment);

int numberOfNeighbors(const Image2D<char>& mask_original, Int2 loc);

Float2 curveDirectionVector(const std::vector<Int2>& curve, int index);

float curveLocalLength(const std::vector<Int2>& curve, int index);

struct AugmentedPoint {
    Float2 point = 0;
    Float2 direction = 0;
    float amount = 0;
    float amountPerLength = 0;
    float orthoVariance = 0;
};

float couplingDistanceAugmented(const std::vector<AugmentedPoint>& curve1, const std::vector<AugmentedPoint>& curve2, float vf, float af, bool enableMirror = false);

float couplingDistanceAugmentedSquared(const std::vector<AugmentedPoint>& curve1, const std::vector<AugmentedPoint>& curve2, float vf, float af, bool enableMirror = false);

float sharpenssOfAugmentedCurve(const std::vector<AugmentedPoint>& curve);

std::vector<AugmentedPoint> augmentCurve(const std::vector<Int2>& curve, const Image2D<unsigned>& plot, int smoothCount);

#endif
