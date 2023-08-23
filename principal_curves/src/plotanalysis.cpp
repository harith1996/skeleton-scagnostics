
#include "plotanalysis.hpp"
#include "math.h"
#include <iostream>

float CovarianceMatrix::varianceAlong(Float2 v) {
    return v[0] * (v[0] * variance[0] + v[1] * covariance) + v[1] * (v[0] * covariance + v[1] * variance[1]);
}

bool RegionStat::empty() const {
    return pointStat.empty();
}

void RegionStat::put(Float2 f) {
    pointStat.put(f);
    productStat.put(f[0] * f[1]);
}

void RegionStat::put(size_t n, Float2 f) {
    pointStat.put(n, f);
    productStat.put(n, f[0] * f[1]);
}

float RegionStat::covariance() const {
    return productStat.average - pointStat.average[0] * pointStat.average[1];
}

RegionStat operator+(RegionStat a, RegionStat b) {
    RegionStat result;
    result.pointStat = a.pointStat + b.pointStat;
    result.productStat = a.productStat + b.productStat;
    return result;
}

CovarianceMatrix RegionStat::covarianceMatrix() const {
    CovarianceMatrix matrix;
    matrix.variance = pointStat.variance();
    matrix.covariance = covariance();
    return matrix;
}

std::vector<RegionStat> regionStats(const Image2D<unsigned>& plot, const Image2D<int>& label, int labelCount) {
    auto w = plot.width();
    auto h = plot.height();

    std::vector<RegionStat> result;

    result.resize(labelCount);

    for(int x = 0; x < w; x++) {
        for(int y = 0; y < h; y++) {
            auto l = label(x, y); // Label should be one of the existent...
            // it may be no label because of mask (?)

            if(l < 0 || l >= labelCount) continue;

            result[l].put(plot(x, y), {x, y});
        }
    }

    return result;
}

std::vector<RegionStat> regionStatsFull(const  Image2D< std::vector<std::pair<float,float>>* >& plot, const Image2D<int>& label, int labelCount) {
    auto w = plot.width();
    auto h = plot.height();

    std::vector<RegionStat> result;

    result.resize(labelCount);

    for(int x = 0; x < w; x++) {
        for(int y = 0; y < h; y++) {
            auto l = label(x, y); // Label should be one of the existent...
            // it may be no label because of mask (?)

            if(l < 0 || l >= labelCount) continue;

            if ( plot(x,y) == nullptr ) continue;

            auto vectorOfPoints = plot(x,y);

            for(int k =0; k < vectorOfPoints->size(); k++){
                auto currentPoint = vectorOfPoints->at(k);
                result[l].put({currentPoint.first, currentPoint.second});
                //result[l].put(plot(x, y), {x, y});
            }
        }
    }

    return result;
}


RegionStat regionStatOfWhole(const Image2D<unsigned>& plot) {
    auto w = plot.width();
    auto h = plot.height();

    RegionStat result;

    for(int x = 0; x < w; x++) {
        for(int y = 0; y < h; y++) {
            result.put(plot(x, y), {x, y});
        }
    }

    return result;
}

std::vector<Int2> principalComponentLine(Float2 centroid, CovarianceMatrix matrix, int width, int height, float step) {
    //do power iteration to find maximum-eigenvalue-eigenvector
    //symmetric positive-semidefinite matrix impliciteyl given by (variance[0], covariance; covariance, variance[1])

    auto covariance = matrix.covariance;
    auto variance = matrix.variance;

    Float2 v {1, 1};

    for(int i = 0; i < 100; i++) {
        v = Float2({v[0] * variance[0] + v[1] * covariance, v[0] * covariance + v[1] * variance[1]});

        //normalize
        float len = std::sqrt(vectorLengthSquared(v));
        v /= len;
    }

    //use step length
    v *= step;

    //move to border in one direction
    auto position = centroid;
    while(position[0] >= 0 && position[0] < width && position[1] >= 0 & position[1] < height) {
        position -= v;
    }

    //move back into area
    position += v;

    //gather points along line into result
    std::vector<Int2> result;
    while(position[0] >= 0 && position[0] < width && position[1] >= 0 & position[1] < height) {
        position += v;

        result.emplace_back((int)std::round(position[0]), (int)std::round(position[1]));
    }


    return result;
}

float averageDistance(const Image2D<unsigned>& plot) {
    auto w = plot.width();
    auto h = plot.height();

    Stat<float> stat;

    for(int x = 0; x < w; x++) {
        for(int y = 0; y < h; y++) {
            if(!plot(x, y)) continue;

            for(int x2 = 0; x2 < w; x2++) {
                for(int y2 = 0; y2 < h; y2++) {
                    if(!plot(x2, y2)) continue;

                    stat.put(plot(x, y) * plot(x2, y2), sqrt((x - x2) * (x - x2) + (y - y2) * (y - y2)));
                }
            }

        }
    }

    return stat.average;
}
