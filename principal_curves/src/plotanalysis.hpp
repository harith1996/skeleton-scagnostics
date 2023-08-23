#ifndef PLOTANALYSIS_HPP
#define PLOTANALYSIS_HPP

#include <vector>

#include "image2d.hpp"
#include "stat.hpp"
#include "multival.hpp"


struct CovarianceMatrix {
    Float2 variance = 0;
    float covariance = 0;

    float varianceAlong(Float2 v);
};

struct RegionStat {
    Stat<Float2> pointStat;
    Stat<float> productStat;

    bool empty() const;
    float covariance() const;
    CovarianceMatrix covarianceMatrix() const;

    void put(Float2 f);
    void put(size_t n, Float2 f);
};

RegionStat operator+(RegionStat a, RegionStat b);

std::vector<RegionStat> regionStats(const Image2D<unsigned>& plot, const Image2D<int>& label, int labelCount);

std::vector<RegionStat> regionStatsFull(const  Image2D< std::vector<std::pair<float,float>>* >& plot, const Image2D<int>& label, int labelCount);

RegionStat regionStatOfWhole(const Image2D<unsigned>& plot);

std::vector<Int2> principalComponentLine(Float2 centroid, CovarianceMatrix matrix, int width, int height, float step);

float averageDistance(const Image2D<unsigned>& plot);

#endif
