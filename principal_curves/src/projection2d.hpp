#ifndef PROJECTION2D_HPP
#define PROJECTION2D_HPP

#include <utility>
#include <functional>
#include <vector>

#include "dataset.hpp"
#include "image2d.hpp"
#include "multival.hpp"

using Projection2D = std::function<Float2(const DataSet&, size_t)>;

class SimpleAxisSelection2D {
    size_t cx, cy;

    public:

    SimpleAxisSelection2D(size_t x, size_t y);

    Float2 operator()(const DataSet& dataSet, size_t index);
};

class AxisSelection2D {
    size_t cx, cy;


    public:

    AxisSelection2D(size_t x, size_t y);

    size_t GetFirstAxis(){return cx;}
    size_t GetSecondAxis(){return cy;}

    Float2 operator()(const DataSet& dataSet, size_t index);
};

class ProjectionMatrix2D {
    std::vector<Float2> entries;

    public:

    ProjectionMatrix2D(size_t n);

    Float2& entry(size_t index);

    Float2 entry(size_t index) const;

    size_t dimension() const;
};

struct CalibratedMatrix2D {
    ProjectionMatrix2D *matrix;
    Float2 min, max;

    public:

    Float2 operator()(const DataSet& dataSet, size_t index);
};

CalibratedMatrix2D calibrateMatrix(ProjectionMatrix2D *matrix, const DataSet& dataSet);


Image2D<unsigned> scatterPlot(const DataSet& dataSet, Projection2D projection, size_t width, size_t height);

Image2D<unsigned> scatterPlotWithAccumulator(const DataSet& dataSet, Projection2D projection, size_t width, size_t height, Image2D<std::vector<std::pair<float, float> > *> *accumulator);

 Image2D<float> gaussedPlot(const DataSet& dataSet, int firstAxis, int secondAxis , size_t width, size_t height, float sigma);

Image2D<float> gaussedPlotWithItk(DataSet dataSet, int firstAxis, int secondAxis , size_t width, size_t height, float sigma, float* sum);

#endif
