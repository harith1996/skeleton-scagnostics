#ifndef PRINCIPALCURVE_HPP
#define PRINCIPALCURVE_HPP

#include <vector>
#include <functional>

#include "multival.hpp"
#include "image2d.hpp"
#include "dataset.hpp"
#include "pipeline.hpp"



int principalCurve(std::vector<Int2>& points, const Image2D<unsigned>& plot, int iterationCount, float convergenceDistance, const float arcSegment, int smoothCount,
                   const std::function<void(int iteration, const std::vector<Int2>&)> eachIteration, bool debug=false);


int principalCurveHS(std::vector<Int2>& points, const DataSet &dataSet, const Projection2D &projection, int sz,const PipelineParameters& p);

#endif
