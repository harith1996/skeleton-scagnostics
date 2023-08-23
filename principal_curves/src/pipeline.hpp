#ifndef PIPELINE_HPP
#define PIPELINE_HPP

#include <vector>

#include "curve2d.hpp"
#include "dataset.hpp"
#include "projection2d.hpp"
#include "image2d.hpp"
#include "skeletongenerator.h"

struct PipelineParameters {
    int width = 64;
    int height = 64;
    int segment = 6;
    int smallPlotFactor = 4;
    float skeletonThreshold = 2.0f;
    int maxIterations = 20;
    float convergenceDistance = 1.5;
    int smoothCount = 2;
    int smoothCountVariance = 0;
    bool resample = false;
    bool asGraph = false;
    bool usePrincipalComponent = false;
    bool useGPUSkeletonization = false;
    float threshold = 10.0;
    float sigma = 3.0;
    float dominationWeight = 0.2f;
    float dtThresholdWeight = 0.5f;
    signed int firstAxis;
    signed int secondAxis;
};


struct ComponentResult {

    std::vector<AugmentedPoint> augmentedCurve;
    std::vector<SkeletonNode*> graphNodes;
    Image2D<unsigned> plot {0, 0};
    Image2D<unsigned> maskedPlot {0, 0};
    Image2D< std::vector<std::pair<float,float>>* > accumulator{0,0};
    Image2D<char> skeleton {0, 0};
    Image2D<char> thresholded {0,0};
    Image2D<float> dt {0,0};
    vector<Int2> controlPoints;
    bool hasSkeleton;
    int iterations;
    bool converged;
    int totalPoints;

};




struct PipelineResult2 {
    Image2D<unsigned> plot {0, 0};
    Image2D<char> thresholded {0,0};
    Image2D< std::vector<std::pair<float,float>>* > accumulator{0,0};
    Image2D<float> plotGauss {0, 0};
    Image2D<float> normalizedPlotGauss {0, 0};

    Image2D<short> components {0, 0};
    vector<Int2> controlPoints;
    double percentageFilled;
    int numberOfComponents;
    int totalNodes;
    int totalAnalysedPoints;
    std::vector<ComponentResult*> results;
};


struct PipelineResult {
    std::vector<AugmentedPoint> augmentedCurve;
    std::vector<SkeletonNode*> graphNodes;
    Image2D<unsigned> plot {0, 0};
    Image2D<unsigned> maskedPlot {0, 0};
    Image2D< std::vector<std::pair<float,float>>* > accumulator{0,0};
    Image2D<char> skeleton {0, 0};
    Image2D<char> thresholded {0,0};
    vector<Int2> controlPoints;
    int iterations;
    bool converged;
    Image2D<float> plotGauss {0, 0};
    Image2D<float> nonNormalizedPlotGauss {0, 0};

    int numberOfComponents;
    int totalNodes;
    std::vector<ComponentResult> results;
};

PipelineResult processPipeline(const PipelineParameters& p, const DataSet& dataSet, const Projection2D& projection, SkeletonGenerator* generator, double* time, bool debug = false);

PipelineResult processPipelineOriginal(const PipelineParameters& p, const DataSet& dataSet, const Projection2D& projection, bool debug=false);

PipelineResult processPipelineTimeWrapper(const PipelineParameters& p, const DataSet& dataSet, const Projection2D& projection, double* completeTime, bool debug=false);

PipelineResult processGPUPipelineTimeWrapper(const PipelineParameters& p, const DataSet& dataSet, const Projection2D& projection, double* completeTime,
                                             SkeletonGenerator* generator, bool debug=false) ;

PipelineResult processGPUPipelineGraphTimeWrapper(const PipelineParameters& p, const DataSet& dataSet, const Projection2D& projection, double* completeTime,
                                             SkeletonGenerator* generator, int threshold, bool debug);
bool processPipelineGraphTimeWrapped(const PipelineParameters& p, const DataSet& dataSet, const Projection2D& projection, SkeletonGenerator* generator, double *completeTime, bool debug, PipelineResult2 *result);

void GetFullRegions(const PipelineParameters& p, ComponentResult* currentResult, int currentComponent);

#endif
