
#include "principalcurve.hpp"

#include <cmath>
#include "transform2d.hpp"
#include "propagate2d.hpp"
#include "plotanalysis.hpp"
#include "curve2d.hpp"
#include "draw.hpp"
#include "dataset.hpp"
#include "propagate2d.hpp"
#include <sstream>
#include <iostream>

void DrawVoronoi(int iter,  Image2D<int> labels, std::vector<Int2>& points, const Image2D<unsigned>& plot){
    DrawContext dc { labels.width(), labels.height(), 600, 600 };
    dc.drawVoronoiRegions(labels, plot, points.size());


    std::stringstream ss;
    ss << "voronoi_" << iter << ".png";
    //dc.setColor(1.0,1.0,0.0);
    //dc.drawMainLine(points);
    dc.setColor(0.0, 1.0,1.0);
    dc.drawControlPoints(points);
    dc.writeToFile(ss.str());
}


//Original HS Approach
int principalCurveHS(std::vector<Int2>& points, const DataSet& dataSet, const Projection2D& projection, int sz,const PipelineParameters& p) {

    Image2D<unsigned> plot = scatterPlot(dataSet, projection, sz,sz);
    auto plotStat = regionStatOfWhole(plot);
    points = principalComponentLine(plotStat.pointStat.average, plotStat.covarianceMatrix(), p.width, p.height, p.segment);



}


/*
    Points... Int -> Pixel Coordinates in the longest path
    Plot  ... Image Containing the Density of the Points per Pixel
    Iteration Count... max number of iterations
    convergence Distance ...
    arcSegment>
    smoothCount>
    eachIteration > ... Function to apply to curve per iteration
*/


int principalCurve(std::vector<Int2>& points, const Image2D<unsigned>& plot, int iterationCount,
                   float convergenceDistance, const float arcSegment, int smoothCount,
                   const std::function<void(int iteration, const std::vector<Int2>&)> eachIteration, bool debug) {



    auto width = plot.width();
    auto height = plot.height();

    Image2D<int> pointsImage { width, height };
    Image2D<float> pointsDistance { width, height };
    Image2D<char> dummyMask { width, height };
    dummyMask.fill(1);

    std::vector<Float2> smoothedPoints;
    std::vector<Int2> pointsBefore, pointsBefore2;

    bool convergence = false;
    int iteration = 0;
    for(iteration = 0; !convergence && iteration < iterationCount; iteration++) {
        pointsBefore2 = pointsBefore;
        pointsBefore = points;

        pointsImage.fill(-1);
        pointsDistance.fill(INFINITY);

        // First we set the guess of principal curve.. f(x) and label the image
        // -1 no point, otherwise index of point
        // Store the distance to the points in principal curve to the image
        // 0 as it is in the index of the point, INFINITY to initialize as non-calculated distance
        for(int i = 0; i < points.size(); i++) {
            pointsImage(points[i][0], points[i][1]) = i;
            pointsDistance(points[i][0], points[i][1]) = 0;
        }

        //Each one starts first with its own origin... propagate
        //the centerline locations so that it saves where is the
        // origin...
        auto origins = propagate2D(dummyMask, pointsDistance);
        // before it was the coordinates.. now labels has
        // for each location the label it belongs to..
        auto labels = indexImage(pointsImage, origins);
        if (debug)
              DrawVoronoi(iteration,labels,points, plot);


        // Plot is the image(density), labels are where each point belongs to, and the number of labels
        auto regions = regionStats(plot, labels, points.size());
        // That labels basically define the voronoi regions of each point...
        smoothedPoints.resize(0);
        smoothedPoints.reserve(regions.size());

        //smooth scatterplot.. look at neighborhoods... it is not weighted so far

        for(int i = 0; i < regions.size(); i++) {
            Stat<Float2> smoothedStat;
            smoothedStat.count = 0;
            for(int j = i - smoothCount; j <= i + smoothCount; j++) {
                smoothedStat = smoothedStat + regions[std::max(0, std::min(j, (int) regions.size() - 1))].pointStat;
            }

            if(smoothedStat.count > 0) {
                //Basically get the centroid of the points in the voronoi regions alongside and
                // then move the point
                smoothedPoints.emplace_back((smoothedStat.average + points[i]) / 2);
            } else {
                //the regions contain no data points and therefore have no average
                //reuse old point
                smoothedPoints.emplace_back(points[i]);
            }
        }


        // Whether we should we sample or not
        if(arcSegment != 0) {
            //extract points to based on arc length
            points.resize(0);
            //not using arc = 0 gives the beginning of the curve "mobility"
            //otherwise it can get stuck in zero-data regions
            float arc = 0;

            for(int i = 1; i < smoothedPoints.size(); i++) {

                auto move = smoothedPoints[i-1] - smoothedPoints[i];
                auto length = sqrt(vectorLengthSquared(move));

                arc += length;

                while(arc >= 0 && length > 0) {
                    auto arcPoint = smoothedPoints[i] + (move * (arc / length));
                    //round to int
                    points.emplace_back(int(round(arcPoint[0])), int(round(arcPoint[1])));
                    arc -= arcSegment;
                }
            }

        } else {
            for(int i = 0; i < smoothedPoints.size(); i++) {
                points[i] = Int2({int(round(smoothedPoints[i][0])), int(round(smoothedPoints[i][1]))});
            }
        }

        if(eachIteration != nullptr) {
            eachIteration(iteration + 1, points);
        }

        //convergent lines should have a fixed amount of points
        //regardless of resampling

        if(points == pointsBefore || points == pointsBefore2) {
            convergence = true;
        } else if(convergenceDistance > 0 && points.size() == pointsBefore.size()) {
            convergence = true;

            for(int i = 0; i < points.size(); i++) {
                auto diff = points[i] - pointsBefore[i];
                if(vectorLengthSquared(diff) > convergenceDistance * convergenceDistance) {
                    convergence = false;
                    break;
                }
            }
        }
    }

    return convergence ? iteration : -1;
}
