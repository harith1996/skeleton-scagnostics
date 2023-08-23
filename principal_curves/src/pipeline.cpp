#include "pipeline.hpp"

#include "transform2d.hpp"
#include "skeleton.hpp"
#include "multival.hpp"
#include "principalcurve.hpp"
#include "curve2d.hpp"
#include "plotanalysis.hpp"
#include "skeletongenerator.h"
#include "draw.hpp"
#include "principalgraph.hpp"
#include "propagate2d.hpp"

#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>

#include <chrono>
#include <iostream>
#include <queue>

PipelineResult processPipelineTimeWrapper(const PipelineParameters& p, const DataSet& dataSet, const Projection2D& projection, double* completeTime, bool debug) {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    auto res = processPipelineOriginal(p, dataSet, projection, debug);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    *completeTime = (std::chrono::duration_cast<std::chrono::microseconds>(end - start).count())/1000.0f; // micro to ms
    return res;
}


PipelineResult processGPUPipelineTimeWrapper(const PipelineParameters& p, const DataSet& dataSet, const Projection2D& projection, double* completeTime,
                                             SkeletonGenerator* generator, bool debug){

    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    PipelineResult result;
    std::vector<Int2> points;

    auto plot = generator->CreateScatterplot(dataSet, projection, p.width,p.height);
    //Image2D<unsigned> plot = scatterPlot(dataSet, projection, p.width, p.height);

    Image2D<char> skeleton {plot.width(), plot.height()};
    // Get as well the mask of the plot... and use that instead for the plot
    // TODO
    Image2D<unsigned> maskedPlot{plot.width(), plot.height()};
    vector<Int2>  skeletonPoints;

    generator->GetSkeleton(plot, p.width, skeleton, maskedPlot,&skeletonPoints,4);

    result.skeleton = skeleton;
    float skeletonCurveLength = 0;
    int numberChecked = 0;
    for(int x = 0; x < p.width; x++) {
        for(int y = 0; y < p.height; y++) {

            // Instead of calculating everything... first calculate whether it should be explored
            if(skeleton(x, y)) {
                // check the neighbors
                if ( numberOfNeighbors(skeleton, {x,y} ) <= 2){
                    tryLongestPathFrom(skeleton, {x, y}, skeletonCurveLength, points, p.segment);
                    numberChecked++;
                }
            }
        }
    }


    result.iterations = principalCurve(points, maskedPlot, p.maxIterations, p.convergenceDistance, p.resample ? p.segment : 0, p.smoothCount, nullptr, debug);

    if(result.iterations != -1) {
        result.converged = true;
    } else {
        result.iterations = p.maxIterations;
        result.converged = false;
    }
    // and finally augment the values

    result.augmentedCurve = augmentCurve(points, plot, p.smoothCountVariance);
    result.plot = plot;

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    *completeTime = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()/1000.0f;
    return result;
}


PipelineResult processGPUPipelineGraphTimeWrapper(const PipelineParameters& p, const DataSet& dataSet, const Projection2D& projection, double* completeTime,
                                             SkeletonGenerator* generator, int threshold, bool debug){

    std::cout << "GPU Pipeline" << std::endl;
    PipelineResult result;

    auto plot = generator->CreateScatterplot(dataSet, projection, p.width,p.height);
    //Image2D<unsigned> plot = scatterPlot(dataSet, projection, p.width, p.height);

    result.plot = plot;
    Image2D<char> skeleton {plot.width(), plot.height()};
    // Get as well the mask of the plot... and use that instead for the plot
    // TODO
    Image2D<unsigned> maskedPlot{plot.width(), plot.height()};
    Image2D<char> newSkeleton{static_cast<unsigned long>(p.width), static_cast<unsigned long>(p.height)};


    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    vector<Int2> tmpSkeletonPoints, skeletonPoints;

    generator->GetSkeleton(plot, p.width, skeleton, maskedPlot, &tmpSkeletonPoints, threshold);


    if ( tmpSkeletonPoints.empty()){
        return result;
    }
    result.skeleton = skeleton;
    result.maskedPlot = maskedPlot;
    generator->CreateOnePixelWidthSkeleton(&tmpSkeletonPoints, &skeletonPoints, plot.width());

    vector<int> endPoints, bifurcations;
    generator->GetEndPointsAndBifurcations(&skeletonPoints, &endPoints, &bifurcations);

    CreateSkeletonGraph(p.segment, endPoints, bifurcations, &skeletonPoints, &(result.graphNodes));

    tmpSkeletonPoints.clear();
    for(unsigned int j = 0; j< result.graphNodes.size(); j++)
        tmpSkeletonPoints.push_back(result.graphNodes.at(j)->location);

    result.controlPoints = tmpSkeletonPoints;
    // Now from the skeleton do the graph...

    result.iterations = principalGraph(&(result.graphNodes), plot, p.maxIterations, p.convergenceDistance,p.resample ? p.segment : 0, p.smoothCount, debug);


    if(result.iterations != -1) {
        result.converged = true;
    } else {
        result.iterations = p.maxIterations;
        result.converged = false;
    }
    // and finally augment the values
    augmentGraph(&(result.graphNodes), plot, p.smoothCountVariance);

    result.plot = plot;

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    *completeTime = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()/1000.0f;
    return result;
}


double Gaussian2DDistributionMaxV( double stddev){
      // mean 0 and std dev as set
      double inverse = 1.0 /( 2.0 * 3.14159 * stddev* stddev);
      double diff =  (pow(0,2.0) + pow(0,2.0)) / (2.0* stddev * stddev);
      return inverse * exp( -diff );
}


bool processPipelineGraphTimeWrapped(const PipelineParameters& p, const DataSet& dataSet, const Projection2D& projection, SkeletonGenerator* generator,
                                     double* completeTime, bool debug, PipelineResult2* result) {


    //First step is to create the compact shape. It might generate
    //multiple shapes, therefore we analyse each component individually
    Image2D<unsigned> maskedplot { static_cast<unsigned long>(p.width), static_cast<unsigned long>(p.height)};
    Image2D<short> componentsPlot { static_cast<unsigned long>(p.width), static_cast<unsigned long>(p.height)};

    //Image2D<unsigned> smallPlot = scatterPlot(dataSet, projection, p.width / p.smallPlotFactor, p.height / p.smallPlotFactor);
    Image2D<char> newSkeleton{static_cast<unsigned long>(p.width), static_cast<unsigned long>(p.height)};
    Image2D<float> blob_grey {static_cast<unsigned long>(p.width), static_cast<unsigned long>(p.height) };

    std::vector<Int2> points;

    //Previous Version... the sigma is based on the resolution
    //It will affect the resolution...

    Image2D< std::vector<std::pair<float,float>>* > accumulator{static_cast<size_t>(p.width),
                                                                static_cast<size_t>(p.height)};

    Image2D<unsigned> plot = scatterPlotWithAccumulator(dataSet, projection, p.width, p.height,&accumulator);
    //The plot itself stores them in truncated locations ...
    result->plot = plot;
    //double sigma = std::max(1.0f, std::min(p.width / 25.0f, 10.0f));
    //blob_grey = gaussianBlob(plot, p.sigma);
    //result.plotGauss = blob_grey;
    //threshold by average
    //auto blob = threshold(blob_grey, density);
    //result.thresholded = blob;
    result->accumulator = accumulator;
    //********
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    //result->plotGauss = gaussedPlot(dataSet, p.firstAxis,p.secondAxis, p.width, p.height, p.sigma);
    //std::cout << "Sigma to be used " << p.sigma <<std::endl;
    //std::cout << "Spacing " << 1.0f/p.width << std::endl;
    float spacing = 1.0f/p.width;
    //std::cout << "Dev in pixels " <<  p.sigma /(spacing);

    float sum= 0;



   result->plotGauss  = gaussedPlotWithItk(dataSet, p.firstAxis,p.secondAxis, p.width, p.height, p.sigma ,&sum);
   result->normalizedPlotGauss = normalizeImageValues(result->plotGauss);

   // const float density = dataSet.size() * 1.0f / (p.width * p.height);

     float density2= (1.0f*dataSet.size())/(p.width * p.height);

    // std::cout << "Max Computed  V " <<  Gaussian2DDistributionMaxV(p.sigma) << std::endl;
    // std::cout << "Max v?" << sum << std::endl;

     //std::cout << "First thresholding " << sum << "..." << density2 << std::endl;

    float alpha =  (p.sigma /spacing);//*(p.sigma/spacing);//ceil((p.sigma /(spacing)));
    // std::cout << "?? max? " << sum << std::endl;



    float computedThreshold = density2*alpha;
    double  percentageFilled = 0;
    result->thresholded = threshold(result->plotGauss, computedThreshold,&percentageFilled);
    auto blob = result->thresholded;
    result->percentageFilled = percentageFilled;
    //Let's find the number of connected components...
    //for each we can apply the skeletonization


    vector<int> sizes;
    int totalObjects =  GetNumberOfComponents(blob, componentsPlot, &sizes);

    // The blob should be first modified....


    result->components = componentsPlot;
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
   // std::cout << "Compute components  " << (std::chrono::duration_cast<std::chrono::microseconds>(end - start).count())/1000.0f << std::endl;


    // std::cout << "Total objects? " << totalObjects << std::endl << std::endl;
    result->totalNodes = 0;
    result->totalAnalysedPoints = 0;
    result->numberOfComponents = totalObjects;
    // The blob should be first modified....trying to get the sizes...

   // std::cout << "Original :" ;


    for(int currentlyAnalysingBlob  = 1; currentlyAnalysingBlob <= totalObjects; currentlyAnalysingBlob++ ){
        ComponentResult* newComponentResult = new ComponentResult;

        int sizeBlob = 0;

        Image2D< std::vector<std::pair<float,float>>* > resAccumulator{static_cast<size_t>(p.width), static_cast<size_t>(p.height)};
        newComponentResult->accumulator = resAccumulator;


        Image2D<unsigned> resPlot{static_cast<size_t>(p.width), static_cast<size_t>(p.height)};
        Image2D<char> blob2{static_cast<size_t>(p.width), static_cast<size_t>(p.height)};

        int rX[2] = {p.width,0};
        int rY[2] = {p.height,0};

        //Getting a single component ...
          for(int x = 0; x < p.width; x++) {
           for(int y = 0; y < p.height; y++) {
               if (componentsPlot(x,y) != currentlyAnalysingBlob){
                  blob2(x,y) = 0;
                  newComponentResult->accumulator(x,y) = nullptr; //new std::vector<std::pair<float,float>>();
                  resPlot(x,y) = 0;

               }
               else {
                   sizeBlob++;
                   blob2(x,y) = 1;

                   if ( x < rX[0] ) rX[0] = x;
                   if ( x > rX[1] ) rX[1] = x;

                   if ( y < rY[0] ) rY[0] = y;
                   if ( y > rY[1] ) rY[1] = y;

                   newComponentResult->accumulator(x,y) = accumulator(x,y);
                   resPlot(x,y) = plot(x,y);
               }
           }
        }

        newComponentResult->plot = resPlot;
        newComponentResult->thresholded = blob2;

        // std::cout << " Size Blob " <<  sizeBlob << std::endl;
        //int minR = min((rX[1]-rX[0]),  (rY[1]-rY[0])  );

        //Here is actually the threshold..... 0.03 times the size of blob
        //float thres = sizes.at(currentlyAnalysingBlob-1)*0.02;

        float maxDt = 0;

        //.... The blob2 is correct

        vector<int> endPoints, bifurcations;
        auto skeleton_grey = computeSkeleton2_Grad(blob2,&maxDt);
        //.....dt(?)
        //use the threshold based on the maximum distance ...

        newComponentResult->dt = skeleton_grey;

        double percentageFilled = 0;
        float thres = maxDt*p.dtThresholdWeight;
        auto skeleton = threshold(skeleton_grey, thres,&percentageFilled);
        newComponentResult->skeleton = skeleton;

        vector<Int2> tmpSkeletonPoints, skeletonPoints;
        for(int x = 0; x < p.width; x++) {
           for(int y = 0; y < p.height; y++) {
           // Instead of calculating everything... first calculate whether it should be explored
              if(skeleton(x, y)) {
                   Int2 pt{x,y};
                   tmpSkeletonPoints.push_back(pt);

               }
               if ( blob2(x,y) ==1 ){
                   maskedplot(x,y) = plot(x,y);
               }
               else  maskedplot(x,y) = 0;
           }
        }



        //std::cout << "Original skeleton points? " << tmpSkeletonPoints.size() << std::endl;
        //std::cout << "Thres? " << thres <<  " ... maxDt " << maxDt <<  std::endl;

        //****************
        // Simplifying the skeleton...
        vector<Int2>*a,*b,*c;
        a = &tmpSkeletonPoints; b = &skeletonPoints;
        int iters = 0;
        while(a->size() != b->size()){
            generator->CreateOnePixelWidthSkeleton(a, b, plot.width());
            //a is used to generate b...
            c = a;
            a = b;
            b = c;
            iters ++;
        }

         if ((&skeletonPoints) != b){
              skeletonPoints.clear();
              for(int k = 0; k < b->size(); k++){
                  skeletonPoints.push_back( b->at(k));
              }
         }


         generator->GetEndPointsAndBifurcations(&skeletonPoints, &endPoints, &bifurcations);
         for(int j = 0; j < skeletonPoints.size(); j++){
               newSkeleton(skeletonPoints.at(j)[0],skeletonPoints.at(j)[1]) = 1;
         }
         //*****************************


         newComponentResult->maskedPlot = maskedplot;

         if ( !skeletonPoints.empty()){

              //   std::cout << "Skeleton points > " << skeletonPoints.size() << "\t";
              //   std::cout << "End Points >" <<  endPoints.size() << "\tBifurcations >" <<  bifurcations.size() <<  std::endl;

             // int controlPointSamplingSize = p.segment;
                  //This... is the function... that needs to be changed?




                 CreateSkeletonGraphAutomaticPoints(endPoints,bifurcations, &skeletonPoints,
                                                  &(newComponentResult->graphNodes), skeleton_grey,p.dominationWeight);

                 /*CreateSkeletonGraph(controlPointSamplingSize, endPoints,
                                      bifurcations, &skeletonPoints,
                                      &(newComponentResult->graphNodes));*/


                 //std::cout << "Before: Nodes in graph " << newComponentResult->graphNodes.size() << std::endl;

                 // Once this is done, use it to calculate the principal curve ....
                  newComponentResult->iterations = principalGraph(&(newComponentResult->graphNodes), maskedplot,
                                                                  p.maxIterations, p.convergenceDistance,
                                                                  p.resample ? p.segment : 0, p.smoothCount, debug);

                 //std::cout << "After: Nodes in graph " << newComponentResult->graphNodes.size() << std::endl;

                  GetFullRegions(p, newComponentResult, currentlyAnalysingBlob);

                  tmpSkeletonPoints.clear();
                  for(unsigned int j = 0; j< newComponentResult->graphNodes.size(); j++){
                      tmpSkeletonPoints.push_back(newComponentResult->graphNodes.at(j)->location);
                      result->controlPoints.push_back(newComponentResult->graphNodes.at(j)->location);
                  }
                  newComponentResult->controlPoints = tmpSkeletonPoints;

                  if(newComponentResult->iterations != -1) {
                     newComponentResult->converged = true;
                  } else {
                      newComponentResult->iterations = p.maxIterations;
                      newComponentResult->converged = false;
                  }
                  augmentGraph(&(newComponentResult->graphNodes), resPlot, p.smoothCountVariance);
                  newComponentResult->hasSkeleton = true;

                  int totPts = 0;
                  for(int k =0;k < newComponentResult->graphNodes.size(); k++){
                      totPts += newComponentResult->graphNodes.at(k)->stat.pointStat.count;
                  }
                  newComponentResult->totalPoints = totPts;
                  result->totalAnalysedPoints += totPts;
                  result->totalNodes += newComponentResult->graphNodes.size();



          }
          else {
             newComponentResult->hasSkeleton = false;
          }


         //std::cout << std::endl << std::endl;
          result->results.push_back(newComponentResult);
    }


   // std::cout << std::endl;
    result->plot = plot;


    end = std::chrono::steady_clock::now();
    //std::cout << "Total time  " << (std::chrono::duration_cast<std::chrono::microseconds>(end - start).count())/1000.0f << std::endl;
    *completeTime = (std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    //return result;
    return true;
}

void GetFullRegions(const PipelineParameters& p, ComponentResult* currentResult, int currentComponent){
    auto width = p.width;
    auto height = p.height;

    Image2D<int> pointsImage { static_cast<size_t>(width), static_cast<size_t>(height) };
    Image2D<float> pointsDistance  { static_cast<size_t>(width), static_cast<size_t>(height) };
    Image2D<char> dummyMask  { static_cast<size_t>(width), static_cast<size_t>(height) };

    pointsImage.fill(-1);
    pointsDistance.fill(INFINITY);
    dummyMask.fill(1);

    auto nodes = &(currentResult->graphNodes);
    for(int i = 0; i < nodes->size(); i++) {
        auto pt = nodes->at(i)->location;
        pointsImage(pt[0], pt[1]) = i;
        pointsDistance(pt[0], pt[1]) = 0;
    }

    auto origins = propagate2D(dummyMask, pointsDistance);
    auto labels = indexImage(pointsImage, origins);

    auto regions = regionStatsFull( currentResult->accumulator, labels, nodes->size());

    for(int i = 0; i < nodes->size(); i++) {
        nodes->at(i)->stat2 = regions.at(i);
    }


}

PipelineResult processPipelineOriginal(const PipelineParameters& p, const DataSet& dataSet, const Projection2D& projection, bool debug) {
    PipelineResult result;

    Image2D<unsigned> plot = scatterPlot(dataSet, projection, p.width, p.height);
    //Image2D<unsigned> smallPlot = scatterPlot(dataSet, projection, p.width / p.smallPlotFactor, p.height / p.smallPlotFactor);

    std::vector<Int2> points;

    // In order to process the pipeline. .. first we calculate either the
    // principal component line... or define the longest path
    //
    if(!p.usePrincipalComponent && !p.useGPUSkeletonization) {
        Image2D<float> blob_grey { static_cast<size_t>(p.width), p.height };
        /*
        auto sigma = p.smallPlotFactor * averageDistance(smallPlot);
        blob_grey = gaussianBlob(plot, std::max(1.0f, std::min(sigma / 4, 10.0f)));
        */
        blob_grey = gaussianBlob(plot, std::max(1.0f, std::min(p.width / 25.0f, 10.0f)));
        result.plotGauss = blob_grey;
        //threshold by average
        const float density = dataSet.size() * 1.0f / (p.width * p.height);
        double percentageFilled = 0;

        auto blob = threshold(blob_grey, density, &percentageFilled);

        float maxDt = 0;
        auto skeleton_grey = computeSkeleton2_Grad(blob,&maxDt);

        auto skeleton = threshold(skeleton_grey, 4.0f, &percentageFilled);
        result.skeleton = skeleton;

        float skeletonCurveLength = 0;

        for(int x = 0; x < p.width; x++) {
            for(int y = 0; y < p.height; y++) {

                // Instead of calculating everything... first calculate whether it should be explored
                if(skeleton(x, y)) {
                    // check the neighbors
                    //if ( numberOfNeighbors(skeleton, {x,y} ) <= 2)
                        tryLongestPathFrom(skeleton, {x, y}, skeletonCurveLength, points, p.segment);
                }
            }
        }       
    }
    else if (p.usePrincipalComponent){
        auto plotStat = regionStatOfWhole(plot);
        points = principalComponentLine(plotStat.pointStat.average, plotStat.covarianceMatrix(), p.width, p.height, p.segment);
    }
    // Once this is done, use it to calculate the principal curve ....
    result.iterations = principalCurve(points, plot, p.maxIterations, p.convergenceDistance, p.resample ? p.segment : 0, p.smoothCount, nullptr, debug);

    if(result.iterations != -1) {
        result.converged = true;
    } else {
        result.iterations = p.maxIterations;
        result.converged = false;
    }

    // and finally augment the values

    result.augmentedCurve = augmentCurve(points, plot, p.smoothCountVariance);
    result.plot = plot;

    return result;
}



PipelineResult processPipeline(const PipelineParameters& p, const DataSet& dataSet, const Projection2D& projection, SkeletonGenerator* generator, double* time, bool debug){


    if (p.asGraph){

        if (p.useGPUSkeletonization)
            return processGPUPipelineGraphTimeWrapper(p, dataSet, projection, time, generator, p.threshold , debug);
        //else
        //    return processPipelineGraphTimeWrapped(p, dataSet, projection, generator,time, debug);
    }
    else {

        if (p.useGPUSkeletonization)
            return processGPUPipelineTimeWrapper(p, dataSet, projection, time, generator, debug);
        else
            return processPipelineTimeWrapper(p, dataSet, projection,time, debug);
    }

}


