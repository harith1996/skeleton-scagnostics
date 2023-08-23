#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <string>
#include <functional>

#include "image2d.hpp"
#include "image2dio.hpp"
#include "dataset.hpp"
#include "parsedata.hpp"
#include "projection2d.hpp"
#include "draw.hpp"
#include "pipeline.hpp"
#include "transform2d.hpp"
#include "curve2d.hpp"
#include "skeletongenerator.h"
#include "principalgraph.hpp"
#include "propagate2d.hpp"
#include "stdlib.h"
#include <vtkFeatureEdges.h>
#include <vtkLine.h>
#include <dirent.h>

#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <stack>
#include <queue>
bool debugAll = false;

#include <cfloat>
#include <limits>

struct TestResults {
    float meanTime;
    float devTime;
    float meanIters;
    float meanDistance;
    float devDistance;

    void PrintResults(bool asTable){

        if (asTable)
            printf("%f\t%f\t%f\n", meanTime, devTime, meanDistance);
         else
            printf("Time: %f\t%f  Iters : %f, Distance: %f,%f\n", meanTime, devTime, meanIters, meanDistance, devDistance);

    }
};

void SaveResultsGraph( PipelineParameters parameters,PipelineResult* res, std::vector<AugmentedPoint> curve, const char* filename){

    {
       std::cout << "Drawing original " << std::endl;
        DrawContext dc { static_cast<int>(res->plot.width()), static_cast<int>(res->plot.height()), 600, 600 };
        dc.setColor(0, 0, 1);
        dc.drawPlot(res->plot);
        //dc.drawPlotBackground(normalizeImageValues(res->plot));
        dc.setColor(1,0,0);
        //dc.drawSkeleton(res->skeleton);

        //dc.drawMainLine(curve);
        //dc.drawControlPoints(curve, 2);
        dc.writeToFile("original.png");

    }



    {
        DrawContext dc { static_cast<int>(res->plot.width()), static_cast<int>(res->plot.height()), 300, 300 };
        //dc.drawNor (res.maskedPlot);
        //
        auto tmp = GetBands(&(res->graphNodes), res->plot.width());
        //dc.drawPlot(res->plot);
        dc.setColor(1, 1, 0);
        //dc.drawPlotBackground(tmp);
        dc.setColor(0.3, 0.3, 1);
        //dc.drawPlotBackground(normalizeImageValues(res->plot));
        dc.drawPlot(res->plot);
        //dc.drawSkeleton(res->thresholded);
        //dc.setColor(0,1,0);
        //dc.drawControlPoints(curve, 2);

        /*
        dc.setColor(1, 0, 0);
        dc.drawVoronoiGraph(res);
        dc.setColor(1,1,0);
        dc.drawControlPoints(res->controlPoints);
        */


       // dc.setColor(1,1,1);
       // dc.drawSkeleton(res->skeleton);
        //dc.setColor(0,1,0);
        //dc.drawMainLine(curve);

        dc.writeToFile(filename);

       // dc.SaveBitmapTemplate(res->plot, "triple.bm");
        //dc.SaveBitmapTemplate(res->maskedPlot, "tripleMs.bm");


    }
}




void SaveResultsGraph2( PipelineParameters parameters,PipelineResult2* res,  const char* filename,
                                                       const char* filename2, const char* filename3,
                                                       const char* filename4){

    int saveRes = 200;
    {
        DrawContext dc { static_cast<int>(res->plot.width()), static_cast<int>(res->plot.height()), saveRes, saveRes };
        dc.setColor(0, 0, 1);
        dc.drawPlot(res->plot);
        //dc.drawPlotBackground(normalizeImageValues(res->plot));

        dc.setColor(1,0,0);
        //dc.drawVoronoiGraph(res);

        /*for(int k = 0; k < res->numberOfComponents; k++){

           dc.drawSkeleton( res->results.at(k).skeleton);
        }*/
        //dc.drawMainLine(curve);
        //dc.drawControlPoints(curve, 2);
        dc.writeToFile(filename2);

    }


    {
        DrawContext dc { static_cast<int>(res->plot.width()), static_cast<int>(res->plot.height()), saveRes, saveRes };
        dc.setColor(1, 1, 1);
        //dc.drawPlot(res->plot);
        dc.drawPlotBackground(normalizeImageValues(res->plotGauss));
        //dc.drawPlotBackground(normalizeImageValues(res->plot));

        //dc.setColor(1,0,0);
        //dc.drawVoronoiGraph(res);

        /*for(int k = 0; k < res->numberOfComponents; k++){

           dc.drawSkeleton( res->results.at(k).skeleton);
        }*/
        //dc.drawMainLine(curve);
        //dc.drawControlPoints(curve, 2);
        dc.writeToFile(filename4);

    }


    {
        DrawContext dc { static_cast<int>(res->plot.width()), static_cast<int>(res->plot.height()),  saveRes, saveRes };
        dc.setColor(1, 1, 1);

        for(int k = 0; k < res->numberOfComponents; k++){

           dc.drawPlotBackground(normalizeImageValues(res->results.at(k)->dt));
        }
        dc.writeToFile(filename3);

    }



    {
        DrawContext dc { static_cast<int>(res->plot.width()), static_cast<int>(res->plot.height()),  saveRes, saveRes };

        //dc.drawNor (res.maskedPlot);
        //
        //auto tmp = GetBands(&(res->graphNodes), res->plot.width());
        //dc.drawPlot(res->plot);
        dc.setColor(1, 1, 0);
        //dc.drawPlotBackground(tmp);
        //dc.drawPlotBackground(normalizeImageValues(res->plot));
        //dc.drawPlot(res->maskedPlot);
       dc.setColor(0.3, 0.3, 1.0);

        for(int k = 0; k < res->numberOfComponents; k++){        
            dc.drawSkeleton(res->results.at(k)->thresholded);

        }
        //dc.setColor(0,1,0);
        //dc.drawControlPoints(curve, 2);

        dc.setColor(1,0,0);
        for(int k = 0; k < res->numberOfComponents; k++){

           dc.drawSkeleton( res->results.at(k)->skeleton);
        }
        dc.setColor(1, 1, 0);


        for(int i = 0; i < res->numberOfComponents; i++ ){

            std::vector<SkeletonNode*>* nodes = &(res->results.at(i)->graphNodes);
            dc.drawVoronoiGraph(nodes);

        }



        dc.setColor(1,1,0);
        dc.drawControlPoints(res->controlPoints);

        //dc.setColor(0,1,0);
        //dc.drawMainLine(curve);

        dc.writeToFile(filename);

        //dc.SaveBitmapTemplate(res->plot, "triple.bm");
        // dc.SaveBitmapTemplate(res->maskedPlot, "tripleMs.bm");


    }
}


void SaveResultsNormal(PipelineParameters parameters,PipelineResult* result, std::vector<AugmentedPoint> curve, const char* filename){
    {
        DrawContext dc { static_cast<int>(result->plot.width()), static_cast<int>(result->plot.height()), 600, 600 };
        dc.setColor(0, 0, 1);
        dc.drawPlotBackground(normalizeImageValues(result->plot));
        dc.setColor(1, 1, 0);
        //dc.drawDashLines(result->augmentedCurve);
        dc.drawBand(result->augmentedCurve);
        dc.drawMainLine(result->augmentedCurve);
        dc.setColor(1,0,0);
        if (!parameters.usePrincipalComponent)
              dc.drawSkeleton(result->skeleton);
        dc.setColor(0,1,0);
        dc.drawMainLine(curve);
        dc.writeToFile(filename);
        dc.SaveBitmapTemplate(result->plot, "normal.bm");

    }
}

void SaveResults(PipelineParameters parameters,PipelineResult* result, std::vector<AugmentedPoint> curve){

    if ( parameters.asGraph){

        if (parameters.useGPUSkeletonization)
            SaveResultsGraph(parameters,result, curve,"gpuGraph.png");
        else{

            if ( parameters.usePrincipalComponent)
                SaveResultsGraph(parameters, result,curve, "compGraph.png");
            else
                SaveResultsGraph(parameters,result, curve, "skelGraph.png");
        }
    }
    else {

        if (parameters.useGPUSkeletonization)
            SaveResultsNormal(parameters,result, curve,"gpu.png");
        else{

            if ( parameters.usePrincipalComponent)
                SaveResultsNormal(parameters, result,curve, "comp.png");
            else
                SaveResultsNormal(parameters,result, curve, "skel.png");
        }
    }
}


TestResults ApplyTest(PipelineParameters parameters, DataSet dataSet, Projection2D projection, std::vector<AugmentedPoint> curve, SkeletonGenerator* generator, int numberOfAttempts,
                      bool debug = false){
    vector<float> times;
    vector<float> iterations;
    vector<float> distances;

    float meanTime = 0;
    float meanIters  = 0;
    float meanDistance = 0;

    auto plot = generator->CreateScatterplot(dataSet, projection, parameters.width,parameters.height);

    {
        DrawContext dc { static_cast<int>(plot.width()), static_cast<int>(plot.height()), 600, 600 };
        dc.setColor(0, 0, 1);
        dc.drawPlotBackground(normalizeImageValues(plot));
        dc.writeToFile("drawPlot.png");

    }

    for(int i =0 ; i < numberOfAttempts; i++){
        double time;
        PipelineResult result;
        result  = processPipeline(parameters, dataSet, projection, generator, &time, debug);

        float diff = 0;



        if ( parameters.asGraph){
            // Hausdorff distance ...
            std::vector<AugmentedPoint> resCurve;
            std::vector<SkeletonNode*>* nodes = &(result.graphNodes);

            for(int i = 0; i < result.graphNodes.size(); i++){
                auto node = result.graphNodes.at(i);
                AugmentedPoint pt;
                pt.amountPerLength = node->amountPerLength;
                pt.orthoVariance = node->orthoVariance;
                pt.point = Int2{ node->location[0], node->location[1]};            
                resCurve.push_back(pt);
                // now get as well the neighbor's locations ...

                for(int j = 0; j < node->indexNeighborNode.size(); j++ ){
                    int o = node->indexNeighborNode.at(j);
                    if ( o <  nodes->size()){
                        auto other = nodes->at(o);
                        Int2 otherLoc = other->location;
                        float totalDistance =  sqrt(pow(otherLoc[0] -pt.point[0],2.0) +pow(otherLoc[1] -pt.point[1],2.0));
                        if ( totalDistance > 1.0){
                            float currentDistance = 1;

                            while( currentDistance < totalDistance){
                                // Linearly interpolate the positions
                                float percent = currentDistance/totalDistance;
                                float newX =  (1.0 -percent)*pt.point[0] + (percent)*otherLoc[0];
                                float newY =  (1.0 -percent)*pt.point[1] + (percent)*otherLoc[1];
                                AugmentedPoint midPoint;
                                midPoint.amountPerLength = node->amountPerLength*(1.0 - percent) + (percent)*other->amountPerLength;
                                midPoint.point[0] = newX;
                                midPoint.point[1] = newY;
                                midPoint.orthoVariance =  node->orthoVariance*(1.0 - percent) + (percent)*other->orthoVariance;

                                resCurve.push_back(midPoint);
                                currentDistance += 1;
                            }
                        }
                    }
                }
            }
            diff = WeightedHausdorffDistance(resCurve, curve, 0.0, 0.0, plot.width() );
        }
        else
            diff = couplingDistanceAugmented(result.augmentedCurve, curve, 0.0, 0.0);

        times.push_back(time);
        iterations.push_back(result.iterations);
        distances.push_back(diff);
        meanTime += time;   meanIters += result.iterations; meanDistance += diff;


        if (  i == 0 ){
               SaveResults(parameters, &result, curve);
        }
    }
    meanTime /= numberOfAttempts;
    meanIters /= numberOfAttempts;
    meanDistance /= numberOfAttempts;

    TestResults result;
    result.meanTime = meanTime;
    result.meanIters = meanIters;
    result.meanDistance = meanDistance;


    // Calculate Averages & Dev.
    float varTime = 0;
    float varDistance = 0;
    for(int j = 0; j < numberOfAttempts; j++){
        varTime     += pow( times.at(j) - meanTime,2.0);
        varDistance += pow( distances.at(j) -meanDistance, 2.0);
    }
    varTime /= numberOfAttempts;
    varDistance /= numberOfAttempts;
    result.devTime = sqrt(varTime);
    result.devDistance = sqrt(varDistance);
    return result;
}


void ApplyRobustnessTest(PipelineParameters parameters, SkeletonGenerator* generator){

   Projection2D projection = AxisSelection2D(0,1);
   float avgDistance = 0;
   float devDistance = 0;
   vector<float> robust;
   for(int i = 2; i <= 30; i++){

       float q =  static_cast<float>(i)/2.0;
       int sz = 256;

       std::vector<float> es{-0.5,0,2,0};

       DataSet dataSet = fakeData(5000,es,q+10.0);
       //DataSet dataSet = fakeSpiral(2000,q);
       //DataSet dataSet = fakeCircle(5000, q);
       auto curve = generatingCurve(50, es, sz, dataSet, 0,1);
       // auto curve = generatingCircle(50, sz, dataSet,0,1);
       //auto curve = generatingSpiral(50,sz,dataSet, 0,1);

       auto res = ApplyTest(parameters,dataSet, projection, curve, generator, 1, false);

       avgDistance += res.meanDistance;
       robust.push_back(res.meanDistance);
   }

   avgDistance /= 30.0;

   for(int k = 0; k < robust.size(); k++)
       devDistance += (robust.at(k) - avgDistance)*(robust.at(k) - avgDistance);


   devDistance /= robust.size();
   std::cout << "Robustness Test Distance " << avgDistance << "," <<  devDistance << std::endl;

}


int NumberOfGraphs(std::vector<SkeletonNode*>* graphNodes, std::vector<int>* elems);

double Clumpy(PipelineResult2* pgraph);
double Clumpy(std::vector<SkeletonNode*>* graphNodes);
///********************
double SpearmanCorrelation(PipelineResult2 *result);
double Pearson(PipelineResult2* result);

//double Convex(Image2D<char> thresholded);
pair<double, double> ConvexAndSkinny(Image2D<char> thresholded);



double Skewed(PipelineResult2* pgraph, bool bowley=true);
double SkewedHelper(std::vector<SkeletonNode*>* graphNodes, double maxRange);
double SkewedHelperBowley(std::vector<SkeletonNode*>* graphNodes, double maxRange);


double Straight(PipelineResult2* result);
double Straight(std::vector<SkeletonNode*>* graphNodes);


double Disjoint(PipelineResult2* result);
double Scattered(PipelineResult2* result);


double Stringy(PipelineResult2* result);
pair<double, double> StringyHelper(std::vector<SkeletonNode*>* graphNodes);


double Bendy(PipelineResult2* result);
pair<double, double> BendyHelper(std::vector<SkeletonNode*>* graphNodes);

//*****************

double Thick(std::vector<SkeletonNode*>* graphNodes, double maxRange);

//Redefine these values ...
double UniformThickness(std::vector<SkeletonNode*>* graphNodes, double maxRange);



void ReadScatterplotsRScag(std::string cpath, std::vector<std::string>* plots ){

    const char* dir_path = cpath.c_str();

    DIR *dir;
    DIR *dir2;

    struct dirent *ent;
    if ((dir = opendir (dir_path)) != nullptr) {
       /* print all the files and directories within directory */
       while ((ent = readdir (dir)) != nullptr) {
            std::string name = ent->d_name;
            if (name.compare(".") == 0) continue;
            if (name.compare("..") == 0) continue;
            std::string path = dir_path;

            std::string fullPath = path + "/" + name;
            if  ((dir2 = opendir (fullPath.c_str() )) != nullptr) continue;
            std::string suffix = name.substr( name.length()-3, name.length());

             if ( suffix.compare("csv") == 0)
                 plots->push_back(fullPath);

            //check the suffix
            /*if (plots->size() == 43108){
                std::cout << suffix << "--" << std::endl;
           }*/
            //->InsertNextValue(fullPath);

            //printf ("%s\n", ent->d_name);
        }
        closedir (dir);
    }
}


void LoadRScagDataFromScatterplot(std::string dataPath, DataSet* dataset){

    std::ifstream file;
    file.open(dataPath);
    std::string line;
    while(std::getline(file, line)){

        std::string del = ",";
        std::vector<std::string> strings;
        boost::split(strings, line, boost::is_any_of(del));
        float v1 = atof( strings.at(0).c_str());
        float v2 = atof( strings.at(1).c_str());

        std::vector<float> values;
        values.push_back(v1);
        values.push_back(v2);
        dataset->append(values);
    }
}

void ReadScatterplots(std::string path, std::vector<std::string>* plots){
    std::ifstream file;
    file.open(path);

    if (!file.is_open()){
        std::cerr << "Error in Reader::ReadDescription() - Error opening the file in the specified directory" << std::endl;
        exit(1);
    }
    std::string line;
    while(std::getline(file, line)){
        plots->push_back(line);
    }

}

void ReadSkip(std::string path, std::set<std::string>* plots){
    std::ifstream file;
    file.open(path);

    if (!file.is_open()){
        std::cerr << "Error in Reader::ReadDescription() - Error opening the file in the specified directory" << std::endl;
        exit(1);
    }
    std::string line;
    int i = 0;
    while(std::getline(file, line)){

        std::vector<std::string> info;
        boost::split(info, line, boost::is_any_of("\t"));
        if (i > 0)
            plots->insert(info.at(1));
        i++;
    }

}

pair<std::string, std::string> LoadDataFromScatterplot(std::string data, DataSet* dataset, bool rotate=false, float angle=0){
    std::string del = ",";
    std::vector<std::string> strings;
    boost::split(strings, data, boost::is_any_of(del));

    std::string firstInfo = strings.at(0);
    std::vector<std::string> info;
    boost::split(info, firstInfo, boost::is_any_of(";"));


    std::string numberOfElemString = info.at(0);
    int numItems =atoi(numberOfElemString.c_str());
    //<< "Total Elements " << numItems <<" .. " << strings.size() <<std::endl;
    for(int i = 0; i < numItems; i++){
        int idx1 = i*2 + 0;
        int idx2 = i*2 + 1;

        float v1 = atof( strings.at(idx1+1).c_str());
        float v2 = atof( strings.at(idx2+1).c_str());
        if (rotate){
            float newV1 = v1*cos(angle) - v2*sin(angle);
            float newV2 = v1*sin(angle) + v2*cos(angle);
            v1 = newV1;
            v2 = newV2;
        }
        std::vector<float> values;
        values.push_back(v1);
        values.push_back(v2);
        dataset->append(values);
    }

    return std::make_pair(info.at(1), info.at(2));
}


void RedefineRanges(DataSet* dataset){

    /*
    << "Redefine  ranges " <<  dataset->size() <<std::endl;
    //Attribute and then index

    << "Current X Range " << std::endl;
    << dataset->getStat(0).minValue <<"," <<  dataset->getStat(0).maxValue << std::endl;


   << "Current Y Range " << std::endl;
    << dataset->getStat(1).minValue <<"," <<  dataset->getStat(1).maxValue << std::endl;*/

    float xRange[2] = { dataset->getStat(0).minValue, dataset->getStat(0).maxValue  };
    float yRange[2] = { dataset->getStat(1).minValue, dataset->getStat(1).maxValue  };


   for(int i = 0; i < dataset->size(); i++){
       float v1 = (*dataset)(0,i);
       float v2 = (*dataset)(1,i);

       for(int k = 0; k < 360; k++){
           float angle = (k/180.0)*3.14159;
           float newV1 = v1*cos(angle) - v2*sin(angle);
           float newV2 = v1*sin(angle) + v2*cos(angle);

           if (newV1 < xRange[0]) xRange[0] = newV1;
           if (newV1 > xRange[1]) xRange[1] = newV1;

           if (newV2 < yRange[0]) yRange[0] = newV2;
           if (newV2 > yRange[1]) yRange[1] = newV2;

       }
   }



}

void CallToRScagnostics(int index, int numBins, int rotation){
    std::stringstream ss;
    ss << index << " " << numBins << " " << rotation;
    std::string call ="Rscript RTest.r " + ss.str();
    system(call.c_str());
}

void Merge(std::string dtName, std::string resName, std::string plot, std::string gauss, std::string finalName){



   //sepVert.png
    std::string call ="convert +append " + dtName + " sepVert.png tmp1.png";
    system(call.c_str());

   call ="convert +append tmp1.png " + resName + " tmp2.png";
   system(call.c_str());

   call ="convert +append tmp2.png sepVert.png tmp3.png";
   system(call.c_str());

   call ="convert +append tmp3.png " + plot + " tmp4.png";
   system(call.c_str());

   call ="convert +append tmp4.png sepVert.png tmp2.png";
   system(call.c_str());

   call ="convert +append tmp2.png " + gauss + " " + finalName;
   system(call.c_str());

    call = "rm tmp1.png tmp2.png tmp3.png " + gauss + " " + dtName  + " " + resName + " " + plot;
    system(call.c_str());

}

void CleanResult(PipelineResult2* result){
    //clean accumulator
    //<< "Cleaning accumulator " << std::endl;
    int width = result->accumulator.width();
    int height = result->accumulator.height();
    for(int x = 0; x < width; x++) {
        for(int y = 0; y < height; y++) {

            (result->accumulator)(x,y)->clear();
            delete (result->accumulator)(x,y);
        }
    }

    for(int i = 0; i < result->results.size(); i++){
        ComponentResult* res = (result->results.at(i));

        for(int j = 0; j < res->graphNodes.size();j++){
            res->graphNodes.at(j)->Clean();
            SkeletonNode* cNode =res->graphNodes.at(j);
            delete cNode;
        }


        res->graphNodes.clear();
        delete res;
    }
    result->results.clear();

}

double ComputeSigma2(const DataSet& dataSet, Projection2D projection, size_t width, size_t height) {
    Image2D<unsigned> img {width, height};


    for(int i = 0; i < dataSet.size(); i++) {
        auto pair = projection(dataSet, i);
        float v1 = pair[0];
        float v2 = pair[1];

        int x = v1 * width;
        x = std::max(0, std::min(x, (int)(width - 1)));
        int y = v2 * height;
        y = std::max(0, std::min(y, (int)(height - 1)));
        img(x, y)++;
    }

    vector< pair<double, double> > pts;

    double spacing = 1.0/width;

    for(int x = 0; x < width; x++) {
        for(int y = 0; y < height; y++) {

            if (img(x,y) > 0){
                pts.push_back(  std::make_pair( x*spacing+ 0.5*spacing, y*spacing+ 0.5*spacing ));
            }
        }
    }



    double sigma = 0;

    double minDs[pts.size()];
    for(int i = 0; i < pts.size(); i++) {
         minDs[i] = DBL_MAX;
    }

    for(int i = 0; i < pts.size(); i++) {


        //find the nearest neighbor... this can be done faster in a much different waz...
        auto pair1  =pts.at(i);

        for(int j = i+1; j < pts.size(); j++) {
            auto pair2 = pts.at(j);


            double d = sqrt(pow(pair1.first- pair2.first,2.0) +pow(pair1.second- pair2.second,2.0));
            minDs[i] = min(minDs[i], d);
            minDs[j] = min(minDs[j], d);
        }

    }
    for(int i = 0; i <  pts.size(); i++) {
         sigma += minDs[i];// = DBL_MAX;
    }

    sigma = sigma / (pts.size()*2.0);

    return sigma;



}

double ComputeSigmaQuadTree(const DataSet& dataSet, const Projection2D& projection);


double ComputeSigma(const DataSet& dataSet, const Projection2D& projection){
    //The radius R is set to the average distance
    // δ of a point in S to its nearest- non-zero neighbor.

    double sigma = 0;

    double minDs[dataSet.size()];
    for(int i = 0; i < dataSet.size(); i++) {
         minDs[i] = DBL_MAX;
    }

    for(int i = 0; i < dataSet.size(); i++) {


        //find the nearest neighbor... this can be done faster in a much different waz...
        auto pair1 = projection(dataSet, i);
        for(int j = i+1; j < dataSet.size(); j++) {
            auto pair2 = projection(dataSet, j);

            double d = sqrt(pow(pair1[0]- pair2[0],2.0) +pow(pair1[1]- pair2[1],2.0));

            //if ( d > 0.0000001){
                minDs[i] = min(minDs[i], d);
                minDs[j] = min(minDs[j], d);
           // }
        }

    }
    for(int i = 0; i < dataSet.size(); i++) {
         sigma += minDs[i];// = DBL_MAX;
    }

    sigma /= dataSet.size();

    if ( sigma < 0.00001 )
        sigma = 0.004;
    return sigma;
}

int main(int argc, char** argv) {
    //std::ifstream infile(argv[1]);
    //
 //   Projection2D projection = AxisSelection2D(atoi(argv[2]), atoi(argv[3]));

    SkeletonGenerator* generator = new SkeletonGenerator;


    std::string path = argv[1];

    std::string skip = argv[2];

    std::cout << "Path?" << path << std::endl;
    std::cout << "Skip?" << skip << std::endl;

    std::vector<std::string> scatterplots;
    std::set<std::string> skipping;
    //ReadScatterplots(path, &scatterplots);
    ReadScatterplotsRScag(path, &scatterplots);
    ReadSkip(skip, &skipping);

    int sz = 192;
    Projection2D projection = AxisSelection2D(0,1);

    PipelineParameters parameters;
                     parameters.maxIterations = 50;
                     parameters.width = parameters.height = sz;
                     parameters.segment = 4;
                     parameters.smoothCount = 1;
                     parameters.useGPUSkeletonization = false;
                     parameters.asGraph = true;
                     parameters.smoothCountVariance = 1;
                     parameters.sigma = 4;
                     parameters.threshold = 20.0f;
                     parameters.firstAxis = 0;
                     parameters.secondAxis = 1;


    const int amtSizes = 10;
    int sizes[amtSizes];
    for(int i = 0; i < amtSizes; i++)
        sizes[i] = 64 + i*16;


    std::stringstream fullRes;


    std::cout <<"Scatterplot\t";    
    std::cout <<"Angle\t";
    std::cout <<"Size\t";
    std::cout <<"Num Nodes\t";
    std::cout <<"Num Graphs\t";

    std::cout <<"Clumpy\t";//1
    std::cout <<"Skewed\t";//2
    std::cout <<"String\t";//3
    std::cout <<"Convex\t";//4
    std::cout <<"Straight\t";//5
    std::cout <<"Bendy\t"; //6
    std::cout <<"Skinny\t"; //7
    std::cout <<"Time\t\n"; //7



    const int numMeasures = 8;
    bool justR = false;

    std::cout << "Total Scatterplots " << scatterplots.size() << " .." <<  skipping.size() << std::endl;

    vector<int> scatterplotsToCheck;
    if ( scatterplots.empty()) return -1;


    enum RobustnessTest { BIN, ADD, REMOVE, JITTER , LIST };


    RobustnessTest currentTest = RobustnessTest::BIN;

    if ( argc > 4){
        if (  strcmp(argv[4],"ADD") == 0 )
            currentTest = RobustnessTest::ADD;
        if ( strcmp(argv[4],"BIN") == 0)
            currentTest = RobustnessTest::BIN;
        if ( strcmp(argv[4],"REMOVE") == 0)
            currentTest = RobustnessTest::REMOVE;
        if ( strcmp(argv[4],"JITTER") == 0)
            currentTest = RobustnessTest::JITTER;
        if ( strcmp(argv[4],"LIST") == 0)
            currentTest = RobustnessTest::LIST;

    }

    if ( currentTest == RobustnessTest::LIST){
        fullRes << "Scatterplot\t";
        fullRes << "Plot name\t";
        fullRes << "#Points\n";

    }
    else {
        fullRes << "Scatterplot\t";
        fullRes << "Size\t";
        fullRes << "Percentage Filled\t";
        fullRes <<"Clumpy\t";//1
        fullRes <<"Skewed\t";//2
        fullRes <<"String\t";//3
        fullRes  <<"Convex\t";//4
        fullRes  <<"Straight\t";//5
        fullRes  <<"Bendy\t"; //6
        fullRes  <<"Skinny\t"; //7
        fullRes  <<"Time\t\n"; //8
    }




    std::cout << "Current Test? " << currentTest << std::endl;

    int N =  scatterplots.size(); //scatterplotsToCheck.size();
    double maxDeviations[numMeasures];
    for(int i = 0; i < numMeasures; i++)
        maxDeviations[i] = 0;


    float percent = 4.5f/100.0f;//

    int incomplete = 0;
    for(int i = 0; i <  N; i++){
        int currentlyTesting = i;//scatterplotsToCheck.at(i);



        std::ifstream infile(scatterplots.at(currentlyTesting));

        std::vector<std::string> info;
        boost::split(info,  scatterplots.at(currentlyTesting) , boost::is_any_of("/"));

        DataSet tmpDataset = parseData(infile);
        bool jumped = false;
        if ( count( skipping.begin(), skipping.end(), info.back() ) > 0)continue;

        if ( tmpDataset.size() < 100) continue;

        std::cout <<"Processing " << currentlyTesting << "," << info.back() <<  "-" <<  tmpDataset.size() << " - " << tmpDataset.dimension() << std::endl;

        if ( currentTest == RobustnessTest::LIST){
            fullRes << i << "\t" <<  info.back() << "\t" << tmpDataset.size() <<  "\n";
            continue;
        }
        // Over 100 attributes
        int angle = 0;
        double measures[numMeasures][amtSizes];
        double deviations[numMeasures];
        double averages[numMeasures];
        double maxs[numMeasures];
        double mins[numMeasures];

        float averagePercent = 0;
        for(int k = 0; k < numMeasures; k++){
           averages[k] = 0;
           deviations[k] = 0;
           for(int l = 0; l < amtSizes; l++)
           measures[k][l] = 0;
        }


            if (!justR){

             std::vector<std::string> finalNames;

             if ( currentTest == RobustnessTest::BIN){
                 //do nothing
             }
             if ( currentTest == RobustnessTest::ADD || currentTest == RobustnessTest::REMOVE || currentTest == RobustnessTest::JITTER){
                 //
                 for(int i = 0; i < amtSizes; i++)
                     sizes[i] = 64;
             }




             tmpDataset.CheckRange(0);
             tmpDataset.CheckRange(1);
             std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

             tmpDataset.DoubleRangeInStat(0, true);
             tmpDataset.DoubleRangeInStat(1, true);
             auto sigma1 = 0.025; //ComputeSigmaQuadTree(tmpDataset, projection);
             auto end = std::chrono::steady_clock::now();
             //std::cout << "Second Sigma Computation " << t1 << ".."  << sigma1 << std::endl;

             for(int k = 0; k < amtSizes; k++){


                 if ( currentTest == RobustnessTest::ADD){//Re-read and compute
                     std::ifstream infile(scatterplots.at(currentlyTesting));
                     tmpDataset = parseData(infile);
                     tmpDataset.addRandomPoints(percent);
                     tmpDataset.CheckRange(0);
                     tmpDataset.CheckRange(1);
                     tmpDataset.DoubleRangeInStat(0, true);
                     tmpDataset.DoubleRangeInStat(1, true);
                 }

                 if ( currentTest == RobustnessTest::REMOVE){//Re-read and compute
                     std::ifstream infile(scatterplots.at(currentlyTesting));
                     tmpDataset = parseData(infile);
                     tmpDataset.removeRandomPoints(percent);
                     tmpDataset.CheckRange(0);
                     tmpDataset.CheckRange(1);
                     tmpDataset.DoubleRangeInStat(0, true);
                     tmpDataset.DoubleRangeInStat(1, true);
                 }


                 if ( currentTest == RobustnessTest::JITTER){//Re-read and compute
                     std::ifstream infile(scatterplots.at(currentlyTesting));
                     tmpDataset = parseData(infile);
                     tmpDataset.jitterRandomPoints(percent);
                     tmpDataset.CheckRange(0);
                     tmpDataset.CheckRange(1);
                     tmpDataset.DoubleRangeInStat(0, true);
                     tmpDataset.DoubleRangeInStat(1, true);
                 }
                //std::cout << "//////////////////////////////////////////////////////////////////// "<< std::endl;
                PipelineResult2 result;
                //std::cout << "CSigma: " << sigma1 << "   "  << sizes[k] <<  " .." << tmpDataset.size() << std::endl;
                parameters.sigma = sigma1; //Automatic...
                parameters.width = parameters.height = sizes[k];
                parameters.dtThresholdWeight = 0.5;
                parameters.dominationWeight = 0.05f;


                std::stringstream ss,ss2,ss3,ss4,ss5;
                ss << "plot_" << currentlyTesting << "_" << sizes[k] << "_" << angle << ".png";
                ss2 <<"res_"  << currentlyTesting << "_" << sizes[k] << "_" << angle << ".png";
                ss3 <<"dt_"  << currentlyTesting << "_" << sizes[k] << "_" << angle << ".png";
                ss4 << "final_"<< currentlyTesting << "_" << sizes[k] << "_" << angle << ".png";
                ss5 <<"gauss_"  << currentlyTesting << "_" << sizes[k] << "_" << angle << ".png";
                double time = 0;

                bool complete = processPipelineGraphTimeWrapped(parameters, tmpDataset, projection, generator, &time,debugAll,&result);

                if ( !complete ) {

                    incomplete++;

                    continue;
                }

                double percentageFilled = result.percentageFilled;
                debugAll = false;
                if (debugAll)
                {   SaveResultsGraph2(parameters,&result, ss.str().c_str(), ss2.str().c_str(),
                                                         ss3.str().c_str(), ss5.str().c_str());
                     Merge(ss3.str(), ss2.str(), ss.str(), ss5.str(), ss4.str());
                     finalNames.push_back(ss4.str());
                 }

                debugAll = false;

                averagePercent += percentageFilled;
                std::vector<int> diffs;
                //Six features in the end...
                // Clumpy, Skewed, Stringy, Convex,  Straight, Bendy, Skinny,
                double clumpy = Clumpy(&result);
                measures[0][k] = clumpy;

                debugAll = false;
                double skewed = Skewed(&result);
                debugAll = false;
                measures[1][k] = skewed;
                double stringy = Stringy(&result);
                measures[2][k] = stringy;

                auto convexAndSkinny = ConvexAndSkinny(result.thresholded);
                measures[3][k] = convexAndSkinny.first;
                measures[4][k] = Straight(&result);
                measures[5][k] = Bendy(&result);
                measures[6][k] = convexAndSkinny.second;
                measures[7][k] = time;

                fullRes << currentlyTesting << "\t" << sizes[k] <<"\t" << percentageFilled << "\t" <<  clumpy  << "\t";
                fullRes << skewed << "\t" << stringy << "\t" << measures[3][k] <<"\t" << measures[4][k] << "\t";
                fullRes << measures[5][k] << "\t" << measures[6][k] << "\t" << measures[7][k] << "\t\n";

                CleanResult(&result);
          }


             if ( finalNames.size() > 1){
                 std::stringstream ss;
                 ss <<   "URes_" <<  currentlyTesting << ".png";
                 std::string call ="convert -append " + finalNames.at(0) + " separation.png " + ss.str();
                 system(call.c_str());
                 call = "rm " + finalNames.at(0); // + " " + finalNames.at(1);
                 system(call.c_str());

                 for(int  g = 1; g < finalNames.size(); g++){
                      call ="convert -append " + ss.str() + " " + finalNames.at(g) + " "  + ss.str();
                      system(call.c_str());
                      call ="convert -append " + ss.str() + " separation.png "  + ss.str();
                      system(call.c_str());
                      call = "rm " + finalNames.at(g);// + " " + finalNames.at(1);
                      system(call.c_str());
                 }
             }

            for(int h = 0; h < numMeasures; h++){
                maxs[h] = measures[h][0];
                mins[h] = measures[h][0];
                for(int k = 0; k < amtSizes; k++){
                     averages[h] += measures[h][k];

                     maxs[h] = max(maxs[h],measures[h][k]);
                     mins[h] = min(mins[h],measures[h][k]);

                }

                averages[h] /= amtSizes;
            }


            for(int h = 0; h < numMeasures; h++){
                for(int k = 0; k < amtSizes; k++){
                     deviations[h] += pow(averages[h]-  measures[h][k],2.0);
                }

                deviations[h] /= amtSizes;
                deviations[h] = sqrt(deviations[h]);
                if ( deviations[h] > maxDeviations[h]) maxDeviations[h] = deviations[h];
            }

            if (!jumped){

                averagePercent /= amtSizes;

                fullRes << " \t \t \t";
                std::cout << currentlyTesting << "\t";
                std::cout << tmpDataset.size() << "\t";
                std::cout<< "\t\t\t";
                for(int h = 0; h < numMeasures; h++){
                    std::cout << averages[h];
                    fullRes << averages[h] << "\t";

                    if( h != numMeasures-1) std::cout << "\t";
                }

                fullRes << "\nDev\t" << currentlyTesting <<  "\t" << averagePercent << "\t";

                std::cout << std::endl;
                std::cout<< "Dev\t\t\t\t\t";
                for(int h = 0; h < numMeasures; h++){
                    std::cout << deviations[h];
                    fullRes << deviations[h] << "\t";
                    if( h != numMeasures-1) std::cout << "\t";
                }
                std::cout << std::endl;
                fullRes << "\n";

                std::cout<< "(Mx-Mn)/Mx\t\t\t\t\t";
                for(int h = 0; h < numMeasures; h++){

                    double relative = (maxs[h] - mins[h])/(maxs[h]);
                    if (isnan(relative)) relative = 0;
                    std::cout << relative;
                    if( h != numMeasures-1) std::cout << "\t";
                }
                std::cout << std::endl;

                std::cout<< "Min\t\t\t\t\t";
                for(int h = 0; h < numMeasures; h++){

                    double relative = mins[h];
                    if (isnan(relative)) relative = 0;
                    std::cout << relative;
                    if( h != numMeasures-1) std::cout << "\t";
                }
                std::cout << std::endl;
                std::cout << std::endl;
            }         
        }
    }

    std::cout << "Incomplete? " << incomplete/amtSizes << std::endl;
    std::cout<< "\t\t\t\t\t";
    for(int h = 0; h < numMeasures; h++){
        std::cout << maxDeviations[h];
        if( h != numMeasures-1) std::cout << "\t";
    }
    std::cout << std::endl;

    std::ofstream out(argv[3]);
    out << fullRes.str();
    out.close();


    delete generator;
    return 0;
}

/*
           Principal Graphs Scagnostics Measures ....


The original approach attempted to look into five aspects of the points...


Changed to Spare ...>

Outliers
1 -> outlying
Shape
2 -> Convex  *
3 -> Skinny  *  (?) Redefine
4 -> Stringy *
5 -> Straight * 1.0 - SOAM
10-> Uniform Thickness * (Redefine)
Trend
6 ->Monotonic
Density
7 -> Skewed   *
8 -> Clumpy   *
coherence
9 ->Striated

*/

double Clumpy(PipelineResult2* pgraph){

    // Create a vector of graphNodes, joining everything...

    std::vector<SkeletonNode*>* totalGraphNodes = new std::vector<SkeletonNode*>();

    for(int i = 0; i < pgraph->results.size(); i++){

        std::vector<SkeletonNode*>* localGraphNodes = &pgraph->results.at(i)->graphNodes;

        for(int k = 0; k < localGraphNodes->size(); k++){
            totalGraphNodes->push_back( localGraphNodes->at(k));
        }
    }

    double val = Clumpy(totalGraphNodes);
    return val;
}

double Clumpy(std::vector<SkeletonNode*>* graphNodes){
    //We defined clumpy as the agglomeration of points in a single area
    //if it is uniformly distributed or has a peak of points in a single point
    int totalNodes =graphNodes->size();
    if (totalNodes == 0)
        return 0;
    if (totalNodes == 1)
        return 1;

    int totalPoints = 0;
    double probAtLoc[totalNodes];
    double probMax[totalNodes];
    for(int i = 0; i < totalNodes; i++){
        probAtLoc[i] = graphNodes->at(i)->stat.pointStat.count;
        totalPoints += graphNodes->at(i)->stat.pointStat.count;
        probMax[i] = 0;
    }

    double divergence = 0;
    double max_divergence = 0;
    probMax[0] = 1.0;

    double uniform[totalNodes];
    for(int i = 0; i < totalNodes; i++){
        probAtLoc[i] /= static_cast<double>(totalPoints);
        uniform[i] = (static_cast<double>(totalPoints)/totalNodes)/static_cast<double>(totalPoints);
         if ( probAtLoc[i] > 0.000001 &&   uniform[i] > 0.00000001){
             divergence += probAtLoc[i]*log2( uniform[i]/probAtLoc[i]);
         }
         if ( probMax[i] > 0.00001)
           max_divergence += probMax[i]*log2(probMax[i]/uniform[i]);
    }

    divergence *= -1;

    double clumpy = divergence/max_divergence;
    return clumpy;
}

double Disjoint(PipelineResult2* result){
    double graphNodes = result->totalNodes;
    std::vector<int> diffs;
    double numberOfGraphs = result->results.size();
    double disjoint = (numberOfGraphs -1)/(graphNodes -1);
    return disjoint;
}

double CalculateTprojectionPointSegment(double *pos1, double *pos2, double *x0)
{
    double t = -((pos2[0] - pos1[0])*(pos1[0] - x0[0]) +
        (pos2[1] - pos1[1])*(pos1[1] - x0[1]) +
        (pos2[2] - pos1[2])*(pos1[2] - x0[2])) /
       ((pos2[0] - pos1[0])*(pos2[0] - pos1[0]) +
        (pos2[1] - pos1[1])*(pos2[1] - pos1[1]) +
        (pos2[2] - pos1[2])*(pos2[2] - pos1[2]));

   return t;
}

void GetProjectedPointSegment(double *pos1, double *pos2, double *x0, double projectedPoint[]){

    double t = CalculateTprojectionPointSegment(pos1, pos2, x0);
    double dv[3] = { pos2[0] - pos1[0], pos2[1] -pos1[1], pos2[2] -pos1[2]};

    for(int i = 0; i < 3; i++) projectedPoint[i] = pos1[i] + t*dv[i];

}

double Skewed(PipelineResult2* pgraph, bool bowley){
    std::vector<SkeletonNode*>* totalGraphNodes = new std::vector<SkeletonNode*>();

    for(int i = 0; i < pgraph->results.size(); i++){

        std::vector<SkeletonNode*>* localGraphNodes = &pgraph->results.at(i)->graphNodes;

        for(int k = 0; k < localGraphNodes->size(); k++){
            totalGraphNodes->push_back( localGraphNodes->at(k));
        }
    }

    double val = 0;
    if ( bowley)
        val = SkewedHelperBowley( totalGraphNodes, 10.0);
    else
        val = SkewedHelper(totalGraphNodes, 10.0);

    return val;
}

double Stringy(PipelineResult2* pgraph){

    double geoDistance = 0;
    double maxDist = 0;
    for(int i = 0; i < pgraph->results.size(); i++){
        std::vector<SkeletonNode*>* localGraphNodes = &pgraph->results.at(i)->graphNodes;
        auto local = StringyHelper(localGraphNodes);
        maxDist += local.first;
        geoDistance += local.second;

    }

    //std::cout << "Max Dist " << maxDist << " Geo Dist" << geoDistance << std::endl;
    if ( fabs(maxDist - geoDistance) < 0.0001)
        return 1.0;
    return maxDist/geoDistance;

}

double Scattered(PipelineResult2* pgraph){
    //If it's sparse there won't be a very full mesh.. so a lot of bifurcations, and endpoints

    int tot = 0;
    int special = 0; //inner nodes are connecting two points otherwise it is multiply
    for(int i = 0; i < pgraph->results.size(); i++){

        std::vector<SkeletonNode*>* localGraphNodes = &pgraph->results.at(i)->graphNodes;

        for(int k = 0; k < localGraphNodes->size(); k++){
            SkeletonNode* node = localGraphNodes->at(k);
              tot++;
              if ( node->indexNeighborNode.size() != 2 )
                  special += 1;
        }


   }

    if ( pgraph->results.size() == 1){
        //if there's only one plot... we don't count two end points....
        special -= 2;
        tot -= 2;

    }


    double sparse = static_cast<double>(special)/static_cast<double>(tot);
   return sparse;
}


double SkewedHelper(std::vector<SkeletonNode*>* graphNodes, double maxRange){
    // We define skewness in relationship with the how the points lie at each direction of the sample points...
    // sample skewness for each node ...
    int totalNodes = graphNodes->size();

    int totalUseNodes = 0;
    double skewed = 0;
    int sumOfCounts = 0;
    int sumOfValues = 0;
    for(int i = 0; i < totalNodes; i++){
        SkeletonNode* node = graphNodes->at(i);
        double pos1[3] = {node->stat2.pointStat.average.data.at(0), node->stat2.pointStat.average.data.at(1), 0};
        double dir[2] ={node->direction.data.at(0), node->direction.data.at(1)};


        bool zero = ( dir[0] < 0.001 && dir[1] < 0.0001);
        bool nand = isnan(node->direction.data.at(0)) ||  isnan(node->direction.data.at(1));
        if ( zero  ||nand){ // dir is 0,0,

            // Find the closest node + use the direction to the closest node
            double minDist = 999999;
            for(int k = 0; k < totalNodes; k++){
                if ( i==k) continue;
                SkeletonNode* tnode = graphNodes->at(k);

                double tmpPos[3] = {tnode->stat2.pointStat.average.data.at(0), tnode->stat2.pointStat.average.data.at(1), 0};

                double d = sqrt(pow(tmpPos[0] - pos1[0],2.0) + pow(tmpPos[1] - pos1[1],2.0) );

                if ( d < minDist){
                    minDist = d;
                    dir[1] = tmpPos[1] - pos1[1];
                    dir[0] = tmpPos[0] - pos1[0];

                }

            }

        }

        double ninety = 3.14159*0.5;

        //90° in one direction
        double oneDir[2], otherDir[2];
        oneDir[0] = dir[0]*cos(ninety) - dir[1]*sin(ninety);
        oneDir[1] = dir[0]*sin(ninety) + dir[1]*cos(ninety);
        double pos2[3] = { pos1[0] + oneDir[0]*maxRange, pos1[1] + oneDir[1]*maxRange, 0 };


        int tot1 = 0, tot2 = 0;
        vector<double> values;
        double avg = 0;
        for(int k = 0; k < node->stat2.pointStat.values.size(); k++){
            double p[3] = { node->stat2.pointStat.values.at(k)[0], node->stat2.pointStat.values.at(k)[1], 0 };
            double t = CalculateTprojectionPointSegment(pos1, pos2, p);
            if ( t >= 0){
                double projectedPoint[3];
                GetProjectedPointSegment(pos1, pos2, p, projectedPoint);
                double d = sqrt(pow(pos1[0] - projectedPoint[0],2.0) + pow(pos1[1] - projectedPoint[1],2.0));
                values.push_back(d);

                //if (fabs(d) < 0.0001)  std::cout << "D is 0 " <<  projectedPoint[0] << " / " <<  projectedPoint[1] << std::endl;

                avg += d;
                tot1 ++;
            }
            if (isnan(t) && debugAll){
                std::cout << "Nan t" << std::endl;

                std::cout << "Not A Number T Point " << p[0]<<"," <<p[1]<<"."<<   std::endl;
                std::cout << "Dir? "<< dir[0] <<"," << dir[1] <<std::endl;
                std::cout << "Position " << pos1[0 ]<< "," << pos1[1] <<std::endl;
                std::cout << "One Dir " << oneDir[0] << "," << oneDir[1] <<std::endl;
                std::cout << "Other End " << pos2[0 ]<< "," << pos2[1] <<std::endl;
                //std::cout << "Total neighbors? " << node->neighborsNodes.size() <<std::endl;
                std::cout << "..............." << std::endl;

            }
        }

        //90° in the other direction
        ninety = -ninety;
        otherDir[0] = dir[0]*cos(ninety) - dir[1]*sin(ninety);
        otherDir[1] = dir[0]*sin(ninety) + dir[1]*cos(ninety);
        double pos3[3] = { pos1[0] + otherDir[0]*maxRange, pos1[1] + otherDir[1]*maxRange, 0 };
        for(int k = 0; k < node->stat2.pointStat.values.size(); k++){
            double p[3] = { node->stat2.pointStat.values.at(k)[0], node->stat2.pointStat.values.at(k)[1], 0 };
            double t = CalculateTprojectionPointSegment(pos1, pos3, p);
            if ( t >= 0){
                double projectedPoint[3];
                GetProjectedPointSegment(pos1, pos3, p, projectedPoint);

                double d =  sqrt(pow(pos1[0] - projectedPoint[0],2.0) + pow(pos1[1] - projectedPoint[1],2.0));

                //if (fabs(d) < 0.0001)  std::cout << "D is 0 " <<  projectedPoint[0] << " / " <<  projectedPoint[1] << std::endl;
                values.push_back(-d);
                avg += (-d);
                tot2++;
            }
            if (isnan(t) && debugAll){
                 std::cout << "Nan t" << std::endl;
            }
        }

        //The samples are the distance to the projected location
        avg /= (values.size()-1);
        double num = 0;
        double den = 0;
        for(unsigned int i = 0; i < values.size(); i++){
            num += pow( values.at(i) - avg,3.0);
            den += pow( values.at(i) - avg,2.0);
        }
        double pnum = num;
        num /= values.size();
        double pden = den;
        den = pow(den/values.size(), 1.5);

        double coef = sqrt( values.size()*(values.size()-1)) /( values.size() - 2);

        double sampleSkewness = coef* (num/den);

        if ( values.empty() &&  node->stat2.pointStat.count != 0 ){
            std::cout << "Empty when should be  "<< node->stat2.pointStat.count << "\n";
            std::cout << "Dir? "<< dir[0] <<"," << dir[1] <<std::endl;
            //std::cout << "Total neighbors? " << node->neighborsNodes.size() <<std::endl;
        }

        sumOfCounts += node->stat2.pointStat.count;
        sumOfValues +=  node->stat2.pointStat.values.size() ;
        if ( node->stat2.pointStat.values.size() != node->stat2.pointStat.count){
            std::cout <<"Not the same size? " << values.size() << ";" << node->stat2.pointStat.values.size() << "," <<  node->stat2.pointStat.count << std::endl;
        }


        if (isinf(sampleSkewness) && values.size() > 2 ){
            std::cout << "Infinte samples skewness? " << coef <<".." << num << ".. " << den <<std::endl;
            std::cout << values.size() <<" .." << node->stat2.pointStat.count  << "/" << node->stat2.pointStat.values.size() << std::endl;
            std::cout << "Position " << pos1[0 ]<< "," << pos1[1] <<std::endl;
            std::cout << "Dir? "<< dir[0] <<"," << dir[1] <<std::endl;
        }


        if (values.size() > 2){

            if (isnan(sampleSkewness)){

            }
            else {
               skewed += fabs(sampleSkewness);
               totalUseNodes ++;
            }
           // std::cout << "Size "<< values.size() <<"/ "<< fabs(sampleSkewness) <<".."<< skewed <<"/"<< totalNodes << ":"<< skewed/totalNodes <<" , " << tot1 <<"," << tot2 << std::endl;
        }
    }

    //std::cout << "Counts / Values " <<  sumOfCounts <<", "<< sumOfValues <<std::endl;

    //::cout << "Total usable nodes? " << totalUseNodes << std::endl;
    //std::cout << "total nodes " << totalNodes << " .. " <<  skewed << std::endl;
    if (totalUseNodes > 0)
         skewed /= totalUseNodes;
    return skewed;
}


double SkewedHelperBowley(std::vector<SkeletonNode*>* graphNodes, double maxRange){
    // We define skewness in relationship with the how the points lie at each direction of the sample points...
    // sample skewness for each node ...
    int totalNodes = graphNodes->size();

       if ( debugAll) std::cout << "trying to go with bowley " << std::endl;
    int totalUseNodes = 0;
    double skewed = 0;

    int totalPoints = 0;
    int actTotalPoints = 0;

    double weights[totalNodes];
    double skewM[totalNodes];

    for(int i = 0; i < totalNodes; i++){
        SkeletonNode* node = graphNodes->at(i);

        if ( debugAll){
            std::cout << "////////////////////////////////////////////////////"  << std::endl;
           std::cout << "Node " << i << " PTS> " <<  node->stat2.pointStat.values.size() << std::endl;

        }
        weights[i] = 0;
        skewM[i] = 0;

        double pos1[3] = {node->stat2.pointStat.average.data.at(0), node->stat2.pointStat.average.data.at(1), 0};
        double dir[2] ={node->direction.data.at(0), node->direction.data.at(1)};
        //int o = currentNode->indexNeighborNode.at(j);
        actTotalPoints += node->stat2.pointStat.values.size();

        bool zero = ( dir[0] < 0.001 && dir[1] < 0.0001);
        bool nand = isnan(node->direction.data.at(0)) ||  isnan(node->direction.data.at(1));
        if ( zero ||nand){ // dir is 0,0,

            // Find the closest node + use the direction to the closest node
            double minDist = 999999;
            for(int k = 0; k < totalNodes; k++){
                if ( i==k) continue;
                SkeletonNode* tnode = graphNodes->at(k);
                double tmpPos[3] = {tnode->stat2.pointStat.average.data.at(0), tnode->stat2.pointStat.average.data.at(1), 0};
                double d = sqrt(pow(tmpPos[0] - pos1[0],2.0) + pow(tmpPos[1] - pos1[1],2.0) );
                if ( d < minDist){
                    minDist = d;
                    dir[1] = tmpPos[1] - pos1[1];
                    dir[0] = tmpPos[0] - pos1[0];

                }
            }
        }

        double ninety = 3.14159*0.5;

        //90° in one direction
        double oneDir[2], otherDir[2];
        oneDir[0] = dir[0]*cos(ninety) - dir[1]*sin(ninety);
        oneDir[1] = dir[0]*sin(ninety) + dir[1]*cos(ninety);
        double n1 = sqrt(oneDir[0]*oneDir[0] + oneDir[1]*oneDir[1]);
        oneDir[0] /= n1;
        oneDir[1] /= n1;

        double pos2[3] = { pos1[0] + oneDir[0]*maxRange, pos1[1] + oneDir[1]*maxRange, 0 };

        //90° in the other direction
        ninety = -ninety;
        otherDir[0] = dir[0]*cos(ninety) - dir[1]*sin(ninety);
        otherDir[1] = dir[0]*sin(ninety) + dir[1]*cos(ninety);
        double n2 = sqrt(otherDir[0]*otherDir[0] + otherDir[1]*otherDir[1]);
        otherDir[0] /= n2;
        otherDir[1] /= n2;

        double pos3[3] = { pos1[0] + otherDir[0]*maxRange, pos1[1] + otherDir[1]*maxRange, 0 };

        if (debugAll){

            std::cout << "Original Dir " << dir[0]<< ", " << dir[1] << std::endl;
            std::cout << "Directions 1.. (" << oneDir[0] << "," << oneDir[1] << ")" << std::endl;
            std::cout << "Directions 2.. (" << otherDir[0] << "," << otherDir[1] << ")" << std::endl;
            std::cout << "One End " << pos2[0 ]<< "," << pos2[1] <<std::endl;
            std::cout << "Other End " << pos3[0 ]<< "," << pos3[1] <<std::endl;


        }


        int tot1 = 0; //, tot2 = 0;
        vector<double> values;
        double avg = 0;
        for(int k = 0; k < node->stat2.pointStat.values.size(); k++){
            double p[3] = { node->stat2.pointStat.values.at(k)[0], node->stat2.pointStat.values.at(k)[1], 0 };
            double t = CalculateTprojectionPointSegment(pos2, pos3, p);

            if ( debugAll){

                std::cout << "Point? " << k << " :: (" << p[0] << ","  << p[1] << ") T :" << t <<   std::endl;
            }

            if ( t >= 0.0){
                double projectedPoint[3];
                GetProjectedPointSegment(pos2, pos3, p, projectedPoint);
                double d = sqrt(pow(pos2[0] - projectedPoint[0],2.0) + pow(pos2[1] - projectedPoint[1],2.0));
                values.push_back(d);


                //if (fabs(d) < 0.0001)  std::cout << "D is 0 " <<  projectedPoint[0] << " / " <<  projectedPoint[1] << std::endl;

                avg += d;
                tot1 ++;
            }

            if ( debugAll){

                values.size();
            }
            if (isnan(t) && debugAll ){
                std::cout << "Nan t" << std::endl;

                std::cout << "Not A Number T Point " << p[0]<<"," <<p[1]<<"."<<   std::endl;
                std::cout << "Dir? "<< dir[0] <<"," << dir[1] <<std::endl;
                std::cout << "Position " << pos1[0 ]<< "," << pos1[1] <<std::endl;
                std::cout << "One Dir " << oneDir[0] << "," << oneDir[1] <<std::endl;
                std::cout << "Other End " << pos2[0 ]<< "," << pos2[1] <<std::endl;
                //std::cout << "Total neighbors? " << node->neighborsNodes.size() <<std::endl;
                std::cout << "..............." << std::endl;

            }
        }


        //The samples are the distance to the projected location
        sort(values.begin(), values.end());

        int Q1Idx = values.size()*0.25;
        int Q2Idx = values.size()*0.5;
        int Q3Idx = values.size()*0.75;


        if ( values.size() > 0 && Q3Idx != Q1Idx && fabs(values.at(Q3Idx) -values.at(Q1Idx)) > 0.000001  ){
            double Bowley = (values.at(Q3Idx) + values.at(Q1Idx) - 2*values.at(Q2Idx))/( values.at(Q3Idx) - values.at(Q1Idx));
            if (isnan(Bowley) && debugAll){
               std::cout << "Bowley? " << Bowley << " .. " << skewed << " ... Q3 " << values.at(Q3Idx) << " ... Q1 " <<  values.at(Q1Idx)<<    std::endl;
               std::cout << Q3Idx << " , " << Q1Idx << std::endl;

            }

            weights[i] = values.size();

            skewM[i] = Bowley;
            totalPoints += values.size();
            //skewed += Bowley;
            totalUseNodes++;
        }//if there are values but the range is 0, then skewness is 0
        else if (values.size() > 0 &&  fabs(values.back() - values.at(0)) < 0.0001){
            weights[i] = values.size();

            skewM[i] = 0;
            totalPoints += values.size();
            //skewed += Bowley;
            totalUseNodes++;
        }

        if ( debugAll){
            std::cout << "////////////////////////////////////////////////////"  << std::endl;
        }

    }


    for(int i = 0; i < totalNodes; i++){
        if (totalPoints > 0){

            // std::cout << weights[i]<< ", " << skewM[i] << std::endl;
            skewed += (weights[i]/totalPoints)*skewM[i];

        }

    }

    //std::cout << "Counts / Values " <<  sumOfCounts <<", "<< sumOfValues <<std::endl;


       if ( debugAll) std::cout << totalPoints << "/" << actTotalPoints << std::endl;
       if ( debugAll) std::cout << "Total usable nodes? " << totalUseNodes << std::endl;
    //std::cout << "total nodes " << totalNodes << " .. " <<  skewed << std::endl;
    if (totalUseNodes > 0)
         skewed /= totalUseNodes;
    return fabs(skewed);

}

#include <boost/bind.hpp>

double SpearmanCorrelation(PipelineResult2* result){
    // first we need to define the ranks of x and y
    std::vector<SkeletonNode*>* graphNodes = new std::vector<SkeletonNode*>();

    for(int i = 0; i < result->results.size(); i++){
        std::vector<SkeletonNode*>* localGraphNodes = &result->results.at(i)->graphNodes;

        for(int k = 0; k < localGraphNodes->size(); k++){
            graphNodes->push_back( localGraphNodes->at(k));
        }
    }



    //auto graphNodes = &(result->graphNodes);
    std::vector<pair<double,int>> rankX;
    std::vector<pair<double,int>> rankY;

    for(unsigned int i = 0; i < graphNodes->size(); i++){
       auto node = graphNodes->at(i);
       double x = node->stat2.pointStat.average.data.at(0);
       double y = node->stat2.pointStat.average.data.at(1);

       rankX.push_back(make_pair(x,i));
       rankY.push_back(make_pair(y,i));

    }

    std::sort( rankX.begin(), rankX.end(),
           boost::bind(&std::pair<double, int>::first,_1) >
           boost::bind(&std::pair<double, int>::first, _2 ));

    std::sort( rankY.begin(), rankY.end(),
           boost::bind(&std::pair<double, int>::first,_1) >
           boost::bind(&std::pair<double, int>::first, _2 ));

    double Xs[graphNodes->size()];
    double Ys[graphNodes->size()];

    for(int i = 0; i < graphNodes->size();i++){
        Xs[rankX.at(i).second] = i+1;
        Ys[rankY.at(i).second] = i+1;
    }

    double meanX = 0, meanY = 0;

    float total = 0;
    for(int j = 0; j < graphNodes->size(); j++){
        double xi = Xs[j];
        double yi = Ys[j];
        meanX += xi;
        meanY += yi;
        total++;

    }

    meanX /= total;
    meanY /= total;

    double up = 0;
    double den1 = 0;
    double den2 = 0;

    for(int j = 0; j <  graphNodes->size(); j++){
        double xi = Xs[j];
        double yi = Ys[j];

         up += (xi - meanX)*(yi - meanY);
         den1 += (xi - meanX)*(xi- meanX);
         den2 += (yi - meanY)*(yi- meanY);
    }

    double den = sqrt(den1)*sqrt(den2);

    double r = up/ den;

    if ( den < 0.0001)
        r = 0;

    return r;
}

double  UniformThickness(std::vector<SkeletonNode*>* graphNodes, double maxRange){
    // We define skewness in relationship with the how the points lie at each direction of the sample points...
    // sample skewness for each node ...
    int totalNodes = graphNodes->size();

    double maxThickness = 0;

    double thickness[totalNodes];
    double sumOfThickness = 0;
    for(int i = 0; i < totalNodes; i++){
        thickness[i] = 0;

        SkeletonNode* node = graphNodes->at(i);
        double pos1[3] = {node->stat.pointStat.average.data.at(0), node->stat.pointStat.average.data.at(1), 0};
        double dir[2] ={node->direction.data.at(0), node->direction.data.at(1)};
        double ninety = 3.14159*0.5;

        //90° in one direction
        double oneDir[2], otherDir[2];
        oneDir[0] = dir[0]*cos(ninety) - dir[1]*sin(ninety);
        oneDir[1] = dir[0]*sin(ninety) + dir[1]*cos(ninety);
        double pos2[3] = { pos1[0] + oneDir[0]*maxRange, pos1[1] + oneDir[1]*maxRange, 0 };

        double max1 = 0, max2 = 0;
        double avg = 0;
        for(int k = 0; k < node->stat.pointStat.values.size(); k++){
            double p[3] = { node->stat.pointStat.values.at(k)[0], node->stat.pointStat.values.at(k)[1], 0 };
            double t = CalculateTprojectionPointSegment(pos1, pos2, p);
            if ( t > 0){
                double projectedPoint[3];
                GetProjectedPointSegment(pos1, pos2, p, projectedPoint);
                double d = sqrt(pow(pos1[0] - projectedPoint[0],2.0) + pow(pos1[1] - projectedPoint[1],2.0));
                if ( d> max1)
                    max1 =d;
            }
        }

        //90° in the other direction
        ninety = -ninety;
        otherDir[0] = dir[0]*cos(ninety) - dir[1]*sin(ninety);
        otherDir[1] = dir[0]*sin(ninety) + dir[1]*cos(ninety);
        double pos3[3] = { pos1[0] + otherDir[0]*maxRange, pos1[1] + otherDir[1]*maxRange, 0 };
        for(int k = 0; k < node->stat.pointStat.values.size(); k++){
            double p[3] = { node->stat.pointStat.values.at(k)[0], node->stat.pointStat.values.at(k)[1], 0 };
            double t = CalculateTprojectionPointSegment(pos1, pos3, p);
            if ( t > 0){
                double projectedPoint[3];
                GetProjectedPointSegment(pos1, pos3, p, projectedPoint);
                double d = sqrt(pow(pos1[0] - projectedPoint[0],2.0) + pow(pos1[1] - projectedPoint[1],2.0));

                if ( d > max2)
                    max2 = d;
            }
        }

        thickness[i] = max1 + max2;

        sumOfThickness += thickness[i];
        if (max1 + max2 > maxThickness)
            maxThickness = max1 + max2;
    }

    //now we can divide
    double avgThickness = sumOfThickness/totalNodes;

    double variance = 0;
    for(int i =0 ; i < totalNodes; i++){

           variance += pow(thickness[i]/maxThickness - avgThickness/maxThickness,2.0);
    }
    variance /= (totalNodes-1);


    return 1.0 - variance;
}

double DFS(int i,std::vector<SkeletonNode*>* graphNodes){
    //
    std::stack<pair<int,double >> currentStack;
    std::set<int> seen;
    currentStack.push(make_pair(i, 0));
    double maxDistance = 0;

    while(!currentStack.empty()){
        auto top =  currentStack.top();
        double currentDistance = top.second;
        int currentId = top.first;
        currentStack.pop();

        if ( currentDistance > maxDistance)
            maxDistance = currentDistance;

        SkeletonNode* currentNode = graphNodes->at(currentId);
        for(int k = 0; k < currentNode->indexNeighborNode.size(); k++){
            int oIdx = currentNode->indexNeighborNode.at(k);
            if (  count(seen.begin(),seen.end(),oIdx) == 0){
                SkeletonNode* node2 = graphNodes->at(oIdx);
                double p1[2] = {currentNode->stat2.pointStat.average.data.at(0), currentNode->stat2.pointStat.average.data.at(1)};
                double p2[2] = {node2->stat2.pointStat.average.data.at(0), node2->stat2.pointStat.average.data.at(1)};
                double d =  sqrt( pow(p1[0]-p2[0],2.0) + pow(p1[1]-p2[1],2.0));
                currentStack.push(std::make_pair(oIdx, d+ currentDistance));
                seen.insert(currentId);
            }
        }
    }

    return maxDistance;
}




double DFSWithHistory(int i,std::vector<SkeletonNode*>* graphNodes,std::set<int>* seen){
    //
    std::stack<pair<int,double >> currentStack;
    currentStack.push(make_pair(i, 0));
    double maxDistance = 0;

    while(!currentStack.empty()){
        auto top =  currentStack.top();
        double currentDistance = top.second;
        int currentId = top.first;
        currentStack.pop();
        seen->insert(currentId);

        if ( currentDistance > maxDistance)
            maxDistance = currentDistance;

        SkeletonNode* currentNode = graphNodes->at(currentId);
        for(int k = 0; k < currentNode->indexNeighborNode.size(); k++){
            int oIdx = currentNode->indexNeighborNode.at(k);
            if (  count(seen->begin(),seen->end(),oIdx) == 0){
                SkeletonNode* node2 = graphNodes->at(oIdx);
                double p1[2] = {currentNode->stat2.pointStat.average.data.at(0), currentNode->stat2.pointStat.average.data.at(1)};
                double p2[2] = {node2->stat2.pointStat.average.data.at(0), node2->stat2.pointStat.average.data.at(1)};
                double d =  sqrt( pow(p1[0]-p2[0],2.0) + pow(p1[1]-p2[1],2.0));
                currentStack.push(std::make_pair(oIdx, d+ currentDistance));
                seen->insert(currentId);
            }
        }


    }

    return maxDistance;
}


int NumberOfGraphs(std::vector<SkeletonNode*>* graphNodes, std::vector<int>* elems){
    std::set<int> seen;

    int iters = 0;
    while( seen.size() != graphNodes->size()){

        //Search for all the graph nodes, which one
        for(int k = 0; k < graphNodes->size(); k++){
            if (count(seen.begin(), seen.end(), k) == 0){
                //haven't seen the node yet
                DFSWithHistory(k,graphNodes, &seen);
                elems->push_back(k);
            }
        }

        iters++;
    }

    return iters;
}


pair<double, double> StringyHelper(std::vector<SkeletonNode*>* graphNodes){
    //Largest geodesic path vs Sum of Geodesic distance
   // A Stringy shape is a skinny shape with no branches

    // First we can get all the edges...
    std::vector<int> leafs;
    std::set< pair<int,int> > edges;
    for(int i = 0; i < graphNodes->size(); i++){
        SkeletonNode* currentNode = graphNodes->at(i);

        if ( currentNode->indexNeighborNode.size() == 1)
            leafs.push_back(i);
        for(int k = 0; k < currentNode->indexNeighborNode.size(); k++){
            int oIdx = currentNode->indexNeighborNode.at(k);
            edges.insert( std::make_pair( min(i,oIdx),max(i, oIdx)));
        }
    }

    double sumDistance = 0;
    for(auto edge:edges){
        SkeletonNode* node1 = graphNodes->at(edge.first);
        SkeletonNode* node2 = graphNodes->at(edge.second);
        double p1[2] = {node1->stat2.pointStat.average.data.at(0), node1->stat2.pointStat.average.data.at(1)};
        double p2[2] = {node2->stat2.pointStat.average.data.at(0), node2->stat2.pointStat.average.data.at(1)};
        sumDistance += sqrt( pow(p1[0]-p2[0],2.0) + pow(p1[1]-p2[1],2.0));
    }

    //For each leaf we perform a DFS... and get the path that has the major distance
    double maxDist = 0;
    for(int i =0 ; i < leafs.size(); i++){
        double d = DFS( leafs.at(i), graphNodes);
        if ( d > maxDist)
            maxDist = d;
    }


    return std::make_pair(maxDist, sumDistance);
}


double Bendy(PipelineResult2* pgraph){
    double inflectionPoints = 0;
    double totalPossible = 0;

   for(int i = 0; i < pgraph->results.size(); i++){
         std::vector<SkeletonNode*>* localGraphNodes = &pgraph->results.at(i)->graphNodes;
        auto local = BendyHelper(localGraphNodes);

        inflectionPoints += local.first;

        if ( local.second == 0)
            totalPossible += 1;
        totalPossible += local.second;
    }

   //std::cout << "Total Bendy Possible " << totalPossible << std::endl;
   if ( totalPossible == 0)
       return 0;
    return inflectionPoints/totalPossible;
}

pair<double,double> BendyHelper(std::vector<SkeletonNode*>* graphNodes){

    //All those locations that have two edges....
    int totalPossible = 0;
    int inflections = 0;
    for(int i = 0; i < graphNodes->size(); i++){
          SkeletonNode* currentNode = graphNodes->at(i);

          if ( currentNode->indexNeighborNode.size()== 2){


              //Let's check the inflections..
              //First let's get the points...
              int idx1 = currentNode->indexNeighborNode.at(0);
              int idx2 = currentNode->indexNeighborNode.at(1);

              SkeletonNode* n1 = graphNodes->at(idx1);
              SkeletonNode* n2 = graphNodes->at(idx2);

              //Create the 5 points
              double pts[7][2];
              double derivative[7][2];
              double secondDerivative[7][2];


              pts[0][0] = n1->stat2.pointStat.average.data.at(0);
              pts[0][1] = n1->stat2.pointStat.average.data.at(1);

              pts[3][0] = currentNode->stat2.pointStat.average.data.at(0);
              pts[3][1] = currentNode->stat2.pointStat.average.data.at(1);

              pts[6][0] = n2->stat2.pointStat.average.data.at(0);
              pts[6][1] = n2->stat2.pointStat.average.data.at(1);

              // midpoint and a 0.90 percent point
              pts[1][0] = pts[0][0]*(0.5)  + (0.5)*pts[3][0];
              pts[1][1] = pts[0][1]*(0.5)  + (0.5)*pts[3][1];

              pts[2][0] = pts[0][0]*(0.1)  + (0.9)*pts[3][0];
              pts[2][1] = pts[0][1]*(0.1)  + (0.9)*pts[3][1];

              pts[5][0] = pts[6][0]*(0.1)  + (0.9)*pts[3][0];
              pts[5][1] = pts[6][1]*(0.1)  + (0.9)*pts[3][1];

              pts[5][0] = pts[6][0]*(0.5)  + (0.5)*pts[3][0];
              pts[5][1] = pts[6][1]*(0.5)  + (0.5)*pts[3][1];

              //std::cout << "******************************" <<std::endl;

              //
              for(int i = 1; i < 6; i++){
                  derivative[i][0] = pts[i][0];

                  double deltaY = pts[i-1][1] - pts[i+1][1];
                  double deltaX = pts[i-1][0] - pts[i+1][0];
                  double slope = deltaY/deltaX;
                  derivative[i][1] = slope;

                  //std::cout << "First derivative  " << i << " <<" << slope << "..." << deltaY << "," << deltaX <<  std::endl;
              }

              for(int i = 2; i < 5; i++){
                  //second derivative
                  secondDerivative[i][0] = derivative[i][0];

                  double deltaY = derivative[i-1][1] - derivative[i+1][1];
                  double deltaX = derivative[i-1][0] - derivative[i+1][0];
                  double slope = deltaY/deltaX;
                  secondDerivative[i][1] =  slope;
                  //std::cout << "Second derivative  " << i << " <<" << slope << std::endl;

              }
              //The necessary condition is that the second derivative is 0
              // or that the signs of f''(x+ e) and f''(x-e) are different

              //std::cout << "******************************" <<std::endl;
              double sign1 = secondDerivative[2][1]/fabs(secondDerivative[2][1]);
              double sign2 = secondDerivative[4][1]/fabs(secondDerivative[4][1]);


              double sum = fabs(sign1+sign2);//same sign will make 2,different will make 0
              if ( sum < 1){
                  //Possible an inflection...


                  double tmpPt[2];
                  tmpPt[0] = pts[0][0]*(0.5)  + (0.5)*pts[6][0];
                  tmpPt[1] = pts[0][1]*(0.5)  + (0.5)*pts[6][1];

                  double d = sqrt(pow(tmpPt[0]-pts[3][0],2.0) + pow(tmpPt[1]-pts[3][1],2.0));

                  if (d > 0.05)
                     inflections++;
              }


              totalPossible += 1;
          }
          else if ( currentNode->indexNeighborNode.size() > 2){
              inflections+= 1;
              totalPossible += 1;
          }
    }
    return make_pair(inflections,totalPossible);
}



double Thick(std::vector<SkeletonNode*>* graphNodes, double maxRange){
    //Original Def .. The ratio of perimeter to area of a polygon measures... how skinny it is...
    int totalNodes = graphNodes->size();
    std::vector<int> leafs;

    //We have the distances... we can calculate the mayor distance to the point segment
    // and then  compare to the largest length of a path...
    double mayorDistance = 0;
    for(int i = 0; i < totalNodes; i++){
        SkeletonNode* node = graphNodes->at(i);
        if ( node->indexNeighborNode.size() == 1)
            leafs.push_back(i);


        double pos1[3] = {node->stat2.pointStat.average.data.at(0), node->stat2.pointStat.average.data.at(1), 0};
        double dir[2] ={node->direction.data.at(0), node->direction.data.at(1)};
        double ninety = 3.14159*0.5;

        //90° in one direction
        double oneDir[2], otherDir[2];
        oneDir[0] = dir[0]*cos(ninety) - dir[1]*sin(ninety);
        oneDir[1] = dir[0]*sin(ninety) + dir[1]*cos(ninety);
        double pos2[3] = { pos1[0] + oneDir[0]*maxRange, pos1[1] + oneDir[1]*maxRange, 0 };

        double mayor1 = 0;
        for(int k = 0; k < node->stat2.pointStat.values.size(); k++){
            double p[3] = { node->stat2.pointStat.values.at(k)[0], node->stat2.pointStat.values.at(k)[1], 0 };
            double t = CalculateTprojectionPointSegment(pos1, pos2, p);
            if ( t > 0){
                double projectedPoint[3];
                GetProjectedPointSegment(pos1, pos2, p, projectedPoint);
                double d = sqrt(pow(pos1[0] - projectedPoint[0],2.0) + pow(pos1[1] - projectedPoint[1],2.0));
                if ( d > mayor1)
                    mayor1 = d;
            }
        }

        //90° in the other direction
        ninety = -ninety;
        otherDir[0] = dir[0]*cos(ninety) - dir[1]*sin(ninety);
        otherDir[1] = dir[0]*sin(ninety) + dir[1]*cos(ninety);
        double pos3[3] = { pos1[0] + otherDir[0]*maxRange, pos1[1] + otherDir[1]*maxRange, 0 };
        double mayor2 = 0;

        for(int k = 0; k < node->stat2.pointStat.values.size(); k++){
            double p[3] = { node->stat2.pointStat.values.at(k)[0], node->stat2.pointStat.values.at(k)[1], 0 };
            double t = CalculateTprojectionPointSegment(pos1, pos3, p);
            if ( t > 0){
                double projectedPoint[3];
                GetProjectedPointSegment(pos1, pos3, p, projectedPoint);
                double d = sqrt(pow(pos1[0] - projectedPoint[0],2.0) + pow(pos1[1] - projectedPoint[1],2.0));
                if ( d > mayor2)
                    mayor2 = d;
            }
        }

        if (mayor1 + mayor2 > mayorDistance)
            mayorDistance = mayor1 + mayor2;
    }

    //Highest projected distances (i.e. thickness)

    double maxDist = 0;
    for(int i =0 ; i < leafs.size(); i++){
        double d = DFS( leafs.at(i), graphNodes);
        if ( d > maxDist)
            maxDist = d;
    }

   return  mayorDistance/maxDist;

}


#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkConvexHull2D.h>
#include <vtkDelaunay2D.h>
#include <vtkTriangle.h>

vtkSmartPointer<vtkPolyData> GetConvexHull(vector<Float2> *points){
    vtkSmartPointer<vtkPoints> vtpoints =  vtkSmartPointer<vtkPoints>::New();

     for(unsigned int i = 0; i < points->size(); i++){
         vtpoints->InsertNextPoint( points->at(i)[0], points->at(i)[1],0);
     }

    // Add the grid points to a polydata object
    vtkSmartPointer<vtkPolyData> polydata =
      vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(vtpoints);

    // Triangulate the grid points
    vtkSmartPointer<vtkConvexHull2D> delaunay =vtkSmartPointer<vtkConvexHull2D>::New();

    //

    delaunay->SetMinHullSizeInWorld(0);

    delaunay->SetHullShape( vtkConvexHull2D::ConvexHull);
    delaunay->SetInputData(polydata);

    delaunay->OutlineOn();
    delaunay->Update();

    vtkSmartPointer<vtkPolyData> AlphaShape = delaunay->GetOutput(1);

    return AlphaShape;
}


double Norm(double  vector[], int n){
    float size = 0;
    for( int i = 0;i < n; i++) size += vector[i]*vector[i];
    size = sqrt(size);
    return size;
}

double DistanceBtwPoints(double p1[],double p2[]){
    double dv[3] = {p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]};
    return Norm(dv, 3);
}



vtkSmartPointer<vtkPolyData>GetAlphaShape(float alphaValue, vector<Float2> *points){
  // We use the sparness measure to generate the alpha-shape as SparsenessMeasure is q90

    vtkSmartPointer<vtkPoints> vtpoints =
      vtkSmartPointer<vtkPoints>::New();

     for(unsigned int i = 0; i < points->size(); i++){
         vtpoints->InsertNextPoint( points->at(i)[0], points->at(i)[1],0);
     }

    // Add the grid points to a polydata object
    vtkSmartPointer<vtkPolyData> polydata =
      vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(vtpoints);

    // Triangulate the grid points
    vtkSmartPointer<vtkDelaunay2D> delaunay =
    vtkSmartPointer<vtkDelaunay2D>::New();
    delaunay->SetAlpha(alphaValue);
    delaunay->SetInputData(polydata);
    delaunay->Update();

    vtkSmartPointer<vtkPolyData> AlphaShape = delaunay->GetOutput();
    return AlphaShape;
}


double GetPerimeter(vtkSmartPointer<vtkPolyData> polydata){
    // The vtkSelectEnclosed surface checks whether points are inside or not
    // if it is on the border then the it sets as not inside
    // we then use the alpha shapes of the border points which should only generate lines

    vtkSmartPointer<vtkFeatureEdges> featureEdges = vtkSmartPointer<vtkFeatureEdges>::New();
      featureEdges->SetInputData(polydata);
      featureEdges->BoundaryEdgesOn();
      featureEdges->FeatureEdgesOff();
      featureEdges->ManifoldEdgesOff();
      featureEdges->NonManifoldEdgesOff();
      featureEdges->Update();

    vtkSmartPointer<vtkPolyData> boundaryEdges = featureEdges->GetOutput();
    vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
    cellArray = boundaryEdges->GetPolys();

    double totalPerim = 0;
    for(int i = 0; i < boundaryEdges->GetNumberOfCells(); i++){

        vtkCell* cell = boundaryEdges->GetCell(i);

        if ( cell->GetNumberOfPoints() == 2 ){
            vtkLine* line = dynamic_cast<vtkLine*>(cell);
            double p0[3];
            double p1[3];
            line->GetPoints()->GetPoint(0, p0);
            line->GetPoints()->GetPoint(1, p1);

            totalPerim += DistanceBtwPoints(p0,p1);
        }

  }
   return totalPerim;
}


double GetArea(vtkSmartPointer<vtkPolyData> polydata, bool debug=false){
    vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
    cellArray = polydata->GetPolys();


    double totalArea = 0;
    int totalTriangles = 0;
    for(int i = 0; i < polydata->GetNumberOfCells(); i++){

        vtkCell* cell = polydata->GetCell(i);

        if ( cell->GetNumberOfPoints() == 3 ){
            if(debug)std::cout << "New Triangle " <<std::endl;
                vtkTriangle* triangle = dynamic_cast<vtkTriangle*>(cell);
                double p0[3];
                double p1[3];
                double p2[3];
                triangle->GetPoints()->GetPoint(0, p0);
                triangle->GetPoints()->GetPoint(1, p1);
                triangle->GetPoints()->GetPoint(2, p2);
                double area = vtkTriangle::TriangleArea(p0, p1, p2);
                totalArea += area;
                totalTriangles++;
        }
        if (debug) std::cout << cell->GetNumberOfPoints() <<std::endl;

    }

    if (debug) std::cout << "Total triangles " << totalTriangles <<std::endl;
    if ( debug) std::cout << "Number of polys" << polydata->GetNumberOfCells() <<std::endl;
    return totalArea;
}

double GetAreaConvexHull(vtkSmartPointer<vtkPolyData> polydata){
    //Assume it's the outline of the convex hull ... centroid should be inside
    double area = 0;
    double centroid[3] = {0,0,0};
    vtkSmartPointer<vtkPoints> points = polydata->GetPoints();

    for(int i = 0; i < points->GetNumberOfPoints();i++ ){
        double p0[3];
        points->GetPoint(i,p0);
        centroid[0] += p0[0];
        centroid[1] += p0[1];

    }

    centroid[0] /= points->GetNumberOfPoints();
    centroid[1] /= points->GetNumberOfPoints();

    for(int i = 0; i < points->GetNumberOfPoints();i++ ){
        double p0[3], p1[3];
        points->GetPoint(i,p0);

        int next = (i +1)% points->GetNumberOfPoints();
        points->GetPoint(next,p1);

        double currentArea = vtkTriangle::TriangleArea(centroid, p0,p1);
        area += currentArea;
    }
    return area;
}


double ConvexMeasure(vtkSmartPointer<vtkPolyData> alphaShape, vtkSmartPointer<vtkPolyData> convexHull){

    double areaA = GetArea(alphaShape);
    double areaH = GetArea(convexHull);
    double convexity = areaA/areaH;
    return convexity;
}

pair<double, double> ConvexAndSkinny(Image2D<char> thresholded){
    //Use the mask instead...
    std::vector<Float2> points;
    auto blob = thresholded;
    for(int x = 0; x < blob.width(); x++) {
        for(int y = 0; y < blob.height(); y++) {
            if ( blob(x,y) ==1 ){
                Float2 pt{x,y};
                points.push_back(pt);
            }
        }
    }

    /* private double computeAlphaValue() {
        int length = sortedOriginalMSTLengths.length;
        if (length == 0) return 100.;
        int n90 = (9 * length) / 10;
        double alpha = sortedOriginalMSTLengths[n90];
        //The 90th percentile of the edge lenghts is the alpha
        return Math.min(alpha, 100.);
    }
   */
    vtkSmartPointer<vtkPolyData> alpha = GetAlphaShape( 100, &points);
    //TODO- convex hull computation is missing (?)

    if ( points.empty())
        return std::make_pair(0,0);
    double areaA = GetArea(alpha);
    vtkSmartPointer<vtkPolyData> convex = GetConvexHull(&points);
    double areaH = GetAreaConvexHull(convex);
    double convexity = areaA/areaH;
    double skinny = 1.0 - sqrt( (4*3.14159*areaH)) /GetPerimeter(alpha);
    return std::make_pair(convexity, skinny);
}

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

bool ArePointsCollinear(float x1[], float x2[], float x3[]){

    double triangleArea = x1[0]*( x2[1] - x3[1]) +  x2[0]*(x3[1] - x1[1]) + x3[0]*(x1[1] - x2[1]);
    return ( fabs(triangleArea) < 0.001);
}

bool  AreEqual(float p1[], float p2[]){
    double error = 0;
    for(int i = 0; i < 3; i++) error += fabs( p1[i] - p2[i]);

    return (error < 0.01);
}

float AngleBwVectors( float* v1, float *v2, bool degrees){

   gsl_vector *a = gsl_vector_alloc(3);
   gsl_vector *b = gsl_vector_alloc(3);
   //MathHelper::Normalize(v1);
   //MathHelper::Normalize(v2);
   for( int i = 0; i < 3; i++){
      gsl_vector_set(a, i, v1[i]);
      gsl_vector_set(b, i, v2[i]);
   }

   double c = 0;
   gsl_blas_ddot(a,b, &c);
   gsl_vector_free(a);
   gsl_vector_free(b);

   // the dot product is 0....

   float origin[3] = {0,0,0};

   float angle =  acos(c);

   if ( ArePointsCollinear(origin, v1,v2)){

       // if the points are collinear, then either they are
       // going in the same direction and the angle is then 0
       // or they are not equal and they are going in different directions....
      if ( AreEqual(v1,v2))
          angle = 0;
      else
          angle = 3.1415;

   }
   // if the two values are really close, then the result will be 1 and epsilon higher
   // so, if isnan then set to 0 the angle.
   else if ( std::isnan(angle)){
       /*std::cout <<" angle is nan " <<  std::endl;
       std::cout << "v1:: " << v1[0] << ", " << v1[1] << ", "  << v1[2] << std::endl;
       std::cout << "v2:: " << v2[0] << ", " << v2[1] << ", "  << v2[2] << std::endl;
       */
       angle = 0;
   }
   if ( degrees )
       angle = angle * (180.0 / 3.1415);
   return angle;
}


void GetLeafNodes(vector<int>* leafsIdx, std::vector<SkeletonNode*>* localGraphNodes){
    for(int i = 0; i < localGraphNodes->size(); i++){
        SkeletonNode* currentNode = localGraphNodes->at(i);
        if ( currentNode->indexNeighborNode.size() == 1 || currentNode->indexNeighborNode.size() == 0){
            //it is a leaf
            leafsIdx->push_back(i);
        }
    }

}
double Straight(PipelineResult2* result){

   if ( result->numberOfComponents == 1)
       return Straight(&(result->results.at(0)->graphNodes));

   //There are multiply components ...
   //what we need is all the connections between the different components....
   double SOAM = 0;

   int totalEdges = 0;
   int extraEdges = 0;
   //std::cout << "Computing Straight with multiple components " <<  result->numberOfComponents <<  std::endl;
   //For each component,
   int withSkeleton = 0;

   std::set< pair<int,int> > graphConnections;


   for(int cmpt = 0; cmpt < result->numberOfComponents; cmpt++){
       if ( result->results.at(cmpt)->hasSkeleton){
           withSkeleton++;
           // First we search for the leafs
           std::vector<SkeletonNode*>* localGraphNodes = &result->results.at(cmpt)->graphNodes;
           vector<int> localLeafs;
           GetLeafNodes(&localLeafs, localGraphNodes);
           // Now we search all the other graphs... to find which one is the closest....

           double minDistance = 99999999;
           int node1Leaf = -1;
           int graph2 = -1;
           int node2Leaf = -1;

           for(int cmpt2 = 0; cmpt2 < result->numberOfComponents; cmpt2++){
               if (cmpt == cmpt2 ) continue;

               if ( result->results.at(cmpt2)->hasSkeleton){
                   std::vector<SkeletonNode*>* otherGraphNodes = &result->results.at(cmpt2)->graphNodes;
                   vector<int> otherLeafs;
                   GetLeafNodes(&otherLeafs, otherGraphNodes);



                   for(int k = 0; k < localLeafs.size(); k++){
                       SkeletonNode* lNode = localGraphNodes->at(localLeafs.at(k));

                       double x1 = lNode->stat2.pointStat.average.data.at(0);
                       double y1 = lNode->stat2.pointStat.average.data.at(1);


                       for(int j = 0; j < otherLeafs.size(); j++){
                           SkeletonNode* oNode = otherGraphNodes->at(otherLeafs.at(j));

                           double x2 = oNode->stat2.pointStat.average.data.at(0);
                           double y2 = oNode->stat2.pointStat.average.data.at(1);

                           double d = sqrt(pow(x1 - x2, 2.0) + pow(y1-y2,2.0));

                           if ( d < minDistance){
                               node1Leaf = localLeafs.at(k);
                               node2Leaf = otherLeafs.at(j);
                               graph2 = cmpt2;
                               minDistance = d;
                           }

                       }

                   }
               }
           }
           //std::cout << cmpt << " ... " << graph2 << " > local leafs? " << localLeafs.size() <<  ".. total nodes " << localGraphNodes->size() <<  std::endl;


           if (graph2 != -1){
               auto connection =  std::make_pair( min(cmpt,graph2),max(cmpt, graph2));
               if (  graphConnections.count(connection) == 0){
                   graphConnections.insert(connection);
             //      std::cout << "Added a connection "<<  minDistance << std::endl;


                   SkeletonNode* node1 = localGraphNodes->at(node1Leaf);
                   SkeletonNode* node2 = (&result->results.at(graph2)->graphNodes)->at(node2Leaf);

                   double x1 = node1->stat2.pointStat.average.data.at(0);
                   double y1 = node1->stat2.pointStat.average.data.at(1);

                   double x2 = node2->stat2.pointStat.average.data.at(0);
                   double y2 = node2->stat2.pointStat.average.data.at(1);

                   float v1[3] = {node1->direction.data[0], node1->direction.data[1],0};
                   float v2[3] = {node2->direction.data[0], node2->direction.data[1],0};

                   if ( v1[0] < 0.00001 && v1[1] < 0.00001){

                       v1[0] = x2- x1;
                       v1[1] = y2 -y1;

                   }
                   if ( v2[0] < 0.00001 && v2[1] < 0.00001){
                       v2[0] = x1- x2;
                       v2[1] = y1 -y2;

                   }

                   //std::cout << "Dir1 " << v1[0] << "," << v1[1] << std::endl;
                   //std::cout << "Dir2 " << v2[0] << "," << v2[1] << std::endl;

                   float angle = AngleBwVectors(v1,v2,false);
                   //std::cout << "Edge " << i << "/" <<  edges.size() << ":"<<  angle << std::endl;
                   SOAM += (angle/(3.14159));//


                   extraEdges++;

               }
           }
       }
   }


  // std::cout << "edges b/w components "<< extraEdges << " .. with skeleton " << withSkeleton << std::endl;
   //within components ....
   //
   for(int cmpt = 0; cmpt < result->numberOfComponents; cmpt++){
       if ( result->results.at(cmpt)->hasSkeleton){
           std::vector<SkeletonNode*>* graphNodes = &result->results.at(cmpt)->graphNodes;

           std::set< pair<int,int> > edges;
           for(int i = 0; i < graphNodes->size(); i++){
               SkeletonNode* currentNode = graphNodes->at(i);

               for(int k = 0; k < currentNode->indexNeighborNode.size(); k++){
                   int oIdx = currentNode->indexNeighborNode.at(k);
                   edges.insert( std::make_pair( min(i,oIdx),max(i, oIdx)));
               }
           }
          for(auto edge:edges){
               SkeletonNode* node1 = graphNodes->at(edge.first);
               SkeletonNode* node2 = graphNodes->at(edge.second);

               float v1[3] = {node1->direction.data[0], node1->direction.data[1],0};
               float v2[3] = {node2->direction.data[0], node2->direction.data[1],0};
               //std::cout << "V1 " << v1[0] <<"," << v1[1] << std::endl;
               //std::cout << "V2 " << v2[0] <<"," << v2[1] << std::endl;

               float angle = AngleBwVectors(v1,v2,false);
               //std::cout << "Edge " << i << "/" <<  edges.size() << ":"<<  angle << std::endl;
               SOAM += (angle/(3.14159));//
               totalEdges++;
           }
       }
   }

   if ( (totalEdges + extraEdges) ==0 )   return 0;

   SOAM /= (totalEdges + extraEdges);
   //How to normalize...
   return (1.0 - SOAM);
}



double Straight(std::vector<SkeletonNode*>* graphNodes){
    double SOAM = 0;

    std::set< pair<int,int> > edges;
    for(int i = 0; i < graphNodes->size(); i++){
        SkeletonNode* currentNode = graphNodes->at(i);

        for(int k = 0; k < currentNode->indexNeighborNode.size(); k++){
            int oIdx = currentNode->indexNeighborNode.at(k);
            edges.insert( std::make_pair( min(i,oIdx),max(i, oIdx)));
        }
    }

    int i = 0;
    for(auto edge:edges){
        SkeletonNode* node1 = graphNodes->at(edge.first);
        SkeletonNode* node2 = graphNodes->at(edge.second);

        float v1[3] = {node1->direction.data[0], node1->direction.data[1],0};
        float v2[3] = {node2->direction.data[0], node2->direction.data[1],0};
        //std::cout << "V1 " << v1[0] <<"," << v1[1] << std::endl;
        //std::cout << "V2 " << v2[0] <<"," << v2[1] << std::endl;

        float angle = AngleBwVectors(v1,v2,false);
        //std::cout << "Edge " << i << "/" <<  edges.size() << ":"<<  angle << std::endl;
        SOAM += (angle/(3.14159));//
        i++;
    }

    if ( edges.size() == 0) return 0;
    SOAM /= edges.size();
    //How to normalize...

    return (1.0 - SOAM);
}


double Pearson(PipelineResult2* result){
    std::vector<SkeletonNode*>* graphNodes = new std::vector<SkeletonNode*>();
    for(int i = 0; i < result->results.size(); i++){
        std::vector<SkeletonNode*>* localGraphNodes = &result->results.at(i)->graphNodes;

        for(int k = 0; k < localGraphNodes->size(); k++){
            graphNodes->push_back( localGraphNodes->at(k));
        }
    }

    double meanX = 0, meanY = 0;
    float total = 0;
    for(int j = 0; j < graphNodes->size(); j++){
        auto node = graphNodes->at(j);
        double xi = node->stat2.pointStat.average.data.at(0);
        double yi = node->stat2.pointStat.average.data.at(1);
        meanX += xi;
        meanY += yi;
        total++;

    }

    meanX /= total;
    meanY /= total;

    double up = 0;
    double den1 = 0;
    double den2 = 0;

    for(int j = 0; j <  graphNodes->size(); j++){
        auto node = graphNodes->at(j);

        double xi = node->stat2.pointStat.average.data.at(0);
        double yi = node->stat2.pointStat.average.data.at(1);

         up += (xi - meanX)*(yi - meanY);
         den1 += (xi - meanX)*(xi- meanX);
         den2 += (yi - meanY)*(yi- meanY);
    }

    double den = sqrt(den1)*sqrt(den2);

    double r = up/ den;

    if ( den < 0.0001)
        r = 0;

    return r;

}

int compareLocalPoint(const void *a, const void * b){

    float* p1 = (float*)a;
    float* p2 = (float*)b;

    // p1 has x,y, color, bucket
    float v1 = p1[3];// We have already made so we know the bucket
    float v2 = p2[3];
    if ( v1 < v2) return -1;
    if ( v1 == v2) return 0;

    return 1;
}


double QueryQuadTree(float px, float py, float* quadTreeMetadata, float* mappedPoints,int maxLevel, int currentIdx){

    // qx,qy is the point to see whether it is close by
    // px,py is the point to query
    const int PER_POINT_INFO = 6;

    int startIndex = 0;
    int endIndex = 0;
    int currentPos = 0;
    int actualPos = currentPos;

    for(int k = 0; k < maxLevel - 1; k++){
         //  We start with the first block, and find which block it belongs to, or
         // whether there are no points left
        int start = (1- pow(4, k)) / (1-4);
        actualPos = currentPos + start;
        float xminP = quadTreeMetadata[actualPos*6 + 0];
        float xmaxP = quadTreeMetadata[actualPos*6 + 1];
        float yminP = quadTreeMetadata[actualPos*6 + 2];
        float ymaxP = quadTreeMetadata[actualPos*6 + 3];
        int parentStartIndex = quadTreeMetadata[actualPos*6 + 4];
        int parentEndIndex = quadTreeMetadata[actualPos*6 + 5];

        int xIndex = 2.0*(px - xminP)/(xmaxP-xminP);
        int yIndex = 2.0*(py - yminP)/(ymaxP-yminP);
        int bucketIndex = 2*yIndex + xIndex;
        currentPos = currentPos*4 + bucketIndex;
        if ( parentStartIndex == parentEndIndex){
            startIndex = parentStartIndex;
            endIndex = parentEndIndex;
            break;
        }
        start = (1- pow(4, k+1)) / (1-4);
        startIndex = quadTreeMetadata[(currentPos + start)*6 + 4];
        endIndex = quadTreeMetadata[(currentPos+ start)*6 + 5];
    }

    // No points in the area,
    // Only if there are points we add color, otherwise we leave it as it is...

    // In the no blending, first it finds the one in the current cell,
    // and then it goes to the others cell, that's the reason to the overlap
    // we need a way to define on how we select the points, such that it always select the same...

    double closestDistance = DBL_MAX;
    if ( startIndex != endIndex) {
        //if this is the case, then there are points in the level...

         // We went to the last level, and now we look at the points
        for(int k = startIndex ; k < endIndex; k++){
              float x = mappedPoints[k*PER_POINT_INFO + 0];
              float y = mappedPoints[k*PER_POINT_INFO + 1];
              //int val = mappedPoints[k*6 + 2];
              int idx = mappedPoints[k*PER_POINT_INFO + 4];

              if ( currentIdx != idx){
                  double d = sqrt(pow(px - x,2.0) + pow(py-y,2.0));
                  closestDistance = std::min(closestDistance,d );
              }
        }
    }

    return closestDistance;
}






double ComputeSigmaQuadTree(const DataSet& dataSet, const Projection2D& projection){
    //The radius R is set to the average distance
    // δ of a point in S to its nearest- non-zero neighbor.
    // use a quadtree to get the points

    int totalNotCulled = 0;
    int totalElements = dataSet.size();
    const int PER_POINT_INFO = 6;
    float* mappedPoints = (float*)malloc(sizeof(float)*totalElements*PER_POINT_INFO); // x, y, valueIndex, bucketIndex, initialIndex, colorIndex
    float xmin =  0,ymin = 0, xmax = 1.0, ymax = 1.0;

    float DesiredRadius = 0.01;

    for(int k = 0; k < totalElements; k++){

        auto pair = projection(dataSet, k);

        float v1 = pair[0];
        float v2 = pair[1];

        float mx = v1;
        float my = v2;

        mappedPoints[totalNotCulled*PER_POINT_INFO + 0] = mx;
        mappedPoints[totalNotCulled*PER_POINT_INFO + 1] = my;
        mappedPoints[totalNotCulled*PER_POINT_INFO +2] = 0;
        mappedPoints[totalNotCulled*PER_POINT_INFO +5] = 0;

        int xIndex = 2.0*(mx - xmin)/(xmax-xmin);
        int yIndex = 2.0*(my - ymin)/(ymax-ymin);

        mappedPoints[totalNotCulled*PER_POINT_INFO + 3] =  2*yIndex + xIndex; // bucket index
        mappedPoints[totalNotCulled*PER_POINT_INFO + 4] = k; // initial index
        totalNotCulled++;
    }

    int maxLevel = 5;
    int totalInLastLevel = sqrt(pow(4, maxLevel));
    float midDistance = (1.0/totalInLastLevel); // half the size of each division...
    //float minD = 3.0/128.0;

    while ( DesiredRadius*2 > midDistance){
        maxLevel -= 1;
        totalInLastLevel = sqrt(pow(4, maxLevel));
        midDistance = (1.0/totalInLastLevel);
    }

    int totalSize =  (1 - pow(4,maxLevel+1))/(1-4);
    float* quadTreeMetadata = (float*) malloc(sizeof(float)*6*totalSize);

    //std::cout << "Start of metadata " << xmin << ", " << ymin << " -> " << xmax << " , " << ymax << std::endl;

    quadTreeMetadata[0] = xmin;           // Complete xmin
    quadTreeMetadata[1] = xmax;           // Complete xmax
    quadTreeMetadata[2] = ymin;           // Complete ymin
    quadTreeMetadata[3] = ymax;           // Complete ymax
    quadTreeMetadata[4] = 0;              // starting index of points
    quadTreeMetadata[5] = totalNotCulled; // end index of points :: total points not culled are the ones that we are drawing


    qsort(mappedPoints,totalNotCulled,sizeof(float)*PER_POINT_INFO, compareLocalPoint);
    //Then we create each level of the quadtree and set the points

    for(int level = 1; level < (maxLevel+1); level++){
        //For each level...
        //We work based on the previous level...

        int totalInLevel = pow(4, level);
        int prevStart = (1- pow(4, level-1)) / (1-4);
        int start = (1- pow(4, level)) / (1-4);

        for(int i = 0; i < totalInLevel; i++){
            // This gives the total items in that level...
            // First we need to figure out where it belongs...
            int cellInPrevLevel = i / 4;
            int loc = cellInPrevLevel * 4 + start +  i%4; // Where to save it...
            int parentLoc = prevStart + cellInPrevLevel;

            float xminParent = quadTreeMetadata[parentLoc*6 + 0];
            float xmaxParent = quadTreeMetadata[parentLoc*6 + 1];
            float yminParent = quadTreeMetadata[parentLoc*6 + 2];
            float ymaxParent = quadTreeMetadata[parentLoc*6 + 3];

            int parentStartIndex = quadTreeMetadata[parentLoc*6 + 4];
            int parentEndIndex = quadTreeMetadata[parentLoc*6 + 5];

            int currentStart = parentStartIndex;
            while( mappedPoints[currentStart*PER_POINT_INFO +3] < static_cast<float>(i%4)){
                currentStart++;
                if ( currentStart > parentEndIndex){ currentStart = parentEndIndex; break; }
            }

            int currentEnd = currentStart;
            while( mappedPoints[currentEnd*PER_POINT_INFO + 3] < static_cast<float>(i%4 + 1)){
                currentEnd++;
                if ( currentEnd > parentEndIndex){ currentEnd = parentEndIndex;  break; }
            }



            float midX = 0.5*(xmaxParent - xminParent) + xminParent;
            float midY = 0.5*(ymaxParent - yminParent) + yminParent;

            float valuesX[3] = {xminParent, midX, xmaxParent};
            float valuesY[3] = {yminParent, midY, ymaxParent};
                                  //xlow, xhigh, ylow, yigh
            int positions[4][4] = { {0,1,0,1},
                                    {1,2,0,1},
                                    {0,1,1,2},
                                    {1,2,1,2}};


           quadTreeMetadata[loc*6 + 0] = valuesX[ positions[i%4][0]];
           quadTreeMetadata[loc*6 + 1] = valuesX[ positions[i%4][1]];
           quadTreeMetadata[loc*6 + 2] = valuesY[ positions[i%4][2]];
           quadTreeMetadata[loc*6 + 3] = valuesY[ positions[i%4][3]];
           quadTreeMetadata[loc*6 + 4] = currentStart;
           quadTreeMetadata[loc*6 + 5] = currentEnd;
        }

        // Now we sort the area
        for(int i = 0; i <  totalInLevel; i++){
            int cellInPrevLevel = i / 4;
            int loc = cellInPrevLevel * 4 + start +  i%4; // Where to save it...

            int currentStart = quadTreeMetadata[loc*6 + 4];
            int currentEnd = quadTreeMetadata[loc*6 + 5];
            for(int k = currentStart; k < currentEnd; k++ ){
                int xIndex = 2.0*(mappedPoints[k*PER_POINT_INFO +0]  - quadTreeMetadata[loc*6 +0])/(quadTreeMetadata[loc*6 +1]-quadTreeMetadata[loc*6 +0]);
                int yIndex = 2.0*(mappedPoints[k*PER_POINT_INFO +1]  - quadTreeMetadata[loc*6 +2])/(quadTreeMetadata[loc*6 +3]-quadTreeMetadata[loc*6 +2]);

                mappedPoints[k*PER_POINT_INFO + 3] =  2*yIndex + xIndex; // bucket index
            }
            qsort(&(mappedPoints[currentStart*PER_POINT_INFO]), currentEnd-currentStart,sizeof(float)*PER_POINT_INFO, compareLocalPoint );
        }
    }



    double sigma = 0;

    int totalDBLMax = 0;
    //The quadtree is generated... now we want the closest point distance . ...
    for(int k = 0; k < totalElements; k++){
        auto pair = projection(dataSet, k);
        double closestD = QueryQuadTree(pair[0], pair[1], quadTreeMetadata, mappedPoints,maxLevel, k );

        if (closestD != DBL_MAX)
             sigma += closestD;
        else
            totalDBLMax++;
    }

    //std::cout << "Total DBL Max " << totalDBLMax << " / " << totalElements << std::endl;
    sigma = (sigma)/(totalElements-totalDBLMax);

    double defaultValue = 0.025;
    if (isnan(sigma))
        sigma = defaultValue;

    free(quadTreeMetadata);
    free(mappedPoints);
    if ( sigma < defaultValue)
        sigma = defaultValue;

    if ( sigma > 0.1)
        sigma = 0.1;

    return sigma;
}




