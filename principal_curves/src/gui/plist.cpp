#include "plist.hpp"

#include "global.hpp"
#include "projectionmatrix.hpp"
#include "../mds.hpp"
#include "../curve2d.hpp"
#include "../draw.hpp"
#include "../tsne/tsne.hpp"
#include "../principalgraph.hpp"

#include <cstdlib>
#include <cmath>

void projectionListReset() {
    global.projectionList.clear();
    global.projectionListCountLabel->set_text("0");
}

void projectionListExtend() {
    if(global.projectionList.size() == 0) {
        global.projectionListCardinalCount = 0;
        // if list empty
        // add carndinal projections
        for(int i = 0; i < global.dataSet.dimension(); i++) {
            if(global.dataSet.getStat(i).minValue == global.dataSet.getStat(i).maxValue) continue;

            for(int j = i + 1; j < global.dataSet.dimension(); j++) {
                if(global.dataSet.getStat(j).minValue == global.dataSet.getStat(j).maxValue) continue;

                ProjectionMatrix2D matrix { global.dataSet.dimension() };
                matrix.entry(i) = Float2({1, 0});
                matrix.entry(j) = Float2({0, 1});
                //auto result = processPipelineNormal(global.parameters, global.dataSet, calibrateMatrix(&matrix, global.dataSet));
                double time;
                auto result = processPipeline(global.parameters, global.dataSet, calibrateMatrix(&matrix, global.dataSet), global.generator, &time);

                ProjectionEntry entry;
                entry.projectionMatrix = std::move(matrix);

                if ( global.parameters.asGraph){
                    std::vector<AugmentedPoint> resCurve;
                    for(int i = 0; i < result.graphNodes.size(); i++){
                        auto node = result.graphNodes.at(i);

                        AugmentedPoint pt;
                        pt.amountPerLength = node->amountPerLength;
                        pt.orthoVariance = node->orthoVariance;
                        pt.point = Int2{ node->location[0], node->location[1]};
                        resCurve.push_back(pt);
                    }
                    entry.curve = std::move(resCurve);

                }
                else {
                    entry.curve = std::move(result.augmentedCurve);
                }
                global.projectionList.emplace_back(std::move(entry));

                global.projectionListCardinalCount++;
            }
        }
    }

    for(int i = 0; i < 10; i++) {
        ProjectionMatrix2D matrix = global.projectionMatrixGenerator.generate(global.dataSet);
        double time;
        auto result = processPipeline(global.parameters, global.dataSet, calibrateMatrix(&matrix, global.dataSet), global.generator, &time);
        //auto result = processPipelineNormal(global.parameters, global.dataSet, calibrateMatrix(&matrix, global.dataSet));
        ProjectionEntry entry;
        entry.projectionMatrix = std::move(matrix);

        if ( global.parameters.asGraph){
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
            entry.curve = std::move(resCurve);
        }
        else {
             entry.curve = std::move(result.augmentedCurve);
        }
        global.projectionList.emplace_back(std::move(entry));
    }
    global.projectionListCountLabel->set_text(std::to_string(global.projectionList.size()));
}

double *projectionListMap_SM() {
    int n = global.projectionList.size();
    double *similarityMatrix = (double*) malloc(n * n * sizeof(double));

    std::cout << "Projection List Computation? " << std::endl;

    bool useFrechet = !global.parameters.asGraph;
    for(int x = 0; x < n; x++) {
        for(int y = x; y < n; y++) {
            float val = 0;
            if(x != y) {

                    if (useFrechet){
                        std::cout <<"Frechet" << std::endl;
                        val = couplingDistanceAugmented(global.projectionList[x].curve, global.projectionList[y].curve, global.adjustmentVarianceFactor->get_value(), global.adjustmentAmountFactor->get_value(), global.checkMirror->get_active());
                    }
                    else {
                        std::cout << "Hausdorff " <<  global.parameters.height <<  std::endl;
                        val = WeightedHausdorffDistance(global.projectionList[x].curve, global.projectionList[y].curve, global.adjustmentAmountFactor->get_value(),  global.adjustmentVarianceFactor->get_value(), global.parameters.height );
                    }

                }

            if(std::isnan(val)) val = 1;

            similarityMatrix[x + y * n] = similarityMatrix[y + x * n] = val;
        }
    }

    return similarityMatrix;
}

void projectionListMap() {
    int n = global.projectionList.size();
    double *similarityMatrix = projectionListMap_SM();

    std::vector<std::pair<float, float>> dots;
    CalculateMDS(similarityMatrix, n, &dots);
    free(similarityMatrix);

    global.dots.resize(0);
    global.dots.reserve(dots.size());
    for(int i = 0; i < dots.size(); i++) {
        global.dots.emplace_back(dots[i].first, dots[i].second);
    }

    projectionListMapDrawOn(&global.projectionListImageController);
    global.projectionListWindow->show();
}

void projectionListMapTSNE() {
    int n = global.projectionList.size();
    double *similarityMatrix = projectionListMap_SM();

    TSNE tsne;
    double *result;

    tsne.Q_run(similarityMatrix, n, result, 2, 10, 42, false);

    global.dots.resize(0);
    global.dots.reserve(n);

    Stat<Float2> stat;
    for(int i = 0; i < n; i++) {
        Float2 val { result[2*i], result[2*i + 1] };
        stat.put(val);
        global.dots.push_back(val);
    }
    for(int i = 0; i < n; i++) {
        global.dots[i] = (global.dots[i] - stat.minValue) / (stat.maxValue - stat.minValue);
    }

    free(result);
    free(similarityMatrix);

    projectionListMapDrawOn(&global.projectionListImageController);
    global.projectionListWindow->show();
}

void projectionListMapDrawOn(ImageController *imgc) {
    auto zoom = imgc->getZoom();

    auto surface = drawDots(global.dots, 20 * zoom, 20 * zoom, global.dotsHighlight, [](int index) { 
        return index < global.projectionListCardinalCount;
    });

    imgc->setImage(surface);
    cairo_surface_destroy(surface);
}

void projectionListMapClick(float x, float y) {
    Float2 click = Float2({x, y});
    int min_i = 0;
    float minDistanceSquared = INFINITY;

    for(int i = 0; i < global.dots.size(); i++) {
        Float2 pos = global.dots[i];
        float distanceSq = vectorLengthSquared(pos - click);
        if(minDistanceSquared > distanceSq) {
            min_i = i;
            minDistanceSquared = distanceSq;
        }
    }

    if(min_i > global.projectionList.size()) return;

    projectionMatrixSet(global.projectionList.at(min_i).projectionMatrix);
    global.dotsHighlight = min_i;

    projectionListMapDrawOn(&global.projectionListImageController);
}

