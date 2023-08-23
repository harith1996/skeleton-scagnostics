#include "transform2d.hpp"
#include "skeleton.hpp"
#include "multival.hpp"
#include "principalcurve.hpp"
#include "curve2d.hpp"
#include "plotanalysis.hpp"
#include "skeletongenerator.h"
#include "draw.hpp"


#include <boost/bind.hpp>

#ifndef PRINCIPALGRAPH_HPP
#define PRINCIPALGRAPH_HPP




int GetNumberOfComponents(Image2D<char> blob,Image2D<short int>& components,vector<int>* sizes);

void GetNonVistedNeighbors(int idx, vector<int>* neighbors,  vector<Int2>* skeletonPoints, bool visited[]);

void GetNeighbors(int idx, vector<int>* neighbors,  vector<Int2>* skeletonPoints);

void GetParentsInGraph(SkeletonNode* node, vector<SkeletonNode*>* principalGraph, set<int>* parents);

void GetChildrenInGraph(SkeletonNode* node, vector<SkeletonNode*>* principalGraph, set<int>* children, set<int>* visited);

bool NodeAlreadyInParents(SkeletonNode* node, SkeletonNode* other);

bool NodeAlreadyInChildren(SkeletonNode* node, SkeletonNode* other);

void CreateSkeletonGraph(int controlPointSamplingSize,  vector<int>& endPoints,  vector<int>& bifurcations,  vector<Int2>* skeletonPoints,
                           vector<SkeletonNode*>* principalGraph);

void CreateSkeletonGraphAutomaticPoints(vector<int>& endPoints,  vector<int>& bifurcations,  vector<Int2>* skeletonPoints,
                                          vector<SkeletonNode*>* principalGraph, Image2D<float> &dt, double w = 0.2);




double HausdorffDistance(std::vector<Int2> pts1, std::vector<Int2> pts2, int sz);

double HausdorffDistance(PipelineResult *result, std::vector<AugmentedPoint>* generatingCurve, int sz);

int principalGraph(vector<SkeletonNode *>* nodes, const Image2D<unsigned>& plot, int iterationCount, float convergenceDistance, const float arcSegment, int smoothCount, bool debug=false );

void augmentGraph(vector<SkeletonNode *> *nodes, const Image2D<unsigned>& plot, int smoothCount);

Image2D<float> GetBands(vector<SkeletonNode *> *nodes, int sz, bool resample= true);


double WeightedHausdorffDistance(const std::vector<AugmentedPoint>& curve1, const std::vector<AugmentedPoint>& curve2, float densityWeight, float varianceWeight, int sz);

double WeightedHausdorffDistanceBase(const std::vector<AugmentedPoint>& curve1, const std::vector<AugmentedPoint>& curve2, float densityWeight, float varianceWeight, int sz);

#endif // PRINCIPALGRAPH_HPP
