#ifndef MEMOIZATION_H
#define MEMOIZATION_H

/*
 * The scatterplot is actually fast enough... let's just worry about the
 * Centerlines ...
 * I need to create a memory class, which basically stores ...
   1. resolution used,
   2. clustering attribute.
   3. filtering attribute
   4. which attributes were drawn
   If we are dealing with centerline then we don't need the typeOfBlend...
*/

#include <chrono>
#include <vector>
#include <iostream>
#include <utility>
#include <cstring>
#include "dataset.h"



class SkeletonSegment {
    public:
       // The points until a bifurcation, with their index
       // of the modified centerline
       vector<pair<LocalPoint, int> > pointsInSegment;

       //  the index in the modified centerline
       vector<int> parentIdx;
       // Image coordinate points
       vector<LocalPoint> imageCoorPoints;
       //
       vector<LocalPoint> tangents;
       //
       vector<double> densities;
       //
       vector<double> distances;

       vector<LocalPoint> debugQueryPoints;
       int indexOfSegment;
       int indexOfParentSegment;
};


class Memoization
{
public:
    Memoization(int attr1, int attr2, int attr, int filt, int size, short plane[], vector<LocalPoint> points, float _thresInfo[]);

    ~Memoization();
    bool IsSame(int attr1, int attr2, int attr, int filt, int size, float _thresInfo[]);

    void GetSkeleton(short texture[]);
    void GetCenterlinePoints(vector<LocalPoint>* points);

    void GetSimplifiedCenterlinePoints(vector<LocalPoint>* points);
    void GetLargestCenterlinePoints(vector<LocalPoint>* points);

    void GetSkeletonSegments(vector<SkeletonSegment*>* segments, double *maxD, double *minD, double *maxDistance);

    void SetSimplification(vector<SkeletonSegment*>* skeleton, vector<LocalPoint>* simplified, vector<int>* largest, double maxD, double minD, double maxDistance);

private:
    double lastTimeUsed; // In order to remove it, we shall use the memory in a removing the ones not used

    int usedResolution;
    int clusteringAttribute;
    int filteringAttribute;
    std::pair<int, int> attributesDrawn;
    short* data; // The actual data

    vector<LocalPoint> centerlinesPoints;
    vector<LocalPoint> simplifiedCenterlines;
    vector<SkeletonSegment*> segmentedCenterlines;
    vector<int> largestCenterline;
    double maxDensity;
    double minDensity;
    double maxDistance;

    float thresInfo[5];
};

#endif // MEMOIZATION_H
