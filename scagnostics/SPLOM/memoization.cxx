#include "memoization.h"

Memoization::Memoization(int attr1, int attr2, int attr, int filt, int size, short plane[], vector<LocalPoint> points, float _thresInfo[])
{

    attributesDrawn = std::make_pair(attr1, attr2);
    clusteringAttribute = attr;
    filteringAttribute = filt;
    usedResolution = size;

    data = (short*) malloc(usedResolution*usedResolution*3*sizeof(short));

    memcpy(&(data[0]), &(plane[0]), sizeof(short)*usedResolution*usedResolution*3);

    for(unsigned int i = 0; i < points.size(); i++)
        centerlinesPoints.push_back(points.at(i));

    for(int i = 0; i< 5; i++)
        thresInfo[i] = _thresInfo[i];
}

void Memoization::GetCenterlinePoints(vector<LocalPoint>* points){

    for(unsigned int i = 0; i < centerlinesPoints.size(); i++)
        points->push_back(centerlinesPoints.at(i));
}

void Memoization::GetSimplifiedCenterlinePoints(vector<LocalPoint>* points){
    for(unsigned int i = 0; i < simplifiedCenterlines.size(); i++)
        points->push_back(simplifiedCenterlines.at(i));
}

void Memoization::GetLargestCenterlinePoints(vector<LocalPoint>* points){
    for(unsigned int i = 0; i <  largestCenterline.size(); i++){

        LocalPoint c = simplifiedCenterlines.at(largestCenterline.at(i));
        LocalPoint imageP;
        imageP.p[0] = c.p[1] / usedResolution;
        imageP.p[1] = 1.0 - c.p[0] / usedResolution;
        points->push_back(imageP);
    }
}

void Memoization::GetSkeletonSegments(vector<SkeletonSegment*>* segments, double* maxD, double* minD, double *maxDistance){

    for(unsigned int i = 0; i < segmentedCenterlines.size(); i++)
        segments->push_back(segmentedCenterlines.at(i));
    *maxD = maxDensity;
    *minD = minDensity;
    *maxDistance = this->maxDistance;
}

Memoization::~Memoization(){
   centerlinesPoints.clear();
   simplifiedCenterlines.clear();
   largestCenterline.clear();

   for(unsigned int i = 0; i < segmentedCenterlines.size(); i++){
       SkeletonSegment* current = segmentedCenterlines.at(i);
       delete current;
   }
   segmentedCenterlines.clear();
   free(data);
}

void Memoization::GetSkeleton(short texture[]){


    memcpy(&(texture[0]), data, sizeof(short)*usedResolution*usedResolution*3);
    /*for(int i= 0; i < usedResolution*usedResolution*3; i++){
        texture[i] = data[i];
    }*/

}



bool Memoization::IsSame(int attr1, int attr2, int attr, int filt, int size, float _thresInfo[]){

    if ( size != usedResolution)
        return false;

    if ( attr != clusteringAttribute)
        return false;

    if ( filt != filteringAttribute)
        return false;

    if ( attributesDrawn.first != attr1)
        return false;
    if ( attributesDrawn.second != attr2)
        return false;

    for(int i = 0; i < 5; i++){
         if ( fabs(_thresInfo[i] - thresInfo[i]) > 0.001 )
             return false;
    }

    return true;
}

void Memoization::SetSimplification(vector<SkeletonSegment *> *skeleton, vector<LocalPoint>* simplified, vector<int> *largest, double maxD, double minD, double maxDistance){

    simplifiedCenterlines.clear();
    segmentedCenterlines.clear();
    largestCenterline.clear();

    for(unsigned int i = 0; i < simplified->size(); i++)
        simplifiedCenterlines.push_back(simplified->at(i));

    for(unsigned int i = 0; i < skeleton->size(); i++)
        segmentedCenterlines.push_back(skeleton->at(i));

    for(unsigned int i = 0; i <largest->size(); i++)
        largestCenterline.push_back(largest->at(i));
    maxDensity = maxD;
    minDensity = minD;
    this->maxDistance = maxDistance;
}
