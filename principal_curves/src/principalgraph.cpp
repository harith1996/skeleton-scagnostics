#include "principalgraph.hpp"
#include "propagate2d.hpp"


void GetNonVistedNeighbors(int idx, vector<int>* neighbors,  vector<Int2>* skeletonPoints, bool visited[]){
    neighbors->clear();
    auto currentLocation = skeletonPoints->at(idx);
    for(unsigned int j = 0; j < skeletonPoints->size(); j++){
        if (idx == static_cast<int>(j)) continue;

        int xDistance = abs(currentLocation[0]- skeletonPoints->at(j)[0]);
        int yDistance = abs(currentLocation[1]- skeletonPoints->at(j)[1]);

        if ( xDistance <= 1 && yDistance <= 1 && !visited[j]){
           neighbors->push_back(j);
        }
    }
}

void GetNeighbors(int idx, vector<int>* neighbors,  vector<Int2>* skeletonPoints){
    neighbors->clear();
    auto currentLocation = skeletonPoints->at(idx);
    for(unsigned int j = 0; j < skeletonPoints->size(); j++){
        if (idx == static_cast<int>(j)) continue;

        int xDistance = abs(currentLocation[0]- skeletonPoints->at(j)[0]);
        int yDistance = abs(currentLocation[1]- skeletonPoints->at(j)[1]);

        if ( xDistance <= 1 && yDistance <= 1 ){
           neighbors->push_back(j);
        }
    }
}

void DrawVoronoiGraph(int iter,  Image2D<int> labels, std::vector<Int2>& points, const Image2D<unsigned>& plot, vector<SkeletonNode *> *nodes){
    DrawContext dc { static_cast<int>(labels.width()), static_cast<int>(labels.height()), 600, 600 };
    dc.drawVoronoiRegions(labels, plot, nodes->size());


    std::stringstream ss;
    ss << "voronoi_" << iter << ".png";
    dc.setColor(1.0,1.0,0.0);
    dc.drawVoronoiGraph(nodes);
    //dc.drawMainLine(points);
    dc.setColor(0.0, 1.0,1.0);
    dc.drawControlPoints(points);
    dc.writeToFile(ss.str());

}

int GetNumberOfComponents(Image2D<char> blob, Image2D<short>& components, vector<int>* sizes){
    typedef itk::Image<short int, 2 > HImageType;
    HImageType::Pointer img = HImageType::New();
    HImageType::RegionType region;
    HImageType::IndexType index;
    index[0] = 0;      index[1] = 0;
    HImageType::SizeType sizeT;

    sizeT[0] =  blob.width();
    sizeT[1] =  blob.height();
    region.SetIndex(index);
    region.SetSize(sizeT);
    img->SetRegions(region);
    img->Allocate(true);
    img->FillBuffer(itk::NumericTraits< HImageType::PixelType >::Zero);
    img->Update();
    itk::ImageRegionIteratorWithIndex<HImageType> it(img, region);

    while(!it.IsAtEnd()){

        auto idx = it.GetIndex();
        int x = idx[0];
        int y = idx[1];

        it.Set( blob(x,y) );
        ++it;
    }


    typedef itk::ConnectedComponentImageFilter <HImageType,HImageType >    ConnectedComponentImageFilterType;
    ConnectedComponentImageFilterType::Pointer connected =   ConnectedComponentImageFilterType::New ();
    connected->SetInput(img);
    connected->Update();

    typedef itk::RelabelComponentImageFilter<HImageType,HImageType> RelabelFilterType;

    RelabelFilterType::Pointer relabel = RelabelFilterType::New();
    relabel->SetInput(connected->GetOutput());
    relabel->SortByObjectSizeOn();
    relabel->Update();

    HImageType::Pointer labeled = relabel->GetOutput();
    itk::ImageRegionIteratorWithIndex<HImageType> it2(labeled, region);

    while(!it2.IsAtEnd()){

        auto idx = it2.GetIndex();
        int x = idx[0];
        int y = idx[1];

        components(x,y) = it2.Get();
       ++it2;
    }

    for(int j = 0; j< relabel->GetNumberOfObjects(); j++){
        auto size = relabel->GetSizeOfObjectInPixels(j+1);
        sizes->push_back(size);
    }

    img = NULL;
    int numComponents = relabel->GetNumberOfObjects();
    return numComponents;
}

void GetParentsInGraph(SkeletonNode* node, vector<SkeletonNode*>* principalGraph, set<int>* parents){

    if ( node->idxInGraph != -1){
        parents->insert(node->idxInGraph);
        return;
    }
    else {
        for(int k = 0; k < node->parents.size(); k++){
            GetParentsInGraph(node->parents.at(k), principalGraph, parents);
        }
    }
}

void GetChildrenInGraph(SkeletonNode* node, vector<SkeletonNode*>* principalGraph, set<int>* children, set<int>* visited){

    visited->insert(node->label);
    if ( node->idxInGraph != -1){
        children->insert(node->idxInGraph);
        return;
    }
    else {
        for(int k = 0; k < node->children.size(); k++){

            if ( count(visited->begin(), visited->end(),  node->children.at(k)->label ) == 0 )
               GetChildrenInGraph(node->children.at(k), principalGraph, children, visited);
        }
    }
}

bool NodeAlreadyInParents(SkeletonNode* node, SkeletonNode* other){
   bool exist = false;

   for(int i = 0; i <node->parents.size(); i++)
       if ( other->label == node->parents.at(i)->label)
           exist = true;
   return true;
}


bool NodeAlreadyInChildren(SkeletonNode* node, SkeletonNode* other){
   bool exist = false;

   for(int i = 0; i <node->children.size(); i++)
       if ( other->label == node->children.at(i)->label)
           exist = true;
   return true;
}


void CreateSkeletonGraphAutomaticPoints(vector<int>& endPoints,  vector<int>& bifurcations,  vector<Int2>* skeletonPoints,
                                          vector<SkeletonNode*>* principalGraph,  Image2D<float> &dt, double w ){

    // > EndPoints and Bifurcations
    // > Sort the rest of skeleton points according to their distance transform

    // > two points (?) are connected iff the maximally-inscribed balls
    //   centered at x and y intersect ( ||x-y|| <  a(DT(x) + DT(y)) )

    float maxDt = 0;
    vector<SkeletonNode*> allNodes;
    vector<int> neighbors;
    vector< pair<float, int>> dtValues;
    // Create a new node for each
    for(unsigned int i = 0; i < skeletonPoints->size(); i++){
        SkeletonNode* newNode = new SkeletonNode;
        newNode->label = i;
        newNode->location = skeletonPoints->at(i);
        newNode->initialDt = dt(newNode->location[0],newNode->location[1]);
        if ( newNode->initialDt > maxDt)
            maxDt = newNode->initialDt;
        dtValues.push_back(make_pair( newNode->initialDt, i));
        allNodes.push_back(newNode);
    }



    // Set all the children
    for(unsigned int i = 0; i < skeletonPoints->size(); i++){
        int cIdx = i;
        GetNeighbors(cIdx, &neighbors, skeletonPoints);
        SkeletonNode* currentNode = allNodes.at(i);

        if ( neighbors.size() > 2){
            currentNode->bifurcation = true;
        }
        if ( neighbors.size() < 2){
            currentNode->endPoint = true;
        }
        for(int j = 0; j < neighbors.size(); j++){
             int idNeighbor = neighbors.at(j);
             currentNode->children.push_back(allNodes.at(idNeighbor));
             allNodes.at(idNeighbor)->parents.push_back(currentNode);
        }
    }

    /*std::sort( dtValues.begin(), dtValues.end(),
           boost::bind(&std::pair<float, int>::first,_1) >
           boost::bind(&std::pair<float, int>::first,_2 ));

    std::cout << dtValues.at(0).first << "," << dtValues.at(0).second << std::endl;
    std::cout << dtValues.at(1).first << "," << dtValues.at(1).second << std::endl;*/


    //We need to apply a greedy approach...
    //first... insert the end points and bifurcations...

    //remove all the points that are dominated by them (w*dt)
    //select the highest point not eliminated ...

    set<int> removed;
    set<int> tobeTested;
    set<int> selected;
    vector<int> tmpSelected;

    for(unsigned  int i = 0; i < allNodes.size(); i++){
        if ( allNodes.at(i)->endPoint || allNodes.at(i)->bifurcation){

            // std::cout << "Init. " << allNodes.at(i)->initialDt << "/" <<  maxDt <<  std::endl;
             selected.insert(i);
        }
        else {
            tobeTested.insert(i);
        }
    }


   // std::cout << "Selected ? " << selected.size() << std::endl;
    /*
        Check the selected ones... they should be over one pixel from distance to each other...
        because currently they have a non-domination free pass
    */


    while(!tobeTested.empty()){
        //if there's still possibility of a point
        set<int> dominated;

        for(auto already: selected){
            //for each point that has already been selected, we remove those
            //that are under a certain value...
            SkeletonNode* currentNode = allNodes.at(already);
            for(auto possible: tobeTested){
                SkeletonNode* otherNode = allNodes.at(possible);
                auto diff = (currentNode->location - otherNode->location);
                double d = sqrt(diff[0]*diff[0] + diff[1]*diff[1]);
                double r1 = currentNode->initialDt*w;
                double r2 = otherNode->initialDt*w;

                if (d < (r1+r2)){
                    dominated.insert(possible);
                }
            }
        }

        for(auto dom: dominated){
            removed.insert(dom);
            tobeTested.erase( tobeTested.find(dom));
        }
        //now, find the one in the to be tested with the maximum dt
        int maxId = -1;
        float maxDtRest = 0;
        for(auto possible: tobeTested){
            if (dtValues.at(possible).first > maxDtRest ){
                maxDtRest = dtValues.at(possible).first;
                maxId = possible;
            }
        }

//        std::cout << "Max DT rest? " << maxDtRest << "/" << maxDtRest/maxDt <<  ".." << maxId << std::endl;
        if (maxId !=-1){
            selected.insert(maxId);
            tobeTested.erase( tobeTested.find(maxId));
        }
    }

    for(auto select: selected){
        allNodes.at(select)->idxInGraph = principalGraph->size();
        principalGraph->push_back(allNodes.at(select));
    }

    // Now that we have the nodes, we
    // connect the points for the neighbors...

    for(int i = 0; i < principalGraph->size(); i++){
        SkeletonNode* currentNode = principalGraph->at(i);
        //Navigate all parents to get the idx
        set<int> neighbors;
        set<int> visited;
        visited.insert(currentNode->label);

        //Navigate all children to get the idx
        //GetChildrenInGraph is DFS search..
        //until it reaches another node in the graph
        //if it's a node, add it to the neighbors list in the graph
        for(int k = 0; k < currentNode->children.size(); k++){
            GetChildrenInGraph(currentNode->children.at(k), principalGraph,
                               &neighbors,&visited);
        }
        //for each neighbor that exist in the graph, make a undirected graph
        // connections both ways..
         for(auto it = neighbors.begin(); it != neighbors.end(); ++it){
             currentNode->indexNeighborNode.push_back(*it);
         }

         if ( currentNode->indexNeighborNode.size() <= 1)
             currentNode->endPoint = true;
    }

}

void CreateSkeletonGraph(int controlPointSamplingSize,  vector<int>& endPoints,  vector<int>& bifurcations,  vector<Int2>* skeletonPoints,
                           vector<SkeletonNode*>* principalGraph){



    //The selection of the points is done via the sampling rate.
    //an automatic version is done via the dt transform
    // -- CreateSkeletonGraphAutomaticPoints
    //-------------------------------------------------------


    //std::cout << "Generating skeleton graph with " << controlPointSamplingSize << std::endl;
    vector<SkeletonNode*> allNodes;
    vector<int> neighbors;
    // Create a new node for each
    for(unsigned int i = 0; i < skeletonPoints->size(); i++){
        SkeletonNode* newNode = new SkeletonNode;
        newNode->label = i;
        newNode->location = skeletonPoints->at(i);
        bool endPoint = false;
        int idx = i;
        if ( count(endPoints.begin(), endPoints.end(),idx) > 0) endPoint = true;
        newNode->endPoint = endPoint;
        allNodes.push_back(newNode);
    }

    // Set all the children
    for(unsigned int i = 0; i < skeletonPoints->size(); i++){
        int cIdx = i;
        GetNeighbors(cIdx, &neighbors, skeletonPoints);
        SkeletonNode* currentNode = allNodes.at(i);

        if ( neighbors.size() > 2){
            currentNode->bifurcation = true;
        }
        for(int j = 0; j < neighbors.size(); j++){
             int idNeighbor = neighbors.at(j);
             currentNode->children.push_back(allNodes.at(idNeighbor));
             allNodes.at(idNeighbor)->parents.push_back(currentNode);
        }
    }
    //From here all skeleton points exist, with their neighbors....



    // Navigate from an endpoint, and define a count number

    queue<  pair<SkeletonNode*, int> > queryPoints;

    //selecting points....
    SkeletonNode* root;
    if ( endPoints.empty())
        root = allNodes.at(0);
    else
        root = allNodes.at( endPoints.at(0));
    int distance = 0;



    //queue of points looked at
    //start from the root.. i.e. sampling = 0
    queryPoints.push(make_pair(root,distance));
    bool visited[skeletonPoints->size()];
    bool added[skeletonPoints->size()];

    for(unsigned int i = 0; i < skeletonPoints->size(); i++){
        visited[i] = false;
        added[i] = false;
    }

    //add children if they haven't been seen, define the distance
    //in image space
    while(!queryPoints.empty()){
        auto topN = queryPoints.front();
        queryPoints.pop();
        SkeletonNode* top = topN.first;
        top->countInPath = topN.second;
        visited[top->label] = true;
        added[top->label] = true;

        for(int k = 0; k < top->children.size(); k++){
              SkeletonNode* child = top->children.at(k);

              if(!visited[child->label] && !added[child->label])
              {
                  added[child->label] = true;
                  queryPoints.push(make_pair(child, topN.second +1));
              }
        }
    }

    //if its a bifurcation, end point or the distance is in the control size,
    // add the point to the graph
    for(unsigned  int i = 0; i < allNodes.size(); i++){
        if ( allNodes.at(i)->countInPath % controlPointSamplingSize == 0 ||  allNodes.at(i)->endPoint || allNodes.at(i)->bifurcation){
            allNodes.at(i)->idxInGraph = principalGraph->size();
            principalGraph->push_back(allNodes.at(i));
        }
    }


    // Now that we have the nodes, we
    // connect the points for the neighbors...
    for(int i = 0; i < principalGraph->size(); i++){
        SkeletonNode* currentNode = principalGraph->at(i);
        //Navigate all parents to get the idx
        set<int> neighbors;
        set<int> visited;
        visited.insert(currentNode->label);

        //Navigate all children to get the idx
        //GetChildrenInGraph is DFS search..
        //until it reaches another node in the graph
        //if it's a node, add it to the neighbors list in the graph

        for(int k = 0; k < currentNode->children.size(); k++){
            GetChildrenInGraph(currentNode->children.at(k), principalGraph,
                               &neighbors,&visited);
        }

        //for each neighbor that exist in the graph, make a undirected graph
        // connections both ways..
         for(auto it = neighbors.begin(); it != neighbors.end(); ++it){

             currentNode->indexNeighborNode.push_back(*it);
             //currentNode->neighborsNodes.push_back(   principalGraph->at(*it)->label);
         }

    }
}

#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"

typedef itk::Image< float, 2 > HImageType;

void CreateITKImageCurve(int sz, bool DT, HImageType::Pointer img, std::vector<Int2>* points){
    HImageType::RegionType region;
    HImageType::IndexType index;
    index[0] = 0;      index[1] = 0;
    HImageType::SizeType sizeT;
    sizeT[0] = sz;
    sizeT[1] = sz;
    region.SetIndex(index);
    region.SetSize(sizeT);
    img->SetRegions(region);
    img->Allocate(true);
    img->FillBuffer(itk::NumericTraits< HImageType::PixelType >::Zero);
    img->Update();

    for(int i = 0; i < points->size(); i++){
        HImageType::IndexType idx;
        idx[0] = points->at(i)[0];
        idx[1] = points->at(i)[1];
        img->SetPixel(idx,itk::NumericTraits< HImageType::PixelType >::One );
    }
    if (DT){
        typedef itk::SignedMaurerDistanceMapImageFilter< HImageType, HImageType >  FilterType;
        typename FilterType::Pointer filter = FilterType::New();

        filter->SetInput( img );
        filter->SetBackgroundValue(itk::NumericTraits< HImageType::PixelType >::Zero);
        filter->SetSquaredDistance(false);
        filter->SetUseImageSpacing(false);
        filter->Update();

        HImageType::Pointer tmp = filter->GetOutput();

        itk::ImageRegionIterator< HImageType>  it(img, img->GetLargestPossibleRegion());
        itk::ImageRegionIterator< HImageType>  dit(tmp, tmp->GetLargestPossibleRegion());

        while(!it.IsAtEnd()){

            it.Set(dit.Get());
            ++it;
            ++dit;
        }

    }

}

double maxInS(HImageType::Pointer skeleton, HImageType::Pointer dt){
   //max_{x in S1}( DT2(x)
   //one simply walk over the points of one of the skeletons and measure the DT
   //of the other skeleton there, to find how close you are to the other skeleton

    double maxVal = 0;

    itk::ImageRegionIterator< HImageType>  it(skeleton, skeleton->GetLargestPossibleRegion());
    itk::ImageRegionIterator< HImageType>  dit(dt, dt->GetLargestPossibleRegion());

    while(!it.IsAtEnd()){

        if(it.Get() != itk::NumericTraits<HImageType::PixelType>::Zero){
           float val = dit.Get();
           if ( val >  maxVal){
               maxVal = val;
           }
        }
        ++it;
        ++dit;
    }

    return maxVal;
}

void RotateCurve(const std::vector<AugmentedPoint>& curve, std::vector<AugmentedPoint>& newCurve, int sz, float angle){
    // We assume the range is [0, sz]
    // so in order to perform a rotation, first we change it to
    // [-sz/2, sz/2 ] so the center of the image is the (0,0) coordinate
    // then apply the rotation matrix
    //
    // R(\theta) = [ cos(theta)   -sin(theta) ]
    //             [ sin(theta)    cos(theta) ]
    float radians = (angle * 3.14159)/180.0;
    for(int i = 0; i < curve.size(); i++){
        AugmentedPoint currentPoint = curve.at(i);


        AugmentedPoint newPoint;
        newPoint.amount = currentPoint.amount;
        newPoint.amountPerLength = currentPoint.amountPerLength;
        newPoint.orthoVariance = currentPoint.orthoVariance;

        // shift to create a center
        double x_prime = currentPoint.point[0] - sz/2.0;
        double y_prime = currentPoint.point[1] - sz/2.0;
        // apply rotation
        double newX = x_prime*cos(radians) - y_prime*cos(radians);
        double newY = x_prime*sin(radians) + y_prime*cos(radians);

        newPoint.point[0] = newX + sz/2.0;
        newPoint.point[1] = newY + sz/2.0;

         newCurve.push_back(newPoint);
    }
}



double WeightedHausdorffDistance(const std::vector<AugmentedPoint>& curve1, const std::vector<AugmentedPoint>& curve2, float densityWeight, float varianceWeight, int sz){

    // Either rotate or mirror the curve2

    // First we can align origins...
    float minLocCurve1[2] = {curve1.at(0).point[0], curve1.at(0).point[1]};
    float minLocCurve2[2] = {curve2.at(0).point[0], curve2.at(0).point[1]};

    for(int i = 0; i < curve1.size();i++){
        minLocCurve1[0] =  min(curve1.at(i).point[0],  minLocCurve1[0] );
        minLocCurve1[1] =  min(curve1.at(i).point[1],  minLocCurve1[1] );
    }
    //
    for(int i = 0; i < curve2.size();i++){
        minLocCurve2[0] =  min(curve2.at(i).point[0],  minLocCurve2[0] );
        minLocCurve2[1] =  min(curve2.at(i).point[1],  minLocCurve2[1] );
    }
    std::vector<AugmentedPoint> alignedToOriginC1;
    for(int i = 0; i< curve1.size(); i++){
        AugmentedPoint currentPoint = curve1.at(i);
        AugmentedPoint newPoint;
        newPoint.amount = currentPoint.amount;
        newPoint.amountPerLength = currentPoint.amountPerLength;
        newPoint.orthoVariance = currentPoint.orthoVariance;

        newPoint.point[0] = currentPoint.point[0] - minLocCurve1[0];
        newPoint.point[1] = currentPoint.point[1] - minLocCurve1[1];
        alignedToOriginC1.push_back(newPoint);
    }

    std::vector<AugmentedPoint> alignedToOriginC2;
    for(int i = 0; i< curve2.size(); i++){
        AugmentedPoint currentPoint = curve2.at(i);
        AugmentedPoint newPoint;
        newPoint.amount = currentPoint.amount;
        newPoint.amountPerLength = currentPoint.amountPerLength;
        newPoint.orthoVariance = currentPoint.orthoVariance;

        newPoint.point[0] = currentPoint.point[0] - minLocCurve2[0];
        newPoint.point[1] = currentPoint.point[1] - minLocCurve2[1];
        alignedToOriginC2.push_back(newPoint);
    }
    // So they are placed to the same corner ...

    float oDistance = WeightedHausdorffDistanceBase(alignedToOriginC1, alignedToOriginC2, densityWeight, varianceWeight, sz);
    std::cout << "Original distance " << oDistance;
    double minDistance= oDistance;

    int stepAngle = 90;
    for(int i =stepAngle; i < 360; i += stepAngle){ //


        std::vector<AugmentedPoint> curve2Rotated;
        RotateCurve(alignedToOriginC2, curve2Rotated, sz, stepAngle);
        double newVal = WeightedHausdorffDistanceBase(alignedToOriginC1, curve2Rotated, densityWeight, varianceWeight, sz);
        minDistance = std::min(newVal, minDistance);
    }
    std::cout << " min " << minDistance << std::endl;
    return minDistance;
}




double WeightedHausdorffDistanceBase(const std::vector<AugmentedPoint>& curve1, const std::vector<AugmentedPoint>& curve2, float densityWeight, float varianceWeight, int sz){

    std::cout << "Variance Weights " << varianceWeight << " , Density Weight " << densityWeight << std::endl;
    // Assume both curves are scaled to the width and height....
    int width = sz;
    int height = sz;


    Float2 way1 = curve1[curve1.size() - 1].point - curve1[0].point;
    Float2 way2 = curve2[curve2.size() - 1].point - curve2[0].point;
    float scale1 = sqrt(vectorLengthSquared(way1));
    float scale2 = sqrt(vectorLengthSquared(way2));

    // Create the DT from the first curve and the labels of where they are located
    Image2D<int> pointsImage1 { static_cast<size_t>(width), static_cast<size_t>(height) };
    Image2D<float> DT1 { static_cast<size_t>(width), static_cast<size_t>(height) };
    Image2D<char> dummyMask1 { static_cast<size_t>(width), static_cast<size_t>(height) };

    pointsImage1.fill(-1);
    DT1.fill(INFINITY);
    dummyMask1.fill(1);

    for(int i = 0; i <  curve1.size(); i++) {
        auto pt = curve1.at(i).point;
        pointsImage1(int(pt[0]), int(pt[1])) = i;
        DT1(pt[0], pt[1]) = 0;
    }

    auto origins1 = propagate2D(dummyMask1, DT1);
    auto labels1 = indexImage(pointsImage1, origins1);

    // Create the DT from the second curve and the labels of where they are located

    Image2D<int> pointsImage2 { static_cast<size_t>(width), static_cast<size_t>(height) };
    Image2D<float> DT2 { static_cast<size_t>(width), static_cast<size_t>(height) };
    Image2D<char> dummyMask2 { static_cast<size_t>(width), static_cast<size_t>(height) };

    pointsImage2.fill(-1);
    DT2.fill(INFINITY);
    dummyMask2.fill(1);

    for(int i = 0; i <  curve2.size(); i++) {
        auto pt = curve2.at(i).point;
        pointsImage2(int(pt[0]), int(pt[1])) = i;
        DT2(pt[0], pt[1]) = 0;
    }

    auto origins2 = propagate2D(dummyMask2, DT2);
    auto labels2 = indexImage(pointsImage2, origins2);


    // Now walk the first curve and get the DT from the second
    double maxS1 = 0;
    for(unsigned int i = 0; i < curve1.size(); i++){
        auto pt = curve1.at(i);
        Int2 loc{pt.point[0], pt.point[1]};
        float distance = DT2(loc[0],loc[1])/static_cast<float>(scale2);
        int idx2 = labels2(loc[0], loc[1]);
        auto secondPoint = curve2.at(idx2);

        auto v1 = pt.orthoVariance / scale1;
        auto v2 = secondPoint.orthoVariance / scale2;
        auto varianceDistanceSquared =  (v1 - v2) * (v1 - v2);

        auto a1 = pt.amountPerLength * scale1;
        auto a2 = secondPoint.amountPerLength * scale2;
        auto amountDistanceSquared = (a1 - a2) * (a1 - a2);
        auto value = densityWeight*amountDistanceSquared + varianceWeight*varianceDistanceSquared + (1.0 - densityWeight - varianceWeight)*distance;

        if ( value > maxS1)
            maxS1 = value;
    }
    // Now the other way around
    double maxS2 = 0;
    for(unsigned int i = 0; i < curve2.size(); i++){
        auto pt = curve2.at(i);
        Int2 loc{pt.point[0], pt.point[1]};
        float distance = DT1(loc[0],loc[1])/static_cast<float>(scale1);
        int idx1 = labels1(loc[0], loc[1]);
        auto secondPoint = curve1.at(idx1);

        auto v1 = pt.orthoVariance / scale2;
        auto v2 = secondPoint.orthoVariance / scale1;
        auto varianceDistanceSquared =  (v1 - v2) * (v1 - v2);

        auto a1 = pt.amountPerLength * scale2;
        auto a2 = secondPoint.amountPerLength * scale1;
        auto amountDistanceSquared = (a1 - a2) * (a1 - a2);
        auto value = densityWeight*amountDistanceSquared + varianceWeight*varianceDistanceSquared + (1.0 - densityWeight - varianceWeight)*distance;

        if ( value > maxS2)
            maxS2 = value;
    }

    double d = maxS1;
    if (maxS2 > d) d = maxS2;

    return d;
}

double HausdorffDistance(std::vector<Int2> pts1, std::vector<Int2> pts2, int sz){

    HImageType::Pointer dtimg1 = HImageType::New();
    HImageType::Pointer skimg1 = HImageType::New();
    CreateITKImageCurve(sz, true, dtimg1, &pts1);
    CreateITKImageCurve(sz, false, skimg1, &pts1);

    HImageType::Pointer dtimg2 = HImageType::New();
    HImageType::Pointer skimg2 = HImageType::New();
    CreateITKImageCurve(sz, true, dtimg2, &pts2);
    CreateITKImageCurve(sz, false, skimg2, &pts2);
    double d1 = maxInS(skimg1,dtimg2);
    double d2 = maxInS(skimg2,dtimg1);

    double d = d1;
    if (d2 > d) d = d2;
    return d;
}

double HausdorffDistance(PipelineResult* result, std::vector<AugmentedPoint>* generatingCurve, int sz){
    //
    std::vector<Int2> graphPoints;
    std::vector<Int2> generatingPoints;
    for(unsigned int i = 0; i < result->graphNodes.size(); i++)
        graphPoints.push_back(  result->graphNodes.at(i)->location);

    for(unsigned int i = 0; i < generatingCurve->size(); i++){
        Float2 pt = generatingCurve->at(i).point;
        Int2 newPt = {pt[0], pt[1]};
        generatingPoints.push_back(newPt);
    }
    return HausdorffDistance(graphPoints, generatingPoints, sz);
}

float graphLocalLength(int index, vector<SkeletonNode*>* nodes, bool useMaximum){
    float length = 0;
    int numberDescendants = 0;
    int numberAntecesors = 0;

    vector<int> antecesors;
    vector<int> descendants;

    SkeletonNode* current = nodes->at(index);
    for(int i =0; i  < current->indexNeighborNode.size(); i++){

        SkeletonNode* other = nodes->at( current->indexNeighborNode.at(i));

        if (current->countInPath <= other->countInPath  ){
            descendants.push_back( nodes->at(index)->indexNeighborNode.at(i));
            numberDescendants++;
        }
        else {
            antecesors.push_back( nodes->at(index)->indexNeighborNode.at(i) );
            numberAntecesors++;
        }
    }

    Float2 result;
    Int2 currentLoc = nodes->at(index)->location;

    float currentMaximum = 0;

    // For the graph we use the average (?)
    if (numberAntecesors == 0){
        // Start

        for(int j = 0; j < numberDescendants; j++){
            Float2 a =   nodes->at(descendants.at(j))->location - currentLoc;
            //a /= ;
            float v = sqrt(vectorLengthSquared(a));
            if (v > currentMaximum)
                currentMaximum = v;
            result += a;
        }
        result /= numberDescendants;
        //
        length = sqrt(vectorLengthSquared(result)); //
    }
    else if ( numberDescendants == 0){
        // Go back
        for(int j = 0; j < numberAntecesors; j++){
            Float2 a = currentLoc - nodes->at(antecesors.at(j))->location ;
            //a /= sqrt(vectorLengthSquared(a));
            float v = sqrt(vectorLengthSquared(a));
            if (v > currentMaximum)
                currentMaximum = v;
            result += a;
        }
        result /= numberAntecesors;
        length = sqrt(vectorLengthSquared(result)); //
    }
    else { // has both antecesors and descendants

        Float2 avgForward;
        for(int j = 0; j < numberDescendants; j++){
            Float2 a =   nodes->at(descendants.at(j))->location - currentLoc;
            //a /= sqrt(vectorLengthSquared(a));
            float v = sqrt(vectorLengthSquared(a));
            if (v > currentMaximum)
                currentMaximum = v;
            avgForward += a;
        }
        avgForward /= numberDescendants;

        Float2 avgBackward;
        for(int j = 0; j < numberAntecesors; j++){
            Float2 a = currentLoc - nodes->at(antecesors.at(j))->location ;
            //a /= sqrt(vectorLengthSquared(a));
            float v = sqrt(vectorLengthSquared(a));
            if (v > currentMaximum)
                currentMaximum = v;
            avgBackward += a;
        }
        avgBackward /= numberAntecesors;
        result = avgForward + avgBackward;
        length = sqrt(vectorLengthSquared(result)); //
    }

    if ( useMaximum)
        return currentMaximum;
    return length;
}

Float2 DirectionGraph(int index, vector<SkeletonNode*>* nodes) {
    int numberDescendants = 0;
    int numberAntecesors = 0;

    vector<int> antecesors;
    vector<int> descendants;
    SkeletonNode* current = nodes->at(index);
    for(int i =0; i  < current->indexNeighborNode.size(); i++){

        SkeletonNode* other = nodes->at( current->indexNeighborNode.at(i));

        if (current->countInPath <= other->countInPath  ){
            descendants.push_back( nodes->at(index)->indexNeighborNode.at(i));
            numberDescendants++;
        }
        else {
            antecesors.push_back( nodes->at(index)->indexNeighborNode.at(i) );
            numberAntecesors++;
        }
    }


    Float2 result;
    //Float2 currentLoc = nodes->at(index)->location;
    Float2 currentLoc = nodes->at(index)->stat2.pointStat.average;

    if (numberAntecesors == 0){
        // Start
        for(int j = 0; j < numberDescendants; j++){
            Float2 a =   nodes->at(descendants.at(j))->stat2.pointStat.average - currentLoc;
            a /= sqrt(vectorLengthSquared(a));
            result += a;
        }

        result /= numberDescendants;
    }
    else if ( numberDescendants == 0){
        // Go back
        for(int j = 0; j < numberAntecesors; j++){
            Float2 a = currentLoc - nodes->at(antecesors.at(j))->stat2.pointStat.average ;
            a /= sqrt(vectorLengthSquared(a));
            result += a;
        }
        result /= numberAntecesors;
    }
    else { // has both antecesors and descendants

        Float2 avgForward;
        for(int j = 0; j < numberDescendants; j++){
            Float2 a =   nodes->at(descendants.at(j))->stat2.pointStat.average - currentLoc;
            a /= sqrt(vectorLengthSquared(a));
            avgForward += a;
        }
        avgForward /= numberDescendants;

        Float2 avgBackward;
        for(int j = 0; j < numberAntecesors; j++){
            Float2 a = currentLoc - nodes->at(antecesors.at(j))->stat2.pointStat.average ;
            a /= sqrt(vectorLengthSquared(a));
            avgBackward += a;
        }
        avgBackward /= numberAntecesors;

        result = avgForward + avgBackward;
    }

    double size= sqrt(vectorLengthSquared(result));
    result /= size;

    if (isnan(result[0])|| isnan(result[1])){
       /* std::cout << std::endl;
        std::cout <<"Size of vector?`" << size <<" .. " << numberAntecesors<<","<< numberDescendants<<std::endl;
        std::cout <<"Neighbors? " <<  current->indexNeighborNode.size() <<std::endl;*/
    }

    if (current->indexNeighborNode.size() == 0){
        result[0] = 0;
        result[1] = 0;
    }
    return result;
}

void augmentGraph(vector<SkeletonNode *> *nodes, const Image2D<unsigned> &plot, int smoothCount){

    auto width = plot.width();
    auto height = plot.height();

    unsigned total = 0;
    for(int x = 0; x < width; x++) {
        for(int y = 0; y < height; y++) {
            total += plot(x, y);
        }
    }

    Image2D<int> pointsImage { width, height };
    Image2D<float> pointsDistance { width, height };
    Image2D<char> dummyMask { width, height };

    pointsImage.fill(-1);
    pointsDistance.fill(INFINITY);
    dummyMask.fill(1);

    for(int i = 0; i < nodes->size(); i++) {
        auto pt = nodes->at(i)->location;
        pointsImage(pt[0], pt[1]) = i;
        pointsDistance(pt[0], pt[1]) = 0;
    }

    auto origins = propagate2D(dummyMask, pointsDistance);
    auto labels = indexImage(pointsImage, origins);
    auto regions = regionStats(plot, labels, nodes->size());

    for(int i = 0; i < nodes->size(); i++) {
        auto direction = DirectionGraph(i, nodes);
        //auto direction = curveDirectionVector(curve, i);
        auto ortho = Float2({direction[1], -direction[0]});

        int levels = smoothCount;
        int startLevel = 0;
        set<int> neighboringRegions;
        neighboringRegions.insert(i);
        while( startLevel <= levels){
            vector<int> newRegions;
            for(auto it = neighboringRegions.begin(); it != neighboringRegions.end(); it++){
                // Get the neighboring regions in the first level and add them to the regions
                SkeletonNode* node = nodes->at(*it);
                for(int h = 0; h < node->indexNeighborNode.size(); h++)
                    newRegions.push_back( node->indexNeighborNode.at(h) );
            }

            // Now that all the close regions are added... we add them to the set
            for(int h = 0; h < newRegions.size(); h++)
                neighboringRegions.insert( newRegions.at(h));
            startLevel++;
        }

        // Location Itsel >  ap.point = Float2({ curve[i][0] * 1.0f, curve[i][1] * 1.0f });


        RegionStat smoothedStat;
        smoothedStat.pointStat.count = 0;
        smoothedStat.productStat.count = 0;
        for(auto it = neighboringRegions.begin(); it != neighboringRegions.end(); it++) {
            int j = *it;
            smoothedStat = smoothedStat + regions[j];
        }

        SkeletonNode* node = nodes->at(i);
        node->amount = regions[i].pointStat.count * 1.0f / total;  // ap.amount = regions[i].pointStat.count * 1.0f / total;

        float graphLength = graphLocalLength(i, nodes, true);
        node->amountPerLength = node->amount / graphLength;
        //  ap.amountPerLength = ap.amount / curveLocalLength(curve, i);

        if(!smoothedStat.empty()) {
             node->orthoVariance = smoothedStat.covarianceMatrix().varianceAlong(ortho); //            ap.orthoVariance = smoothedStat.covarianceMatrix().varianceAlong(ortho);
        }
        node->direction = direction; // ap.direction = direction;
    }
}


void MergeNodesInGraph(vector<SkeletonNode *> *nodes, vector<int> LabelstoJoin){

     int firstLabel = LabelstoJoin.at(0);
     set<int> redirect;


     vector<SkeletonNode*> tmpNodes;

     int indexOfFirstLabel = -1;
     for(int k =0; k < nodes->size(); k++){
         int currentLabel = nodes->at(k)->label;
         //std::cout <<  "(" <<   nodes->at(k)->location[0] << "," <<  nodes->at(k)->location[1] << ")  ";
         if ( currentLabel != firstLabel ){
            // Either merge or redirect ...
             if ( count(LabelstoJoin.begin(), LabelstoJoin.end(), currentLabel ) == 0)
                 tmpNodes.push_back( nodes->at(k));
         }
         else if ( currentLabel == firstLabel){
               //desired one
             indexOfFirstLabel = tmpNodes.size();
             tmpNodes.push_back( nodes->at(k));
         }
     }
     //std::cout << std::endl;
     //tmpNodes have the new indices ....




     for(int k =0; k < nodes->size(); k++){
         int currentLabel = nodes->at(k)->label;



         std::vector<int> neighborsLabels;
         for(int h = 0; h < nodes->at(k)->indexNeighborNode.size(); h++){
             neighborsLabels.push_back( nodes->at(nodes->at(k)->indexNeighborNode.at(h))->label);
         }

       if ( currentLabel ==  firstLabel){ // add all the neighbors of the ones removed
             for(int j = 1; j < LabelstoJoin.size(); j++){
                 for(int i = 0; i < nodes->size(); i++){
                     if ( count(LabelstoJoin.begin(), LabelstoJoin.end(), nodes->at(i)->label) > 0){
                         SkeletonNode* currentNode = nodes->at(i);
                         for(int h = 0; h < currentNode->indexNeighborNode.size(); h++){
                             neighborsLabels.push_back( nodes->at(currentNode->indexNeighborNode.at(h))->label);
                         }
                     }
                 }
             }
         }
         //------
         //std::cout << "A " <<  nodes->at(k)->indexNeighborNode.size() << std::endl;
         nodes->at(k)->indexNeighborNode.clear();
         for(int i = 0; i < neighborsLabels.size(); i++){
             int neighborLabel = neighborsLabels.at(i);
             //find in tmp nodes the index
             int newIndex = -1;
             for(int h = 0; h < tmpNodes.size(); h++){
                 if (tmpNodes.at(h)->label == neighborLabel){
                     newIndex = h;
                     break;
                 }
             }


             if ( newIndex != -1){

                 //std::cout << "enters " << std::endl;
                 if ( count(LabelstoJoin.begin(), LabelstoJoin.end(), neighborLabel ) > 0)
                 {
                     if ( neighborLabel == firstLabel){

                         if ( currentLabel != firstLabel){
                             if ( count(  nodes->at(k)->indexNeighborNode.begin(),
                                          nodes->at(k)->indexNeighborNode.end(), newIndex) == 0)
                                 nodes->at(k)->indexNeighborNode.push_back(newIndex);

                         }
                     }
                     else {
                           //need to check here ...
                         if (  count(  nodes->at(k)->indexNeighborNode.begin(),
                                       nodes->at(k)->indexNeighborNode.end(), indexOfFirstLabel) == 0 && currentLabel != firstLabel)
                             nodes->at(k)->indexNeighborNode.push_back(indexOfFirstLabel);
                  }
                 }
                 else {

                    if ( count(  nodes->at(k)->indexNeighborNode.begin(),
                                 nodes->at(k)->indexNeighborNode.end(), newIndex) == 0)
                             nodes->at(k)->indexNeighborNode.push_back(newIndex);
                 }
             }
             else {
                 //it was removed?...
                 if ( count(LabelstoJoin.begin(), LabelstoJoin.end(), neighborLabel ) > 0 && currentLabel != firstLabel)
                 {
                     //std::cout << "Removed but still a neighbor label? " << std::endl;
                     if (  count(  nodes->at(k)->indexNeighborNode.begin(),
                                   nodes->at(k)->indexNeighborNode.end(), indexOfFirstLabel) == 0){
                         nodes->at(k)->indexNeighborNode.push_back(indexOfFirstLabel);
                     }
                 }

             }
         }
         //std::cout << "B " <<  nodes->at(k)->indexNeighborNode.size() << std::endl;



         if (nodes->at(k)->indexNeighborNode.size() <= 1)
             nodes->at(k)->endPoint = true;
         else nodes->at(k)->endPoint = false;

         //Need to change, index Neighbornode, children and parents...
         std::vector<SkeletonNode*> newChildren;
         for(int h = 0; h < nodes->at(k)->children.size(); h++){
             int childLabel = nodes->at(k)->children.at(h)->label;
             if ( count(LabelstoJoin.begin(), LabelstoJoin.end(), childLabel ) > 0)
             {
           //if it's in the places to merge..
                 if ( childLabel == firstLabel && currentLabel != firstLabel)
                     newChildren.push_back( nodes->at(k)->children.at(h));
             }
             else {
                 newChildren.push_back( nodes->at(k)->children.at(h));
             }
         }
         nodes->at(k)->children.clear();
         for(int h = 0; h < newChildren.size(); h++){
              nodes->at(k)->children.push_back( newChildren.at(h));
         }


         std::vector<SkeletonNode*> newParents;
         for(int h = 0; h < nodes->at(k)->parents.size(); h++){
             int parentLabel = nodes->at(k)->parents.at(h)->label;
             if ( count(LabelstoJoin.begin(), LabelstoJoin.end(), parentLabel ) > 0)
             {
                 //if it's in the places to merge..
                 if (parentLabel == firstLabel && currentLabel != firstLabel)
                   newParents.push_back( nodes->at(k)->parents.at(h));
             }
             else {
                   newParents.push_back( nodes->at(k)->parents.at(h));
             }
         }

     }
     //****
     nodes->clear();
     for(int i = 0; i < tmpNodes.size(); i++){

         if (tmpNodes.at(i)->indexNeighborNode.size() < 2)
             tmpNodes.at(i)->endPoint = true;
         else
              tmpNodes.at(i)->endPoint = false;
         nodes->push_back(tmpNodes.at(i));
     }
     // std::cout << "End... "<< std::endl;
    // for(int k =0; k < nodes->size(); k++){
         //int currentLabel = nodes->at(k)->label;
      //   std::cout <<  "(" <<   nodes->at(k)->location[0] << "," <<  nodes->at(k)->location[1] << ")  ";

     //}

    // std::cout << std::endl;

     //std::cout << "************************************" << std::endl;


}

int principalGraph(vector<SkeletonNode *> *nodes, const Image2D<unsigned>& plot, int iterationCount, float convergenceDistance, const float arcSegment, int smoothCount, bool debug ){

    auto width = plot.width();
    auto height = plot.height();

    Image2D<int> pointsImage { width, height };
    Image2D<float> pointsDistance { width, height };
    Image2D<char> dummyMask { width, height };
    dummyMask.fill(1);


    vector<Int2> points;
    int nodesSize = nodes->size();

    // std::cout << "INITIAL NODES SIZE " <<nodesSize <<std::endl;
    for(int i = 0; i < nodesSize; i++) {
        auto pt = nodes->at(i)->location;
        points.push_back(pt);
    }
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
        for(int i = 0; i <  nodesSize; i++) {
            auto pt = nodes->at(i)->location;
            pointsImage(pt[0], pt[1]) = i;
            pointsDistance(pt[0], pt[1]) = 0;
        }

        //Each one starts first with its own origin... propagate
        //the centerline locations so that it saves where is the
        // origin...
        auto origins = propagate2D(dummyMask, pointsDistance);
        // before it was the coordinates.. now labels has
        // for each location the label it belongs to..
        auto labels = indexImage(pointsImage, origins);
        if (debug)
              DrawVoronoiGraph(iteration,labels,points, plot, nodes);


        // Plot is the image(density), labels are where each point belongs to, and the number of labels
        auto regions = regionStats(plot, labels, points.size());
        // That labels basically define the voronoi regions of each point...
        // Now for each region basically store the covariance ( needed for other steps )


        smoothedPoints.resize(0);
        smoothedPoints.reserve(regions.size());

        for(int i = 0; i < regions.size(); i++){
            Stat<Float2> smoothedStat;
            smoothedStat.count = 0;
            // The label is the location in the principal graph...
            int levels = smoothCount;
            int startLevel = 0;
            if ( nodes->at(i)-> endPoint)
                levels = 0;

            set<int> neighboringRegions;
            neighboringRegions.insert(i);
            while( startLevel <= levels){
                vector<int> newRegions;
                for(auto it = neighboringRegions.begin(); it != neighboringRegions.end(); it++){
                    // Get the neighboring regions in the first level and add them to the regions
                    SkeletonNode* node = nodes->at(*it);
                    for(int h = 0; h < node->indexNeighborNode.size(); h++)
                        newRegions.push_back( node->indexNeighborNode.at(h) );
                }

                // Now that all the close regions are added... we add them to the set
                for(int h = 0; h < newRegions.size(); h++)
                    neighboringRegions.insert( newRegions.at(h));
                startLevel++;
            }

            for(auto it = neighboringRegions.begin(); it != neighboringRegions.end(); it++) {
                int j = *it;
                smoothedStat = smoothedStat +  regions[j].pointStat; //regions[std::max(0, std::min(j, (int) regions.size() - 1))].pointStat;
            }

            // The smooth count is how far along the neighborhood we look... levels ....
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

        //---- if we should re-sample or not is here...  the skeleton nodes should be merged...
        std::vector<Int2> tmpSmoothedPoints;
        for(int i = 0; i < smoothedPoints.size(); i++) {
            auto tmp = Int2({int(round(smoothedPoints[i][0])), int(round(smoothedPoints[i][1]))});
            tmpSmoothedPoints.push_back(tmp);
        }

        bool applyMerging = false;

        int groups[tmpSmoothedPoints.size()];
        for(int k = 0; k < tmpSmoothedPoints.size(); k++){
             groups[k] = k;
        }

        double maxDistanceThreshold = 1.5;

        if ( tmpSmoothedPoints.size() == 2){
            maxDistanceThreshold = 0.05*width;
        }
        for(int k = 0; k < tmpSmoothedPoints.size(); k++){
            auto pt1 = tmpSmoothedPoints.at(k);
            for(int h = k+1; h < tmpSmoothedPoints.size(); h++){
                auto pt2 = tmpSmoothedPoints.at(h);
                double d = sqrt((pt1[0] - pt2[0])*(pt1[0] - pt2[0]) + (pt1[1] - pt2[1])*(pt1[1] - pt2[1]));

                if ( d < maxDistanceThreshold){
                    groups[k] = std::min(groups[k],groups[h]);
                    groups[h] = std::min(groups[k],groups[h]);
                    applyMerging = true;
                }
            }
        }

        if ( applyMerging){
           // std::cout << "Merging? Iter:" << iteration << std::endl;
            std::set<int> diffGroups;
            std::map<int, int> currentIndices;


            for(int k = 0; k < tmpSmoothedPoints.size(); k++){
                     diffGroups.insert( groups[k]);
                     currentIndices[k] = nodes->at(k)->label;
            }

            vector<vector<int>> elementsInGroups;
            for(auto group: diffGroups){
               vector<int> elementsInGroup;
               for(int k = 0; k < tmpSmoothedPoints.size(); k++){
                   if (groups[k] == group) elementsInGroup.push_back(nodes->at(k)->label);
               }
               elementsInGroups.push_back(elementsInGroup);
            }



            for(auto elementsInGroup: elementsInGroups){
               if ( elementsInGroup.size() > 1)
                        MergeNodesInGraph(nodes, elementsInGroup);

            }
           nodesSize = nodes->size();


           points.clear();
            for(int i = 0; i < smoothedPoints.size(); i++) {
                //i is the original index...
                int labelOfOriginalIndex = currentIndices[i];
                int newIndex = -1;
                for(int h = 0;h <nodes->size(); h++){

                    if (nodes->at(h)->label == labelOfOriginalIndex){
                        newIndex = h;
                    }
                }

                if (newIndex != -1 ){
                    nodes->at(newIndex)->location = Int2({int(round(smoothedPoints[i][0])), int(round(smoothedPoints[i][1]))});
                    nodes->at(newIndex)->stat = regions.at(i);
                    points.push_back(nodes->at(newIndex)->location );

                }
            }
        }
        else {

            for(int i = 0; i < smoothedPoints.size(); i++) {
                points[i] = Int2({int(round(smoothedPoints[i][0])), int(round(smoothedPoints[i][1]))});
                nodes->at(i)->location = points[i];
                nodes->at(i)->stat = regions.at(i);
            }

        }



        //convergent lines should have a fixed amount of points
        //regardless of resampling

        if(points == pointsBefore || points == pointsBefore2) {
            convergence = true;
            //Save the invo of the regions

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

    //std::cout << "iteration " << iteration << "...." << iterationCount << std::endl;

   // std::cout << "END NODES SIZE " <<nodesSize <<std::endl;


   return convergence ? iteration : -1;
}




Image2D<float> GetBands(vector<SkeletonNode *> *nodes, int sz, bool resample){
    auto width = sz;
    auto height = sz;

    Image2D<int> pointsImage { static_cast<size_t>(width), static_cast<size_t>(height) };
    Image2D<float> pointsDistance { static_cast<size_t>(width), static_cast<size_t>(height) };
    Image2D<char> dummyMask { static_cast<size_t>(width), static_cast<size_t>(height) };
    dummyMask.fill(1);
    pointsImage.fill(-1);
    pointsDistance.fill(INFINITY);
    float curveLength = 0;


    vector<SkeletonNode*> resampled;
    vector< pair<int,int> > seenEdges;

    // Values that the nodes should have
    // amount, amountPerLength, orthoVriance

     if ( resample){

         for(int i = 0; i < nodes->size(); i++){
             //std::vector<int> indexNeighborNode;
             SkeletonNode *cNode = nodes->at(i);
             resampled.push_back(cNode);

             for(int j = 0; j < cNode->indexNeighborNode.size(); j++){
                 int o = cNode->indexNeighborNode.at(j);
                 auto key1 = make_pair(i, o);
                 auto key2 = make_pair(o, i);
                 bool seen = false;

                seen =( count(seenEdges.begin(), seenEdges.end(), key1) > 0 ||  count(seenEdges.begin(), seenEdges.end(), key2) > 0);

                if (!seen){
                    // Create a new node by interpolating
                    seenEdges.push_back(key1);
                    seenEdges.push_back(key2);
                    SkeletonNode *oNode = nodes->at(o);


                    SkeletonNode* newNode = new SkeletonNode();
                    newNode->amount = (cNode->amount + oNode->amount)/2.0;
                    newNode->amountPerLength = ( cNode->amountPerLength + oNode->amountPerLength)/2.0;
                    newNode->location = (cNode->location + oNode->location)/2.0;
                    newNode->orthoVariance = (cNode->orthoVariance + oNode->orthoVariance)/2.0;

                    resampled.push_back(newNode);

                }

             }

         }


     }
     else {
         for(int i = 0; i < nodes->size(); i++) {
             resampled.push_back(nodes->at(i));
         }
     }



    for(int i = 0; i < resampled.size(); i++) {
        auto pt = resampled.at(i)->location;

        if ( resampled.at(i)->amount > 0)
            curveLength +=  (resampled.at(i)->amount / resampled.at(i)->amountPerLength);
        pointsImage(pt[0], pt[1]) = i;
        pointsDistance(pt[0], pt[1]) = 0;
    }

    auto origins = propagate2D(dummyMask, pointsDistance);
    auto labels = indexImage(pointsImage, origins);

    Image2D<float> test{static_cast<size_t>(width), static_cast<size_t>(height)};

    const float sigmaFactor = 1.2;
    for(unsigned int x = 0; x < test.width(); x++){
          for(unsigned int y = 0; y < test.height(); y++){

               int idx = labels(x,y);
               float distance = pointsDistance(x,y);

               float sigma = sqrt(resampled.at(idx)->orthoVariance);
               float diff = sigma*sigmaFactor - distance;
               // if distance is higher than variance diff will be negative
               const float lr = 0.92, lf = 50;
               float amountPerLength = resampled.at(idx)->amountPerLength;

               float l1 = std::atan(lf * (amountPerLength - 1 / curveLength)) / (0.5 * M_PI) + 0.5;
               l1 = 1 - lr + lr * l1;

               if ( resampled.at(idx)->endPoint)
                   diff = -1;

               if ( diff < 0){
                   test(x,y) = 0;
               }
               else{
                   test(x,y) = l1;
               }

          }
    }
    return test;
}
