#include "wilkinsonscagnostics.h"

WilkinsonScagnostics::WilkinsonScagnostics()
{
   wrapper = NULL;
}

void WilkinsonScagnostics::SetSplomGenerator(SPLOMWrapper *_wrapper){
    wrapper = _wrapper;
}

WilkinsonScagnostics::~WilkinsonScagnostics(){

}

vtkSmartPointer<vtkTree> WilkinsonScagnostics::CreateMST(int attr1, int attr2, int clusteringAttribute, int valueToFilter,vector<float> *edgeDistances,
                                                         vector<LocalPoint> *points, bool withBinning){

    vtkSmartPointer<vtkTree> minimumSpanningTree = vtkSmartPointer<vtkTree>::New();
    if ( wrapper == NULL) return minimumSpanningTree;

    if ( withBinning){
        return CreateMSTWithBinning(attr1, attr2, clusteringAttribute, valueToFilter, edgeDistances, points);
    }
    else {
        return CreateMSTWithPoints(attr1, attr2, clusteringAttribute, valueToFilter, edgeDistances, points);
    }
}

double WilkinsonScagnostics::GetRVM(int attr1, int attr2, int clusteringAttribute, int valueToFilter){
  // Andrada Tatu et al. measure as defined in Combining automated analysis and visualization techniques for effective
  // exploration of high dimensional data
  if (wrapper->IsRunning()) return -1;
  int size = wrapper->GetSizeTexture();
  short singlePlane[size*size*3];

  wrapper->GetImage(attr1, attr2, clusteringAttribute, valueToFilter, singlePlane);
  vector<LocalPoint> points;
  // We need to generate a fully connected
  for(int k = 0; k < size*size; k++){
      if ( singlePlane[k*3 + 0] != 0 && singlePlane[k*3 + 0] == singlePlane[k*3 +1] && singlePlane[k*3 + 0] == singlePlane[k*3 +2] ){

          LocalPoint p;
          int row = k / size;
          int col = k % size;
          p.p[0] = static_cast<float>(col)/ static_cast<float>(size);
          p.p[1] = static_cast<float>(row)/ static_cast<float>(size);
          points.push_back(p);
      }
  }
  // got the drawn points

  float* cudaPoints = (float*) malloc(sizeof(float)*points.size()*2);

  for(unsigned int i = 0; i < points.size(); i++){
      cudaPoints[i*2 + 0] = points.at(i).p[0];
      cudaPoints[i*2 + 1] = points.at(i).p[1];
  }
  double val = RotatingVarianceMeasure(size, cudaPoints, points.size(), 10);

  free(cudaPoints);
  return val;
}

vtkSmartPointer<vtkTree> WilkinsonScagnostics::CreateMSTWithBinning(int attr1, int attr2, int clusteringAttribute, int valueToFilter,
                                                                    vector<float>* edgeDistances, vector<LocalPoint>* binPoints){
    vtkSmartPointer<vtkTree> minimumSpanningTree = vtkSmartPointer<vtkTree>::New();
    if (wrapper->IsRunning()) return minimumSpanningTree;
    edgeDistances->clear();

    int size = wrapper->GetSizeTexture();
    short singlePlane[size*size*3];

    wrapper->GetImage(attr1, attr2, clusteringAttribute, valueToFilter, singlePlane);

    // We need to generate a fully connected
    for(int k = 0; k < size*size; k++){
        if ( singlePlane[k*3 + 0] != 0 && singlePlane[k*3 + 0] == singlePlane[k*3 +1] && singlePlane[k*3 + 0] == singlePlane[k*3 +2] ){

            LocalPoint p;
            int row = k / size;
            int col = k % size;
            p.p[0] = static_cast<float>(col)/ static_cast<float>(size);
            p.p[1] = static_cast<float>(row)/ static_cast<float>(size);
            binPoints->push_back(p);
        }
    }
    //
    vtkSmartPointer<vtkMutableUndirectedGraph> g =
        vtkSmartPointer<vtkMutableUndirectedGraph>::New();

    // added the amount of points as vertices
    vtkIdType ids[binPoints->size()];
    for(unsigned int i = 0; i < binPoints->size(); i++){
        ids[i] = g->AddVertex();
    }
    // Create a fully connected graph
    vtkSmartPointer<vtkDoubleArray> weights = vtkSmartPointer<vtkDoubleArray>::New();
    weights->SetNumberOfComponents(1);
    weights->SetName("Weights");
    //
    for(unsigned int i = 0; i < binPoints->size(); i++){
        // since it is undirected, we don't need to go from 0
        LocalPoint pi = binPoints->at(i);

        for(unsigned int j = i+1; j < binPoints->size(); j++)
        {
               g->AddEdge( ids[i], ids[j]);
               LocalPoint pj = binPoints->at(j);

               double distance = sqrt( (pi.p[0] - pj.p[0])*(pi.p[0] - pj.p[0])  + (pi.p[1] - pj.p[1])*(pi.p[1] - pj.p[1]) );
               weights->InsertNextValue(distance);
        }
    }
    // The weights are the distance between the two points
     g->GetEdgeData()->AddArray(weights);
     // Setup the minimum spanning tree filter
  vtkSmartPointer<vtkBoostPrimMinimumSpanningTree> minimumSpanningTreeFilter = vtkSmartPointer<vtkBoostPrimMinimumSpanningTree>::New();
  minimumSpanningTreeFilter->SetOriginVertex(ids[0]);
  minimumSpanningTreeFilter->SetInputData(g);
  minimumSpanningTreeFilter->SetEdgeWeightArrayName("Weights");

  // Compute the minimum spanning tree
  minimumSpanningTreeFilter->CreateGraphVertexIdArrayOn();
  minimumSpanningTreeFilter->Update();
  minimumSpanningTree->ShallowCopy(minimumSpanningTreeFilter->GetOutput());

  vtkSmartPointer<vtkEdgeListIterator> it = vtkSmartPointer<vtkEdgeListIterator>::New();
  minimumSpanningTree->GetEdges(it);

  while(it->HasNext()){
      vtkEdgeType e = it->Next();
      LocalPoint pi = binPoints->at(e.Source);
      LocalPoint pj = binPoints->at(e.Target);
      double distance = sqrt( (pi.p[0] - pj.p[0])*(pi.p[0] - pj.p[0])  + (pi.p[1] - pj.p[1])*(pi.p[1] - pj.p[1]) );
      edgeDistances->push_back(distance);
  }

  std::sort(edgeDistances->begin(), edgeDistances->end());
  return minimumSpanningTree;
}

double WilkinsonScagnostics::GetW(vector<float>* edgeDistances){

    float size = edgeDistances->size();
    float q75 = edgeDistances->at( size*0.75 );
    float q25 = edgeDistances->at( size*0.25 );

    float w = q75 + 1.5*(q75- q25);
    return w;
}

void WilkinsonScagnostics::GetOutliers(vtkSmartPointer<vtkTree> mst, vector<int>* idsOfOutliers, double w, vector<LocalPoint> *points){
   // We consider an outlier to be a vertex whose adjacent edges in the MST all have a weight greater than w
   vtkSmartPointer<vtkVertexListIterator> it = vtkSmartPointer<vtkVertexListIterator>::New();

   mst->GetVertices( it);

   vtkSmartPointer<vtkInEdgeIterator> inIt = vtkSmartPointer<vtkInEdgeIterator>::New();
   vtkSmartPointer<vtkOutEdgeIterator> outIt = vtkSmartPointer<vtkOutEdgeIterator>::New();

   while(it->HasNext()){
       vtkIdType id = it->Next();
       mst->GetInEdges(id, inIt);
       bool lowerThanW = false;
       LocalPoint pi = points->at(id);
       while(inIt->HasNext()){
          vtkInEdgeType inEdge = inIt->Next();
          LocalPoint pj = points->at(inEdge.Source);
          double distance = sqrt( (pi.p[0] - pj.p[0])*(pi.p[0] - pj.p[0])  + (pi.p[1] - pj.p[1])*(pi.p[1] - pj.p[1]) );

          if ( distance <= w) lowerThanW = true;
       }

       mst->GetOutEdges(id, outIt);
       while(outIt->HasNext()){
           vtkOutEdgeType outEdge = outIt->Next();
           LocalPoint pj = points->at(outEdge.Target);
           double distance = sqrt( (pi.p[0] - pj.p[0])*(pi.p[0] - pj.p[0])  + (pi.p[1] - pj.p[1])*(pi.p[1] - pj.p[1]) );
           if ( distance <= w) lowerThanW = true;
       }

       if (!lowerThanW){
           idsOfOutliers->push_back(id);
       }
   }
}

double WilkinsonScagnostics::SkewednessMeasure(vtkSmartPointer<vtkTree> mst,vector<int>* idsOfOutliers, vector<LocalPoint> *points ){
   double skewed = 0;

   vector<float> distances;
   vtkSmartPointer<vtkEdgeListIterator> it = vtkSmartPointer<vtkEdgeListIterator>::New();
   mst->GetEdges(it);

   while(it->HasNext()){
       vtkEdgeType e = it->Next();
       if (count(idsOfOutliers->begin(), idsOfOutliers->end(), e.Source ) > 0) continue;
       if (count(idsOfOutliers->begin(), idsOfOutliers->end(), e.Target ) > 0) continue;
       // removed the outliers
       LocalPoint pi = points->at(e.Source);
       LocalPoint pj = points->at(e.Target);
       double distance = sqrt( (pi.p[0] - pj.p[0])*(pi.p[0] - pj.p[0])  + (pi.p[1] - pj.p[1])*(pi.p[1] - pj.p[1]) );
       distances.push_back(distance);
   }

   std::sort(distances.begin(), distances.end());

   float size = distances.size();
   float q90  = distances.at( size*0.90 );
   float q50  = distances.at( size*0.50 );
   float q10  = distances.at( size*0.10 );

   skewed = (q90 - q50) / (q90 - q10);


   if ( fabs(fabs(q90-50) - fabs(q90-q10)) < 0.000001)
       skewed = 1;
   std::cout << "Skewed " << q90 << "? " << q50 << " , " << q10 <<  ":: " << skewed << std::endl;

   return skewed;
}

double WilkinsonScagnostics::ClumpyMeasure(vtkSmartPointer<vtkTree> mst,vector<int>* idsOfOutliers, vector<LocalPoint> *points ){
    vector<double> distances;
    vtkSmartPointer<vtkEdgeListIterator> it = vtkSmartPointer<vtkEdgeListIterator>::New();
    mst->GetEdges(it);

    while(it->HasNext()){
        vtkEdgeType e = it->Next();
        if (count(idsOfOutliers->begin(), idsOfOutliers->end(), e.Source ) > 0) continue;
        if (count(idsOfOutliers->begin(), idsOfOutliers->end(), e.Target ) > 0) continue;
        // removed the outliers
        LocalPoint pi = points->at(e.Source);
        LocalPoint pj = points->at(e.Target);
        double distance = sqrt( (pi.p[0] - pj.p[0])*(pi.p[0] - pj.p[0])  + (pi.p[1] - pj.p[1])*(pi.p[1] - pj.p[1]) );
        distances.push_back(distance);
    }

    std::sort(distances.begin(), distances.end());


    // The runt graph corresponding to each edge is the smaller of the two subsets
    // that are still connected to each of the two vertices in ej after deleting edges in the
    // mst with lenghts less than ej

    double c_clumpy = 0;
    for(unsigned int j = 0; j < distances.size(); j++){
        double length_ej = distances.at(j);

        float max_ek = -1;
        for(unsigned int k =0; k < distances.size(); k++){
            if (distances.at(k) < length_ej) continue; // we remove edges in the MST with leghts less than length(ej)
            if (distances.at(k) > max_ek) max_ek = distances.at(k);
        }

        float c_val = 1.0 - max_ek/length_ej;
         if ( c_val > c_clumpy) c_clumpy = c_val;
    }

    return c_clumpy;
}

double WilkinsonScagnostics::SparsenessMeasure(vtkSmartPointer<vtkTree> mst,vector<int>* idsOfOutliers, vector<LocalPoint> *points ){

    vector<float> distances;
    vtkSmartPointer<vtkEdgeListIterator> it = vtkSmartPointer<vtkEdgeListIterator>::New();
    mst->GetEdges(it);

    while(it->HasNext()){
        vtkEdgeType e = it->Next();
        if (count(idsOfOutliers->begin(), idsOfOutliers->end(), e.Source ) > 0) continue;
        if (count(idsOfOutliers->begin(), idsOfOutliers->end(), e.Target ) > 0) continue;
        // removed the outliers
        LocalPoint pi = points->at(e.Source);
        LocalPoint pj = points->at(e.Target);
        double distance = sqrt( (pi.p[0] - pj.p[0])*(pi.p[0] - pj.p[0])  + (pi.p[1] - pj.p[1])*(pi.p[1] - pj.p[1]) );
        distances.push_back(distance);
    }

    std::sort(distances.begin(), distances.end());

    float size = distances.size();
    float q90  = distances.at( size*0.90 );
    return q90;
}

double WilkinsonScagnostics::OutlyingMeasure(vtkSmartPointer<vtkTree> mst,vector<int>* idsOfOutliers, vector<LocalPoint> *points ){
    /*
     *   Gesamtlänge der langen Kanten im MST
     *   ------------------------------------
     *      Gesamtlänge aller Kanten im MST
    */
    float gesamtLangeLangenKanten = 0;
    float gesamtLangeAllerKanten = 0;

    vtkSmartPointer<vtkEdgeListIterator> it = vtkSmartPointer<vtkEdgeListIterator>::New();
    mst->GetEdges(it);

    while(it->HasNext()){
        vtkEdgeType e = it->Next();
        LocalPoint pi = points->at(e.Source);
        LocalPoint pj = points->at(e.Target);
        double distance = sqrt( (pi.p[0] - pj.p[0])*(pi.p[0] - pj.p[0])  + (pi.p[1] - pj.p[1])*(pi.p[1] - pj.p[1]) );

        if (count(idsOfOutliers->begin(), idsOfOutliers->end(), e.Source ) > 0){
            gesamtLangeLangenKanten += distance;
        }
        else  if (count(idsOfOutliers->begin(), idsOfOutliers->end(), e.Target ) > 0){
            gesamtLangeLangenKanten += distance;
        }

        gesamtLangeAllerKanten += distance;
    }

    double c_outlying = gesamtLangeLangenKanten / gesamtLangeAllerKanten;

    return c_outlying;
}

double WilkinsonScagnostics::StringyMeasure(vtkSmartPointer<vtkTree> mst,vector<int>* idsOfOutliers ){


    vtkSmartPointer<vtkVertexListIterator> it = vtkSmartPointer<vtkVertexListIterator>::New();
    mst->GetVertices( it);

    vtkSmartPointer<vtkInEdgeIterator> inIt = vtkSmartPointer<vtkInEdgeIterator>::New();
    vtkSmartPointer<vtkOutEdgeIterator> outIt = vtkSmartPointer<vtkOutEdgeIterator>::New();

    float numVertices= 0;
    float numSingleDegree = 0;
    float numDoubleDegree = 0;

    while(it->HasNext()){
        vtkIdType id = it->Next();
        if (count(idsOfOutliers->begin(), idsOfOutliers->end(), id) ) continue; //if its outlier then

        mst->GetInEdges(id, inIt);
        mst->GetOutEdges(id, outIt);

        int totalDegree = 0;

        while(inIt->HasNext()){
           vtkInEdgeType inEdge = inIt->Next();
           totalDegree++;
        }

        while(outIt->HasNext()){
           vtkOutEdgeType outEdge = outIt->Next();
           totalDegree++;
        }

        numVertices += 1;
        if ( totalDegree == 1) numSingleDegree += 1;
        if ( totalDegree == 2) numDoubleDegree += 1;
    }

    double c_stringy =  numDoubleDegree / (numVertices - numSingleDegree );

    return c_stringy;
}

double WilkinsonScagnostics::StriateMeasure(vtkSmartPointer<vtkTree> mst,vector<int>* idsOfOutliers, vector<LocalPoint> *points ){

    vtkSmartPointer<vtkVertexListIterator> it = vtkSmartPointer<vtkVertexListIterator>::New();
    mst->GetVertices( it);

    vtkSmartPointer<vtkInEdgeIterator> inIt = vtkSmartPointer<vtkInEdgeIterator>::New();
    vtkSmartPointer<vtkOutEdgeIterator> outIt = vtkSmartPointer<vtkOutEdgeIterator>::New();


    int totalVerts = 0;
    int totalStriate = 0;
    while(it->HasNext()){
        vtkIdType id = it->Next();
        if (count(idsOfOutliers->begin(), idsOfOutliers->end(), id) ) continue; //if its outlier then
        mst->GetInEdges(id, inIt);
        mst->GetOutEdges(id, outIt);

        LocalPoint p1 = points->at(id);
        vector<int> others;

        int sizeVEdges = 0;
        while( inIt->HasNext() ) {
            vtkInEdgeType inEdge = inIt->Next();
            others.push_back(inEdge.Source);
            sizeVEdges++;
        }


        while( outIt->HasNext() ) {
            vtkOutEdgeType outEdge = outIt->Next();
            others.push_back(outEdge.Target);
            sizeVEdges++;
        }

        if ( sizeVEdges == 2){
           // cos e(v,a) e(v,b)
            LocalPoint a = points->at(others.at(0));
            LocalPoint b = points->at(others.at(1));

            float v1[3] = {a.p[0] - p1.p[0], a.p[1] - p1.p[1], 0};
            float v2[3] = {b.p[0] - p1.p[0], b.p[1] - p1.p[1], 0};

            float c = MathHelper::CosineAngleOfVectors(v1,v2);

            if ( c <-0.75 ){
                totalStriate++;
            }
        }
        totalVerts++;
    }

    float c_striate = static_cast<float>(totalStriate) / static_cast<float>(totalVerts);
    return c_striate;
}

vtkSmartPointer<vtkPolyData> WilkinsonScagnostics::GetAlphaShape(float alphaValue,vector<int>* idsOfOutliers, vector<LocalPoint> *points){
  // We use the sparness measure to generate the alpha-shape as SparsenessMeasure is q90

    vtkSmartPointer<vtkPoints> vtpoints =
      vtkSmartPointer<vtkPoints>::New();

     for(unsigned int i = 0; i < points->size(); i++){
         if (count(idsOfOutliers->begin(), idsOfOutliers->end(), i) ) continue;

         vtpoints->InsertNextPoint( points->at(i).p[0], points->at(i).p[1],0);
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

vtkSmartPointer<vtkPolyData> WilkinsonScagnostics::GetConvexHull(vector<int> *idsOfOutliers, vector<LocalPoint> *points){
    vtkSmartPointer<vtkPoints> vtpoints =
      vtkSmartPointer<vtkPoints>::New();

     for(unsigned int i = 0; i < points->size(); i++){

          if (count(idsOfOutliers->begin(), idsOfOutliers->end(), i) ) continue;

         vtpoints->InsertNextPoint( points->at(i).p[0], points->at(i).p[1],0);
     }

    // Add the grid points to a polydata object
    vtkSmartPointer<vtkPolyData> polydata =
      vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(vtpoints);

    // Triangulate the grid points
    vtkSmartPointer<vtkDelaunay2D> delaunay =
    vtkSmartPointer<vtkDelaunay2D>::New();
    delaunay->SetInputData(polydata);
    delaunay->Update();

    vtkSmartPointer<vtkPolyData> AlphaShape = delaunay->GetOutput();
    return AlphaShape;
}

vtkSmartPointer<vtkPolyData> WilkinsonScagnostics::GetConvexHull(vector<LocalPoint> *points){
    vtkSmartPointer<vtkPoints> vtpoints =
      vtkSmartPointer<vtkPoints>::New();

     for(unsigned int i = 0; i < points->size(); i++){
         vtpoints->InsertNextPoint( points->at(i).p[0], points->at(i).p[1],0);
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




double WilkinsonScagnostics::ConvexMeasure(vtkSmartPointer<vtkPolyData> alphaShape, vtkSmartPointer<vtkPolyData> convexHull){

    double areaA = GetArea(alphaShape);
    double areaH = GetArea(convexHull);
    double convexity = areaA/areaH;

    return convexity;
}
double WilkinsonScagnostics::SkinnyMeasure(vtkSmartPointer<vtkPolyData> alphaShape){
   // We need area and perimeter
   double totalArea = GetArea(alphaShape);
   double totalPerimeter = GetPerimeter(alphaShape);

   std::cout << "Area? " <<  totalArea <<  ", " << totalPerimeter << std::endl;
   double skinny = 1.0 - sqrt(4*3.14159*totalArea)/ totalPerimeter;

   return skinny;
}

double WilkinsonScagnostics::GetArea(vtkSmartPointer<vtkPolyData> polydata){
    vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
    cellArray = polydata->GetPolys();

    double totalArea = 0;
    int totalTriangles = 0;
    for(int i = 0; i < polydata->GetNumberOfCells(); i++){

        vtkCell* cell = polydata->GetCell(i);

        if ( cell->GetNumberOfPoints() == 3 ){
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
    }

    return totalArea;
}


double WilkinsonScagnostics::GetPerimeter(vtkSmartPointer<vtkPolyData> polydata){
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

            totalPerim += MathHelper::DistanceBtwPoints(p0,p1);
        }

  }
   return totalPerim;
}




double WilkinsonScagnostics::MonotonicMeasure(vector<int> *idsOfOutliers, vector<LocalPoint> *points){

    // first we need the mean of x and y
    double meanX = 0, meanY = 0;
    double total = 0;
    for(unsigned int i = 0; i < points->size(); i++){
       if (count(idsOfOutliers->begin(), idsOfOutliers->end(), i) ) continue; // if 0 then false,
       meanX += points->at(i).p[0];
       meanY += points->at(i).p[1];
       total +=1;
    }
    meanX /= total;
    meanY /= total;


    double up = 0;
    double den1 = 0;
    double den2 = 0;

    for(unsigned int i = 0; i < points->size(); i++){
         if (count(idsOfOutliers->begin(), idsOfOutliers->end(), i) ) continue;
         double xi = points->at(i).p[0];
         double yi = points->at(i).p[1];

         up += (xi - meanX)*(yi - meanY);
         den1 += (xi - meanX)*(xi- meanX);
         den2 += (yi - meanY)*(yi- meanY);
    }

    double den = sqrt(den1)*sqrt(den2);

    double r = up/ den;
    double r2 = r*r;

    return r2;
}



vtkSmartPointer<vtkTree> WilkinsonScagnostics::CreateMSTWithPoints(int attr1, int attr2, int clusteringAttribute, int valueToFilter, vector<float> *edgeDistances, vector<LocalPoint> *points){
    vtkSmartPointer<vtkTree> minimumSpanningTree = vtkSmartPointer<vtkTree>::New();

    Statistics* _stats = wrapper->GetStats();

    float minAttr1 = _stats->GetMinimumInAttribute(attr1), maxAttr1 = _stats->GetMaximumInAttribute(attr1);
    float range1 = maxAttr1 - minAttr1;

    float minAttr2 = _stats->GetMinimumInAttribute(attr2), maxAttr2 = _stats->GetMaximumInAttribute(attr2);
    float range2 = maxAttr2 - minAttr2;

    int totalElements = wrapper->GetData()->GetTotalNumberOfElements();

    for(int k = 0; k < totalElements; k++){

        float v1 = wrapper->GetData()->GetElementValue(k, attr1);
        float v2 = wrapper->GetData()->GetElementValue(k, attr2);
        if ( clusteringAttribute != -1){

            int value = wrapper->GetData()->GetElementValue(k, clusteringAttribute);
            if ( value != valueToFilter) continue;
        }

        double mx = (v1 -minAttr1)/range1;
        double my = (v2 -minAttr2)/range2;
        LocalPoint p;

        p.p[0] = mx;
        p.p[1] = my;
        points->push_back(p);

    }


    vtkSmartPointer<vtkMutableUndirectedGraph> g =
        vtkSmartPointer<vtkMutableUndirectedGraph>::New();

    // added the amount of points as vertices
    vtkIdType ids[points->size()];
    for(unsigned int i = 0; i < points->size(); i++){
        ids[i] = g->AddVertex();
    }
    // Create a fully connected graph
    vtkSmartPointer<vtkDoubleArray> weights = vtkSmartPointer<vtkDoubleArray>::New();
    weights->SetNumberOfComponents(1);
    weights->SetName("Weights");
    //
    for(unsigned int i = 0; i < points->size(); i++){
        // since it is undirected, we don't need to go from 0
        LocalPoint pi = points->at(i);

        for(unsigned int j = i+1; j < points->size(); j++)
        {
               g->AddEdge( ids[i], ids[j]);
               LocalPoint pj = points->at(j);

               double distance = sqrt( (pi.p[0] - pj.p[0])*(pi.p[0] - pj.p[0])  + (pi.p[1] - pj.p[1])*(pi.p[1] - pj.p[1]) );
               weights->InsertNextValue(distance);
        }
    }
    // The weights are the distance between the two points
     g->GetEdgeData()->AddArray(weights);
     // Setup the minimum spanning tree filter
  vtkSmartPointer<vtkBoostPrimMinimumSpanningTree> minimumSpanningTreeFilter = vtkSmartPointer<vtkBoostPrimMinimumSpanningTree>::New();
  minimumSpanningTreeFilter->SetOriginVertex(ids[0]);
  minimumSpanningTreeFilter->SetInputData(g);
  minimumSpanningTreeFilter->SetEdgeWeightArrayName("Weights");

  // Compute the minimum spanning tree
  minimumSpanningTreeFilter->CreateGraphVertexIdArrayOn();
  minimumSpanningTreeFilter->Update();
  minimumSpanningTree->ShallowCopy(minimumSpanningTreeFilter->GetOutput());

  vtkSmartPointer<vtkEdgeListIterator> it = vtkSmartPointer<vtkEdgeListIterator>::New();
  minimumSpanningTree->GetEdges(it);

  while(it->HasNext()){
      vtkEdgeType e = it->Next();
      LocalPoint pi = points->at(e.Source);
      LocalPoint pj = points->at(e.Target);
      double distance = sqrt( (pi.p[0] - pj.p[0])*(pi.p[0] - pj.p[0])  + (pi.p[1] - pj.p[1])*(pi.p[1] - pj.p[1]) );
      edgeDistances->push_back(distance);
  }

  std::sort(edgeDistances->begin(), edgeDistances->end());
  return minimumSpanningTree;
}




