#ifndef WILKINSONSCAGNOSTICS_H
#define WILKINSONSCAGNOSTICS_H


#include "splomwrapper.h"

/*
   Tukeys orginally proposed the Scagnostics approach. About 30! years ago.
   Wilkinson et al. [1] summarizes Tukeys approach.

   The feature measures depend on geometric graphs.
        * The convex hull, alpha shape & minimum spanning tree need to be calculated.
          VTK has implementation for all this 3 graphs.


  First we have to detect outliers, and remove them
     ** An outlier to be a vertex whose adjacent edges in the MST all have a weight (lenght) greater than w
     ** w = q_75 + 1.5 (q_75 - q_25) where
     **  q_75 is the 75th percentile of the MST edge lenghts
     ** and q_75 - q_25 is the interquartile range of the edge lenghts

  Measure 1. Outlying : c_outlying = length(T_outliers)/ length(T)
         where length(A) is the sum of the length of it edges

  Measure 2. Skew:
         c_skew = (q_90 - q50) / (q_90 - q10)

  Measure 3. Clumpy:
         c_clumpy = max_j [ 1 - max[length(ek)]/ length[e_j]]  // This one needs the euclidean distance btw all points
             length(e) = Euclidean distance between its vertices

  Measure 4. Sparce:
       c_sparse = q_90
       Also the alpha value for the alpha shape...

  Measure 5: Striated:



  Measure 6: Convex:

  Measure 7: Skinny:

  Measure 8: Stingy:

  Measure 9: Monotonic:
          squared Spearman correlation coefficient

  IDEA 1: Use all the points
  IDEA 2: Create from the scatterplot the points needed to draw
  [1] High-Dimensional Visual Analytics: Interactive Exploration Guided by Pairwise Views of Point Distributions


  // The convex hull based on http://www.paraview.org/Wiki/VTK/Examples/Cxx/PolyData/PointsProjectedHull
*/

#include <vtkTree.h>
#include <vtkSmartPointer.h>
#include <vtkBoostPrimMinimumSpanningTree.h>
#include <vtkEdgeListIterator.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkDataSetAttributes.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkDataArray.h>
#include <vtkVertexListIterator.h>
#include <vtkInEdgeIterator.h>
#include <vtkOutEdgeIterator.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkConvexHull2D.h>

#include <vtkDelaunay2D.h>
#include <vtkCell.h>
#include <vtkTriangle.h>
#include <vtkFeatureEdges.h>
#include <vtkLine.h>

#include "cuda/density.h"

class WilkinsonScagnostics
{
public:
    WilkinsonScagnostics();
    ~WilkinsonScagnostics();
    void SetSplomGenerator(SPLOMWrapper* _wrapper);


    vtkSmartPointer<vtkTree> CreateMST(int attr1, int attr2, int clusteringAttribute, int valueToFilter, vector<float> *edgeDistances, vector<LocalPoint> *points, bool withBinning = true);

    // Weight defined as threshold were a vertex is eliminated from
    double GetW(vector<float>* edgeDistances);
    // "We consider an outlier to be a vertex whose adjacent edges in the MST all have a weight ( length ) greater than w
    void GetOutliers(vtkSmartPointer<vtkTree> mst, vector<int>* idsOfOutliers, double w, vector<LocalPoint>* points);


    double GetRVM(int attr1, int attr2, int clusteringAttribute, int valueToFilter);

    vtkSmartPointer<vtkTree> CreateMSTWithBinning(int attr1, int attr2, int clusteringAttribute, int valueToFilter, vector<float> *edgeDistances, vector<LocalPoint> *binPoints);
    vtkSmartPointer<vtkTree> CreateMSTWithPoints(int attr1, int attr2, int clusteringAttribute, int valueToFilter, vector<float> *edgeDistances, vector<LocalPoint> *points);

    double SkewednessMeasure(vtkSmartPointer<vtkTree> mst,vector<int>* idsOfOutliers, vector<LocalPoint> *points );
    double SparsenessMeasure(vtkSmartPointer<vtkTree> mst,vector<int>* idsOfOutliers, vector<LocalPoint> *points );
    double OutlyingMeasure(vtkSmartPointer<vtkTree> mst,vector<int>* idsOfOutliers, vector<LocalPoint> *points );

    double StringyMeasure(vtkSmartPointer<vtkTree> mst,vector<int>* idsOfOutliers );
    double StriateMeasure(vtkSmartPointer<vtkTree> mst,vector<int>* idsOfOutliers, vector<LocalPoint> *points );
    double ClumpyMeasure(vtkSmartPointer<vtkTree> mst,vector<int>* idsOfOutliers, vector<LocalPoint> *points );

    // Convex, Skinny & Monotonic

    double ConvexMeasure(vtkSmartPointer<vtkPolyData> alphaShape, vtkSmartPointer<vtkPolyData> convexHull);
    double SkinnyMeasure(vtkSmartPointer<vtkPolyData> alphaShape);
    double MonotonicMeasure(vector<int> *idsOfOutliers, vector<LocalPoint> *points);

    double GetArea(vtkSmartPointer<vtkPolyData> polydata);
    double GetPerimeter(vtkSmartPointer<vtkPolyData> polydata);
    vtkSmartPointer<vtkPolyData> GetAlphaShape(float alphaValue, vector<int> *idsOfOutliers, vector<LocalPoint> *points);
    vtkSmartPointer<vtkPolyData> GetConvexHull(vector<LocalPoint> *points);
    vtkSmartPointer<vtkPolyData> GetConvexHull(vector<int>* idsOfOutliers, vector<LocalPoint> *points);

private:
    SPLOMWrapper* wrapper;


};

#endif // WILKINSONSCAGNOSTICS_H
