#ifndef SKELETONVIS_H
#define SKELETONVIS_H


/*
   The following class has 3 modes:

   Skeleton Vis:
      SPLOM & Skeleton Visualization

   Categorical Histogram:
      As the name implies Categorical Histogram
            categorical attributes could be shown using a 2D histogram. The columns of the histogram represent the individ-
            ual categorical dimensions, while the rows represent the individual categories within each category.
            Columns shall be sorted by correlation between dimensions and shall be interactively adjustable.
            Rows shall be sorted by size (number of subjects falling into the categories) or by maximally match-
            ing neighboring columns (in terms of correlation) and, again, interactively adjustable. Number of
            occurences/subjects per category shall be color coded using a suitable transfer function.


   Skeleton Explorer:
      Similarity between different skeletons & grouping ..


   // Methods that need to be modified accordingly:
       -> Interaction Methods
             * Mouse Press Event
             * Mouse Double Click Event
             * Mouse Move Event
             * Mouse Release Event
             * wheel event
       -> Painting methods
             * PaintGL
        -> Change MultiVis:
             * Reset scaling & translation

*/

#include <QGLWidget>
#include <QMouseEvent>
#include <QKeyEvent>
#include <iostream>
#include <QTime>
#include <QWidget>
#include <GL/glut.h>
#include <math.h>
#include <utility>      // std::pair
#include <cstring>
#include "../skeletongenerator.h"
#include "../splomwrapper.h"
#include "../memoization.h"
#include "../mathhelper.h"
#include "../cooccurrencecalculator.h"
#include "../statistics.h"
#include "../wilkinsonscagnostics.h"
#include "transferfunction.h"

#include "transferfunctiondialog.h"
#include <itkDirectedHausdorffDistanceImageFilter.h>
#include <itkFilterWatcher.h>
#include <itkHausdorffDistanceImageFilter.h>
#include "stack"
#include "itkBinaryThinningImageFilter.h"
#include "itkLaplacianSharpeningImageFilter.h"
#include "itkNeighborhoodIterator.h"
#include "itkImageToVTKImageFilter.h"

#include "vtkPolyData.h"
#include "vtkContourFilter.h"
#include "vtkMutableDirectedGraph.h"
#include "vtkDijkstraImageGeodesicPath.h"
#include "vtkGraphToPolyData.h"
#include "density.h"

#include "freeformwidget.h"

using namespace std;


typedef itk::Image< short, 2 > DensityImageType;

class SkeletonVis : public QGLWidget
{
    Q_OBJECT

public:
     SkeletonVis(QWidget *parent = 0);
    ~SkeletonVis();

     void CreateInitialHistogramOrdering();
     void ReorderAccordingToColumn(int col);
     void SetThreshParams(float start, float stride, float limit);

    void ChangeMultiVisMethod(int method);

    void SetSplomGenerator(SPLOMWrapper* wrapper, SPLOMWrapper *drawer);
    void SetCoocurrenceInfo(CooccurrenceCalculator* cooccurCalc);

    void Generate3DTexture(int clusteringVariable, int filterVariable, bool blendingChange);
    void GenerateDensityContours(int clusteringVariable, int filterVariable );

    // Generates the Centerline Texture for a single variable in a single class
    void Generate3DCenterlineTexture(int clusteringVariable, int filterVariable);

    //
    void Generate3DCenterlineTextureGrouped(int clusteringVariable, int totalVariables);
    void ChangeContourOpacity(float val);

    void SetMainTF(TransferFunction* tf){ mainTf = tf; }

    void GenerateSkeleton(int attr1, int attr2, int clusteringAttribute, int valuetoFilter, int type, short texture[]);

    void SetMappingProperties(int map){
        this->mapProperties = map;

        Redraw();
    }
    void SkeletonPaint();
    void CategoricalHistPaint();
    void MDSPaint();
    void ScagnosticsPaint();

    void ChangeMultiVisMeasure(int idx);


    void CompareModel(bool compare, float a, float b, float c );
    void CreateMDSProjection();

    void CreateMDSProjectionFrechet();
    void CreateMDSProjectionScagnostics();
    void CreateMDSProjectionTatu();
    void CreateMDSProjectionHausdorff();
    void DisplayMDSCenterlineWithInfo();


    void ChangeSimilarityThreshold(double value);
    void ClearMemory();
    int GetDesiredResolution(){ return latestDesiredSize; }
    void GetGridsPairToGen(  vector< pair<int, int> >* res){
        for(unsigned int k = 0; k < gridsToGen.size(); k++)
            res->push_back(gridsToGen.at(k));
    }
    void GetCenterlinesPairToGen(  vector< pair<int, int> >* res){
        for(unsigned int k = 0; k < centerlineToGen.size(); k++)
            res->push_back(centerlineToGen.at(k));
    }

    bool* GetHighlightedSamples(){ return highlightedPoints;}

    void MovingOnSkeletonGrid(double pX, double pY);
    void MovingOnSkeletonCenterline(float pX, float pY);
    void GetWorldCoordinates(double x, double y, double *xAtPress, double *yAtPress);

    void PassAsGeneratedGrid();
    void SetNames(vector<string>* attributes);
    void Redraw();

    void ChangeSkeletonLineWidth(float newWidth){
        centerlineWidth = newWidth;
    }
    void ZoomSkeleton();
    void MovementSkeleton();

    void ShowTF();
    float alpha;
    float beta;
    bool useSubset;
    bool pause;

    void ChangeContourColor(QColor val);
    void ChangeCenterlineColor(QColor val);

    void ToggleDensityContour(bool densityBool);
    void SetDensityContourValue(float value);

    void ChangeFromSimilarityColor(QColor c);
    void ChangeToSimilarityColor(QColor c);

    void ChangeDoubleClickMode(int mode){ this->doubleClickMode = mode; Redraw(); }
    void SetMultiRes(bool multi){ this->multiResolution = multi;}

    void ToggleFreeForm();
public slots:
    void NewFreeFormModel();

signals:
    void NewSkeletonMade();
    void ChangeResolution();
    void SamplesOnGridSelected();
    void SamplesOnHistogramSelected();
protected:
   void paintGL();
   void initializeGL();
   void resizeGL(int width, int height);
   virtual void mousePressEvent(QMouseEvent* event);
   virtual void mouseDoubleClickEvent(QMouseEvent* event);
   virtual void mouseMoveEvent(QMouseEvent *event);
   virtual void mouseReleaseEvent(QMouseEvent *event);
   //virtual void keyPressEvent(QKeyEvent *event);
   virtual void wheelEvent(QWheelEvent *event);

private:

   enum MultiDimVisMethod { SPLOM, CategoricalHist, MDS, Scagnostics };
   enum MultiDimVisMeasure { Wilkinson, Frechet, Tatu, Hausdorff };

   int currentMultiVisMethod;

   bool InRange(float px, float py, float startX, float endX, float startY, float endY);
   bool GetCenterlinePoints(int attr1, int attr2, int clusteringAttribute, int valuetoFilter, vector<LocalPoint>* points);
   bool GetSimplifiedCenterlinePoints(int attr1, int attr2, int clusteringAttribute, int valuetoFilter, vector<LocalPoint>* points);
   bool GetLargestCenterlinePoints(int attr1, int attr2, int clusteringAttribute, int valuetoFilter, vector<LocalPoint>* points);

   bool GetSkeletonSegments(int attr1, int attr2, int clusteringAttribute, int valuetoFilter, vector<SkeletonSegment*>* segments, double *maxD, double *minD, double *maxDistance);

   void CreateOnePixelWidthSkeleton(vector<LocalPoint>* inputSkeleton, vector<LocalPoint>* outputSkeleton);
   void MapPoints(vector<LocalPoint>* mappedPoints, int attr1, int attr2, int clusteringAttribute, int valueToFilter);

   void GetEndPointsAndBifurcations(vector<LocalPoint>* centerlinePoints, vector<int>* endPoints, vector<int>* bifurcations);
   int GetStartingPoint(vector<LocalPoint>* centerlinePoints, vector<int>* endPoints, bool seen[]);

   void CreateParentsAndTangents(vector<SkeletonSegment*>* currentTree, int* totalCenterlinePoints);

   void CreateTreeDensity(vector<SkeletonSegment*>* currentTree , double *maxDensity, double *minDensity, double *maxDistance,
                          int totalCenterlinePoints, int attr1, int attr2, int clusteringAttribute, int valueToFilter);


   double GetGeodesicDistance( vector< pair<LocalPoint,int> >* points, vector<int>* ids );


   void GetLargestContinuousCenterline(vector<SkeletonSegment *> *currentTree, vector<int> *largestCenterline);

   void AddToVector( vector< pair<int, int> >* vec, pair<int, int> loc);
   void GetPointsInGrid();

   void CreateDensityOnCenterline(int attr1, int attr2, int clusteringAttribute, int valueToFilter,
                                    vector<LocalPoint> *simplifiedCenterline, vector<SkeletonSegment *> *currentTree, double *maxDensity, double *minDensity, double *maxDistance, vector<int> *largestCenterline);


   double GetRightAngle(LocalPoint firstTangent, LocalPoint firstPoint);

   void OrderSkeletonsAccordingToFrechetDistance();
   void OrderSkeletonsAccordingToTatuDistance();
   void OrderSkeletonsAccordingToHausDorffDistance();



   void mouseSkeletonMovement(QMouseEvent *event);
   void mouseHistogramMovement(QMouseEvent *event);
   void mouseScagnosticsMovement(QMouseEvent *event);
   void mouseMDSMovement(QMouseEvent* event);

   void mouseSkeletonDoubleClick(QMouseEvent *event);
   void mouseHistogramDoubleClick(QMouseEvent *event);
   void mouseScagnosticsDoubleClick(QMouseEvent *event);

   void CreateITKImage(int attr1, int attr2, HImageType::Pointer img, bool DT);
   void CreateITKModelImage(HImageType::Pointer img, bool DT);

   double maxInS(HImageType::Pointer skeleton, HImageType::Pointer dt);

   void mouseSkeletonPress(QMouseEvent *event);
   void mouseHistogramPress(QMouseEvent *event);
   void mouseScagnosticsPress(QMouseEvent *event);
   void mouseMDSPress(QMouseEvent* event);

   void mouseSkeletonRelease(QMouseEvent *event);
   void mouseHistogramRelease(QMouseEvent *event);
   void mouseScagnosticsRelease(QMouseEvent *event);
   void mouseMDSRelease(QMouseEvent* event);

   void SwapRowsTopHistogram();
   void SwapRowsBottomHistogram();

   double Gaussian2D(LocalPoint center, LocalPoint point, double amplitude, const double stddev);

   double GetCurvatureWeight(  vector< pair<LocalPoint,int> >* points, vector<int>* ids);

   void SwapColsRightHistogram();
   void SwapColsLeftHistogram();

   void ChangeSPLOMOrderingTatu();
   void ChangeSPLOMOrderingHausdorff();
   void ChangeSPLOMOrderingFrechet();

  void DisplayCenterlineWithInfo();
  void HighlightSimilarCenterlines();

  void CreateScagnosticsMeasures();
  void CreateTatuMeasures();

  //int GetStartingPoint(vector<LocalPoint>* centerlinePoints, vector<int>* endPoints, bool seen[]);

  int  CreateTree2(int startingPoint, vector<SkeletonSegment*>* currentTree, bool seenEndPoints[],
                               vector<LocalPoint> centerlinePoints, vector<int> endPointsIndices,bool visited[], bool debugFunction = false);


   GLint viewport[4];

   GLubyte* tex3DCombined;
   GLuint texture3DCombined;


   float  scale, transX, transY;
   bool   isLeftMouseActive, isRightMouseActive;
   int    oldMouseX, oldMouseY;
   int    winSize;
   int    imgSize;

   SPLOMWrapper* _splomGenerator;
   SPLOMWrapper* _splomDrawer;

   CooccurrenceCalculator* _cooccurCalc;

   bool changeInTexture;
   //bool changeInTextureSkel;

   vector<string> names;
   SkeletonGenerator generator;
   //***********************************
   pair<int,int> selectedGrid;
   pair<int,int> selectedCenterline;
   pair<int,int> selectedCenterlineTree;


   // In order to check whether we have to change the resolution or what to generate...
   // then look at the memoization ....
   vector< pair<int,int> > gridsToGen;
   vector< pair<int,int> > centerlineToGen;
   int latestDesiredSize;
   vector< pair<int, int> > currentlyDrawn;

   set< pair<int,int> > gridsGen;
   set< pair<int,int> > centerlinesGen;


   vector<Memoization*>  centerlineMemory;

   int totalSeenBefore;
   // 0 means no selection, 1 means selection in grid, 2 means selection in centerline
   int selectingOnSkeletonType;
   bool runningCenterline;
   //*****************

   //***************************
   bool areaSelectionInGrid;
   vector<LocalPoint> areaSelectionOutlineInGrid;
   float currentMovingDirectionInGrid[3];
   bool* highlightedPoints;

   TransferFunction* mainTf;
   // For re-ordering
   vector< int > columnOrder;
   vector< pair<int,int> > rowOrder;
   vector< pair<int,int> > rowBelonging;
   int currentRowSelected;
   int currentColSelected;
   int actionOnHistogram;

   int hoverRow;
   int hoverColumn;
   //
   float hoverSphereLoc[2];
   float hoverRowLoc[2];

   vector<LocalPoint> brushingArea;
   float brushHistSize;
   float brushSkeletonSize;
   TransferFunctionDialog* tfDialog;
   FreeFormWidget* freeFromDialog;

  //***********

   vector< pair<LocalPoint,LocalPoint> > selectedCenterlinePoints;
   //vector<LocalPoint> selectedCenterlinePoints;
   int mapProperties;


   float thresStart;
   float thresStride;
   float thresLimit;


   bool useDensityContour;
   float densityContourPercentage;
   vector< vtkSmartPointer<vtkPolyData> > contours;

   vector< pair<int,int> > similarCenterlines;

   float similarityTreshold;
   float contourOpacity;
   QColor colorContour;
   QColor colorCenterline;
   //******
   QColor fromBackground;
   QColor toBackground;
   vector<float> similarityDistancesFromSelected;
   //****

   bool displayCenterlineLine;
   float centerlineWidth;

   //
   vector<LocalPoint> mdsPoints;
   bool useMDSCenterlines;
   //
   int typeOfMeasureMD; // 0 for wilkinson, 1 for frechet, 2 for Tatu , 3 for Yates Grouping

   //
   double *scagnosticsValues;
   double *tatuValues;
   double scagnosticsMaxsMins[18];

   vector<  vector<pair<int, float> > > orderedMeasures;

   //bool changeSPLOMOrdering;
   int doubleClickMode;

   vector<int> reOrderedSplom;
   vector<int> orderAccordingToSimilarity;


   bool modelComparison;
   // a x^2 +  b^x +  c
   float a_model;
   float b_model;
   float c_model;

   // Multi Resolution is basically
   // when we zoom in, to recalculate the Scatterplots to a better
   // resolution. ... i.e. focus and context
   bool multiResolution;

   bool useFreeFormModel;


};

#endif // SKELETONVIS_H
