#ifndef SKELETONGENERATOR_H
#define SKELETONGENERATOR_H


#include "cuda/skelft.h"
#include "./skeleton/field.h"

#include <cuda_runtime_api.h>

#include "itkImage.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkConnectedComponentAlgorithm.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileWriter.h"
#include "itkVTKImageIO.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkNumericTraits.h"
#include "image2d.hpp"
#include "multival.hpp"
#include "dataset.hpp"
#include "projection2d.hpp"

#include <chrono>


using namespace std::chrono;
typedef itk::Image< short, 2> SPLOMImageType;

class SkeletonGenerator
{
public:
    SkeletonGenerator();
    ~SkeletonGenerator();

    //Loads an RGB Image, masks it if it has color and generates the centerline out of it
    Image2D<unsigned> CreateScatterplot(const DataSet& dataSet, Projection2D projection, size_t width, size_t height);
    void GetColorFromQuadTree(float qx, float qy, float px, float py, float* quadTreeMetadata, float* mappedPoints, int maxLevel, float color[], int qts[],
                                           int blocksSeen[], double minD);
    void GetSkeleton(Image2D<unsigned> loadedImage, int size, Image2D<char> &skelImg, Image2D<unsigned>& maskedPlot, vector<Int2> *skeletonPoints, int thres=1); //vector<Int2 > *centerlinePoints
    void CreateOnePixelWidthSkeleton(vector<Int2>* inputSkeleton, vector<Int2>* outputSkeleton, int sz);
    void GetEndPointsAndBifurcations(vector<Int2>* centerlinePoints, vector<int>* endPoints, vector<int>* bifurcations);

    void SetThreshParams(float start, float stride, float limit);

private:

    bool initialized;
    void InitializeDisplay(int texsize, short* FT, unsigned char* skel, float* siteParam, int endpoints, short* nendpoints,
                           float* skelDT, float length, float thr, int xm, int ym, int xM, int yM);

    void CreateLargestConnectedComponent(short* image, int sz, short *mask);

    void ChangeOfThreshold(float threshold);

    void ProcessImageInMemory(float threshold);

    void GetSkeletonInformation(int type, short data[], int sizeTexture);
    double GetR(const DataSet& dataSet, Projection2D projection, int sizeT);

    void SetOriginalImage(short* image, short* scatter);

    void MorphoErode(SPLOMImageType::Pointer img, int size);
    void MorphoDilate(SPLOMImageType::Pointer img, int size);

    //////////////////////////////////////////////

    /*Variables used to maintain The Telea approach */
    void InitializeCUDAMemory(float*& siteParam, short*& outputFT, unsigned char*& outputSkeleton, short*& outputTopo, float*& skelDT, int size);

    void DeinitializationCudaMemory(float* siteParam, short* outputFT, unsigned char* outputSkeleton, short* outputTopo, float* skelDT);

    //Boundary parameterization: fboSize x fboSize. value(i,j) = boundary-param of (i,j) if on boundary, else unused
    float* siteParam;

    // The original image with the fboSize x fboSize
    short int* originalScatterplot;
    short int* originalImage;


    //1-point FT output: fboSize x fboSize x 2. value(i,j) = the two coords of closest site to (i,j)
    short* outputFT;


    //Skeleton: fboSize x fboSize. value(i,j) = skeleton at (i,j) (i.e. zero or one)
    unsigned char* outputSkeleton;

    short* outputTopo;
    // Distance Transform from the skeleton, different than the Distance Transform from Border used for the
    // centerline calculation
    float* skelDT;

    // minimums and maximums of x and y, so we have a bounding box
    short xm,ym,xM,yM;
    int fboSize;

    short* orig;
    short* splot;

    int    imgSize;												//Effective size (pow 2) of texture images to be displayed

    short* FT;
    unsigned char* skel;

    int    show_what;
    float  threshold;
    short* endpoints;
    int    nendpoints;
    bool   tex_interp;
    float length;
    //////////////////////////////////////////////

    float thresStart;
    float thresStride;
    float thresLimit;

    short* computed;
    int* aggregated;

};

#endif // SKELETONGENERATOR_H
