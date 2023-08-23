#ifndef SKELETONGENERATOR_H
#define SKELETONGENERATOR_H


#include "../cuda/skelft.h"
#include "../include/field.h"

#include <cuda_runtime_api.h>

#include "dataset.h"
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
#include <chrono>


using namespace std::chrono;
typedef itk::Image< short, 2> SPLOMImageType;

class SkeletonGenerator
{
public:
    SkeletonGenerator();
    ~SkeletonGenerator();

    //Loads an RGB Image, masks it if it has color and generates the centerline out of it
    void GetSkeleton(short loadedImage[], int size, short finalImage[], int type, vector<LocalPoint>* centerlinePoints, string prefix);
    void SetThreshParams(float start, float stride, float limit);

private:
    void InitializeDisplay(int texsize, short* FT, unsigned char* skel, float* siteParam, int endpoints, short* nendpoints,
                           float* skelDT, float length, float thr, int xm, int ym, int xM, int yM);

    void CreateLargestConnectedComponent(short* image, int sz, short *mask, string prefix);

    void ChangeOfThreshold(float threshold);

    void ProcessImageInMemory(float threshold, string suffix);

    void GetSkeletonInformation(int type, short data[], int sizeTexture);

    void SetOriginalImage(short* image, short* scatter);



    void MorphoClose(SPLOMImageType::Pointer img, int size);
    void MorphoErode(SPLOMImageType::Pointer img, int size);
    void MorphoDilate(SPLOMImageType::Pointer img, int size);
    void MorphoOpen(SPLOMImageType::Pointer img, int size);

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
};

#endif // SKELETONGENERATOR_H
