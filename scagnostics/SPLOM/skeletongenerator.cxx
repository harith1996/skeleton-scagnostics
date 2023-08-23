#include "skeletongenerator.h"

SkeletonGenerator::SkeletonGenerator()
{
    //std::cout<< "Constructor " << std::endl;
    originalImage = (short*) malloc(sizeof(short)*512*512);
    originalScatterplot = (short*) malloc(sizeof(short)*512*512);

    thresStart = 0.1;
    thresStride = 0.05;
    thresLimit = 1.0;
    length = -99999;
    InitializeCUDAMemory(this->siteParam, this->outputFT, this->outputSkeleton,
                         this->outputTopo, this->skelDT, 1024);
}


SkeletonGenerator::~SkeletonGenerator(){
    std::cout << "Destructor " << std::endl;
    DeinitializationCudaMemory(this->siteParam, this->outputFT, this->outputSkeleton,
                                      this->outputTopo, this->skelDT);
}

void SkeletonGenerator::SetThreshParams(float start, float stride, float limit){
    this->thresStart = start;
    this->thresStride = stride;
    this->thresLimit = limit;
}

void SkeletonGenerator::GetSkeleton(short loadedImage[], int size, short texture[], int type, vector<LocalPoint> *centerlinePoints, string prefix){


    //Create the scatter plot
    // draw is 0
    // not draw is 255
    int sz = size;
    int nx = sz, ny = sz;
    short int* image = (short int*) malloc(sizeof(short int)*sz*sz);//[sz*sz];
    short int* mask = (short int*) malloc(sizeof(short int)*sz*sz);//[sz*sz];

    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            bool hasColor = false;

            for(int k = 0; k < 3; k++)
                if ( loadedImage[ (i+j*sz)*3 + k] != 0 )
                    hasColor = true;

            image[i*sz +j] = (hasColor)?0 :255;
       }
    }


    CreateLargestConnectedComponent(image, sz, mask, prefix);

    fboSize = skelft2DSize(nx,ny);

    int maxTexSize = 1024;
    memset(siteParam,0,maxTexSize*maxTexSize*sizeof(float));
    memset(outputFT,0,maxTexSize*maxTexSize*2*sizeof(short));

    xm=ym=nx; xM=yM=0;

    for(int i=0;i<nx;++i)
      for(int j=0;j<ny;++j)
      {
        originalScatterplot[j*fboSize+i] = image[i + j*nx];
        if (!image[i + j*nx])
        {
           originalImage[j*fboSize+i] = 0;
           siteParam[j*fboSize+i] = 1.0f;
           xm = min(xm,i); ym = min(ym,j);
           xM = max(xM,i); yM = max(yM,j);
        }
        else {
            originalImage[j*fboSize+i] = 255;
            siteParam[j*fboSize+i] = 0.0f;
        }
    }

    xM = nx-1; yM = ny-1;

    short computed[fboSize*fboSize*3];
    memset(computed, 0,fboSize*fboSize*3*sizeof(short));

    int aggregated[size*size*3]; // 3 because of rgb

    if ( length == -99999){
        // not initialized value,
        length = 200;
    }
    length = 1000;

    float startOfThreshold = this->length*this->thresStart;

    //std::cout << "length " << length << std::endl;
    //std::cout << "start of thres? " << startOfThreshold << std::endl;
    int n = 0;
    float currentTreshold = startOfThreshold;
    centerlinePoints->clear();
    // How to select a proper threshold (?)

    // Image In Memory is processed
    // and then run each time..

    //std::cout << "current "  << currentTreshold << "... ? " << this->length << std::endl;

    //currentTreshold = 20;
    while( currentTreshold < this->thresLimit*this->length){
    //for(int k = 2; k < 2 + n; k++){
       if ( currentTreshold == startOfThreshold){
           ProcessImageInMemory(startOfThreshold, prefix);
        }
        else {
            // The length is the maximum threshold that we can get
            // but it seems that there is a lot that eliminates the centerline
            ChangeOfThreshold(currentTreshold);
        }
        GetSkeletonInformation(type, computed, sz);
        for(int i=0;i<nx;++i)
          for(int j=0;j<ny;++j)
          {
              int id1 = j*size + i;
              int id2 = j*fboSize + i;

              if ( currentTreshold == startOfThreshold)
              {
                  aggregated[id1*3 + 0] = 0;
                  aggregated[id1*3 + 1] = 0;
                  aggregated[id1*3 + 2] = 0;
              }

              // the two options, either skeleton or
              if ( computed[id2*3 + 0] == 255 && computed[id2*3 + 1] == 255  && computed[id2*3 + 1] == 255 ){

                  if ( originalImage[id2] == 0){
                      aggregated[id1*3 + 0] += computed[id2*3 + 0];
                      aggregated[id1*3 + 1] += computed[id2*3 + 1];
                      aggregated[id1*3 + 2] += computed[id2*3 + 2];
                  }


              }
              else if ( computed[id2*3 + 0] == 255 && computed[id2*3 + 1] == 0  && computed[id2*3 + 1] == 0 ){
                  aggregated[id1*3 + 0] += computed[id2*3 + 0];
                  aggregated[id1*3 + 1] += computed[id2*3 + 1];
                  aggregated[id1*3 + 2] += computed[id2*3 + 2];
              }

          }

        this->length = 1000;

        currentTreshold += this->length*this->thresStride;
        n++;
    }

    //std::cout << n << std::endl;
    // Once the image is processed, we have two possible cases that are drawn
    // 1) Is the border information, that is the red(?), on where the area was created
    // 2) The actual centerline info...
    //
    if (n == 0) n = 1; // there is an issue when no centerline is found, and therefore the length is set to -999

    //std::cout << "n " << n << std::endl;
    int filled_pixels = 0;
    for(int i=0;i<nx;++i)
      for(int j=0;j<ny;++j)
      {
          int id1 = j*size + i;


          texture[id1*3 + 0] = aggregated[id1*3 + 0]/n ;
          texture[id1*3 + 1] = aggregated[id1*3 + 1]/n ;
          texture[id1*3 + 2] = aggregated[id1*3 + 2]/n ;

          if ( texture[id1*3 + 0] == texture[id1*3 +1] && texture[id1*3 + 0] > 0){
                   LocalPoint p;
                   p.p[0] = i;    p.p[1] = j;
                   centerlinePoints->push_back(p);

          }

          if ( texture[id1*3 + 0 ] != 0 ) filled_pixels++;
      }

    //std::cout << "Filled pixels are " << filled_pixels << " , " << nx*ny << std::endl;
    free(image);
    free(mask);
}

void SkeletonGenerator::InitializeDisplay(int texsize, short* FT, unsigned char* skel, float* siteParam, int endpoints, short* nendpoints,
                       float* skelDT, float length, float thr, int xm, int ym, int xM, int yM){

    this->imgSize = texsize;
    this->FT = FT;
    this->skel = skel;
    this->siteParam = siteParam;
    this->endpoints = nendpoints;
    this->nendpoints = endpoints;
    this->skelDT = skelDT;
    this->length = length;
    this->threshold = thr;
    this->xm = xm;
    this->xM = xM;
    this->ym = ym;
    this->yM = yM;

}


void SkeletonGenerator::SetOriginalImage(short *image, short* scatter){
    orig = image;
    splot = scatter;
}

void SkeletonGenerator::GetSkeletonInformation(int type, short data[], int sizeTexture){
    int show_what = type;
    //int dataSize = sizeTexture;
    for (int i = 0; i < imgSize; ++i)															// Generate visualization texture
        for (int j = 0; j < imgSize; ++j)
        {
            int id = j * imgSize + i;															// Position of (i,j) in the various 1D arrays
            //int id2 = j * dataSize + i;
            float p = siteParam[id];
            short r,g,b; r=g=b=0;
            if (p)																				// Boundary (site): mark as red
            { r = 255; g = b = 0; }
            else																				// Regular (off-boundary): show data
            {
              if (show_what==0)																	// Show either simplified skeleton or FT
              {
                float value = 0;
                if (FT)
                {
                 int ox = FT[id * 2];															// Coords of the closest site to crt pixel
                 int oy = FT[id * 2 + 1];
                 int vid = oy * imgSize + ox;													// Idx in image arrays of site with coords (ox,oy)
                 value = siteParam[vid]/length;													// FT value (color-coded) at current pixel (normalized 0..1 i.e. b/w)
                }
                 r = g = b = static_cast<short>(255*value);
              }
              else if (show_what==1)
              {																					// Show skeleton branches
                  if (skel[id]) r=g=b=255;														// This is a branch non-endpoint
              }
              else if (show_what == 2)//show_what==2
              {
                float value = (skelDT)? skelDT[id] : 0;
                r = g = b = static_cast<short>(255*pow(1-value,0.2));
              }
              else if (show_what == 3){

                // Show the original for debugging purposes
                // Image that was used to create the centerlines.
                short v = orig[id];
                r = g = b = static_cast<short>(v);
              }
              else {
                // Show the scatter plot
                  short v = splot[id];
                  r = g = b = static_cast<short>(v);
              }
            }
            data[id * 3 + 0] = r;
            data[id * 3 + 1] = g;
            data[id * 3 + 2] = b;
        }

    if (show_what==1)
        for(int i=0;i<nendpoints;i++)
        {
           short x = endpoints[2*i], y = endpoints[2*i+1];
           int  id = y * imgSize + x;

           data[id * 3 + 0] = 0;
           data[id * 3 + 1] = 255;
           data[id * 3 + 2] = 0;
        }

}


// Initialization
void SkeletonGenerator::InitializeCUDAMemory(float*& siteParam, short*& outputFT, unsigned char*& outputSkeleton, short*& outputTopo, float*& skelDT, int size)
{
    skelft2DInitialization(size);

    cudaMallocHost((void**)&outputFT,size*size*2*sizeof(short));
    cudaMallocHost((void**)&outputSkeleton,size*size*sizeof(unsigned char));
    cudaMallocHost((void**)&siteParam,size*size*sizeof(float));
    cudaMallocHost((void**)&outputTopo,size*size*2*sizeof(short));
    cudaMallocHost((void**)&skelDT,size*size*2*sizeof(float));
    // Let's initialize everything
    cudaMemset(outputFT, 0, size*size*2*sizeof(short));
    cudaMemset(outputSkeleton, 0, size*size*sizeof(unsigned char));
    cudaMemset(siteParam,0, size*size*sizeof(float));
    cudaMemset(outputTopo,0, size*size*2*sizeof(short));
    cudaMemset(skelDT,0, size*size*2*sizeof(float));
    //*************************************
}


// Deinitialization
void SkeletonGenerator::DeinitializationCudaMemory(float* siteParam, short* outputFT, unsigned char* outputSkeleton, short* outputTopo, float* skelDT)
{
    skelft2DDeinitialization();

    cudaFreeHost(outputFT);
    cudaFreeHost(outputSkeleton);
    cudaFreeHost(siteParam);
    cudaFreeHost(outputTopo);
    cudaFreeHost(skelDT);
}

//oid SkeletonGenerator::CreateLargestConnectedComponent(short *image, int sz, short* mask){
void SkeletonGenerator::CreateLargestConnectedComponent(short* image, int sz, short *mask, string prefix){

    //int radius = 1;
     SPLOMImageType::Pointer imageITK = SPLOMImageType::New();
    SPLOMImageType::RegionType region;
    SPLOMImageType::IndexType start;

    start[0] = 0; start[1] = 0;
    SPLOMImageType::SizeType size;
    size[0] = sz; size[1] = sz;

    region.SetIndex(start);
    region.SetSize(size);

    imageITK->SetRegions(region);
    imageITK->Allocate();

    itk::ImageRegionIteratorWithIndex< SPLOMImageType>  it(imageITK, imageITK->GetLargestPossibleRegion());

    int drawn = 0;
    while(!it.IsAtEnd()){

         SPLOMImageType::IndexType idx = it.GetIndex();
        int pos = idx[0] + idx[1]*sz;

        if (image[pos] == 0){
           it.Set(1);
           drawn++;
        }
        else it.Set(0);
        ++it;
    }


    typedef itk::Image< unsigned short, 2> OutputImageType;
    typedef itk::ConnectedComponentImageFilter < SPLOMImageType, OutputImageType >    ConnectedComponentImageFilterType;

    ConnectedComponentImageFilterType::Pointer connected =   ConnectedComponentImageFilterType::New ();
    connected->SetInput(imageITK);
    connected->Update();

    typedef itk::RelabelComponentImageFilter<OutputImageType,OutputImageType> RelabelFilterType;

    RelabelFilterType::Pointer relabel = RelabelFilterType::New();
    relabel->SetInput(connected->GetOutput());
    relabel->Update();

    OutputImageType::Pointer CCMask = relabel->GetOutput();

    itk::ImageRegionIteratorWithIndex<OutputImageType>  ccit(CCMask, CCMask->GetLargestPossibleRegion());
    itk::ImageRegionIteratorWithIndex< SPLOMImageType>  it2(imageITK, imageITK->GetLargestPossibleRegion());
    while(!ccit.IsAtEnd()){
        if ( ccit.Get() != 1){
              it2.Set(itk::NumericTraits<SPLOMImageType::PixelType>::ZeroValue());
        }
        else {
            it2.Set( itk::NumericTraits<SPLOMImageType::PixelType>::OneValue());
        }
        ++it2;
        ++ccit;
    }
    //*****************************************

    MorphoDilate(imageITK,2);
    MorphoErode(imageITK ,2);

    itk::ImageRegionIteratorWithIndex< SPLOMImageType>  oit(imageITK, imageITK->GetLargestPossibleRegion());
    while(!oit.IsAtEnd()){
         SPLOMImageType::IndexType idx = oit.GetIndex();
        if ( oit.Get() == 0){
             image[ idx[0] + idx[1]*sz] = 255;
             mask[idx[0] + idx[1]*sz] = 0;
        }
        else {
            image[idx[0]+ idx[1]*sz] = 0;
            mask[idx[0] + idx[1]*sz] = 1;
        }
        ++oit;
    }

}




void SkeletonGenerator::ProcessImageInMemory(float threshold, string suffix) {


    int nendpoints = 2000; //TODO- why this value
    //Skeleton lower threshold
    int inflation_dist = 0;
    //float threshold = 50;  //Default initial value
    /**/
    short int FTtemporal[fboSize*fboSize*2];// = (short int*)malloc(fboSize*fboSize*sizeof(short int));

    //3. Compute FT of the sites (for inflation)
    skelft2DFT(FTtemporal,siteParam,0,0,fboSize,fboSize,fboSize);
    //4. Inflate input shape with some distance (by thresholding its DT). Reuse outputFT to save thresholded-DT
    //void skelft2DDT(short* outputDT, float threshold,								//Compute (thresholded) DT (into pbaTextures[2]) from resident FT (in pbaTextures[1])
    //                short xm, short ym, short xM, short yM)

    //string name = "preDDT_" +suffix + ".ppm";
    //skelft2DSave(FTtemporal, fboSize, fboSize, name.c_str());


    //skelft2DDT(outputFT,inflation_dist,xm,ym,xM,yM);
    skelft2DDTshort(outputFT,inflation_dist,xm,ym,xM,yM);

    //string name1 = "afterDDT_" +suffix + ".ppm";
    //skelft2DSave(outputFT, fboSize, fboSize, name1.c_str());

    //Adjust bounds with inflation_dist; careful not to exceed img size (not done!!)
    //Currently inflation distance is set to zero, so no worries about it...
    //The first two steps take 1/3 of the time ...
    //Each step takes about 5 ms..
    xm -= inflation_dist; ym -= inflation_dist;
    xM += inflation_dist; yM += inflation_dist;


    //skelft2DFillHoles((unsigned char*)outputFT,xm+1,ym+1,0,1);		//Slow..
    //Fill background in shape with value 128 (i.e. sth which isn't foreground or background)

    //int skelft2DFill(unsigned char* outputFill, short sx, short sy, short xm, short ym, short xM, short yM, unsigned char fill_value, int limit)
    skelft2DFill((unsigned char*)outputFT,xm+1,ym+1,xm-2,ym-2,xM+2,yM+2,  128 , 100);

    //string name2 = "afterFill_" +suffix + ".ppm";

    //std::cout << "iters " << iter << std::endl;
    //skelft2DSave(outputFT, fboSize, fboSize, name2.c_str());

    // It takes half the time until here...around 30 ms for a 64x64 image, it used 15 until here

    //5. Parameterize boundary of inflated shape into 'siteParam' for skeleton computation


    float length = skelft2DMakeBoundary((unsigned char*)outputFT,xm,ym,xM,yM,siteParam,fboSize,127,false);
    //-----------

    //6. Compute FT of 'siteParam'
    skelft2DFT(outputFT,siteParam,xm,ym,xM,yM,fboSize);

    //std::cout << " Step 6 " << std::endl;
    //7. Skeletonize the FT into 'outputSkeleton'
    skelft2DSkeleton(outputSkeleton,length,threshold,xm,ym,xM,yM);
    //8. Detect endpoints of the skeleton, put them in outputTopo[]
    skelft2DTopology(0,&nendpoints,outputTopo,xm,ym,xM,yM);

    // Compute DT of in-CUDA-memory skeleton
    // There are two dif. DT, one is the DTB and DTS....
    skel2DSkeletonDT(skelDT,xm,ym,xM,yM);
    // only need the original scatterplot
    // saved into "splot"
    // this can be done before....
    SetOriginalImage(originalImage, originalScatterplot);
    InitializeDisplay(fboSize,outputFT,outputSkeleton,siteParam,
                               nendpoints,outputTopo,skelDT,length,threshold,xm,ym,xM,yM);

    // display->GenerateTexture();
    //free(FTtemporal);
}



void SkeletonGenerator::ChangeOfThreshold(float threshold){
    skelft2DSkeleton(skel, length, threshold, xm, ym, xM, yM);
    nendpoints = 1000;
    skelft2DTopology(0,&nendpoints,endpoints, xm, ym, xM, yM);
}

void SkeletonGenerator::MorphoErode(SPLOMImageType::Pointer img, int size){
    typedef itk::BinaryBallStructuringElement< SPLOMImageType::PixelType,2>                  StructuringElementType;
    StructuringElementType erodeSt;
    erodeSt.SetRadius(size);
    erodeSt.CreateStructuringElement();

    typedef itk::BinaryDilateImageFilter <SPLOMImageType, SPLOMImageType, StructuringElementType>
            BinaryDilateImageFilterType;

    BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
    dilateFilter->SetInput(img);
    dilateFilter->SetKernel(erodeSt);
    dilateFilter->SetBackgroundValue(itk::NumericTraits<SPLOMImageType::PixelType>::OneValue());
    dilateFilter->SetForegroundValue(itk::NumericTraits<SPLOMImageType::PixelType>::ZeroValue());
    dilateFilter->Update();

    SPLOMImageType::Pointer out = dilateFilter->GetOutput();

    itk::ImageRegionIterator< SPLOMImageType>  it(img, img->GetLargestPossibleRegion());
    itk::ImageRegionIterator< SPLOMImageType>  git(out, out->GetLargestPossibleRegion());

    while(!it.IsAtEnd()){

        it.Set(git.Get());
        ++it;
        ++git;
    }
}

void SkeletonGenerator::MorphoDilate(SPLOMImageType::Pointer img, int size){
    typedef itk::BinaryBallStructuringElement< SPLOMImageType::PixelType,2>                  StructuringElementType;
    StructuringElementType erodeSt;
    erodeSt.SetRadius(size);
    erodeSt.CreateStructuringElement();

    typedef itk::BinaryDilateImageFilter <SPLOMImageType, SPLOMImageType, StructuringElementType>
            BinaryDilateImageFilterType;

    BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
    dilateFilter->SetInput(img);
    dilateFilter->SetKernel(erodeSt);
    dilateFilter->SetBackgroundValue(itk::NumericTraits<SPLOMImageType::PixelType>::ZeroValue());
    dilateFilter->SetForegroundValue(itk::NumericTraits<SPLOMImageType::PixelType>::OneValue());
    dilateFilter->Update();

    SPLOMImageType::Pointer out = dilateFilter->GetOutput();

    itk::ImageRegionIterator< SPLOMImageType>  it(img, img->GetLargestPossibleRegion());
    itk::ImageRegionIterator< SPLOMImageType>  git(out, out->GetLargestPossibleRegion());

    while(!it.IsAtEnd()){

        it.Set(git.Get());
        ++it;
        ++git;
    }
}

void SkeletonGenerator::MorphoClose(SPLOMImageType::Pointer img, int size){


    typedef itk::BinaryBallStructuringElement< SPLOMImageType::PixelType,2>                  StructuringElementType;
    StructuringElementType erodeSt;
    erodeSt.SetRadius(size);
    erodeSt.CreateStructuringElement();

    typedef itk::BinaryDilateImageFilter <SPLOMImageType, SPLOMImageType, StructuringElementType>
            BinaryDilateImageFilterType;

    BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
    dilateFilter->SetInput(img);
    dilateFilter->SetKernel(erodeSt);
    dilateFilter->SetBackgroundValue(itk::NumericTraits<SPLOMImageType::PixelType>::ZeroValue());
    dilateFilter->SetForegroundValue(itk::NumericTraits<SPLOMImageType::PixelType>::OneValue());
    dilateFilter->Update();

    BinaryDilateImageFilterType::Pointer erodeFilter = BinaryDilateImageFilterType::New();
    erodeFilter->SetInput(dilateFilter->GetOutput());
    erodeFilter->SetKernel(erodeSt);
    erodeFilter->SetBackgroundValue(itk::NumericTraits<SPLOMImageType::PixelType>::OneValue());
    erodeFilter->SetForegroundValue(itk::NumericTraits<SPLOMImageType::PixelType>::ZeroValue());
    erodeFilter->Update();

    SPLOMImageType::Pointer out = erodeFilter->GetOutput();

    itk::ImageRegionIterator< SPLOMImageType>  it(img, img->GetLargestPossibleRegion());
    itk::ImageRegionIterator< SPLOMImageType>  git(out, out->GetLargestPossibleRegion());

    while(!it.IsAtEnd()){

        it.Set(git.Get());
        ++it;
        ++git;
    }
}

void SkeletonGenerator::MorphoOpen(SPLOMImageType::Pointer img, int size){


    typedef itk::BinaryBallStructuringElement< SPLOMImageType::PixelType,2>                  StructuringElementType;
    StructuringElementType erodeSt;
    erodeSt.SetRadius(size);
    erodeSt.CreateStructuringElement();

    typedef itk::BinaryDilateImageFilter <SPLOMImageType, SPLOMImageType, StructuringElementType>
            BinaryDilateImageFilterType;

    BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
    dilateFilter->SetInput(img);
    dilateFilter->SetKernel(erodeSt);
    dilateFilter->SetBackgroundValue(itk::NumericTraits<SPLOMImageType::PixelType>::OneValue());
    dilateFilter->SetForegroundValue(itk::NumericTraits<SPLOMImageType::PixelType>::ZeroValue());
    dilateFilter->Update();

    BinaryDilateImageFilterType::Pointer erodeFilter = BinaryDilateImageFilterType::New();
    erodeFilter->SetInput(dilateFilter->GetOutput());
    erodeFilter->SetKernel(erodeSt);
    erodeFilter->SetBackgroundValue(itk::NumericTraits<SPLOMImageType::PixelType>::ZeroValue());
    erodeFilter->SetForegroundValue(itk::NumericTraits<SPLOMImageType::PixelType>::OneValue());
    erodeFilter->Update();

    SPLOMImageType::Pointer out = erodeFilter->GetOutput();

    itk::ImageRegionIterator< SPLOMImageType>  it(img, img->GetLargestPossibleRegion());
    itk::ImageRegionIterator< SPLOMImageType>  git(out, out->GetLargestPossibleRegion());

    while(!it.IsAtEnd()){

        it.Set(git.Get());
        ++it;
        ++git;
    }
}
