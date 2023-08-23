

#include "skeletongenerator.h"
#include "cuda/skelft.h"

int comparePoint(const void *a, const void * b){

    float* p1 = (float*)a;
    float* p2 = (float*)b;

    // p1 has x,y, color, bucket
    float v1 = p1[3];// We have already made so we know the bucket
    float v2 = p2[3];
    if ( v1 < v2) return -1;
    if ( v1 == v2) return 0;

    return 1;
}


SkeletonGenerator::SkeletonGenerator()
{
    //std::cout<< "Constructor " << std::endl;
    originalImage = (short*) malloc(sizeof(short)*1024*1024);
    originalScatterplot = (short*) malloc(sizeof(short)*1024*1024);
    computed = (short*) malloc(sizeof(short)*1024*1024);
    aggregated = (int*)malloc(sizeof(int)*1024*1024);

    initialized = false;

    thresStart = 0.1;
    thresStride = 1.0;
    thresLimit = 1.0;
    length = -99999;

}


SkeletonGenerator::~SkeletonGenerator(){
    if (initialized)
    DeinitializationCudaMemory(this->siteParam, this->outputFT, this->outputSkeleton,
                                      this->outputTopo, this->skelDT);
    if (computed != nullptr)free(computed);
    if (aggregated != nullptr) free(aggregated);
}

void SkeletonGenerator::SetThreshParams(float start, float stride, float limit){
    this->thresStart = start;
    this->thresStride = stride;
    this->thresLimit = limit;
}

Image2D<unsigned> SkeletonGenerator::CreateScatterplot(const DataSet& dataSet, Projection2D projection, size_t width, size_t height){

    //  {"Lin","Std","Log","Mnl-Lin","Mnl-Std","Mnl-Log"};

    //int typeScaleAttr1 = _reader->GetScale(attr1);
    //int typeScaleAttr2 = _reader->GetScale(attr2);
    float minD = (1.0/128.0);
    float xmin =  0,ymin = 0, xmax = 1.0, ymax = 1.0;
    //First step is to get the points that are not culled.

    //float r_tmp = GetR(dataSet, projection,width);
    //std::cout << "Min D "<< minD << " r_ tmp? " << r_tmp << std::endl;

    //minD = r_tmp;*/


    int totalNotCulled = 0;
    int totalElements = dataSet.size();
    const int PER_POINT_INFO = 6;
    float* mappedPoints = (float*)malloc(sizeof(float)*totalElements*PER_POINT_INFO); // x, y, valueIndex, bucketIndex, initialIndex, colorIndex

    float DesiredRadius = minD;

    for(int k = 0; k < totalElements; k++){

        auto pair = projection(dataSet, k);

        float v1 = pair[0];
        float v2 = pair[1];

        float mx = v1;
        float my = v2;

        mappedPoints[totalNotCulled*PER_POINT_INFO + 0] = mx;
        mappedPoints[totalNotCulled*PER_POINT_INFO + 1] = my;
        mappedPoints[totalNotCulled*PER_POINT_INFO +2] = 0;
        mappedPoints[totalNotCulled*PER_POINT_INFO +5] = 0;

        int xIndex = 2.0*(mx - xmin)/(xmax-xmin);
        int yIndex = 2.0*(my - ymin)/(ymax-ymin);

        mappedPoints[totalNotCulled*PER_POINT_INFO + 3] =  2*yIndex + xIndex; // bucket index
        mappedPoints[totalNotCulled*PER_POINT_INFO + 4] = k; // initial index
        totalNotCulled++;
    }

    //Lets do a break down up to 5 levels

    /*if ( !useAsDrawer){
        //DesiredRadius = (1.0/sizeT)*3;
        DesiredRadius *= 2;
        //DesiredRadius = GetR(mappedPoints, totalNotCulled );
        //std::cout << "Before? " << minD << " after " << DesiredRadius << std::endl;
    }*/

    int maxLevel = 5;
    int totalInLastLevel = sqrt(pow(4, maxLevel));
    float midDistance = (1.0/totalInLastLevel); // half the size of each division...
    //float minD = 3.0/128.0;

    while ( DesiredRadius*2 > midDistance){
        maxLevel -= 1;
        totalInLastLevel = sqrt(pow(4, maxLevel));
        midDistance = (1.0/totalInLastLevel);
    }

    int totalSize =  (1 - pow(4,maxLevel+1))/(1-4);
    float* quadTreeMetadata = (float*) malloc(sizeof(float)*6*totalSize);

    //std::cout << "Start of metadata " << xmin << ", " << ymin << " -> " << xmax << " , " << ymax << std::endl;

    quadTreeMetadata[0] = xmin;           // Complete xmin
    quadTreeMetadata[1] = xmax;           // Complete xmax
    quadTreeMetadata[2] = ymin;           // Complete ymin
    quadTreeMetadata[3] = ymax;           // Complete ymax
    quadTreeMetadata[4] = 0;              // starting index of points
    quadTreeMetadata[5] = totalNotCulled; // end index of points :: total points not culled are the ones that we are drawing


    // We sort the points at the first level...

    qsort(mappedPoints,totalNotCulled,sizeof(float)*PER_POINT_INFO, comparePoint);
    //Then we create each level of the quadtree and set the points

    for(int level = 1; level < (maxLevel+1); level++){
        //For each level...
        //We work based on the previous level...

        int totalInLevel = pow(4, level);
        int prevStart = (1- pow(4, level-1)) / (1-4);
        int start = (1- pow(4, level)) / (1-4);

        for(int i = 0; i < totalInLevel; i++){
            // This gives the total items in that level...
            // First we need to figure out where it belongs...
            int cellInPrevLevel = i / 4;
            int loc = cellInPrevLevel * 4 + start +  i%4; // Where to save it...
            int parentLoc = prevStart + cellInPrevLevel;

            float xminParent = quadTreeMetadata[parentLoc*6 + 0];
            float xmaxParent = quadTreeMetadata[parentLoc*6 + 1];
            float yminParent = quadTreeMetadata[parentLoc*6 + 2];
            float ymaxParent = quadTreeMetadata[parentLoc*6 + 3];

            int parentStartIndex = quadTreeMetadata[parentLoc*6 + 4];
            int parentEndIndex = quadTreeMetadata[parentLoc*6 + 5];

            int currentStart = parentStartIndex;
            while( mappedPoints[currentStart*PER_POINT_INFO +3] < static_cast<float>(i%4)){
                currentStart++;
                if ( currentStart > parentEndIndex){ currentStart = parentEndIndex; break; }
            }

            int currentEnd = currentStart;
            while( mappedPoints[currentEnd*PER_POINT_INFO + 3] < static_cast<float>(i%4 + 1)){
                currentEnd++;
                if ( currentEnd > parentEndIndex){ currentEnd = parentEndIndex;  break; }
            }



            float midX = 0.5*(xmaxParent - xminParent) + xminParent;
            float midY = 0.5*(ymaxParent - yminParent) + yminParent;

            float valuesX[3] = {xminParent, midX, xmaxParent};
            float valuesY[3] = {yminParent, midY, ymaxParent};
                                  //xlow, xhigh, ylow, yigh
            int positions[4][4] = { {0,1,0,1},
                                    {1,2,0,1},
                                    {0,1,1,2},
                                    {1,2,1,2}};


           quadTreeMetadata[loc*6 + 0] = valuesX[ positions[i%4][0]];
           quadTreeMetadata[loc*6 + 1] = valuesX[ positions[i%4][1]];
           quadTreeMetadata[loc*6 + 2] = valuesY[ positions[i%4][2]];
           quadTreeMetadata[loc*6 + 3] = valuesY[ positions[i%4][3]];
           quadTreeMetadata[loc*6 + 4] = currentStart;
           quadTreeMetadata[loc*6 + 5] = currentEnd;
        }

        // Now we sort the area
        for(int i = 0; i <  totalInLevel; i++){
            int cellInPrevLevel = i / 4;
            int loc = cellInPrevLevel * 4 + start +  i%4; // Where to save it...

            int currentStart = quadTreeMetadata[loc*6 + 4];
            int currentEnd = quadTreeMetadata[loc*6 + 5];
            for(int k = currentStart; k < currentEnd; k++ ){
                int xIndex = 2.0*(mappedPoints[k*PER_POINT_INFO +0]  - quadTreeMetadata[loc*6 +0])/(quadTreeMetadata[loc*6 +1]-quadTreeMetadata[loc*6 +0]);
                int yIndex = 2.0*(mappedPoints[k*PER_POINT_INFO +1]  - quadTreeMetadata[loc*6 +2])/(quadTreeMetadata[loc*6 +3]-quadTreeMetadata[loc*6 +2]);

                mappedPoints[k*PER_POINT_INFO + 3] =  2*yIndex + xIndex; // bucket index
            }
            qsort(&(mappedPoints[currentStart*PER_POINT_INFO]), currentEnd-currentStart,sizeof(float)*PER_POINT_INFO, comparePoint );
        }
    }

    float spacing_x = 1.0 / width;
    float spacing_y = 1.0 / height;
    float upper_cornerX =  0;
    float upper_cornerY =  0.0;

    Image2D<unsigned> img {width, height};

    for(unsigned int  j = 0; j < height; j++) {

        float py = upper_cornerY + spacing_y*j + spacing_y/2.0f;
        for(int i = 0; i < width; i++) {
            float px = upper_cornerX + spacing_x*i + spacing_x/2.0f;

            // For each level find which bucket it belongs to
            //unsigned int color[4] = {0,0,0,0};
            float color[4] = {0,0,0,0};

            int qts[1];
            for(int k = 0; k < 1; k++){
                qts[k] = 0;
            }
            float other = DesiredRadius*1.01;

            int negX =  static_cast<int>(floor((px - other - xmin)/(midDistance)))  - static_cast<int>(floor((px - xmin)/(midDistance)));
            int posX =  static_cast<int>(floor((px + other - xmin)/(midDistance)))  - static_cast<int>(floor((px - xmin)/(midDistance)));
            int difX = posX + negX;

            int negY =  static_cast<int>(floor((py - other - ymin)/(midDistance)))  - static_cast<int>(floor((py - ymin)/(midDistance)));
            int posY =  static_cast<int>(floor((py + other - ymin)/(midDistance)))  - static_cast<int>(floor((py - ymin)/(midDistance)));
            int difY = posY + negY;
            int startOfBlocks[8] = {-1,-1,-1,-1, -1,-1,-1,-1};

            int movX = 0, movY = 0;

            while( abs(movX) <= abs(difX)){

                while(abs(movY) <= abs(difY)){
                    GetColorFromQuadTree(px, py, px + movX*DesiredRadius, py + movY*DesiredRadius, quadTreeMetadata,mappedPoints,
                                         maxLevel, color, qts, startOfBlocks, DesiredRadius);

                    if (difY == 0) break;

                    movY += difY;
                }
                movY  = 0;
                if (difX == 0) break;
                movX += difX;
            }

            img(i,j) += color[3];        
        }
    }


    free(quadTreeMetadata);
    free(mappedPoints);
    return img;
}




void SkeletonGenerator::GetColorFromQuadTree(float qx, float qy,float px, float py, float* quadTreeMetadata, float* mappedPoints,
                                       int maxLevel, float color[],  int qts[],  int blocksSeen[], double minD){

    // qx,qy is the point to see whether it is close by
    // px,py is the point to query
    const int PER_POINT_INFO = 6;

    int startIndex = 0;
    int endIndex = 0;
    int currentPos = 0;
    int actualPos = currentPos;

    for(int k = 0; k < maxLevel - 1; k++){
         //  We start with the first block, and find which block it belongs to, or
         // whether there are no points left
        int start = (1- pow(4, k)) / (1-4);
        actualPos = currentPos + start;
        float xminP = quadTreeMetadata[actualPos*6 + 0];
        float xmaxP = quadTreeMetadata[actualPos*6 + 1];
        float yminP = quadTreeMetadata[actualPos*6 + 2];
        float ymaxP = quadTreeMetadata[actualPos*6 + 3];
        int parentStartIndex = quadTreeMetadata[actualPos*6 + 4];
        int parentEndIndex = quadTreeMetadata[actualPos*6 + 5];

        int xIndex = 2.0*(px - xminP)/(xmaxP-xminP);
        int yIndex = 2.0*(py - yminP)/(ymaxP-yminP);
        int bucketIndex = 2*yIndex + xIndex;
        currentPos = currentPos*4 + bucketIndex;
        if ( parentStartIndex == parentEndIndex){
            startIndex = parentStartIndex;
            endIndex = parentEndIndex;
            break;
        }
        start = (1- pow(4, k+1)) / (1-4);
        startIndex = quadTreeMetadata[(currentPos + start)*6 + 4];
        endIndex = quadTreeMetadata[(currentPos+ start)*6 + 5];
    }

    // No points in the area,
    // Only if there are points we add color, otherwise we leave it as it is...

    // In the no blending, first it finds the one in the current cell,
    // and then it goes to the others cell, that's the reason to the overlap
    // we need a way to define on how we select the points, such that it always select the same...

    if ( startIndex != endIndex) {
        //if this is the case, then there are points in the level...
        for(int k = 0; k < 4; k++){
            if ( blocksSeen[2*k] == startIndex && blocksSeen[2*k +1] == endIndex)
            {
                // then we just return, because we have seen it ...
                return;
            }
        }
         // We went to the last level, and now we look at the points
        for(int k = startIndex ; k < endIndex; k++){
              float x = mappedPoints[k*PER_POINT_INFO + 0];
              float y = mappedPoints[k*PER_POINT_INFO + 1];
              //int val = mappedPoints[k*6 + 2];
              int idx = mappedPoints[k*PER_POINT_INFO + 5];

              float d = pow(qx - x,2.0) + pow(qy-y,2.0);
              if ( d < minD*minD){
                  color[3] += 1;
                  qts[idx] += 1;
                  int tmpColor[3] ={255,255,255};
                  color[0] += tmpColor[0];
                  color[1] += tmpColor[1];
                  color[2] += tmpColor[2];
              }
        }


        // Given that we processed it, now we find the first -1 and replace in the seen blocks
          for(int k = 0; k < 4; k++){
                if (blocksSeen[2*k] == -1){
                    blocksSeen[2*k] = startIndex;
                    blocksSeen[2*k+1] = endIndex;

                    break;
                }
            }
    }
}


void SkeletonGenerator::GetSkeleton(Image2D<unsigned> loadedImage, int size, Image2D<char>& skelImg,  Image2D<unsigned>& maskedPlot, vector<Int2>* skeletonPoints, int thres){

    InitializeCUDAMemory(this->siteParam, this->outputFT, this->outputSkeleton,
                         this->outputTopo, this->skelDT, 1024);
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

            //for(int k = 0; k < 3; k++)
                //if ( loadedImage[ (i+j*sz)*3 + k] != 0 )
            if ( loadedImage(j,i) >= thres )
                  hasColor = true;
            image[i*sz +j] = (hasColor)?0 :255;
       }
    }

    int type = 1; //Get Skeleton Information
    CreateLargestConnectedComponent(image, sz, mask);

    skeletonPoints->clear();
    fboSize = skelft2DSize(nx,ny);

    int maxTexSize = 1024;
    memset(siteParam,0,maxTexSize*maxTexSize*sizeof(float));
    memset(outputFT,0,maxTexSize*maxTexSize*2*sizeof(short));

    xm=ym=nx; xM=yM=0;

    for(int i=0;i<nx;++i)
      for(int j=0;j<ny;++j)
      {
        originalScatterplot[j*fboSize+i] = image[i + j*nx];

        if( mask[j*sz + i]){
            maskedPlot(i,j) = loadedImage(i,j);
        }

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

    //short computed[fboSize*fboSize*3];
    //memset(computed, 0,fboSize*fboSize*3*sizeof(short));
    for(int k = 0; k < fboSize*fboSize*3; k++)
        computed[k] = 0;

    std::cout << "Computed size? " << fboSize*fboSize*3 << std::endl;
    //int aggregated[size*size*3]; // 3 because of rgb
    for(int k = 0; k < size*size*3;k++)
        aggregated[k] = 0;

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
    //centerlinePoints->clear();

    std::cout << "Current threshold " << currentTreshold << std::endl;
    ProcessImageInMemory(currentTreshold);

    GetSkeletonInformation(type, computed, sz);
    for(int i=0;i<nx;++i)
      for(int j=0;j<ny;++j)
      {
          int id1 = j*size + i;
          int id2 = j*fboSize + i;

          if ( currentTreshold == startOfThreshold)
          //if ( tV == 20)
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

    n = 1;

    //std::cout << n << std::endl;
    // Once the image is processed, we have two possible cases that are drawn
    // 1) Is the border information, that is the red(?), on where the area was created
    // 2) The actual centerline info...
    //
    if (n == 0) n = 1; // there is an issue when no centerline is found, and therefore the length is set to -999

    int filled_pixels = 0;
    for(int i=0;i<nx;++i)
      for(int j=0;j<ny;++j)
      {
          int id1 = j*size + i;


          //texture[id1*3 + 0] = aggregated[id1*3 + 0]/n ;
          //texture[id1*3 + 1] = aggregated[id1*3 + 1]/n ;
          //texture[id1*3 + 2] = aggregated[id1*3 + 2]/n ;
          int v1 = aggregated[id1*3 + 0]/n;
          int v2 = aggregated[id1*3 + 1]/n;


          if ( v1 == v2 && v1 > 0){
                   Int2 p;
                   p[0] = i; p[1] = j;
                   //p.p[0] = i;    p.p[1] = j;
                   //centerlinePoints->push_back(p);
                  skeletonPoints->push_back(p);
                  skelImg(i,j) = aggregated[id1*3 + 0]/n ;
          }


          if ( v1 > 0 && v1 == v2 ) filled_pixels++;
      }

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
    int totalSkeletons = 0;

    for (int i = 0; i < imgSize; ++i)															// Generate visualization texture
        for (int j = 0; j < imgSize; ++j)
        {
            int id = j * imgSize + i;															// Position of (i,j) in the various 1D arrays
            //int id2 = j * dataSize + i;
            float p = siteParam[id];
            short r,g,b; r=g=b=0;
             if (skel[id]) totalSkeletons++;
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

    if (show_what==1 && totalSkeletons > 0){
        for(int i=0;i<nendpoints;i++)
        {
           short x = endpoints[2*i], y = endpoints[2*i+1];
           int  id = y * imgSize + x;
           if ( id >= imgSize*imgSize) continue;

           data[id * 3 + 0] = 0;
           data[id * 3 + 1] = 255;
           data[id * 3 + 2] = 0;
        }

    }

}


// Initialization
void SkeletonGenerator::InitializeCUDAMemory(float*& siteParam, short*& outputFT, unsigned char*& outputSkeleton, short*& outputTopo, float*& skelDT, int size)
{
    initialized = true;
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

void SkeletonGenerator::CreateLargestConnectedComponent(short *image, int sz, short* mask){

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




void SkeletonGenerator::ProcessImageInMemory(float threshold) {


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


    skelft2DDT(outputFT,inflation_dist,xm,ym,xM,yM);
    //skelft2DDTshort(outputFT,inflation_dist,xm,ym,xM,yM);

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

    //int skelft2DFill(unsigned char* outputFill, short seedx, short seedy, short xm, short ym, short xM, short yM, unsigned char foreground);

    skelft2DFill((unsigned char*)outputFT,xm+1,ym+1,xm-2,ym-2,xM+2,yM+2,  128);

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
    std::cout << "After topology? " << nendpoints << std::endl;

    std::cout << "A? " << std::endl;
    // Compute DT of in-CUDA-memory skeleton
    // There are two dif. DT, one is the DTB and DTS....
    skel2DSkeletonDT(skelDT,xm,ym,xM,yM);
    // only need the original scatterplot
    // saved into "splot"
    // this can be done before....
    SetOriginalImage(originalImage, originalScatterplot);

    InitializeDisplay(fboSize,outputFT,outputSkeleton,siteParam,
                               nendpoints,outputTopo,skelDT,length,threshold,xm,ym,xM,yM);


    std::cout << "B---" << std::endl;
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


double SkeletonGenerator::GetR(const DataSet& dataSet, Projection2D projection, int sizeT){
   // x, y, valueIndex, bucketIndex, initialIndex, colorIndex points are in this form
    int totalElements = dataSet.size();

    float sumOfClosest = 0;
    for(int j = 0; j < totalElements; j++){
        auto pair = projection(dataSet, j);

        float cx = pair[0];
        float cy = pair[1];

        float closestPointDistance = 9999;

        for(int k = 0; k < totalElements; k++){
            if (j == k ) continue;
            auto pair2 = projection(dataSet, k);

            float ox = pair2[0];
            float oy = pair2[1];
            float d = sqrt(  pow(cx -ox,2.0) + pow(cy -oy,2.0));
            if ( d  < closestPointDistance && d > (2.0/static_cast<float>(sizeT))) // that are not sampled in the same resolution
                closestPointDistance = d;
        }

        sumOfClosest += closestPointDistance;
    }

    sumOfClosest /= static_cast<float>(totalElements);

    return sumOfClosest;
}




void SkeletonGenerator::CreateOnePixelWidthSkeleton(vector<Int2>* inputSkeleton, vector<Int2>* outputSkeleton, int sz){

    // Create an empty image
    // Create an empty image
    typedef itk::Image<unsigned char, 2> ImageType;
    ImageType::Pointer image = ImageType::New();
    ImageType::IndexType start;
    start.Fill(0);

    ImageType::SizeType size;
    size.Fill(sz);

    ImageType::RegionType region(start, size);
    image->SetRegions(region);
    image->Allocate();
    image->FillBuffer(0);
    // load the points
    for(unsigned int i = 0; i < inputSkeleton->size(); i++){

        ImageType::IndexType idx;
        idx[0] = inputSkeleton->at(i)[0];
        idx[1] = inputSkeleton->at(i)[1];
        image->SetPixel(idx, itk::NumericTraits<SPLOMImageType::PixelType>::OneValue() );
    }

    //****************************
    // apply filter
    /*typedef itk::BinaryThinningImageFilter < ImageType, ImageType> BinaryThinningImageFilterType;
    BinaryThinningImageFilterType::Pointer binaryFilter = BinaryThinningImageFilterType::New();
    binaryFilter->SetInput(image);
    binaryFilter->Update();
   // get the points

    ImageType::Pointer thinned = binaryFilter->GetOutput(); //laplacianSharpeningImageFilter->GetOutput();*/

    // The image is thinned, but there are still issues with it...
    // There are certain patterns, corners wise to simplify the issues
    ImageType::SizeType radius;
    radius[0] = 1;
    radius[1] = 1;

    // Possible cases, 1 is filled, 0 is unfilled, and 2 is doesn't matter
    int cases[20][9] = {{0,0,0,1,1,0,0,1,2},   {0,0,0,0,1,1,2,1,0},
                        {0,1,2,1,1,0,0,0,0},   {2,1,0,0,1,1,0,0,0},
                        {2,0,2,1,1,1,0,1,0},   {2,1,0,0,1,1,2,1,0},
                        {0,1,0,1,1,1,2,0,2},   {0,1,2,1,1,0,0,1,2},
                        {0,0,0,0,1,1,0,1,0},   {0,1,0,1,1,0,0,0,0},
                        {0,1,0,0,1,1,0,0,0},   {0,0,0,1,1,0,0,1,0},
                        {0,0,1,0,1,1,0,1,0},   {0,0,0,1,1,0,0,1,1},
                        {0,1,0,1,1,0,1,0,0},   {1,1,0,0,1,1,0,0,0},
                        {0,0,0,0,1,1,1,1,0},   {1,0,0,1,1,0,0,1,0},
                        {0,1,0,0,1,1,0,0,1},   {0,1,1,1,1,0,0,0,0}};

    int specialCases[4][9] = { {1,1,1,1,1,0,1,0,0},
                               {1,1,1,0,1,1,0,0,1},
                               {1,0,0,1,1,0,1,1,1},
                               {0,0,1,0,1,1,1,1,1}};

    int fillInCase[4][2] = { {1,3},{1,5} , {3,7},{5,7}};

    itk::NeighborhoodIterator<ImageType> niterator(radius, image, image->GetLargestPossibleRegion());

    int changes = 0;

    while(!niterator.IsAtEnd()){

        if ( niterator.GetCenterPixel() != 0){
           // Center pixel as something...

            for(int tCase = 0; tCase< 20; tCase++){
                int fulfillsProperty = 0;

                for(int i = 0; i < 9; i++){

                    if ( cases[tCase][i] == 2) {fulfillsProperty++; continue; }

                    bool IsInBounds;
                    int val = niterator.GetPixel(i, IsInBounds);
                    if ( IsInBounds){
                        if ( cases[tCase][i] == val)
                            fulfillsProperty++;
                    }
                }

                if ( fulfillsProperty == 9){
                    changes++;
                    niterator.SetCenterPixel(0);
                }
            }


            for(int sCase = 0; sCase< 4; sCase++){
                int fulfillsProperty = 0;

                for(int i = 0; i < 9; i++){

                    if ( cases[sCase][i] == 2) {fulfillsProperty++; continue; }

                    bool IsInBounds;
                    int val = niterator.GetPixel(i, IsInBounds);
                    if ( IsInBounds){
                        if ( specialCases[sCase][i] == val)
                            fulfillsProperty++;
                    }
                }

                if ( fulfillsProperty == 9){
                    changes++;

                    niterator.SetPixel(fillInCase[sCase][0],0);
                    niterator.SetPixel(fillInCase[sCase][1],0);
                    //niterator.SetCenterPixel(0);
                }
            }

        }
        ++niterator;
    }

    //
    outputSkeleton->clear();
    itk::ImageRegionIteratorWithIndex< ImageType>  it(image, image->GetLargestPossibleRegion());
    while(!it.IsAtEnd()){
        ImageType::IndexType idx = it.GetIndex();

        if ( it.Get() != 0){
           Int2 p;
           p[0] = idx[0];
           p[1] = idx[1];
           outputSkeleton->push_back(p);
        }
        ++it;
    }

}



void SkeletonGenerator::GetEndPointsAndBifurcations(vector<Int2>* centerlinePoints, vector<int>* endPoints, vector<int>* bifurcations){


    vector<int> tmpEndPoints;
    vector<int> tmpBifurcationPoints;

    for(unsigned int i = 0; i < centerlinePoints->size(); i++){

        Int2 c = centerlinePoints->at(i);
        int numNeighbors =0;
        for(unsigned int j = 0; j < centerlinePoints->size(); j++){
            if (i ==j) continue;

            int xDistance = abs(c[0]- centerlinePoints->at(j)[0]);
            int yDistance = abs(c[1]- centerlinePoints->at(j)[1]);

            if ( xDistance <= 1 && yDistance <= 1){      numNeighbors++;   }
            //if (numNeighbors > 1) break;
        }
        if (numNeighbors == 1)
           tmpEndPoints.push_back(i);
        if (numNeighbors > 2)
             tmpBifurcationPoints.push_back(i);
    }

    //if distance between end point is at most one diagonal, change to a single point.
    //
    set<int> remove;
    for(int i = 0; i < tmpEndPoints.size();  i++){
        Int2 c = centerlinePoints->at( tmpEndPoints.at(i));
        for(int j = i+1; j < tmpEndPoints.size();  j++){
            Int2 o = centerlinePoints->at( tmpEndPoints.at(j));

            int xDistance = abs(c[0]- o[0]);
            int yDistance = abs(c[1]- o[1]);
            float totalDistance = sqrt(xDistance*xDistance +  yDistance*yDistance );

            if ( totalDistance < 2.01f) //Merge
            {
                remove.insert(tmpEndPoints.at(j));
            }
        }
    }

    for(int i = 0; i < tmpEndPoints.size(); i++){
        if (  count(remove.begin(), remove.end(),tmpEndPoints.at(i)) == 0 ){
            endPoints->push_back(i);
        }
    }

    //*** do the same for the bifurcations

    remove.clear();
    for(int i = 0; i < tmpBifurcationPoints.size();  i++){
        Int2 c = centerlinePoints->at( tmpBifurcationPoints.at(i));
        for(int j = i+1; j < tmpBifurcationPoints.size();  j++){
            Int2 o = centerlinePoints->at( tmpBifurcationPoints.at(j));

            int xDistance = abs(c[0]- o[0]);
            int yDistance = abs(c[1]- o[1]);
            float totalDistance = sqrt(xDistance*xDistance +  yDistance*yDistance );

            if ( totalDistance < 2.01f) //Merge
            {
                remove.insert(tmpBifurcationPoints.at(j));
            }
        }
    }

    for(int i = 0; i < tmpBifurcationPoints.size(); i++){
        if (  count(remove.begin(), remove.end(),tmpBifurcationPoints.at(i)) == 0 ){
            bifurcations->push_back(i);
        }
    }

}
