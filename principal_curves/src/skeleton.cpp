#include "skeleton.hpp"

#include <cmath>
//#include <cstdio>
#include <algorithm>

#include "skeleton/mfmm.h"
#include "skeleton/flags.h"

#include "image2d.hpp"
#include "transform2d.hpp"
#include "propagate2d.hpp"

inline float distance(int cd, float length)
//Computes distance along boundary in a wrap-around fashion...
{ return min(fabs(float(cd)),fabs(length-fabs(float(cd)))); }

Image2D<char> computeSkeleton(const Image2D<char>& img) {
    int w = int(img.width());
    int h = int(img.height());

    FIELD<float> field { w, h };
    //copy data
    for(int x = 0; x < w; x++) {
        for(int y = 0; y < h; y++) {
            field(x, y) = img(x, y) ? 0 : INFINITY;
        }
    }

    //compute flags
    float k = -1;
    FLAGS flags { field, k };

    //alloc other data
    FIELD<float> count;
    FLAGS flags_copy { flags };

    //execute
    ModifiedFastMarchingMethod fmm { &field, &flags_copy, &count, /*origs*/ nullptr };

    fmm.setScanDir(0);											//Set its various parameters
    fmm.setMethod(ModifiedFastMarchingMethod::AFMM_STAR);
    fmm.setBoundaryLabeling(ModifiedFastMarchingMethod::ARC_LENGTH);

    int nfail,nextr;
    int iter = fmm.execute(nfail,nextr);
    float length = fmm.getLength();

    /*
    //print field
    for(int x = 0; x < w; x++) {
        for(int y = 0; y < h; y++) {
            printf("%1.1f ", field(x, y));
        }
        printf("\n");
    }
    */

    //compute grad
    FIELD<float> grad { w, h };

    int i,j;
    for(j=0;j<grad.dimY();j++)              //Compute grad in a special way, i.e. on a 2-pixel
        for(i=0;i<grad.dimX();i++)          //neighbourhood - this ensures pixel-size skeletons!
        {
            float ux = count.value(i+1,j) - count.value(i,j);
            float uy = count.value(i,j+1) - count.value(i,j);
            grad.value(i,j) = flags.faraway(i,j)? max(distance(ux, length),distance(uy, length)) : 0;
        }

    //copy to result
    float level = 1;
    Image2D<char> result { static_cast<size_t>(w), static_cast<size_t>(h) };

    for(int x = 0; x < w; x++) {
        for(int y = 0; y < h; y++) {
            auto g = grad(x, y);
            result(x, y) = ((g == INFINITY) || (g > level));
        }
    }

    return result;
}

Image2D<float> computeSkeleton2_Grad(const Image2D<char>& img, float* maxDt) {

    //  Setting the distance 0 in border locations.. everywhere else infinity;
    Image2D<float> distance = boundary(img, 0.0f, INFINITY);

    auto w = img.width();
    auto h = img.height();

    float maxLength = 0;
    auto labelOrigin = arcLengthLabel(distance, 0.0f, maxLength);

    //std::cout << std::endl;
    //std::cout << "max length " << maxLength << std::endl;
    /*----- TODO: one of this two methods is not working greatly */
    auto origin = propagate2D(img, distance);
    auto label = indexImage(labelOrigin, origin);

    /*----- Fix THIS*/
    for(int x = 0; x < w; x++) {
        for(int y = 0; y < h; y++) {

            //?inf label...  then it makes label -1,-1
            // std::cout << label(x,y) << "  ";
            if( label(x,y) == INFINITY) distance(x,y) = 0;
        }
    }

    //std::cout << std::endl << std::endl;
    //std::cout << "Result " << std::endl;

    Image2D<float> result {w, h};
    float mx = 0;

    for(int x = 0; x < w; x++) {
        for(int y = 0; y < h; y++) {
            if(distance(x, y) == 0) continue;
            auto diffX = std::abs(label(x, y) - label(x+1, y));
            diffX = std::min(diffX, maxLength - diffX);
            auto diffY = std::abs(label(x, y) - label(x, y+1));
            diffY = std::min(diffY, maxLength - diffY);
            // diffX and diffY get to be nan..
            //if(diffX == INFINITY) diffX = 0;//Previously commented, why is this happening (?)
            //if(diffY == INFINITY) diffY = 0;//
            result(x, y) = std::max(diffX, diffY);
            mx = std::max(result(x,y), mx);
        }
    }

    //std::cout << std::endl;
    *maxDt = mx;

    return result;
}
