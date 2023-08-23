#ifndef SKELETON_HPP
#define SKELETON_HPP

#include "image2d.hpp"

Image2D<char> computeSkeleton(const Image2D<char>& img);

Image2D<float> computeSkeleton2_Grad(const Image2D<char>& img, float *maxDt);

#endif
