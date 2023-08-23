#ifndef PROPAGATE2D_HPP
#define PROPAGATE2D_HPP

#include <map>
#include <cstdio>
#include <cmath>

#include "image2d.hpp"
#include "multival.hpp"


Image2D<Int2> propagate2D(const Image2D<char>& mask, Image2D<float>& distance);

#endif
