#ifndef TRANSFORM2D_HPP
#define TRANSFORM2D_HPP

#include <algorithm>
#include <utility>
#include <cmath>
#include <vector>

#include "image2d.hpp"
#include "multival.hpp"
#include <iostream>

template<class T>
Image2D<char> distanceThreshold(const Image2D<T>& img, float radius) {
    auto w = img.width();
    auto h = img.height();

    Image2D<char> result { w, h };

    for(int x = 0; x < w; x++) {
        for(int y = 0; y < h; y++) {
            if(img(x, y) != 0) {
                for(int dx = int(-radius); dx <= radius; dx++) {
                    for(int dy = int(-radius); dy <= radius; dy++) {
                        if(dx * dx + dy * dy < radius * radius) {
                            auto cx = x + dx;
                            auto cy = y + dy;
                            if(cx >= 0 && cx < w && cy >= 0 && cy < h) {
                                result(cx, cy) = 1;
                            }
                        }
                    }
                }
            }
        }
    }
    return result;
}

template<class T>
Image2D<float> inverseSquareBlob(const Image2D<T>& img, float radius) {
    auto w = img.width();
    auto h = img.height();

    Image2D<float> result { w, h };
    result.fill(0);

    for(int x = 0; x < w; x++) {
        for(int y = 0; y < h; y++) {
            for(int dx = int(-radius); dx <= radius; dx++) {
                for(int dy = int(-radius); dy <= radius; dy++) {
                    auto cx = x + dx;
                    auto cy = y + dy;
                    if(cx >= 0 && cx < w && cy >= 0 && cy < h) {
                        result(cx, cy) += img(x, y) * 1.0 / std::max(1, dx * dx + dy * dy);
                    }
                }
            }
        }
    }
    return result;
}

template<class T>
Image2D<float> gaussianBlob(const Image2D<T>& img, int sigma) {
    auto w = img.width();
    auto h = img.height();

    int radius = 3 * sigma;

    Image2D<float> kernel { 2 * radius + 1, 2 * radius + 1 };

    for(int x = 0; x < kernel.width(); x++) {
        kernel(x, radius) = (1.0 / sqrt(2 * M_PI * sigma * sigma)) * exp(- (x - radius) * (x - radius) / (2.0 * sigma * sigma));
    }

    for(int x = 0; x < kernel.width(); x++) {
        for(int y = 0; y < kernel.height(); y++) {
            if(y != radius) {
                kernel(x, y) = kernel(x, radius) * kernel(y, radius);
            }
        }
    }
    float q = kernel(radius, radius);
    for(int x = 0; x < kernel.width(); x++) {
        kernel(x, radius) *= q;
    }

    Image2D<float> result { w, h };
    result.fill(0);

    for(int x = 0; x < w; x++) {
        for(int y = 0; y < h; y++) {
            for(int dx = -radius; dx <= radius; dx++) {
                for(int dy = -radius; dy <= radius; dy++) {
                    auto cx = x + dx;
                    auto cy = y + dy;
                    if(cx >= 0 && cx < w && cy >= 0 && cy < h) {
                        result(cx, cy) += kernel(dx + radius, dy + radius) * img(x, y);
                    }

                }
            }
        }
    }

    return result;
}

template<class T>
Image2D<char> threshold(const Image2D<T>& img, T t, double* percentage) {
    auto w = img.width();
    auto h = img.height();

    int tot = 0;
    int comp = 0;
    Image2D<char> result { w, h };

    for(int x = 0; x < w; x++) {
        for(int y = 0; y < h; y++) {

            if ( x > w*0.25 && x < w*0.75){
                if ( y > h*0.25 && y < h*0.75)
                    comp++;
            }

            if(img(x, y) >= t) {
                result(x, y) = 1;
                if ( x > w*0.25 && x < w*0.75){
                    if ( y > h*0.25 && y < h*0.75)
                        tot++;
                }
            }
        }
    }

   *percentage = static_cast<float>(tot)/comp;
    return result;
}

template<class T, class S>
Image2D<T> boundary(const Image2D<S>& img, T valueBoundary, T valueOther) {
    auto w = img.width();
    auto h = img.height();
    Image2D<T> result { w, h };

    //boundary is outside of 1_img
    //but 4-connected to it

    for(int x = 0; x < w; x++) {
        for(int y = 0; y < h; y++) {
            result(x, y) = valueOther;

            if(img(x, y)) {
                if(x == 0 || y == 0 || x == w - 1 || y == h - 1) {
                    result(x, y) = valueBoundary;
                }
            } else {

                for(int x2 = std::max(x-1, 0); x2 <= x+1 && x2 < w; x2++) {
                    for(int y2 = std::max(y-1, 0); y2 <= y+1 && y2 < h; y2++) {

                        //check for 4-neighbourhood
                        auto xd = x2 - x;
                        auto yd = y2 - y;

                        if( ((xd == 0) ^ (yd == 0)) && img(x2, y2)) {
                            result(x, y) = valueBoundary;
                        }
                    }
                }
            }
        }
    }

    return result;
}

template<class T>
void arcLengthLabel_checkNeighbour(Int2 pos, Int2 diff, std::vector<Int2>& stack, const Image2D<T>& img, T marker, Image2D<float>& result) {
    int w = int(img.width());
    int h = int(img.height());

    auto newPos = pos + diff;

    if(newPos[0] < 0 || newPos[1] >= w) return;
    if(newPos[1] < 0 || newPos[1] >= h) return;

    if(img(newPos[0], newPos[1]) == marker && result(newPos[0], newPos[1]) == INFINITY) {
        result(newPos[0], newPos[1]) = result(pos[0], pos[1]) + sqrt(vectorLengthSquared(diff));
        stack.push_back(newPos);
    }
}

template<class T>
Image2D<float> arcLengthLabel(const Image2D<T>& img, T marker, float& maxLength) {
    int w = int(img.width());
    int h = int(img.height());

    Image2D<float> result { static_cast<size_t>(w), static_cast<size_t>(h) };
    result.fill(INFINITY);

    maxLength = 0;

    std::vector<Int2> stack;

    //find a marker. The marker in this case is the border ....
    int reaches = 0;
    for(int x = 0; x < w; x++) {
        for(int y = 0; y < h; y++) {
            //skip non-markers
            if(img(x, y) != marker) continue;
            reaches++;
            //skip already visited
            if(result(x, y) != INFINITY) continue;


            stack.emplace_back(x, y);
            result(x, y) = 0;

            //while label not already set
            do {
                auto pos = stack.back();
                stack.pop_back();

                //update maxLength
                maxLength = std::max(result(pos[0], pos[1]), maxLength);

                //find all marker neighbours where result is infinity

                arcLengthLabel_checkNeighbour(pos, {+1,  0}, stack, img, marker, result);
                arcLengthLabel_checkNeighbour(pos, { 0, +1}, stack, img, marker, result);
                arcLengthLabel_checkNeighbour(pos, {-1,  0}, stack, img, marker, result);
                arcLengthLabel_checkNeighbour(pos, { 0, -1}, stack, img, marker, result);
                arcLengthLabel_checkNeighbour(pos, {+1, +1}, stack, img, marker, result);
                arcLengthLabel_checkNeighbour(pos, {+1, -1}, stack, img, marker, result);
                arcLengthLabel_checkNeighbour(pos, {-1, +1}, stack, img, marker, result);
                arcLengthLabel_checkNeighbour(pos, {-1, -1}, stack, img, marker, result);
            } while(!stack.empty());
        }
    }

    //std::cout << "Markers? " << reaches <<std::endl;
    return result;
}

Image2D<Int2> identityImage(size_t width, size_t height);

template<class T>
Image2D<T> indexImage(const Image2D<T>& img, const Image2D<Int2>& index) {
    auto w = index.width();
    auto h = index.height();

    Image2D<T> result { w, h };

    for(int x = 0; x < w; x++) {
        for(int y = 0; y < h; y++) {
            auto ind = index(x, y);
            result(x, y) = img(ind[0], ind[1]);
        }
    }
    return result;
}

template<class T>
Image2D<float> normalizeImageValues(const Image2D<T>& img) {
    auto w = img.width();
    auto h = img.height();

    Image2D<float> result { w, h };

    if(w * h == 0) return result;

    auto max = img(0, 0);

    for(int x = 0; x < w; x++) {
        for(int y = 0; y < h; y++) {
            max = std::max(max, img(x, y));
        }
    }
    for(int x = 0; x < w; x++) {
        for(int y = 0; y < h; y++) {
            result(x, y) = img(x, y) * 1.0 / max;
        }
    }

    return result;
}

#endif
