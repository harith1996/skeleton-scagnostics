#include "propagate2d.hpp"

#include <algorithm>

#include "transform2d.hpp"

float propagate2D_diagonal(float a, float b) {
    float d = 2 - (a-b) * (a-b);

    if(d < 0) return INFINITY;

    float r = sqrt(d);
    float solution = ((a + b) - r) / 2;

    float result =  std::max(solution, solution + r);

    //diagonal cannot cause "free moves"
    if(result <= a || result <= b) return INFINITY;

    return result;
}

void propagate2D_updateCell(int x, int y, const Image2D<char>& mask, Image2D<float>& distance, std::multimap<float, Int2>& queue) {
    size_t width = mask.width();
    size_t height = mask.height();

    //nothing to do for non-existing cells
    if(x < 0 || x >= width || y < 0 || y >= height) return;

    //skip cells outside of mask
    if(!mask(x, y)) return;

    float oldDistance = distance(x, y);

    float newDistance = INFINITY;

    //try simple straight moves
    // [maybe need check for finalised]
    if(x >= 1 && distance(x - 1, y) != INFINITY) {
        newDistance = std::min(newDistance, 1 + distance(x - 1, y));
    }
    if(y >= 1 && distance(x, y - 1) != INFINITY) {
        newDistance = std::min(newDistance, 1 + distance(x, y - 1));
    }
    if(x < width - 1 && distance(x + 1, y) != INFINITY) {
        newDistance = std::min(newDistance, 1 + distance(x + 1, y));
    }
    if(y < height - 1 && distance(x, y + 1) != INFINITY) {
        newDistance = std::min(newDistance, 1 + distance(x, y + 1));
    }

    //try diagonal moves
    for(int xd = -1; xd <= 1; xd += 2) {
        for(int yd = -1; yd <= 1; yd += 2) {
            if(x+xd >= 0 && x+xd < width && y+yd >= 0 && y+yd < height
                    && distance(x+xd, y) != INFINITY && distance(x, y+yd) != INFINITY) {
                newDistance = std::min(newDistance, propagate2D_diagonal( distance(x+xd, y), distance(x, y+yd) ));
            }
        }
    }

    //return if no better move to this cell found
    if(newDistance == INFINITY || newDistance >= oldDistance) return;

    //add new entry
    distance(x, y) = newDistance;
    queue.insert({newDistance, {x, y}});
}
#include <iostream>

Image2D<Int2> propagate2D(const Image2D<char>& mask, Image2D<float>& distance) {
    size_t width = mask.width();
    size_t height = mask.height();

    //construct origin image (result)
    auto origin = identityImage(width, height);

    //construct binary image `finalised`
    Image2D<char> finalised { width, height };

    //construct priority queue
    std::multimap<float, Int2> queue;

    //init data... get all points with 0 distance i.e. border or centerline
    for(int x = 0; x < width; x++) {
        for(int y = 0; y < height; y++) {
            if(distance(x, y) == 0) {
                queue.insert( {0, {x, y} } );
            } else {
                distance(x, y) = INFINITY;
            }
        }
    }

    while(!queue.empty()) {
        //pop element
        Int2 position;
        {
            auto begin_iterator = queue.begin();
            position =  begin_iterator->second;
            queue.erase(begin_iterator);
        }
        int x = position[0], y = position[1];

        // basically if the position has been propagated,
        // we set the value to true and not query again
        if(finalised(x, y)) continue;

        //get closest origin from finalised 8-neighbours
        if(distance(x,y) > 0) {
            int minX = x, minY = y;
            float minDistance = INFINITY;

            for(int xd = -1; xd <= 1; xd++) {
                for(int yd = -1; yd <= 1; yd++) {
                    if(x+xd >= 0 && x+xd < width && y+yd >= 0 && y+yd < height
                            && finalised(x+xd, y+yd)) {

                        auto o = origin(x+xd, y+yd);

                        //distance between (x,y) and neighbour's origin (squared)
                        float dist = vectorLengthSquared(o - Int2({x, y}));

                        if(dist < minDistance) {
                            minX = x+xd;
                            minY = y+yd;
                            minDistance = dist;
                        }
                    }
                }
            }
            //Save the label of closest distance
            origin(x, y) = origin(minX, minY);
        }

        //set finalised
        //the distance of this cell is final
        finalised(x, y) = 1;


        //update neighbours... initially the mask is set to 1s completely
        //mask can be changed...modify the distance... if no better
        //distance is found from the current, then do nothing otherwise
        //return ...
        propagate2D_updateCell(x-1, y, mask, distance, queue);
        propagate2D_updateCell(x+1, y, mask, distance, queue);
        propagate2D_updateCell(x, y-1, mask, distance, queue);
        propagate2D_updateCell(x, y+1, mask, distance, queue);
    }

    return origin;
}
