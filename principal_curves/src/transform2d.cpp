
#include "transform2d.hpp"

Image2D<Int2> identityImage(size_t width, size_t height) {
    Image2D<Int2> result { width, height };

    for(int x = 0; x < width; x++) {
        for(int y = 0; y < height; y++) {
            result(x, y) = {int(x), int(y)};
        }
    }

    return result;
}
