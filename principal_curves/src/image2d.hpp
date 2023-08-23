#ifndef IMAGE2D_HPP
#define IMAGE2D_HPP

#include <vector>
#include <cstring>

template<class T>
class Image2D {
    std::vector<T> buffer;
    size_t _width, _height;

    public:

    Image2D(size_t w, size_t h) : buffer { }, _width { w }, _height { h } {
        buffer.resize(w * h);
    }

    T& operator()(size_t x, size_t y) {
        return buffer[x + y * _width];
    }

    const T& operator()(size_t x, size_t y) const {
        return buffer[x + y * _width];
    }

    size_t width() const {
        return _width;
    }

    size_t height() const {
        return _height;
    }

    void fill(const T& t) {
        for(int x = 0; x < _width; x++) {
            for(int y = 0; y < _height; y++) {
                (*this)(x, y) = t;
            }
        }
    }
};



#endif
