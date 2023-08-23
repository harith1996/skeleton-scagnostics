#ifndef IMAGE2DIO_HPP
#define IMAGE2DIO_HPP

#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>

#include "image2d.hpp"

template<class T>
void writePBM(std::ostream& stream, const Image2D<T>& img) {
    auto w = img.width();
    auto h = img.height();

    stream << "P1\n" << w << " " << h << "\n";

    for(int j = 0; j < img.height(); j++) {
        for(int i = 0; i < img.width(); i++) {
            stream << (img(i, j) != 0) << " ";
        }
        stream << "\n";
    }
}

template<class T>
void writePBM(std::string fileName, const Image2D<T>& img) {
    std::ofstream file { fileName };
    writePBM(file, img);
}

template<class T>
void writePGM(std::ostream& stream, const Image2D<T>& img, T min, T max) {
    auto w = img.width();
    auto h = img.height();

    stream << "P2\n" << w << " " << h << " 255\n";

    for(int j = 0; j < img.height(); j++) {
        for(int i = 0; i < img.width(); i++) {
            auto val = img(i, j);
            val = std::min(max, std::max(min, val));
            stream << static_cast<int>(((val - min) * 255 / (max - min))) << " ";
        }
        stream << "\n";
    }
}

template<class T>
void writePGM(std::string fileName, const Image2D<T>& img, T min, T max) {
    std::ofstream file { fileName };
    writePGM(file, img, min, max);
}


#endif
