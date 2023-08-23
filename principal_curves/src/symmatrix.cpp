#include "symmatrix.hpp"

SymMatrix::SymMatrix(size_t n) : data {}, dim { n } {
    data.resize(n * (n + 1) / 2);
    for(int x = 0; x < n; x++) {
        for(int y = x; y < n; y++) {
            (*this)(x, y) = 0;
        }
    }
}

float& SymMatrix::operator()(int x, int y) {
    if(y > x) return (*this)(y, x);
    return data.at(x * (x + 1) / 2 + y);
}

float SymMatrix::operator()(int x, int y) const {
    if(y > x) return (*this)(y, x);
    return data.at(x * (x + 1) / 2 + y);
}

void SymMatrix::resize(size_t n) {
    data.resize(n * (n + 1) / 2);
    for(int i = dim * (dim + 1) / 2; i < data.size(); i++) {
        data[i] = 0;
    }
    dim = n;
}

size_t SymMatrix::dimension() const {
    return dim;
}
