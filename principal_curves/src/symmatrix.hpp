#ifndef SYMMATRIX_HPP
#define SYMMATRIX_HPP

#include <vector>

class SymMatrix {
    std::vector<float> data;
    size_t dim;

    public:

    SymMatrix(size_t n);

    float& operator()(int x, int y);
    float operator()(int x, int y) const;

    size_t dimension() const;

    void resize(size_t n);
};


#endif
