#include "randomprojection.hpp"

#include <cmath>

void normalizeProjectionMatrix(ProjectionMatrix2D& matrix) {
    for(int j = 0; j < 2; j++) {
        //normalize by sum (not sum-of-squares)
        float sum = 0;
        for(int i = 0; i < matrix.dimension(); i++) {
            sum += std::abs(matrix.entry(i)[j]) * std::abs(matrix.entry(i)[j]);
        }

        if(sum != 0) {
            for(int i = 0; i < matrix.dimension(); i++) {
                matrix.entry(i)[j] /= sum;
            }
        }
    }
}

ProjectionMatrix2D ProjectionGenerator::generate(const DataSet& dataSet) {
    ProjectionMatrix2D result { dataSet.dimension() };

    for(int j = 0; j < 2; j++) {
        for(int i = 0; i < dataSet.dimension(); i++) {
            float value = 0;
            //skip "flat" dimensions
            if(dataSet.getStat(i).minValue != dataSet.getStat(i).maxValue) {
                //use RNG
                value = distribution(generator);
            }

            result.entry(i)[j] = value;
        }
    }

    normalizeProjectionMatrix(result);

    return result;
}
