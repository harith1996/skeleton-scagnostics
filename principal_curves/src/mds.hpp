#ifndef MDS_HPP
#define MDS_HPP

#include <vector>
#include <utility>

void CalculateMDS(double *SimilarityMatrix, int n, std::vector<std::pair<float, float> > *normalizedMDSpoints);

void CalculateMDS_Squared(double* similarityMatrixSquared, int n, std::vector<std::pair<float,float>>* normalizedMDSpoints);

#endif

