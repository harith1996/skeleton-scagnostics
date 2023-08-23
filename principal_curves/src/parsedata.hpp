#ifndef PARSEDATA_HPP
#define PARSEDATA_HPP

#include <fstream>

#include "dataset.hpp"
#include "curve2d.hpp"

DataSet parseData(std::istream& infile, char delimiter=',', int multiplier=1, bool debug=false);

DataSet fakeData(const int numPoints,const std::vector<float>& generatingCoefficients, float noise);

DataSet fakeCircle(const int numPoints, float noise);

DataSet fakeX(const int numPoints, float noise);

DataSet fakeSpiral(const int numPoints, float noise);

DataSet fakeTripleCircle(const int numPoints, float noise);

void fakeLinear(const int numPoints, float noise,  DataSet& normal, DataSet& up, DataSet& down);

DataSet growingDNA(const int numPoints, float noise);

DataSet growingDNA2(const int numPoints, float noise);


std::vector<AugmentedPoint> generatingLine(const int numPoints, int sz, DataSet originalDataset, int axis1, int axis2);

std::vector<AugmentedPoint> generatingSquareRoot(const int numPoints, int sz, DataSet originalDataset, int axis1, int axis2);

std::vector<AugmentedPoint> generatingSpiral(const int numPoints, int sz, DataSet originalDataset, int axis1, int axis2);

std::vector<AugmentedPoint> generatingCurve(const int numPoints, const std::vector<float>& generatingCoefficients, int sz, DataSet originalDataset, int axis1, int axis2);

std::vector<AugmentedPoint> generatingCircle(const int numPoints, int sz, DataSet originalDataset, int axis1, int axis2);

std::vector<AugmentedPoint> generatingX(const int numPoints, int sz, DataSet originalDataset, int axis1, int axis2);

std::vector<AugmentedPoint> generatingDNA(const int numPoints, int sz, DataSet originalDataset, int axis1, int axis2);

std::vector<AugmentedPoint> generatingTripleCircle(const int numPoints, int sz,  DataSet originalDataset,int axis1, int axis2);

#endif
