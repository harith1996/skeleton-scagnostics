#ifndef STATISTICS_H
#define STATISTICS_H

#include "math.h"
#include "mathhelper.h"
#include "dataset.h"
#include "Reader.h"
#include "ordering.h"
#include <QString>
#include <map>
#include <stdlib.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <utility>
#include <functional>
#include <iostream>
#include <stdio.h>
#include "cuda/projection.h"
// lets use blas for mds

using namespace std;


class Statistics
{
public:
    Statistics();
    ~Statistics();

    void SetDataset(Dataset* data, Reader* mapper);

    bool HasDataset(){ return !dataNotSet;}
    void SetToNoDataset(){ dataNotSet = true; }
    int GetNumberElements(){ return dataset->GetTotalNumberOfElements(); }
    float GetFrequency(int attributeIndex, int valueIndex, bool force = false);
    float GetPDF(int attributeIndex, int valueIndex, int order = 0, bool force = false, float *freq = nullptr);
    float GetMidPointPDF(int attributeIndex, int valueIndex, int order = 0, bool force = false);


    void CreateOrderings();
    void CreateRanges();
    void SetToEqualSpace();
    // Given a probability value get which is the value that it belongs to
    int GetValueIndex(int attributeIndex, float frequency, int order = 0 );

    void ToggleRepresentationToEqualSpaced();
    bool IsRepresentationEquallySpaced(){ return this->equalSpaceRepresentation; }

    int GetTotalOfCategoricalValue(int attributeIndex, int valueIndex);

    float GetElementValue(int elementIndex, int attributeIndex ){ return dataset->GetElementValue(elementIndex, attributeIndex);}

    float GetDistance(vector<float> dataPoint, int j);
    float GetSimilarity(int i, int j, int attributeToIgnore = -1);
    float GetSimilarityFromMem(int i, int j);

    float GetCategoricalSimilarity(int i, int j, int attributeIndex);
    float GetCategoricalSimilarity(vector<float> dataPoint, int j, int attributeIndex);
    float GetCategoricalSimilarity(vector<float> dataPointA, vector<float> dataPointB, int attributeIndex);

    float GetNumericalDistance(vector<float> &dataPoint, int j);

    void CalculateMDSWithPointIndices(vector<int>* indices, vector<pair<float, float> > *normalizedMDSpoints);

    float GetCategoricalSimilarityFromMemory(int i, int j, int attributeIndex);
    double CalculatePearsonCoefficientWithOrdinal(int attrib1, int ordinalAttrib2, bool skip[]);

    float GetNumericalSimilarity(int i, int j);
    void CreateNumericalSimilarity();
    float maxNumericalDistance;

    void CreateCategoricalSimilarity();
    void CreateCategoricalSimilarityForAttribute(int attributeIndex);

    void GetDensityDistributionViaBinning(int attributeIndex, int numberOfBins, vector<double>* bins, bool skip[]);

    void CreateAllSimilarity(){
        CreateNumericalSimilarity();
        CreateCategoricalSimilarity();
        CreateSimilarityMatrix();
    }
    void CreateSimilarityMatrix();


    gsl_matrix* GetCovarianceMatrix(){ return covarianceMatrix; }
    void CreateCPUSimilarityMatrix(int attributeToIgnore = -1);
    void CreatePCAProjectionMatrix(float* prevProjectionMatrix, float *newProjectionMatrix, pair<int, int> eigenVectorsToUse, gsl_matrix* covarianceMatrix = nullptr);

    void CreateLDAProjectionMatrix(float* prevProjectionMatrix, float *newProjectionMatrix, int clusteringAttribute, float regularizationParameter, vector<pair<float, float> > numericalClassesRanges);
    void CreateSkippableLDAProjectionMatrix(float* prevProjectionMatrix, float *newProjectionMatrix, int clusteringAttribute, float regularizationParameter,
                                               vector<pair<float, float> > numericalClassesRanges, bool skip[]);
    void CreateLongProjectionMatrix(float* prevProjectionMatrix, float *newProjectionMatrix, int clusteringAttribute, vector<pair<float, float> > numericalClassesRanges);


    void CreateLaplacianProjectionMatrix(float* prevProjectionMatrix, float *newProjectionMatrix, bool normalize);
    void CreateSkippableLaplacianProjectionMatrix(float* prevProjectionMatrix, float *newProjectionMatrix, bool normalize, bool skip[]);

    void GetTextInfo(vector<QString>* info,int elementIndex);
    void GetHeader(vector<QString>* info);
    float GetMaximumInAttribute(int attributeIndex);
    float GetMinimumInAttribute(int attributeIndex);
    float GetMeanInAttribute(int attributeIndex);
    AttributeRanges GetAttributeRange(int attributeIndex){ return ranges.at(attributeIndex); }

    void SetSimilarityMeasure(int v){ currentSimilarityMeasure = v;}
    int GetSimilarityMeasure(){ return currentSimilarityMeasure; }
    enum SimilarityMeasure { Gower, Lin, Goodall ,Eskin,  Own };

    void SetBestOrdering(bool val){ useBestOrdering = val;}
    int GetBestOrdering(int attributeIndex);
    int GetCustomOrdering(int attributeIndex);
    void SetCustomOrdering(int attributeIndex, int newOrdering);

    int GetPrevValueInOrder(int attributeIndex, int valueIndex, int order);
    int GetNextValueInOrder(int attributeIndex, int valueIndex, int order);
    int GetValueInOrder(int attributeIndex, int valueIndex, int order);

    Attribute GetAttr(int attributeIndex){ return _mapper->GetAttribute(attributeIndex); }
    void SetUseCustomOrdering(bool useCustom){ useCustomOrdering = useCustom; }

    void GetOrder(int orderIndex, int attributeIndex, vector<int>* order );

    int GetDesiredOrder(int desiredOrder[], int attributeIndex);

    void SetElementValue(int elementIndex, int attributeIndex, float newVal );

    void ForceRecalculationOrdering(int attributeIndex);

    void ClearMaps();

    static void CalculateMDS(double *similarityMatrixSquared, int n, vector<pair<float, float> > *normalizedMDSpoints);
    void CalculateMDS(vector<pair<float, float> > *normalizedMDSpoints);


    void CreateCovarianceMatrix();
    void CreateCovarianceMatrix(int attributeToJump, vector<int>* indices, gsl_matrix* covarianceMatrix, bool skip[]);
    void CreateWeightedCovarianceMatrix(int attributeToJump, vector<int>* indices, float weights[], gsl_matrix* covarianceMatrix, bool skip[]);

    AttributeRanges GetNumericalVariableRange(int index);
    AttributeRanges GetNumericalVariableRange(int index, vector<int>* elementIndices, bool skip[]);

    AttributeRanges GetCategoricalVariableRange(int index);
    AttributeRanges GetCategoricalVariableRange(int index, vector<int>* elementIndices, bool skip[]);

    // Dimension similarities now...
    double GetDimensionSimilarity(int attrib1, int attrib2, bool skip[]);
    double CalculateCramersV(int attrib1, int attrib2, bool skip[]);
    double GetCorrelation(int attrib1, int attrib2, bool skip[]);
    double CalculatePearsonCoefficient(int attrib1, int attrib2, bool skip[]);

    int GetTotalOfCategoricalValueSkippable(int attributeIndex, int valueIndex,bool skip[]);
    int GetTotalOfTwoCategoricalValueSkippable(int attributeIndex1, int valueIndex1,int attributeIndex2, int valueIndex2, bool skip[]);

    double SpearmanCorrelation(vector<LocalPoint> *points);

    double MixedNormalizedMutualInformation(int categoricalAttribute, int numericalAttribute, bool skip[]);
    double NumericalNormalizedMutualInformation(int numericalAttribute1, int numericalAttribute2, bool skip[]);
    double CategoricalNormalizedMutualInformation(int categoricalAttribute1, int categoricalAttribute2, bool skip[]);


    double CategoricalEntropy(int attrib, bool skip[]);
    double NumericalEntropy(int attrib, bool skip[]);
    double NumericalKernelEntropy(int attrib, bool skip[], int numberOfBins);

    double MixedMutualInformation(int categoricalAttribute, int numericalAttribute, bool skip[], int k);
    double CategoricalMutualInformation(int categoricalAttribute1, int categoricalAttribute2, bool skip[]);
    double NumericalMutualInformation(int numericalAttribute1, int numericalAttribute2, bool skip[], int k);

    //Kullback Leibler Divergence for categorical variables
    // Dkl (P||Q)... divergence from Q to P
    double CategoricalKullBackDivergence(int attributeIndex, vector<int>& P, vector<int>& Q);
    double NumericalKullBackDivergence(int attributeIndex, vector<int>& P, vector<int>& Q, int numBins);
    //
    double CategoricalKullBackDivergenceGivenValues(int attributeIndex, vector<int>& P, vector<int>& Q);
    double NumericalKullBackDivergenceGivenValues(int attributeIndex, vector<float>& P, vector<float>& Q);

private:


    void PrintMatrix(std::string str, gsl_matrix* A){
        cout << endl << "*********************************" << endl;
        cout << str << endl;
        cout << A->size1 << "x" << A->size2 << endl;
       for(unsigned long i = 0; i < A->size1; i++){
           for(unsigned long j = 0; j < A->size2; j++)
               cout << gsl_matrix_get(A,i,j) << ", ";
           cout << endl;
       }
    }


    bool dataNotSet;

    Dataset* dataset;
    Reader* _mapper;


    gsl_matrix* covarianceMatrix;

    /// Similarity calculations
    /// The unweighted is the similarities in each attribute by itself.
    /// The weighted is already the sum of all into one...
    /// -1 is the numerical
    /// and then the other values are the categorical and ordinal...
    map< const int, float*> attributeSimilarity;

    float* SimilarityMatrix;


    int currentSimilarityMeasure;

    /// Variables and methods for setting the order
    void CalculateBestOrdering(bool force);
    map< pair<int, int>, float> frequencies;
    map< pair<int, int>, float> categoricalValueTotals;

    map< std::string, float> pdfWithOrdering;
    map< std::string, float> midpointpdfWithOrdering;
    // We only need as many orderings as the maximal to minimal sizes
    vector<Ordering> possibleOrderings;

    int* bestOrdering;
    int* customOrdering;
    bool useBestOrdering;
    bool useCustomOrdering;


    bool nonBestOrdering;

    bool equalSpaceRepresentation; // whether ordinal data or categorical use a representation of equal spacing in between

    /// Variables for re-ranging the data
    //vector<float> maximumInAttribute;
    //vector<float> minimumInAttribute;
    //vector<float> meanInAttribute;
    //vector<float> varianceInAttribute;
    vector<AttributeRanges> ranges;
};

#endif // STATISTICS_H
