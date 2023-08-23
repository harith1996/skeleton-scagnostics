#ifndef MATHHELPER_H
#define MATHHELPER_H

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <boost/dynamic_bitset.hpp>

#include <math.h>
#include <utility>
#include <functional>
#include <iostream>
#include <algorithm>    // std::sort
#include <iomanip>
#include <locale>
#include <sstream>
#include <string>
#include <bitset>
#include "dataset.h"



using namespace std;

#ifndef MATH_PI
#define MATH_PI 3.14159
#endif


class MathHelper
{
public:
    MathHelper();

    static bool  ArePointsCollinear(float x1[], float x2[], float x3[]);
    static bool  ArePointsCollinear(double x1[],double x2[], double x3[]);

    static float AngleBwVectors( float* v1, float *v2, bool degrees = true);
    static bool  AreEqual(float p1[], float p2[]);

    template <typename T>
      static bool  IsInTriangle(T p[], T a[], T b[], T c[]);

    template <typename T>
      static float Dot(T u[], T v[]);

    template <typename T>
      static void Normalize(T vector[]);

    template <typename T>
      static double Norm(T vector[], int n);

    static double CalculateDistancePointSegment(double *pos1, double *pos2, double *x0);
    static double CalculateTprojectionPointSegment(double *pos1, double *pos2, double *x0);

    //static void Normalize(double vector[]);
    static void Normalize(double vector[], int n);

    static double DistanceBtwPoints(double p1[],double p2[]);

    static int StringDistance(string a, string b);

    static int HammingDistance(string A, string B);


    static int LevDamDist(std::string s1,  std::string s2);

    static double NormalDistribution(double x, double stddev){
          // mean 0 and std dev as set
          double inverse = 1.0 / sqrt(2.0 * MATH_PI * stddev );
          double diff =  pow(x,2.0) / (2.0* stddev * stddev);
          return inverse * exp( -diff );
    }

    //static double Norm(T vector[], int n);

    static string FormattedNumber(double value, int precision);

    static void RotatePointAroundAxis( double angle, double axis[], double point[], double populate[]);
    static void RotationAroundAxis( double angle, gsl_vector* axis, gsl_matrix* R);

    static void SkewSymmetric(const gsl_vector* v,gsl_matrix* Vx);
    static float CosineAngleOfVectors(float v1[], float v2[]);

    static void GetTangent(vector<LocalPoint> *points, double populate[], int position);
    static void GetCurrentPosition(vector<LocalPoint>* points ,double populate[], int position );
    static void IntersectionBetweenTwoLines( std::pair<LocalPoint,LocalPoint> l1,
                                             std::pair<LocalPoint,LocalPoint> l2, double* x ,double *y);


    static double GetDistanceBetweenPoints(LocalPoint a, LocalPoint b);
    static float objectiveFrechet(int method, float* curve1, float* curve2, int numberOfPointsCurve1, int numberOfPointsCurve2, int indexCurve1, int indexCurve2, bool debug = false);

    static float FrechetDistanceBtwTwoLines(float* curve1, float* curve2,
                                            int numberOfPointsCurve1 , int numberOfPointsCurve2, int dx, int dy, int method = 0, bool debug = false);


    static void GetTangentFromCurve(float* curve, float* tangent, int position, int numberSamples);

    static float RateOfChangeAngle(float* curve, int index, int numberOfPoints, bool backwards);
    static void GetCurrentPosition(float* point, int position, float* curve);

    static bool InverseMatrixGSL2(gsl_matrix* A, gsl_matrix* inv, int size){
        // Here, we have to check if the matrix is invertible



        gsl_permutation * p = gsl_permutation_alloc (size);
        int s;
        gsl_linalg_LU_decomp(A, p, &s);

        double det = gsl_linalg_LU_det(A, s);
        std::cout << "Determinant? " << det << std::endl;
        if ( fabs(det) < DBL_MIN*2){
            std::cout << "Non-invertible? " << std::endl;
            return false;
        }
        gsl_linalg_LU_invert(A, p, inv);


        return true;
    }


    static double Gaussian2D(float px, float py,  float qx, float qy, double amplitude, const double stddev2){

        double t1 =  pow(px - qx,2)/ (stddev2);
        double t2 =  pow(py - qy,2)/ (stddev2);
        double gaussianValue = amplitude*exp(-(t1+t2));
        if ( isnan(gaussianValue)){
            std::cout << "Nan " << std::endl;
            return 0;
        }
        return gaussianValue;
    }

   static float HueToRgb(float p, float q, float t);

    static void HSLToRGB(float h, float s, float l, float color[]);

   static double  CalculateDistancePointSegment2(double *pos1, double *pos2, double *x0);

    static double DigammaFunction(double x){
            // en.wikipedia.org/wiki/Digamma_function#computation and approximation  ---
            //
            double v =  log(x) - 0.5*pow(x,-1.0) - (1.0/12.0)*pow(x,-2) + (1.0/120.0)*pow(x,-4);
            double smallV = -(1.0/252.0)*pow(x,-6) + (1.0/240)*pow(x, -8) - (5.0/660)*pow(x, -10) + (691.0/32760)*pow(x, -12) - (1.0/12)*pow(x,-14);
            v += smallV;
            return v;
    }

    static constexpr double EulerMascheroniConstant =  0.57721566;
};

template <typename T>
double MathHelper::Norm(T vector[], int n){
    float size = 0;
    for( int i = 0;i < n; i++) size += vector[i]*vector[i];
    size = sqrt(size);
    return size;
}


template <typename T>
void MathHelper::Normalize(T vector[]){
  double size = Norm(vector, 3);
  for (int i = 0; i < 3; i++) vector[i] /= size;
}

template <typename T>
bool MathHelper::IsInTriangle(T p[], T a[], T b[], T c[])
{
    T v0[3] = {b[0] - a[0], b[1] - a[1], b[2] -a[2]};
    T v1[3] = {c[0] - a[0], c[1] - a[1], c[2] -a[2]};
    T v2[3] = {p[0] - a[0], p[1] - a[1], p[2] -a[2]};

    T d00 = Dot(v0, v0);
    T d01 = Dot(v0, v1);
    T d11 = Dot(v1, v1);
    T d20 = Dot(v2, v0);
    T d21 = Dot(v2, v1);
    T denom = d00 * d11 - d01 * d01;
    T v = (d11 * d20 - d01 * d21) / denom;
    T w = (d00 * d21 - d01 * d20) / denom;
    T u = 1.0f - v - w;

    if ( u < 0.0 || u > 1.0 || v < 0.0 || v > 1.0 || w < 0.0 || w > 1.0){
       return false;
    }
    return true;
}

template <typename T>
float  MathHelper::Dot(T u[], T v[]){
    // Simple Dot Product
    float res = 0.0;
    for(int i = 0; i < 3; i++) res += u[i]*v[i];
    return res;
}

#endif // MATHHELPER_H
