#include "mathhelper.h"

MathHelper::MathHelper()
{

}



bool MathHelper::ArePointsCollinear(float x1[], float x2[], float x3[]){

    double triangleArea = x1[0]*( x2[1] - x3[1]) +  x2[0]*(x3[1] - x1[1]) + x3[0]*(x1[1] - x2[1]);
    return ( fabs(triangleArea) < 0.001);
}


bool MathHelper::ArePointsCollinear(double x1[], double x2[], double x3[]){

    double triangleArea = x1[0]*( x2[1] - x3[1]) +  x2[0]*(x3[1] - x1[1]) + x3[0]*(x1[1] - x2[1]);
    return ( fabs(triangleArea) < 0.001);
}

int MathHelper::StringDistance(string a, string b){

    return HammingDistance(a,b);//LevDamDist(a,b);
}

int MathHelper::LevDamDist(std::string s1,  std::string s2)
{
  size_t size1 = s1.size();
  size_t size2 = s2.size();
  size_t d[size1 + 1][size2 + 1];
  for (int i = 0; i <= size1; i ++)
    d[i][0] = i;
  for (int i = 0; i <= size2; i ++)
    d[0][i] = i;

  int cost = 0;
  for (int i = 1; i <= size1; i ++)
    for (int j = 1; j <= size2; j ++)
    {
      cost = (s1[i - 1] == s2[j - 1]) ? 0 : 1 ;
      if ( (i > 1) && (j > 1) && (s1[i] == s2[j - 1]) && (s1[i - 1] == s2[j]))
      {
        size_t a = std::min(d[i - 1][j], d[i][j - 1] + 1);
        size_t b = std::min(d[i][j] + cost, d[i - 2][j - 2]);
        d[i][j] = std::min(a, b);
      }
      else
      {
        d[i][j] = std::min(std::min(d[i][j -1] + 1, d[i - 1][j] + 1), d[i - 1][j - 1] + cost);
      }
    }
  return d[size1][size2];
}


int MathHelper::HammingDistance(string A, string B){
    const int N = A.length();

    boost::dynamic_bitset<> a(N);
    boost::dynamic_bitset<> b(N);
    for(int i = 0 ; i < N; i++){
        char nowA = A.at(i);
        char nowB = B.at(i);

        a[i] = (nowA == '1');
        b[i] = (nowB == '1');

    }


    auto C = a^b; // if binary, then we xor the value and count the bits set
    return C.count();
}

string MathHelper::FormattedNumber(double value, int precision){
    ostringstream Convert;


    Convert.imbue(locale());       // Imbue the custom locale to the stringstream

    Convert << fixed << setprecision(precision) << value; // Use some manipulators

    string Result = Convert.str();

    return Result;

}

double MathHelper::GetDistanceBetweenPoints(LocalPoint a, LocalPoint b){

    double d = sqrt( (a.p[0] - b.p[0])*(a.p[0] - b.p[0]) + (a.p[1] - b.p[1])*(a.p[1] - b.p[1])) ;
    return d;
}


float MathHelper::CosineAngleOfVectors(float v1[], float v2[]){

    MathHelper::Normalize(v1);
    MathHelper::Normalize(v2);

    gsl_vector *a = gsl_vector_alloc(3);
    gsl_vector *b = gsl_vector_alloc(3);
    //MathHelper::Normalize(v1);
    //MathHelper::Normalize(v2);
    for( int i = 0; i < 3; i++){
       gsl_vector_set(a, i, v1[i]);
       gsl_vector_set(b, i, v2[i]);
    }

    double c = 0;
    gsl_blas_ddot(a,b, &c);
    gsl_vector_free(a);
    gsl_vector_free(b);
    return c;
}
float MathHelper::AngleBwVectors( float* v1, float *v2, bool degrees){

   gsl_vector *a = gsl_vector_alloc(3);
   gsl_vector *b = gsl_vector_alloc(3);
   //MathHelper::Normalize(v1);
   //MathHelper::Normalize(v2);
   for( int i = 0; i < 3; i++){
      gsl_vector_set(a, i, v1[i]);
      gsl_vector_set(b, i, v2[i]);
   }

   double c = 0;
   gsl_blas_ddot(a,b, &c);
   gsl_vector_free(a);
   gsl_vector_free(b);

   // the dot product is 0....

   float origin[3] = {0,0,0};

   float angle =  acos(c);

   if ( ArePointsCollinear(origin, v1,v2)){

       // if the points are collinear, then either they are
       // going in the same direction and the angle is then 0
       // or they are not equal and they are going in different directions....
      if ( AreEqual(v1,v2))
          angle = 0;
      else
          angle = 3.1415;

   }
   // if the two values are really close, then the result will be 1 and epsilon higher
   // so, if isnan then set to 0 the angle.
   else if ( std::isnan(angle)){
       /*std::cout <<" angle is nan " <<  std::endl;
       std::cout << "v1:: " << v1[0] << ", " << v1[1] << ", "  << v1[2] << std::endl;
       std::cout << "v2:: " << v2[0] << ", " << v2[1] << ", "  << v2[2] << std::endl;
       */
       angle = 0;
   }
   if ( degrees )
       angle = angle * (180.0 / 3.1415);
   return angle;
}

bool  MathHelper::AreEqual(float p1[], float p2[]){
    double error = 0;
    for(int i = 0; i < 3; i++) error += fabs( p1[i] - p2[i]);

    return (error < 0.00001);
}





double MathHelper::CalculateDistancePointSegment(double *pos1, double *pos2, double *x0){

    double t = MathHelper::CalculateTprojectionPointSegment(pos1, pos2, x0);


    if ( t < 0 || t > 1){
        return 999999999999;
    }


    double dv[3] = { pos2[0] - pos1[0], pos2[1] -pos1[1], pos2[2] -pos1[2]};

    //? Should I really normalize this...
    //MathHelper::Normalize(dv);

    double projectedPoint[3];
    for(int i = 0; i < 3; i++) projectedPoint[i] = pos1[i] + t*dv[i];


    return MathHelper::DistanceBtwPoints(projectedPoint, x0);
}

double MathHelper::CalculateTprojectionPointSegment(double *pos1, double *pos2, double *x0)
{
    double t = -((pos2[0] - pos1[0])*(pos1[0] - x0[0]) +
        (pos2[1] - pos1[1])*(pos1[1] - x0[1]) +
        (pos2[2] - pos1[2])*(pos1[2] - x0[2])) /
       ((pos2[0] - pos1[0])*(pos2[0] - pos1[0]) +
        (pos2[1] - pos1[1])*(pos2[1] - pos1[1]) +
        (pos2[2] - pos1[2])*(pos2[2] - pos1[2]));

   return t;
}

/*void MathHelper::Normalize(double vector[]){
   Normalize(vector, 3);
}*/

void MathHelper::Normalize(double vector[], int n){
  float size = Norm(vector, n);
  for (int i = 0; i < n; i++) vector[i] /= size;
}

double MathHelper::DistanceBtwPoints(double p1[],double p2[]){
    double dv[3] = {p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]};
    return Norm(dv, 3);
}



void MathHelper::RotatePointAroundAxis( double angle, double axis[], double point[], double populate[])
{
    MathHelper::Normalize(axis);// in case it is not normalized.

    gsl_matrix* R = gsl_matrix_alloc(3, 3);
    gsl_vector *a = gsl_vector_alloc(3);
    gsl_vector *p = gsl_vector_alloc(3);

    for(int i = 0; i < 3 ; i++){
       gsl_vector_set(a, i, axis[i]);
       gsl_vector_set(p, i, point[i]);
    }

    MathHelper::RotationAroundAxis( angle, a, R);

    gsl_vector *result = gsl_vector_alloc(3);

    gsl_blas_dgemv(CblasNoTrans, 1.0, R, p, 0, result);

    for(int i = 0; i < 3; i++)
      populate[i] = gsl_vector_get(result,i);

    gsl_matrix_free(R);
    gsl_vector_free(a);
    gsl_vector_free(p);
    gsl_vector_free(result);
}

void MathHelper::RotationAroundAxis( double angle, gsl_vector* axis, gsl_matrix* R){
    // Create the rotation matrix according to Euler-Rodrigues formula.
    double radians = angle * (MATH_PI/180.0);
    double cosine = cos(radians);
    double sine = sin(radians);
    int degree = 3;

    gsl_matrix_set_zero(R);

    gsl_matrix* I = gsl_matrix_alloc(degree, degree);
    gsl_matrix_set_identity(I);


    gsl_matrix* CrossU = gsl_matrix_alloc(degree, degree);
    gsl_matrix* TensorU = gsl_matrix_alloc(degree, degree);

    MathHelper::SkewSymmetric(axis,  CrossU);

    gsl_matrix_set_zero(TensorU);
    for( int i = 0; i < degree ; i++){
      for(int j = 0; j < degree; j++){
          gsl_matrix_set(TensorU, i, j, gsl_vector_get(axis, i)* gsl_vector_get(axis,j) );
       }
    }

   gsl_matrix_scale(I, cosine);
   gsl_matrix_scale(CrossU, sine  );
   gsl_matrix_scale(TensorU, (1.0 - cosine));


   gsl_matrix_add(R, I);
   gsl_matrix_add(R, CrossU);
   gsl_matrix_add(R, TensorU);

   gsl_matrix_free(I);
   gsl_matrix_free(CrossU);
   gsl_matrix_free(TensorU);
}

void MathHelper::SkewSymmetric(const gsl_vector* v, gsl_matrix* Vx){
   gsl_matrix_set_zero(Vx);

   gsl_matrix_set(Vx, 0, 1, -gsl_vector_get(v, 2) );
   gsl_matrix_set(Vx, 1, 0, gsl_vector_get(v, 2) );

   gsl_matrix_set(Vx, 0, 2, gsl_vector_get(v, 1) );
   gsl_matrix_set(Vx, 2, 0, -gsl_vector_get(v, 1) );


   gsl_matrix_set(Vx, 1, 2, -gsl_vector_get(v, 0) );
   gsl_matrix_set(Vx, 2, 1, gsl_vector_get(v, 0) );
}

void MathHelper::GetCurrentPosition(vector<LocalPoint>* points ,double populate[], int position ){
    populate[0] = points->at(position).p[0];
    populate[1] = points->at(position).p[1];
    populate[2] = 0;
}

void MathHelper::GetTangent(vector<LocalPoint>* points,double populate[], int position){
    //Lets define the tangent depending on the amount of points available.
    // In the case of being 2 <= i <= n -2 we use a method to get an error of O(h^5)
    // Otherwise we use only the difference.
    int numberSamples = points->size();
    if (numberSamples <= 1) { populate[0] = 0; populate[1] = 0; populate[2] = 0;  return;}
    if ( numberSamples < 5){
        // forward difference   f(x+h) - f(x)
        // backward difference  f(x) -  f(x-h)
        // central difference   (f(x+h)-f(x-h))/2.0
        // if number samples < 2, check forward or backward integration
        // well, actually, if they are in the borders...

       if (position == 0){
           // forward
           double p_ph[3], p[3];
           GetCurrentPosition(points, p_ph, position+1);
           GetCurrentPosition(points, p   , position);
           for(int i  =0 ; i < 3; i++){ populate[i]= p_ph[i] - p[i]; }
           MathHelper::Normalize( populate);
       }
       else if ( position == numberSamples -1){
           // backward
           double p[3], p_mh[3];
           GetCurrentPosition(points, p   , position);
           GetCurrentPosition(points, p_mh, position -1);
           for(int i  =0 ; i < 3; i++){ populate[i]= p[i] - p_mh[i]; }
           MathHelper::Normalize( populate);
       }
       else {
           // central diff.
           double p_mh[3], p_ph[3];
           GetCurrentPosition(points, p_ph, position +1);
           GetCurrentPosition(points, p_mh, position -1);
           for(int i  =0 ; i < 3; i++){ populate[i]= (p_ph[i] - p_mh[i])/2.0; }
           MathHelper::Normalize( populate);
       }

    }
    else if ( position >= 2 && position < numberSamples - 2){
        double x_minus2[3], x_minus1[3], x_plus1[3], x_plus2[3];

        GetCurrentPosition(points, x_minus2, position -2);
        GetCurrentPosition(points, x_minus1, position -1);
        GetCurrentPosition(points, x_plus1, position + 1);
        GetCurrentPosition(points, x_plus2, position + 2);
        for(int i  =0 ; i < 3; i++){ populate[i]= x_minus2[i] - 8.0f*x_minus1[i] + 8.0f*x_plus1[i] - x_plus2[i]; }
        MathHelper::Normalize( populate);
    }
    else {
      // boundary conditions.
      if ( position < 0 ) position = 0;

      if ( position == 0){
          double x0[3], x1[3], x2[3], x3[3], x4[3];
          GetCurrentPosition(points, x0, 0);
          GetCurrentPosition(points, x1, 1);
          GetCurrentPosition(points, x2, 2);
          GetCurrentPosition(points, x3, 3);
          GetCurrentPosition(points, x4, 4);
          for(int i  =0 ; i < 3; i++){ populate[i]=   -25.0*x0[i]   + 48.0*x1[i] - 36.0*x2[i]  + 16*x3[i] - 3.0*x4[i]; }
          MathHelper::Normalize( populate);
      }
      else if ( position == 1){
          double x0[3], x1[3], x2[3], x3[3], x4[3];
          GetCurrentPosition(points, x0, 0);
          GetCurrentPosition(points, x1, 1);
          GetCurrentPosition(points, x2, 2);
          GetCurrentPosition(points, x3, 3);
          GetCurrentPosition(points, x4, 4);
          for(int i  =0 ; i < 3; i++){ populate[i]=   -3.0*x0[i]   -10.0*x1[i] +18.0*x2[i]  - 6*x3[i] + x4[i]; }
          MathHelper::Normalize( populate);
      }
      else if ( position == numberSamples - 2){ // tn-1
          double x0[3], x1[3], x2[3], x3[3], x4[3];
          GetCurrentPosition(points, x0, numberSamples - 1);
          GetCurrentPosition(points, x1, numberSamples - 2);
          GetCurrentPosition(points, x2, numberSamples - 3);
          GetCurrentPosition(points, x3, numberSamples - 4);
          GetCurrentPosition(points, x4, numberSamples - 5);
          for(int i  =0 ; i < 3; i++){ populate[i]=   3.0*x0[i]  + 10.0*x1[i] - 18.0*x2[i]  + 6*x3[i] - x4[i]; }
          MathHelper::Normalize( populate);
      }
      else {
          // position == numberSamples -1 or above.
          double x0[3], x1[3], x2[3], x3[3], x4[3];
          GetCurrentPosition(points, x0, numberSamples - 1);
          GetCurrentPosition(points, x1, numberSamples - 2);
          GetCurrentPosition(points, x2, numberSamples - 3);
          GetCurrentPosition(points, x3, numberSamples - 4);
          GetCurrentPosition(points, x4, numberSamples - 5);
          for(int i  =0 ; i < 3; i++){ populate[i]=   25.0*x0[i]  - 48.0*x1[i] + 36.0*x2[i]  - 16.0*x3[i] + 3.0f*x4[i]; }
          MathHelper::Normalize( populate);
      }
    }
}


void MathHelper::IntersectionBetweenTwoLines(std::pair<LocalPoint, LocalPoint> l1, std::pair<LocalPoint, LocalPoint> l2, double* x , double *y){
    //Formula taken from http://mathworld.wolfram.com/Line-LineIntersection.html
    float denS1 = l1.first.p[0] - l1.second.p[0];// x1 - x2
    float denS2 = l1.first.p[1] - l1.second.p[1];// y1 - y2
    float denS3 = l2.first.p[0] - l2.second.p[0];// x3 - x4
    float denS4 = l2.first.p[1] - l2.second.p[1];// y3 - y4

    float den = denS1*denS4 - denS2*denS3;
    // The same denominator for both

    // | x1 y1 |   | x3  y3 |
    // | x2 y2 |   | x4  y4 |
    float numS1 = l1.first.p[0]*l1.second.p[1]  - l1.second.p[0]*l1.first.p[1];
    float numS2 = l2.first.p[0]*l2.second.p[1]  - l2.second.p[0]*l2.first.p[1];

    float numX = numS1*(l2.first.p[0] - l2.second.p[0]) -numS2*(l1.first.p[0] - l1.second.p[0]);
    float numY = numS1*(l2.first.p[1] - l2.second.p[1]) -numS2*(l1.first.p[1] - l1.second.p[1]);

    double newX = numX / den;
    double newY = numY / den;
    *x = newX;
    *y = newY;
}




float MathHelper::FrechetDistanceBtwTwoLines(float* curve1, float* curve2, int numberOfPointsCurve1, int numberOfPointsCurve2,
                                             int dx, int dy, int method, bool debug){

    //float* frechetDistancesArray = (float*) malloc( sizeof(float) * numberOfPointsCurve1 *numberOfPointsCurve2);
    float frechetDistancesArray[numberOfPointsCurve1*numberOfPointsCurve2];
    for(int i = 0; i < numberOfPointsCurve1*numberOfPointsCurve2; i++)
        frechetDistancesArray[i] = -2;
    // i is for the projected
    // j is for the original centerline
    int ti, tj;
    ti = 0; tj = 0;
    // Frechet distance is variant to euclidean transformations ... variant to rotation is good
    // Scaling, somewhat...
    // Translation is what we need to take care...


    double tx = curve1[0] - curve2[0];
    double ty = curve1[1] - curve2[1];

    for(int i = 0; i < numberOfPointsCurve2; i++){
        curve2[i*3 + 0] += tx;
        curve2[i*3 + 1] += ty;
    }

    // This is the value that we need to calculate  the maximazing and minimizing properties....
    // float current = sqrt( (p_q[0] -xyz[0])*(p_q[0]-xyz[0]) +  (p_q[1] -xyz[1])*(p_q[1]-xyz[1]) + (p_q[2] -xyz[2])*(p_q[2]-xyz[2]) );  /** 1 change **/
    float current = objectiveFrechet(method, curve1, curve2, numberOfPointsCurve1, numberOfPointsCurve2,  ti, tj , debug);

    if (debug) std::cout << "Debug starting points " << current << std::endl;
    // The first position, 0,0 is the distance between the starting points of the two curves....
    frechetDistancesArray[0] = current;

    // Now we can do the first row and the first column of the frechet distance array
    // ti = row,  tj = col
    // we start for each row

    if ( debug ) std::cout << "Filling row " << std::endl;
    for(ti = 1; ti < numberOfPointsCurve1; ti++){

          // get the point in the projection for comparison for with the first centreline point
          // get the distance from the first centerline point and the projection
          current = objectiveFrechet(method, curve1, curve2, numberOfPointsCurve1, numberOfPointsCurve2, ti, tj );
          // now let's set the max to the location, comparing previous and current.
          //frechetDistancesArray[ti*numberOfPointsCurve1 + tj] = fmaxf(current, frechetDistancesArray[(ti-1)*numberOfPointsCurve1 + tj]);
          frechetDistancesArray[ti*numberOfPointsCurve2 + tj] = fmaxf(current, frechetDistancesArray[(ti-1)*numberOfPointsCurve2 + tj]);

          if ( debug ) std::cout <<  frechetDistancesArray[ti*numberOfPointsCurve2 + tj]<< " , ";
    }
    // we need to do the same with the original, only need to project the first point again
    ti = 0;


    if ( debug ) std::cout << "Filling column  " << std::endl;

    for(tj = 1; tj < numberOfPointsCurve2; tj++){
        current = objectiveFrechet(method, curve1, curve2, numberOfPointsCurve1,numberOfPointsCurve2, ti, tj );
        //frechetDistancesArray[ti*numberOfPointsCurve1 + tj] = fmaxf(current, frechetDistancesArray[ti*numberOfPointsCurve1 + tj-1]);
        frechetDistancesArray[ti*numberOfPointsCurve2 + tj] = fmaxf(current, frechetDistancesArray[ti*numberOfPointsCurve2 + tj-1]);
    }


    // now that we dealt with the first three cases, now we deal with the last one.
    for(ti = 1; ti < numberOfPointsCurve1; ti++){
        // ti is the projected
        for(tj = 1; tj < numberOfPointsCurve2; tj++){
               current = objectiveFrechet(method, curve1, curve2, numberOfPointsCurve1, numberOfPointsCurve2, ti, tj );
               // we need to look at the 3 previous points.
               //float minv  = frechetDistancesArray[ti*numberOfPointsCurve1 + tj -1];
               //minv = fminf(minv, frechetDistancesArray[(ti-1)*numberOfPointsCurve1 + tj -1]);
               //minv = fminf(minv, frechetDistancesArray[(ti-1)*numberOfPointsCurve1 + tj]);
               //frechetDistancesArray[ti*numberOfPointsCurve1 + tj] = fmaxf( minv, current);
               float minv  = frechetDistancesArray[ti*numberOfPointsCurve2 + tj -1];
               minv = fminf(minv, frechetDistancesArray[(ti-1)*numberOfPointsCurve2 + tj -1]);
               minv = fminf(minv, frechetDistancesArray[(ti-1)*numberOfPointsCurve2 + tj]);
               frechetDistancesArray[ti*numberOfPointsCurve2 + tj] = fmaxf( minv, current);
        }
    }

    return   frechetDistancesArray[dx*numberOfPointsCurve2 + dy];
}


float MathHelper::objectiveFrechet(int method, float* curve1, float* curve2, int numberOfPointsCurve1, int numberOfPointsCurve2, int indexCurve1, int indexCurve2, bool debug){


    // methods to change the objective value to be used in the frechet distance
    // in the original frechet distance calculation this is done using the euclidean distance
    // the second flow based method is used the angle between two vectors
    // the thirds one is based on the second where the integral is used (?) ... need to be certain how to implement it
    // the fourth is one proposed, look at the difference in continuous orientation between points...

    if ( method == 0){
         float p[3] = {curve1[indexCurve1*3 + 0],curve1[indexCurve1*3 + 1],curve1[indexCurve1*3 + 2] };
         float q[3] = {curve2[indexCurve2*3 + 0],curve2[indexCurve2*3 + 1],curve2[indexCurve2*3 + 2] };

         if (debug){
             std::cout << "p " << p[0] << " , " << p[1] << " , " << p[2] << " :: q " << q[0] << ", " << q[1] << ", " << q[2] << std::endl;
         }
         float distance =  sqrt( (p[0] -q[0])*(p[0]-q[0]) +  (p[1] -q[1])*(p[1]-q[1]) + (p[2] -q[2])*(p[2]-q[2]) );

         return distance;
    }
    if ( method == 1){
        float tangentC1[3];
        float tangentC2[3];

        GetTangentFromCurve(curve1, tangentC1, indexCurve1, numberOfPointsCurve1);
        GetTangentFromCurve(curve2, tangentC2, indexCurve2, numberOfPointsCurve2);

        float angle = AngleBwVectors(tangentC1, tangentC2);
        return angle;

    }
    if ( method == 2){
        // Curve integral.... need to figure this one out...

    }

    if ( method == 3){
        // Change of tangent
        float angle1 = RateOfChangeAngle(curve1, indexCurve1, numberOfPointsCurve1, true);
        float angle2 = RateOfChangeAngle(curve2, indexCurve2, numberOfPointsCurve2, true);
        return fabs(angle1- angle2);
    }
    if ( method == 4){
        float angle1 = RateOfChangeAngle(curve1, indexCurve1, numberOfPointsCurve1, false);
        float angle2 = RateOfChangeAngle(curve2, indexCurve2, numberOfPointsCurve2, false);
        return fabs(angle1- angle2);

    }
    return 99999999;
}

void MathHelper::GetCurrentPosition(float* point, int position, float* curve){

    point[0] = curve[position*3 + 0];
    point[1] = curve[position*3 + 1];
    point[2] = curve[position*3 + 2];
}

void MathHelper::GetTangentFromCurve(float* curve, float* tangent, int position, int numberSamples){
    //Lets define the tangent depending on the amount of points available.
    // In the case of being 2 <= i <= n -2 we use a method to get an error of O(h^5)
    // Otherwise we use only the difference.

    // TODO - check this function.
    if ( position >= 2 && position < numberSamples - 2){
        float x_minus2[3], x_minus1[3], x_plus1[3], x_plus2[3];

        GetCurrentPosition(x_minus2, position -2, curve);
        GetCurrentPosition(x_minus1, position -1, curve);
        GetCurrentPosition(x_plus1, position + 1, curve);
        GetCurrentPosition(x_plus2, position + 2, curve);
        for(int i  =0 ; i < 3; i++){ tangent[i]= x_minus2[i] - 8.0f*x_minus1[i] + 8.0f*x_plus1[i] - x_plus2[i]; }
        Normalize( tangent);

    }
    else {
      // boundary conditions.
      if ( position < 0 ) position = 0;

      if ( position == 0){
          float x0[3], x1[3], x2[3], x3[3], x4[3];
          GetCurrentPosition(x0, 0, curve);
          GetCurrentPosition(x1, 1, curve);
          GetCurrentPosition(x2, 2, curve);
          GetCurrentPosition(x3, 3, curve);
          GetCurrentPosition(x4, 4, curve);
          for(int i  =0 ; i < 3; i++){ tangent[i]=   -25.0*x0[i]   + 48.0*x1[i] - 36.0*x2[i]  + 16*x3[i] - 3.0*x4[i]; }
          Normalize( tangent);


      }
      else if ( position == 1){
          float x0[3], x1[3], x2[3], x3[3], x4[3];
          GetCurrentPosition(x0, 0, curve);
          GetCurrentPosition(x1, 1, curve);
          GetCurrentPosition(x2, 2, curve);
          GetCurrentPosition(x3, 3, curve);
          GetCurrentPosition(x4, 4, curve);
          for(int i  =0 ; i < 3; i++){ tangent[i]=   -3.0*x0[i]   -10.0*x1[i] +18.0*x2[i]  - 6*x3[i] + x4[i]; }
          Normalize( tangent);

      }
      else if ( position == numberSamples - 2){ // tn-1
          float x0[3], x1[3], x2[3], x3[3], x4[3];
          GetCurrentPosition(x0, numberSamples - 1, curve);
          GetCurrentPosition(x1, numberSamples - 2, curve);
          GetCurrentPosition(x2, numberSamples - 3, curve);
          GetCurrentPosition(x3, numberSamples - 4, curve);
          GetCurrentPosition(x4, numberSamples - 5, curve);
          for(int i  =0 ; i < 3; i++){ tangent[i]=   3.0*x0[i]  + 10.0*x1[i] - 18.0*x2[i]  + 6*x3[i] - x4[i]; }
          Normalize( tangent);
      }
      else {
          // position == numberSamples -1 or above.

          float x0[3], x1[3], x2[3], x3[3], x4[3];
          GetCurrentPosition(x0, numberSamples - 1, curve);
          GetCurrentPosition(x1, numberSamples - 2, curve);
          GetCurrentPosition(x2, numberSamples - 3, curve);
          GetCurrentPosition(x3, numberSamples - 4, curve);
          GetCurrentPosition(x4, numberSamples - 5, curve);
          for(int i  =0 ; i < 3; i++){ tangent[i]=   25.0*x0[i]  - 48.0*x1[i] + 36.0*x2[i]  - 16.0*x3[i] + 3.0f*x4[i]; }
          Normalize( tangent);
      }
    }
}

float MathHelper::RateOfChangeAngle(float* curve, int index, int numberOfPoints, bool backwards){

    float tangent[3], tangentPrev[3];
    int prev = index -1;
    if ( prev < 0) prev = 0;

    if (!backwards){
        prev = index +1;
        if (prev > numberOfPoints -1) prev = numberOfPoints;
    }

    GetTangentFromCurve(curve, tangent, index, numberOfPoints);
    GetTangentFromCurve(curve, tangentPrev, prev, numberOfPoints);

    float angle = AngleBwVectors(tangent, tangentPrev);
    return angle;
}

float MathHelper::HueToRgb(float p, float q, float t){
        if(t < 0.0f) t += 1.0;
        if(t > 1.0f) t -= 1.0;
        if(t < 1.0/6.0f) return p + (q - p) * 6.0 * t;
        if(t < 1.0f/2.0f) return q;
        if(t < 2.0f/3.0f) return p + (q - p) * (2.0/3.0 - t) * 6.0;
        return p;
}

void   MathHelper::HSLToRGB(float h, float s, float l, float color[]){
   float r, g, b;

    if(s == 0){
        r = g = b = l; // achromatic
    }else{


        float q = (l < 0.5 ) ? l * (1 + s) : l + s - l * s;
        float p = 2 * l - q;
        r = HueToRgb(p, q, h + 1.0f/3.0f);
        g = HueToRgb(p, q, h);
        b = HueToRgb(p, q, h - 1.0f/3.0f);
    }

    color[0] = r;
    color[1] = g;
    color[2] = b;
}




double  MathHelper::CalculateDistancePointSegment2(double *pos1, double *pos2, double *x0){

    double t = MathHelper::CalculateTprojectionPointSegment(pos1, pos2, x0);


    if ( t < 0 || t > 1.0 ) return 99999;

    double dv[3] = { pos2[0] - pos1[0], pos2[1] -pos1[1], pos2[2] -pos1[2]};

    double dvSize = 0;
    for(int i = 0; i < 3; i++) dvSize +=  (dv[i])*(dv[i]);
    dvSize = sqrt(dvSize);
    for(int i = 0; i < 3; i++) dv[i] /= dvSize;

    double projectedPoint[3]; // dv is normalized
    for(int i = 0; i < 3; i++) projectedPoint[i] = pos1[i] + t*dvSize*dv[i]; // t is only from 0.. to 1,

    float distance = 0;
    for(int i = 0; i < 3; i++) distance +=  (projectedPoint[i] - x0[i])*(projectedPoint[i] - x0[i]);

    distance = sqrt(distance);
    return distance;
}
