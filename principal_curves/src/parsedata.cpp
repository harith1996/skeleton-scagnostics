
#include "parsedata.hpp"

#include <locale>
#include <string>
#include <sstream>
#include <iostream>
#include <random>
#include <limits>

float parseFloat(const std::string& str, size_t start, size_t end) {
    std::stringstream stream { str.substr(start) };
    float f;
    stream.imbue(std::locale::classic());
    stream >> f;
    return f;
}



// The generating coefficients create a generating curve... this curve is used to create the points
// and to test the robustness of the data ...
float runPolynomial(const float x,const std::vector<float>& generatingCoefficients){
    float r = 0;
    for (auto c:generatingCoefficients) {
        r = r*x +c;
    }
    return r;
}





std::vector<AugmentedPoint> generatingCircle(const int numPoints, int sz, DataSet originalDataset, int axis1, int axis2){
    std::vector<AugmentedPoint> curve;

    Stat<float> xStat = originalDataset.getStat(axis1);
    Stat<float> yStat = originalDataset.getStat(axis2);

    for(int i = 0; i < numPoints; i++)
    {
        AugmentedPoint pt;
        float percent = (static_cast<float>(i)/static_cast<float>(numPoints));
        float rad = (1.5*3.14159)*percent;
        float x = cos(rad);
        float y = sin(rad);
        pt.point[0] = ((x - xStat.minValue)/( xStat.maxValue-xStat.minValue))*sz;
        pt.point[1] = ((y - yStat.minValue)/( yStat.maxValue-yStat.minValue))*sz;
        curve.push_back(pt);
    }
    return curve;
}


std::vector<AugmentedPoint> generatingCurve(const int numPoints, const std::vector<float>& generatingCoefficients, int sz, DataSet originalDataset, int axis1, int axis2){

    std::vector<AugmentedPoint> curve;

    Stat<float> xStat = originalDataset.getStat(axis1);
    Stat<float> yStat = originalDataset.getStat(axis2);

    for(int i = 0; i < numPoints; i++)
    {
        AugmentedPoint pt;

        float x = (static_cast<float>(i)/static_cast<float>(numPoints))*4 -2.0; // 2-2 range
        float y = runPolynomial(x, generatingCoefficients);

        pt.point[0] = ((x - xStat.minValue)/( xStat.maxValue-xStat.minValue))*sz;
        pt.point[1] = ((y - yStat.minValue)/( yStat.maxValue-yStat.minValue))*sz;
        curve.push_back(pt);

    }
    return curve;
}

DataSet fakeTripleCircle(const int numPoints, float noise){
    DataSet dataSet { 0 };
    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)


   float ranges[3][2] = { {0.0, 2.0}, {0.4, 1.66}, {0.66, 1.5}};

    for(int k = 1; k <= 3; k++ ){

            float r = k*1.5;

            float st = ranges[k-1][0]*3.14159;
            float end = ranges[k-1][1]*3.14159;

            std::uniform_real_distribution<float> uni(st, end);// guaranteed unbiased

            std::normal_distribution<float> distribution(0, 0.01*noise*k);

            for(int i = 0; i < numPoints; i++)
            {

               std::vector<float> row;

               float t = uni(rng);
               float x = 0;
               float y = 0;

               x = r*cos(t);      y = r*sin(t);

               x += distribution(rng);
               y += distribution(rng);

               row.push_back(x);
               row.push_back(y);
               if(dataSet.dimension() == 0) {
                    dataSet = DataSet(row.size());
               }
               //TODO check row size
               dataSet.append(row);
            }
    }

    std::uniform_real_distribution<float> uni2(0, 3.0);// guaranteed unbiased
    std::uniform_real_distribution<float> uni3(0, 3.5);// guaranteed unbiased

    std::normal_distribution<float> distribution(0, 0.015*noise);

    for(int i = 0; i < numPoints; i++)
    {

       std::vector<float> row, row2;

       float r = 1.0 + uni2(rng);
       float r2 = 1.0 + uni3(rng);

       float t = (0.25)*3.14159;
       float t2 = (1.25)*3.14159;

       float x = 0;
       float y = 0;

       x = r*cos(t);      y = r*sin(t);

       x += distribution(rng);
       y += distribution(rng);

       row.push_back(x);
       row.push_back(y);


       x = r2*cos(t2);      y = r2*sin(t2);

       x += distribution(rng);
       y += distribution(rng);

       row2.push_back(x);
       row2.push_back(y);

       if(dataSet.dimension() == 0) {
            dataSet = DataSet(row.size());
       }
       //TODO check row size
       dataSet.append(row);
       dataSet.append(row2);

    }
    return dataSet;


}


std::vector<AugmentedPoint> generatingTripleCircle(const int numPoints, int sz, DataSet originalDataset, int axis1, int axis2){
     std::vector<AugmentedPoint> curve;

    Stat<float> xStat = originalDataset.getStat(axis1);
    Stat<float> yStat = originalDataset.getStat(axis2);


   float ranges[3][2] = { {0.0, 2.0}, {0.4, 1.66}, {0.66, 1.5}};

    for(int k = 1; k <= 3; k++ ){

            float r = k*1.5;
            float st = ranges[k-1][0]*3.14159;
            float end = ranges[k-1][1]*3.14159;
            for(int i = 0; i < numPoints/2; i++)
            {


               float t = (end - st)*(static_cast<float>(i)/(static_cast<float>(numPoints/2))) + st;
               float x = 0;
               float y = 0;

               AugmentedPoint pt;

               x = r*cos(t);      y = r*sin(t);
               pt.point[0] = ((x - xStat.minValue)/( xStat.maxValue-xStat.minValue))*sz;
               pt.point[1] = ((y - yStat.minValue)/( yStat.maxValue-yStat.minValue))*sz;
               curve.push_back(pt);
               //x += distribution(rng);
               //y += distribution(rng);
            }
    }

    for(int i = 0; i < numPoints/10; i++)
    {


       float percent = (static_cast<float>(i)/(static_cast<float>(numPoints/10.0)));
       float r = 1.0 + 3.0*percent;
       float r2 = 1.0 + 3.5*percent;

       float t = (0.25)*3.14159;
       float t2 = (1.25)*3.14159;

       float x = 0;
       float y = 0;

       AugmentedPoint pt;

       x = r*cos(t);      y = r*sin(t);

       //x += distribution(rng);
       //y += distribution(rng);
       pt.point[0] = ((x - xStat.minValue)/( xStat.maxValue-xStat.minValue))*sz;
       pt.point[1] = ((y - yStat.minValue)/( yStat.maxValue-yStat.minValue))*sz;
       curve.push_back(pt);



       x = r2*cos(t2);      y = r2*sin(t2);

       //x += distribution(rng);
       //y += distribution(rng);

       //row2.push_back(x);
       //row2.push_back(y);
       AugmentedPoint pt2;

       pt2.point[0] = ((x - xStat.minValue)/( xStat.maxValue-xStat.minValue))*sz;
       pt2.point[1] = ((y - yStat.minValue)/( yStat.maxValue-yStat.minValue))*sz;
       curve.push_back(pt2);


    }
    return curve;
}

std::vector<AugmentedPoint> generatingSquareRoot(const int numPoints, int sz, DataSet originalDataset, int axis1, int axis2){

    std::vector<AugmentedPoint> curve;

    Stat<float> xStat = originalDataset.getStat(axis1);
    Stat<float> yStat = originalDataset.getStat(axis2);

    for(int i = 0; i < numPoints; i++)
    {
        AugmentedPoint pt;

        float x = (static_cast<float>(i)/static_cast<float>(numPoints))*5.0;
        float y = sqrt(x);

        pt.point[0] = ((x - xStat.minValue)/( xStat.maxValue-xStat.minValue))*sz;
        pt.point[1] = ((y - yStat.minValue)/( yStat.maxValue-yStat.minValue))*sz;
        curve.push_back(pt);

    }
    return curve;



}

std::vector<AugmentedPoint> generatingLine(const int numPoints, int sz, DataSet originalDataset, int axis1, int axis2){
    std::vector<AugmentedPoint> curve;

    Stat<float> xStat = originalDataset.getStat(axis1);
    Stat<float> yStat = originalDataset.getStat(axis2);

    for(int i = 0; i < numPoints; i++)
    {
        AugmentedPoint pt;

        float x = (static_cast<float>(i)/static_cast<float>(numPoints))*5.0;
        float y = x;

        pt.point[0] = ((x - xStat.minValue)/( xStat.maxValue-xStat.minValue))*sz;
        pt.point[1] = ((y - yStat.minValue)/( yStat.maxValue-yStat.minValue))*sz;
        curve.push_back(pt);

    }
    return curve;
}

void fakeLinear(const int numPoints, float noise,  DataSet& normal, DataSet& up, DataSet& down){


    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_real_distribution<float> uni(0, 5.0);// guaranteed unbiased
    std::normal_distribution<float> distribution(0, 0.015*(noise));

    int offsided = 3;

    for(int i = 0; i < numPoints; i++)
    {

       std::vector<float> row, row2;
       float x =  uni(rng);
       float y =  x;

       float noise_y =  fabs(distribution(rng));
       //x += distribution(rng);
       float y1 = y + noise_y;
       float y2 = y - noise_y;

       row.push_back(x);
       row.push_back(y1);

       row2.push_back(x);
       row2.push_back(y2);


       if(normal.dimension() == 0) {
            normal = DataSet(row.size());
            up = DataSet(row.size());
            down = DataSet(row.size());
       }

           for(int k =0; k < offsided;k++ )
              normal.append(row);
           for(int k =0; k < offsided
               ;k++ )
              normal.append(row2);


           if ( noise_y > 0.015*(noise*0.75) && noise_y < 0.015*(noise*0.25) ){ //if the noise is too big, lets not put it, as not to change the shape...

               while(normal.size() != up.size()){
                  up.append(row);
                  up.append(row2);
                  down.append(row);
                  down.append(row2);
               }
           }
           else{
               down.append(row); // just skewed ...
               up.append(row2);
               while(normal.size() != down.size()){
                   down.append(row2);
                   up.append(row);
               }
           }


    }

    // now add skewness
}


DataSet fakeSpiral(const int numPoints, float noise){

    DataSet dataSet { 0 };
    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_real_distribution<float> uni(0, 2*3.14159);// guaranteed unbiased
    std::normal_distribution<float> distribution(0, 0.015*noise);

    std::cout << "Noise " << 0.015*noise << std::endl;
    for(int i = 0; i < numPoints; i++)
    {

       std::vector<float> row, row2;
       float t =  uni(rng);

       float x = -cos(t)*t + distribution(rng);
       float y = sin(t)*t + distribution(rng);

       row.push_back(x);
       row.push_back(y);
       row2.push_back(-x);
       row2.push_back(-y);

       if(dataSet.dimension() == 0) {
            dataSet = DataSet(row.size());
       }

       //TODO check row size
       dataSet.append(row);
       dataSet.append(row2);
    }


    return dataSet;
}

std::vector<AugmentedPoint> generatingSpiral(const int numPoints, int sz, DataSet originalDataset, int axis1, int axis2){
    std::vector<AugmentedPoint> curve;

    Stat<float> xStat = originalDataset.getStat(axis1);
    Stat<float> yStat = originalDataset.getStat(axis2);

    for(int i = 0; i < numPoints; i++)
    {
        AugmentedPoint pt;
        float t = ((static_cast<float>(numPoints-i)/static_cast<float>(numPoints)))*2*3.14159;
        float x = -cos(t)*t;
        float y = sin(t)*t;
        pt.point[0] = ((x - xStat.minValue)/( xStat.maxValue-xStat.minValue))*sz;
        pt.point[1] = ((y - yStat.minValue)/( yStat.maxValue-yStat.minValue))*sz;
        curve.push_back(pt);
    }

    for(int i = 0; i < numPoints; i++)
    {
        AugmentedPoint pt;
        float t = ((static_cast<float>(i)/static_cast<float>(numPoints)))*2*3.14159;
        float x = cos(t)*t;
        float y = -sin(t)*t;
        pt.point[0] = ((x - xStat.minValue)/( xStat.maxValue-xStat.minValue))*sz;
        pt.point[1] = ((y - yStat.minValue)/( yStat.maxValue-yStat.minValue))*sz;
        curve.push_back(pt);
    }

    return curve;
}

DataSet fakeX(const int numPoints, float noise){
    DataSet dataSet { 0 };
    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_real_distribution<float> uni(0.25*3.14159, 1.75*3.14159);// guaranteed unbiased
    std::normal_distribution<float> distribution(0, 0.015*noise);

    for(int i = 0; i < numPoints; i++)
    {

       std::vector<float> row, row2;

       float t=  uni(rng);
       float x = 0;
       float y = 0;

       x = cos(t);
       y = sin(t);

       x += distribution(rng);
       y += distribution(rng);

       row.push_back(x+1);
       row.push_back(y);

       row2.push_back(-x-1);
       row2.push_back(y);


       if(dataSet.dimension() == 0) {
            dataSet = DataSet(row.size());
       }


       //TODO check row size
       dataSet.append(row);
       dataSet.append(row2);

    }


    return dataSet;

}

std::vector<AugmentedPoint> generatingX(const int numPoints, int sz, DataSet originalDataset, int axis1, int axis2){
    std::vector<AugmentedPoint> curve;
    Stat<float> xStat = originalDataset.getStat(axis1);
    Stat<float> yStat = originalDataset.getStat(axis2);

    for(int i = 0; i < numPoints; i++)
    {
       AugmentedPoint pt, pt2;

       float t = ((static_cast<float>(numPoints-i)/static_cast<float>(numPoints)))*1.75*3.14159  +0.25*3.14159;
       float x = 0;
       float y = 0;

       x = cos(t);
       y = sin(t);


       pt.point[0] = ((x+1 - xStat.minValue)/( xStat.maxValue-xStat.minValue))*sz;
       pt.point[1] = ((y - yStat.minValue)/( yStat.maxValue-yStat.minValue))*sz;
       curve.push_back(pt);


       pt2.point[0] = ((-x-1 - xStat.minValue)/( xStat.maxValue-xStat.minValue))*sz;
       pt2.point[1] = ((y - yStat.minValue)/( yStat.maxValue-yStat.minValue))*sz;
       curve.push_back(pt2);
    }
    return curve;
}

std::vector<AugmentedPoint> generatingDNA(const int numPoints, int sz, DataSet originalDataset, int axis1, int axis2){
    std::vector<AugmentedPoint> curve;
    Stat<float> xStat = originalDataset.getStat(axis1);
    Stat<float> yStat = originalDataset.getStat(axis2);

    for(int i = 0; i < numPoints; i++)
    {
       AugmentedPoint pt, pt2;

       float x = ((static_cast<float>(numPoints-i)/static_cast<float>(numPoints)))*3*3.14159;
       float y1 = (1.0 +x)*cos(x); // increase the amplitude along the distance
       float y2 = -(1.0 + x)*cos(x);

       pt.point[0] = ((x - xStat.minValue)/( xStat.maxValue-xStat.minValue))*sz;
       pt.point[1] = ((y1 - yStat.minValue)/( yStat.maxValue-yStat.minValue))*sz;
       curve.push_back(pt);

       pt2.point[0] = ((x - xStat.minValue)/( xStat.maxValue-xStat.minValue))*sz;
       pt2.point[1] = ((y2 - yStat.minValue)/( yStat.maxValue-yStat.minValue))*sz;
       curve.push_back(pt2);
    }
    return curve;
}

DataSet growingDNA(const int numPoints, float noise){
    DataSet dataSet { 0 };
    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)

    float st = 0*3.14159;
    float end = 3*3.14159;

    std::uniform_real_distribution<float> uni(st, end);// guaranteed unbiased
    std::normal_distribution<float> distribution(0, 0.015*noise);

    for(int i = 0; i < numPoints; i++)
    {


        std::vector<float> row, row2;

        float x = uni(rng);
        float y1 = (1.0 +x)*cos(x); // increase the amplitude along the distance
        float y2 = -(1.0 + x)*cos(x);

        x += distribution(rng);

        y1 += distribution(rng);
        y2 += distribution(rng);


        row.push_back(x);
        row.push_back(y1);
        if(dataSet.dimension() == 0) {
             dataSet = DataSet(row.size());
        }


        row2.push_back(x);
        row2.push_back(y2);

        //TODO check row size
        dataSet.append(row);
        dataSet.append(row2);

    }
    return dataSet;
}

DataSet growingDNA2(const int numPoints, float noise){
    DataSet dataSet { 0 };
    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)



    std::uniform_real_distribution<float> uni1(0, 1*3.14159);// guaranteed unbiased
    std::normal_distribution<float> distribution1(0, 0.015*noise*1);


    std::uniform_real_distribution<float> uni2(1*3.14159, 2*3.14159);// guaranteed unbiased
    std::normal_distribution<float> distribution2(0, 0.015*noise*2);


    std::uniform_real_distribution<float> uni3(2*3.14159, 3*3.14159);// guaranteed unbiased
    std::normal_distribution<float> distribution3(0, 0.015*noise*3);


    for(int i = 0; i < numPoints/3; i++)
    {
        std::vector<float> row, row2;
        float x = uni1(rng);
        float y1 = (1.0 +x)*cos(x); // increase the amplitude along the distance
        float y2 = -(1.0 + x)*cos(x);
        x += distribution1(rng);        y1 += distribution1(rng);        y2 += distribution1(rng);

        row.push_back(x);        row.push_back(y1);
        if(dataSet.dimension() == 0) {
             dataSet = DataSet(row.size());
        }
        row2.push_back(x);        row2.push_back(y2);
        dataSet.append(row);        dataSet.append(row2);
    }

    for(int i = 0; i < numPoints/3; i++)
    {
        std::vector<float> row, row2;
        float x = uni2(rng);
        float y1 = (1.0 +x)*cos(x); // increase the amplitude along the distance
        float y2 = -(1.0 + x)*cos(x);
        x += distribution2(rng);        y1 += distribution2(rng);        y2 += distribution2(rng);

        row.push_back(x);        row.push_back(y1);
        if(dataSet.dimension() == 0) {
             dataSet = DataSet(row.size());
        }
        row2.push_back(x);        row2.push_back(y2);
        dataSet.append(row);        dataSet.append(row2);
    }

    for(int i = 0; i < numPoints/3; i++)
    {
        std::vector<float> row, row2;
        float x = uni3(rng);
        float y1 = (1.0 +x )*cos(x); // increase the amplitude along the distance
        float y2 = -(1.0 + x)*cos(x);
        x += distribution3(rng);        y1 += distribution3(rng);        y2 += distribution3(rng);

        row.push_back(x);        row.push_back(y1);
        if(dataSet.dimension() == 0) {
             dataSet = DataSet(row.size());
        }
        row2.push_back(x);        row2.push_back(y2);
        dataSet.append(row);        dataSet.append(row2);
    }


    return dataSet;
}

DataSet fakeCircle(const int numPoints, float noise){
    DataSet dataSet { 0 };
    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)

    float st = 0;
    float end = 1.5*3.14159;

    std::uniform_real_distribution<float> uni(st, end);// guaranteed unbiased
    std::normal_distribution<float> distribution(0, 0.015*noise);

    for(int i = 0; i < numPoints; i++)
    {

       std::vector<float> row;

       float t=  uni(rng);
       float x = 0;
       float y = 0;

       x = cos(t);
       y = sin(t);

       x += distribution(rng);
       y += distribution(rng);

       row.push_back(x);
       row.push_back(y);
       if(dataSet.dimension() == 0) {
            dataSet = DataSet(row.size());
       }


       //TODO check row size
       dataSet.append(row);
    }

    return dataSet;
}

DataSet fakeData(const int numPoints,const std::vector<float>& generatingCoefficients, float noise){

   DataSet dataSet { 0 };


   std::random_device rd;     // only used once to initialise (seed) engine
   std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
   std::uniform_real_distribution<float> uni(-2,2);// guaranteed unbiased
   std::normal_distribution<float> distribution(0, 0.015*noise);

   for(int i = 0; i < numPoints; i++)
   {

      std::vector<float> row;

      float x = 0;
      float y = 0;

      x = uni(rng);
      y = runPolynomial(x, generatingCoefficients);

      x += distribution(rng);
      y += distribution(rng);

      row.push_back(x);
      row.push_back(y);
      if(dataSet.dimension() == 0) {
           dataSet = DataSet(row.size());
      }

      //TODO check row size
      dataSet.append(row);
   }




   return dataSet;
}

DataSet parseData(std::istream& infile, char delimiter, int multiplier, bool debug) {
    infile.imbue(std::locale::classic());

    DataSet     dataSet { 0 };

    std::vector<float> row;

    for(std::string lineString; std::getline(infile, lineString); ) {
        row.resize(0);

        std::istringstream line { lineString };

        //munch whitespace
        line >> std::ws;

        while(!line.eof()) {
            float value;
            line >> value;
            if(!line.fail()) {
                row.push_back(value);
            } else {
                row.push_back(0);
                //unset erro flag
                line.clear();
            }

            line >> std::ws;

            char c;
            do {
                if(line.eof()) break;
                c = line.get();
            } while(c != delimiter);
            line >> std::ws;
        }

        //skip empty rows
        if(row.empty()) continue;

        if (debug){

            for(int k = 0; k < row.size(); k++){
                std::cout << row.at(k) << "_ ";
            }
            std::cout << std::endl;
        }
        //init dataSet when with first non-empty row
        if(dataSet.dimension() == 0) {
            dataSet = DataSet(row.size());
        }

        //Append the same row several times, just to get a higher sized
        //dataset if multiplier > 1
        for(int j = 0; j < multiplier; j++)
            dataSet.append(row);
    }

    return dataSet;
}

