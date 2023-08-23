
#include "projection2d.hpp"
#include <iostream>

#include "itkScatterplotDrawer.h"

SimpleAxisSelection2D::SimpleAxisSelection2D(size_t x, size_t y) : cx { x }, cy { y } {
}

Float2 SimpleAxisSelection2D::operator()(const DataSet& dataSet, size_t index) {
    return { dataSet(cx, index), dataSet(cy, index) };
}

AxisSelection2D::AxisSelection2D(size_t x, size_t y) : cx { x }, cy { y } {
}

Float2 AxisSelection2D::operator()(const DataSet& dataSet, size_t index) {
    Float2 point = { dataSet(cx, index), dataSet(cy, index) };
    Float2 min = { dataSet.getStat(cx).minValue, dataSet.getStat(cy).minValue };
    Float2 max = { dataSet.getStat(cx).maxValue, dataSet.getStat(cy).maxValue };

    return (point - min) / (max - min);
}

ProjectionMatrix2D::ProjectionMatrix2D(size_t n) : entries {} {
    entries.resize(n);
    for(int i = 0; i < n; i++) {
        entries[0] = 0;
    }
}

size_t ProjectionMatrix2D::dimension() const {
    return entries.size();
}

Float2& ProjectionMatrix2D::entry(size_t index) {
    return entries.at(index);
}

Float2 ProjectionMatrix2D::entry(size_t index) const {
    return entries.at(index);
}

Float2 CalibratedMatrix2D::operator()(const DataSet& dataSet, size_t index) {
    Float2 result = 0;

    for(int i = 0; i < dataSet.dimension() && i < matrix->dimension(); i++) {
        auto value = dataSet(i, index);
        auto minValue = dataSet.getStat(i).minValue;
        auto maxValue = dataSet.getStat(i).maxValue;
        if(maxValue != minValue) {
            result += matrix->entry(i) * dataSet(i, index);
        }
    }

    result = (result - min) / (max - min);

    return result;
}

CalibratedMatrix2D calibrateMatrix(ProjectionMatrix2D *matrix, const DataSet& dataSet) {
    CalibratedMatrix2D result;

    result.matrix = matrix;
    result.min = Float2({0, 0});
    result.max = Float2({1, 1});

    //compute normalization
    Stat<Float2> stat;
    for(int i = 0; i < dataSet.size(); i++) {
        //use 0-1-calibrated matrix to calculate min and max
        //so it can it be sensibly normalized
        stat.put(result(dataSet, i));
    }
    result.min = stat.minValue;
    result.max = stat.maxValue;

    return result;
}

double NormalDistribution(double x, double stddev){
      // mean 0 and std dev as set
      double inverse = 1.0 / sqrt(2.0 * 3.14159 * stddev );
      double diff =  pow(x,2.0) / (2.0* stddev * stddev);
      return inverse * exp( -diff );
}

double Gaussian2DDistribution(double x, double y, double stddev){
      // mean 0 and std dev as set
      double inverse = 1.0 /( 2.0 * 3.14159 * stddev* stddev);
      double diff =  (pow(x,2.0) + pow(y,2.0)) / (2.0* stddev * stddev);
      return inverse * exp( -diff );
}

Image2D<float> gaussedPlotWithItk(DataSet dataSet, int firstAxis, int secondAxis ,
                                  size_t width, size_t height, float sigma, float* sum){

    auto xStat = dataSet.getStat( firstAxis );
    auto yStat = dataSet.getStat( secondAxis );
    Image2D<float> fnl {width, height};


    typedef itk::Image< float, 2 > HImageType;
    HImageType::Pointer img = HImageType::New();
    HImageType::RegionType region;
    HImageType::IndexType index;
    index[0] = 0;      index[1] = 0;
    HImageType::SizeType sizeT;
    sizeT[0] = width;
    sizeT[1] = height;
    region.SetIndex(index);
    region.SetSize(sizeT);
    img->SetRegions(region);
    img->Allocate(true);
    img->FillBuffer(itk::NumericTraits< HImageType::PixelType >::Zero);
    img->Update();

   typedef itk::GaussianScatterplotDrawerFilter<HImageType> drawerFilter;
   drawerFilter::Pointer drawer = drawerFilter::New();



   drawer->SetHelperInfo(&dataSet, std::make_pair(firstAxis,secondAxis), sigma, std::make_pair(xStat.minValue,xStat.maxValue),
                         std::make_pair(yStat.minValue, yStat.maxValue), width, height);
   drawer->SetInput(img);
   drawer->Update();

   HImageType::Pointer res = drawer->GetOutput();
   itk::ImageRegionIteratorWithIndex<HImageType> it(res, region);
   //double minValueNonZero = 999999;
   float avg = 0;
   while(!it.IsAtEnd()){

       auto idx = it.GetIndex();
       int x = idx[0];
       int y = idx[1];
        fnl(x,y) = it.Get();

        avg =  std::max(it.Get(),avg); //it.Get();
       ++it;
   }

   *sum = avg;
   return fnl;
}

Image2D<float> gaussedPlot(const DataSet& dataSet, int firstAxis, int secondAxis , size_t width, size_t height, float sigma){
    Image2D<float> img {width, height};

    //Let's make a gauss version of this ..
    //Let's see how it works ...
    auto xStat = dataSet.getStat( firstAxis );
    auto yStat = dataSet.getStat( secondAxis );

    float xRange = xStat.maxValue - xStat.minValue;
    float yRange = yStat.maxValue - yStat.minValue;

    float xSpacing = 1.0/(width);
    float ySpacing = 1.0/(height);

    for(int x = 0; x < width; x++) {
        float o1 = x*xSpacing + xSpacing*0.5; //Center of pixel


        for(int y = 0; y < height; y++) {
            img(x,y) = 0;
            float v = 0;
            float o2 = y*ySpacing + ySpacing*0.5; //Center of pixel

            for(int k = 0; k < dataSet.size(); k++){

               float v1 = (dataSet(firstAxis,k) - xStat.minValue)/xRange;
               float v2 = (dataSet(secondAxis,k) -yStat.minValue)/yRange;


               // if (x==0&& y ==0)
               //    std::cout << "(" << v1 <<","<< v2 <<")  ";
               v +=  Gaussian2DDistribution(o1-v1,o2-v2,sigma);
            }
            img(x,y) = v/dataSet.size();

            //std::cout << std::endl;
            //std::cout << "(x,y) "<< x <<"," << y <<".."  << img(x,y) <<"::>" << o1 <<"," << o2 <<std::endl;
        }
    }

    return img;
}


Image2D<unsigned> scatterPlotWithAccumulator(const DataSet& dataSet, Projection2D projection, size_t width, size_t height,Image2D< std::vector<std::pair<float,float>>* >* accumulator) {
    Image2D<unsigned> img {width, height};
    for(int x = 0; x < width; x++) {
        for(int y = 0; y < height; y++) {
             (*accumulator)(x,y) = new std::vector<std::pair<float,float>>();
        }
    }


    int maxAccum = 0;

    for(int i = 0; i < dataSet.size(); i++) {
        auto pair = projection(dataSet, i);
        float v1 = pair[0];
        float v2 = pair[1];

        int x = v1 * width;
        x = std::max(0, std::min(x, (int)(width - 1)));
        int y = v2 * height;
        y = std::max(0, std::min(y, (int)(height - 1)));
        (*accumulator)(x,y)->push_back( std::make_pair(v1,v2));

        img(x, y)++;
        if ( img(x,y) > maxAccum) maxAccum = img(x,y);
    }

    //std::cout << "??? " << maxAccum << std::endl;
    return img;
}

Image2D<unsigned> scatterPlot(const DataSet& dataSet, Projection2D projection, size_t width, size_t height) {
    Image2D<unsigned> img {width, height};


    for(int i = 0; i < dataSet.size(); i++) {
        auto pair = projection(dataSet, i);
        float v1 = pair[0];
        float v2 = pair[1];

        int x = v1 * width;
        x = std::max(0, std::min(x, (int)(width - 1)));
        int y = v2 * height;
        y = std::max(0, std::min(y, (int)(height - 1)));
        img(x, y)++;
    }

    return img;
}
