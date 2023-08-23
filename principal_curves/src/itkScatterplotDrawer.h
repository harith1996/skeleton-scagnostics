#ifndef __itkGaussianScatterplotDrawer_h
#define __itkGaussianScatterplotDrawer_h

#include "itkImageToImageFilter.h"
#include <vector>
#include "dataset.hpp"


namespace itk
{
    template< class TImage>
    class GaussianScatterplotDrawerFilter : public ImageToImageFilter< TImage, TImage>
    {
        public:
            typedef  GaussianScatterplotDrawerFilter Self;
            typedef ImageToImageFilter<TImage, TImage> Superclass;
            typedef SmartPointer< Self> Pointer;
            
            typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
            typedef typename TImage::IndexType IndexType;
            typedef typename TImage::SpacingType SpacingType;

            itkNewMacro(Self);
            
            itkTypeMacro( GaussianScatterplotDrawerFilter, ImageToImageFilter);
        
            void SetHelperInfo(DataSet* dataSet, std::pair<int,int> indices, float sigma, std::pair<float,float> xRange,
                               std::pair<float,float> yRange, int width, int height)
            {
                    _dataset = dataSet;
                    _indices = indices;
                    _sigma = sigma;
                    _xRange = xRange;
                    _yRange = yRange;
                    _width = width;
                    _height = height;
            }


            double Gaussian2DDistribution(double x, double y, double stddev){
                  // mean 0 and std dev as set
                  double inverse = 1.0 /( 2.0 * 3.14159 * stddev* stddev);
                  double diff =  (pow(x,2.0) + pow(y,2.0)) / (2.0* stddev * stddev);
                  return inverse * exp( -diff );
            }

    protected:
            GaussianScatterplotDrawerFilter(){}
           ~ GaussianScatterplotDrawerFilter(){}
           virtual void ThreadedGenerateData(const OutputImageRegionType &, ThreadIdType);
        private:
             GaussianScatterplotDrawerFilter( const Self &);
            void operator=(const Self&); 
            DataSet* _dataset;
            std::pair<int,int> _indices;
            float _sigma;
            std::pair<float,float> _xRange;
            std::pair<float,float> _yRange;
            int _height;
            int _width;
    };
}

#ifndef ITK_MANUAL_INSTALLATION
#include "itkScatterplotDrawer.hxx"
#endif

#endif
