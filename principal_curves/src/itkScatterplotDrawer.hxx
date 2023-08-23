#ifndef __itkQuadTreeBinaryDrawFilter_hxx
#define __itkQuadTreeBinaryDrawFilter_hxx

#include "itkScatterplotDrawer.h"

#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkPixelTraits.h"
#include "itkImageAlgorithm.h"




namespace itk
{

  template<class TImage>
  void GaussianScatterplotDrawerFilter<TImage>
  ::ThreadedGenerateData( const OutputImageRegionType & region, ThreadIdType threadId)
  {
     typename TImage::ConstPointer input = this->GetInput();
     typename TImage::Pointer output = this->GetOutput();
     
     ImageAlgorithm::Copy(input.GetPointer(), output.GetPointer(), region, region);
     itk::ImageRegionIteratorWithIndex<TImage> it(output, region);

     float xSpacing = 1.0/(_width);
     float ySpacing = 1.0/(_height);
     float xRange = _xRange.second - _xRange.first;
     float yRange = _yRange.second - _yRange.first;
     int n = _dataset->size();
     while( ! it.IsAtEnd()){

         auto idx = it.GetIndex();
         int x = idx[0];
         int y = idx[1];
         float o1 = x*xSpacing + xSpacing*0.5;
         float o2 = y*ySpacing + ySpacing*0.5;


         double v=0;
         for(int i =0; i < n; i++){
             float v1 = ((*_dataset)(_indices.first,i) - _xRange.first)/xRange;
             float v2 = ((*_dataset)(_indices.second,i) -_yRange.first)/yRange;
             v += Gaussian2DDistribution(o1-v1,o2-v2,_sigma);
         }

         v /= n;
         it.Set(v);

         ++it;
     }
  }




}

#endif
