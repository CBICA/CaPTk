#pragma once

#include "itkThresholdImageFilter.h"
#include "itkMaskedImageToHistogramFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkDivideImageFilter.h"

#include "cbicaITKUtilities.h"
#include "cbicaStatistics.h"

template< class TImageType = itk::Image< float, 3 > >
class P1P2Normalizer
{
  using TIteratorType = itk::ImageRegionIteratorWithIndex< TImageType >;
  using TConstIteratorType = itk::ImageRegionConstIterator< TImageType >;

public:
  P1P2Normalizer() {};
  ~P1P2Normalizer() {};

  //! Sets the input image
  void SetInputImage(typename TImageType::Pointer image);

  //! Actual computer
  void Update();

  typename TImageType::Pointer GetOutput();

private:

  //! Get basic statistics from image
  std::map< std::string, double > GetStatisticsForImage(const typename TImageType::Pointer m_inputImage, bool considerMask = true);

  typename TImageType::Pointer m_inputImage, m_mask, m_output;
  float m_quantLower = 0.02, m_quantUpper = 0.95;
  bool m_algorithmDone = false;
};

#include "P1P2Normalizer.hxx"