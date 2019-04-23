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

  //! Set the quantile locations - this should not be changed by the developer
  void SetQuantiles(float lower, float upper);

  //! Get basic statistics from image
  std::map< std::string, double > GetStatisticsForImage(const typename TImageType::Pointer m_inputImage, bool considerMask = true);

  typename TImageType::Pointer m_inputImage, m_mask, m_inputImageMasked, m_output;
  float m_quantLower = 0.05, m_quantUpper = 0.95;
  float m_cutoffLower = 3, m_cutoffUpper = 3;
  bool m_wholeImageMeanThreshold = false;
  bool m_doSanityCheck = true; //! not needed when called from CaPTk since it does those checks already
  bool m_algorithmDone = false;
};

#include "P1P2Normalizer.hxx"