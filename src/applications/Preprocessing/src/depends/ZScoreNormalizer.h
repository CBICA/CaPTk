#pragma once

#include "itkConstantPadImageFilter.h"
#include "itkNeighborhoodIterator.h"
#include "itkShrinkImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkMaskedImageToHistogramFilter.h"
#include "itkImageToListSampleAdaptor.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkDivideImageFilter.h"

#include "cbicaITKUtilities.h"
#include "cbicaStatistics.h"

template< class TImageType = itk::Image< float, 3 > >
class ZScoreNormalizer
{
  using TIteratorType = itk::ImageRegionIteratorWithIndex< TImageType >;
  using TConstIteratorType = itk::ImageRegionConstIterator< TImageType >;

public:
  ZScoreNormalizer() {};
  ~ZScoreNormalizer() {};

  //! Sets the input image
  void SetInputImage(typename TImageType::Pointer image);

  /**
  \brief Sets the input mask

  If no mask is provided, it is assumed that all non-zero pixel/voxel values are to be normalized
  */
  void SetInputMask(typename TImageType::Pointer image);

  //! Set the quantile locations
  void SetQuantiles(float lower, float upper);

  //! Set the cutoff locations
  void SetCutoffs(float lower, float upper);

  //! Actual computer
  void Update();

  typename TImageType::Pointer GetOutput();

private:

  //! Get basic statistics from image
  std::map< std::string, double > GetStatisticsForImage(const typename TImageType::Pointer m_inputImage, bool considerMask = true);

  typename TImageType::Pointer m_inputImage, m_mask, m_inputImageMasked, m_output;
  float m_quantLower = 0.05, m_quantUpper = 0.95;
  float m_cutoffLower = 3, m_cutoffUpper = 3;
  bool m_wholeImageMeanThreshold = false;
  bool m_doSanityCheck = true; //! not needed when called from CaPTk since it does those checks already
  bool m_algorithmDone = false;
};

#include "ZScoreNormalizer.hxx"