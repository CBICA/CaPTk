/**
\file TextureFeatureBase.h

This file contains the base for all Feature Extraction classes.

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html
*/
#pragma once

#include <typeinfo>
#include <map>
#include <cmath>

#include "itkImage.h"
#include "itkVectorContainer.h"
#include "itkHistogram.h"

#include "cbicaLogging.h"

//! the histogram binning type
enum HistogramBinningType
{
  FixedBinNumber,
  FixedBinSize,
  Equal
};

/**
\class TextureFeatureBase

\brief This is the base class for all histogram-based feature extraction classes that are not ITK-based

*/
template< class TImageType = itk::Image< float, 3 > >
class TextureFeatureBase
{
public:
  using OffsetVectorPointer = typename itk::VectorContainer< unsigned char, typename TImageType::OffsetType >::Pointer;
  
  //! Default Constructor
  TextureFeatureBase() {};

  //! Default Destructor
  ~TextureFeatureBase() {};

  //! Actual algorithm runner - this is to be changed for each feature family
  virtual void Update() {};

  //! Set the pixel minimum
  void SetMinimum(typename TImageType::PixelType min) { m_minimum = min; };

  //! Set the pixel maximum
  void SetMaximum(typename TImageType::PixelType max) { m_maximum = max; };

  //! Set the offsets
  void SetOffsets(OffsetVectorPointer offsets)
  {
    m_offsets = offsets;
  }

  /**
  \brief Set the number of bins for quantization of the input image; defaults to 10 unless overridden by user

  \param numBinValue integer value for the number of bins you want for image quantization
  **/
  void SetNumBins(unsigned int numBinValue)
  {
    m_Bins = numBinValue;
  }

  //! Set the histogram binning type
  void SetHistogramBinningType(int type)
  {
    m_histogramBinningType = type;
  }

protected:

  unsigned int m_Bins = 10; //! the binning information

  OffsetVectorPointer m_offsets; //! the offsets to consider

  typename TImageType::PixelType 
    m_minimum = 0, //! the minimum to consider during binning
    m_maximum = 0; //! the maximum to consider during binning

  itk::Statistics::Histogram< double >::Pointer m_histogram; //! the actual histogram 

  int m_histogramBinningType; //! the default binning type
};

#include "TextureFeatureBase.hxx"