/**
\file FeatureBase.h

This file contains the base for all Feature Extraction classes.

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html
*/
#pragma once

#include <typeinfo>
#include <map>
#include <cmath>

#include "itkImage.h"
#include "itkNeighborhoodIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkMath.h"

#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"
#include "cbicaLogging.h"

auto const PI = itk::Math::pi; //! initialize PI

/**
\class FeatureBase

\brief This is the base class for all feature extraction classes that are not ITK-based

The idea is that all classes are to inherit from this base class so that all common functions are populated.
*/
template< class TImageType = itk::Image< float, 3 > >
class FeatureBase
{
public:
  
  using TConstIteratorType = itk::ImageRegionConstIterator< TImageType >;
  using TIteratorType = itk::ImageIteratorWithIndex< TImageType >;
  using TRegionIteratorType = itk::ImageRegionIterator< TImageType >;

  //! Default Constructor
  FeatureBase() {};

  //! Default Destructor
  ~FeatureBase() {};

  //! Sets the input image for the class
  void SetInputImage(const typename TImageType::Pointer image);

  //! Sets the input mask for the class
  void SetInputMask(const typename TImageType::Pointer image);

  //! Set the non-zero indeces for the input image
  void SetNonZeroIndeces(std::vector< typename TImageType::IndexType > &nonZeroIndeces);

  //! Set the current image starting index
  void SetStartingIndex(const typename TImageType::IndexType index) { m_startIndex = index; };

  //! Enables flag, actual writing is done in inherited classes, dependent on the type of features
  void WriteIntermediateFiles() { m_writeIntermediateFiles = true; };

  //! Returns the features (names and respective values)
  std::map< std::string, double > GetOutput()
  {
    if (!m_algorithmDone)
    {
      Update();
    }
    return m_features;
  }

  //! Actual algorithm runner - this is to be changed for each feature family
  virtual void Update() {};

  //! Enable debug mode
  void EnableDebugMode() { m_debugMode = true; };

protected:

  typename TImageType::Pointer m_inputImage, //! the image on which the computations need to happen
    m_Mask; //! the mask on which the computations need to happen
  std::vector< typename TImageType::IndexType > m_nonZeroIndeces; //! non-zero indeces from mask // unused since image is always in mask
  typename TImageType::IndexType m_startIndex; //! starting index of the image (useful for lattice computation)
  bool m_algorithmDone = false; //! whether the algorithm has finished computing or not
  bool m_writeIntermediateFiles = false; //! used for debugging only
  std::map< std::string, double > m_features; //! the output with feature names and their respective values
  bool m_debugMode = false;
};

#include "FeatureBase.hxx"