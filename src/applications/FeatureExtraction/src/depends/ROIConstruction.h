/**
\file ROIConstruction.h

This file holds the declaration of the class ROIConstruction.

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html
*/
#pragma once

#include <typeinfo>
#include <chrono>

#include "itkImage.h"
#include "itkNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"

#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"
#include "cbicaLogging.h"

//using DefaultImageType = itk::Image< float, 3 >;

/**
\brief This class constructs the ROIs based on the parameters provided

If lattice window size & step are provided, then lattice points and patches are also computed.
*/
template< class TImageType = itk::Image< float, 3 > >
class ROIConstruction
{
public:
  //! Some useful aliases
  using TIteratorType = itk::ImageRegionIteratorWithIndex< TImageType >;
  using TConstIteratorType = itk::ImageRegionConstIterator< TImageType >;

  /**
  \brief Helper structure to store the ROI properties for Feature Extraction
  */
  struct ROIProperties
  {
    int value; //! the value of the ROI coming from the complete mask image
    std::string label; //! label corresponding to above value
    float weight = 1.0; //! the weight assigned to the ROI based on how much of the mask is present
    std::vector< typename TImageType::IndexType > nonZeroIndeces; // store the ROI indeces for computation of intensity features
    typename TImageType::IndexType centerIndex; // save the center index for lattice computation
    bool latticeGridPoint = false; //! whether the current ROI is part of a lattice grid point or not
  };

  //! Default Constructor
  ROIConstruction() {};

  //! Default Destructor
  ~ROIConstruction() {};

  //! Get Output
  std::vector< ROIProperties > GetOutput() 
  {
    if (!m_algorithmDone)
    {
      Update();
    }
    return m_output;
  };

  //! Run
  void Update();

  //! Set the input mask
  void SetInputMask(typename TImageType::Pointer input);

  //! Set the input images
  //void SetInputImages(std::vector< typename TImageType::Pointer > images);

  /*
  \brief Set the requested ROIs and Labels

  \param roi The vector of ROIs as integers (computed in FeatureExtraction class)
  \param roi_labels The vector of ROI labels as strings
  */
  void SetSelectedROIsAndLabels(std::vector< int > roi, std::vector< std::string > roi_labels);

  /**
  \brief Sets the lattice window size (radius around each grid node)

  \param window Input Lattice Window in mm
  */
  void SetLatticeWindowSize(float window) { m_latticeWindow = window; };

  /**
  \brief This sets whether the fluxNeumann Boundary condition is to be used or not.

  Ref: https://itk.org/Doxygen/html/classitk_1_1ZeroFluxNeumannBoundaryCondition.html
  */
  void SetBoundaryCondition(bool fluxNeumannCondition) { m_fluxNeumannEnabled = fluxNeumannCondition; };

  /**
  \brief Sets whether the entire patch is taken into consideration or just the part in the selected ROI
  */
  void SetPatchConstructionConditionROI(bool patchConstructionROI) { m_patchOnRoiEnabled = patchConstructionROI; };

  /**
  \brief Sets whether the patch gets diregarded when the weight is less than 1
  */
  void SetPatchConstructionConditionNone(bool patchConstructionNone) { m_patchBoundaryDisregarded = patchConstructionNone; };

  /**
  \brief Sets the lattice grid step (distance between 2 grid nodes)

  \param step Input Lattice Grid Step in mm
  */
  void SetLatticeGridStep(float step) 
  { 
    m_latticeStep = step; 
    if (!m_mask.IsNull())
    {
      auto tempSpacing = m_mask->GetSpacing();
      for (size_t d = 0; d < TImageType::ImageDimension; d++)
      {
        auto temp = m_latticeStep / tempSpacing[d];
        if (temp < 1) // condition where the spacing is larger than the lattice step
        {
          m_latticeStepImage[d] = 1;
        }
        else
        {
          m_latticeStepImage[d] = std::round(temp);
        }
      }
    }
  };

  //! Set new logger file
  void SetNewLogFile(const std::string &logFile);

  //! Check if Lattice computation has been enabled or not
  bool IsLatticeEnabled()
  {
    if (!m_algorithmDone)
    {
      Update();
    }
    return m_latticeEnabled;
  }

  //! Get the number of valid lattice points (as defined by the properties in the parameter file)
  size_t GetNumberOfValidLatticePoints() { return m_latticePoints; };

  //! Get the lattice patch radius in image coordinates
  itk::Size< TImageType::ImageDimension > GetLatticeRadius() { return m_latticeRadius; };

  //! Get the lattice step in image coordinates
  typename TImageType::IndexType GetLatticeStepImage() { return m_latticeStepImage; };

private:

  /**
  \brief Populate the output for non-lattice computation

  \param inputLatticeGrid The lattice grid from which patches are to be extracted; defined in world coordiantes
  \param inputLatticeRadius The lattice patch radius; defined in image coordiantes
  */
  void ROIComputation(std::vector< typename TImageType::IndexType > &inputLatticeGrid, itk::Size< TImageType::ImageDimension > &inputLatticeRadius);

  size_t m_latticePoints = 0; // number of valid lattice points
  bool m_latticeEnabled = false; //! flag to check if lattice-based computation has been enabled or not
  itk::Size< TImageType::ImageDimension > m_latticeRadius; //! lattice radius in image space
  bool m_algorithmDone = false; // check if processing is finished or not
  typename TImageType::Pointer m_mask; //! the full mask image
  std::vector< int > m_roi; //! the ROI values on which calculation needs to happen
  std::vector< std::string > m_roiLabels; //! ROI label names 
  float m_latticeWindow = 0, m_latticeStep = 0; //! these are defined in mm
  typename TImageType::IndexType m_latticeStepImage; //! lattice step as defined in the image coordinates
  bool m_fluxNeumannEnabled = false, m_zeroPaddingEnabled = false; //! zero-padding is the default behavior for boundary condition
  bool m_patchOnRoiEnabled = false; //! whether to pull the entire patch or only along the ROI
  bool m_patchBoundaryDisregarded = false; //! only considers patches with all pixels != 0
  bool m_maskBackgroundWeighting = false; //! whether the mask background needs to be weighted or not
  std::vector< ROIProperties > m_output; //! the output of the class
  cbica::Logging m_logger;
};

#include "ROIConstruction.hxx"