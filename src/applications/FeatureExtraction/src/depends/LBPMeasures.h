/**
\file LBPMeasures.h

This file calculates the LBP Features of an input image and mask.

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html
*/
#pragma once

#include <typeinfo>
#include <chrono>
#include <map>
#include <cmath>

#include "itkImage.h"
#include "itkNeighborhoodIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkMirrorPadImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkPowImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkMaskImageFilter.h"
#include "itkImageDuplicator.h"

#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"
#include "cbicaLogging.h"

#include "FeatureBase.h"

template< class TImageType = itk::Image< float, 3 > >
class LBPMeasures : public FeatureBase < TImageType >
{
  using TConstIteratorType = itk::ImageRegionConstIterator< TImageType >;
  using TConstNeighborhoodIteratorType = itk::NeighborhoodIterator< TImageType >;
public:
  //! Default Constructor
  LBPMeasures() {};

  //! Default Destructor
  ~LBPMeasures() {};

  //! Set the radius for computation
  void SetRadius(int radius) { m_radius = radius; };
  
  //! Set the radius for computation in world coordinates
  void SetRadius(float radius) { m_radius_float = radius; };

  //! Set number of neighbors to consider for computation
  void SetNeighbors(size_t neighbors) { m_neighbors = neighbors; };

  //! Set the LBP style for computation
  void SetLBPStyle(int style) { m_style = style; };

  //! Actual algorithm runner
  void Update()
  {
    if (!this->m_algorithmDone)
    {
      typename TConstNeighborhoodIteratorType::RadiusType radius;
      // finalize radius if it has been defined in world coordinates
      if (m_radius == -1)
      {
        auto spacing = this->m_inputImage->GetSpacing();
        for (size_t d = 0; d < TImageType::ImageDimension; d++)
        {
          auto temp = m_radius_float / spacing[d];
          if ((temp < 1) && (temp > 0)) // this is a contingency in cases where the radius has been initialized to be less than the pixel spacing
          {
            radius[d] = 1;
          }
          else
          {
            radius[d] = std::round(temp);
          }
        }
      }
      else
      {
        for (unsigned int i = 0; i < TImageType::ImageDimension; i++)
          radius[i] = m_radius;
      }

      // mask the input mask
      auto masker = itk::MaskImageFilter< TImageType, TImageType >::New();
      masker->SetInput(this->m_inputImage);
      masker->SetMaskImage(this->m_Mask);
      masker->SetOutsideValue(0);
      masker->Update();

      auto maskedInput = masker->GetOutput();

      auto outputImage = cbica::CreateImage< TImageType >(maskedInput);

      TConstNeighborhoodIteratorType ImageIterator(radius, maskedInput, maskedInput->GetLargestPossibleRegion());
      itk::ImageRegionIterator< TImageType > outputIterator(outputImage, outputImage->GetLargestPossibleRegion());
      
      for (ImageIterator.GoToBegin(); !ImageIterator.IsAtEnd(); ++ImageIterator)
      {
        auto currentindex = ImageIterator.GetIndex();
        float CenterIntensityValue = ImageIterator.GetCenterPixel();

        // make parallel here, perhaps?
        double lbp_pattern = 0;
        int neighborhoodNumber = 0;
        for (unsigned int LocalNeighborhoodIterator = 0; LocalNeighborhoodIterator < ImageIterator.Size(); ++LocalNeighborhoodIterator)
        {
          auto offsetType1 = ImageIterator.ComputeInternalIndex(LocalNeighborhoodIterator);
          float NeighborIntensityValue = ImageIterator.GetPixel(LocalNeighborhoodIterator);

          lbp_pattern += (NeighborIntensityValue > CenterIntensityValue) << neighborhoodNumber;
          neighborhoodNumber++;
        }

        outputIterator.SetIndex(currentindex);
        outputIterator.Set(lbp_pattern);
      }


      itk::ImageRegionConstIterator< TImageType > outputConstIterator(outputImage, outputImage->GetLargestPossibleRegion());
      double sum = 0;
      int counter = 0;
      for (outputConstIterator.GoToBegin(); !outputConstIterator.IsAtEnd(); ++outputConstIterator)
      {
        sum += outputConstIterator.Get();
        counter++;
      }

      double mean_lbp_value = sum / counter;
      this->m_features["LBP"] = mean_lbp_value;
 
      this->m_algorithmDone = true;

    }
  }

private:

  int m_radius = -1; //! radius around which features are to be extracted
  float m_radius_float = -1; //! radius around which features are to be extracted
  size_t m_neighbors; //! neighbors for LBP computation
  int m_style = 2; //!case 0: original LBP; case 1: uniform LBP; case 2 [DEFAULT]: rotation invariant LBP; case 3: uniform + rotation invariant LBP
};