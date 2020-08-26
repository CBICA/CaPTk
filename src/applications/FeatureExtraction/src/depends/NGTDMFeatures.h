/**
\file  NGTDMFeatures.h

\brief NGTDM feature calculation (Neighbourhood Grey Tone Difference Matrix based features)

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html

*/
#pragma once

#include "itkImage.h"
#include "itkHistogram.h"
#include "itkNumericTraits.h"
#include "itkVectorContainer.h"
#include "itkNeighborhoodIterator.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkMaskImageFilter.h"

#include <map>
#include <string>
#include <numeric>
#include <vector>

#include "FeatureBase.h"
#include "TextureFeatureBase.h"

template< typename TImageType >
class NGTDMFeatures : 
  public FeatureBase < TImageType >, public TextureFeatureBase< TImageType >
{
public:
  //! Default constructor
  NGTDMFeatures() { };

  //! Default destructor
  ~NGTDMFeatures() { };

  /**
  \brief Set the range (how far from the center index of interest do you want to calculate neighborhood tone difference); defaults to 1

  \param rangeValue integer value for the value you want for range
  **/
  void SetRange(int rangeValue)
  {
    this->m_range = rangeValue;
  }

  /**
  \brief Update calculate five feature values
  **/
  void Update()
  {
    if (!this->m_algorithmDone)
    {
      // do the computation that is needed here
      if ((this->m_minimum == 0) || (this->m_maximum == 0))
      {
        auto maskFilter = itk::MaskImageFilter< TImageType, TImageType, TImageType >::New();
        maskFilter->SetInput(this->m_inputImage); //full input image
        maskFilter->SetMaskImage(this->m_Mask);   //Assumed this is single label binary mask (already extracted from multi-label segmentation ROI)
        maskFilter->SetOutsideValue(0);
        maskFilter->Update();
        //Masked out regions outside of the mask (ROI region, single label of multi label ROi).
        //Question: Is the mask being read in only one label of multi-label (e.g. edema of multi label segmentation?)

        auto minMaxComputer = itk::MinimumMaximumImageCalculator< TImageType >::New();
        minMaxComputer->SetImage(maskFilter->GetOutput());
        minMaxComputer->Compute();
        //Calculated min and max of the input image masked out for the given ROI single label

        if (this->m_minimum == 0)
        {
          this->m_minimum = minMaxComputer->GetMinimum();
        }
        if (this->m_maximum == 0)
        {
          this->m_maximum = minMaxComputer->GetMaximum();
        }
      }

      //std::cout << "\n[DEBUG] NGTDMFeatures.h - Update() - this->m_minimum = " << this->m_minimum << std::endl;
      //std::cout << "\n[DEBUG] NGTDMFeatures.h - Update() - this->m_maximum = " << this->m_maximum << std::endl;

      /// histogram calculation from ITK -- for texture feature pipeline
      using TMaskImageType = itk::Image< int, TImageType::ImageDimension >;
      auto caster = itk::CastImageFilter< TImageType, TMaskImageType >::New();
      caster->SetInput(this->m_Mask); //original input binary mask, not the masked image
      caster->Update();

      auto stats = itk::LabelStatisticsImageFilter< TImageType, TMaskImageType >::New();
      stats->SetLabelInput(caster->GetOutput());  //Masked Image
      stats->SetInput(this->m_inputImage);  //input, unmasked full image
      stats->Update();
      stats->SetUseHistograms(true);
      stats->SetHistogramParameters(this->m_Bins, this->m_minimum, this->m_maximum);
      stats->Update();


      this->m_histogram = stats->GetHistogram(1); //Get Histogram for the label value one

      m_radius.Fill(m_range);

      if (this->m_debugMode)
      {
        std::cout << "\n[DEBUG] NGTDMFeatures.h - Update() : m_minimum = " << this->m_minimum << std::endl;
        std::cout << "\n[DEBUG] NGTDMFeatures.h - Update() : m_maximum = " << this->m_maximum << std::endl;
        std::cout << "\n[DEBUG] NGTDMFeatures.h - Update() : itk::LabelStatisticsImageFilter->SetHistogramParameters: (m_Bins=" << this->m_Bins << " | m_minimum=" << this->m_minimum << " | m_maximum=" << this->m_maximum << ")" << std::endl;
        std::cout << "\n[DEBUG] NGTDMFeatures.h - Update() : m_range = " << m_range << std::endl;
      }
      
      using NeighborhoodType = itk::NeighborhoodIterator< TImageType, itk::ConstantBoundaryCondition< TImageType > >;
      NeighborhoodType iter(m_radius, this->m_inputImage, this->m_inputImage->GetBufferedRegion()),
        iterMask(m_radius, this->m_Mask, this->m_Mask->GetBufferedRegion());

      pVector.resize(this->m_Bins, 0);
      sVector.resize(this->m_Bins, 0);

      int count = 0;
      iter.GoToBegin();
      iterMask.GoToBegin();
      while (!iter.IsAtEnd())
      {
        if (iterMask.GetCenterPixel() > 0)
        {
          int localCount = 0;
          double localMean = 0;
          unsigned int localIndex = this->IntensityToIndex(iter.GetCenterPixel());
          for (itk::SizeValueType i = 0; i < iter.Size(); ++i)
          {
            if (i == (iter.Size() / 2))
              continue;
            if (iterMask.GetPixel(i) > 0)
            {
              ++localCount;
              localMean += this->IntensityToIndex(iter.GetPixel(i)) + 1;
            }
          }
          if (localCount > 0)
          {
            localMean /= localCount;
          }
          localMean = std::abs<double>(localIndex + 1 - localMean);

          pVector[localIndex] += 1;
          sVector[localIndex] += localMean;
          ++count;

        }
        ++iterMask;
        ++iter;
      }

      unsigned int Ngp = 0;
      for (unsigned int i = 0; i < this->m_Bins; ++i)
      {
        if (pVector[i] > 0.1)
        {
          ++Ngp;
        }
        pVector[i] /= count;
      }

      double sumS = 0;
      double sumStimesP = 0;

      double contrastA = 0;
      double busynessA = 0;
      double complexity = 0;
      double strengthA = 0;

      for (unsigned int i = 0; i < this->m_Bins; ++i)
      {
        sumS += sVector[i];
        sumStimesP += pVector[i] * sVector[i];

        //TBD - for debugging NGTDM matrix - values were correct for phantom as of 2019-04-26, do not uncomment unless debugging for more complex test cases
        // if (this->m_debugMode){
        //   std::cout << "\n[DEBUG] NGTDMFeatures.h - Update() - bin[" << i << "] | p[" << i << "] * s[" << i << "] = " << pVector[i] << " * " << sVector[i] << " = " << (pVector[i] * sVector[i]) << std::endl;
        // }
        //TBD - for debugging NGTDM matrix

        for (unsigned int j = 0; j < this->m_Bins; ++j)
        {
          double iMinusj = 1.0*i - 1.0*j;
          contrastA += pVector[i] * pVector[j] * iMinusj*iMinusj;
          if ((pVector[i] > 0) && (pVector[j] > 0))
          {
            busynessA += std::abs<double>((i + 1.0)*pVector[i] - (j + 1.0)*pVector[j]);
            complexity += std::abs<double>(iMinusj)*(pVector[i] * sVector[i] + pVector[j] * sVector[j]) / (pVector[i] + pVector[j]);
            strengthA += (pVector[i] + pVector[j])*iMinusj*iMinusj;
          }
        }
      }
      double coarsness_temp = 1.0 / sumStimesP;
      double coarsness = std::min<double>(coarsness_temp, 1000000); // SH: The MAX possible courseness is 10^6

      // SH - orignilly: double coarsness = 1.0 / std::min<double>(sumStimesP, 1000000); // TBD (SP): does it make sense to use "std::numeric_limits< double >::max()" here?
      double contrast = 0;
      double busyness = 0;
      if (Ngp > 1)
      {
        contrast = contrastA / Ngp / (Ngp - 1) / count * sumS;
        busyness = sumStimesP / busynessA;
      }
      complexity /= count;
      double strength = 0;
      if (sumS > 0)
      {
        strength = strengthA / sumS;
      }


      this->m_features["Coarsness"] = coarsness;
      this->m_features["Contrast"] = contrast;
      this->m_features["Busyness"] = busyness;
      this->m_features["Complexity"] = complexity;
      this->m_features["Strength"] = strength;

      //TBD - for debugging NGTDM matrix - values were correct for phantom as of 2019-04-26, do not uncomment unless debugging for more complex test cases
      // if (this->m_debugMode){
      //   std::cout << "\n[DEBUG] NGTDMFeatures.h - Update() - Coarsness = " << coarsness << std::endl;
      //   std::cout << "\n[DEBUG] NGTDMFeatures.h - Update() - Contrast = " << contrast << std::endl;
      //   std::cout << "\n[DEBUG] NGTDMFeatures.h - Update() - Busyness = " << busyness << std::endl;
      //   std::cout << "\n[DEBUG] NGTDMFeatures.h - Update() - Complexity = " << complexity << std::endl;
      //   std::cout << "\n[DEBUG] NGTDMFeatures.h - Update() - Strength = " << strength << std::endl;
      // }

      this->m_algorithmDone = true;
    }
  };

private:

  /**
  \brief IntensityQuantifier::IntensityToIndex(double intensity)
  **/
  unsigned int IntensityToIndex(double intensity)
  {
    itk::Statistics::Histogram< double >::IndexType temp_idx(1);

    itk::Statistics::Histogram< double >::MeasurementVectorType temp_mv(1);
    temp_mv.Fill(intensity);

    if (!this->m_histogram->GetIndex(temp_mv, temp_idx))
    {
      //std::cerr << "Couldn't find index for intensity value '" << intensity << "' in histogram for NGTDM.\n";
    }
    {
      double index = temp_idx[0]; // for uniform histogram, use "std::floor((intensity - this->m_minimum) / this->m_Bins)"
      return std::max<double>(0, std::min<double>(index, this->m_Bins - 1));
    }
  }

  //using HistogramType = THistogramFrequencyContainer;

  unsigned int m_range = 1;

  typename TImageType::SizeType m_radius;

  std::vector< double > pVector;
  std::vector< double > sVector;

};
