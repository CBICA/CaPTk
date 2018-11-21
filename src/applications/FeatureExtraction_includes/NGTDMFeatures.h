/**
\file  NGTDMFeatures.h

\brief NGTDM feature calculation

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

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

template< typename TImageType >
class NGTDMFeatures : public FeatureBase < TImageType >
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
    m_range = rangeValue;
  }

  /**
  \brief Get the range (how far from the center index of interest do you want to calculate neighborhood tone difference); returns int value
  **/
  int GetRange() {
    return m_range;
  }

  /**
  \brief Get the number of bins for quantization of the input image
  **/
  int GetNumBins() {
    return m_bins;
  }

  /**
  \brief Set the number of bins for quantization of the input image; defaults to 10 unless overridden by user

  \param numBinValue integer value for the number of bins you want for image quantization
  **/
  void SetNumBins(int numBinValue)
  {
    m_bins = numBinValue;
  }

  /**
  \brief Set the minimum
  */
  void SetMinimum(int minimumInput)
  {
    m_minimum = minimumInput;
  }

  /**
  \brief Set the maximum
  */
  void SetMaximum(int maximumInput)
  {
    m_maximum = maximumInput;
  }

  /**
  \brief Update calculate five feature values
  **/
  void Update()
  {
    if (!this->m_algorithmDone)
    {
      // do the computation that is needed here
      if ((m_minimum == 0) || (m_maximum == 0))
      {
        auto maskFilter = itk::MaskImageFilter< TImageType, TImageType, TImageType >::New();
        maskFilter->SetInput(this->m_inputImage);
        maskFilter->SetMaskImage(this->m_Mask);
        maskFilter->SetOutsideValue(0);
        maskFilter->Update();

        auto minMaxComputer = itk::MinimumMaximumImageCalculator< TImageType >::New();
        minMaxComputer->SetImage(maskFilter->GetOutput());
        minMaxComputer->Compute();

        if (m_minimum == 0)
        {
          m_minimum = minMaxComputer->GetMinimum();
        }
        if (m_maximum == 0)
        {
          m_maximum = minMaxComputer->GetMaximum();
        }
      }

      /// histogram calculation from ITK -- for texture feature pipeline
      using TMaskImageType = itk::Image< int, TImageType::ImageDimension >;
      auto caster = itk::CastImageFilter< TImageType, TMaskImageType >::New();
      caster->SetInput(this->m_Mask);
      caster->Update();

      auto stats = itk::LabelStatisticsImageFilter< TImageType, TMaskImageType >::New();
      stats->SetLabelInput(caster->GetOutput());
      stats->SetInput(this->m_inputImage);
      stats->Update();
      stats->SetUseHistograms(true);
      stats->SetHistogramParameters(m_bins, m_minimum, m_maximum);
      stats->Update();

      m_histogram = stats->GetHistogram(1);

      m_radius.Fill(m_range);

      using NeighborhoodType = itk::NeighborhoodIterator< TImageType, itk::ConstantBoundaryCondition< TImageType > >;
      //using NeighborhoodType = itk::NeighborhoodIterator< TImageType >;
      NeighborhoodType iter(m_radius, this->m_inputImage, this->m_inputImage->GetBufferedRegion()),
        iterMask(m_radius, this->m_Mask, this->m_Mask->GetBufferedRegion());

      pVector.resize(m_bins, 0);
      sVector.resize(m_bins, 0);

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
      for (unsigned int i = 0; i < m_bins; ++i)
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

      for (unsigned int i = 0; i < m_bins; ++i)
      {
        sumS += sVector[i];
        sumStimesP += pVector[i] * sVector[i];
        //TBD - for debugging NGTDM matrix
        //std::cout << "\n bin[" << i << "] | p[" << i << "] * s[" << i << "] = " << pVector[i] << " * " << sVector[i] << " = " << (pVector[i] * sVector[i]) << std::endl;
        //TBD - for debugging NGTDM matrix
        for (unsigned int j = 0; j < m_bins; ++j)
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

      //TBD - for debugging NGTDM features
      //std::cout << "\n coursness_temp = " << coarsness_temp << std::endl;
      //TBD - for debugging NGTDM features

      double coarsness = std::min<double>(coarsness_temp, 1000000); // SH: The MAX possible courseness is 10^6

      //TBD - for debugging NGTDM matrix
      //std::cout << "\n coarsness = " << coarsness << std::endl;
      //TBD - for debugging NGTDM matrix

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

      this->m_algorithmDone = true;
    }
  };

  /**
  \brief return the map of feature names and feature values
  **/
  std::map< std::string, double > GetOutput()
  {
    if (!this->m_algorithmDone)
    {
      Update();
    }
    return this->m_features;
  };

  /**
  \brief return the Coarsness
  **/
  double GetCoarsness()
  {
    if (!this->m_algorithmDone)
    {
      Update();
    }
    return this->m_features["Coarsness"];
  }

  /**
  \brief return the Contrast
  **/
  double GetContrast()
  {
    if (!this->m_algorithmDone)
    {
      Update();
    }
    return this->m_features["Contrast"];
  }

  /**
  \brief return the Busyness
  **/
  double GetBusyness()
  {
    if (!this->m_algorithmDone)
    {
      Update();
    }
    return this->m_features["Busyness"];
  }

  /**
  \brief return the Complexity
  **/
  double GetComplexity()
  {
    if (!this->m_algorithmDone)
    {
      Update();
    }
    return this->m_features["Complexity"];
  }

  /**
  \brief return the Strength
  **/
  double GetStrength()
  {
    if (!this->m_algorithmDone)
    {
      Update();
    }
    return this->m_features["Strength"];
  }

private:

  /**
  \brief IntensityQuantifier::IntensityToIndex(double intensity)
  **/
  unsigned int IntensityToIndex(double intensity)
  {
    itk::Statistics::Histogram< double >::IndexType temp_idx(1);

    itk::Statistics::Histogram< double >::MeasurementVectorType temp_mv(1);
    temp_mv.Fill(intensity);

    if (!m_histogram->GetIndex(temp_mv, temp_idx))
    {
      std::cerr << "Couldn't find index for intensity value '" << intensity << "' in histogram for NGTDM.\n";
    }
    {
      double index = temp_idx[0]; // for uniform histogram, use "std::floor((intensity - m_minimum) / m_bins)"
      return std::max<double>(0, std::min<double>(index, m_bins - 1));
    }
  }

  //using HistogramType = THistogramFrequencyContainer;

  unsigned int m_range = 1;
  unsigned int m_bins = 10;
  int m_minimum = 0;
  int m_maximum = 0;

  typename TImageType::SizeType m_radius;

  std::vector< double > pVector;
  std::vector< double > sVector;

  itk::Statistics::Histogram< double >::Pointer m_histogram;
};
