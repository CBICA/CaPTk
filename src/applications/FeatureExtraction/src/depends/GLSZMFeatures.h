/**
\file  GLSZMFeatures.h

\brief GLSZM feature calculation

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
#include "itkVariableSizeMatrix.h"
#include "itkOffset.h"

//#include "Eigen/src/Core/Matrix.h"
// #include "Eigen/src/Core/Array.h"
#include "Eigen/Dense"

#include <map>
#include <string>
#include <numeric>
#include <vector>

#include "FeatureBase.h"

//! Helper Class
struct GreyLevelSizeZoneMatrixHolder
{
public:
  GreyLevelSizeZoneMatrixHolder(double min, double max, int number, int maxSize) :
    m_MinimumRange(min),
    m_MaximumRange(max),
    m_NumberOfBins(number),
    m_MaximumSize(maxSize)
  {
    m_Matrix.resize(m_NumberOfBins, m_MaximumSize);
    m_Matrix.fill(0);
    m_Stepsize = (m_MaximumRange - m_MinimumRange) / (m_NumberOfBins);
  }

  int IntensityToIndex(double intensity)
  {
    return std::floor((intensity - m_MinimumRange) / m_Stepsize);
  }

  double IndexToMinIntensity(int index)
  {
    return m_MinimumRange + index * m_Stepsize;
  }

  double IndexToMeanIntensity(int index)
  {
    return m_MinimumRange + (index + 0.5) * m_Stepsize;
  }

  double IndexToMaxIntensity(int index)
  {
    return m_MinimumRange + (index + 1) * m_Stepsize;
  }

  double m_MinimumRange;
  double m_MaximumRange;
  double m_Stepsize;
  int m_NumberOfBins;
  int m_MaximumSize;
  Eigen::MatrixXd m_Matrix;
};

//! Helper Class
struct GreyLevelSizeZoneFeatures
{
  GreyLevelSizeZoneFeatures() :
    SmallZoneEmphasis(0),
    LargeZoneEmphasis(0),
    LowGreyLevelEmphasis(0),
    HighGreyLevelEmphasis(0),
    SmallZoneLowGreyLevelEmphasis(0),
    SmallZoneHighGreyLevelEmphasis(0),
    LargeZoneLowGreyLevelEmphasis(0),
    LargeZoneHighGreyLevelEmphasis(0),
    GreyLevelNonUniformity(0),
    GreyLevelNonUniformityNormalized(0),
    ZoneSizeNonUniformity(0),
    ZoneSizeNoneUniformityNormalized(0),
    ZonePercentage(0),
    GreyLevelMean(0),
    GreyLevelVariance(0),
    ZoneSizeMean(0),
    ZoneSizeVariance(0),
    ZoneSizeEntropy(0)
  {
  }

public:
  double SmallZoneEmphasis;
  double LargeZoneEmphasis;
  double LowGreyLevelEmphasis;
  double HighGreyLevelEmphasis;
  double SmallZoneLowGreyLevelEmphasis;
  double SmallZoneHighGreyLevelEmphasis;
  double LargeZoneLowGreyLevelEmphasis;
  double LargeZoneHighGreyLevelEmphasis;
  double GreyLevelNonUniformity;
  double GreyLevelNonUniformityNormalized;
  double ZoneSizeNonUniformity;
  double ZoneSizeNoneUniformityNormalized;
  double ZonePercentage;
  double GreyLevelMean;
  double GreyLevelVariance;
  double ZoneSizeMean;
  double ZoneSizeVariance;
  double ZoneSizeEntropy;
};

template< typename TImageType >
class GLSZMFeatures : public FeatureBase < TImageType >
{
public:
  //! Default constructor
  GLSZMFeatures() { };

  //! Default destructor
  ~GLSZMFeatures() { };

  //! Set the radius for computation
  void SetRadius(int radius) { m_radius = radius; };

  //! Set the radius for computation in world coordinates
  void SetRadius(float radius) { m_radius_float = radius; };

  /**
  \brief Set the range (how far from the center index of interest do you want to calculate neighborhood tone difference); defaults to 1

  \param rangeValue integer value for the value you want for range
  **/
  void SetMaxSize(int rangeValue)
  {
    m_maxSize = rangeValue;
  }

  /**
  \brief Get the range (how far from the center index of interest do you want to calculate neighborhood tone difference); returns int value
  **/
  int GetMaxSize()
  {
    return m_maxSize;
  }

  /**
  \brief Get the number of bins for quantization of the input image
  **/
  int GetNumBins()
  {
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
  Set the offsets
  */
  void SetOffsets(itk::VectorContainer< unsigned char, typename TImageType::OffsetType > * offsets)
  {
    m_offsets = offsets;
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

      std::cout << "\n[DEBUG] GLSZMFeatures.h - Update() - m_minimum = " << m_minimum << std::endl;
      std::cout << "\n[DEBUG] GLSZMFeatures.h - Update() - m_maximum = " << m_maximum << std::endl;
      std::cout << "\n[DEBUG] GLSZMFeatures.h - Update() - m_bins = " << m_bins << std::endl;

      GreyLevelSizeZoneMatrixHolder tmpHolder(m_minimum, m_maximum, m_bins, TImageType::ImageDimension);
      int largestRegion = CalculateGlSZMatrix(true, tmpHolder);

      std::cout << "\n[DEBUG] GLSZMFeatures.h - Update() - largestRegion = " << largestRegion << std::endl;

      GreyLevelSizeZoneMatrixHolder holderOverall(m_minimum, m_maximum, m_bins, largestRegion);
      CalculateGlSZMatrix(false, holderOverall);
      CalculateFeatures(holderOverall);

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

  void CalculateFeatures(GreyLevelSizeZoneMatrixHolder &holder)
  {
    auto SgzMatrix = holder.m_Matrix;
    auto pgzMatrix = holder.m_Matrix;
    auto pgMatrix = holder.m_Matrix;
    auto pzMatrix = holder.m_Matrix;

    double Ns = pgzMatrix.sum();
    pgzMatrix /= Ns;
    pgMatrix.rowwise().normalize();
    pzMatrix.colwise().normalize();

    for (int i = 0; i < holder.m_NumberOfBins; ++i)
      for (int j = 0; j < holder.m_NumberOfBins; ++j)
      {
        if (pgzMatrix(i, j) != pgzMatrix(i, j))
          pgzMatrix(i, j) = 0;
        if (pgMatrix(i, j) != pgMatrix(i, j))
          pgMatrix(i, j) = 0;
        if (pzMatrix(i, j) != pzMatrix(i, j))
          pzMatrix(i, j) = 0;
      }

    Eigen::VectorXd SgVector = SgzMatrix.rowwise().sum();
    Eigen::VectorXd SzVector = SgzMatrix.colwise().sum();

    GreyLevelSizeZoneFeatures results;

    for (int j = 0; j < SzVector.size(); ++j)
    {
      results.SmallZoneEmphasis += SzVector(j) / (j + 1) / (j + 1);
      results.LargeZoneEmphasis += SzVector(j) * (j + 1.0) * (j + 1.0);
      results.ZoneSizeNonUniformity += SzVector(j) * SzVector(j);
      results.ZoneSizeNoneUniformityNormalized += SzVector(j) * SzVector(j);
    }
    for (int i = 0; i < SgVector.size(); ++i)
    {
      results.LowGreyLevelEmphasis += SgVector(i) / (i + 1) / (i + 1);
      results.HighGreyLevelEmphasis += SgVector(i) * (i + 1) * (i + 1);
      results.GreyLevelNonUniformity += SgVector(i)*SgVector(i);
      results.GreyLevelNonUniformityNormalized += SgVector(i)*SgVector(i);
    }

    for (int i = 0; i < SgzMatrix.rows(); ++i)
    {
      for (int j = 0; j < SgzMatrix.cols(); ++j)
      {
        results.SmallZoneLowGreyLevelEmphasis += SgzMatrix(i, j) / (i + 1) / (i + 1) / (j + 1) / (j + 1);
        results.SmallZoneHighGreyLevelEmphasis += SgzMatrix(i, j) * (i + 1) * (i + 1) / (j + 1) / (j + 1);
        results.LargeZoneLowGreyLevelEmphasis += SgzMatrix(i, j) / (i + 1) / (i + 1) * (j + 1.0) * (j + 1.0);
        results.LargeZoneHighGreyLevelEmphasis += SgzMatrix(i, j) * (i + 1) * (i + 1) * (j + 1.0) * (j + 1.0);
        results.ZonePercentage += SgzMatrix(i, j)*(j + 1);

        results.GreyLevelMean += (i + 1)*pgzMatrix(i, j);
        results.ZoneSizeMean += (j + 1)*pgzMatrix(i, j);
        if (pgzMatrix(i, j) > 0)
          results.ZoneSizeEntropy -= pgzMatrix(i, j) * std::log(pgzMatrix(i, j)) / std::log(2);
      }
    }

    for (int i = 0; i < SgzMatrix.rows(); ++i)
    {
      for (int j = 0; j < SgzMatrix.cols(); ++j)
      {
        results.GreyLevelVariance += (i + 1 - results.GreyLevelMean)*(i + 1 - results.GreyLevelMean)*pgzMatrix(i, j);
        results.ZoneSizeVariance += (j + 1 - results.ZoneSizeMean)*(j + 1 - results.ZoneSizeMean)*pgzMatrix(i, j);
      }
    }

    results.SmallZoneEmphasis /= Ns;
    results.LargeZoneEmphasis /= Ns;
    results.LowGreyLevelEmphasis /= Ns;
    results.HighGreyLevelEmphasis /= Ns;

    results.SmallZoneLowGreyLevelEmphasis /= Ns;
    results.SmallZoneHighGreyLevelEmphasis /= Ns;
    results.LargeZoneLowGreyLevelEmphasis /= Ns;
    results.LargeZoneHighGreyLevelEmphasis /= Ns;
    results.GreyLevelNonUniformity /= Ns;
    results.GreyLevelNonUniformityNormalized /= Ns * Ns;
    results.ZoneSizeNonUniformity /= Ns;
    results.ZoneSizeNoneUniformityNormalized /= Ns * Ns;
    results.ZonePercentage = Ns / results.ZonePercentage;

    this->m_features["SmallZoneEmphasis"] = results.SmallZoneEmphasis;
    this->m_features["LargeZoneEmphasis"] = results.LargeZoneEmphasis;
    this->m_features["LowGreyLevelEmphasis"] = results.LowGreyLevelEmphasis;
    this->m_features["HighGreyLevelEmphasis"] = results.HighGreyLevelEmphasis;
    this->m_features["SmallZoneLowGreyLevelEmphasis"] = results.SmallZoneLowGreyLevelEmphasis;
    this->m_features["SmallZoneHighGreyLevelEmphasis"] = results.SmallZoneHighGreyLevelEmphasis;
    this->m_features["LargeZoneLowGreyLevelEmphasis"] = results.LargeZoneLowGreyLevelEmphasis;
    this->m_features["LargeZoneHighGreyLevelEmphasis"] = results.LargeZoneHighGreyLevelEmphasis;
    this->m_features["GreyLevelNonUniformity"] = results.GreyLevelNonUniformity;
    this->m_features["GreyLevelNonUniformityNormalized"] = results.GreyLevelNonUniformityNormalized;
    this->m_features["ZoneSizeNonUniformity"] = results.ZoneSizeNonUniformity;
    this->m_features["ZoneSizeNoneUniformityNormalized"] = results.ZoneSizeNoneUniformityNormalized;
    this->m_features["ZonePercentage"] = results.ZonePercentage;
    this->m_features["GreyLevelMean"] = results.GreyLevelMean;
    this->m_features["GreyLevelVariance"] = results.GreyLevelVariance;
    this->m_features["ZoneSizeMean"] = results.ZoneSizeMean;
    this->m_features["ZoneSizeVariance"] = results.ZoneSizeVariance;
    this->m_features["ZoneSizeEntropy"] = results.ZoneSizeEntropy;

    std::cout << "\n[DEBUG] GLSZMFeatures.h - CalculateFeatures - SmallZoneEmphasis = " << results.SmallZoneEmphasis << std::endl;
    std::cout << "\n[DEBUG] GLSZMFeatures.h - CalculateFeatures - LargeZoneEmphasis = " << results.LargeZoneEmphasis << std::endl;
    std::cout << "\n[DEBUG] GLSZMFeatures.h - CalculateFeatures - LowGreyLevelEmphasis = " << results.LowGreyLevelEmphasis << std::endl;
    std::cout << "\n[DEBUG] GLSZMFeatures.h - CalculateFeatures - HighGreyLevelEmphasis = " << results.HighGreyLevelEmphasis << std::endl;
    std::cout << "\n[DEBUG] GLSZMFeatures.h - CalculateFeatures - SmallZoneLowGreyLevelEmphasis = " << results.SmallZoneLowGreyLevelEmphasis << std::endl;
    std::cout << "\n[DEBUG] GLSZMFeatures.h - CalculateFeatures - SmallZoneHighGreyLevelEmphasis = " << results.SmallZoneHighGreyLevelEmphasis << std::endl;
    std::cout << "\n[DEBUG] GLSZMFeatures.h - CalculateFeatures - LargeZoneLowGreyLevelEmphasis = " << results.LargeZoneLowGreyLevelEmphasis << std::endl;
    std::cout << "\n[DEBUG] GLSZMFeatures.h - CalculateFeatures - LargeZoneHighGreyLevelEmphasis = " << results.LargeZoneHighGreyLevelEmphasis << std::endl;
    std::cout << "\n[DEBUG] GLSZMFeatures.h - CalculateFeatures - GreyLevelNonUniformity = " << results.GreyLevelNonUniformity << std::endl;
    std::cout << "\n[DEBUG] GLSZMFeatures.h - CalculateFeatures - GreyLevelNonUniformityNormalized = " << results.GreyLevelNonUniformityNormalized << std::endl;
    std::cout << "\n[DEBUG] GLSZMFeatures.h - CalculateFeatures - ZoneSizeNonUniformity = " << results.ZoneSizeNonUniformity << std::endl;
    std::cout << "\n[DEBUG] GLSZMFeatures.h - CalculateFeatures - ZoneSizeNoneUniformityNormalized = " << results.ZoneSizeNoneUniformityNormalized << std::endl;
    std::cout << "\n[DEBUG] GLSZMFeatures.h - CalculateFeatures - ZonePercentage = " << results.ZonePercentage << std::endl;
    std::cout << "\n[DEBUG] GLSZMFeatures.h - CalculateFeatures - GreyLevelVariance = " << results.GreyLevelVariance << std::endl;
    std::cout << "\n[DEBUG] GLSZMFeatures.h - CalculateFeatures - ZoneSizeMean = " << results.ZoneSizeMean << std::endl;
    std::cout << "\n[DEBUG] GLSZMFeatures.h - CalculateFeatures - ZoneSizeVariance = " << results.ZoneSizeVariance << std::endl;
    std::cout << "\n[DEBUG] GLSZMFeatures.h - CalculateFeatures - ZoneSizeEntropy = " << results.ZoneSizeEntropy << std::endl;
  }

  int CalculateGlSZMatrix(bool estimateLargestRegion, GreyLevelSizeZoneMatrixHolder &holder)
  {
    typedef typename TImageType::IndexType IndexType;

    typedef itk::ImageRegionIteratorWithIndex<TImageType> ConstIterType;
    typedef itk::ImageRegionIteratorWithIndex<TImageType> ConstMaskIterType;

    auto region = this->m_Mask->GetLargestPossibleRegion();
    typename TImageType::RegionType newRegion;
    newRegion.SetSize(region.GetSize());
    newRegion.SetIndex(region.GetIndex());

    ConstIterType imageIter(this->m_inputImage, this->m_inputImage->GetLargestPossibleRegion());
    ConstMaskIterType maskIter(this->m_Mask, this->m_Mask->GetLargestPossibleRegion());

    typename TImageType::Pointer visitedImage = TImageType::New();
    visitedImage->SetRegions(newRegion);
    visitedImage->Allocate();
    visitedImage->FillBuffer(0);

    int largestRegion = 0;

    while (!maskIter.IsAtEnd())
    {
      if (maskIter.Value() > 0)
      {
        auto startIntensityIndex = holder.IntensityToIndex(imageIter.Value());
        std::vector<IndexType> indices;
        indices.push_back(maskIter.GetIndex());
        unsigned int steps = 0;

        while (indices.size() > 0)
        {
          auto currentIndex = indices.back();
          indices.pop_back();

          if (!region.IsInside(currentIndex))
          {
            continue;
          }

          auto wasVisited = visitedImage->GetPixel(currentIndex);
          auto newIntensityIndex = holder.IntensityToIndex(this->m_inputImage->GetPixel(currentIndex));
          auto isInMask = this->m_Mask->GetPixel(currentIndex);

          if ((isInMask > 0) &&
            (newIntensityIndex == startIntensityIndex) &&
            (wasVisited < 1))
          {
            ++steps;
            visitedImage->SetPixel(currentIndex, 1);
            for (size_t i = 0; i < m_offsets->size(); i++)
            {
              auto newIndex = currentIndex + m_offsets->at(i);
              indices.push_back(newIndex);
              newIndex = currentIndex - m_offsets->at(i);
              indices.push_back(newIndex);
            }
          }
        }
        if ((steps > 0) /*&& (startIntensityIndex < holder.m_Matrix.rows())*/)
        {
          //std::cout << "[DEBUG] GLSZMFeatures.h - steps = " << steps << std::endl;
          largestRegion = std::max<int>(steps, largestRegion);
          steps = std::min<unsigned int>(steps, holder.m_MaximumSize);

          if (!estimateLargestRegion)
          {
            if (startIntensityIndex >= holder.m_Matrix.rows())
            {
              startIntensityIndex = holder.m_Matrix.rows() - 1;
            }
            holder.m_Matrix(startIntensityIndex, steps - 1) += 1;
          }

        }
      }
      ++imageIter;
      ++maskIter;
    }

    //TBD - for debugging GLSZM matrix
    if (!estimateLargestRegion){
      std::cout << "\n[DEBUG] NGTDMFeatures.h - CalculateGLSZMatrix() - holder.m_Matrix = \n" << holder.m_Matrix << std::endl;
    }
    //TBD - for debugging GLSZM matrix

    return largestRegion;
  }

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

  unsigned int m_maxSize = 1;
  unsigned int m_bins = 10;
  int m_minimum = 0;
  int m_maximum = 0;

  int m_radius = -1;
  float m_radius_float = -1;

  std::vector< double > pVector;
  std::vector< double > sVector;

  itk::Statistics::Histogram< double >::Pointer m_histogram;

  itk::VectorContainer< unsigned char, typename TImageType::OffsetType > * m_offsets;
};
