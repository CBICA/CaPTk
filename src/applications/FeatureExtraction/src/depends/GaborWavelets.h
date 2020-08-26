/**
\file GaborWavelets.h

This file calculates the Gabor Wavelet features of an input image and mask.

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html
*/
#pragma once

#include <typeinfo>
#include <chrono>
#include <map>
#include <cmath>
#include <math.h>

#include "itkImage.h"
#include "itkNeighborhoodIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkMirrorPadImageFilter.h"
#include "itkStatisticsImageFilter.h"

#include "itkOpenCVImageBridge.h"

#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"
#include "cbicaLogging.h"

#include "FeatureBase.h"

template< class TImageType = itk::Image< float, 3 > >
class GaborWavelets : public FeatureBase < TImageType >
{
public:
  //! Default Constructor
  GaborWavelets() {};

  //! Default Destructor
  ~GaborWavelets() {};

  //! Set the radius for computation
  void SetRadius(int radius) { m_radius = radius; };
  void SetRadius(float radius) { m_radius_float = radius; };

  //! Set the directions for computation
  void SetDirections(int direction) { m_Direction = direction; };

  //! Set the level for computation
  void SetLevel(int level) { m_gaborLevel = level; };

  //! Set the gamma for computation
  void SetGamma(float gamma) { m_gaborGamma = gamma; };

  //! Set the FMax for computation
  void SetFMax(float fmax) { m_gaborFMax = fmax; };

  //! Actual algorithm runner
  void Update()
  {
    if (!this->m_algorithmDone)
    {
      if (TImageType::ImageDimension == 2)
      {
        auto cvImage_input = itk::OpenCVImageBridge::ITKImageToCVMat< TImageType >(this->m_inputImage);

        size_t count = 0;
        // finalize radius if it has been defined in world coordinates
        if (m_radius == -1)
        {
          auto spacing = this->m_inputImage->GetSpacing();
          for (size_t d = 0; d < TImageType::ImageDimension; d++)
          {
            auto temp = m_radius_float / spacing[d];
            if ((temp < 1) && (temp > 0)) // this is a contingency in cases where the radius has been initialized to be less than the pixel spacing
            {
              m_radius = 1;
            }
            else if (m_radius < static_cast< int >(temp))
            {
              m_radius = std::round(temp);
            }
            else
            {
              m_radius = std::round(temp);
            }
          }
        }
        for (int u = 1; u < m_gaborLevel + 1; u++)
        {
          for (int v = 1; v < m_Direction + 1; v++)
          {
            //auto gaborImage = GaborConstruction(m_radius, m_gaborFMax, m_gaborGamma, u, v);

            //typename FeatureBase< TImageType >::TConstIteratorType inputIterator(m_inputImage, m_inputImage->GetBufferedRegion());
            //TRegionIteratorType gaborIterator(gaborImage, gaborImage->GetBufferedRegion());

            //for (gaborIterator.GoToBegin(); !gaborIterator.IsAtEnd(); ++gaborIterator)
            //{
            //  double sum = 0;
            //  auto gaborValue = gaborIterator.Get();
            //  for (inputIterator.GoToBegin(); !inputIterator.IsAtEnd(); ++inputIterator)
            //  {
            //    sum += gaborValue + inputIterator.Get();
            //  }
            //  gaborIterator.Set(sum);
            //}

            auto result = gabor_wavelet(m_radius, m_gaborFMax, m_gaborGamma, u, v);
            auto size = this->m_inputImage->GetLargestPossibleRegion().GetSize();
            //result image
            auto result_matrix = TImageType::New();
            typename TImageType::RegionType region_new;
            typename TImageType::IndexType start;
            start.Fill(0);
            //try to ignore the boundary effect, we can do the convolution without padding the image
            typename TImageType::SizeType size_new_matrix;
            size_new_matrix[0] = size[0] - m_radius + 1;
            size_new_matrix[1] = size[1] - m_radius + 1;
            region_new.SetSize(size_new_matrix);
            result_matrix->SetRegions(region_new);
            result_matrix->Allocate();
            result_matrix->FillBuffer(itk::NumericTraits< double >::Zero);

            std::vector< std::vector< double > > image_matrix(size[1]);

            for (uint j = 0; j < size[1]; j++)
            {
              image_matrix[j].resize(size[0]);
              for (uint i = 0; i < size[0]; i++)
              {
                image_matrix[j][i] = cvImage_input.template ptr< typename TImageType::PixelType >(i)[j];
              }
            }

            typename TImageType::IndexType pixelIndex;

            for (uint i = 0; i < size_new_matrix[1]; i++)
            {
              for (uint j = 0; j < size_new_matrix[0]; j++)
              {
                double sum = 0;
                for (int k = 0; k < m_radius; k++)
                {
                  for (int l = 0; l < m_radius; l++)
                  {
                    sum = sum + result[k][l] * image_matrix[i + k][j + l];
                  }
                }
                pixelIndex[0] = j;
                pixelIndex[1] = i;
                result_matrix->SetPixel(pixelIndex, sum);
              }
            }

            auto statisticsImageFilter = itk::StatisticsImageFilter< TImageType >::New();
            statisticsImageFilter->SetInput(/*gaborImage*/result_matrix);
            statisticsImageFilter->Update();
            //lattice_features.push_back(statisticsImageFilter->GetMean());
            // TBD: check what is happening here
            this->m_features["Scale_" + std::to_string(u) + "_Orientation_" + std::to_string(v) + "_Mean"] = statisticsImageFilter->GetMean();
            this->m_features["Scale_" + std::to_string(u) + "_Orientation_" + std::to_string(v) + "_STD"] = statisticsImageFilter->GetSigma();
            this->m_features["Scale_" + std::to_string(u) + "_Orientation_" + std::to_string(v) + "_Variance"] = statisticsImageFilter->GetVariance();
            this->m_features["Scale_" + std::to_string(u) + "_Orientation_" + std::to_string(v) + "_Max"] = statisticsImageFilter->GetMaximum();
            this->m_features["Scale_" + std::to_string(u) + "_Orientation_" + std::to_string(v) + "_Sum"] = statisticsImageFilter->GetSum();
            count++;
          }
        }
      }

      this->m_algorithmDone = true;
    }
  }

private:

  /**
  \brief Compute the Gabor wavelet image
  */
  //typename TImageType::Pointer GaborConstruction(double R0, double fmax0, double gamma_gabor0, int level0, int orientation0)
  //{
  //  double fu, alpha, tetav;
  //  fu = fmax0 / std::pow(gamma_gabor0, level0 - 1);
  //  alpha = fu / gamma_gabor0;
  //  tetav = (static_cast<double>(orientation0) - 1) / 8 * PI;
  //  auto tetav_sin = std::sin(tetav);
  //  auto tetav_cos = std::cos(tetav);
  //  double xprime, yprime;

  //  auto returnImage = TImageType::New();
  //  typename TImageType::SizeType size;
  //  size.Fill(static_cast< size_t >(R0)); // ensure a static cast so that there is no buffer overflow
  //  typename TImageType::RegionType region(size);
  //  returnImage->SetRegions(region);
  //  returnImage->Allocate();
  //  returnImage->FillBuffer(0);
  //  //vector<vector <double>> realGW;
  //  //vector <double> realtemp;

  //  TRegionIteratorType iterator(returnImage, returnImage->GetBufferedRegion());
  //  iterator.GoToBegin();
  //  for (double i = -R0 / 2 + 1; i <= R0 / 2; i++)
  //  {
  //    //realtemp.clear();
  //    double realpart/*, imagpart*/;
  //    for (double j = -R0 / 2 + 1; j <= R0 / 2; j++)
  //    {
  //      xprime = (i - 0.5)*tetav_cos + (j - 0.5)*tetav_sin;
  //      yprime = -(i - 0.5)*tetav_sin + (j - 0.5)*tetav_cos;
  //      realpart = fu * fu /
  //        (PI * std::pow(gamma_gabor0, 2)) * std::exp(-alpha * alpha * (std::pow(xprime, 2) + std::pow(yprime, 2))) *
  //        std::cos(2 * PI*fu*xprime);
  //      //realpart = fu * fu / (PI*gamma_gabor0*gamma_gabor0)*exp(-alpha * alpha*(xprime*xprime + yprime * yprime))*cos(2 * PI*fu*xprime);
  //      //imagpart = fu * fu /                                                              
  //      //  (itk::Math::pi * std::pow(gamma_gabor0, 2)) * std::exp(-alpha * alpha * (std::pow(xprime, 2) + std::pow(yprime, 2))) * 
  //      //  std::sin(2 * itk::Math::pi*fu*xprime);
  //      //realtemp.push_back(realpart);
  //      iterator.Set(realpart);
  //      ++iterator;
  //    }
  //    //realGW.push_back(realtemp);
  //  }
  //  //return realGW;
  //  return returnImage;
  //}
  std::vector< std::vector< double > > gabor_wavelet(double R0, double fmax0, double gamma_gabor0, int level0, int orientation0)
  {
    double fu, alpha, tetav;
    fu = fmax0 / pow(gamma_gabor0, level0 - 1);
    alpha = fu / gamma_gabor0;
    tetav = (orientation0 - 1 + 0.0) / 8 * PI;
    auto tetav_cos = std::cos(tetav);
    auto tetav_sin = std::sin(tetav);
    double xprime, yprime;
    std::vector< std::vector< double > > realGW;
    std::vector <double> realtemp;
    for (double i = -R0 / 2 + 1; i <= R0 / 2; i++)
    {
      realtemp.clear();
      double realpart;
      for (double j = -R0 / 2 + 1; j <= R0 / 2; j++)
      {
        xprime = (i - 0.5)*tetav_cos + (j - 0.5)*tetav_sin;
        yprime = -(i - 0.5)*tetav_sin + (j - 0.5)*tetav_cos;
        realpart = fu * fu / (PI*gamma_gabor0*gamma_gabor0)*std::exp(-alpha * alpha*(xprime*xprime + yprime * yprime))*std::cos(2 * PI*fu*xprime);
        //imagpart = fu * fu / (PI*gamma_gabor0*gamma_gabor0)*exp(-alpha * alpha*(xprime*xprime + yprime * yprime))*sin(2 * PI*fu*xprime);
        realtemp.push_back(realpart);
      }
      realGW.push_back(realtemp);;
    }
    return realGW;
  }

  /**
  \brief Computes the log of a number by considering a special case for '0'

  When fitting the regression line, we have log function on the expression, the log is defined as this in order to avoid inf value.
  */
  double log_with_zero(double x)
  {
    if (x == 0)
      return -8;
    else
      return std::log(x);
  }

  int m_radius = -1; //! radius around which features are to be extracted
  float m_radius_float = -1; //! radius around which features are to be extracted
  int m_Direction; //! direction around which features are to be extracted
  float m_gaborFMax = 0.25; //! TBD: what is the description of this?
  float m_gaborGamma = sqrtf(2); //! TBD: what is the description of this?
  int m_gaborLevel = 4; //! TBD: what is the description of this?
};