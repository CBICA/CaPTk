/**
\file EdgeEnhancement.h

This file calculates the Laws Measures Features of an input image and mask.

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

#include "itkImage.h"
#include "itkNeighborhoodIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkMirrorPadImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkPowImageFilter.h"

#include "itkOpenCVImageBridge.h"

#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"
#include "cbicaLogging.h"
#include "cbicaStatistics.h"

#include "FeatureBase.h"

template< class TImageType = itk::Image< float, 3 > >
class EdgeEnhancement : public FeatureBase < TImageType >
{
public:
  //! Default Constructor
  EdgeEnhancement() {};

  //! Default Destructor
  ~EdgeEnhancement() {};

  //! Sets the ETA
  void SetETA(float eta) { m_Eta = eta; };

  //! Sets the Epsilon
  void SetEpsilon(float eps) { m_epsilon = eps; };

  //! Set the radius for computation
  void SetRadius(int radius) { m_radius = radius; };
  void SetRadius(float radius) { m_radius_float = radius; };

  //! Actual algorithm runner
  void Update()
  {
    if (!this->m_algorithmDone)
    {
      auto imageSize = this->m_inputImage->GetLargestPossibleRegion().GetSize();

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

      auto cvImage_input = itk::OpenCVImageBridge::ITKImageToCVMat< TImageType >(this->m_inputImage);

      /// this is only valid for 2D images... need to generalize for 3D
      if (TImageType::ImageDimension == 2)
      {
        std::vector< std::vector< double > > image_matrix(imageSize[1]);

        for (uint j = 0; j < imageSize[1]; j++)
        {
          image_matrix[j].resize(imageSize[0]);
          for (uint i = 0; i < imageSize[0]; i++)
          {
            image_matrix[j][i] = cvImage_input.template ptr< typename TImageType::PixelType >(i)[j];
          }
        }

        //for (int j = 0; j < imageSize[1]; j++)
        //{
        //  image_matrix[j].resize(imageSize[0]);
        //  for (int i = 0; i < imageSize[0]; i++)
        //  {
        //    typename TImageType::IndexType index;
        //    index[0] = i;
        //    index[1] = j;
        //    image_matrix[j][i] = m_inputImage->GetPixel(index);
        //  }
        //}

        //double Eta, eps;
        auto fsx = gD(image_matrix, m_radius, 1, 0);
        auto fsy = gD(image_matrix, m_radius, 0, 1);
        std::vector< double > inputVector(fsx.size() * fsx[0].size());
        for (uint i = 0; i < fsx.size(); i++)
        {
          for (uint j = 0; j < fsx[0].size(); j++)
          {
            double beta2 = std::exp(-(fsx[i][j] * fsx[i][j] + fsy[i][j] * fsy[i][j]) / (m_Eta * m_Eta));
            double beta1 = 0.2 * beta2;
            double factor = std::sqrt(fsx[i][j] * fsx[i][j] + fsy[i][j] * fsy[i][j]);
            double A = (beta1 * fsx[i][j] * fsx[i][j] + beta2 * fsy[i][j] * fsy[i][j]) / factor;
            double B = ((beta2 - beta1) * fsx[i][j] * fsy[i][j]) / factor;
            double C = B;
            double D = (beta2 * fsx[i][j] * fsx[i][j] + beta1 * fsy[i][j] * fsy[i][j]) / factor;
            inputVector[i] = ((A - D)*(A - D) + 4 * B*C) / ((A + D + m_epsilon)*(A + D + m_epsilon));
          }
        }

        cbica::Statistics< double > statsCalculator;
        statsCalculator.SetInput(inputVector);

        this->m_features["Minimum"] = statsCalculator.GetMinimum();
        this->m_features["Maximum"] = statsCalculator.GetMaximum();
        this->m_features["Mean"] = statsCalculator.GetMean();
        this->m_features["Sum"] = statsCalculator.GetSum();
        this->m_features["Mode"] = statsCalculator.GetMode();
        this->m_features["Median"] = statsCalculator.GetMedian();
        this->m_features["Variance"] = statsCalculator.GetVariance();
        this->m_features["StandardDeviation"] = statsCalculator.GetStandardDeviation();
        this->m_features["Skewness"] = statsCalculator.GetSkewness();
        this->m_features["Kurtosis"] = statsCalculator.GetKurtosis();
      } // end 2D 
      else if (TImageType::ImageDimension == 3) ///TBD: 3D calculation is not defined, yet
      {

      } // end 3D 
      else
      {
        std::cerr << "Image Dimensions is not supported.\n";
        exit(EXIT_FAILURE);
      }

      this->m_algorithmDone = true;

    }
  }

private:

  /**
  \brief As far as I can tell, this is doing convolution

  /// TBD: DO NOT USE RAW INDECES, this gets messed up with world+image coordinate information is present

  This needs to be converted to an itk::Image based computation
  */
  std::vector< std::vector< double > > gD(std::vector< std::vector< double > > &f, int scale, int ox, int oy)
  {
    int J = std::ceil(3 * scale);
    std::vector< double > x, Gs;
    double sum = 0;
    for (int i = -J; i < J + 1; i++)
    {
      x.push_back(i);
      Gs.push_back(exp(-i * i / (2 * scale*scale)));
      sum += exp(-i * i / (2 * scale*scale));
    }
    for (int i = 0; i < 2 * J + 1; i++)
    {
      Gs[i] = Gs[i] / sum;
    }
    std::vector< double > Gsx = gDerivative(ox, x, Gs, scale);
    std::vector< double > Gsy = gDerivative(oy, x, Gs, scale);
    int N = f.size();
    int M = f[0].size();
    int K = (Gsx.size() - 1) / 2;
    int L = (Gsy.size() - 1) / 2;

    std::vector< int > iind;
    iind.resize(N + 2 * K);
    for (int i = 0; i < N + 2 * K; i++)
    {
      iind[i] = std::min(std::max(i + 1 - K, 0), N - 1);
    }
    
    std::vector< int > jind;
    jind.resize(M + 2 * L);
    for (int i = 0; i < M + 2 * L; i++)
    {
      jind[i] = std::min(std::max(i + 1 - L, 0), M - 1);
    }

    std::vector< std::vector< double > > fwb;
    fwb.resize(iind.size());
    for (uint i = 0; i < iind.size(); i++)
    {
      fwb[i].resize(jind.size());
      for (uint j = 0; j < jind.size(); j++)
      {
        fwb[i][j] = f[iind[i]][jind[j]];
      }
    }

    std::vector< std::vector< double > > filter;
    filter.resize(Gsx.size());
    for (uint i = 0; i < Gsx.size(); i++)
    {
      filter[i].resize(Gsy.size());
      for (uint j = 0; j < Gsy.size(); j++)
      {
        filter[i][j] = Gsx[i] * Gsy[j];
      }
    }

    std::vector< std::vector< double > > conv_final;
    conv_final.resize(iind.size() - Gsx.size() + 1);
    for (uint i = 0; i < iind.size() - Gsx.size() + 1; i++)
    {
      conv_final[i].resize(jind.size() - Gsy.size() + 1);
      for (uint j = 0; j < jind.size() - Gsy.size() + 1; j++)
      {
        double sum = 0;
        for (uint k = 0; k < Gsx.size(); k++)
        {
          for (uint l = 0; l < Gsy.size(); l++)
          {
            sum += filter[k][l] * fwb[i + k][j + l];
          }
        }
        conv_final[i][j] = sum;
      }
    }
    return conv_final;
  }


  std::vector< double > gDerivative(int a, std::vector< double > &x, std::vector< double > &Gs, int scale)
  {
    std::vector< double > r;
    r.resize(x.size());
    for (uint i = 0; i < x.size(); i++)
    {
      if (a == 0)
      {
        r[i] = Gs[i];
      }
      else
      {
        r[i] = -x[i] / (scale*scale)*Gs[i];
      }
    }
    return r;
  }

  float m_Eta = 10, 
    m_epsilon = 10;
  int m_radius = -1; //! radius around which features are to be extracted
  float m_radius_float = -1; //! radius around which features are to be extracted
};