/**
\file FractalBoxCount.h

This file calculates the Fractal Dimension features of an input image and mask.

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

#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"
#include "cbicaLogging.h"

#include "FeatureBase.h"

template< class TImageType = itk::Image< float, 3 > >
class FractalBoxCount : public FeatureBase < TImageType >
{
public:
  //! Default Constructor
  FractalBoxCount() {};

  //! Default Destructor
  ~FractalBoxCount() {};

  //! Set the radius for computation
  void SetRadius(int radius) { m_radius = radius; };

  //! Whether the current image is part of a lattice point or not
  void SetLatticePointStatus(bool flag) { m_latticePoint = flag; };

  //! Get the Box Count
  typename TImageType::PixelType GetBoxCount()
  {
    if (!this->m_algorithmDone)
    {
      Update();
    }
    return m_boxCount;
  }

  //! Get the Minkovski dimension
  typename TImageType::PixelType GetMinkovski()
  {
    if (!this->m_algorithmDone)
    {
      Update();
    }
    return m_minkovski;
  }

  //! Actual algorithm runner
  void Update()
  {
    if (!this->m_algorithmDone)
    {
      //std::vector< double > box_counting_dimension;
      //std::vector< double > minkovski_dimension;

      //int maxDiameter = m_radius * 2 + 1; // this is not being used anywhere

      auto imageSize = this->m_inputImage->GetBufferedRegion().GetSize();

      /// TBD: need clarification from Yifan about this; Michael says this width needs to be the window size
      int width = 0;
      for (size_t d = 0; d < TImageType::ImageDimension; d++)
      {
        if (static_cast< int >(imageSize[d]) > width)
        {
          width = imageSize[d];
        }
      }

      float p = std::log(width) / std::log(2);
      p = std::ceil(p);
      width = std::pow(2, p);

      /// TBD: need clarification from Yifan about this; why would we use a mirror pad filter at all?
      // set up the lower and upper regions for the mirror pad filter
      typename TImageType::SizeType lowerExtendedRegion, upperExtendedRegion;
      for (size_t d = 0; d < TImageType::ImageDimension; d++)
      {
        lowerExtendedRegion[d] = 0;
        upperExtendedRegion[d] = width - imageSize[d];
      }
      if (!m_latticePoint)
      {
        upperExtendedRegion = lowerExtendedRegion;
      }

      // initialize and update the pad filter
      auto padFilter = itk::MirrorPadImageFilter< TImageType, TImageType >::New();
      padFilter->SetInput(this->m_inputImage);
      padFilter->SetPadLowerBound(lowerExtendedRegion);
      padFilter->SetPadUpperBound(upperExtendedRegion);
      padFilter->Update();

      auto paddedImage = padFilter->GetOutput();
      auto newSize = paddedImage->GetBufferedRegion().GetSize();

      typename FeatureBase< TImageType >::TConstIteratorType paddedImageIterator(paddedImage, paddedImage->GetBufferedRegion());

      // start the computation

      // initialize variables that are to be populated
      std::vector< double > area_sum;
      std::vector< double > minkovski_diff;
      std::vector< double > border_length;
      std::vector< double > border_length_minkov;

      /// this is only valid for 2D images... need to generalize for 3D
      if (TImageType::ImageDimension == 2)
      {
        std::vector< std::vector< double > > image_matrix;
        for (size_t j = 0; j < newSize[1]; j++)
        {
          std::vector< double > temp;
          for (size_t i = 0; i < newSize[0]; i++)
          {
            typename TImageType::IndexType index;
            index[0] = i;
            index[1] = j;
            temp.push_back(paddedImage->GetPixel(index));
          }
          image_matrix.push_back(temp);
        }

        // Definition of box-counting and Minkowski fractal dimension is in this paper 
        // "Fractal Analysis of Mammographic Parenchymal Patterns in Breast Cancer Risk Assessment", 
        // equation 2 for box-counting; equation 3 and 4 for minkowski

        for (float g = 0; g < p; g++)
        {
          // sum is for box-counting; sum1 is for minkowski
          double scale = std::pow(2, g); 
          double sum = 0;
          double sum1 = 0;
          double square_sum1 = 0, square_sum2 = 0;
          for (int i = 0; i < width - scale; i = i + scale)
          {
            for (int j = 0; j < width - scale; j = j + scale)
            {
              // box counting
              for (int k = 0; k < scale; k++)
              {
                for (int l = 0; l < scale; l++)
                {
                  // square_sum1 and square_sum2 are part of box-counting equation 2 in the reference
                  square_sum1 += std::abs(image_matrix[i + k][j + l] - image_matrix[i + k + scale][j + l]);
                  square_sum2 += std::abs(image_matrix[i + k][j + l] - image_matrix[i + k][j + l + scale]);
                }
              }
            }
          }

          // Minkovski
          if (scale > 1)
          {
            for (int i = 0; i < width; i++)
            {
              for (int j = 0; j < width; j++)
              {

                double temp_min = image_matrix[i][j], temp_max = image_matrix[i][j];
                for (int k = std::max(0, i - static_cast<int>(scale) / 2); k < std::min(i + static_cast<int>(scale) / 2, width); k++)
                {
                  for (int l = std::max(0, j - static_cast<int>(scale) / 2); l < std::min(j + static_cast<int>(scale) / 2, width); l++)
                  {
                    // temp_min and temp_max is from equation 4 in the reference
                    temp_min = temp_min < image_matrix[k][l] ? temp_min : image_matrix[k][l];
                    temp_max = temp_max > image_matrix[k][l] ? temp_max : image_matrix[k][l];
                  }
                }
                // equation 4
                sum1 += temp_max - temp_min;
              }
            }
          }
          // equation 2
          sum = sum + scale * scale + (std::abs(square_sum1) + std::abs(square_sum2)) * scale;

          // box counting
          area_sum.push_back(sum);
          border_length.push_back(scale);

          // Minkovski
          if (scale > 1)
          {
            // equation 3
            minkovski_diff.push_back(sum1 / (scale*scale*scale));
            border_length_minkov.push_back(1.0 / scale);
          }
        }

        double beta0 = 0, beta1 = 0;
        // box_counting
        double sum_y = 0;
        double sum_x = 0;
        double sum_y_x = 0;
        double sum_x_x = 0;

        // fit the line and get the slope which is the fractal dimension
        for (size_t i = 0; i < border_length.size(); i++)
        {
          sum_y_x += log_with_zero(area_sum[i]) * log_with_zero(border_length[i]);
          sum_y += log_with_zero(area_sum[i]);
          sum_x_x += log_with_zero(border_length[i])*log_with_zero(border_length[i]);
          sum_x += log_with_zero(border_length[i]);
        }

        //if (!border_length.empty())
        {
          beta0 = (sum_y_x - sum_y * sum_x / border_length.size()) / (sum_x_x - sum_x * sum_x / border_length.size());
        }

        //box_counting_dimension.push_back(2 - beta0);
        //cout << 2-beta0 << endl;
        //second time for minkowski, variable with same meaning
        sum_y = 0;
        sum_x = 0;
        sum_y_x = 0;
        sum_x_x = 0;

        for (size_t i = 0; i < border_length_minkov.size(); i++)
        {
          sum_y_x += log_with_zero(minkovski_diff[i]) * log_with_zero(border_length_minkov[i]);
          sum_y += log_with_zero(minkovski_diff[i]);
          sum_x_x += log_with_zero(border_length_minkov[i])*log_with_zero(border_length_minkov[i]);
          sum_x += log_with_zero(border_length_minkov[i]);
        }

        //if (!border_length_minkov.empty())
        {
          beta1 = (sum_y_x - sum_y * sum_x / border_length_minkov.size()) / (sum_x_x - sum_x * sum_x / border_length_minkov.size());
        }

        m_boxCount = 2 - beta0;
        m_minkovski = beta1;
        this->m_features["BoxCount"] = m_boxCount; /// TBD: why "2 - beta0"? I was not able to find the reference for this in the paper
        this->m_features["Minkovski"] = m_minkovski;
        //minkovski_dimension.push_back(beta1);
      } // end 2D 
      else if (TImageType::ImageDimension == 3) ///TBD: 3D calculation is not defined, yet
      {
        //std::vector< std::vector< std::vector< double > > > image_matrix;
        //for (int k = 0; k < newSize[2]; k++)
        //{
        //  std::vector< std::vector< double > > temp2;
        //  for (int j = 0; j < newSize[1]; j++)
        //  {
        //    std::vector< double > temp1;
        //    for (int i = 0; i < newSize[0]; i++)
        //    {
        //      typename TImageType::IndexType index;
        //      index[0] = i;
        //      index[1] = j;
        //      index[2] = k;
        //      temp1.push_back(paddedImage->GetPixel(index));
        //    }
        //    temp2.push_back(temp1);
        //  }
        //  image_matrix.push_back(temp2);
        //}

        //// Definition of box-counting and Minkowski fractal dimension is in this paper 
        //// "Fractal Analysis of Mammographic Parenchymal Patterns in Breast Cancer Risk Assessment", 
        //// equation 2 for box-counting; equation 3 and 4 for minkowski

        //for (int g = 0; g < p; g++)
        //{
        //  // sum is for box-counting; sum1 is for minkowski
        //  /// TBD: Michael: what is this scale?
        //  int scale = std::pow(2, g);
        //  float sum = 0;
        //  float sum1 = 0;
        //  float square_sum1 = 0, square_sum2 = 0, square_sum3 = 0;
        //  for (int i = 0; i < width - scale; i = i + scale)
        //  {
        //    for (int j = 0; j < width - scale; j = j + scale)
        //    {
        //      for (int m = 0; m < width - scale; m = m + scale)
        //      {
        //        // box counting
        //        for (int k = 0; k < scale; k++)
        //        {
        //          for (int l = 0; l < scale; l++)
        //          {
        //            for (int n = 0; n < scale; n++)
        //            {
        //              // square_sum1 and square_sum2 are part of box-counting equation 2 in the reference
        //              square_sum1 += abs(image_matrix[i + k][j + l][m + n] - image_matrix[i + k + scale][j + l][m + n]);
        //              square_sum2 += abs(image_matrix[i + k][j + l][m + n] - image_matrix[i + k][j + l + scale][m + n]);
        //              square_sum3 += abs(image_matrix[i + k][j + l][m + n] - image_matrix[i + k][j + l][m + n + scale]);
        //            }
        //          }
        //        }
        //      }
        //    }
        //  }

        //  // Minkovski
        //  if (scale > 1)
        //  {
        //    for (int i = 0; i < width; i++)
        //    {
        //      for (int j = 0; j < width; j++)
        //      {
        //        for (int m = 0; m < width; m++)
        //        {
        //          float temp_min = image_matrix[i][j][m]; float temp_max = image_matrix[i][j][m];
        //          for (int k = std::max(0, i - scale / 2); k < std::min(i + scale / 2, width); k++)
        //          {
        //            for (int l = std::max(0, j - scale / 2); l < std::min(j + scale / 2, width); l++)
        //            {
        //              for (int n = std::max(0, m - scale / 2); n < std::min(m + scale / 2, width); n++)
        //              {
        //                // temp_min and temp_max is from equation 4 in the reference
        //                temp_min = temp_min < image_matrix[k][l][n] ? temp_min : image_matrix[k][l][n];
        //                temp_max = temp_max > image_matrix[k][l][n] ? temp_max : image_matrix[k][l][n];
        //              }
        //            }
        //          }
        //          // equation 4
        //          sum1 += temp_max - temp_min;
        //        }
        //      }
        //    }
        //  }
        //  // equation 2
        //  sum = sum + scale * scale + (std::abs(square_sum1) + std::abs(square_sum2) + std::abs(square_sum3)) * scale;

        //  // box counting
        //  area_sum.push_back(sum);
        //  border_length.push_back(scale);

        //  // Minkovski
        //  if (scale > 1)
        //  {
        //    // equation 3
        //    minkovski_diff.push_back(sum1 / (scale*scale*scale));
        //    border_length_minkov.push_back(1.0 / scale);
        //  }
        //}

        //// box_counting
        //double sum_y = 0;
        //double sum_x = 0;
        //double sum_y_x = 0;
        //double sum_x_x = 0;

        //// fit the line and get the slope which is the fractal dimension
        //for (int i = 0; i < border_length.size(); i++)
        //{
        //  sum_y_x += log_with_zero(area_sum[i]) * log_with_zero(border_length[i]);
        //  sum_y += log_with_zero(area_sum[i]);
        //  sum_x_x += log_with_zero(border_length[i])*log_with_zero(border_length[i]);
        //  sum_x += log_with_zero(border_length[i]);
        //}

        //double beta0;
        //beta0 = (sum_y_x - sum_y * sum_x / border_length.size()) / (sum_x_x - sum_x * sum_x / border_length.size());


        ////box_counting_dimension.push_back(2 - beta0);
        ////cout << 2-beta0 << endl;
        ////second time for minkowski, variable with same meaning
        //sum_y = 0;
        //sum_x = 0;
        //sum_y_x = 0;
        //sum_x_x = 0;

        //for (int i = 0; i < border_length_minkov.size(); i++)
        //{
        //  sum_y_x += log_with_zero(minkovski_diff[i]) * log_with_zero(border_length_minkov[i]);
        //  sum_y += log_with_zero(minkovski_diff[i]);
        //  sum_x_x += log_with_zero(border_length_minkov[i])*log_with_zero(border_length_minkov[i]);
        //  sum_x += log_with_zero(border_length_minkov[i]);
        //}

        //double beta1;
        //beta1 = (sum_y_x - sum_y * sum_x / border_length_minkov.size()) / (sum_x_x - sum_x * sum_x / border_length_minkov.size());

        //m_boxCount = 2 - beta0;
        //m_minkovski = beta1;
        //m_features["BoxCount"] = m_boxCount; /// TBD: why "2 - beta0"? I was not able to find the reference for this in the paper
        //m_features["Minkovski"] = m_minkovski;
        //minkovski_dimension.push_back(beta1);
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

  int m_radius; //! radius around which features are to be extracted
  bool m_latticePoint = false; //! whether the current image is part of a lattice patch or not
  double m_boxCount; //! the first output
  double m_minkovski; //! the first output
};