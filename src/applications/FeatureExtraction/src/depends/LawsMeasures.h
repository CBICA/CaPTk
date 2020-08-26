/**
\file LawsMeasures.h

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

#include "FeatureBase.h"

template< class TImageType = itk::Image< float, 3 > >
class LawsMeasures : public FeatureBase < TImageType >
{
  using TConstIteratorType = itk::ImageRegionConstIterator< TImageType >;
public:
  //! Default Constructor
  LawsMeasures() {};

  //! Default Destructor
  ~LawsMeasures() {};

  //! Actual algorithm runner
  void Update()
  {
    if (!this->m_algorithmDone)
    {
      auto imageSize = this->m_inputImage->GetBufferedRegion().GetSize();

      /// this is only valid for 2D images... need to generalize for 3D
      if (TImageType::ImageDimension == 2)
      {
        auto cvImage_input = itk::OpenCVImageBridge::ITKImageToCVMat< TImageType >(this->m_inputImage);

        std::vector< std::vector< double > > image_matrix;
        image_matrix.resize(imageSize[1]);
        for (size_t j = 0; j < imageSize[1]; j++)
        {
          image_matrix[j].resize(imageSize[0]);
          for (size_t i = 0; i < imageSize[0]; i++)
          {
            //auto Mi = cvImage_input.ptr< typename TImageType::PixelType >(i);
            //typename TImageType::IndexType index;
            //index[0] = i;
            //index[1] = j;
            //temp2.push_back(cvImage_input.at< typename TImageType::PixelType >(i, j));
            image_matrix[j][i] = cvImage_input.template ptr< typename TImageType::PixelType >(i)[j];
          }
          //image_matrix.push_back(temp2);
        }

        typename TImageType::RegionType region_new;
        typename TImageType::IndexType start;
        start.Fill(0);

        typename TImageType::SizeType size_new_matrix;
        size_new_matrix[0] = imageSize[0] - 5 + 1;
        size_new_matrix[1] = imageSize[1] - 5 + 1;
        region_new.SetSize(size_new_matrix);

        // here we can also choose mask size=3 by using laws_masks_3 function
        auto mask_matrix = laws_masks();
        int mask_num = mask_matrix.size();
        typename TImageType::IndexType pixelIndex;

        //25 masks(9 masks when using mask size=3)
        //feature for each mask
        for (int m = 0; m < mask_num; m++)
        {
          auto result_matrix = TImageType::New();
          result_matrix->SetRegions(region_new);
          result_matrix->Allocate();
          result_matrix->FillBuffer(itk::NumericTraits< double >::Zero);

          //convolution step
          for (size_t i = 0; i < size_new_matrix[1]; i++) ///TBD: why does the iteration start from [1] but the variable refers to the rows?
          {
            for (size_t j = 0; j < size_new_matrix[0]; j++)
            {
              double sum = 0;
              for (int k = 0; k < 5; k++)
              {
                for (int l = 0; l < 5; l++)
                {
                  sum += mask_matrix[m][k][l] * image_matrix[i + k][j + l];
                }
              }
              pixelIndex[0] = j;
              pixelIndex[1] = i;
              if (sum != 0)
              {
                result_matrix->SetPixel(pixelIndex, sum);
              }
            }
          }

          auto statisticsImageFilter = itk::StatisticsImageFilter< TImageType >::New();
          statisticsImageFilter->SetInput(result_matrix);
          statisticsImageFilter->Update();
          double origin_mean = statisticsImageFilter->GetMean();
          double origin_sigma = statisticsImageFilter->GetSigma();
          //mean and standard deviation
          this->m_features["Mean_" + std::to_string(m + 1)] = origin_mean;
          this->m_features["STD_" + std::to_string(m + 1)] = origin_sigma;

          auto addImageFilter = itk::AddImageFilter< TImageType, TImageType, TImageType >::New();
          addImageFilter->SetInput(result_matrix);
          addImageFilter->SetConstant2(-origin_mean);
          addImageFilter->Update();
          auto result_matrix1 = addImageFilter->GetOutput();

          //power of value: 3 and 4
          using PowerFilter = itk::PowImageFilter< TImageType, TImageType, TImageType >;
          auto powImageFilter1 = PowerFilter::New();
          auto powImageFilter2 = PowerFilter::New();
          //skewness
          powImageFilter1->SetInput(result_matrix1);
          powImageFilter1->SetConstant2(3);
          powImageFilter1->Update();
          auto skewness_matrix = powImageFilter1->GetOutput();
          statisticsImageFilter->SetInput(skewness_matrix);
          statisticsImageFilter->Update();
          if (origin_sigma == 0)
            this->m_features["Skewness_" + std::to_string(m + 1)] = 0;
          else
            this->m_features["Skewness_" + std::to_string(m + 1)] = (statisticsImageFilter->GetMean()) / origin_sigma;
          //kurtosis
          powImageFilter2->SetInput(result_matrix1);
          powImageFilter2->SetConstant2(4);
          powImageFilter2->Update();
          auto kurtosis_matrix = powImageFilter2->GetOutput();
          statisticsImageFilter->SetInput(kurtosis_matrix);
          statisticsImageFilter->Update();
          if (origin_sigma == 0)
            this->m_features["Kurtosis_" + std::to_string(m + 1)] = 0;
          else
            this->m_features["Kurtosis_" + std::to_string(m + 1)] = (statisticsImageFilter->GetMean()) / pow(origin_sigma, 4) - 3;
          //entropy
          auto powImageFilter3 = PowerFilter::New();
          powImageFilter3->SetInput(result_matrix);
          powImageFilter3->SetConstant2(2);
          powImageFilter3->Update();
          statisticsImageFilter->SetInput(powImageFilter3->GetOutput());
          statisticsImageFilter->Update();
          this->m_features["Entropy_" + std::to_string(m + 1)] = statisticsImageFilter->GetMean();
        }

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
  \brief Get masks
  */
  std::vector< std::vector< std::vector <int> > > laws_masks()
  {
    int all_vectors[5][5] = { { 1, 4, 6, 4, 1 },        //Level detection
    { -1, -2, 0, 2, 1 },                                //Edge detection
    { -1, 0, 2, 0, -1 },                                //Spot detection
    { -1, 2, 0, -2, 1 },                                //Wave detection
    { 1, -4, 6, -4, 1 } };                              //Ripple detection

                                                        //mask by level and edge could be considered as gradient
    std::vector< std::vector< std::vector <int> > > all_masks;

    for (int i = 0; i < 5; i++)
    {
      for (int j = 0; j < 5; j++)
      {
        std::vector< std::vector <int> > temp;
        for (int row = 0; row < 5; row++)
        {
          std::vector< int > temp_vector;
          for (int col = 0; col < 5; col++)
          {
            temp_vector.push_back(all_vectors[i][row] * all_vectors[j][col]);
          }
          temp.push_back(temp_vector);
        }
        all_masks.push_back(temp);
      }
    }
    return all_masks;
  }

  /**
  \brief Get masks
  */
  std::vector< std::vector< std::vector <int> > > laws_masks_3()
  {
    int all_vectors[3][3] = { { 1, 2, 1 },        //Level detection
    { -1, 0, 1 },                                //Edge detection
    { -1, 2, -1 }, };                            //Spot detection

                                                 //mask by level and edge could be considered as gradient
    std::vector< std::vector< std::vector <int> > > all_masks;

    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        std::vector< std::vector <int> > temp;
        for (int row = 0; row < 3; row++)
        {
          std::vector< int > temp_vector;
          for (int col = 0; col < 3; col++)
          {
            temp_vector.push_back(all_vectors[i][row] * all_vectors[j][col]);
          }
          temp.push_back(temp_vector);
          temp_vector.clear();
        }
        all_masks.push_back(temp);
      }
    }
    return all_masks;
  }

};