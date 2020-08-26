/**
\file PowerSpectrum.h

This file calculates the Laws Measures Features of an input image and mask.

https://www.med.upenn.edu/cbica/captk/ <br>
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
#include "itkComplexToModulusImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkFFTShiftImageFilter.h"
#include "itkForwardFFTImageFilter.h"
#include "itkWrapPadImageFilter.h"
#include "itkCastImageFilter.h"

#include "itkOpenCVImageBridge.h"

#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"
#include "cbicaLogging.h"
#include "cbicaStatistics.h"

#include "FeatureBase.h"
//using TImageType = itk::Image< float, 3 >;

template< class TImageType = itk::Image< float, 3 > >
class PowerSpectrum : public FeatureBase < TImageType >
{
public:
  //! Default Constructor
  PowerSpectrum() {};

  //! Default Destructor
  ~PowerSpectrum() {};

  void SetCenter(const std::string &center) 
  {
    m_centerLocation = center;
    m_centerLocation = cbica::replaceString(m_centerLocation, "(", "");
    m_centerLocation = cbica::replaceString(m_centerLocation, ")", "");
    m_centerLocation = cbica::replaceString(m_centerLocation, "|", "x");
  }

  //! Actual algorithm runner
  void Update()
  {
    if (!this->m_algorithmDone)
    {
      auto imageSize = this->m_inputImage->GetLargestPossibleRegion().GetSize();
      
      using ITKOpenCVBridgeType = itk::OpenCVImageBridge;

      /// this is only valid for 2D images... need to generalize for 3D
      if (TImageType::ImageDimension == 2)
      {
        auto minSize = imageSize[0];
        for (size_t i = 1; i < TImageType::ImageDimension; i++)
        {
          if (minSize > imageSize[i])
          {
            minSize = imageSize[i];
          }
        }

        /// add image watch visualization for m_inputImage

        using TRealImageType = itk::Image< float, TImageType::ImageDimension >;
        auto caster = itk::CastImageFilter< TImageType, TRealImageType >::New();
        caster->SetInput(this->m_inputImage);
        caster->Update();
        auto padFilter = itk::WrapPadImageFilter< TRealImageType, TRealImageType >::New();
        padFilter->SetInput(caster->GetOutput());
        auto padding = imageSize;
        // Input size is [48, 62, 42].  Pad to [48, 64, 48].

        std::vector< int > reference;
        //create an array with numbers are a multiplication by 2,3 and 5 because the itk filter can only computed on a window size with multiplication by 2,3 and 5; it depends on our lattice size;
        for (int i = 0; i < 10; i++)
        {
          for (int j = 0; j < 7; j++)
          {
            for (int k = 0; k < 5; k++)
            {
              int s = std::pow(2, i) * std::pow(3, j) * std::pow(5, k);
              if (s < 1000)
              {
                reference.push_back(s);
              }
            }
          }
        }
        std::sort(reference.begin(), reference.end());
        //padding the window into the size in the array
        int kk;
        for (uint i = 0; i < reference.size(); i++)
        {
          kk = reference[i] - /*imageSize[0]*/minSize; /// TBD: why only the first dimension of the size and not any of the others or maybe use maxSize instead?
          if (kk >= 0)
            break;
        }
        padding.Fill(kk);

        padFilter->SetPadUpperBound(padding);
        padFilter->Update();

        //auto cvImage_padded = ITKOpenCVBridgeType::ITKImageToCVMat< TRealImageType >(padFilter->GetOutput());

        //use itk filter to calculate power spectrum: https://itk.org/ITKExamples/src/Filtering/FFT/ComputeImageSpectralDensity/Documentation.html
        using ForwardFFTFilterType = itk::ForwardFFTImageFilter< TRealImageType >;
        using ComplexImageType = typename ForwardFFTFilterType::OutputImageType;
        auto forwardFFTFilter = ForwardFFTFilterType::New();
        forwardFFTFilter->SetInput(padFilter->GetOutput());
        forwardFFTFilter->Update();

        auto complexToModulusFilter = itk::ComplexToModulusImageFilter< ComplexImageType, TRealImageType >::New();
        complexToModulusFilter->SetInput(forwardFFTFilter->GetOutput());
        complexToModulusFilter->Update();

        //auto cvImage_modulus = ITKOpenCVBridgeType::ITKImageToCVMat< TRealImageType >(complexToModulusFilter->GetOutput());

        using OutputImageType = itk::Image< unsigned short, TRealImageType::ImageDimension >;

        auto windowingFilter = itk::IntensityWindowingImageFilter< TRealImageType, OutputImageType >::New();
        windowingFilter->SetInput(complexToModulusFilter->GetOutput());
        windowingFilter->SetWindowMinimum(0);
        windowingFilter->SetWindowMaximum(100000);
        windowingFilter->Update();

        //auto cvImage_windowed = ITKOpenCVBridgeType::ITKImageToCVMat< OutputImageType >(windowingFilter->GetOutput());

        auto fftShiftFilter = itk::FFTShiftImageFilter< OutputImageType, OutputImageType >::New();
        fftShiftFilter->SetInput(windowingFilter->GetOutput());
        fftShiftFilter->Update();

        auto fftShiftFilterOutput = fftShiftFilter->GetOutput();

        auto cvImage_fftShift = ITKOpenCVBridgeType::ITKImageToCVMat< OutputImageType >(fftShiftFilterOutput);

        //cbica::WriteImage< OutputImageType >(fftShiftFilterOutput, "C:/Projects/CaPTk-dev/branches/FeatureExtraction/data/fromMichael/fe_output/fftShift.nii.gz");

        itk::ImageRegionConstIteratorWithIndex< OutputImageType > fftShiftFilterOutputIterator(fftShiftFilterOutput, fftShiftFilterOutput->GetBufferedRegion());

        auto size = fftShiftFilterOutput->GetBufferedRegion().GetSize();

        int power_number = minSize / 2;
        std::vector< typename OutputImageType::PixelType > count(power_number + 1);
        std::vector< typename OutputImageType::PixelType > sum_power_density(power_number + 1);

        auto refer_matrix = index_matrix(size[0], size[1]);

        /// this code does not work because the currentIndex starts at 143x143
        //for (fftShiftFilterOutputIterator.GoToBegin(); !fftShiftFilterOutputIterator.IsAtEnd(); ++fftShiftFilterOutputIterator)
        //{
        //  auto currentIndex = fftShiftFilterOutputIterator.GetIndex();
        //  auto temp = fftShiftFilterOutputIterator.Get();

        //  auto indexToAccumulate = std::min(power_number, refer_matrix[currentIndex[0]][currentIndex[1]]);

        //  count[indexToAccumulate]++;
        //  sum_power_density[indexToAccumulate] += fftShiftFilterOutputIterator.Get();
        //}

        //compute average power correspond to different radius; paper: Modelling the Power Spectra of Natural Images: Statistics and Information
        for (uint i = 0; i < size[0]; i++)
        {
          /// TBD: DO NOT USE RAW INDECES, this gets messed up with world+image coordinate information is present
          auto Mi = cvImage_fftShift.template ptr< typename OutputImageType::PixelType >(i);
          for (uint j = 0; j < size[1]; j++)
          {
            //typename TRealImageType::IndexType index;
            //index[0] = i;
            //index[1] = j;
            //cout << fftShiftFilter->GetOutput()->GetPixel(index) << endl;
            //fftShiftFilterOutputIterator.SetIndex(index);
            //if (fftShiftFilterOutputIterator.Get() == 0)
            //{
            //  std::cerr << "Zero detected at [" << i << "," << j << "]\n";
            //}

            //auto temp1 = fftShiftFilterOutputIterator.Get();
            //double temp2 = fftShiftFilterOutputIterator.Get();
            //int temp3 = fftShiftFilterOutputIterator.Get();
            //unsigned short temp4 = fftShiftFilterOutputIterator.Get();
            //auto temp5 = fftShiftFilterOutput->GetPixel(index);
            //auto temp6 = Mi[j];
            
            auto indexToAccumulate = std::min(power_number, refer_matrix[i][j]);
            //std::cout << "[DEBUG] fftShiftFilterOutputIterator.Get(" << i << "," << j << ") = " << fftShiftFilterOutputIterator.Get() << "\n";
            //if (refer_matrix[i][j] > power_number)
            //{
            //  count[power_number]++;
            //  sum_power_density[power_number] += fftShiftFilterOutput->GetPixel(index)/* fftShiftFilterOutputIterator.Get()*/;
            //}
            //else
            //{
              count[indexToAccumulate]++;
              sum_power_density[indexToAccumulate] += /*fftShiftFilterOutput->GetPixel(index)*//*fftShiftFilterOutputIterator.Get()*/Mi[j];
            //}
          }
        }
        //compute the average on different radius
        double sum_y = 0;
        double sum_x = 0;
        double sum_y_x = 0;
        double sum_x_x = 0;
        double n = static_cast<double>(power_number);
        for (int i = 1; i < power_number + 1; i++)
        {
          //double power_y = (sum_power_density[i] + 0.0) / (count[i] + 0.0);
          //aver_power.push_back((sum_power_density[i] + 0.0) / (count[i] + 0.0));
          //dist.push_back(i);
          auto log_i = std::log(i);
          auto log_powerY = std::log((static_cast<double>(sum_power_density[i])) / (static_cast<double>(count[i])));
          sum_y_x += log_powerY * log_i;
          sum_y += log_powerY;
          sum_x_x += std::pow(log_i, 2);
          sum_x += log_i;
        }
        //fit the line and get the slope
        double beta = (n * sum_y_x - sum_y * sum_x) / (n * sum_x_x - sum_x * sum_x);
        this->m_features["Beta"] = -beta;

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

  std::string m_centerLocation;

  std::vector< std::vector< int > > index_matrix(int n1, int n2)
  {
    std::vector< std::vector< int > > index(n1);
    std::vector< int > center = { static_cast<int>(std::floor(n1 / 2)) + 1, static_cast<int>(std::floor(n2 / 2)) + 1 };
    for (int i = 0; i < n1; i = i + 1)
    {
      index[i].resize(n2);
      for (int j = 0; j < n2; j = j + 1)
      {
        index[i][j] = std::max(std::abs(i - center[0]), std::abs(j - center[1]));
      }
    }
    //std::vector< std::vector< int > > index;
    //int center[2] = { static_cast<int>(std::floor(n1 / 2) + 1), static_cast<int>(std::floor(n2 / 2) + 1) };
    //for (int i = 0; i < n1; i = i + 1)
    //{
    //  std::vector< int > temp;
    //  for (int j = 0; j < n2; j = j + 1)
    //  {
    //    temp.push_back(std::max(std::abs(i - center[0]), std::abs(j - center[1])));
    //  }
    //  index.push_back(temp);
    //}
    return index;
  }

  //typename TImageType::Pointer index_matrix(typename TImageType::SizeType size)
  //{
  //  auto returnImage = typename TImageType::New();
  //  typename TImageType::RegionType region_new;
  //  typename TImageType::IndexType start;
  //  start.Fill(0);
  //  region_new.SetSize(size);
  //  returnImage->SetRegions(region_new);
  //  returnImage->Allocate();

  //  typename TImageType::IndexType center;
  //  for (size_t d = 0; d < TImageType::ImageDimension; d++)
  //  {
  //    center[d] = std::floor(size[d] / 2) + 1;
  //  }
  //  itk::ImageIteratorWithIndex< TImageType > iterator(returnImage, returnImage->GetLargestPossibleRegion());
  //  for (iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator)
  //  {
  //    auto temp = iterator.Get();
  //    auto currentIndex = iterator.GetIndex();
  //    for (size_t d = 0; d < TImageType::ImageDimension; d++)
  //    {
  //      temp = std::max(std::abs(currentIndex[d] - center[d]), temp);
  //    }
  //  }
  //  return returnImage;
  //}

};