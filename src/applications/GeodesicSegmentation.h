/**
\file  GeodesicSegmentation.h

\brief The header file containing the Geodesic segmentation class, used to apply an adaptive geodesic transform

Library Dependecies: ITK 4.7+ <br>
Header Dependencies: cbicaUtilities.h, cbicaLogging.h

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html

*/
#pragma once
//
//#include <iostream>
//#include <limits.h>

//#include "itkImage.h"
//#include "itkConnectedThresholdImageFilter.h"
//#include "itkImageRegionIterator.h"
//#include "itkBSplineControlPointImageFilter.h"
//#include "itkExpImageFilter.h"
//#include "itkImageRegionIterator.h"
//#include "itkOtsuThresholdImageFilter.h"
//#include "itkShrinkImageFilter.h"
//#include "itkMedianImageFunction.h"
//#include "itkNeighborhoodIterator.h"
//#include "itkMinimumMaximumImageCalculator.h"
//#include "itkConnectedComponentImageFilter.h"
//#include "itkBinaryThresholdImageFilter.h"
//#include "itkThresholdImageFilter.h"

//#include "fProgressDialog.h"

//#include "PreprocessingPipelineClass.h"
#ifdef APP_BASE_CAPTK_H
#include "ApplicationBase.h"
#endif
//#include "CAPTk.h"
#include "CaPTkDefines.h"
#include "cbicaLogging.h"
#include "itkNeighborhoodIterator.h"
//#include <vector>
//#include "cbicaLogging.h"
//using VectorVectorDouble = std::vector< std::vector < double > >;

const int MAX_VAL = /*static_cast<int>(_I32_MAX)*/100000; // maximum possible Geodesic distance -- this value doesn't work for float type images


/**
\class GeodesicSegmentation

\brief Applies an adaptive Geodesic filter to image

Reference:

@inproceedings{gaonkar2014adaptive,
title={Adaptive geodesic transform for segmentation of vertebrae on CT images},
author={Gaonkar, Bilwaj and Shu, Liao and Hermosillo, Gerardo and Zhan, Yiqiang},
booktitle={SPIE Medical Imaging},
pages={903516--903516},
year={2014},
organization={International Society for Optics and Photonics}
}
*/
template< class ImageTypeGeodesic = ImageTypeShort3D >
class GeodesicSegmentation 
#ifdef APP_BASE_CAPTK_H
  : public ApplicationBase
#endif
{
public:
  explicit GeodesicSegmentation()
  {
    m_maxVal = /*std::numeric_limits< typename ImageTypeGeodesic::PixelType >::max()*/100000;
  }

  virtual ~GeodesicSegmentation()
  {

  }

  void cleanUp()
  {
    
  }

  //template<class ImageTypeGeodesic = ImageTypeShort3D >
  typename ImageTypeGeodesic::Pointer Run(typename ImageTypeGeodesic::Pointer Inp, VectorVectorDouble &tumorPoints)
  {
    typename ImageTypeGeodesic::Pointer mask = ImageTypeGeodesic::New();
    mask->CopyInformation(Inp);
    mask->SetRequestedRegion(Inp->GetLargestPossibleRegion());
    mask->SetBufferedRegion(Inp->GetBufferedRegion());
    mask->Allocate();
    mask->FillBuffer(2);

    return Run/*< ImageTypeGeodesic >*/(Inp, mask, tumorPoints);
  }

  //template<class ImageTypeGeodesic = ImageTypeShort3D >
  typename ImageTypeGeodesic::Pointer Run(typename ImageTypeGeodesic::Pointer Inp, typename ImageTypeGeodesic::Pointer MaskImage, VectorVectorDouble &tumorPoints)
  {
    //--------------allocate a few images ------------------------
#ifdef APP_BASE_CAPTK_H
    messageUpdate("Geodesic Segmentation");
    progressUpdate(0);
#endif

    //typename ImageTypeGeodesic::Pointer Init = ImageTypeGeodesic::New();
    //Init->CopyInformation(Inp);
    //Init->SetRequestedRegion(Inp->GetLargestPossibleRegion());
    //Init->SetBufferedRegion(Inp->GetBufferedRegion());
    //Init->Allocate();
    //Init->FillBuffer(0);

    typename ImageTypeGeodesic::Pointer Geos = ImageTypeGeodesic::New();
    Geos->CopyInformation(Inp);
    Geos->SetRequestedRegion(Inp->GetLargestPossibleRegion());
    Geos->SetBufferedRegion(Inp->GetBufferedRegion());
    Geos->Allocate();
    Geos->FillBuffer(0);

    typename ImageTypeGeodesic::Pointer Gamma = ImageTypeGeodesic::New();
    Gamma->CopyInformation(Inp);
    Gamma->SetRequestedRegion(Inp->GetLargestPossibleRegion());
    Gamma->SetBufferedRegion(Inp->GetBufferedRegion());
    Gamma->Allocate();
    Gamma->FillBuffer(1);

    typename ImageTypeGeodesic::Pointer tumorMask = ImageTypeGeodesic::New();
    tumorMask->CopyInformation(Inp);
    tumorMask->SetRequestedRegion(Inp->GetLargestPossibleRegion());
    tumorMask->SetBufferedRegion(Inp->GetBufferedRegion());
    tumorMask->Allocate();
    tumorMask->FillBuffer(0);

#ifdef APP_BASE_CAPTK_H
    progressUpdate(10);
    qApp->processEvents();
#endif


    //---------------calculation of initial mask--------------------------
    typedef itk::ImageRegionIteratorWithIndex <ImageTypeGeodesic> IteratorType;
    for (unsigned int i = 0; i < tumorPoints.size(); i++)
    {
      // get index from the input points
      typename ImageTypeGeodesic::IndexType index;
      index[0] = tumorPoints[i][0];
      index[1] = tumorPoints[i][1];
      index[2] = tumorPoints[i][2];

      // initialize the mask and geodesic images 
      //Init->SetPixel(index, static_cast<typename ImageTypeGeodesic::PixelType>(255));
      Geos->SetPixel(index, static_cast<typename ImageTypeGeodesic::PixelType>(m_maxVal));
    }
#ifdef APP_BASE_CAPTK_H
    progressUpdate(20);
    qApp->processEvents();
#endif

    //---------------------------actual geodesic segmentation--------------------------
    typedef itk::ImageRegionIteratorWithIndex <ImageTypeGeodesic> IteratorType;
    IteratorType GeosIt(Geos, Geos->GetLargestPossibleRegion());
    IteratorType GamIt(Gamma, Gamma->GetLargestPossibleRegion());
    IteratorType MaskIt(MaskImage, MaskImage->GetLargestPossibleRegion());

#ifdef APP_BASE_CAPTK_H
    progressUpdate(22);
    qApp->processEvents();
#endif

    //Setting up the neighborhood iterator
    typename ImageTypeGeodesic::SizeType radius;
    radius[0] = 1;
    radius[1] = 1;
    radius[2] = 1;
    itk::NeighborhoodIterator<ImageTypeGeodesic> ResNIt(radius, Geos, Geos->GetLargestPossibleRegion());
    itk::NeighborhoodIterator<ImageTypeGeodesic> InpNIt(radius, Inp, Inp->GetLargestPossibleRegion());

    //The main loops
    cbica::Logging(loggerFile, "Main loops execution : Forward pass");

    MaskIt.GoToBegin();

#ifdef APP_BASE_CAPTK_H
    progressUpdate(25);
    qApp->processEvents();
#endif

    while (!MaskIt.IsAtEnd())
    {
      if (MaskIt.Get() > 1)
      {
        GamIt.SetIndex(MaskIt.GetIndex());
        GeosIt.SetIndex(MaskIt.GetIndex());

        double C_f_arr[14];

        // forward pass
        C_f_arr[13] = GeosIt.Get();
        C_f_arr[0] = ResNIt.GetPixel(4) + sqrt(1.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(4))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(4))));
        C_f_arr[1] = ResNIt.GetPixel(10) + sqrt(1.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(10))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(10))));
        C_f_arr[2] = ResNIt.GetPixel(12) + sqrt(1.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(12))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(12))));
        C_f_arr[3] = ResNIt.GetPixel(1) + sqrt(2.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(1))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(1))));
        C_f_arr[4] = ResNIt.GetPixel(3) + sqrt(2.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(3))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(3))));
        C_f_arr[5] = ResNIt.GetPixel(9) + sqrt(2.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(9))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(9))));
        C_f_arr[6] = ResNIt.GetPixel(0) + sqrt(3.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(0))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(0))));
        C_f_arr[7] = ResNIt.GetPixel(7) + sqrt(2.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(7))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(7))));
        C_f_arr[8] = ResNIt.GetPixel(6) + sqrt(3.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(6))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(6))));
        C_f_arr[9] = ResNIt.GetPixel(15) + sqrt(2.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(15))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(15))));
        C_f_arr[10] = ResNIt.GetPixel(24) + sqrt(3.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(24))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(24))));
        C_f_arr[11] = ResNIt.GetPixel(21) + sqrt(2.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(21))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(21))));
        C_f_arr[12] = ResNIt.GetPixel(18) + sqrt(3.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(18))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(18))));

        double minval = m_maxVal * 200;
        for (int i = 0; i < 14; ++i)
        {
          if (C_f_arr[i] < minval)
            minval = C_f_arr[i];
        }

        GeosIt.Set(minval);
      }
      ++MaskIt;
      ++ResNIt;
      ++InpNIt;
    }
#ifdef APP_BASE_CAPTK_H
    progressUpdate(55);
    qApp->processEvents();
#endif

    MaskIt.GoToReverseBegin();
    ResNIt.GoToEnd();
    InpNIt.GoToEnd();
    --ResNIt;
    --InpNIt;
    cbica::Logging(loggerFile, "Main loops execution : Backward pass");

    while (!MaskIt.IsAtReverseEnd())
    {
      if (MaskIt.Get() > 1)
      {
        GamIt.SetIndex(MaskIt.GetIndex());
        GeosIt.SetIndex(MaskIt.GetIndex());

        double C_b_arr[14];
        C_b_arr[13] = GeosIt.Get();
        C_b_arr[0] = ResNIt.GetPixel(22) + sqrt(1.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(22))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(22))));
        C_b_arr[1] = ResNIt.GetPixel(16) + sqrt(1.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(16))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(16))));
        C_b_arr[2] = ResNIt.GetPixel(14) + sqrt(1.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(14))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(14))));
        C_b_arr[3] = ResNIt.GetPixel(25) + sqrt(2.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(25))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(25))));
        C_b_arr[4] = ResNIt.GetPixel(23) + sqrt(2.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(23))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(23))));
        C_b_arr[5] = ResNIt.GetPixel(17) + sqrt(2.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(17))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(17))));
        C_b_arr[6] = ResNIt.GetPixel(26) + sqrt(3.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(26))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(26))));
        C_b_arr[7] = ResNIt.GetPixel(19) + sqrt(2.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(19))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(19))));
        C_b_arr[8] = ResNIt.GetPixel(20) + sqrt(3.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(20))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(20))));
        C_b_arr[9] = ResNIt.GetPixel(11) + sqrt(2.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(11))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(11))));
        C_b_arr[10] = ResNIt.GetPixel(2) + sqrt(3.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(2))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(2))));
        C_b_arr[11] = ResNIt.GetPixel(5) + sqrt(2.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(5))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(5))));
        C_b_arr[12] = ResNIt.GetPixel(8) + sqrt(3.0 + GamIt.Get()*((InpNIt.GetPixel(13) - InpNIt.GetPixel(8))*(InpNIt.GetPixel(13) - InpNIt.GetPixel(8))));
        double minval = m_maxVal * 200;
        for (int i = 0; i < 14; ++i)
        {
          if (C_b_arr[i] < minval)
            minval = C_b_arr[i];
        }

        if (minval >= m_maxVal)
        {
          GeosIt.Set(0);
        }
        else
        {
          GeosIt.Set(minval);
        }
      }
      --MaskIt;
      --ResNIt;
      --InpNIt;
    }

#ifdef APP_BASE_CAPTK_H
    progressUpdate(90);
    qApp->processEvents();
#endif
    //------------------------------geodesic thresholding-----------------------------
    //typedef itk::MinimumMaximumImageCalculator <ImageTypeGeodesic> ImageCalculatorFilterType;
    //typename ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
    //imageCalculatorFilter->SetImage(Geos);
    //imageCalculatorFilter->Compute();
    //typename ImageTypeGeodesic::PixelType minValue = imageCalculatorFilter->GetMinimum();
    //typename ImageTypeGeodesic::PixelType maxValue = imageCalculatorFilter->GetMaximum();
    //typename ImageTypeGeodesic::PixelType average = (minValue + maxValue) / 2;
    //typename ImageTypeGeodesic::PixelType average_sqrt = std::sqrt(std::abs(average));
    //typename ImageTypeGeodesic::PixelType lowerLimit, upperLimit;

    //typedef itk::RescaleIntensityImageFilter<ImageTypeGeodesic> RescalerType;
    //typename RescalerType::Pointer rescaler = RescalerType::New();
    //rescaler->SetOutputMaximum(255);
    //rescaler->SetOutputMinimum(0);
    //rescaler->SetInput(Geos);
    //rescaler->Update();

    //if (minValue < 0)
    //{
    //  lowerLimit = minValue + average_sqrt / 2;
    //  upperLimit = minValue - average_sqrt;
    //}
    //else
    //{
    //  lowerLimit = minValue - average_sqrt / 2;
    //  upperLimit = minValue + average_sqrt;
    //}
    //typename ImageTypeGeodesic::PixelType lowerLimit = minValue + average_sqrt / 2;
    //typename ImageTypeGeodesic::PixelType upperLimit = minValue - average_sqrt;

    //typedef itk::BinaryThresholdImageFilter< ImageTypeGeodesic, ImageTypeGeodesic > BinaryThresholderType;
    //typename BinaryThresholderType::Pointer binaryThresholdFilter = BinaryThresholderType::New();
    //binaryThresholdFilter->SetLowerThreshold(0);
    //binaryThresholdFilter->SetUpperThreshold(25);
    //binaryThresholdFilter->SetInput(rescaler->GetOutput());
    //binaryThresholdFilter->SetOutsideValue(static_cast<typename ImageTypeGeodesic::PixelType>(0));
    //binaryThresholdFilter->SetInsideValue(static_cast<typename ImageTypeGeodesic::PixelType>(1));
    //binaryThresholdFilter->Update();

    //typedef itk::ConnectedComponentImageFilter <ImageTypeGeodesic, ImageTypeGeodesic > ConnectedComponentImageFilterType;
    //typename ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
    //connected->SetInput(binaryThresholdFilter->GetOutput());
    //connected->Update();

#ifdef APP_BASE_CAPTK_H
    progressUpdate(100);
    qApp->processEvents();
#endif

    cleanUp();
    return Geos;
  }


private:
  inline void SetLongRunning(bool longRunning = false)
  {

  }

  double m_maxVal;

  //typename ImageTypeGeodesic::Pointer Init, Geos, Gamma, tumorMask

};
