/**
\file SusanDenoising.h

This file holds the declaration of the class SusanDenoising.

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html
*/

#pragma once

//#include <QApplication>
#include "iostream"
#include "itkImage.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkMedianImageFunction.h"
#include "itkNeighborhoodIterator.h"
#include "itkImageDuplicator.h"
//#include "ApplicationBase.h"

/**
\class SusanDenoising

\brief This class implements the SUSAN denoising algorithm.

Reference:

@inproceedings{SmithBrady1997,
title={SUSAN - a new approach to low level image processing},
author={Smith, S.M., and Brady, J.M.},
journal={International Journal of Computer Vision},
pages={45�78},
year={1997},
organization={Springer}
}
*/
class SusanDenoising /*: public ApplicationBase*/
{

public:
  SusanDenoising() {};
  ~SusanDenoising() {};

  template<class ImageType>
  typename ImageType::Pointer Run(typename ImageType::Pointer image);

  void SetSigma(float input)
  {
    m_sigma = input;
  }

  void SetIntensityVariationThreshold(float input)
  {
    m_intensityVariationThreshold = input;
  }

  void SetRadius(size_t input)
  {
    m_radius = input;
  }


private:
  inline void SetLongRunning(bool longRunning);

  float m_sigma = 0.5;
  float m_intensityVariationThreshold = 80;
  size_t m_radius = 1;
};

template<class ImageType>
typename ImageType::Pointer SusanDenoising::Run(const typename ImageType::Pointer image)
{
	
  //messageUpdate("SUSAN Denoising");
  //progressUpdate(0);
  //qApp->processEvents();

  typedef itk::ImageDuplicator<ImageType> DuplicatorType;
  typename DuplicatorType::Pointer outputImageFilter = DuplicatorType::New();
  outputImageFilter->SetInputImage(image);
  outputImageFilter->Update();

  double sigma = m_sigma;
  double intensityVariationThreshold = m_intensityVariationThreshold;

  // typename ImageType::SizeType imageSize = image->GetLargestPossibleRegion().GetSize();

  typedef typename itk::MedianImageFunction< ImageType > MedianImageFunctionType;
  typename MedianImageFunctionType::Pointer medianImageFunction = MedianImageFunctionType::New();
  medianImageFunction->SetInputImage(outputImageFilter->GetOutput());


  typedef typename itk::NeighborhoodIterator<ImageType> NeighborhoodIteratorType;
  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill(m_radius);
  NeighborhoodIteratorType ImageIterator(radius, outputImageFilter->GetOutput(), outputImageFilter->GetOutput()->GetLargestPossibleRegion());

  //progressUpdate(10);
  //qApp->processEvents();

  for (ImageIterator.GoToBegin(); !ImageIterator.IsAtEnd(); ++ImageIterator)
  {
    double NumeratorSum = 0;
    double DenominatorSum = 0;

    typename itk::NeighborhoodIterator<ImageType>::IndexType currentindex = ImageIterator.GetIndex();
    float CenterIntensityValue = ImageIterator.GetCenterPixel();

    // make parallel here, perhaps?
    for (unsigned int LocalNeighborhoodIterator = 0; LocalNeighborhoodIterator < ImageIterator.Size(); ++LocalNeighborhoodIterator)
    {
      typename ImageType::OffsetType offsetType1 = ImageIterator.ComputeInternalIndex(LocalNeighborhoodIterator);

      float NeighborIntensityValue = ImageIterator.GetPixel(LocalNeighborhoodIterator);

      typename itk::NeighborhoodIterator<ImageType>::IndexType neighborindex;
      for (unsigned int index = 0; index < 3; index++)
        neighborindex[index] = (currentindex[index] - radius[index]) + offsetType1[index];

      typename ImageType::IndexType LocalizedVoxelIndex;

      for (unsigned int index = 0; index < LocalizedVoxelIndex.GetIndexDimension(); index++)
        LocalizedVoxelIndex[index] = neighborindex[index] - currentindex[index];

      if (LocalizedVoxelIndex == currentindex)
        continue;

      double radius = 0;
      for (unsigned int index = 0; index < LocalizedVoxelIndex.GetIndexDimension(); index++)
        radius = radius + pow(LocalizedVoxelIndex[index], 2);
      radius = sqrt(radius);

      double weightfactor = exp(-pow(radius, 2) / (2 * pow(sigma, 2)) - pow((NeighborIntensityValue - CenterIntensityValue), 2) / pow(intensityVariationThreshold, 2));
      DenominatorSum = DenominatorSum + weightfactor;
      NumeratorSum = NumeratorSum + (NeighborIntensityValue* weightfactor);
    }
    if (DenominatorSum == 0)
    {
      //medianImageFunction->EvaluateAtIndex(currentindex);
      //ImageIterator->SetPixel()
      ImageIterator.SetPixel(ImageIterator.Size() - 1, medianImageFunction->EvaluateAtIndex(currentindex));
    }
    else
    {
      // double previousval = image->GetPixel(currentindex);
      double newval = std::round(NumeratorSum / DenominatorSum);
      outputImageFilter->GetOutput()->SetPixel(currentindex, newval);
      //ImageIterator->SetPixel(currentindex, newval);
    }
  }
  //progressUpdate(100);
  return outputImageFilter->GetOutput();
}