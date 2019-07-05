/**
\file N3BiasCorrection.h

This file holds the declaration of the class N3BiasCorrection.

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html
*/
#pragma once

#include "iostream"
//#include "vtkImageData.h"
//#include "itkImage.h"
//#include "itkImageFileWriter.h"
//#include "itkConnectedThresholdImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkExpImageFilter.h"
//#include "itkImageRegionIterator.h"
#include "itkN3MRIBiasFieldCorrectionImageFilter.h"
//#include "itkN4BiasFieldCorrectionImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkShrinkImageFilter.h"
//#include "itkMedianImageFunction.h"
//#include "itkNeighborhoodIterator.h"
//#include "itkMinimumMaximumImageCalculator.h"
//#include "itkFlatStructuringElement.h"
//#include "itkBinaryDilateImageFilter.h"
//#include "itkRegularStepGradientDescentOptimizer.h"
//#include "qmessagebox.h"
//#include "qstring.h"
#include "cbicaUtilities.h"
#include "cbicaLogging.h"
#include "CaPTkDefines.h"
#include "ApplicationBase.h"

//template<class TFilter>
//class CommandIterationUpdate : public itk::Command
//{
//public:
//  typedef CommandIterationUpdate   Self;
//  typedef itk::Command             Superclass;
//  typedef itk::SmartPointer<Self>  Pointer;
//  itkNewMacro(Self);
//protected:
//  CommandIterationUpdate() {};
//public:
//
//  void Execute(itk::Object *caller, const itk::EventObject & event)
//  {
//    Execute((const itk::Object *) caller, event);
//  }
//
//  void Execute(const itk::Object * object, const itk::EventObject & event)
//  {
//    const TFilter * filter =
//      dynamic_cast< const TFilter * >(object);
//    if (typeid(event) != typeid(itk::IterationEvent))
//    {
//      return;
//    }
//
//    std::string msg = "Iteration " + std::to_string(filter->GetElapsedIterations()) + " (of " + std::to_string(filter->GetMaximumNumberOfIterations()) + ").  ";
//    msg = msg + " Current convergence value = " + std::to_string(filter->GetCurrentConvergenceMeasurement()) + " (threshold = " + std::to_string(filter->GetConvergenceThreshold()) + ")";
//
//    QMessageBox::warning(NULL, "Error", QString::fromStdString(msg), QMessageBox::Ok, NULL);
//  }
//
//};

/**
\class N3BiasCorrection

\brief This class implements the N3BiasCorrection algorithm.

}
*/
class N3BiasCorrection : public ApplicationBase
{

public:
  N3BiasCorrection() {};
  ~N3BiasCorrection() {};
  
  template<class ImageType>
typename ImageType::Pointer Run(typename ImageType::Pointer inputimage);

private:

};

template<class ImageType>
typename ImageType::Pointer N3BiasCorrection::Run(typename ImageType::Pointer inputimage)
{
  messageUpdate("N4 Bias Correction");
  progressUpdate(0);
  //typedef signed short RealType;
  //typedef itk::Image<RealType, 3> RealImageType;
  typedef itk::Image<unsigned char, 3> MaskImageType;

  typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkerType;
  typename ShrinkerType::Pointer shrinker = ShrinkerType::New();
  shrinker->SetInput(inputimage);
  shrinker->SetShrinkFactors(4);
  //Parameter1: Shrink factor
  typename MaskImageType::Pointer maskImage = NULL;
  if (!maskImage)
  {
    typedef itk::OtsuThresholdImageFilter<ImageType, MaskImageType> ThresholderType;
    typename ThresholderType::Pointer otsu = ThresholderType::New();
    otsu->SetInput(inputimage);
    otsu->SetNumberOfHistogramBins(200);
    otsu->SetInsideValue(0);
    otsu->SetOutsideValue(1);
    otsu->Update();
    maskImage = otsu->GetOutput();
  }
  //parameter2: Number of Bins
  
  progressUpdate(25);

  typedef itk::ShrinkImageFilter<MaskImageType, MaskImageType> MaskShrinkerType;
  typename MaskShrinkerType::Pointer maskshrinker = MaskShrinkerType::New();
  maskshrinker->SetInput(maskImage);
  maskshrinker->SetShrinkFactors(4);
  shrinker->Update();
  maskshrinker->Update();

  typedef itk::N3MRIBiasFieldCorrectionImageFilter<ImageType, MaskImageType, ImageType> CorrecterType;
  //typedef itk::N4BiasFieldCorrectionImageFilter<ImageType, MaskImageType, ImageType> CorrecterType;
  typename CorrecterType::Pointer correcter = CorrecterType::New();
  correcter->SetInput(shrinker->GetOutput());
  correcter->SetMaskImage(maskshrinker->GetOutput());
  //CorrecterType::VariableSizeArrayType max_iterations;
  //max_iterations.Fill(100);
  correcter->SetMaximumNumberOfIterations(/*max_iterations*/100);
  //correcter->SetNumberOfFittingLevels(50);

  //parameter3: Number of iterations

  //typedef CommandIterationUpdate<CorrecterType> CommandType; 
  //typename CommandType::Pointer observer = CommandType::New();
  //correcter->AddObserver(itk::IterationEvent(), observer);
  try
  {
    correcter->Update();
  }
  catch (std::exception &e)
  {
    cbica::Logging(loggerFile, "Error caught: " + std::string(e.what()));
  }

  progressUpdate(50);

  typedef itk::BSplineControlPointImageFilter<typename CorrecterType::BiasFieldControlPointLatticeType, typename CorrecterType::ScalarImageType> BSplinerType;
  typename BSplinerType::Pointer bspliner = BSplinerType::New();
  bspliner->SetInput(correcter->GetLogBiasFieldControlPointLattice());
  bspliner->SetSplineOrder(correcter->GetSplineOrder());
  bspliner->SetSize(inputimage->GetLargestPossibleRegion().GetSize());
  bspliner->SetOrigin(inputimage->GetOrigin());
  bspliner->SetDirection(inputimage->GetDirection());
  bspliner->SetSpacing(inputimage->GetSpacing());
  bspliner->Update();

  typename ImageType::Pointer logField = ImageType::New();
  logField->SetOrigin(bspliner->GetOutput()->GetOrigin());
  logField->SetSpacing(bspliner->GetOutput()->GetSpacing());
  logField->SetRegions(bspliner->GetOutput()->GetLargestPossibleRegion().GetSize());
  logField->SetDirection(bspliner->GetOutput()->GetDirection());
  logField->Allocate();

  itk::ImageRegionIterator<typename CorrecterType::ScalarImageType> ItB(bspliner->GetOutput(), bspliner->GetOutput()->GetLargestPossibleRegion());
  itk::ImageRegionIterator<ImageType> ItF(logField, logField->GetLargestPossibleRegion());
  for (ItB.GoToBegin(), ItF.GoToBegin(); !ItB.IsAtEnd(); ++ItB, ++ItF)
  {
    ItF.Set(ItB.Get()[0]);
  }

  progressUpdate(75);

  typedef itk::ExpImageFilter<ImageType, ImageType> ExpFilterType;
  typename ExpFilterType::Pointer expFilter = ExpFilterType::New();
  expFilter->SetInput(logField);
  expFilter->Update();

  typedef itk::DivideImageFilter<ImageType, ImageType, ImageType> DividerType;
  typename DividerType::Pointer divider = DividerType::New();
  divider->SetInput1(inputimage);
  divider->SetInput2(expFilter->GetOutput());
  divider->Update();


  const typename ImageType::SpacingType spacing = inputimage->GetSpacing();
  const typename ImageType::PointType origin = inputimage->GetOrigin();
  const typename ImageType::DirectionType direction = inputimage->GetDirection();
  const typename ImageType::RegionType region = inputimage->GetLargestPossibleRegion();
  typename ImageType::Pointer outputImage = ImageType::New();

  outputImage->SetSpacing(spacing);
  outputImage->SetLargestPossibleRegion(region);
  outputImage->SetBufferedRegion(region);
  outputImage->Allocate();
  outputImage->SetOrigin(origin);
  outputImage->SetDirection(direction);

  itk::ImageRegionIterator<ImageType> outputIterator(outputImage, outputImage->GetLargestPossibleRegion());
  itk::ImageRegionIterator<ImageType> dividerIterator(divider->GetOutput(), divider->GetOutput()->GetLargestPossibleRegion());

  for (outputIterator.GoToBegin(), dividerIterator.GoToBegin(); !dividerIterator.IsAtEnd(); ++dividerIterator, ++outputIterator)
  {
    outputIterator.Set(dividerIterator.Get());
  }

  progressUpdate(100);
  return outputImage;
}


