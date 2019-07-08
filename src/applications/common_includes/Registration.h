/**
\file Registration.h

This file holds the declaration of the class Registration.

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html
*/
#pragma once

#include <fstream>
#include "itkImageRegistrationMethodv4.h"
#include "itkAffineTransform.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizerv4.h"
#include "itkConjugateGradientLineSearchOptimizerv4.hxx"
#include "itkHistogramMatchingImageFilter.h"
#include "ApplicationBase.h"


template <typename TRegistration>
class RegistrationInterfaceCommand : public itk::Command
{
public:
  typedef  RegistrationInterfaceCommand   Self;
  typedef  itk::Command                   Superclass;
  typedef  itk::SmartPointer<Self>        Pointer;
  itkNewMacro(Self);
protected:
  RegistrationInterfaceCommand() {};
  // Software Guide : EndCodeSnippet
  // Software Guide : BeginLatex
  //
  // For convenience, we declare types useful for converting pointers
  // in the \code{Execute()} method.
  //
  // Software Guide : EndLatex
  // Software Guide : BeginCodeSnippet
public:
  typedef   TRegistration      RegistrationType;
  void Execute(itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
  {
    if (!(itk::MultiResolutionIterationEvent().CheckEvent(&event)))
    {
      return;
    }

    RegistrationType* registration = static_cast<RegistrationType *>(object);
    unsigned int currentLevel = registration->GetCurrentLevel();
    typename RegistrationType::ShrinkFactorsPerDimensionContainerType shrinkFactors =
      registration->GetShrinkFactorsPerDimension(currentLevel);
    typename RegistrationType::SmoothingSigmasArrayType smoothingSigmas =
      registration->GetSmoothingSigmasPerLevel();
    //std::cout << "-------------------------------------" << std::endl;
    //std::cout << " Current multi-resolution level = " << currentLevel << std::endl;
    //std::cout << "    shrink factor = " << shrinkFactors << std::endl;
    //std::cout << "    smoothing sigma = " << smoothingSigmas[currentLevel] << std::endl;
    //std::cout << std::endl;
  }
    void Execute(const itk::Object *, const itk::EventObject &) ITK_OVERRIDE
  { return; }
};

class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef  itk::SmartPointer<Self>  Pointer;
  itkNewMacro(Self);
protected:
  CommandIterationUpdate() : m_CumulativeIterationIndex(0) {};
public:
  typedef   itk::GradientDescentOptimizerv4Template<double>   OptimizerType;
  typedef   const OptimizerType *                               OptimizerPointer;
  void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
  {
    Execute((const itk::Object *)caller, event);
  }
    void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
  {
    OptimizerPointer optimizer = static_cast< OptimizerPointer >(object);
    if (!(itk::IterationEvent().CheckEvent(&event)))
    {
      return;
    }
    std::cout << optimizer->GetCurrentIteration() << "   ";
    std::cout << optimizer->GetValue() << "   ";
    std::cout << optimizer->GetCurrentPosition() << "   ";



    char filename[] = "Iterations.csv";
    std::fstream FileToWorkWith;

    FileToWorkWith.open(filename, std::fstream::in | std::fstream::out | std::fstream::app);


    // If file does not exist, Create new file
    if (!FileToWorkWith)
    {
      std::cout << "Cannot open file, file does not exist. Creating new file..";

      FileToWorkWith.open(filename, std::fstream::in | std::fstream::out | std::fstream::trunc);
      FileToWorkWith << m_CumulativeIterationIndex << "," << optimizer->GetValue() << "\n";
      FileToWorkWith.close();

    }
    else
    {    // use existing file
      FileToWorkWith << m_CumulativeIterationIndex << "," << optimizer->GetValue() << "\n";
      FileToWorkWith.close();

    }
    std::cout << m_CumulativeIterationIndex++ << std::endl;
  }
private:
  unsigned int m_CumulativeIterationIndex;
};


/**
\class Registration

\brief Class that handles affine registration between a fixed and moving image

This class uses standard ITK-based filters/optimizers/etc. for performing this multi-resolution registration.

Optimizer:    RegularStepGradientDescentOptimizer
Interpolator: LinearInterpolateImageFunction
Metric:       MattesMutualInformationImageToImageMetric

*/
const unsigned int  Dimension = 3;
typedef  float           PixelType;
typedef itk::Image< PixelType, Dimension > ImageType;

class Registration : public ApplicationBase
{

public:
  Registration() {};
  ~Registration() {};

  template<class ImageType/*do we need this second template?*/>
  typename itk::CastImageFilter< ImageType, ImageType >::Pointer Run(typename ImageType::Pointer fixedImagePointer,
    typename ImageType::Pointer movingImagePointer);


private:
  inline void SetLongRunning(bool longRunning);

};

template<class ImageType >
typename itk::CastImageFilter< ImageType, ImageType >::Pointer Registration::Run(typename ImageType::Pointer fixedImagePointer,
  typename ImageType::Pointer movingImagePointer)
{

  ///**preprocessing input images**/
  const int numHistogramLevels = 45, numMatchPoints = 30;

  typedef itk::HistogramMatchingImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer normalizeFilter = FilterType::New();
  normalizeFilter->SetInput(movingImagePointer);
  normalizeFilter->SetReferenceImage(fixedImagePointer);
  normalizeFilter->ThresholdAtMeanIntensityOn();
  normalizeFilter->SetNumberOfHistogramLevels(numHistogramLevels);
  normalizeFilter->SetNumberOfMatchPoints(numMatchPoints);
  normalizeFilter->Update();

  /***************************************************************/
  /****Initial translation transform ****/
  typedef itk::TranslationTransform< double, Dimension >      TTransformType;
  typedef itk::AffineTransform< double, Dimension >  ATransformType;
  typedef itk::RegularStepGradientDescentOptimizerv4<double>    TOptimizerType;
  typedef itk::MattesMutualInformationImageToImageMetricv4 <ImageType, ImageType > MetricType;
  typedef itk::ImageRegistrationMethodv4 <ImageType, ImageType, TTransformType > TRegistrationType;
  


  TOptimizerType::Pointer      transoptimizer = TOptimizerType::New();
  typename MetricType::Pointer         transmetric = MetricType::New();
  typename TRegistrationType::Pointer   transregistration = TRegistrationType::New();
  TTransformType::Pointer     ttransform = TTransformType::New();

  transregistration->SetOptimizer(transoptimizer);
  transregistration->SetMetric(transmetric);
  typedef TOptimizerType::ParametersType ParametersType;
  ParametersType initialParameters(ttransform->GetNumberOfParameters());
  initialParameters[0] = 0;
  initialParameters[1] = 0;
  initialParameters[2] = 0;

  ttransform->SetParameters(initialParameters);
  transregistration->SetMovingInitialTransform(ttransform);

  typedef itk::CompositeTransform< double, Dimension >  CompositeTransformType;
  CompositeTransformType::Pointer  compositeTransform = CompositeTransformType::New();
  compositeTransform->AddTransform(ttransform);
  transregistration->SetFixedImage(normalizeFilter->GetReferenceImage());
  transregistration->SetMovingImage(normalizeFilter->GetSourceImage());

  const unsigned int numberOfLevels1 = 1;
  typename TRegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel1;
  shrinkFactorsPerLevel1.SetSize(numberOfLevels1);
  shrinkFactorsPerLevel1[0] = 3;
  typename TRegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel1;
  smoothingSigmasPerLevel1.SetSize(numberOfLevels1);
  smoothingSigmasPerLevel1[0] = 2;
  transregistration->SetNumberOfLevels(numberOfLevels1);
  transregistration->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel1);
  transregistration->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel1);

  transmetric->SetNumberOfHistogramBins(24);

  transoptimizer->SetNumberOfIterations(200);
  transoptimizer->SetRelaxationFactor(0.5);
  transoptimizer->SetLearningRate(2);
  transoptimizer->SetMinimumStepLength(0.001);

  CommandIterationUpdate::Pointer observer1 = CommandIterationUpdate::New();
  transoptimizer->AddObserver(itk::IterationEvent(), observer1);

  typedef RegistrationInterfaceCommand<TRegistrationType> TranslationCommandType;
  typename TranslationCommandType::Pointer command1 = TranslationCommandType::New();
  transregistration->AddObserver(itk::MultiResolutionIterationEvent(), command1);

  try
  {
    transregistration->Update();
    std::cout << "Optimizer stop condition: "
      << transregistration->GetOptimizer()->GetStopConditionDescription()
      << std::endl;
  }
  catch (itk::ExceptionObject & err)
  {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return EXIT_FAILURE;
  }
  compositeTransform->AddTransform(
    transregistration->GetModifiableTransform());
  /**End of initial transform**/


  /** Affine transform**/
  using FixedImageType = ImageType;
  using MovingImageType = ImageType;
  typedef itk::ConjugateGradientLineSearchOptimizerv4Template<double > AOptimizerType;
  typedef itk::ImageRegistrationMethodv4<FixedImageType, MovingImageType, ATransformType > ARegistrationType;
  AOptimizerType::Pointer      affineoptimizer = AOptimizerType::New();
  typename MetricType::Pointer          affinemetric = MetricType::New();
  typename ARegistrationType::Pointer   affineregistration = ARegistrationType::New();
  ATransformType::Pointer affinetransform = ATransformType::New();

  affineregistration->SetOptimizer(affineoptimizer);
  affineregistration->SetMetric(affinemetric);
  affineregistration->SetMovingInitialTransform(compositeTransform);

  affineregistration->SetFixedImage(normalizeFilter->GetReferenceImage());
  affineregistration->SetMovingImage(normalizeFilter->GetSourceImage());
  affinemetric->SetNumberOfHistogramBins(24);
  
  using SpacingType = typename FixedImageType::SpacingType;
  using OriginType = typename FixedImageType::PointType;
  using RegionType = typename FixedImageType::RegionType;
  using SizeType =  typename FixedImageType::SizeType;
  
  typename FixedImageType::Pointer fixedImage = (FixedImageType*)normalizeFilter->GetReferenceImage();

  const SpacingType fixedSpacing = fixedImage->GetSpacing();
  const OriginType  fixedOrigin = fixedImage->GetOrigin();
  const RegionType  fixedRegion = fixedImage->GetLargestPossibleRegion();
  const SizeType    fixedSize = fixedRegion.GetSize();
  ATransformType::InputPointType centerFixed;
  centerFixed[0] = fixedOrigin[0] + fixedSpacing[0] * fixedSize[0] / 2.0;
  centerFixed[1] = fixedOrigin[1] + fixedSpacing[1] * fixedSize[1] / 2.0;
  centerFixed[2] = fixedOrigin[2] + fixedSpacing[2] * fixedSize[2] / 2.0;

  const unsigned int numberOfFixedParameters = affinetransform->GetFixedParameters().Size();
  ATransformType::ParametersType fixedParameters(numberOfFixedParameters);
  for (unsigned int i = 0; i < numberOfFixedParameters; ++i)
  {
    fixedParameters[i] = centerFixed[i];
  }
  affinetransform->SetFixedParameters(fixedParameters);

  affineregistration->SetInitialTransform(affinetransform);

  typedef itk::RegistrationParameterScalesFromPhysicalShift<MetricType> ScalesEstimatorType;
  typename ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();
  scalesEstimator->SetMetric(affinemetric);
  scalesEstimator->SetTransformForward(true);
  affineoptimizer->SetScalesEstimator(scalesEstimator);

  affineoptimizer->SetDoEstimateLearningRateOnce(true);
  affineoptimizer->SetDoEstimateLearningRateAtEachIteration(false);
  affineoptimizer->SetLowerLimit(0);
  affineoptimizer->SetUpperLimit(2);
  affineoptimizer->SetEpsilon(0.2);
  affineoptimizer->SetNumberOfIterations(200);
  affineoptimizer->SetMinimumConvergenceValue(1e-6);
  affineoptimizer->SetConvergenceWindowSize(5);

  CommandIterationUpdate::Pointer observer2 = CommandIterationUpdate::New();
  affineoptimizer->AddObserver(itk::IterationEvent(), observer2);



  const unsigned int numberOfLevels2 = 2;
  typename ARegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel2;
  shrinkFactorsPerLevel2.SetSize(numberOfLevels2);
  shrinkFactorsPerLevel2[0] = 2;
  shrinkFactorsPerLevel2[1] = 1;
  typename ARegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel2;
  smoothingSigmasPerLevel2.SetSize(numberOfLevels2);
  smoothingSigmasPerLevel2[0] = 1;
  smoothingSigmasPerLevel2[1] = 0;
  affineregistration->SetNumberOfLevels(numberOfLevels2);
  affineregistration->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel2);
  affineregistration->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel2);

  typedef RegistrationInterfaceCommand<ARegistrationType> AffineCommandType;
  typename AffineCommandType::Pointer command2 = AffineCommandType::New();
  affineregistration->AddObserver(itk::MultiResolutionIterationEvent(), command2);
  //itk::TimeProbesCollectorBase chronometer;
  try
  {
    //chronometer.Start("Registration");
    affineregistration->Update();
    //chronometer.Stop("Registration");
    //std::cout << "Optimizer stop condition: "
      //<< affineregistration->GetOptimizer()->GetStopConditionDescription()
      //<< std::endl;
  }
  catch (itk::ExceptionObject & err)
  {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return EXIT_FAILURE;
  }
  //chronometer.Report(std::cout);


  compositeTransform->AddTransform(affineregistration->GetModifiableTransform());
  //std::cout << "\nInitial parameters of the registration process: " << std::endl
  //  << ttransform->GetParameters() << std::endl;
  //std::cout << "\nTranslation parameters after registration: " << std::endl
  //  << transoptimizer->GetCurrentPosition() << std::endl
  //  << " Last LearningRate: " << transoptimizer->GetCurrentStepLength() << std::endl;
  //std::cout << "\nAffine parameters after registration: " << std::endl
  //  << affineoptimizer->GetCurrentPosition() << std::endl
  //  << " Last LearningRate: " << affineoptimizer->GetLearningRate() << std::endl;

  /**End of affine registration**/

  /*** Resample and write back registered image***/
  typedef itk::ResampleImageFilter<MovingImageType, FixedImageType >ResampleFilterType;
  typename ResampleFilterType::Pointer resample = ResampleFilterType::New();
  resample->SetTransform(compositeTransform);
  resample->SetInput(normalizeFilter->GetSourceImage()/*movingImageReader->GetOutput()*/);
  resample->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
  resample->SetOutputOrigin(fixedImage->GetOrigin());
  resample->SetOutputSpacing(fixedImage->GetSpacing());
  resample->SetOutputDirection(fixedImage->GetDirection());
  resample->SetDefaultPixelValue(100);

  typedef itk::CastImageFilter< ImageType, ImageType > CastFilterType;

  typename CastFilterType::Pointer  caster = CastFilterType::New();
  caster->SetInput(resample->GetOutput());
  caster->Update();
  return caster->GetOutput();
}
