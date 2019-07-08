/**
\file  PreprocessingPipelineClass.h

\brief File that holds the PreprocessingPipelineClass class

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/
#pragma once

//#include "iostream"
//#include "vtkImageData.h"
//#include "itkImage.h"
//#include "itkImageFileWriter.h"
//#include "itkConnectedThresholdImageFilter.h"
//#include "itkImageRegionIterator.h"
#include "itkBSplineControlPointImageFilter.h"
//#include "itkExpImageFilter.h"
//#include "itkImageRegionIterator.h"
#include "itkN3MRIBiasFieldCorrectionImageFilter.h"
//#include "itkN4BiasFieldCorrectionImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
//#include "itkShrinkImageFilter.h"
#include "itkMedianImageFunction.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkConnectedThresholdImageFilter.h"
//#include "itkNeighborhoodIterator.h"
//#include "itkMinimumMaximumImageCalculator.h"
//#include "itkAffineTransform.h"
//#include "itkRegularStepGradientDescentOptimizer.h"
//#include "itkMattesMutualInformationImageToImageMetric.h"
//#include "itkMultiResolutionImageRegistrationMethod.h"
//#include "itkNormalizedMutualInformationHistogramImageToImageMetric.h"
#include "itkFlatStructuringElement.h"
//#include "itkBinaryDilateImageFilter.h"
#include "itkHistogramMatchingImageFilter.h"
//
//#include "qmessagebox.h"
//#include "qstring.h"

#include "cbicaUtilities.h"

//#include "CAPTk.h"
#include "CaPTkDefines.h"
#include "CaPTkEnums.h"
#include <QMessageBox>
#include <QString>


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
public:
  typedef   TRegistration                              RegistrationType;
  typedef   RegistrationType *                         RegistrationPointer;
  typedef   itk::RegularStepGradientDescentOptimizer   OptimizerType;
  typedef   OptimizerType *                            OptimizerPointer;

  void Execute(itk::Object * object, const itk::EventObject & event)
  {
    if (!(itk::IterationEvent().CheckEvent(&event)))
    {
      return;
    }
    RegistrationPointer registration = dynamic_cast<RegistrationPointer>(object);

    OptimizerPointer optimizer = dynamic_cast<OptimizerPointer>(registration->GetModifiableOptimizer());
    /*std::cout << "-------------------------------------" << std::endl;
    std::cout << "MultiResolution Level : " << registration->GetCurrentLevel() << std::endl;
    std::cout << "mutual information value" << optimizer->GetValue()<<std::endl;*/
    std::cout << std::endl;
    if (registration->GetCurrentLevel() == 0)
    {
      optimizer->SetMaximumStepLength(0.0525);
      optimizer->SetMinimumStepLength(0.00001);
    }
    else
    {
      optimizer->SetMaximumStepLength(optimizer->GetMaximumStepLength() * 0.25);
      optimizer->SetMinimumStepLength(optimizer->GetMinimumStepLength() * 0.1);
    }
  }
  void Execute(const itk::Object *, const itk::EventObject &)
  {
    return;
  }
};

template<class TFilter>
class CommandIterationUpdate : public itk::Command
{
public:
  typedef CommandIterationUpdate   Self;
  typedef itk::Command             Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  itkNewMacro(Self);
protected:
  CommandIterationUpdate() {};
public:

  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute((const itk::Object *) caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    const TFilter * filter =
      dynamic_cast< const TFilter * >(object);
    if (typeid(event) != typeid(itk::IterationEvent))
    {
      return;
    }

    std::string msg = "Iteration " + std::to_string(filter->GetElapsedIterations()) + " (of " + std::to_string(filter->GetMaximumNumberOfIterations()) + ").  ";
    msg = msg + " Current convergence value = " + std::to_string(filter->GetCurrentConvergenceMeasurement()) + " (threshold = " + std::to_string(filter->GetConvergenceThreshold()) + ")";

    QMessageBox::warning(NULL, "Error", QString::fromStdString(msg), QMessageBox::Ok, NULL);
  }

};

class CommandIterationUpdateRegistration : public itk::Command
{
public:
  typedef  CommandIterationUpdateRegistration   Self;
  typedef  itk::Command             Superclass;
  typedef  itk::SmartPointer<Self>  Pointer;
  itkNewMacro(Self);
protected:
  CommandIterationUpdateRegistration() {};
public:
  typedef   itk::RegularStepGradientDescentOptimizer  OptimizerType;
  typedef   const OptimizerType *                     OptimizerPointer;
  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute((const itk::Object *)caller, event);
  }
  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    //OptimizerPointer optimizer = dynamic_cast<OptimizerPointer>(object);
    if (!(itk::IterationEvent().CheckEvent(&event)))
    {
      return;
    }
    //std::cout << optimizer->GetCurrentIteration() << "   ";
    //std::cout << optimizer->GetValue() << "   ";
    //std::cout << optimizer->GetCurrentPosition() << std::endl;
  }
};

class  PreprocessingPipelineClass
{

public:
  /*BiasFieldCorrectionClass biasObject;
  RegistrationClass registrationObject;
  SusanDenoisingClass susanObject;*/


public:
  PreprocessingPipelineClass();
  ~PreprocessingPipelineClass(){ };

  template<class ImageType>
  typename ImageType::Pointer PrepareTumroImageFromPoints(typename ImageType::Pointer inputimage, VectorVectorDouble &allDrawingPoints);

  template<class ImageType>
  typename ImageType::Pointer N4BiasCorrectionMethodCall(typename ImageType::Pointer inputimage);
  
  template<class ImageType>
  typename ImageType::Pointer DenoisingMethodCall(typename ImageType::Pointer image);


  template<class ImageType, class InternalImageType>
  typename itk::MultiResolutionImageRegistrationMethod<ImageType, ImageType>::Pointer Registration(typename ImageType::Pointer fixedImagePointer,
    typename ImageType::Pointer movingImagePointer);

  template<class ImageType>
  typename ImageType::Pointer ResampleTransform(typename itk::MultiResolutionImageRegistrationMethod<ImageType, ImageType>::Pointer registrationPointer, typename ImageType::Pointer fixedImagePointer, typename ImageType::Pointer movingImagePointer);

  template<class ImageType>
  typename ImageType::IndexType ApplyTransformationToFarPoints(typename ImageType::IndexType pointIndex, double transformation[3][3], double shifting[3]);

  template<class ImageType>
  typename ImageType::Pointer Edema3DSegmentationInGivenImage(typename ImageType::Pointer inputimage, VectorVectorDouble &allDrawingPoints);

  template<class ImageType>
  typename ImageType::Pointer Tumor3DSegmentationInGivenImage(typename ImageType::Pointer inputimage, VectorVectorDouble &allDrawingPoints);

  template<class ImageType>
  void GeodesicThresholding(typename ImageType::Pointer Geos, typename ImageType::Pointer Mask);

  template<class ImageType>
  void GetInitialGeodesicMask(typename ImageType::Pointer Geos, VectorVectorDouble &points);

  template<class ImageType>
  typename ImageType::Pointer HistogramMatching(const typename ImageType::Pointer inputImage, const typename ImageType::Pointer templateImage,
    const size_t histLevels = 100, const size_t noOfMatchPoints = 15, bool threshholdAtMeanIntensity = true);

};

template<class ImageType>
typename ImageType::Pointer PreprocessingPipelineClass::DenoisingMethodCall(const typename ImageType::Pointer image)
{
  typedef itk::ImageDuplicator<ImageType> DuplicatorType;
  typename DuplicatorType::Pointer outputImageFilter = DuplicatorType::New();
  outputImageFilter->SetInputImage(image);
  outputImageFilter->Update();

  double sigma = 0;
  double intensityVariationThreshold = 80;

  // typename ImageType::SizeType imageSize = image->GetLargestPossibleRegion().GetSize();

  typedef typename itk::MedianImageFunction< ImageType > MedianImageFunctionType;
  typename MedianImageFunctionType::Pointer medianImageFunction = MedianImageFunctionType::New();
  medianImageFunction->SetInputImage(outputImageFilter->GetOutput());


  typedef typename itk::NeighborhoodIterator<ImageType> NeighborhoodIteratorType;
  typename NeighborhoodIteratorType::RadiusType radius;
  for (unsigned int i = 0; i < ImageType::ImageDimension; i++)
    radius[i] = 1;
  NeighborhoodIteratorType ImageIterator(radius, outputImageFilter->GetOutput(), outputImageFilter->GetOutput()->GetLargestPossibleRegion());

  for (ImageIterator.GoToBegin(); !ImageIterator.IsAtEnd(); ++ImageIterator)
  {
    typename ImageType::IndexType index = ImageIterator.GetIndex();
    double NumeratorSum = 0;
    double DenominatorSum = 0;
    typename itk::NeighborhoodIterator<ImageType>::IndexType currentindex = ImageIterator.GetIndex();
    float CenterIntensityValue = ImageIterator.GetCenterPixel();

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
      medianImageFunction->EvaluateAtIndex(currentindex);
    else
    {
      // double previousval = image->GetPixel(currentindex);
      double newval = std::round(NumeratorSum / DenominatorSum);
      outputImageFilter->GetOutput()->SetPixel(currentindex, newval);
    }
  }
  return outputImageFilter->GetOutput();
}

template<class ImageType>
typename ImageType::Pointer PreprocessingPipelineClass::N4BiasCorrectionMethodCall(typename ImageType::Pointer inputimage)
{
  //typedef signed short RealType;
  //typedef itk::Image<RealType, 3> RealImageType;
  typedef itk::Image<unsigned char, 3> MaskImageType;

  typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkerType;
  typename ShrinkerType::Pointer shrinker = ShrinkerType::New();
  shrinker->SetInput(inputimage);
  shrinker->SetShrinkFactors(4);

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

  //typedef CommandIterationUpdate<CorrecterType> CommandType; 
  //typename CommandType::Pointer observer = CommandType::New();
  //correcter->AddObserver(itk::IterationEvent(), observer);
  try
  {
    correcter->Update();
  }
  catch (std::exception &e)
  {
    std::cerr << "Exception caught." << e.what() << "\n";
  }

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

  typedef itk::ExpImageFilter<ImageType, ImageType> ExpFilterType;
  typename ExpFilterType::Pointer expFilter = ExpFilterType::New();
  expFilter->SetInput(logField);
  expFilter->Update();

  typedef itk::DivideImageFilter<ImageType, ImageType, ImageType> DividerType;
  typename DividerType::Pointer divider = DividerType::New();
  divider->SetInput1(inputimage);
  divider->SetInput2(expFilter->GetOutput());
  divider->Update();

  return divider->GetOutput();
}

template<class ImageType, class InternalImageType>
typename itk::MultiResolutionImageRegistrationMethod<ImageType, ImageType>::Pointer PreprocessingPipelineClass::Registration(typename ImageType::Pointer fixedImagePointer,
  typename ImageType::Pointer movingImagePointer)
{
  typedef itk::AffineTransform <double, 3> TransformType;
  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
  typedef itk::MattesMutualInformationImageToImageMetric<ImageType, ImageType> MetricType;
  // typedef itk::NormalizedMutualInformationHistogramImageToImageMetric<ImageType, ImageType> MetricType2;
  typedef itk::MultiResolutionImageRegistrationMethod<ImageType, ImageType> RegistrationType;

  typedef itk::MultiResolutionPyramidImageFilter<ImageType, ImageType> FixedImagePyramidType;
  typedef itk::MultiResolutionPyramidImageFilter<ImageType, ImageType> MovingImagePyramidType;

  TransformType::Pointer transform = TransformType::New();
  OptimizerType::Pointer optimizer = OptimizerType::New();
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
  typename RegistrationType::Pointer   registrar = RegistrationType::New();
  typename MetricType::Pointer         metric = MetricType::New();
  typename FixedImagePyramidType::Pointer fixedImagePyramid = FixedImagePyramidType::New();
  typename MovingImagePyramidType::Pointer movingImagePyramid = MovingImagePyramidType::New();

  registrar->SetOptimizer(optimizer);
  registrar->SetInterpolator(interpolator);
  registrar->SetMetric(metric);
  registrar->SetTransform(transform);
  registrar->SetFixedImagePyramid(fixedImagePyramid);
  registrar->SetMovingImagePyramid(movingImagePyramid);

  // typedef itk::Image<ImageType> FixedImageType;
  // typedef itk::Image<ImageType> MovingImageType;
  //typedef itk::CastImageFilter<FixedImageType, InternalImageType> FixedCastFilterType;
  ////typedef itk::CastImageFilter<MovingImageType, InternalImageType> MovingCastFilterType;
  //FixedCastFilterType::Pointer fixedCaster = FixedCastFilterType::New();
  //MovingCastFilterType::Pointer movingCaster = MovingCastFilterType::New();

  //fixedCaster->SetInput(fixedImagePointer);
  //movingCaster->SetInput(movingImagePointer);

  registrar->SetFixedImage(fixedImagePointer);
  registrar->SetMovingImage(movingImagePointer);

  //	fixedCaster->Update();
  registrar->SetFixedImageRegion(fixedImagePointer->GetBufferedRegion());
  typedef typename RegistrationType::ParametersType ParametersType;

  ParametersType initialParameters(transform->GetNumberOfParameters());
  initialParameters[0] = 1;
  initialParameters[1] = 0.0;
  initialParameters[2] = 0.0;
  initialParameters[3] = 0.0;
  initialParameters[4] = 1;
  initialParameters[5] = 0.0;
  initialParameters[6] = 0.0;
  initialParameters[7] = 0.0;
  initialParameters[8] = 1;

  initialParameters[9] = 0.0;
  initialParameters[10] = 0.0;
  initialParameters[11] = 0.0;

  registrar->SetInitialTransformParameters(initialParameters);

  //parameters for metric1

  metric->SetNumberOfHistogramBins(128);
  metric->SetNumberOfSpatialSamples(30000);
  metric->ReinitializeSeed(76926294);
  metric->SetUseExplicitPDFDerivatives(false);

  optimizer->SetNumberOfIterations(100);
  optimizer->SetRelaxationFactor(0.9);
  optimizer->MaximizeOn();

  CommandIterationUpdateRegistration::Pointer observer = CommandIterationUpdateRegistration::New();
  optimizer->AddObserver(itk::IterationEvent(), observer);

  typedef RegistrationInterfaceCommand<RegistrationType> CommandType;
  typename CommandType::Pointer command = CommandType::New();
  registrar->AddObserver(itk::IterationEvent(), command);
  registrar->SetNumberOfLevels(5);
  try
  {
    registrar->Update();
    //std::cout << "stop condition" << registrar->GetOptimizer()->GetStopConditionDescription() << std::endl;
  }
  catch (itk::ExceptionObject &ex)
  {
    std::cout << ex << std::endl;
  }

  //---------------------------obtaining the registartion parameters---------------------
  ParametersType finalParameters = registrar->GetLastTransformParameters();
  //transformation[0][0] =  finalParameters[0];
  //transformation[0][1] = finalParameters[1];
  //transformation[0][2] = finalParameters[2];
  //transformation[1][0] = finalParameters[3];
  //transformation[1][1] = finalParameters[4];
  //transformation[1][2] = finalParameters[5];	
  //transformation[2][0] = finalParameters[6];
  //transformation[2][1] = finalParameters[7];
  //transformation[2][2] = finalParameters[8];
  //
  //shifting[0] = finalParameters[9];
  //shifting[1] = finalParameters[10];
  //shifting[2] = finalParameters[11];

  // unsigned int numberOfIterations = optimizer->GetCurrentIteration();
  // double bestValue = optimizer->GetValue();
  // TransformType::Pointer finalTransform = TransformType::New();

  return registrar;
}


template<class ImageType>
typename ImageType::Pointer PreprocessingPipelineClass::ResampleTransform(typename itk::MultiResolutionImageRegistrationMethod<ImageType, ImageType>::Pointer registrationPointer, typename ImageType::Pointer fixedImagePointer, typename ImageType::Pointer movingImagePointer)
{
  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
  typename ResampleFilterType::Pointer resample = ResampleFilterType::New();

  resample->SetTransform(registrationPointer->GetOutput()->Get());
  resample->SetInput(movingImagePointer);
  resample->SetSize(fixedImagePointer->GetLargestPossibleRegion().GetSize());
  resample->SetOutputOrigin(fixedImagePointer->GetOrigin());
  resample->SetOutputSpacing(fixedImagePointer->GetSpacing());
  resample->SetOutputDirection(fixedImagePointer->GetDirection());
  resample->SetDefaultPixelValue(0);
  resample->Update();

  return resample->GetOutput();
}


template<class ImageType>
typename ImageType::IndexType PreprocessingPipelineClass::ApplyTransformationToFarPoints(typename ImageType::IndexType pointIndex, double transformation[3][3], double shifting[3])
{
  typename ImageType::IndexType transformed_point;
  for (int i = 0; i < transformed_point.GetIndexDimension(); i++)
  {
    double sum = 0.0;
    for (int k = 0; k < transformed_point.GetIndexDimension(); k++)
      sum = sum + transformation[i][k] * pointIndex[k];
    transformed_point[i] = sum + shifting[i];
  }
  return transformed_point;
}
template<class ImageType>
typename ImageType::Pointer PreprocessingPipelineClass::Edema3DSegmentationInGivenImage(typename ImageType::Pointer inputimage, VectorVectorDouble &allDrawingPoints)
{
  typedef itk::ConnectedThresholdImageFilter<ImageType, ImageType > ConnectedFilterType; // typedef itk::VariableLengthVector< double > VariableLengthVectorType;
  typename ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();
  connectedThreshold->SetInput(inputimage);
  itk::VariableLengthVector< typename ImageType::PixelType > intensityValue;
  intensityValue.SetSize(allDrawingPoints.size());

  for (unsigned int i = 0; i < allDrawingPoints.size(); i++)
  {
    typename ImageType::IndexType  index;
    index[0] = allDrawingPoints[i][0];
    index[1] = allDrawingPoints[i][1];
    index[2] = allDrawingPoints[i][2];
    connectedThreshold->AddSeed(index);
    intensityValue[i] = inputimage.GetPointer()->GetPixel(index);
  }
  
  double mean = 0;
  for (unsigned int featureNo = 0; featureNo < allDrawingPoints.size(); featureNo++)
    mean = mean + intensityValue[featureNo];
  mean = mean / allDrawingPoints.size();

  // double stdDev = 0;
  double temp = 0;
  for (unsigned int featureNo = 0; featureNo < allDrawingPoints.size(); featureNo++)
    temp = temp + (intensityValue[featureNo] - mean)*(intensityValue[featureNo] - mean);
  // stdDev = std::sqrt(temp / (allDrawingPoints.size() - 1));

  
  typename ImageType::PixelType lowerThreshold = 300;
  typename ImageType::PixelType  upperThreshold = 500;
  connectedThreshold->SetLower(lowerThreshold);
  connectedThreshold->SetUpper(upperThreshold);
  connectedThreshold->SetReplaceValue(CAPTK::VOXEL_STATUS::ON);
  connectedThreshold->Update();

  typedef itk::FlatStructuringElement<3>	StructuringElementType;
  StructuringElementType::RadiusType radius;
  radius.Fill(5);
  StructuringElementType structuringElement = StructuringElementType::Ball(radius);
  typedef itk::BinaryDilateImageFilter< ImageType, ImageType, StructuringElementType > BinaryDilateImageFilterType;
  typename BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
  dilateFilter->SetInput(connectedThreshold->GetOutput());
  dilateFilter->SetKernel(structuringElement);
  dilateFilter->Update();

  typename ImageType::Pointer segmentedEdema = dilateFilter->GetOutput();

  //typedef itk::ImageFileWriter<ImageType> WriterType;
  //typename WriterType::Pointer writer = WriterType::New();
  //std::string filename = "segmentedEdema.nii.gz";
  //writer->SetFileName(filename);
  //writer->SetInput(dilateFilter->GetOutput());
  //writer->Update();

  return dilateFilter->GetOutput();
}
template<class ImageType>
typename ImageType::Pointer PreprocessingPipelineClass::PrepareTumroImageFromPoints(typename ImageType::Pointer inputimage, VectorVectorDouble &allDrawingPoints)
{
	typename ImageType::Pointer TumorImage = ImageType::New();
	TumorImage->CopyInformation(inputimage);
	TumorImage->SetRequestedRegion(inputimage->GetLargestPossibleRegion());
	TumorImage->SetBufferedRegion(inputimage->GetBufferedRegion());
	TumorImage->Allocate();
	TumorImage->FillBuffer(0);

	for (unsigned int i = 0; i < allDrawingPoints.size(); i++)
	{
		typename ImageType::IndexType  index;
		index[0] = allDrawingPoints[i][0];
		index[1] = allDrawingPoints[i][1];
		index[2] = allDrawingPoints[i][2];
		TumorImage->SetPixel(index, 1);
	}
	return TumorImage;
}
template<class ImageType>
typename ImageType::Pointer PreprocessingPipelineClass::Tumor3DSegmentationInGivenImage(typename ImageType::Pointer inputimage, VectorVectorDouble &allDrawingPoints)
{
  typedef itk::ConnectedThresholdImageFilter<ImageType, ImageType > ConnectedFilterType; // typedef itk::VariableLengthVector< double > VariableLengthVectorType;
  typename ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();
  connectedThreshold->SetInput(inputimage);
  itk::VariableLengthVector< typename ImageType::PixelType > intensityValue;
  intensityValue.SetSize(allDrawingPoints.size());

  for (unsigned int i = 0; i < allDrawingPoints.size(); i++)
  {
    typename ImageType::IndexType  index;
    index[0] = allDrawingPoints[i][0];
    index[1] = allDrawingPoints[i][1];
    index[2] = allDrawingPoints[i][2];
    connectedThreshold->AddSeed(index);
    intensityValue[i] = inputimage.GetPointer()->GetPixel(index);
  }
  double mean = 0;
  for (unsigned int featureNo = 0; featureNo < allDrawingPoints.size(); featureNo++)
    mean = mean + intensityValue[featureNo];
  mean = mean / allDrawingPoints.size();

   double stdDev = 0;
  double temp = 0;
  for (unsigned int featureNo = 0; featureNo < allDrawingPoints.size(); featureNo++)
    temp = temp + (intensityValue[featureNo] - mean)*(intensityValue[featureNo] - mean);
   stdDev = std::sqrt(temp / (allDrawingPoints.size() - 1));

  typename ImageType::PixelType lowerThreshold = mean-stdDev;
  typename ImageType::PixelType  upperThreshold = mean+stdDev;
  connectedThreshold->SetLower(lowerThreshold);
  connectedThreshold->SetUpper(upperThreshold);
  connectedThreshold->SetReplaceValue(175);
  connectedThreshold->Update();

  typedef itk::FlatStructuringElement<3>	StructuringElementType;
  StructuringElementType::RadiusType radius;
  radius.Fill(2);
  StructuringElementType structuringElement = StructuringElementType::Ball(radius);
  typedef itk::BinaryDilateImageFilter< ImageType, ImageType, StructuringElementType > BinaryDilateImageFilterType;
  typename BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
  dilateFilter->SetInput(connectedThreshold->GetOutput());
  dilateFilter->SetKernel(structuringElement);
  dilateFilter->Update();

  typename ImageType::Pointer segmentedTumor = dilateFilter->GetOutput();

  //typedef itk::ImageFileWriter<ImageType> WriterType;
  //typename WriterType::Pointer writer = WriterType::New();
  //std::string filename = "segmentedTumor.nii.gz";
  //writer->SetFileName(filename);
  //writer->SetInput(dilateFilter->GetOutput());
  //writer->Update();

  return segmentedTumor;
}


template< class ImageType >
typename ImageType::Pointer PreprocessingPipelineClass::HistogramMatching(const typename ImageType::Pointer inputImage, const typename ImageType::Pointer templateImage,
  const size_t histLevels, const size_t noOfMatchPoints, bool threshholdAtMeanIntensity)
{
  using HEFilterType = itk::HistogramMatchingImageFilter<ImageType, ImageType>;
  typename HEFilterType::Pointer IntensityEqualizeFilter = HEFilterType::New();
  IntensityEqualizeFilter->SetReferenceImage(templateImage);
  IntensityEqualizeFilter->SetInput(inputImage);
  IntensityEqualizeFilter->SetNumberOfHistogramLevels(histLevels);
  IntensityEqualizeFilter->SetNumberOfMatchPoints(noOfMatchPoints);
  if (threshholdAtMeanIntensity)
  {
    IntensityEqualizeFilter->ThresholdAtMeanIntensityOn();
  }
  IntensityEqualizeFilter->Update();

  return IntensityEqualizeFilter->GetOutput();
}