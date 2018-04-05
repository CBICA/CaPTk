/**
\file PerfusionPCA.h

This file holds the declaration of the class PerfusionPCA.

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or http://www.med.upenn.edu/sbia/software/license.html
*/
#pragma once

#include "iostream"
#include "vtkImageData.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkExpImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkN3MRIBiasFieldCorrectionImageFilter.h"
#include "itkN4BiasFieldCorrectionImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "itkMedianImageFunction.h"
#include "itkNeighborhoodIterator.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkFlatStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "qmessagebox.h"
#include "qstring.h"
#include "cbicaUtilities.h"
#include "ApplicationBase.h"
#include "FeatureReductionClass.h"
#include "CAPTk.h"
#include <itkExtractImageFilter.h>

#ifdef APP_BASE_CAPTK_H
#include "ApplicationBase.h"
#endif


class PerfusionPCA
#ifdef APP_BASE_CAPTK_H
	: public ApplicationBase
#endif
{
public:
  PerfusionPCA() {};
  ~PerfusionPCA() {};
  
  template<class ImageTypeFloat4D, class ImageTypeFloat3D>
  std::vector<typename ImageTypeFloat3D::Pointer> Run(typename ImageTypeFloat3D::Pointer maskImagePointerNifti, typename ImageTypeFloat4D::Pointer perfImagePointerNifti);

private:

};

template<class ImageTypeFloat4D,class ImageTypeFloat3D>
std::vector<typename ImageTypeFloat3D::Pointer> PerfusionPCA::Run(typename ImageTypeFloat3D::Pointer maskImagePointerNifti, typename ImageTypeFloat4D::Pointer perfImagePointerNifti)
{
	VectorVectorDouble pIntensities;
	messageUpdate("PCA Calculation");
	progressUpdate(0);
	//----------------------find the voxels of the mask------------------------------
	typedef itk::ImageRegionIteratorWithIndex <ImageTypeFloat3D> IteratorType;
	IteratorType maskIt(maskImagePointerNifti, maskImagePointerNifti->GetLargestPossibleRegion());
	std::vector<typename ImageTypeFloat3D::IndexType> qualifiedIndices;

	maskIt.GoToBegin();
	int mask_counter = 0;
	while (!maskIt.IsAtEnd())
	{
		if (maskIt.Get() == GLISTR_OUTPUT_LABELS::NONENHANCING)
		{
			mask_counter++;
			qualifiedIndices.push_back(maskIt.GetIndex());
		}
		++maskIt;
	}
	//--------------------------populate the covariance matrix --------------------------------
	typename ImageTypeFloat4D::RegionType region = perfImagePointerNifti->GetLargestPossibleRegion();
	vnl_matrix<double> revisedcovariancematrix;
	revisedcovariancematrix.set_size(mask_counter, region.GetSize()[3]);
	for (size_t i = 0; i < qualifiedIndices.size(); i++)
	{
		typename ImageTypeFloat4D::IndexType regionIndex;
		regionIndex[0] = qualifiedIndices[i][0];
		regionIndex[1] = qualifiedIndices[i][1];
		regionIndex[2] = qualifiedIndices[i][2];
		for (size_t j = 0; j < region.GetSize()[3]; j++)
		{
			regionIndex[3] = j;
			revisedcovariancematrix(i, j) = perfImagePointerNifti->GetPixel(regionIndex);
		}
		
	}

	//typename ImageTypeFloat4D::RegionType region = perfImagePointerNifti->GetLargestPossibleRegion();
 //   typename ImageTypeFloat4D::IndexType regionIndex;
	//typename ImageTypeFloat4D::SizeType regionSize;

	//vnl_matrix<double> covarianceMatrix;
	//vnl_matrix<double> revisedcovariancematrix;
	//covarianceMatrix.set_size(region.GetSize()[0] * region.GetSize()[1] * region.GetSize()[2], region.GetSize()[3]);
	//std::vector<typename ImageTypeFloat3D::IndexType> initialIndices;
	//std::vector<typename ImageTypeFloat3D::IndexType> qualifiedIndices;
	//std::vector<int> qualifiedIndicesNumbers;

	//regionSize[0] = region.GetSize()[0];
	//regionSize[1] = region.GetSize()[1];
	//regionSize[2] = region.GetSize()[2];
	//regionSize[3] = 0; // this is 0 because we need a 3D image
	//regionIndex[0] = 0;
	//regionIndex[1] = 0;
	//regionIndex[2] = 0;
	//regionIndex[3] = 0;

	//for (size_t i = 0; i < region.GetSize()[3]; i++)
	//{
	//	regionIndex[3] = i;
	//	typename ImageTypeFloat4D::RegionType desiredRegion(regionIndex, regionSize);
	//	auto filter = itk::ExtractImageFilter< ImageTypeFloat4D, ImageTypeFloat3D >::New();
	//	filter->SetExtractionRegion(desiredRegion);
	//	filter->SetInput(perfImagePointerNifti);

	//	filter->SetDirectionCollapseToIdentity(); // This is required.
	//	filter->Update();
	//	typedef itk::ImageRegionIteratorWithIndex <ImageTypeFloat3D> IteratorType;
	//	IteratorType tumorIt(filter->GetOutput(), filter->GetOutput()->GetLargestPossibleRegion());
	//	tumorIt.GoToBegin();
	//	int counter = 0;
	//	while (!tumorIt.IsAtEnd())
	//	{
	//		covarianceMatrix(counter, i) = tumorIt.Get();
	//		counter++;
	//		++tumorIt;

	//		if (i==0)
	//			initialIndices.push_back(tumorIt.GetIndex());
	//	}
	//}
	////trimming down the number of voxels
	//for (unsigned int index = 0; index < covarianceMatrix.rows(); index++)
	//{
	//	double mNearMean = 0.0;
	//	double mNearStd = 0.0;
	//	double temp = 0.0;
	//	for (unsigned int featureNo = 0; featureNo < covarianceMatrix.columns(); featureNo++)
	//		temp = temp + covarianceMatrix[index][featureNo];
	//	mNearMean = temp / covarianceMatrix.columns();

	//	temp = 0.0;
	//	for (unsigned int featureNo = 0; featureNo < covarianceMatrix.columns(); featureNo++)
	//		temp = temp + (covarianceMatrix[index][featureNo] - mNearMean)*(covarianceMatrix[index][featureNo] - mNearMean);
	//	mNearStd = std::sqrt(temp / (covarianceMatrix.columns() - 1));

	//	if (mNearStd != 0)
	//	{
	//		qualifiedIndices.push_back(initialIndices[index]);
	//		qualifiedIndicesNumbers.push_back(index);
	//	}
	//}
	//revisedcovariancematrix.set_size(qualifiedIndicesNumbers.size(), region.GetSize()[3]);
	//for (unsigned int index = 0; index < qualifiedIndicesNumbers.size(); index++)
	//{
	//	for (unsigned int featureNo = 0; featureNo < covarianceMatrix.columns(); featureNo++)
	//		revisedcovariancematrix(index, featureNo) = covarianceMatrix(qualifiedIndicesNumbers[index], featureNo);
	//}
	FeatureReductionClass obj;
	vtkSmartPointer<vtkTable> principalcomponents = obj.GetDiscerningPerfusionTimePoints(revisedcovariancematrix);

	//--------------------------------
	typename ImageTypeFloat4D::IndexType regionIndex;
	typename ImageTypeFloat4D::SizeType regionSize;
	std::vector<typename ImageTypeFloat3D::Pointer> principalcomponentImages;
	regionSize[0] = region.GetSize()[0];
	regionSize[1] = region.GetSize()[1];
	regionSize[2] = region.GetSize()[2];
	regionSize[3] = 0; 
	regionIndex[0] = 0;
	regionIndex[1] = 0;
	regionIndex[2] = 0;
	regionIndex[3] = 0;

	typename ImageTypeFloat4D::RegionType desiredRegion(regionIndex, regionSize);
	auto filter = itk::ExtractImageFilter< ImageTypeFloat4D, ImageTypeFloat3D >::New();
	filter->SetExtractionRegion(desiredRegion);
	filter->SetInput(perfImagePointerNifti);
	filter->SetDirectionCollapseToIdentity(); // This is required.
	filter->Update();

	for (unsigned int i = 0; i < region.GetSize()[3]; i++)
	{
		typename ImageTypeFloat3D::Pointer currentImage = ImageTypeFloat3D::New();
		currentImage->CopyInformation(filter->GetOutput());
		currentImage->SetRequestedRegion(filter->GetOutput()->GetLargestPossibleRegion());
		currentImage->SetBufferedRegion(filter->GetOutput()->GetBufferedRegion());
		currentImage->Allocate();
		currentImage->FillBuffer(0);

		for (unsigned int index = 0; index < qualifiedIndices.size(); index++)
			currentImage->SetPixel(qualifiedIndices[index], principalcomponents->GetValue(index,i).ToDouble());

		principalcomponentImages.push_back(currentImage);
	}
	return principalcomponentImages;
}


