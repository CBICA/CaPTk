/**
\file  AdaptiveGeodesicDistance.h
\brief The header file containing the Geodesic segmentation class, used to apply an adaptive geodesic transform
Library Dependecies: ITK 4.7+ <br>
Header Dependencies: cbicaUtilities.h, cbicaLogging.h
http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu
Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or http://www.med.upenn.edu/sbia/software/license.html
*/
#ifndef  H_CBICA_ADAPTIVE_GEODESIC_DISTANCE
#define  H_CBICA_ADAPTIVE_GEODESIC_DISTANCE

#include <iostream>
#include <limits.h>

#include "itkImage.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkExpImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "itkMedianImageFunction.h"
#include "itkNeighborhoodIterator.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkThresholdImageFilter.h"

/**
\namespace AdaptiveGeodesicDistance
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
namespace AdaptiveGeodesicDistance
{
	typedef int LabelsPixelType;

	typedef itk::Image< int, 3 > ImageTypeDefault; // Default image type

	template <class TImageType>
	using ImagePointer = typename TImageType::Pointer;

	/*For internal use*/
	template <class TImageType>
	double square(typename TImageType::PixelType x) {
		return static_cast<double>(x) * x;
	}

	/**
	Runs the Adaptive Geodesic Distance algorithm on a single image
	@param input the input image.
	@param labels an image the same size as input. A sample of labels for the input image. Only the pixels with value=labelOfInterest will be used.
	@param labelOfInterest For which label (of the possibly many) from the labels image to perform AGD
	@param limitAt255 if set to true the return image will be [0-255]
	@return the AGD result
	*/
	template<typename PixelType = int, unsigned int Dimensions = 3>
	ImagePointer< itk::Image<PixelType, Dimensions> >
	Run(const ImagePointer< itk::Image<PixelType, Dimensions> >       input,
        const ImagePointer< itk::Image<LabelsPixelType, Dimensions> > labels, 
        const int labelOfInterest = 1, const bool verbose = false, const bool limitAt255 = false)
	{
		typedef itk::Image<PixelType, Dimensions> ImageTypeGeodesic;

		ImagePointer<ImageTypeGeodesic> skipZerosGuideImage = ImageTypeGeodesic::New();
		skipZerosGuideImage->CopyInformation(input);
		skipZerosGuideImage->SetRequestedRegion(input->GetLargestPossibleRegion());
		skipZerosGuideImage->SetBufferedRegion(input->GetBufferedRegion());
		skipZerosGuideImage->Allocate();
		skipZerosGuideImage->FillBuffer(1);

		return Run<ImageTypeGeodesic>(input, skipZerosGuideImage, labels, labelOfInterest, verbose, limitAt255);
	}

	/**
	Runs the Adaptive Geodesic Distance algorithm on a single image
	@param input the input image.
	@param skipZerosGuideImage an image the same size as input. For the pixels that are zero no calculation will be done. (can be the input image itself)
	@param labels an image the same size as input. A sample of labels for the input image. Only the pixels with value=labelOfInterest will be used.
	@param labelOfInterest For which label (of the possibly many) from the labels image to perform AGD
	@param limitAt255 if set to true the return image will be [0-255]
	@return the AGD result
	*/
	template<typename PixelType = int, unsigned int Dimensions = 3>
	ImagePointer< itk::Image<PixelType, Dimensions> >
	Run(const ImagePointer< itk::Image<PixelType, Dimensions> >       input,
		const ImagePointer< itk::Image<PixelType, Dimensions> >       skipZerosGuideImage,
        const ImagePointer< itk::Image<LabelsPixelType, Dimensions> > labels,
        const int labelOfInterest = 1, const bool verbose = false, const bool limitAt255 = false)
	{
		static_assert((Dimensions == 2 || Dimensions == 3), "2D or 3D Images supported");

		typedef itk::Image<PixelType, Dimensions>                    ImageTypeGeodesic;
		typedef itk::Image<LabelsPixelType, Dimensions>              LabelsImageType;
		typedef ImagePointer<ImageTypeGeodesic>                      ImageGeodesicPointer;
		typedef ImagePointer<LabelsImageType>                        LabelsImagePointer;
		typedef itk::ImageRegionIteratorWithIndex<ImageTypeGeodesic> NormalIteratorIndexedGeo;
		typedef itk::ImageRegionIteratorWithIndex<LabelsImageType>   NormalIteratorIndexedLabels;
		typedef itk::NeighborhoodIterator<ImageTypeGeodesic>         NeighborhoodIteratorGeo;

		const double pixelTypeMaxVal = itk::NumericTraits< typename ImageTypeGeodesic::PixelType >::max();

		//-------------- Allocate output image ------------------------

		typename ImageTypeGeodesic::Pointer output = ImageTypeGeodesic::New();
		output->SetRegions(input->GetLargestPossibleRegion());
		output->SetRequestedRegion(input->GetLargestPossibleRegion());
		//output->SetBufferedRegion(input->GetBufferedRegion());
		output->Allocate();
		output->FillBuffer((limitAt255) ? 255 : static_cast<typename ImageTypeGeodesic::PixelType>(pixelTypeMaxVal - 1));
		output->SetDirection(input->GetDirection());
		output->SetOrigin(input->GetOrigin());
		output->SetSpacing(input->GetSpacing());

		//-------------- Initialize values for output image (0 at pixels of label of interest, max elsewhere [set above]) ------------------------

		NormalIteratorIndexedLabels labIter(labels, labels->GetLargestPossibleRegion());
		NormalIteratorIndexedGeo    outIter(output, output->GetLargestPossibleRegion());
		labIter.GoToBegin();
		outIter.GoToBegin();

		while (!labIter.IsAtEnd()) {
			if (labIter.Get() == labelOfInterest) {
				outIter.Set(0);
			}
			++labIter;
			++outIter;
		}

		//-------------- If there is need to use gamma ------------------------

		//typename ImageTypeGeodesic::Pointer Gamma = ImageTypeGeodesic::New();
		//Gamma->CopyInformation(input);
		//Gamma->SetRequestedRegion(input->GetLargestPossibleRegion());
		//Gamma->SetBufferedRegion(input->GetBufferedRegion());
		//Gamma->Allocate();
		//Gamma->FillBuffer(1);

		//--------------------------- For actual AGD --------------------------

		// Iterators that are used in the loops
		NormalIteratorIndexedGeo iterSkipZeros(skipZerosGuideImage, skipZerosGuideImage->GetLargestPossibleRegion());
		//NormalIteratorIndexedGeo iterGamma(Gamma, Gamma->GetLargestPossibleRegion());

		typename ImageTypeGeodesic::SizeType radius;
		
		for (int i = 0; i < ImageTypeGeodesic::ImageDimension; i++) {
			radius[i] = 1;
		}

		NeighborhoodIteratorGeo outNIter(radius, output, output->GetLargestPossibleRegion());
		NeighborhoodIteratorGeo inputNIter(radius, input, input->GetLargestPossibleRegion());

		// Variables that are used in the loops
		typename ImageTypeGeodesic::PixelType inpCenterPixel;
		//ImageTypeGeodesic::PixelType gamPixel;
		double arr[14];
		double minVal;

		//----- Backward pass -----

		if (verbose) {
			std::cout << "AdaptiveGeodesicDistance: \tForward pass\n";
		}
		//cbica::Logging(loggerFile, "AGD Main loops execution : Forward pass");

		outNIter.GoToBegin();
		inputNIter.GoToBegin();
		iterSkipZeros.GoToBegin();
		//iterGamma.GoToBegin();

		while (!iterSkipZeros.IsAtEnd())
		{
			if (iterSkipZeros.Get() != 0)
			{
				// Forward pass
				
				if (ImageTypeGeodesic::ImageDimension == 2)
				{
					// 2D

					inpCenterPixel = inputNIter.GetPixel(4);
					//gamPixel = iterGamma.Get();

					arr[4] = outNIter.GetCenterPixel();
					arr[0] = outNIter.GetPixel(0) + sqrt(2.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(0)));
					arr[1] = outNIter.GetPixel(1) + sqrt(1.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(1)));
					arr[2] = outNIter.GetPixel(2) + sqrt(2.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(2)));
					arr[3] = outNIter.GetPixel(3) + sqrt(1.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(3)));
					
					minVal = arr[4];
					for (int i = 0; i < 4; i++)
					{
						if (arr[i] < minVal) {
							minVal = arr[i];
						}
					}
				}
				else {
					// 3D
					
					inpCenterPixel = inputNIter.GetPixel(13);
					//gamPixel = iterGamma.Get();

					arr[13] = outNIter.GetCenterPixel();
					arr[0]  = outNIter.GetPixel(4)  + sqrt(1.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(4)));
					arr[1]  = outNIter.GetPixel(10) + sqrt(1.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(10)));
					arr[2]  = outNIter.GetPixel(12) + sqrt(1.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(12)));
					arr[3]  = outNIter.GetPixel(1)  + sqrt(2.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(1)));
					arr[4]  = outNIter.GetPixel(3)  + sqrt(2.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(3)));
					arr[5]  = outNIter.GetPixel(9)  + sqrt(2.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(9)));
					arr[6]  = outNIter.GetPixel(0)  + sqrt(3.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(0)));
					arr[7]  = outNIter.GetPixel(7)  + sqrt(2.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(7)));
					arr[8]  = outNIter.GetPixel(6)  + sqrt(3.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(6)));
					arr[9]  = outNIter.GetPixel(15) + sqrt(2.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(15)));
					arr[10] = outNIter.GetPixel(24) + sqrt(3.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(24)));
					arr[11] = outNIter.GetPixel(21) + sqrt(2.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(21)));
					arr[12] = outNIter.GetPixel(18) + sqrt(3.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(18)));

					minVal = arr[13];
					for (int i = 0; i < 13; i++)
					{
						if (arr[i] < minVal) {
							minVal = arr[i];
						}
					}
				}

				outNIter.SetCenterPixel(minVal);
			}
			++outNIter;
			++inputNIter;
			++iterSkipZeros;
			//++iterGamma;
		}

		//----- Backward pass -----

		if (verbose) {
			std::cout << "AdaptiveGeodesicDistance: \tBackward pass\n";
		}
		//cbica::Logging(loggerFile, "AGD Main loops execution : Backward pass");

		outNIter.GoToEnd();
		inputNIter.GoToEnd();
		--outNIter;
		--inputNIter;
		iterSkipZeros.GoToReverseBegin();
		//iterGamma.GoToReverseBegin();

		while (!iterSkipZeros.IsAtReverseEnd())
		{
			if (iterSkipZeros.Get() != 0)
			{
				// Backward pass

				if (ImageTypeGeodesic::ImageDimension == 2) 
				{
					// 2D

					inpCenterPixel = inputNIter.GetPixel(4);
					//gamPixelB = iterGamma.Get();

					arr[4] = outNIter.GetCenterPixel();
					arr[0] = outNIter.GetPixel(5) + sqrt(1.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(5)));
					arr[1] = outNIter.GetPixel(6) + sqrt(2.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(6)));
					arr[2] = outNIter.GetPixel(7) + sqrt(1.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(7)));
					arr[3] = outNIter.GetPixel(8) + sqrt(2.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(8)));
					
					minVal = arr[4];
					for (int i = 0; i < 4; i++)
					{
						if (arr[i] < minVal) {
							minVal = arr[i];
						}
					}
				}
				else
				{
					// 3D

					inpCenterPixel = inputNIter.GetPixel(13);
					//gamPixelB = iterGamma.Get();

					arr[13] = outNIter.GetCenterPixel();
					arr[0] = outNIter.GetPixel(22) + sqrt(1.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(22)));
					arr[1] = outNIter.GetPixel(16) + sqrt(1.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(16)));
					arr[2] = outNIter.GetPixel(14) + sqrt(1.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(14)));
					arr[3] = outNIter.GetPixel(25) + sqrt(2.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(25)));
					arr[4] = outNIter.GetPixel(23) + sqrt(2.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(23)));
					arr[5] = outNIter.GetPixel(17) + sqrt(2.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(17)));
					arr[6] = outNIter.GetPixel(26) + sqrt(3.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(26)));
					arr[7] = outNIter.GetPixel(19) + sqrt(2.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(19)));
					arr[8] = outNIter.GetPixel(20) + sqrt(3.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(20)));
					arr[9] = outNIter.GetPixel(11) + sqrt(2.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(11)));
					arr[10] = outNIter.GetPixel(2) + sqrt(3.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(2)));
					arr[11] = outNIter.GetPixel(5) + sqrt(2.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(5)));
					arr[12] = outNIter.GetPixel(8) + sqrt(3.0 + /*gamPixel * */ square<ImageTypeGeodesic>(inpCenterPixel - inputNIter.GetPixel(8)));
				
					minVal = arr[13];
					for (int i = 0; i < 13; i++)
					{
						if (arr[i] < minVal) {
							minVal = arr[i];
						}
					}
				}

				outNIter.SetCenterPixel(minVal);
			}
			--outNIter;
			--inputNIter;
			--iterSkipZeros;
			//--iterGamma;
		}

		return output;
	}
}

#endif