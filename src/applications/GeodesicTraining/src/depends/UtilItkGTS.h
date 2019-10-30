#ifndef H_CBICA_ITK_UTIL_GTS
#define H_CBICA_ITK_UTIL_GTS

#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkExceptionObject.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkNormalizeImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkResampleImageFilter.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkBilateralImageFilter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkBinaryBallStructuringElement.h"

namespace GeodesicTrainingSegmentation
{
	namespace ItkUtilGTS
	{
		/**Normalizes image values to be inside [0,255]*/
		template<class UImageType>
		typename UImageType::Pointer normalizeImage(typename UImageType::Pointer image, typename UImageType::PixelType max = 255)
		{
			auto filter = itk::RescaleIntensityImageFilter< UImageType, UImageType >::New();
			filter->SetInput(image);
			filter->SetOutputMinimum(0);
			filter->SetOutputMaximum(max);
			filter->Update();
			return filter->GetOutput();
		}

		/**Normalizes image values to be inside [0,255]*/
		template<class UImageType>
		void normalizeImageVector(std::vector<typename UImageType::Pointer>& images)
		{
			for (size_t i = 0; i < images.size(); i++) {
				images[i] = normalizeImage<UImageType>(images[i]);
			}
		}

		/**Normalizes image values to be inside [0,255]*/
		template<typename KeysType, class UImageType>
		void normalizeImageMap(std::map< KeysType, typename UImageType::Pointer>& map)
		{
			if (map.size() == 0) { return; }
			
			// Find a list of the different keys
			std::vector<KeysType> keys;
			keys.reserve(map.size());
			for (auto const& imap : map) {
				keys.push_back(imap.first);
			}

			// Normalize
			for (auto key : keys) {
				map[key] = normalizeImage<UImageType>(map[key]);
			}
		}

		/**Calculates the maximum pixel value of an image*/
		template<class UImageType>
		typename UImageType::PixelType getImageMaximum(typename UImageType::Pointer image)
		{
			auto imageCalculatorFilter = itk::MinimumMaximumImageCalculator<UImageType>::New();
			imageCalculatorFilter->SetImage(image);
			imageCalculatorFilter->ComputeMaximum();
			return imageCalculatorFilter->GetMaximum();
		}

		/**Changes the type of an image and rescales it to be inside [0,255]*/
		template<class TImageType, class UImageType>
		typename UImageType::Pointer castAndRescaleImage(typename TImageType::Pointer input)
		{
			typedef itk::CastImageFilter<TImageType, UImageType> CastFilterType;

			typename CastFilterType::Pointer castFilter = CastFilterType::New();
			castFilter->SetInput(input);

			auto filter = itk::RescaleIntensityImageFilter< UImageType, UImageType >::New();
			filter->SetInput(castFilter->GetOutput());
			filter->SetOutputMinimum(0);
			filter->SetOutputMaximum(255);
			filter->Update();
			return filter->GetOutput();
		}

		/**Create a zero image with exactly the same dimensions as another image*/
		template <class TImageType, class UImageType = TImageType>
		typename UImageType::Pointer initializeOutputImageBasedOn(typename TImageType::Pointer image)
		{
			typename UImageType::Pointer res = UImageType::New();

			res->SetRegions(image->GetLargestPossibleRegion());
			res->SetRequestedRegion(image->GetRequestedRegion());
			res->SetBufferedRegion(image->GetBufferedRegion());
			res->Allocate();
			res->FillBuffer(0);
			res->SetDirection(image->GetDirection());
			res->SetOrigin(image->GetOrigin());
			res->SetSpacing(image->GetSpacing());

			return res;
		}

		/**Normalize an image by setting its mean to 0 and variance to 1*/
		template<class UImageType, class VImageType = UImageType>
		static typename VImageType::Pointer statisticalImageNormalization(typename UImageType::Pointer image, typename VImageType::PixelType normMax = 255)
		{
			auto normalizeFilter = itk::NormalizeImageFilter<UImageType, VImageType>::New();
			normalizeFilter->SetInput(image);
			normalizeFilter->Update();
			return normalizeImage<VImageType>(normalizeFilter->GetOutput(), normMax);
		}

		/**Normalize an image by setting its mean to 0 and variance to 1*/
		template<class UImageType, class VImageType = UImageType>
		static void statisticalImageVectorNormalization(std::vector<typename UImageType::Pointer>& images, typename VImageType::PixelType normMax = 255)
		{
			for (size_t i = 0; i < images.size(); i++) {
				images[i] = statisticalImageNormalization<UImageType, VImageType>(images[i], normMax);
			}
		}

		/**Normalize an image by setting its mean to 0 and variance to 1*/
		template<typename KeysType, class UImageType, class VImageType = UImageType>
		static void statisticalImageMapNormalization(std::map< KeysType, typename UImageType::Pointer>& map, typename VImageType::PixelType normMax = 255)
		{
			if (map.size() == 0) { return; }

			// Find a list of the different keys
			std::vector<KeysType> keys;
			keys.reserve(map.size());
			for (auto const& imap : map) {
				keys.push_back(imap.first);
			}

			// Statistical Normalize
			for (auto key : keys) {
				map[key] = statisticalImageNormalization<UImageType, VImageType>(map[key], normMax);
			}
		}

		template<class TImageType>
		static typename TImageType::Pointer
		changeImageSpacing(typename TImageType::Pointer input, 
		                   typename TImageType::SpacingType outputSpacing,
		                   bool isLabelsImage = false, 
		                   bool forceCertainSize  = false, 
		                   typename TImageType::SizeType forcedSize = typename TImageType::SizeType(), 
		                   bool forceCertainOrigin = false,
		                   typename TImageType::PointType forcedOrigin = typename TImageType::PointType())
		{
			typename TImageType::SizeType inputSize = input->GetLargestPossibleRegion().GetSize();
			typename TImageType::SpacingType inputSpacing = input->GetSpacing();
			const unsigned int dimensions = inputSize.GetSizeDimension();

			bool isAnythingDifferent = false;
			for (unsigned int i=0; i < dimensions; i++) {
				if (inputSpacing[i] != outputSpacing[i]) {
					isAnythingDifferent = true;
					break;
				}
			}
			if (!isAnythingDifferent) { return input; }

			// Print input information
			std::cout << "\tInput image size:    " << inputSize << "\n";
			std::cout << "\tInput image spacing: " << inputSpacing << "\n";

			// Create the new image
			std::cout << "\tOperating...\n";
			typename TImageType::SizeType outputSize;

			for (unsigned int i=0; i < dimensions; i++) {
				outputSize[i]    = (inputSpacing[i] / outputSpacing[i]) * inputSize[i];
			}

			typedef itk::IdentityTransform<double, TImageType::ImageDimension> TransformType;

			typedef itk::ResampleImageFilter<TImageType, TImageType> ResampleImageFilterType;
			typename ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();
			resample->SetInput(input);
			resample->SetOutputParametersFromImage(input);
			resample->SetOutputStartIndex(input->GetLargestPossibleRegion().GetIndex());
			resample->SetOutputSpacing(outputSpacing);
			resample->SetOutputDirection(input->GetDirection());
			resample->SetTransform(TransformType::New());

			if (forceCertainSize) {
				resample->SetSize(forcedSize);
			}
			else {
				resample->SetSize(outputSize);
			}

			if (forceCertainSize) {
				resample->SetOutputOrigin(forcedOrigin);
			}
			else {
				resample->SetOutputOrigin(input->GetOrigin());
			}

			if (isLabelsImage) {
				typedef itk::NearestNeighborInterpolateImageFunction<TImageType, double> Interpolator;
				resample->SetInterpolator(Interpolator::New());
			}

			resample->UpdateLargestPossibleRegion();

			// Print output information
			std::cout << "\tOutput image size:    " << ((forceCertainSize)? forcedSize : outputSize) << "\n";
			std::cout << "\tOutput image spacing: " << outputSpacing << "\n";

			return resample->GetOutput();			
		}

		/**Change the image spacing to be equal in all directions. Picks the the smallest value
		*  from the original spacing for everything. Example [0.2,0.5,2] -> [0.2,0.2,0.2]*/
		template<class TImageType>
		static typename TImageType::Pointer 
		normalizeImageSpacing(typename TImageType::Pointer input, 
                              bool isLabelsImage = false)
		{
			typename TImageType::SizeType inputSize = input->GetLargestPossibleRegion().GetSize();
			typename TImageType::SpacingType inputSpacing = input->GetSpacing();
			const unsigned int dimensions = inputSize.GetSizeDimension();

			// Find the minimum value from input spacing
			float minSpacing = inputSpacing[0];
			for (unsigned int i=1; i < dimensions; i++) {
				if (inputSpacing[i] < minSpacing) {
					minSpacing = inputSpacing[i];
				}
			}

			// Set it as the output spacing
			typename TImageType::SpacingType outputSpacing;

			for (unsigned int i=0; i < dimensions; i++) {
				outputSpacing[i] = minSpacing;
			}

			return changeImageSpacing<TImageType>(input, outputSpacing, isLabelsImage);
		}

		template <class TImageType>
		static typename TImageType::Pointer 
		resizeImage(
			typename TImageType::Pointer input, 
			typename TImageType::SizeType outputSize,
			bool isLabelsImage = false) 
		{
			typename TImageType::SizeType inputSize = input->GetLargestPossibleRegion().GetSize();
			const unsigned int dimensions = inputSize.GetSizeDimension();

			// Print input information
			std::cout << "Input image size:         " << inputSize << "\n";

			// Resize
			typename TImageType::SpacingType outputSpacing;

			//outputSpacing = input->GetSpacing();
			for (unsigned int i=0; i < dimensions; i++) {
				outputSpacing[i] = input->GetSpacing()[i] * 
					(static_cast<double>(inputSize[i]) / static_cast<double>(outputSize[i]));
			}


			std::cout << "Resizing...\n";
			typedef itk::IdentityTransform<double, TImageType::ImageDimension> TransformType;

			typedef itk::ResampleImageFilter<TImageType, TImageType> ResampleImageFilterType;
			typename ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();
			resample->SetInput(input);
			resample->SetOutputParametersFromImage(input);
			resample->SetSize(outputSize);
			resample->SetOutputSpacing(outputSpacing);
			resample->SetOutputOrigin(input->GetOrigin());
			resample->SetOutputStartIndex(input->GetLargestPossibleRegion().GetIndex());
			resample->SetOutputDirection(input->GetDirection());
			resample->SetTransform(TransformType::New());

			if (isLabelsImage) {
				typedef itk::NearestNeighborInterpolateImageFunction<TImageType, double> Interpolator;
				resample->SetInterpolator(Interpolator::New());
			}

			resample->UpdateLargestPossibleRegion();
			return resample->GetOutput();
		}

		/**Resize an image to have about requestedNumberOfPixels pixels*/
		template <class TImageType>
		static typename TImageType::Pointer 
		resizeImageMaximumPixelNumber(
			typename TImageType::Pointer input, 
			unsigned int requestedNumberOfPixels,
			bool isLabelsImage = false) 
		{
			typename TImageType::SizeType inputSize = input->GetLargestPossibleRegion().GetSize();
			const unsigned int dimensions = inputSize.GetSizeDimension();

			unsigned int inputPixelTotal = 1;
			for (unsigned int i=0; i < dimensions; i++) {
				inputPixelTotal *= inputSize[i];
			}

			if (inputPixelTotal <= requestedNumberOfPixels) { return input; }

			// Print input information
			std::cout << "\tInput image size:         " << inputSize << "\n";
			std::cout << "\tInput total # of pixels:  " << inputPixelTotal << "\n";
			std::cout << "\tRequested # of pixels:    " << requestedNumberOfPixels << "\n---\n";

			double ratio = 1.0 * requestedNumberOfPixels / inputPixelTotal;
			double ratioPerDimension = std::pow(ratio, 1.0 / dimensions);

			typename TImageType::SizeType outputSize;
			
			for (unsigned int i=0; i < dimensions; i++) {
				outputSize[i] = std::lround(inputSize[i] * ratioPerDimension);
				if (outputSize[i] == 0) { outputSize[i] = 1; }
			}

			unsigned int outputPixelTotal = 1;
			for (unsigned int i=0; i < dimensions; i++) {
			  outputPixelTotal *= outputSize[i];
			}
			std::cout << "\tOutput image size:        " << outputSize << "\n";
			std::cout << "\tOutput total # of pixels: " << outputPixelTotal << "\n---\n";
			std::cout << "\tProposed ratio per dimensions: " << ratioPerDimension << "\n---\n"; 

			// Resize
			return resizeImage<TImageType>(input, outputSize, isLabelsImage);
		}

		template<class TImageType, class UImageType = TImageType>
		static typename UImageType::Pointer
		curvatureAnisotropicDiffusionImageFilter(typename TImageType::Pointer input, 
		                                         bool useImageSpacing = true,
			                                     bool rescaleAtTheEnd = false, 
			                                     float rescaleMax = 255)
		{
			const unsigned int numberOfIterations = 5;
			const double       timeStep = ((TImageType::ImageDimension == 2) ? 0.12 : 0.05);
			const double       conductance = 3;

			using FilterType = itk::CurvatureAnisotropicDiffusionImageFilter< TImageType, UImageType >;
			typename FilterType::Pointer filter = FilterType::New();
			filter->SetInput( input );
			filter->SetNumberOfIterations( numberOfIterations );
			filter->SetTimeStep( timeStep );
			filter->SetConductanceParameter( conductance );
			
			if (useImageSpacing) 
			{ 
				filter->UseImageSpacingOn(); 
			}

			filter->Update();

			auto output = filter->GetOutput();

			if (rescaleAtTheEnd)
			{
			  using RescaleFilterType = itk::RescaleIntensityImageFilter<UImageType, UImageType>;
			  typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
			  rescaler->SetOutputMinimum(   0 );
			  rescaler->SetOutputMaximum( rescaleMax );
			  rescaler->SetInput( output );
			  rescaler->Update();
			  output = rescaler->GetOutput();
			}

			return output;
		}

		template<class TImageType, class UImageType = TImageType>
		static typename UImageType::Pointer
		gradientAnisotropicDiffusionImageFilter(typename TImageType::Pointer input, 
			                                    bool rescaleAtTheEnd = false,
			                                    float rescaleMax = 255)
		{
			const unsigned int numberOfIterations = 15;
			const double       timeStep = ((TImageType::ImageDimension == 2) ? 0.12 : 0.05);
			const double       conductance = 3;

			using FilterType = itk::GradientAnisotropicDiffusionImageFilter< TImageType, UImageType >;
			typename FilterType::Pointer filter = FilterType::New();
			filter->SetInput( input );
			filter->SetNumberOfIterations( numberOfIterations );
			filter->SetTimeStep( timeStep );
			filter->SetConductanceParameter( conductance );
			
			filter->Update();

			typename UImageType::Pointer output = filter->GetOutput();

			if (rescaleAtTheEnd)
			{
			  using RescaleFilterType = itk::RescaleIntensityImageFilter<UImageType, UImageType>;
			  typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
			  rescaler->SetOutputMinimum(   0 );
			  rescaler->SetOutputMaximum( rescaleMax );
			  rescaler->SetInput( output );
			  rescaler->Update();
			  output = rescaler->GetOutput();
			}

			return output;
		}

		template<class TImageType, class UImageType = TImageType>
		static typename UImageType::Pointer
		BilateralImageFilter(typename TImageType::Pointer input, 
			                 bool rescaleAtTheEnd = false, float rescaleMax = 255)
		{
			using FilterType = itk::BilateralImageFilter< TImageType, UImageType >;
			typename FilterType::Pointer filter = FilterType::New();
			filter->SetInput( input );
			filter->SetDomainSigma( 3.0 );
			filter->SetRangeSigma( 30.0 );

			filter->Update();

			typename UImageType::Pointer output = filter->GetOutput();

			if (rescaleAtTheEnd)
			{
				std::cout << "Rescaling...\n";
				using RescaleFilterType = itk::RescaleIntensityImageFilter<UImageType, UImageType>;
				typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
				rescaler->SetOutputMinimum(   0 );
				rescaler->SetOutputMaximum( rescaleMax );
				rescaler->SetInput( output );
				rescaler->Update();
				output = rescaler->GetOutput();
			}

			return output;
		}


		/** Increase the radius of points in a labels image. 
		*   Obviously useful for hand drawings, rather than an image where labels cover everything.
		*   Think of it like increasing the "marker size" after drawing. 
		@param inputLabels the input labels image
		@param boldRadius an array with size equal to the image dimensions. It has the desired bolding.
		@return the bolded labels image 
		*/
		template <class TImageType>
		static typename TImageType::Pointer 
		boldLabelsImage(const typename TImageType::Pointer& inputLabels, 
		                const int (&boldRadius)[TImageType::ImageDimension])
		{
			typedef itk::ImageRegionIterator<TImageType>        Iterator;
			typedef itk::ShapedNeighborhoodIterator<TImageType> NIterator;
			typedef itk::BinaryBallStructuringElement<typename TImageType::PixelType, 
				TImageType::ImageDimension> StructuringElementType;

			typename StructuringElementType::RadiusType elementRadius;
			for (unsigned int i = 0; i < TImageType::ImageDimension; i++) { 
				elementRadius[i] = boldRadius[i]; 
			}

			// Structuring element is used to give shape to the neighborhood iterator (here a ellipsoid)
			StructuringElementType structuringElement;
			structuringElement.SetRadius(elementRadius);
			structuringElement.CreateStructuringElement();

			// Initialize the output image
			typename TImageType::Pointer outputLabels = initializeOutputImageBasedOn<TImageType>(inputLabels);

			// Neighborhood iterators for input and output
			// NIterator snIterInput(structuringElement.GetRadius(), inputLabels, inputLabels->GetRequestedRegion());
			// snIterInput.CreateActiveListFromNeighborhood(structuringElement);
			// snIterInput.NeedToUseBoundaryConditionOff();
			Iterator iterInput(inputLabels, inputLabels->GetRequestedRegion());
			
			NIterator snIterOutput(structuringElement.GetRadius(), outputLabels, outputLabels->GetRequestedRegion());
			snIterOutput.CreateActiveListFromNeighborhood(structuringElement);
			snIterOutput.NeedToUseBoundaryConditionOff();

			// Helper variables for the loop
			typename TImageType::PixelType val;
			bool statusIgnore;
			// typename NIterator::Iterator nIterInput;
			typename NIterator::Iterator nIterOutput;

			// Iterate throughout the input/output images
			for (/*snIterInput*/iterInput.GoToBegin(), snIterOutput.GoToBegin(); 
				 !iterInput/*snIterInput*/.IsAtEnd();
				 ++iterInput/*snIterInput*/, ++snIterOutput)
			{
				val = iterInput/*nIterInput*/.Get();
				if (val == 0) { continue; } // 0 means no label
				
				// Iterate through the neighborhood
				for (nIterOutput = snIterOutput.Begin(); !nIterOutput.IsAtEnd(); ++nIterOutput)
				{
					nIterOutput.Set(val);
				}
			}

			return outputLabels;
		}

		template <class TImageType>
		static typename TImageType::Pointer 
		resampleNormalImage(const typename TImageType::Pointer& image, 
		                const typename TImageType::SizeType& desiredSize, 
		                const typename TImageType::SpacingType& desiredSpacing)
		{
			typename TImageType::SizeType inputSize = image->GetLargestPossibleRegion().GetSize();
			typename TImageType::SpacingType inputSpacing = image->GetSpacing();

			using TransformType = itk::IdentityTransform<double, TImageType::ImageDimension>;
			using ResampleFilterType = itk::ResampleImageFilter<TImageType, TImageType>;
			
			// Instantiate the transform and specify it should be the identity transform.
			typename TransformType::Pointer transform = TransformType::New();
			transform->SetIdentity();
			
			// Initiate the resampler
			typename ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
			resampleFilter->SetTransform(transform);
			resampleFilter->SetDefaultPixelValue(0);
			resampleFilter->SetOutputParametersFromImage(image);
			resampleFilter->SetOutputStartIndex(image->GetLargestPossibleRegion().GetIndex());
			resampleFilter->SetOutputOrigin(image->GetOrigin());
			resampleFilter->SetOutputDirection(image->GetDirection());

			resampleFilter->SetSize(desiredSize);
			resampleFilter->SetOutputSpacing(desiredSpacing);

			// Execute the filter
			resampleFilter->SetInput(image);
			resampleFilter->UpdateLargestPossibleRegion();

			return resampleFilter->GetOutput();
		}

		template <class TImageType>
		static typename TImageType::Pointer 
		resampleLabelsImage(const typename TImageType::Pointer& inputLabels, 
		                const typename TImageType::SizeType& desiredSize, 
		                const typename TImageType::SpacingType& desiredSpacing, 
		                bool verbose = false)
		{
			typename TImageType::SizeType inputSize = inputLabels->GetLargestPossibleRegion().GetSize();
			typename TImageType::SpacingType inputSpacing = inputLabels->GetSpacing();

			using TransformType = itk::IdentityTransform<double, TImageType::ImageDimension>;
			// using GaussianInterpolatorType = 
			// 	itk::LabelImageGaussianInterpolateImageFunction<TImageType, double>;
			using NearestNeighborInterpolatorType = 
				itk::NearestNeighborInterpolateImageFunction<TImageType, double>;
			using ResampleFilterType = itk::ResampleImageFilter<TImageType, TImageType>;
			
			// Instantiate the transform and specify it should be the identity transform.
			typename TransformType::Pointer transform = TransformType::New();
			transform->SetIdentity();
			
			// Initiate the interpolator
			typename NearestNeighborInterpolatorType::Pointer nearestNeighborInterpolator =
				NearestNeighborInterpolatorType::New();

			// Initiate the resampler
			typename ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
			resampleFilter->SetTransform(transform);
			resampleFilter->SetInterpolator(nearestNeighborInterpolator);
			resampleFilter->SetDefaultPixelValue(0);
			resampleFilter->SetOutputParametersFromImage(inputLabels);
			resampleFilter->SetOutputStartIndex(inputLabels->GetLargestPossibleRegion().GetIndex());
			resampleFilter->SetOutputOrigin(inputLabels->GetOrigin());
			resampleFilter->SetOutputDirection(inputLabels->GetDirection());

			resampleFilter->SetSize(desiredSize);
			resampleFilter->SetOutputSpacing(desiredSpacing);
			
			// The input labels might need to be bolded before rescaling
			// Find the bold radius for each dimension
			int boldRadius[TImageType::ImageDimension];
			bool shouldBold = false;
			if (verbose) { std::cout << "\tboldRadius: ["; }
			for (unsigned int i=0; i<TImageType::ImageDimension; i++) {
				double ratio = 1.0 * inputSize[i] / desiredSize[i];
				boldRadius[i] = std::ceil(ratio / 2);
				if (boldRadius[i] >= 1) { shouldBold = true; }
				if (boldRadius[i] < 1) { boldRadius[i] = 0; }
				if (verbose) { std::cout << boldRadius[i] << ((i+1==TImageType::ImageDimension)?"":", "); }
			} 
			if (verbose) { std::cout << "]\n"; }

			// Bolding to more that 1 is probably pointless, though.
			// It's better to keep the above code though, as it might be useful in the future.
			// So this is a hack to use a maximum of 1 if there is no zeros, or 2 otherwise
			{
				if (shouldBold)
				{
					// Find if there are zeros
					bool areThereZeros = false;
					bool areThereValuesBiggerThanOne = false;
					for (unsigned int i=0; i<TImageType::ImageDimension; i++) {
						if (boldRadius[i] == 0) { 
							areThereZeros = true;
						}
						if (boldRadius[i] > 1) {
							areThereValuesBiggerThanOne = true;
						}
					}

					int maxRadius = ((areThereValuesBiggerThanOne) ? 2 : 1);
					maxRadius = ((areThereZeros) ? maxRadius : 1);
				
					if (verbose) { std::cout << "\tActually used boldRadius: ["; }
					for (unsigned int i=0; i<TImageType::ImageDimension; i++) {
						if (boldRadius[i] > 1) {
							boldRadius[i] = maxRadius;
						}
						if (verbose) { std::cout << boldRadius[i] << ((i+1==TImageType::ImageDimension)?"":", "); }
					}
					if (verbose) { std::cout << "]\n"; }
				}
			}

			// Bold the labels image if necessary
			if (shouldBold)
			{
				if (verbose) { std::cout << "\tBolding...\n"; }
				resampleFilter->SetInput(boldLabelsImage<TImageType>(inputLabels, boldRadius));
			}
			else {
				resampleFilter->SetInput(inputLabels);
			}

			// Execute the filter
			std::cout << "\tUpdating Nearest Neighbor Interpolator Rescaling...\n";
			resampleFilter->UpdateLargestPossibleRegion();

			return resampleFilter->GetOutput();
		}

	}
}

#endif // !H_CBICA_ITK_UTIL_GTS