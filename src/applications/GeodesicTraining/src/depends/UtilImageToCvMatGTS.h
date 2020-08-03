#ifndef H_CBICA_UTIL_IMAGE_TO_CV_MAT_GTS
#define H_CBICA_UTIL_IMAGE_TO_CV_MAT_GTS

#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkExceptionObject.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"

#include <opencv2/ml.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <vector>
#include <map>
#include <unordered_map>

namespace GeodesicTrainingSegmentation 
{
	namespace ParserGTS
	{
		typedef int LabelsPixelType;

		typedef struct Result {
			cv::Mat trainingMat;
			cv::Mat labelsMat;
			cv::Mat weightsMat;
			cv::Mat testingMat;
			cv::Mat skipZerosMat;
			std::set<LabelsPixelType> differentLabels;
		} Result;

		template <class TImageType>
		int getNumberOfPixels(typename TImageType::Pointer image)
		{
			// Get number of pixels of the images
			typename TImageType::RegionType region = image->GetLargestPossibleRegion();
			typename TImageType::SizeType   imageSize = region.GetSize();
			int numberOfPixels = 1;

			for (int i = 0; i < TImageType::ImageDimension; i++) {
				// Multiply with image size in each dimension
				numberOfPixels *= imageSize[i];
			}

			return numberOfPixels;
		}

		template <typename TPixelType, unsigned int TDimensions>
		int getNumberOfNotZerosFromLabels(typename itk::Image<LabelsPixelType, TDimensions>::Pointer labels,
			typename itk::Image<TPixelType, TDimensions>::Pointer oneOfInput, bool considerZeros)
		{
			typedef itk::Image<TPixelType, TDimensions>                        TImageType;
			typedef typename itk::Image<TPixelType, TDimensions>::Pointer      TImagePointer;
			typedef itk::Image<LabelsPixelType, TDimensions>                   LabelsImageType;
			typedef typename LabelsImageType::Pointer                          LabelsImagePointer;

			itk::ImageRegionIteratorWithIndex<LabelsImageType> iter_l(labels, labels->GetRequestedRegion());
			itk::ImageRegionIteratorWithIndex<TImageType> iter_i(oneOfInput, oneOfInput->GetRequestedRegion());
			int labeled = 0;

			for (iter_l.GoToBegin(), iter_i.GoToBegin(); !iter_l.IsAtEnd(); ++iter_l, ++iter_i) {
				if (iter_l.Get() != 0) {
					labeled++;

					if (!considerZeros && iter_i.Get() == 0) {
						//std::cout << "One of the labels is at a voxel out of bounds and will not be counted."
						//	<< " Consider redrawing the mask\n";
						labeled--;
					}
				}
			}

			return labeled;
		}

		/**
		Finds how many samples are there for each label used in an itk::Image
		The PixelType of the input image should always be int
		@param labels the input labels image
		@return unordered_map with the different labels as keys and the counts for each label as values
		*/
		template <class TImageType>
		static std::unordered_map<LabelsPixelType, int> CountsOfEachLabel(const typename TImageType::Pointer labels)
		{
			std::unordered_map<LabelsPixelType, int> labelsCountMap;

			itk::ImageRegionIteratorWithIndex<TImageType> iter_l(labels, labels->GetRequestedRegion());

			LabelsPixelType val;

			for (iter_l.GoToBegin(); !iter_l.IsAtEnd(); ++iter_l) {
				val = iter_l.Get();

				if (val != 0) {
					if (labelsCountMap.find(val) == labelsCountMap.end()) {
						labelsCountMap[val] = 1;
					}
					else {
						labelsCountMap[val] += 1;
					}
				}
			}

			return labelsCountMap;
		}

		template <typename TPixelType, unsigned int TDimensions>
		std::shared_ptr<Result> Parse(const std::vector< typename itk::Image<TPixelType, TDimensions>::Pointer > &input_images,
			const typename itk::Image<LabelsPixelType, TDimensions>::Pointer &input_labels, bool considerZeros = false)
		{
			typedef itk::Image<TPixelType, TDimensions>                        InputImageType;
			typedef typename itk::Image<TPixelType, TDimensions>::Pointer      InputImagePointer;
			typedef itk::Image<LabelsPixelType, TDimensions>                   LabelsImageType;
			typedef typename LabelsImageType::Pointer                          LabelsImagePointer;

			std::shared_ptr<Result> res(new Result());

			// Initialize the labels iterator
			itk::ImageRegionIteratorWithIndex<LabelsImageType> iter_labels(input_labels, input_labels->GetRequestedRegion());

			// Get number of pixels of the images
			int numberOfPixels = getNumberOfPixels<InputImageType>(input_images[0]);

			// Number of voxels in the labels image that are not zero (and are usable labels)
			int labelsNumber = getNumberOfNotZerosFromLabels<TPixelType, TDimensions>(input_labels, input_images[0], considerZeros);

			// Initialize result
			res->testingMat   = cv::Mat::zeros(numberOfPixels, input_images.size(), CV_32F);
			res->skipZerosMat = cv::Mat::zeros(numberOfPixels, input_images.size(), CV_32F);
			res->trainingMat  = cv::Mat::zeros(labelsNumber, input_images.size(), CV_32F);
			res->labelsMat    = cv::Mat::zeros(labelsNumber, 1, CV_32S);

			// If there is a label for this voxel (use 1,2,... for labels, not zero)
			cv::Mat features = cv::Mat::zeros(1, input_images.size(), CV_32F);
			std::unordered_map< LabelsPixelType, int > weightSums;
			unsigned long allWeightsSum = 0;

			float val;
			LabelsPixelType label;
			bool isInputInt = (typeid(TPixelType).name() == typeid(int).name()) ? true : false;
			int inputImagesSize = input_images.size();
			bool skipPixel;
			int row_i = 0, row_testing_i = 0;
			size_t fi;

			for (iter_labels.GoToBegin(); !iter_labels.IsAtEnd(); ++iter_labels)
			{
				skipPixel = false;

				// Set the values of the feature vector
				for (fi = 0; fi < inputImagesSize; fi++) {
					val = input_images[fi]->GetPixel(iter_labels.GetIndex());

					// Only the first image is used to skip zeros
					if ((fi == 0) && (val == 0) && (!considerZeros)) {
						skipPixel = true;
						break;
					}

					features.ptr< TPixelType >(0)[fi] = (isInputInt) ? std::lround(val) : val;
				}

				// Populate matrices
				if (!skipPixel || considerZeros) {
					features.copyTo(res->testingMat.row(row_i++));

					label = iter_labels.Get();
					if (label != 0) {
						// Set training and labels data

						res->labelsMat.ptr< LabelsPixelType >(row_testing_i)[0] = label;
						res->differentLabels.insert(label);

						features.copyTo(res->trainingMat.row(row_testing_i++));

						// For weights calculation
						if (weightSums.find(label) == weightSums.end()) {
							weightSums[label] = 1;
						}
						else {
							weightSums[label] += 1;
						}
						allWeightsSum++;
					}
				}
				else {
					row_i++;
				}
			}

			// Extra for weights
			res->weightsMat = cv::Mat::zeros(res->differentLabels.size(), 1, CV_32F);

			int wi = 0;
			for (LabelsPixelType label : res->differentLabels)
			{
				if (allWeightsSum != 0) {
					res->weightsMat.at< float >(wi++, 0) = static_cast<float>(allWeightsSum - weightSums[label]) / allWeightsSum;
				}
			}

			res->testingMat.copyTo(res->skipZerosMat);

			return res;
		}

		template <typename TPixelType, unsigned int TDimensions>
		std::shared_ptr<Result> NormalizedParse(const std::vector< typename itk::Image<TPixelType, TDimensions>::Pointer > &input_images,
			const typename itk::Image<LabelsPixelType, TDimensions>::Pointer &input_labels, bool considerZeros = false)
		{
			typedef itk::Image<TPixelType, TDimensions>                        InputImageType;
			typedef typename itk::Image<TPixelType, TDimensions>::Pointer      InputImagePointer;
			typedef itk::Image<LabelsPixelType, TDimensions>                   LabelsImageType;
			typedef typename LabelsImageType::Pointer                          LabelsImagePointer;

			std::shared_ptr<Result> res(new Result());

			// Initialize the labels iterator
			itk::ImageRegionIteratorWithIndex<LabelsImageType> iter_labels(input_labels, input_labels->GetRequestedRegion());

			// Get number of pixels of the images
			int numberOfPixels = getNumberOfPixels<InputImageType>(input_images[0]);

			// Number of voxels in the labels image that are not zero (and are usable labels)
			int labelsNumber = getNumberOfNotZerosFromLabels<TPixelType, TDimensions>(input_labels, input_images[0], considerZeros);

			// Initialize result
			res->testingMat   = cv::Mat::zeros(numberOfPixels, input_images.size(), CV_32F);
			res->skipZerosMat = cv::Mat::zeros(numberOfPixels, input_images.size(), CV_32F);
			res->trainingMat  = cv::Mat::zeros(labelsNumber, input_images.size(), CV_32F);
			res->labelsMat    = cv::Mat::zeros(labelsNumber, 1, CV_32S);

			// If there is a label for this voxel (use 1,2,... for labels, not zero)
			cv::Mat features = cv::Mat::zeros(1, input_images.size(), CV_32F);
			std::unordered_map< LabelsPixelType, int > weightSums;
			unsigned long allWeightsSum = 0;

			float val;
			LabelsPixelType label;
			bool isInputInt = (typeid(TPixelType).name() == typeid(int).name()) ? true : false;
			int inputImagesSize = input_images.size();
			bool skipPixel;
			int row_i = 0, row_testing_i = 0;
			size_t fi;

			for (iter_labels.GoToBegin(); !iter_labels.IsAtEnd(); ++iter_labels)
			{
				skipPixel = false;

				// Set the values of the feature vector
				for (fi = 0; fi < inputImagesSize; fi++) {
					val = input_images[fi]->GetPixel(iter_labels.GetIndex());

					// Only the first image is used to skip zeros
					if ((fi == 0) && (val == 0) && (!considerZeros)) {
						skipPixel = true;
						break;
					}

					features.ptr< TPixelType >(0)[fi] = (isInputInt) ? std::lround(val) : val;
				}

				// Populate matrices
				if (!skipPixel || considerZeros) {
					features.copyTo(res->testingMat.row(row_i++));

					label = iter_labels.Get();
					if (label != 0) {
						// Set training and labels data

						res->labelsMat.ptr< LabelsPixelType >(row_testing_i)[0] = label;
						res->differentLabels.insert(label);

						features.copyTo(res->trainingMat.row(row_testing_i++));

						// For weights calculation
						if (weightSums.find(label) == weightSums.end()) {
							weightSums[label] = 1;
						}
						else {
							weightSums[label] += 1;
						}
						allWeightsSum++;
					}
				}
				else {
					row_i++;
				}
			}

			// Extra for weights
			res->weightsMat = cv::Mat::zeros(res->differentLabels.size(), 1, CV_32F);

			int wi = 0;
			for (LabelsPixelType label : res->differentLabels)
			{
				if (allWeightsSum != 0) {
					res->weightsMat.at< float >(wi++, 0) = static_cast<float>(allWeightsSum - weightSums[label]) / allWeightsSum;
				}
			}

			res->testingMat.copyTo(res->skipZerosMat);

			// Normalization

			// Collect only the relevant samples
			cv::Mat onlyRelevantSamples;
			if (considerZeros) {
				res->trainingMat.copyTo(onlyRelevantSamples);
			}
			else {
				for (int row_i = 0; row_i < res->trainingMat.rows; row_i++)
				{
					if (!(res->trainingMat.at<float>(row_i, 0) == 0)) {
						onlyRelevantSamples.push_back(res->trainingMat.row(row_i));
					}
				}
			}

			// Normalize by subtracting mean and dividing by twice the standard deviation
			cv::Scalar meanValue, stdValue;
			double mean, std;

			for (int i_col = 0; i_col < onlyRelevantSamples.cols; i_col++)
			{
				cv::meanStdDev(onlyRelevantSamples.col(i_col), meanValue, stdValue);

				mean = meanValue[0];
				std = (stdValue[0] == 0) ? 0.001 : stdValue[0];

				res->trainingMat.col(i_col).convertTo(res->trainingMat.col(i_col), CV_32F, 1, -1 * mean);
				res->trainingMat.col(i_col).convertTo(res->trainingMat.col(i_col), CV_32F, 1 / (2 * std), 0);

				res->testingMat.col(i_col).convertTo(res->testingMat.col(i_col), CV_32F, 1, -1 * mean);
				res->testingMat.col(i_col).convertTo(res->testingMat.col(i_col), CV_32F, 1 / (2 * std), 0);
			}

			return res;
		}

		void ScaleSomeOfTheColumns(cv::Mat& mat, int colStart, int colEnd, double ratio);
	
	}

}

#endif // !H_CBICA_UTIL_IMAGE_TO_CV_MAT_GTS

