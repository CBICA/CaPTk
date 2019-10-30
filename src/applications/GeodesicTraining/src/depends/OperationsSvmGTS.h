#ifndef H_CBICA_SVM_GTS
#define H_CBICA_SVM_GTS

#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkExceptionObject.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkInvertIntensityImageFilter.h"

#include <opencv2/ml.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "UtilItkGTS.h"
#include "UtilImageToCvMatGTS.h"
#include "UtilCvMatToImageGTS.h"
#include "SvmSuite.h"

#include <unordered_map>
#include <set>
#include <memory>

namespace GeodesicTrainingSegmentation
{
	template<typename TPixelType = float, unsigned int TDimensions = 3>
	class SvmManagerGTS : public SvmSuite::Manager
	{
	public:
		typedef int   LabelsPixelType;
		typedef float PseudoProbPixelType;

		typedef itk::Image< TPixelType, TDimensions >         InputImageType;
		typedef itk::Image< LabelsPixelType, TDimensions >     LabelsImageType;
		typedef itk::Image< PseudoProbPixelType, TDimensions > PseudoProbImageType;

		typedef typename InputImageType::Pointer              InputImagePointer;
		typedef typename LabelsImageType::Pointer             LabelsImagePointer;
		typedef typename PseudoProbImageType::Pointer         PseudoProbImagePointer;

		typedef struct ResultSvmGTS {
			PseudoProbImagePointer posImage;
			PseudoProbImagePointer negImage;
			LabelsPixelType posLabel = 0;
			LabelsPixelType negLabel = 0;
			LabelsImagePointer labelsImage;
		} ResultSvmGTS;

		std::shared_ptr<ResultSvmGTS> TrainAndTestGTS(std::shared_ptr<ParserGTS::Result> data, std::vector<InputImagePointer> images, 
                                                      LabelsImagePointer labels, bool predictFlags/*, cv::Mat sampleIdx = cv::Mat()*/)
		{
			this->SetTrainData(data->trainingMat, data->labelsMat, data->weightsMat/*, sampleIdx*/);

			this->Train();
			auto result = this->Test(data->testingMat, data->skipZerosMat, predictFlags, true);

			std::shared_ptr<ResultSvmGTS> resultGTS(new ResultSvmGTS());

			if (predictFlags) {
				// Pseudoprob
				resultGTS->posImage = ItkUtilGTS::initializeOutputImageBasedOn<LabelsImageType, PseudoProbImageType>(labels);
				resultGTS->negImage = ItkUtilGTS::initializeOutputImageBasedOn<LabelsImageType, PseudoProbImageType>(labels);
				
				std::cout << "Converting output matrices to images...";
				CvMatToImageGTS::FillImage<PseudoProbPixelType, TDimensions>(resultGTS->posImage, result->posMat);
				CvMatToImageGTS::FillImage<PseudoProbPixelType, TDimensions>(resultGTS->negImage, result->negMat);
				std::cout << "finished\n";

				resultGTS->posLabel = result->posLabel;
				resultGTS->negLabel = result->negLabel;
			}
			else {
				// Labels
				resultGTS->labelsImage = ItkUtilGTS::initializeOutputImageBasedOn<LabelsImageType>(labels);
				
				std::cout << "Converting output matrix to image...";
				CvMatToImageGTS::FillImage<LabelsPixelType, TDimensions>(resultGTS->labelsImage, result->labelsMat);
				std::cout << "finished\n";
			}

			return resultGTS;
		}
	};

	/** Balanced Subsampling */
	bool /*cv::Mat*/ CreateBalancedSubsample(std::shared_ptr<ParserGTS::Result>& data, 
		std::string& errorMessageIfApplicable,
		std::unordered_map<int, int> labelsCountMap, int maxSamples)
	{
		typedef int LabelsPixelType;

		// Find how many samples to keep for each label
		std::unordered_map<LabelsPixelType, int> maxForEachLabel;
		double totalLabels = 0;

		for (std::set<LabelsPixelType>::iterator it = data->differentLabels.begin(); it != data->differentLabels.end(); it++)
		{
			// Initialize
			maxForEachLabel[*it] = 0;
		}

		for (std::set<LabelsPixelType>::iterator it = data->differentLabels.begin(); it != data->differentLabels.end(); it++)
		{
			totalLabels += labelsCountMap[*it];
		}

		std::cout << "\n";
		for (std::set<LabelsPixelType>::iterator it = data->differentLabels.begin(); it != data->differentLabels.end(); it++)
		{
			maxForEachLabel[*it] = std::lround(maxSamples * labelsCountMap[*it] / totalLabels);
			std::cout << "MAX FOR " << *it << ": " << maxForEachLabel[*it] << "\n";
			if (maxForEachLabel[*it] < 3) {
				errorMessageIfApplicable = std::string("Cannot create a balanced subsample, ") +
					std::string("because one (or more) labels is outweighted by the others a lot. ") +
					std::string("Draw more for the labels that you haven't drawn a lot for, or remove ") +
					std::string("some samples from the ones that are large.");
				return false;
			}
		}

		// Shuffle training matrix and labels (with the same seed)
		std::vector <int> seeds;
		for (int cont = 0; cont < data->trainingMat.rows; cont++) {
			seeds.push_back(cont);
		}

		cv::randShuffle(seeds);

		cv::Mat randTrainingMat, randLabelsMat;
		for (int cont = 0; cont < data->trainingMat.rows; cont++) {
			randTrainingMat.push_back(data->trainingMat.row(seeds[cont]));
			randLabelsMat.push_back(data->labelsMat.row(seeds[cont]));
		}
		
		data->trainingMat = randTrainingMat;
		data->labelsMat = randLabelsMat;

		// // Create a mask for which samples to keep for training
		// cv::Mat sampleIdx = cv::Mat::zeros(data->trainingMat.rows, 1, CV_8U);

		std::unordered_map<LabelsPixelType, int> keptForEachLabel;

		for (std::set<LabelsPixelType>::iterator it = data->differentLabels.begin(); it != data->differentLabels.end(); it++)
		{
			// Initialize
			keptForEachLabel[*it] = 0;
		}

		cv::Mat trainingMatNew, labelsMatNew;

		for (int row_i = 0; row_i < data->trainingMat.rows; row_i++)
		{
			LabelsPixelType label = data->labelsMat.ptr<LabelsPixelType>(row_i)[0];
			
			if (keptForEachLabel[label] <= maxForEachLabel[label]) {
				//sampleIdx.ptr<uchar>(row_i)[0] = 1;
				keptForEachLabel[label] += 1;
				trainingMatNew.push_back(data->trainingMat.row(row_i).clone());
				labelsMatNew.push_back(data->labelsMat.row(row_i).clone());
			}
		}
		data->trainingMat = trainingMatNew;
		data->labelsMat   = labelsMatNew;

		return true;
		//return sampleIdx;
	}

}

#endif // !H_CBICA_SVM_GTS
