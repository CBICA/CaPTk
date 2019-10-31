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
		std::unordered_map<int, int> labelsCountMap, int maxSamples);

}

#endif // !H_CBICA_SVM_GTS
