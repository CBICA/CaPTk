#ifndef H_CBICA_SVM_SUITE_CONVERTIONS_OPEN_CV
#define H_CBICA_SVM_SUITE_CONVERTIONS_OPEN_CV

#include <opencv2/ml.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <string>

namespace SvmSuiteConvertions {

	// OpenCV to String (and vice versa) Conversions

	cv::ml::SVM::Types TypeFromString(std::string typeString);

	std::string StringFromType(cv::ml::SVM::Types type);

	cv::ml::SVM::KernelTypes KernelTypeFromString(std::string kernelTypeString);

	std::string StringFromKernelType(cv::ml::SVM::KernelTypes kernelType);

	cv::ml::SVM::ParamTypes ParamTypeFromString(std::string paramTypeString);

	std::string StringFromParamType(cv::ml::SVM::ParamTypes paramType);

	int TermCriteriaTypeFromString(std::string termCriteriaTypeString);

	std::string StringFromTermCriteriaType(int termCriteriaType);

}

#endif