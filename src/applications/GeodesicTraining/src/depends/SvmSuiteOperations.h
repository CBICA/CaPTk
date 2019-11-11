#ifndef H_CBICA_SVM_SUITE_OPERATIONS
#define H_CBICA_SVM_SUITE_OPERATIONS

#include <opencv2/ml.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <iostream>
#include <fstream>
#include <vector>

#include "SvmSuiteDescription.h"

namespace SvmSuite
{
	/** Normalizes input for the SVMs
	@param mat the input matrix
	*/
	void NormalizeInput(cv::Mat& mat);

	std::vector< SvmSuite::SvmDescription > getSvmDescriptionsFromConfig(std::string configPath);

	SvmSuite::SvmDescription convertModelToSvmDescription(cv::Ptr<cv::ml::SVM> svm_model, 
		int neighborhood_radius, double importance);

	SvmSuite::SvmDescription convertModelToSvmDescription(std::string model_path, 
		int neighborhood_radius, double importance);

	void generateConfig(std::vector< SvmSuite::SvmDescription > &svm_descriptions, 
		std::string outputFilePath, bool m_save_models = false);

	void generateConfig(SvmSuite::SvmDescription &svm_description, 
		std::string outputFilePath, bool m_save_models = false);

	void generateConfig(cv::Ptr<cv::ml::SVM> model, int neighborhood_radius, 
		double importance, std::string outputFilePath, bool m_save_models = false);

	void generateConfig(std::vector< cv::Ptr<cv::ml::SVM> > multiple_models, std::vector< int > neighborhood_radii,
		std::vector< double > importance_values, std::string outputFilePath, bool m_save_models = false);

	void generateConfig(std::vector< std::string > multiple_models_paths, std::vector< int > neighborhood_radii,
		std::vector< double > importance_values, std::string outputFilePath, bool m_save_models = false);
}

#endif