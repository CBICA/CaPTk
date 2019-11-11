#include "SvmSuiteOperations.h"

#include "ConvertionsYAML.h"

void SvmSuite::NormalizeInput(cv::Mat& mat)
{
	//cv::Mat meanValue = cv::Mat::zeros(1,1,CV_32F), stdValue = cv::Mat::zeros(1, 1, CV_32F);
	cv::Scalar meanValue, stdValue;
	double mean, std;

	for (int i_col = 0; i_col < mat.cols; i_col++)
	{
		cv::meanStdDev(mat.col(i_col), meanValue, stdValue);

		//mean = meanValue.at<float>(0,0);
		//std  = (stdValue.at<float>(0, 0) == 0) ? 0.001 : stdValue.at<float>(0, 0);
		mean = meanValue[0];
		std  = (stdValue[0] == 0) ? 0.001 : stdValue[0];

		mat.col(i_col).convertTo(mat.col(i_col), CV_32F, 1, -1 * mean);
		mat.col(i_col).convertTo(mat.col(i_col), CV_32F, 1 / (2 * std), 0);
	}
}

std::vector< SvmSuite::SvmDescription > SvmSuite::getSvmDescriptionsFromConfig(std::string configPath)
{
	std::vector< SvmSuite::SvmDescription > svm_descriptions;
	YAML::Node root = YAML::LoadFile(configPath);

	// Check if the root node of the document is named as stated in YAML_ROOT_NODE
	if (root[YAML_ROOT_NODE])
	{
		// Iterate through the different SVMs in the config
		for (YAML::const_iterator it_svms = root[YAML_ROOT_NODE].begin(); it_svms != root[YAML_ROOT_NODE].end(); ++it_svms)
		{
			SvmSuite::SvmDescription svm;

			// Iterate through the different options for each SVM in the config
			for (YAML::const_iterator it_single_svm = it_svms->begin(); it_single_svm != it_svms->end(); ++it_single_svm)
			{
				std::string configOption = it_single_svm->first.as<std::string>();

				if (configOption == YAML_KERNEL_TYPE) {
					std::string val = it_single_svm->second.as<std::string>();
					svm.SetKernelType(val);
				}
				else if (configOption == YAML_TYPE) {
					std::string val = it_single_svm->second.as<std::string>();
					svm.SetType(val);
				}
				else if (configOption == YAML_KFOLD) {
					int val = it_single_svm->second.as<int>();
					svm.SetKfold(val);
				}
				else if (configOption == YAML_NEIGHBORHOOD_RADIUS) {
					int val = it_single_svm->second.as<int>();
					svm.SetNeighborhoodRadius(val);
				}
				else if (configOption == YAML_IMPORTANCE) {
					double val = it_single_svm->second.as<double>();
					svm.SetImportance(val);
				}
				else if (configOption == YAML_CONSIDER_WEIGHTS) {
					bool val = it_single_svm->second.as<bool>();
					svm.SetConsiderWeights(val);
				}
				else if (configOption == YAML_MODEL_PATH) {
					std::string val = it_single_svm->second.as<std::string>();
					svm.SetModelPath(val);
				}
				else if (configOption == YAML_TERM_CRITERIA) {
					// Term criteria must have the following variables
					std::string criteria_type = "";
					std::string max_count = "";
					std::string epsilon = "";

					for (YAML::const_iterator it_term = it_single_svm->second.begin(); it_term != it_single_svm->second.end(); ++it_term)
					{
						std::string term_param_name = it_term->first.as<std::string>();
						std::string term_param_val = it_term->second.as<std::string>();

						if (term_param_name == YAML_TERM_CRITERIA_TYPE) {
							criteria_type = term_param_val;
						}
						else if (term_param_name == YAML_TERM_CRITERIA_MAX) {
							max_count = term_param_val;
						}
						else if (term_param_name == YAML_TERM_CRITERIA_EPS) {
							epsilon = term_param_val;
						}
						else {
							//invalid
						}
					}

					svm.SetTermCriteria(criteria_type, max_count, epsilon);
				}
				else if ((configOption == YAML_C) || (configOption == YAML_GAMMA) || (configOption == YAML_P) ||
					(configOption == YAML_NU) || (configOption == YAML_COEF) || (configOption == YAML_DEGREE))
				{
					if (it_single_svm->second.IsScalar())
					{
						std::string val = it_single_svm->second.as<std::string>();
						if (val == YAML_AUTO) {
							// Parameter was set to auto automatically
							svm.SetParameterRangeAuto(configOption);
						}
						else {
							// Specific value for the parameter was given
							svm.SetParameter(configOption, val);
						}
					}
					else
					{
						// Parameter range must have the following variables
						double minVal = 0.0;
						double maxVal = 0.0;
						double logStep = 0.0;

						for (YAML::const_iterator it_param = it_single_svm->second.begin(); it_param != it_single_svm->second.end(); ++it_param)
						{
							std::string range_param_name = it_param->first.as<std::string>();
							double range_param_val = it_param->second.as<double>();

							if (range_param_name == YAML_MIN_VAL) {
								minVal = range_param_val;
							}
							else if (range_param_name == YAML_MAX_VAL) {
								maxVal = range_param_val;
							}
							else if (range_param_name == YAML_LOG_STEP) {
								logStep = range_param_val;
							}
							else {
								//invalid
							}
						}

						svm.SetParameterRange(configOption, minVal, maxVal, logStep);
					}
				}
			}

			svm_descriptions.push_back(svm);
		}
	}
	else {
		//invalid
	}

	return svm_descriptions;
}

SvmSuite::SvmDescription SvmSuite::convertModelToSvmDescription(cv::Ptr<cv::ml::SVM> svm_model, int neighborhood_radius, double importance)
{
	SvmSuite::SvmDescription svm_description;
	svm_description.SetKernelType(static_cast<cv::ml::SVM::KernelTypes>(svm_model->getKernelType()));
	svm_description.SetType(static_cast<cv::ml::SVM::Types>(svm_model->getType()));

	if (svm_model->getC() != 1.0) {
		svm_description.SetC(svm_model->getC());
	}
	if (svm_model->getGamma() != 1.0) {
		svm_description.SetGamma(svm_model->getGamma());
	}
	if (svm_model->getP() != 0.0) {
		svm_description.SetP(svm_model->getP());
	}
	if (svm_model->getNu() != 0.0) {
		svm_description.SetNu(svm_model->getNu());
	}
	if (svm_model->getCoef0() != 0.0) {
		svm_description.SetCoef(svm_model->getCoef0());
	}
	if (svm_model->getDegree() != 0) {
		svm_description.SetDegree(svm_model->getDegree());
	}

	if (neighborhood_radius != 0) {
		svm_description.SetNeighborhoodRadius(neighborhood_radius);
	}
	if (importance != 1.0) {
		svm_description.SetImportance(importance);
	}

	svm_description.SetTermCriteria(svm_model->getTermCriteria().type, svm_model->getTermCriteria().maxCount, svm_model->getTermCriteria().epsilon);

	return svm_description;
}

SvmSuite::SvmDescription SvmSuite::convertModelToSvmDescription(std::string model_path, int neighborhood_radius, double importance)
{
	SvmSuite::SvmDescription svm_desc = convertModelToSvmDescription(cv::ml::SVM::load(model_path), neighborhood_radius, importance);
	svm_desc.SetModelPath(model_path);
	return svm_desc;
}

void SvmSuite::generateConfig(std::vector< SvmSuite::SvmDescription > &svm_descriptions, std::string outputFilePath, bool m_save_models)
{
	YAML::Node node;

	for (SvmSuite::SvmDescription& svm_desc : svm_descriptions) {
		YAML::Node toAdd = SvmSuiteConvertions::yamlConvertSvmDescriptionToNode(svm_desc);
		node[YAML_ROOT_NODE].push_back(toAdd);
	}

	if (m_save_models) {
		std::ofstream out(outputFilePath);
		out << "---\n" << node << "\n...";
		out.close();
	}
}

void SvmSuite::generateConfig(SvmSuite::SvmDescription &svm_description, std::string outputFilePath, bool m_save_models)
{
	std::vector< SvmSuite::SvmDescription > svm_descriptions;
	svm_descriptions.push_back(svm_description);

	generateConfig(svm_descriptions, outputFilePath);
}

void SvmSuite::generateConfig(cv::Ptr<cv::ml::SVM> model, int neighborhood_radius, double importance, std::string outputFilePath, bool m_save_models)
{
	std::vector< SvmSuite::SvmDescription > svm_descriptions;
	svm_descriptions.push_back(convertModelToSvmDescription(model, neighborhood_radius, importance));

	generateConfig(svm_descriptions, outputFilePath);
}

void SvmSuite::generateConfig(std::vector< cv::Ptr<cv::ml::SVM> > multiple_models, std::vector< int > neighborhood_radii,
	std::vector< double > importance_values, std::string outputFilePath, bool m_save_models)
{
	std::vector< SvmSuite::SvmDescription > svm_descriptions;

	for (int i = 0; i < multiple_models.size(); i++) {
		SvmSuite::SvmDescription svm_desc = convertModelToSvmDescription(multiple_models[i], neighborhood_radii[i], importance_values[i]);
		svm_desc.SetNeighborhoodRadius(neighborhood_radii[i]);
		svm_desc.SetImportance(importance_values[i]);
		svm_descriptions.push_back(svm_desc);
	}

	generateConfig(svm_descriptions, outputFilePath);
}

void SvmSuite::generateConfig(std::vector< std::string > multiple_models_paths, std::vector< int > neighborhood_radii,
	std::vector< double > importance_values, std::string outputFilePath, bool m_save_models)
{
	std::vector< SvmSuite::SvmDescription > svm_descriptions;

	for (int i = 0; i < multiple_models_paths.size(); i++) {
		SvmSuite::SvmDescription svm_desc = convertModelToSvmDescription(multiple_models_paths[i], neighborhood_radii[i], importance_values[i]);
		svm_desc.SetNeighborhoodRadius(neighborhood_radii[i]);
		svm_desc.SetImportance(importance_values[i]);
		svm_descriptions.push_back(svm_desc);
	}

	generateConfig(svm_descriptions, outputFilePath);
}
