#include "ConvertionsYAML.h"

#include "ConvertionsOpenCV.h"

YAML::Node SvmSuiteConvertions::yamlConvertSvmDescriptionToNode(SvmSuite::SvmDescription& svm_description)
{
	YAML::Node node;
	node[YAML_KERNEL_TYPE] = svm_description.GetKernelTypeAsString();
	node[YAML_TYPE] = svm_description.GetTypeAsString();

	if (svm_description.GetKfold() != 10) {
		node[YAML_KFOLD] = svm_description.GetKfold();
	}
	if (svm_description.GetNeighborhoodRadius() != 0) {
		node[YAML_NEIGHBORHOOD_RADIUS] = svm_description.GetNeighborhoodRadius();
	}
	if (svm_description.GetConsiderWeights()) {
		node[YAML_CONSIDER_WEIGHTS] = svm_description.GetConsiderWeights();
	}
	if (svm_description.GetImportance() != 1.0) {
		node[YAML_IMPORTANCE] = svm_description.GetImportance();
	}
	if (svm_description.GetModelPath() != "") {
		node[YAML_MODEL_PATH] = svm_description.GetModelPath();
	}

	for (const auto param : { cv::ml::SVM::ParamTypes::C,  cv::ml::SVM::ParamTypes::GAMMA, cv::ml::SVM::ParamTypes::P,
							  cv::ml::SVM::ParamTypes::NU, cv::ml::SVM::ParamTypes::COEF,  cv::ml::SVM::ParamTypes::DEGREE })
	{
		if (svm_description.isParameterSetToSpecificValue(param))
		{
			switch (param)
			{
			case cv::ml::SVM::ParamTypes::C: {
				double c = svm_description.GetC();
				if (c != 1) {
					node[YAML_C] = c;
				}
			} break;
			case cv::ml::SVM::ParamTypes::GAMMA: {
				double gamma = svm_description.GetGamma();
				if (gamma != 1) {
					node[YAML_GAMMA] = gamma;
				}
			} break;
			case cv::ml::SVM::ParamTypes::P: {
				double p = svm_description.GetP();
				if (p != 0) {
					node[YAML_P] = p;
				}
			} break;
			case cv::ml::SVM::ParamTypes::NU: {
				double nu = svm_description.GetNu();
				if (nu != 0) {
					node[YAML_NU] = nu;
				}
			} break;
			case cv::ml::SVM::ParamTypes::COEF: {
				double coef = svm_description.GetCoef();
				if (coef != 0) {
					node[YAML_COEF] = coef;
				}
			} break;
			case cv::ml::SVM::ParamTypes::DEGREE: {
				double degree = svm_description.GetDegree();
				if (degree != 0) {
					node[YAML_DEGREE] = degree;
				}
			} break;
			default:
				//invalid
				break;
			}
		}
		else {
			auto grid = svm_description.GetParamGridForParameter(param);
			auto default_grid = cv::ml::SVM::getDefaultGrid(param);

			if ((grid.minVal == default_grid.minVal) && (grid.maxVal == default_grid.maxVal) && (grid.logStep == default_grid.logStep)) {
				node[StringFromParamType(param)] = YAML_AUTO;
			}
			else {
				YAML::Node map;
				map[YAML_MIN_VAL] = grid.minVal;
				map[YAML_MAX_VAL] = grid.maxVal;
				map[YAML_LOG_STEP] = grid.logStep;

				node[StringFromParamType(param)] = map;
			}
		}
	}

	// For Term Criteria

	std::string termCriteriaType = StringFromTermCriteriaType(svm_description.GetTermCriteria().type);
	int         termCriteriaMax = svm_description.GetTermCriteria().maxCount;
	double      termCriteriaEps = svm_description.GetTermCriteria().epsilon;

	if ((termCriteriaType != "MAX_ITER+EPS") && (termCriteriaMax != 1000) && (termCriteriaEps != FLT_EPSILON)) {
		YAML::Node map;
		map[YAML_TERM_CRITERIA_TYPE] = termCriteriaType;
		map[YAML_TERM_CRITERIA_MAX] = termCriteriaMax;
		map[YAML_TERM_CRITERIA_EPS] = termCriteriaEps;
		node[YAML_TERM_CRITERIA] = map;
	}

	return node;
}