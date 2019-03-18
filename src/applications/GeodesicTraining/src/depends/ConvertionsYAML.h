#ifndef H_CBICA_SVM_SUITE_CONVERTIONS_YAML
#define H_CBICA_SVM_SUITE_CONVERTIONS_YAML

#include <opencv2/ml.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <string>

#include "yaml-cpp/yaml.h"

#include "SvmSuiteDescription.h"

const std::string YAML_KERNEL_TYPE = "kernel_type";
const std::string YAML_TYPE = "type";
const std::string YAML_KFOLD = "kfold";
const std::string YAML_NEIGHBORHOOD_RADIUS = "neighborhood_radius";
const std::string YAML_CONSIDER_WEIGHTS = "consider_weights";
const std::string YAML_IMPORTANCE = "importance";
const std::string YAML_MODEL_PATH = "pretrained";
const std::string YAML_C = "c";
const std::string YAML_GAMMA = "gamma";
const std::string YAML_P = "p";
const std::string YAML_NU = "nu";
const std::string YAML_COEF = "coef";
const std::string YAML_DEGREE = "degree";
const std::string YAML_AUTO = "auto";
const std::string YAML_TERM_CRITERIA = "term_criteria";
const std::string YAML_TERM_CRITERIA_TYPE = "criteria_type";
const std::string YAML_TERM_CRITERIA_MAX = "max_count";
const std::string YAML_TERM_CRITERIA_EPS = "epsilon";
const std::string YAML_MIN_VAL = "min_value";
const std::string YAML_MAX_VAL = "max_value";
const std::string YAML_LOG_STEP = "log_step";
const std::string YAML_ROOT_NODE = "svms";

namespace SvmSuiteConvertions
{
	YAML::Node yamlConvertSvmDescriptionToNode(SvmSuite::SvmDescription& svm_description);
}

#endif