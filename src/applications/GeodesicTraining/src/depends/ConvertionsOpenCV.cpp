#include "ConvertionsOpenCV.h"

cv::ml::SVM::Types SvmSuiteConvertions::TypeFromString(std::string typeString) {
	if (typeString == "C_SVC") {
		return cv::ml::SVM::Types::C_SVC;
	}
	else if (typeString == "NU_SVC") {
		return cv::ml::SVM::Types::NU_SVC;
	}
	else if (typeString == "ONE_CLASS") {
		return cv::ml::SVM::Types::ONE_CLASS;
	}
	else if (typeString == "EPS_SVR") {
		return cv::ml::SVM::Types::EPS_SVR;
	}
	else if (typeString == "NU_SVR") {
		return cv::ml::SVM::Types::NU_SVR;
	}
	else {
		//invalid
		return cv::ml::SVM::Types::C_SVC;
	}
}

std::string SvmSuiteConvertions::StringFromType(cv::ml::SVM::Types type) {
	switch (type) {
	case cv::ml::SVM::Types::C_SVC:
		return "C_SVC";
	case cv::ml::SVM::Types::NU_SVC:
		return "NU_SVC";
	case cv::ml::SVM::Types::ONE_CLASS:
		return "ONE_CLASS";
	case cv::ml::SVM::Types::EPS_SVR:
		return "EPS_SVR";
	case cv::ml::SVM::Types::NU_SVR:
		return "NU_SVR";
	default:
		return "invalid";
	}
}

cv::ml::SVM::KernelTypes SvmSuiteConvertions::KernelTypeFromString(std::string kernelTypeString) {
	if (kernelTypeString == "linear") {
		return cv::ml::SVM::KernelTypes::LINEAR;
	}
	else if (kernelTypeString == "rbf") {
		return cv::ml::SVM::KernelTypes::RBF;
	}
	else if (kernelTypeString == "poly") {
		return cv::ml::SVM::KernelTypes::POLY;
	}
	else if (kernelTypeString == "sigmoid") {
		return cv::ml::SVM::KernelTypes::SIGMOID;
	}
	else if (kernelTypeString == "chi2") {
		return cv::ml::SVM::KernelTypes::CHI2;
	}
	else if (kernelTypeString == "inter") {
		return cv::ml::SVM::KernelTypes::INTER;
	}
	else {
		//invalid
		return cv::ml::SVM::KernelTypes::RBF;
	}
}

std::string SvmSuiteConvertions::StringFromKernelType(cv::ml::SVM::KernelTypes kernelType) {
	switch (kernelType)
	{
	case cv::ml::SVM::LINEAR:
		return "linear";
	case cv::ml::SVM::RBF:
		return "rbf";
	case cv::ml::SVM::POLY:
		return "poly";
	case cv::ml::SVM::SIGMOID:
		return "sigmoid";
	case cv::ml::SVM::CHI2:
		return "chi2";
	case cv::ml::SVM::INTER:
		return "inter";
	default:
		return "invalid";
	}
}

cv::ml::SVM::ParamTypes SvmSuiteConvertions::ParamTypeFromString(std::string paramTypeString) {
	if (paramTypeString == "c") {
		return cv::ml::SVM::ParamTypes::C;
	}
	else if (paramTypeString == "gamma") {
		return cv::ml::SVM::ParamTypes::GAMMA;
	}
	else if (paramTypeString == "p") {
		return cv::ml::SVM::ParamTypes::P;
	}
	else if (paramTypeString == "nu") {
		return cv::ml::SVM::ParamTypes::NU;
	}
	else if (paramTypeString == "coef") {
		return cv::ml::SVM::ParamTypes::COEF;
	}
	else if (paramTypeString == "degree") {
		return cv::ml::SVM::ParamTypes::DEGREE;
	}
	else {
		//invalid
		return cv::ml::SVM::ParamTypes::GAMMA;
	}
}

std::string SvmSuiteConvertions::StringFromParamType(cv::ml::SVM::ParamTypes paramType) {
	switch (paramType)
	{
	case cv::ml::SVM::C:
		return "c";
	case cv::ml::SVM::GAMMA:
		return "gamma";
	case cv::ml::SVM::P:
		return "p";
	case cv::ml::SVM::NU:
		return "nu";
	case cv::ml::SVM::COEF:
		return "coef";
	case cv::ml::SVM::DEGREE:
		return "degree";
	default:
		return "invalid";
	}
}

int SvmSuiteConvertions::TermCriteriaTypeFromString(std::string termCriteriaTypeString) {
	if (termCriteriaTypeString == "MAX_ITER") {
		return cv::TermCriteria::MAX_ITER;
	}
	else if (termCriteriaTypeString == "EPS") {
		return cv::TermCriteria::EPS;
	}
	else if (termCriteriaTypeString == "MAX_ITER+EPS") {
		return cv::TermCriteria::MAX_ITER + cv::TermCriteria::EPS;
	}
	else {
		//invalid
		return cv::TermCriteria::MAX_ITER + cv::TermCriteria::EPS;
	}
}

std::string SvmSuiteConvertions::StringFromTermCriteriaType(int termCriteriaType) {
	switch (termCriteriaType)
	{
	case cv::TermCriteria::MAX_ITER:
		return "MAX_ITER";
	case cv::TermCriteria::EPS:
		return "EPS";
	case (cv::TermCriteria::MAX_ITER + cv::TermCriteria::EPS):
		return "MAX_ITER+EPS";
	default:
		return "invalid";
	}
}