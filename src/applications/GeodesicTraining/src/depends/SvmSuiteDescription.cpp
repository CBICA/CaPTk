#include "SvmSuiteDescription.h"

#include "ConvertionsOpenCV.h"

SvmSuite::SvmDescription::SvmDescription()
{
	// Initialize the ParamGrid for each parameter
	parametersRanges.push_back(cv::ml::ParamGrid(1, 1, 1)); //C
	parametersRanges.push_back(cv::ml::ParamGrid(1, 1, 1)); //GAMMA
	parametersRanges.push_back(cv::ml::ParamGrid(0, 0, 1)); //P
	parametersRanges.push_back(cv::ml::ParamGrid(0, 0, 1)); //NU
	parametersRanges.push_back(cv::ml::ParamGrid(0, 0, 1)); //COEF
	parametersRanges.push_back(cv::ml::ParamGrid(0, 0, 1)); //DEGREE

	// Initialize the flag vector for setting specific value to a parameter
	isParameterSetToSpecificValueVector.resize(6, true);
}

void SvmSuite::SvmDescription::SetParameterRangeAuto(cv::ml::SVM::ParamTypes param)
{
	isParameterSetToSpecificValueVector[param] = false;
	parametersRanges[param] = cv::ml::SVM::getDefaultGrid(param);
}

void SvmSuite::SvmDescription::SetParameterRangeAuto(std::string param)
{
	SvmDescription::SetParameterRangeAuto(SvmSuiteConvertions::ParamTypeFromString(param));
}

void SvmSuite::SvmDescription::SetParameterRange(cv::ml::SVM::ParamTypes param, double minVal, double maxVal, double logStep)
{
	isParameterSetToSpecificValueVector[param] = false;
	parametersRanges[param].minVal = minVal;
	parametersRanges[param].maxVal = maxVal;
	parametersRanges[param].logStep = logStep;
}

void SvmSuite::SvmDescription::SetParameterRange(std::string param, double minVal, double maxVal, double logStep)
{
	SvmDescription::SetParameterRange(SvmSuiteConvertions::ParamTypeFromString(param), minVal, maxVal, logStep);
}

void SvmSuite::SvmDescription::SetTermCriteria(int criteriaType, int maxCount, double epsilon)
{
	termCriteria = cv::TermCriteria(criteriaType, maxCount, epsilon);
}

void SvmSuite::SvmDescription::SetTermCriteria(std::string criteriaType, int maxCount, double epsilon)
{
	termCriteria = cv::TermCriteria(SvmSuiteConvertions::TermCriteriaTypeFromString(criteriaType), maxCount, epsilon);
}

void SvmSuite::SvmDescription::SetTermCriteria(std::string criteriaType, std::string maxCount, std::string epsilon)
{
	termCriteria = cv::TermCriteria(SvmSuiteConvertions::TermCriteriaTypeFromString(criteriaType), std::stoi(maxCount), std::stod(epsilon));
}

cv::ml::SVM::KernelTypes SvmSuite::SvmDescription::GetKernelType()
{
	return kernelType;
}

std::string SvmSuite::SvmDescription::GetKernelTypeAsString() {
	return SvmSuiteConvertions::StringFromKernelType(this->kernelType);
}

cv::ml::SVM::Types SvmSuite::SvmDescription::GetType()
{
	return type;
}

std::string SvmSuite::SvmDescription::GetTypeAsString() {
	return SvmSuiteConvertions::StringFromType(this->type);
}

int SvmSuite::SvmDescription::GetKfold()
{
	return kfold;
}

int SvmSuite::SvmDescription::GetNeighborhoodRadius()
{
	return neighborhoodRadius;
}

bool SvmSuite::SvmDescription::GetConsiderWeights()
{
	return considerWeights;
}

double SvmSuite::SvmDescription::GetImportance()
{
	return importance;
}

std::string SvmSuite::SvmDescription::GetModelPath()
{
	return this->modelPath;
}

cv::Ptr<cv::ml::SVM> SvmSuite::SvmDescription::GetModel()
{
	return this->model;
}

double SvmSuite::SvmDescription::GetC()
{
	return c;
}

double SvmSuite::SvmDescription::GetGamma()
{
	return gamma;
}

double SvmSuite::SvmDescription::GetP()
{
	return p;
}

double SvmSuite::SvmDescription::GetNu()
{
	return nu;
}

double SvmSuite::SvmDescription::GetCoef()
{
	return coef;
}

double SvmSuite::SvmDescription::GetDegree()
{
	return degree;
}

double SvmSuite::SvmDescription::GetParameter(cv::ml::SVM::ParamTypes param)
{
	switch (param)
	{
	case cv::ml::SVM::C:
		return c;
		break;
	case cv::ml::SVM::GAMMA:
		return gamma;
		break;
	case cv::ml::SVM::P:
		return p;
		break;
	case cv::ml::SVM::NU:
		return nu;
		break;
	case cv::ml::SVM::COEF:
		return coef;
		break;
	case cv::ml::SVM::DEGREE:
		return degree;
		break;
	default:
		//invalid
		return c;
		break;
	}
}

double SvmSuite::SvmDescription::GetParameter(std::string param)
{
	return SvmDescription::GetParameter(SvmSuiteConvertions::ParamTypeFromString(param));
}

bool SvmSuite::SvmDescription::isParameterSetToSpecificValue(cv::ml::SVM::ParamTypes param)
{
	return isParameterSetToSpecificValueVector[param];
}

bool SvmSuite::SvmDescription::isParameterSetToSpecificValue(std::string param)
{
	return isParameterSetToSpecificValueVector[SvmSuiteConvertions::ParamTypeFromString(param)];
}

cv::ml::ParamGrid SvmSuite::SvmDescription::GetParamGridForParameter(cv::ml::SVM::ParamTypes param)
{
	return parametersRanges[param];
}

cv::ml::ParamGrid SvmSuite::SvmDescription::GetParamGridForParameter(std::string param)
{
	return SvmDescription::GetParamGridForParameter(SvmSuiteConvertions::ParamTypeFromString(param));
}

cv::TermCriteria SvmSuite::SvmDescription::GetTermCriteria()
{
	return termCriteria;
}

std::string SvmSuite::SvmDescription::GetTermCriteriaTypeAsString() {
	return SvmSuiteConvertions::StringFromTermCriteriaType(this->termCriteria.type);
}

void SvmSuite::SvmDescription::SetKernelType(cv::ml::SVM::KernelTypes kernelType)
{
	this->kernelType = kernelType;
}

void SvmSuite::SvmDescription::SetKernelType(std::string kernelType) {
	this->kernelType = SvmSuiteConvertions::KernelTypeFromString(kernelType);
}

void SvmSuite::SvmDescription::SetType(cv::ml::SVM::Types type)
{
	this->type = type;
}

void SvmSuite::SvmDescription::SetType(std::string type)
{
	this->type = SvmSuiteConvertions::TypeFromString(type);
}

void SvmSuite::SvmDescription::SetKfold(int kfold)
{
	this->kfold = kfold;
}

void SvmSuite::SvmDescription::SetKfold(std::string kfold)
{
	this->kfold = std::stoi(kfold);
}

void SvmSuite::SvmDescription::SetNeighborhoodRadius(int neighborhoodRadius)
{
	this->neighborhoodRadius = neighborhoodRadius;
}

void SvmSuite::SvmDescription::SetNeighborhoodRadius(std::string neighborhoodRadius)
{
	this->neighborhoodRadius = std::stoi(neighborhoodRadius);
}

void SvmSuite::SvmDescription::SetConsiderWeights(bool considerWeights)
{
	this->considerWeights = considerWeights;
}

void SvmSuite::SvmDescription::SetConsiderWeights(std::string considerWeights)
{
	this->considerWeights = (considerWeights == "true") ? true : false;
}

void SvmSuite::SvmDescription::SetImportance(double importance)
{
	this->importance = importance;
}

void SvmSuite::SvmDescription::SetImportance(std::string importance)
{
	this->importance = std::stod(importance);
}

void SvmSuite::SvmDescription::SetModelPath(std::string modelPath)
{
	this->modelPath = modelPath;
}

void SvmSuite::SvmDescription::SetModel(cv::Ptr<cv::ml::SVM> model)
{
	this->model = model;
}

void SvmSuite::SvmDescription::SetParameter(cv::ml::SVM::ParamTypes param, double val)
{
	this->isParameterSetToSpecificValueVector[param] = true;

	// Set the ParamGrid values for trainAuto (don't use SetParameterRanges(), it changes isParameterSetToSpecificValueVector)
	this->parametersRanges[param].minVal = val;
	this->parametersRanges[param].maxVal = val + 1;
	this->parametersRanges[param].logStep = 2;

	switch (param) {
	case cv::ml::SVM::C:
		c = val;
		break;
	case cv::ml::SVM::GAMMA:
		gamma = val;
		break;
	case cv::ml::SVM::P:
		p = val;
		break;
	case cv::ml::SVM::NU:
		nu = val;
		break;
	case cv::ml::SVM::COEF:
		coef = val;
		break;
	case cv::ml::SVM::DEGREE:
		degree = val;
		break;
	default:
		//invalid
		break;
	}
}

void SvmSuite::SvmDescription::SetParameter(std::string param, double val)
{
	this->SetParameter(SvmSuiteConvertions::ParamTypeFromString(param), val);
}

void SvmSuite::SvmDescription::SetParameter(std::string param, std::string val)
{
	this->SetParameter(SvmSuiteConvertions::ParamTypeFromString(param), std::stod(val));
}

void SvmSuite::SvmDescription::SetC(double c)
{
	this->SetParameter(cv::ml::SVM::C, c);
}

void SvmSuite::SvmDescription::SetC(std::string c)
{
	this->SetParameter(cv::ml::SVM::C, std::stod(c));
}

void SvmSuite::SvmDescription::SetGamma(double gamma)
{
	this->SetParameter(cv::ml::SVM::GAMMA, gamma);
}

void SvmSuite::SvmDescription::SetGamma(std::string gamma)
{
	this->SetParameter(cv::ml::SVM::GAMMA, std::stod(gamma));
}

void SvmSuite::SvmDescription::SetP(double p)
{
	this->SetParameter(cv::ml::SVM::P, p);
}

void SvmSuite::SvmDescription::SetP(std::string p)
{
	this->SetParameter(cv::ml::SVM::P, std::stod(p));
}

void SvmSuite::SvmDescription::SetNu(double nu)
{
	this->SetParameter(cv::ml::SVM::NU, nu);
}

void SvmSuite::SvmDescription::SetNu(std::string nu)
{
	this->SetParameter(cv::ml::SVM::NU, std::stod(nu));
}

void SvmSuite::SvmDescription::SetCoef(double coef)
{
	this->SetParameter(cv::ml::SVM::COEF, coef);
}

void SvmSuite::SvmDescription::SetCoef(std::string coef)
{
	this->SetParameter(cv::ml::SVM::COEF, std::stod(coef));
}

void SvmSuite::SvmDescription::SetDegree(double degree)
{
	this->SetParameter(cv::ml::SVM::DEGREE, degree);
}

void SvmSuite::SvmDescription::SetDegree(std::string degree)
{
	this->SetParameter(cv::ml::SVM::DEGREE, std::stod(degree));
}

std::vector<SvmSuite::SvmDescription> SvmSuite::SvmDescription::GetDefaultSvmDescriptions()
{
	SvmDescription rbfDesc;
	rbfDesc.SetKernelType(cv::ml::SVM::KernelTypes::RBF);
	// rbfDesc.SetParameterRangeAuto(cv::ml::SVM::C);
	rbfDesc.SetParameterRange(cv::ml::SVM::C, 1.0, 400, 2.5);
	rbfDesc.SetParameterRangeAuto(cv::ml::SVM::GAMMA);
	rbfDesc.SetConsiderWeights(true);
	rbfDesc.SetImportance(0.4);

	SvmDescription chi2Desc;
	chi2Desc.SetKernelType(cv::ml::SVM::KernelTypes::CHI2);
	chi2Desc.SetParameterRange(cv::ml::SVM::C, 1.0, 400, 2.5);
	// chi2Desc.SetParameterRangeAuto(cv::ml::SVM::C);
	chi2Desc.SetParameterRangeAuto(cv::ml::SVM::GAMMA);
	chi2Desc.SetConsiderWeights(true);
	chi2Desc.SetImportance(0.3);

	SvmDescription interDesc;
	interDesc.SetKernelType(cv::ml::SVM::KernelTypes::INTER);
	interDesc.SetParameterRangeAuto(cv::ml::SVM::C);
	interDesc.SetConsiderWeights(true);
	interDesc.SetImportance(0.3);

	std::vector<SvmDescription> defDescs;
	defDescs.push_back(rbfDesc);
	defDescs.push_back(chi2Desc);
	defDescs.push_back(interDesc);

	return defDescs;
}

void SvmSuite::SvmDescription::errorOccured(std::string msg)
{
	std::cerr << "Manager error: " << msg << std::endl;
}