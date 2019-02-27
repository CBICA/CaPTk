#ifndef H_CBICA_SVM_SUITE_DESCRIPTION
#define H_CBICA_SVM_SUITE_DESCRIPTION

#include <opencv2/ml.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <string>

namespace SvmSuite 
{
	// Each SVM in the configuration (there might be more than one) gets stored in one instance of SvmDescription
	class SvmDescription
	{
	public:
		/**
		Constructor (initializes ParamGrids)
		*/
		SvmDescription();

		~SvmDescription() {}

		// Setters

		void SetKernelType(cv::ml::SVM::KernelTypes kernelType);
		void SetKernelType(std::string kernelType);
		void SetType(cv::ml::SVM::Types type);
		void SetType(std::string type);
		void SetKfold(int kfold);
		void SetKfold(std::string kfold);
		void SetNeighborhoodRadius(int neighborhoodRadius);
		void SetNeighborhoodRadius(std::string neighborhoodRadius);
		void SetConsiderWeights(bool considerWeights);
		void SetConsiderWeights(std::string considerWeights);
		void SetImportance(double importance);
		void SetImportance(std::string importance);
		void SetModelPath(std::string modelPath);
		void SetModel(cv::Ptr<cv::ml::SVM> model);

		// Note: setting a parameters to value X will also change the paramgrid to {X, X+1, 2} (that way only value X is tried in OpenCV's trainAuto)
		
		void SetParameter(cv::ml::SVM::ParamTypes param, double val);
		void SetParameter(std::string param, double val);
		void SetParameter(std::string param, std::string val);

		void SetC(double c);
		void SetC(std::string c);
		void SetGamma(double gamma);
		void SetGamma(std::string gamma);
		void SetP(double p);
		void SetP(std::string p);
		void SetNu(double nu);
		void SetNu(std::string nu);
		void SetCoef(double coef);
		void SetCoef(std::string coef);
		void SetDegree(double degree);
		void SetDegree(std::string degree);

		void SetParameterRangeAuto(cv::ml::SVM::ParamTypes param);
		void SetParameterRangeAuto(std::string param);

		void SetParameterRange(cv::ml::SVM::ParamTypes param, double minVal, double maxVal, double logStep);
		void SetParameterRange(std::string param, double minVal, double maxVal, double logStep);

		void SetTermCriteria(int criteriaType, int maxCount, double epsilon);
		void SetTermCriteria(std::string criteriaType, int maxCount, double epsilon);
		void SetTermCriteria(std::string criteriaType, std::string maxCount, std::string epsilon);

		// Getters

		cv::ml::SVM::KernelTypes GetKernelType();
		std::string GetKernelTypeAsString();
		cv::ml::SVM::Types GetType();
		std::string GetTypeAsString();
		int GetKfold();
		int GetNeighborhoodRadius();
		bool GetConsiderWeights();
		double GetImportance();
		std::string GetModelPath();
		cv::Ptr<cv::ml::SVM> GetModel();

		double GetC();
		double GetGamma();
		double GetP();
		double GetNu();
		double GetCoef();
		double GetDegree();

		double GetParameter(cv::ml::SVM::ParamTypes param);
		double GetParameter(std::string param);
		bool isParameterSetToSpecificValue(cv::ml::SVM::ParamTypes param);
		bool isParameterSetToSpecificValue(std::string param);

		cv::ml::ParamGrid GetParamGridForParameter(cv::ml::SVM::ParamTypes param);
		cv::ml::ParamGrid GetParamGridForParameter(std::string param);
		cv::TermCriteria GetTermCriteria();
		std::string GetTermCriteriaTypeAsString();

		/**
		Get default ensemble of SvmDescriptions
		@return Default vector of SvmDescriptions
		*/
		static std::vector<SvmDescription> GetDefaultSvmDescriptions();

	private:
		// Everything is initialized with OpenCV's default values

		cv::ml::SVM::KernelTypes kernelType = cv::ml::SVM::KernelTypes::RBF;
		cv::ml::SVM::Types type = cv::ml::SVM::Types::C_SVC;
		std::vector< cv::ml::ParamGrid > parametersRanges;
		std::vector< bool > isParameterSetToSpecificValueVector;
		cv::TermCriteria termCriteria = cv::TermCriteria(cv::TermCriteria::MAX_ITER + cv::TermCriteria::EPS, 1000, FLT_EPSILON);
		bool considerWeights = false;
		int    kfold = 10;
		int	   neighborhoodRadius = 0;
		double importance = 1.0;
		std::string modelPath = "";
		cv::Ptr<cv::ml::SVM> model;
		double c = 1.0;
		double gamma = 1.0;
		double p = 0.0;
		double nu = 0.0;
		double coef = 0.0;
		double degree = 0.0;

		void errorOccured(std::string msg);
	};

}

#endif