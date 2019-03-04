#include <opencv2/ml.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <cmath>

#include "ConfigParserRF.h"
#include "RFPrepareTrainData.h"

namespace RFSuite
{
	class Manager
	{
	public:
		const double DEFAULT_TRAINING_SAMPLE_PERCENTAGE = 75.0;
		const int    DEFAULT_MAX_DEPTH = 10;
		const double DEFAULT_MIN_SAMPLE_COUNT_PERCENTAGE = 10.0;
		const int    DEFAULT_MAX_CATEGORIES = 24;
		const int    DEFAULT_ACTIVE_VAR_COUNT = 0; // 0 is the default and it is sqrt(# of features)
		const int    DEFAULT_NUMBER_OF_TREES = 150;

		explicit Manager() {}

		virtual ~Manager() {}

		/**
		Trains while tuning parameters with prespecified ranges
		@return lower training error
		*/
		float TrainAuto();

		/**
		Trains while tuning parameters with ranges specified in the parameters
		Every ParamGrid has fields: minVal, maxVal, logStep(>=1).
		if (logStep>1) then values: minVal, minVal*logStep^1, minVal*logStep^2,... are tried
		while minVal*logStep^n < maxVal (see definition of cv::ml::ParamGrid for more)
		if (logStep=1) then values: minVal, minVal+1, minVal+2,... are tried
		while minVal+n < maxVal
		Note: If there is no need for tuning for a parameter and value X is to be used,
		give a ParamGrid with minVal=X, maxVal=X+1, logStep=2 (only X will be tried)
		@param maxDepthGrid ParamGrid for max depth
		@param minSampleCountPercentageGrid ParamGrid for min sample count percentage
		@param maxCategoriesGrid ParamGrid for max categories
		@param activeVarCountGrid ParamGrid for active variable count
		@return lower training error
		*/
		float TrainAuto(cv::ml::ParamGrid maxDepthGrid,      cv::ml::ParamGrid minSampleCountPercentageGrid,
		                cv::ml::ParamGrid maxCategoriesGrid, cv::ml::ParamGrid activeVarCountGrid);

		/**
		Trains will values specified by setters
		@return training error
		*/
		float Train();

		/**
		Trains with specified parameters 
		Note: Training sample percentage should still be changed from SetTrainingSamplePercentage()
		@param maxDepth the max depth for the decision trees
		@param minSampleCountPercentage percentage of training data to be randomly selected for a tree
		@param maxCategories Max categories in training for each class (Relevant only in n>2 classification)
		@param activeVarCount subset of features for each tree
		@param numberOfTrees size of the forest
		@return training error
		*/
		float Train(int maxDepth, double minSampleCountPercentage, int maxCategories, int activeVarCount, int numberOfTrees);

		/**
		Make predictions
		@param testingMat matrix with one sample per row, each column is a feature
		@return pointer to the responses mat
		*/
		std::shared_ptr<cv::Mat> Test(cv::Mat &testingMat);

		/**
		Make predictions
		@param testingMat matrix with one sample per row, each column is a feature
		@param skipZerosMat matrix with the same number of rows as trainingMat (if val at [i][0] is 0 -> no prediction will happen)
		@param skipZeros whether skipZerosMat will be used
		@return pointer to the responses mat
		*/
		std::shared_ptr<cv::Mat> Test(cv::Mat &testingMat, cv::Mat &skipZerosMat, bool skipZeros = true);

		/**
		Save random forest model to file
		@param filename full path to desired save location
		*/
		void SaveModel(const std::string filename);

		// Setters

		void SetTrainDataFromMats(cv::Mat &trainingMat, cv::Mat &labelsMat);
		
		void SetPriorsMat(cv::Mat &priorsMat);
		
		void SetOutputPath(std::string path);
		
		void SetSaveAll(bool saveAll);

		void SetVerbose(bool verbose);
		
		void SetParametersFromConfig(std::string filePath);

		void SetTrainingSamplePercentage(double trainingSamplePercentage);

		void SetMaxDepth(int maxDepth);

		void SetMinSampleCountPercentage(double minSampleCountPercentage);
		
		void SetMaxCategories(int maxCategories);

		void SetActiveVarCount(int activeVarCount);

		void SetNumberOfTrees(int numberOfTrees);

	private:
		cv::Ptr<cv::ml::RTrees> m_rtrees = cv::ml::RTrees::create();
		cv::Mat m_priors_mat = cv::Mat();
		cv::Ptr<cv::ml::TrainData> m_traindata;
		std::string m_output_path = "./", m_rf_config_file_path = "";
		bool m_save_all = false, m_verbose = false;

		// Parameters
		//    Note on m_training_sample_percentage:
		//    From the training data, only m_training_sample_percentage % will actually be used for training
		//    The rest will be used for checking the accuracy of the model
		double m_training_sample_percentage  = DEFAULT_TRAINING_SAMPLE_PERCENTAGE;
		int    m_max_depth                   = DEFAULT_MAX_DEPTH;
		double m_min_sample_count_percentage = DEFAULT_MIN_SAMPLE_COUNT_PERCENTAGE;
		int    m_max_categories              = DEFAULT_MAX_CATEGORIES;
		int    m_active_var_count            = DEFAULT_ACTIVE_VAR_COUNT;
		int    m_number_of_trees             = DEFAULT_NUMBER_OF_TREES;

		/**
		Trains the classifier
		@return training error
		*/
		float train_and_print_errs(cv::Ptr<cv::ml::StatModel> model, const cv::Ptr<cv::ml::TrainData>& data);
	};
}
