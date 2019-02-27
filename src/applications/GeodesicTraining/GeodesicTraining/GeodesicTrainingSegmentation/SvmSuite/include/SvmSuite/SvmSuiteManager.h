#ifndef H_CBICA_SVM_SUITE_MANAGER
#define H_CBICA_SVM_SUITE_MANAGER

#include <opencv2/ml.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <map>
#include <unordered_map>
#include <thread>

#include "SvmSuiteDescription.h"
#include "SvmSuiteOperations.h"
#include "SvmSuiteUtil.h"

namespace SvmSuite
{
	class Manager
	{
	public:

		typedef int   LabelsType;
		typedef float PseudoProbType;

		explicit Manager() {}

		virtual ~Manager() {}

		typedef struct Result {
			cv::Mat posMat;
			cv::Mat negMat;
			LabelsType posLabel = 0;
			LabelsType negLabel = 0;
			cv::Mat labelsMat;
		} Result;

		void GenerateConfigFromBestValues(std::string outputFileName);

		void GenerateConfigFromBestValues();

		/**
		Train the models specified in the svm descriptions
		*/
		void Train();

		/**
		Test using the ensemble of trained svm models
		@param testingMat cv::Mat where the rows are the samples to test and columns are the features
		@param pseudoProbMapResult true->  the result would be pos and neg pseudoprobability maps (only for n=2 classification)
		                           false-> the result would be the predicted labels (n>=2 classification)
		@return pointer to SvmSuite::Manager::Result object (which contains the result images)
		*/
		std::shared_ptr<Result> Test(cv::Mat &testingMat, bool pseudoProbMapResult = false);

		/**
		Test using the ensemble of trained svm models
		@param testingMat cv::Mat where the rows are the samples to test and columns are the features
		@param skipZerosMat cv::Mat that has as many rows as testingMat, where for each row if skipZerosMat has zero the line will be skipped
		                    and the label 0 will be set. In practice you can pass testingMat two times (extra columns would not matter)
		@param pseudoProbMapResult true->  the result would be pos and neg pseudoprobability maps (only for n=2 classification)
								   false-> the result would be the predicted labels (n>=2 classification)
		@param skipZeros whether to user skipZerosMat
		@return pointer to SvmSuite::Manager::Result object (which contains the result images)
		*/
		std::shared_ptr<Result> Test(cv::Mat &testingMat, cv::Mat &skipZerosMat, bool pseudoProbMapResult = false, bool skipZeros = true);

		// Adders and Setters

		void AddSvmDescriptionToList(SvmDescription svmDesc);
		void AddSvmDescriptions(std::vector< SvmDescription > svmDescs);
		void AddPretrainedModel(std::string pretrainedModelPath, int neighborhoodRadius = 0);
		void AddSvmsFromConfig(std::string configPath);
		void SetTrainData(cv::Mat &trainingMat, cv::Mat &labelsMat, cv::Mat &weightsMat/*, cv::Mat sampleIdx = cv::Mat()*/);
		void SetVerbose(bool verbose);
		void SetOutputPath(std::string path);
		void SetSavingModelsEnabled(bool modelsEnabled);
		void SetTimerEnabled(bool timerEnabled);
		void SetNumberOfThreads(int numberOfThreads);
		void SetSubsampling(bool subsample, int maxSamples = 3000);
		void SetInputNormalization(bool normalize);
		
	private:
		std::vector< SvmDescription > m_svm_descriptions;
		cv::Ptr<cv::ml::TrainData>    m_traindata;
		cv::Mat                       m_weights_mat;
		std::set<LabelsType>          m_different_labels;
		std::string                   m_output_path = "./";
		SvmSuiteUtil::Timer           m_timer;
		int                           m_number_of_threads = 32, m_max_samples = 3000;
		bool m_timer_enabled = false, m_save_models = false, m_verbose = false, m_subsample = false, m_normalize = false;
		
		void testingLabelsThreadJob(cv::Mat &testingMat, int iStart, int interval, cv::Mat &resultLabelsMat,
		                            bool skipZeros, cv::Mat &skipZerosMat, std::vector< double > &importanceValues);

		// For timer

		void startTimer();

		void stopTimerAndReport(std::string desc);

		// Util

		void message(std::string message, bool overdraw = false, bool finished = false, int progress = -1);

		void errorOccured(std::string msg);

		template< typename U, typename V >
		std::vector< U > getMapKeyset(std::map< U, V > map)
		{
			// Find a list of the different keys
			std::vector< U > res;

			res.reserve(map.size());
			for (auto const& imap : map) {
				res.push_back(imap.first);
			}

			return res;
		}

	};

}

#endif