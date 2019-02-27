#include "RFSuiteManager.h"

float RFSuite::Manager::TrainAuto()
{
	cv::ml::ParamGrid maxDepthGrid = cv::ml::ParamGrid(2, 19, 1.75);
	cv::ml::ParamGrid minSampleCountPercentageGrid = cv::ml::ParamGrid(10, 11, 2); //5,27,1.75
	cv::ml::ParamGrid maxCategoriesGrid = cv::ml::ParamGrid(16, 17, 2);
	cv::ml::ParamGrid activeVarCountGrid = cv::ml::ParamGrid(1, m_traindata->getTrainSamples().cols + 1, 1);

	return TrainAuto(maxDepthGrid, minSampleCountPercentageGrid, maxCategoriesGrid, activeVarCountGrid);
}

float RFSuite::Manager::TrainAuto(cv::ml::ParamGrid maxDepthGrid, cv::ml::ParamGrid minSampleCountPercentageGrid,
	cv::ml::ParamGrid maxCategoriesGrid, cv::ml::ParamGrid activeVarCountGrid)
{
	if (m_verbose) {
		std::cout << "RF Manager:\t Training random forest while tuning parameters...\n";
	}

	if (maxDepthGrid.logStep < 1 || minSampleCountPercentageGrid.logStep < 1 ||
		maxCategoriesGrid.logStep < 1 || activeVarCountGrid.logStep < 1)
	{
		std::cerr << "logStep should be >= 1. (if 1 then the values: minVal, minVal+1, minVal+2... are tried, else see ParamGrid definition)\n";
		return -1;
	}

	float minError = 100;
	cv::Ptr<cv::ml::RTrees> best_rforest; // = cv::ml::RTrees::create();

	int i_maxDepth = 0;

	for (double maxDepth = maxDepthGrid.minVal;
		((maxDepthGrid.logStep == 1) ?
		(maxDepth < maxDepthGrid.maxVal) :
			(maxDepthGrid.minVal * std::pow(maxDepthGrid.logStep, i_maxDepth) < maxDepthGrid.maxVal));
		i_maxDepth++)
	{
		int i_minSampleCountPercentage = 0;
		for (double minSampleCountPercentage = minSampleCountPercentageGrid.minVal;
			((minSampleCountPercentageGrid.logStep == 1) ?
			(minSampleCountPercentage < minSampleCountPercentageGrid.maxVal) :
				(minSampleCountPercentageGrid.minVal * std::pow(minSampleCountPercentageGrid.logStep, i_minSampleCountPercentage) <
					minSampleCountPercentageGrid.maxVal));
			i_minSampleCountPercentage++)
		{
			int i_maxCategories = 0;

			for (double maxCategories = maxCategoriesGrid.minVal;
				((maxCategoriesGrid.logStep == 1) ?
				(maxCategories < maxCategoriesGrid.maxVal) :
					(maxCategoriesGrid.minVal * std::pow(maxCategoriesGrid.logStep, i_maxCategories) < maxCategoriesGrid.maxVal));
				i_maxCategories++)
			{
				int i_activeVarCount = 0;

				for (double activeVarCount = activeVarCountGrid.minVal;
					((activeVarCountGrid.logStep == 1) ?
					(activeVarCount < activeVarCountGrid.maxVal) :
						(activeVarCountGrid.minVal * std::pow(activeVarCountGrid.logStep, i_activeVarCount) < activeVarCountGrid.maxVal));
					i_activeVarCount++)
				{
					if (m_verbose) {
						std::cout << "---\n" << maxDepth << "," << minSampleCountPercentage << "," <<
							maxCategories << "," << activeVarCount << "\n---\n";
					}
					if (m_save_all) {
						std::ofstream rfReportFile;
						rfReportFile.open(m_output_path + "/rf_report.txt", std::ios_base::app); //append file

						rfReportFile << "---\n" << maxDepth << "," << minSampleCountPercentage << "," <<
							maxCategories << "," << activeVarCount << "\n---\n";
					}

					m_rtrees = cv::ml::RTrees::create();
					float val = Train(static_cast<int>(std::floor(maxDepth)),
						minSampleCountPercentage,
						static_cast<int>(std::floor(maxCategories)),
						static_cast<int>(std::floor(activeVarCount)),
						static_cast<int>(std::floor(m_number_of_trees)));

					if (val < minError) {
						minError = val;
						best_rforest = m_rtrees;

						if (m_verbose) {
							std::cout << "Last model is current best with error: " << minError << "\n";
						}
						if (m_save_all) {
							std::ofstream rfReportFile;
							rfReportFile.open(m_output_path + "/rf_report.txt", std::ios_base::app); //append file

							rfReportFile << "Last model is current best with error: " << minError << "\n";
						}
					}

					if (activeVarCountGrid.logStep == 1) {
						activeVarCount++;
					}
					else {
						activeVarCount *= activeVarCountGrid.logStep;
					}
				}

				if (maxCategoriesGrid.logStep == 1) {
					maxCategories++;
				}
				else {
					maxCategories *= maxCategoriesGrid.logStep;
				}
			}

			if (minSampleCountPercentageGrid.logStep == 1) {
				minSampleCountPercentage++;
			}
			else {
				minSampleCountPercentage *= minSampleCountPercentageGrid.logStep;
			}
		}

		if (maxDepthGrid.logStep == 1) {
			maxDepth++;
		}
		else {
			maxDepth *= maxDepthGrid.logStep;
		}
	}

	if (m_verbose) {
		std::cout << "Best model has error: " << minError;
		std::cout << "\n\t MAX DEPTH          = " << best_rforest->getMaxDepth();
		std::cout << "\n\t MIN SAMPLE COUNT % = " << 100 * best_rforest->getMinSampleCount() / m_traindata->getTrainSamples().rows;
		std::cout << "\n\t MAX CATEGORIES     = " << best_rforest->getMaxCategories();
		std::cout << "\n\t ACTIVE VAR COUNT   = " << best_rforest->getActiveVarCount() << "\n";
	}
	if (m_save_all) {
		std::ofstream rfReportFile;
		rfReportFile.open(m_output_path + "/rf_report.txt", std::ios_base::app); //append file

		rfReportFile << "\n[Best model has error: " << minError << "]";
		rfReportFile << "\n\t MAX DEPTH          = " << best_rforest->getMaxDepth();
		rfReportFile << "\n\t MIN SAMPLE COUNT % = " << 100 * best_rforest->getMinSampleCount() / m_traindata->getTrainSamples().rows;
		rfReportFile << "\n\t MAX CATEGORIES     = " << best_rforest->getMaxCategories();
		rfReportFile << "\n\t ACTIVE VAR COUNT   = " << best_rforest->getActiveVarCount() << "\n";
	}
	m_rtrees = best_rforest;

	return minError;
}

float RFSuite::Manager::Train()
{
	return Train(m_max_depth, m_min_sample_count_percentage, m_max_categories, m_active_var_count, m_number_of_trees);
}

float RFSuite::Manager::Train(int maxDepth, double minSampleCountPercentage, int maxCategories, int activeVarCount, int numberOfTrees)
{
	/*if (m_verbose) {
		ConfigParserRF::PrintParseResult(maxDepth, 0, true, minSampleCountPercentage,
			maxCategories, m_active_var_count, m_number_of_trees, m_priors_mat);
	}*/

	if (maxDepth != 0) {
		m_rtrees->setMaxDepth(maxDepth);
	}
	m_rtrees->setMinSampleCount(std::lround((minSampleCountPercentage / 100) * m_traindata->getTrainSamples().rows));
	m_rtrees->setRegressionAccuracy(0); // Not useful for classification
	m_rtrees->setUseSurrogates(false);  // Probably not useful
	m_rtrees->setMaxCategories(maxCategories);
	m_rtrees->setPriors(m_priors_mat);
	m_rtrees->setCalculateVarImportance(true);
	m_rtrees->setActiveVarCount(activeVarCount);
	//m_rtrees->setTermCriteria(cv::TermCriteria(cv::TermCriteria::MAX_ITER + cv::TermCriteria::EPS, 100, FLT_EPSILON));
	m_rtrees->setTermCriteria(cv::TermCriteria(cv::TermCriteria::MAX_ITER, numberOfTrees, 0));

	return train_and_print_errs(m_rtrees, m_traindata);
}

std::shared_ptr<cv::Mat> RFSuite::Manager::Test(cv::Mat &testingMat)
{
	cv::Mat empty; // unnecessary but can't provide argument cv::Mat() in linux
	return Test(testingMat, empty, false);
}

std::shared_ptr<cv::Mat> RFSuite::Manager::Test(cv::Mat &testingMat, cv::Mat &skipZerosMat, bool skipZeros)
{
	int testSize = testingMat.rows;

	std::shared_ptr<cv::Mat> res(new cv::Mat());
	//std::shared_ptr<cv::Mat> res(cv::Mat::zeros(testSize, 1, CV_32S));
	*res = cv::Mat::zeros(testSize, 1, CV_32S);
	//cv::Mat predictLabels;

	int progress = 0;
	int realProgress;
	int val;

	if (m_verbose) {
		std::cout << "RF Manager:\t Testing...0%";
	}
	for (int i = 0; i < testSize; i++)
	{
		if (skipZeros && skipZerosMat.ptr<float>(i)[0] != 0)
		{
			//val = m_rtrees->predict(testingMat.row(i), predictLabels);
			val = std::lround(m_rtrees->predict(testingMat.row(i)));
			res->ptr< int >(i)[0] = val;
		}

		if (m_verbose) {
			realProgress = 100 * i / testSize;

			if (realProgress >= progress + 5) {
				progress = realProgress;
				std::cout << "\r" << "RF Manager:\t Testing..." << progress << "%";
			}
		}
	}
	//m_rtrees->predict(testingMat, predictLabels);

	if (m_verbose) {
		std::cout << "\r" << "RF Manager:\t Testing...finished\n";
	}

	if (m_verbose || m_save_all) {
		cv::Mat variable_importance = m_rtrees->getVarImportance();

		if (m_verbose) {
			std::cout << "RF Manager:\t \tEstimated variable importance:\n";

			for (int i = 0; i < variable_importance.rows; i++) {
				std::cout << "RF Manager:\t \t\tVariable " << i + 1 << ": " << variable_importance.at<float>(i, 0) << "\n";
			}
		}
		if (m_save_all) {
			std::ofstream rfReportFile;
			rfReportFile.open(m_output_path + "/rf_report.txt", std::ios_base::app); //append file

			rfReportFile << "Estimated variable importance:\n";

			for (int i = 0; i < variable_importance.rows; i++) {
				rfReportFile << "\tVariable " << i + 1 << ": " << variable_importance.at<float>(i, 0) << "\n";
			}
			rfReportFile << "\n";
		}
	}

	return res;
}

void RFSuite::Manager::SaveModel(const std::string filename)
{
	m_rtrees->save(filename);
}

float RFSuite::Manager::train_and_print_errs(cv::Ptr<cv::ml::StatModel> model, const cv::Ptr<cv::ml::TrainData>& data)
{
	if (m_verbose) {
		std::cout << "RF Manager:\t Training...";
	}

	bool ok = model->train(data);
	if (!ok)
	{
		if (m_verbose) {
			std::cout << "FAILED\n";
		}

		return -1;
	}
	else
	{
		if (m_verbose) {
			std::cout << "finished\n";

			//std::cout << "# of train: " << data->getTrainSamples().rows << "\n";
		}

		float calcErrorTrain, calcErrorTest;

		calcErrorTest = model->calcError(data, true, cv::noArray());

		if (m_verbose || m_save_all) {
			calcErrorTrain = model->calcError(data, false, cv::noArray());

			if (m_verbose) {
				printf("RF Manager:\t \tTrain error (train part): %f\n", calcErrorTrain);
				printf("RF Manager:\t \tTrain error (test part):  %f\n", calcErrorTest);
			}
			if (m_save_all) {
				std::ofstream rfReportFile;
				rfReportFile.open(m_output_path + "/rf_report.txt", std::ios_base::app); //append file
				rfReportFile << "Train error (train part): " << calcErrorTrain << "\n";
				rfReportFile << "Train error (test part):  " << calcErrorTest << "\n\n";
			}
		}

		return calcErrorTest;
	}
}

void RFSuite::Manager::SetTrainDataFromMats(cv::Mat &trainingMat, cv::Mat &labelsMat)
{
	if (m_verbose) {
		std::cout << "RF Manager:\t Number of samples: " << trainingMat.rows << "\n";
	}
	if (m_save_all) {
		std::ofstream rfReportFile;
		rfReportFile.open(m_output_path + "/rf_report.txt", std::ios_base::app); //append file
		rfReportFile << "Number of samples: " << trainingMat.rows << "\n";
	}
	
	m_traindata = cv::ml::TrainData::create(trainingMat, cv::ml::ROW_SAMPLE, labelsMat);
	m_traindata->setTrainTestSplitRatio(m_training_sample_percentage / 100, true);
}

void RFSuite::Manager::SetPriorsMat(cv::Mat &priorsMat)
{
	m_priors_mat = priorsMat;
}

void RFSuite::Manager::SetOutputPath(std::string path) {
	m_output_path = path;
}
void RFSuite::Manager::SetSaveAll(bool saveAll) {
	m_save_all = saveAll;
}
void RFSuite::Manager::SetVerbose(bool verbose) {
	m_verbose = verbose;
}

void RFSuite::Manager::SetParametersFromConfig(std::string filePath) {
	if (filePath != "") {
		ConfigParserRF::Parse(filePath, m_training_sample_percentage, m_max_depth, m_min_sample_count_percentage,
			m_max_categories, m_active_var_count, m_number_of_trees, m_priors_mat);

		if (m_verbose) {
			ConfigParserRF::PrintParseResult(m_training_sample_percentage, m_max_depth, m_min_sample_count_percentage,
				m_max_categories, m_active_var_count, m_number_of_trees, m_priors_mat);
		}
		if (m_save_all) {
			ConfigParserRF::PrintParseResultToFile(m_output_path + "/rf_report.txt",
				m_training_sample_percentage, m_max_depth, m_min_sample_count_percentage,
				m_max_categories, m_active_var_count, m_number_of_trees, m_priors_mat);
		}

		if (m_traindata != nullptr) {
			m_traindata->setTrainTestSplitRatio(m_training_sample_percentage / 100, true);
		}
	}
}

void RFSuite::Manager::SetTrainingSamplePercentage(double trainingSamplePercentage) {
	m_training_sample_percentage = trainingSamplePercentage;

	if (m_traindata != nullptr) {
		m_traindata->setTrainTestSplitRatio(m_training_sample_percentage / 100, true);
	}
}

void RFSuite::Manager::SetMaxDepth(int maxDepth) {
	m_max_depth = maxDepth;
}
void RFSuite::Manager::SetMinSampleCountPercentage(double minSampleCountPercentage) {
	if (minSampleCountPercentage >= 0 && minSampleCountPercentage <= 100) {
		m_min_sample_count_percentage = minSampleCountPercentage;
	}
	else {
		std::cerr << "Min sample percentage should be in [0,100]\n";
	}
}
void RFSuite::Manager::SetMaxCategories(int maxCategories) {
	m_max_categories = maxCategories;
}
void RFSuite::Manager::SetActiveVarCount(int activeVarCount) {
	m_active_var_count = activeVarCount;
}
void RFSuite::Manager::SetNumberOfTrees(int numberOfTrees) {
	m_number_of_trees = numberOfTrees;
}