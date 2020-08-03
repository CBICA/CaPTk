#include "SvmSuiteManager.h"

void SvmSuite::Manager::GenerateConfigFromBestValues(std::string outputFileName) {
	generateConfig(m_svm_descriptions, outputFileName, m_save_models);
}

void SvmSuite::Manager::GenerateConfigFromBestValues() {
	this->GenerateConfigFromBestValues(m_output_path + "/config.yaml");
}

void SvmSuite::Manager::Train()
{
	// Normalize

	if (m_normalize) {
		startTimer();
		cv::Mat trainSamples = m_traindata->getTrainSamples();
		NormalizeInput(trainSamples);
		stopTimerAndReport("Normalizing SVM input (Training)");
	}

	// Subsampling

	if (m_subsample && m_traindata->getTrainSamples().rows > m_max_samples) {
		m_traindata->setTrainTestSplit(m_max_samples);
		message("Number of samples after subsampling: " + std::to_string(m_traindata->getTrainSamples().rows) + "\n");
	}

	// Training
	
	int counterForFileName = 1;

	for (SvmDescription& svm_desc : m_svm_descriptions)
	{
		if (svm_desc.GetModel() != nullptr) {
			// A pretrained model is used
			continue;
		}
		else if (svm_desc.GetModelPath() != "") {
			// A pretrained model is used
			svm_desc.SetModel(cv::ml::SVM::load(svm_desc.GetModelPath()));
		}
		else {
			startTimer();

			// Train
			auto svm = cv::ml::SVM::create();
			svm->setKernel(svm_desc.GetKernelType());
			svm->setType(svm_desc.GetType());
			svm->setTermCriteria(svm_desc.GetTermCriteria());
			if (svm_desc.GetConsiderWeights()) {
				svm->setClassWeights(m_weights_mat);
			}

			bool res = false;
			try
			{
				message("Training " + svm_desc.GetKernelTypeAsString() + std::string(" kernel..."));

				res = svm->trainAuto(
					m_traindata,
					svm_desc.GetKfold(),
					svm_desc.GetParamGridForParameter(cv::ml::SVM::ParamTypes::C),
					svm_desc.GetParamGridForParameter(cv::ml::SVM::ParamTypes::GAMMA),
					svm_desc.GetParamGridForParameter(cv::ml::SVM::ParamTypes::P),
					svm_desc.GetParamGridForParameter(cv::ml::SVM::ParamTypes::NU),
					svm_desc.GetParamGridForParameter(cv::ml::SVM::ParamTypes::COEF),
					svm_desc.GetParamGridForParameter(cv::ml::SVM::ParamTypes::DEGREE),
					false
				);
				message("", false, true);
			}
			catch (cv::Exception ex) {
				errorOccured(ex.what());
			}

			stopTimerAndReport("Training of " + svm_desc.GetKernelTypeAsString()
				+ std::string(" (n=") + std::to_string(m_traindata->getResponses().rows) + std::string(")")
			);

			// Set and Save results of training
			if (res) {
				// Save svm model
				if (m_save_models) {
					// Create directory for this neighborhood radius if it doesn't exist
					//auto extraDir = m_output_path + "/NRadius" + std::to_string(svm_desc.GetNeighborhoodRadius());
					
					//auto extraDir = m_output_path + "/models";
					auto extraDir = m_output_path;
					
					// if (!cbica::isDir(extraDir)) {
					// 	cbica::createDir(extraDir);
					// }

					// Save model
					std::string savePath = extraDir + "/model" + std::to_string(counterForFileName) + "_" 
						+ svm_desc.GetKernelTypeAsString() + ".xml";
					
					svm->save(savePath);
					svm_desc.SetModelPath(savePath);
				}

				// Set svm to SvmDescription
				svm_desc.SetModel(svm);

				// Set svm parameters to SvmDescription
				svm_desc.SetC(svm->getC());
				svm_desc.SetGamma(svm->getGamma());
				svm_desc.SetP(svm->getP());
				svm_desc.SetNu(svm->getNu());
				svm_desc.SetCoef(svm->getCoef0());
				svm_desc.SetDegree(svm->getDegree());
			}
			counterForFileName++;
		}
	}
}

std::shared_ptr<SvmSuite::Manager::Result> SvmSuite::Manager::Test(cv::Mat &testingMat, bool pseudoProbMapResult) {
	cv::Mat empty; // unnecessary but can't provide argument cv::Mat() in linux
	return Test(testingMat, empty, pseudoProbMapResult, false);
}

std::shared_ptr<SvmSuite::Manager::Result> SvmSuite::Manager::Test(cv::Mat &testingMat, cv::Mat &skipZerosMat, bool pseudoProbMapResult, bool skipZeros)
{
	message("Testing...");


	// Normalize

	if (m_normalize) {
		startTimer();
		NormalizeInput(testingMat);
		stopTimerAndReport("Normalizing SVM input (Testing)");
	}
	
	startTimer();

	// Save importance values for each svm used and the sum
	double importanceSum = 0;
	std::vector< double > importanceValues;

	for (SvmDescription& svm_desc : m_svm_descriptions)
	{
		double importance = svm_desc.GetImportance();
		importanceSum += importance;
		importanceValues.push_back(importance);
	}

	// Initialize output
	std::shared_ptr<Result> res(new Result());

	// pseudoProbMapResult==true -> Pseudoprobability maps for 2-class, false -> n-class classification

	if (pseudoProbMapResult) {
		res->posMat = cv::Mat::zeros(testingMat.rows, 1, CV_32F);
		res->negMat = cv::Mat::zeros(testingMat.rows, 1, CV_32F);

		bool posLabelSet = false, negLabelSet = false;

		cv::Mat predicted(1, 1, CV_32F);
		std::vector<PseudoProbType> outputPseudo(m_svm_descriptions.size());

		PseudoProbType decisionAccu, decision, pos, neg;
		float dist;

		// Iterate through the images and make predictions for each pixel/voxel
		int isvm;

		for (int iTest = 0; iTest < testingMat.rows; iTest++)
		{
			if (skipZeros && (skipZerosMat.row(iTest).at<float>(0, 0) == 0)) {
				continue;
			}

			isvm = 0;
			// For each svm
			for (SvmDescription& svm_desc : m_svm_descriptions)
			{
				// Predict
				if (svm_desc.GetModel() == nullptr) {
					errorOccured("Trying to make predictions on untrained SVM");
					return res;
				}
				else {
					svm_desc.GetModel()->predict(testingMat.row(iTest), predicted, pseudoProbMapResult);

					dist = predicted.at<float>(0, 0);
					if ((!posLabelSet) && (dist > 0)) {
						// The label hasn't been set for pos image
						svm_desc.GetModel()->predict(testingMat.row(iTest), predicted, false);
						res->posLabel = std::lround(predicted.at<float>(0, 0));
						posLabelSet = true;
					}
					else if ((!negLabelSet) && (dist < 0)) {
						// The label hasn't been set for neg image
						svm_desc.GetModel()->predict(testingMat.row(iTest), predicted, false);
						res->negLabel = std::lround(predicted.at<float>(0, 0));
						negLabelSet = true;
					}

					// The value is distance to the hyperplane (signed positive for one label, negative for the other)
					outputPseudo[isvm++] = dist;
				}

				decisionAccu = 0;
				decision = 0;

				for (int i = 0; i < outputPseudo.size(); i++) {
					decisionAccu += static_cast<float>(outputPseudo[i] * importanceValues[i]);
				}

				if (importanceSum != 0) {
					// Average of the each different svm's decision
					decision = static_cast<float>(decisionAccu / importanceSum);
				}

				// By using the sigmoid function: f(x) = 1 / (1+e^(-x)) the values are normalized to be doubles between [0,1]
				pos = static_cast<float>(1.0 / (1.0 + std::exp(-decision))); // Accounts for the distances for one label  (the one with pos distances)
				neg = static_cast<float>(1.0 / (1.0 + std::exp(decision))); // Accounts for the distances for the other label (the one with neg distances)

				// Set values to output
				res->posMat.ptr<PseudoProbType>(iTest)[0] = pos;
				res->negMat.ptr<PseudoProbType>(iTest)[0] = neg;
			}
		}
	}
	else {
		//Create output image with the same dimensions as the input images for this subject
		res->labelsMat = cv::Mat::zeros(testingMat.rows, 1, CV_32S);

		int totalThreadsNumber = (m_number_of_threads > testingMat.rows) ? 1 : m_number_of_threads;
		int counterForThreadsVec = 0;
		std::vector<std::thread> threads(totalThreadsNumber);

		for (int iStart = 0; iStart < totalThreadsNumber; iStart++)
		{
			threads[counterForThreadsVec++] = std::thread(&SvmSuite::Manager::testingLabelsThreadJob, this,
				std::ref(testingMat), iStart, totalThreadsNumber,
				std::ref(res->labelsMat), skipZeros, std::ref(skipZerosMat), std::ref(importanceValues)
			);
		}

		for (int i = 0; i < totalThreadsNumber; i++) {
			threads[i].join();
		}
	}

	stopTimerAndReport("SVM predictions");

	message("", false, true);

	return res;
}

void SvmSuite::Manager::AddSvmDescriptionToList(SvmDescription svmDesc) {
	m_svm_descriptions.push_back(svmDesc);
}

void SvmSuite::Manager::AddSvmDescriptions(std::vector< SvmDescription > svmDescs) {
	for (size_t i = 0; i < svmDescs.size(); i++) {
		m_svm_descriptions.push_back(svmDescs[i]);
	}
}

void SvmSuite::Manager::AddPretrainedModel(std::string pretrainedModelPath, int neighborhoodRadius) {
	SvmDescription svm_desc;
	svm_desc.SetModelPath(pretrainedModelPath);
	svm_desc.SetNeighborhoodRadius(neighborhoodRadius);
	m_svm_descriptions.push_back(svm_desc);
}

void SvmSuite::Manager::AddSvmsFromConfig(std::string configPath) {
	this->AddSvmDescriptions(getSvmDescriptionsFromConfig(configPath));
}

void SvmSuite::Manager::SetTrainData(cv::Mat &trainingMat, cv::Mat &labelsMat, cv::Mat &weightsMat/*, cv::Mat sampleIdx*/)
{
	message("Number of samples:  " + std::to_string(trainingMat.rows) + std::string("\n"));
	message("Number of features: " + std::to_string(trainingMat.cols) + std::string("\n"));

	m_traindata = cv::ml::TrainData::create(trainingMat, cv::ml::ROW_SAMPLE, labelsMat/*, cv::Mat(), sampleIdx*/);

	m_weights_mat = weightsMat;

	for (int i = 0; i < m_traindata->getClassLabels().rows; i++) {
		m_different_labels.insert(m_traindata->getClassLabels().at<LabelsType>(i, 0));
	}
}

void SvmSuite::Manager::SetVerbose(bool verbose) {
	m_verbose = verbose;
}

void SvmSuite::Manager::SetOutputPath(std::string path) {
	m_output_path = path;
}

void SvmSuite::Manager::SetSavingModelsEnabled(bool modelsEnabled) {
	m_save_models = modelsEnabled;
}

void SvmSuite::Manager::SetTimerEnabled(bool timerEnabled) {
	m_timer_enabled = timerEnabled;
}

void SvmSuite::Manager::SetNumberOfThreads(int numberOfThreads) {
	m_number_of_threads = numberOfThreads;
}

void SvmSuite::Manager::SetSubsampling(bool subsample, int maxSamples) {
	m_subsample = subsample;
	m_max_samples = maxSamples;
}

void SvmSuite::Manager::SetInputNormalization(bool normalize) {
	m_normalize = normalize;
}

void SvmSuite::Manager::testingLabelsThreadJob(cv::Mat &testingMat, int iStart, int interval, cv::Mat &resultLabelsMat,
	bool skipZeros, cv::Mat &skipZerosMat, std::vector< double > &importanceValues)
{
	std::map< LabelsType, double > decisionImportanceValues;

	cv::Mat predicted(1, 1, CV_32F);
	std::vector<LabelsType> outputLabels(m_svm_descriptions.size());
	int isvm;

	// Iterate through the images and make predictions for each pixel/voxel
	for (int iTest = iStart; iTest < testingMat.rows; iTest += interval)
	{
		if (skipZeros && (skipZerosMat.row(iTest).at<float>(0, 0) == 0)) {
			continue;
		}

		isvm = 0;

		// For each svm
		for (SvmDescription& svm_desc : m_svm_descriptions)
		{
			// Predict
			if (svm_desc.GetModel() == nullptr) {
				errorOccured("Trying to make predictions on untrained SVM");
				//return res;
			}
			else {
				svm_desc.GetModel()->predict(testingMat.row(iTest), predicted, false);

				// The value is the predicted label
				outputLabels[isvm++] = std::lround(predicted.at<float>(0, 0));
				//outputLabels[isvm++] = std::round(predicted.at<float>(0, 0));
			}
		}

		// The values are the predicted labels
		// In case of different svms predicting a different label
		// then the prediction with the most accumulated importance will be used
		// In case of a draw the decision from the first in order svm is used (draws should be avoided in the configuration)
		for (LabelsType label : m_different_labels) {
			// This is done for time optimization (not creating a new object each time)
			decisionImportanceValues[label] = 0;
		}

		// For each different label in the predictions, find the sum of the importance values of the svms that predicted it
		for (int i = 0; i < outputLabels.size(); i++) {
			//if (decisionImportanceValues.find(outputLabels[i]) == decisionImportanceValues.end())
			//{
			//	decisionImportanceValues[outputLabels[i]] = importanceValues[i];
			//}
			//else {
			decisionImportanceValues[outputLabels[i]] += importanceValues[i];
			//}
		}

		// Find the label with the most importance value
		//std::vector< int > keys = getMapKeyset(decisionImportanceValues);
		int    decision = 0;
		double bestDecisionImportance = 0;

		for (int decisionCandidate : m_different_labels)
		{
			if (decisionImportanceValues[decisionCandidate] > bestDecisionImportance)
			{
				decision = decisionCandidate;
				bestDecisionImportance = decisionImportanceValues[decisionCandidate];
			}
		}

		// Set value to output
		resultLabelsMat.ptr<LabelsType>(iTest)[0] = decision;
	}
}

void SvmSuite::Manager::startTimer() {
	if (m_timer_enabled) {
		m_timer.Reset();
	}
}

void SvmSuite::Manager::stopTimerAndReport(std::string desc) {
	if (m_timer_enabled) {
		float diff = m_timer.Diff();

		std::ofstream timerFile;
		timerFile.open(m_output_path + "/time_report.txt", std::ios_base::app); //append file
		timerFile << desc << ": " << diff << "s\n";
	}
}

void SvmSuite::Manager::message(std::string message, bool overdraw, bool finished, int progress)
{
	if (m_verbose) {
		if (overdraw) {
			std::cout << "\r";
		}

		if (message != "") {
			std::cout << "SVM Manager:\t" << message;

			if (progress != -1) {
				std::cout << " [" << progress << "%]";
			}
		}

		if (finished) {
			std::cout << "finished\n";
		}
	}
}

void SvmSuite::Manager::errorOccured(std::string msg) {
	std::cerr << "SVM Manager error: " << msg << std::endl;
}
