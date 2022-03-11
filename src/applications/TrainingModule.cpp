/**
\file  TrainingModule.cpp

\brief Source file containing the implementation for the TrainingModule class.
This class is a refactored framework based on the original TrainingModule by Saima Rathore.
Author: Alexander Getka
Library Dependencies: ITK 4.7+ <br>

https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2020 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/captk/license.html

*/

#include "TrainingModule.h"
//#include "TrainingModuleThreadWorker.h"
#include "opencv2/core/core.hpp"
#include "opencv2/ml.hpp"
#include <cmath>
#include "CaPTkUtils.h"

// Machine Learning / Binary Classification
#include "OpenCVSVMClassifierStrategy.h"
#include "OpenCVRandomForestStrategy.h"
#include "OpenCVSGDSVMStrategy.h"
#include "OpenCVBoostedTreesStrategy.h"

// Feature Selection
#include "BackwardFeatureSelectionStrategy.h"
#include "EffectSizeFeatureSelectionStrategy.h"
#include "ForwardFeatureSelectionStrategy.h"
#include "RandomForestFeatureSelectionStrategy.h"
#include "ReliefFFeatureSelectionStrategy.h"

//etc
#include <fstream>
#include <iostream>


TrainingModule::TrainingModule()
{
	//isCurrentlyRunning = false;
	//lastResult.success = false;
}

TrainingModule::~TrainingModule()
{
	// See note in RunThread about modifying the underlying TrainingModule code,
	//and why terminate is needed right now
	//m_workerThread->terminate();
	//m_workerThread->requestInterruption(); // Advise code to 
	//m_workerThread->quit(); // Instruct worker thread to exit
	//m_workerThread->wait(); // Ensure worker thread is actually terminated before continuing
}

const TrainingModuleResult TrainingModule::Run(const TrainingModuleParameters& params)
{
	std::cout << "Entered the Training Module." << std::endl;
	
	TrainingModuleResult result;
	bool success = false;
	try
	{
		if (params.configurationType == CAPTK::ClassificationConfigurationType::CONF_TYPE_SPLIT_TEST)
		{
			success = RunSplitTesting(params);
		}
		else if (params.configurationType == CAPTK::ClassificationConfigurationType::CONF_TYPE_SPLIT_TRAIN)
		{
			success = RunSplitTraining(params);
		}
		else if (params.configurationType == CAPTK::ClassificationConfigurationType::CONF_TYPE_KFOLD_CV)
		{
			success = RunKFoldCrossValidation(params);
		}
		else 
		{
			success = false;
			auto msg = "Requested an unsupported TrainingModule configuration. Supported configurations are Split Train, Split Test, K-fold CrossValidation.";
			logger.WriteError(msg);
			result.message = msg;
		}
	}
	catch (const std::exception& exc)
	{
		success = false;
		logger.WriteError("Exception caught while running the training module. Error information: " + std::string(exc.what()) );
		logger.WriteError("Please double check the validity of your data and ensure that multiple examples exist of each class.");
		result.success = false;
		result.message = exc.what();
		return result;
	}
	result.success = success;

	// TODO: Return a result struct for passing to CLI/GUI? Needs handling there.
	if (success)
	{
		const auto successMsg = "Training module task completed successfully.";
		std::cout << successMsg << std::endl;
		logger.Write(successMsg);
		result.message = successMsg;
	}
	else
	{
		const auto errMsg = "The training module task failed. See logs for more information.";
		std::cerr << errMsg << std::endl;
		logger.WriteError("The training module task failed.");
		result.message = errMsg;
	}
	return result;

}

bool TrainingModule::RunSplitTraining(const TrainingModuleParameters& params)
{
	//// Load training data into memory 
	std::cout << "Running training." << std::endl;
	VariableSizeMatrixType FeaturesOfAllSubjects;
	VariableLengthVectorType LabelsOfAllSubjects;
	std::vector<std::string> featureRowHeaders;
	std::vector<std::string> featureColHeaders;
	std::vector<std::string> labelRowHeaders;
	std::vector<std::string> labelColHeaders;
	typedef vnl_matrix<double> MatrixType;
	MatrixType dataMatrix;

	auto featureLoadingSucceeded = GetFeatureDataFromFile(params.inputFeaturesFile, FeaturesOfAllSubjects, featureRowHeaders, featureColHeaders);
	auto labelLoadingSucceeded = GetLabelsFromFile(params.inputLabelsFile, LabelsOfAllSubjects, labelRowHeaders, labelColHeaders);

	if (!std::get<0>(featureLoadingSucceeded) || !std::get<0>(labelLoadingSucceeded))
	{
		std::cerr << "Error when loading training data. See logs for more information." << std::endl;
	    return false;
	}
	std::cout << "Finished loading training data." << std::endl;

	//// Generate scaled features, save means/std-devs to disk
	//// NOTE: We are scaling (normalizing) according to the z-scores of entire set of samples here.
	// ALSO NOTE: This is *feature-wise* z-scoring.
	// This is a separate process from z-score normalization of image intensities, which should be done before extraction of features.
	
	// Features are z-scored by their individual means/std-devs on all sample data loaded here.
	// But the z-scores for the whole dataset will not be the same as the z-scores for a smaller sample.
	// We need to determine, if we split this dataset into training and testing/holdout sets (such as for CV),
	// Do we determine means/std devs from the training set, and then scale ALL data according to that?
	// Or do we perform z-scoring *on each split sample set*, against their own means/std-devs?
	// If so, we need to do this scaling inside the code that performs this split.
	FeatureScalingClass mFeatureScaler;
	VariableSizeMatrixType scaledFeatureSet;
	VariableLengthVectorType meanVector;
	VariableLengthVectorType stdVector;
	// Performs z-scoring
	mFeatureScaler.ScaleGivenTrainingFeatures(FeaturesOfAllSubjects, scaledFeatureSet, meanVector, stdVector);
	// Remove nans before saving
	RemoveNans(scaledFeatureSet);
	RemoveNans(meanVector);
	RemoveNans(stdVector);

	// Write scaling information to disk
	WriteCSVFiles(scaledFeatureSet, params.outputDirectory + "/scaled-feature-set.csv");
	WriteCSVFiles(meanVector, params.outputDirectory + "/zscore-mean.csv");
	WriteCSVFiles(stdVector, params.outputDirectory + "/zscore-std.csv");
	std::cout << "Scaling parameters written." << std::endl;


	// Pointers to strategy interfaces -- logic for each approach is implemented in its own class.
	// Any class that publically implements that interface can be used here, allowing easy extension.
	std::shared_ptr<IFeatureSelectionStrategy> featureSelector;
	std::shared_ptr<IClassifierStrategy> classifier;
	
	// Set Feature Selection strategy
	if (params.featureSelectionType == CAPTK::FeatureSelectionType::FS_TYPE_ES)
	{
		featureSelector = std::make_shared<EffectSizeFeatureSelectionStrategy>();
	}
	else if (params.featureSelectionType == CAPTK::FeatureSelectionType::FS_TYPE_FFS)
	{
		featureSelector = std::make_shared<ForwardFeatureSelectionStrategy>();
	}
	else if (params.featureSelectionType == CAPTK::FeatureSelectionType::FS_TYPE_BACKWARDS)
	{
		featureSelector = std::make_shared<BackwardFeatureSelectionStrategy>();
	}
	else if (params.featureSelectionType == CAPTK::FeatureSelectionType::FS_TYPE_RANDOMFOREST)
	{
		featureSelector = std::make_shared<RandomForestFeatureSelectionStrategy>();
	}
	else if (params.featureSelectionType == CAPTK::FeatureSelectionType::FS_TYPE_RELIEF_F)
	{
		featureSelector = std::make_shared<ReliefFFeatureSelectionStrategy>();
	}
	// TODO: Enable Optimal FS when it's fixed.
	/*
	else if (params.featureSelectionType == CAPTK::FeatureSelectionType::FS_TYPE_OPTIMIZE_ANY)
	{
		featureSelector = std::make_shared<OptimalFeatureSelectionStrategy>();
	}
	*/
	else
	{
		// If we reach here, the user interface is broken. Only supported methods should be exposed to the GUI/CLI.
		throw std::logic_error("Requested an unsupported Feature Selection strategy.");
	}


	if (params.classifierType == CAPTK::ClassifierType::CLASS_TYPE_SVM_LINEAR || 
		params.classifierType == CAPTK::ClassifierType::CLASS_TYPE_SVM_RBF ||
		params.classifierType == CAPTK::ClassifierType::CLASS_TYPE_SVM_POLYNOMIAL ||
		params.classifierType == CAPTK::ClassifierType::CLASS_TYPE_SVM_SIGMOID || 
		params.classifierType == CAPTK::ClassifierType::CLASS_TYPE_SVM_CHISQUARED ||
		params.classifierType == CAPTK::ClassifierType::CLASS_TYPE_SVM_INTERSECTION)
	{
		// OpenCVSVMClassifierStrategy covers OpenCV SVM kernel types
		// see https://docs.opencv.org/3.4/d1/d2d/classcv_1_1ml_1_1SVM.html for information
		classifier = std::make_shared<OpenCVSVMClassifierStrategy>();
	}
	else if (params.classifierType == CAPTK::ClassifierType::CLASS_TYPE_RANDOMFOREST)
	{
		classifier = std::make_shared<OpenCVRandomForestStrategy>();
	}
	else if (params.classifierType == CAPTK::ClassifierType::CLASS_TYPE_SGD_SVM)
	{
		classifier = std::make_shared<OpenCVSGDSVMStrategy>();
	}
	else if (params.classifierType == CAPTK::ClassifierType::CLASS_TYPE_BOOSTEDTREES)
	{
		classifier = std::make_shared<OpenCVBoostedTreesStrategy>();
	}
	/*
	else if (params.classifierType == CAPTK::ClassifierType::CLASS_TYPE_OPTIMIZE_ANY)
	{
		classifier = std::make_shared<SelectOptimalClassifierStrategy>();
	} */
	else
	{
		// If we reach here, the user interface is broken. Only supported methods should be exposed to the GUI/CLI.
		throw std::logic_error("Requested an unsupported classifier.");
	}
	
	// copy current execution parameters to the classifier and feature selector
	// Note that the feature selection strategy may choose to modify classifier parameters,
	// But it should generally re-set them when it is done.
	classifier->SetParameters(params);
	featureSelector->SetParameters(params);
	

	// Convert labels to VectorDouble first...
	VectorDouble labels;
	labels.resize(LabelsOfAllSubjects.Size());
	for (int i = 0; i < LabelsOfAllSubjects.Size(); i++)
	{
		labels[i] = LabelsOfAllSubjects[i];
	}

	std::vector<int> selectedFeatureIndices;

	// Perform training with the selected strategies using the previously scaled features.
	selectedFeatureIndices = featureSelector->PerformFeatureSelectionBasedTraining(scaledFeatureSet, labels, classifier.get());
	WriteCSVFiles(selectedFeatureIndices, params.outputDirectory + "/selected-features.csv");
	if (!featureColHeaders.empty())
	{
		std::ofstream myfile;
		myfile.open(params.outputDirectory + "/selected-feature-names.csv");
		std::vector<std::string> selectedFeatureNames;
		for (int i = 0; i < selectedFeatureIndices.size(); i++)
		{
			std::string currentFeatureName = featureColHeaders[selectedFeatureIndices[i]];
			if (i == 0)
			{
				myfile << currentFeatureName;
			}
			else
			{
				myfile << "," << currentFeatureName;
			}
		}
		myfile.close();
	}

	VariableSizeMatrixType selectedScaledFeatureData;
	selectedScaledFeatureData.SetSize(labels.size(), selectedFeatureIndices.size()); // samples x selected features
	// copy data from selected features into final dataset
	for (int sample_idx = 0; sample_idx < labels.size(); sample_idx++)
	{
		for (int feature_idx = 0; feature_idx < selectedFeatureIndices.size(); feature_idx++)
			selectedScaledFeatureData(sample_idx, feature_idx) = scaledFeatureSet(sample_idx, selectedFeatureIndices[feature_idx]);
	}


	classifier->train(selectedScaledFeatureData, labels);
	classifier->saveModel(); // defaults to params.outputDirectory, no need to specify

	// TODO: Output trained model metadata as YML, output version YML (can also do in saveModel?)
	std::cout << "Misc. Information about the final model:" << std::endl;
	std::cout << classifier->getModelInfoAsString() << std::endl;
	auto cvResults = classifier->CalculatePerformanceMeasuresAgainstLabels(selectedScaledFeatureData, labels);
	std::cout << "Overall accuracy: " + std::to_string(std::get<0>(cvResults)) << std::endl;
	std::cout << "Sensitivity: " + std::to_string(std::get<1>(cvResults)) << std::endl;
	std::cout << "Specificity: " + std::to_string(std::get<2>(cvResults)) << std::endl;
	std::cout << "Balanced accuracy: " + std::to_string(std::get<3>(cvResults)) << std::endl;
	std::cout << "Finished training. Results placed in " + params.outputDirectory + "." << std::endl;

	return true;
}

bool TrainingModule::RunSplitTesting(const TrainingModuleParameters& params)
{
	// This can probably be broken up further.

	// Read in feature means and std-devs for feature scaling
	CSVFileReaderType::Pointer reader = CSVFileReaderType::New();

	VariableLengthVectorType featureMeans;
	VariableLengthVectorType featureStdDeviations;
	std::vector<double> selectedFeatureIndices;
	
	try
	{
		MatrixType meanMatrix;
		reader->SetFileName(params.modelDirectory + "/zscore-mean.csv");
		reader->SetFieldDelimiterCharacter(',');
		reader->HasColumnHeadersOff();
		reader->HasRowHeadersOff();
		reader->Parse();
		meanMatrix = reader->GetArray2DDataObject()->GetMatrix();

		featureMeans.SetSize(meanMatrix.size());
		for (unsigned int i = 0; i < meanMatrix.size(); i++)
			featureMeans[i] = meanMatrix(0, i);
	}
	catch (const std::exception& e1)
	{
		throw std::runtime_error("Error in reading the file: " + params.modelDirectory + "/zscore-mean.csv. Error code : " + std::string(e1.what()));
	}

	try
	{
		MatrixType stdMatrix;
		reader->SetFileName(params.modelDirectory + "/zscore-std.csv");
		reader->SetFieldDelimiterCharacter(',');
		reader->HasColumnHeadersOff();
		reader->HasRowHeadersOff();
		reader->Parse();
		stdMatrix = reader->GetArray2DDataObject()->GetMatrix();

		featureStdDeviations.SetSize(stdMatrix.size());
		for (unsigned int i = 0; i < stdMatrix.size(); i++)
			featureStdDeviations[i] = stdMatrix(0, i);
	}
	catch (const std::exception& e1)
	{
		throw std::runtime_error("Error in reading the file: " + params.modelDirectory + "/zscore-std.csv. Error code : " + std::string(e1.what()));
	}

	try
	{
		MatrixType selectionMatrix;
		reader->SetFileName(params.modelDirectory + "/selected-features.csv");
		reader->SetFieldDelimiterCharacter(',');
		reader->HasColumnHeadersOff();
		reader->HasRowHeadersOff();
		reader->Parse();
		selectionMatrix = reader->GetArray2DDataObject()->GetMatrix();

		for (unsigned int i = 0; i < selectionMatrix.cols(); i++)
			selectedFeatureIndices.push_back(selectionMatrix(0, i));
	}
	catch (const std::exception& e1)
	{
		throw std::runtime_error("Error in reading the file: " + params.modelDirectory + "/selected-features.csv. Error code : " + std::string(e1.what()));
	}

	VariableSizeMatrixType allFeaturesMatrix;
	std::vector<std::string> featureRowHeaders;
	std::vector<std::string> featureColHeaders;
	auto featuresLoadingSucceeded = GetFeatureDataFromFile(params.inputFeaturesFile, allFeaturesMatrix, featureRowHeaders, featureColHeaders);

	std::shared_ptr<IClassifierStrategy> predictorPtr; // to be filled with the appropriate classifier for model handling

	// Need to determine what strategy to use for model load/prediction.
	// TBD: Do this based on user specification/parameter or try to auto-detect?
	// For right now, we'll just check for presence of a predetermined filename and use that.
	// In the future, could use some kind of mapping to iterate over, or use the "Factory Method" pattern
	// Can make the filename to look for a static variable of the classifier strategy class in the header
	// then create a mapping between the classifier enum and that string
	int modelsFound = 0;
	QString modelQuery = QString::fromStdString(params.modelDirectory + "/OpenCV_SVM_Model.xml");
	if (QFileInfo::exists(modelQuery))
	{
		modelsFound++;
		predictorPtr = std::make_shared<OpenCVSVMClassifierStrategy>();
	}
	modelQuery = QString::fromStdString(params.modelDirectory + "/OpenCV_RF_Model.xml");
	if (QFileInfo::exists(modelQuery))
	{
		modelsFound++;
		predictorPtr = std::make_shared<OpenCVRandomForestStrategy>();
	}
	modelQuery = QString::fromStdString(params.modelDirectory + "/OpenCV_SVMSGD_Model.xml");
	if (QFileInfo::exists(modelQuery))
	{
		modelsFound++;
		predictorPtr = std::make_shared<OpenCVSGDSVMStrategy>();
	}
	modelQuery = QString::fromStdString(params.modelDirectory + "/OpenCV_BoostedTrees_Model.xml");
	if (QFileInfo::exists(modelQuery))
	{
		modelsFound++;
		predictorPtr = std::make_shared<OpenCVBoostedTreesStrategy>();
	}

	// add additional model-type handling (beyond SVMs and Random Forest) here...

	if (modelsFound == 0) // no models present
	{
		std::string errMsg = "Model not found: model directory ("
			+ params.modelDirectory + ") does not contain a detectable model file.";

		std::cerr << errMsg << std::endl;
		logger.WriteError(errMsg);
		return false;
	}

	// error out if too many models are present
	if (modelsFound > 1) {
		std::string errMsg = "Model type ambiguous: model directory ("
			+ params.modelDirectory + ") contains multiple conflicting model files.";

		std::cerr << errMsg << std::endl;
		logger.WriteError(errMsg);
		return false;
	}

	predictorPtr->SetParameters(params);
	// Now that we know what/how to load...
	predictorPtr->loadModel();


	// Set up storage for results
	VectorDouble predictedDistances;
	VectorDouble predictedLabels;

	// Scale features (z-scoring, given feature means/std-devs)
	FeatureScalingClass mFeatureScaler;
	VariableSizeMatrixType scaledFeaturesMatrix;
	// Perform z-scoring given the information we loaded earlier
	scaledFeaturesMatrix = mFeatureScaler.ScaleGivenTestingFeatures(allFeaturesMatrix, featureMeans, featureStdDeviations);
	RemoveNans(scaledFeaturesMatrix);

	// Select features from loaded selection indices
	VariableSizeMatrixType finalFeatures;
	finalFeatures.SetSize(scaledFeaturesMatrix.Rows(), selectedFeatureIndices.size());
	for (unsigned int j = 0; j < scaledFeaturesMatrix.Rows(); j++)
		for (unsigned int k = 0; k < finalFeatures.Cols(); k++)
			finalFeatures(j, k) = scaledFeaturesMatrix(j, selectedFeatureIndices[k]);

	// Run prediction for both classes and distances
	predictedLabels = predictorPtr->predict(finalFeatures);
	predictedDistances = predictorPtr->predict(finalFeatures, true);
	std::vector<std::string> emptyHeaders;
	WriteLabelFilesWithHeaders(params.outputDirectory + "/predicted-distances.csv", predictedDistances, featureRowHeaders, emptyHeaders);
	WriteLabelFilesWithHeaders(params.outputDirectory + "/predicted-labels.csv", predictedLabels, featureRowHeaders, emptyHeaders);

	if (params.testPredictionsAgainstProvidedLabels)
	{
		VariableLengthVectorType actualLabels;
		std::vector<std::string> labelRowHeaders;
		std::vector<std::string> labelColHeaders;
		auto succeeded = GetLabelsFromFile(params.inputLabelsFile, actualLabels, labelRowHeaders, labelColHeaders);
		VariableLengthVectorType predictedLabelsAsItkVector;
		predictedLabelsAsItkVector.SetSize(predictedLabels.size());
		for (unsigned int i = 0; i < predictedLabels.size(); i++)
		{
			predictedLabelsAsItkVector[i] = predictedLabels[i];
		}


		// TODO: also consider validating with AUC (needs thresholding effect etc). In general, add performance metrics here.
		double balancedAccuracy = GetBinaryClassificationBalancedAccuracy(predictedLabelsAsItkVector, actualLabels);
		std::ofstream performancefile;
		performancefile.open(params.outputDirectory + "/performance.txt");
		performancefile << "Balanced accuracy: " << std::to_string(balancedAccuracy) << std::endl;
		performancefile.close();
		

	}

	return true;
}


bool TrainingModule::RunKFoldCrossValidation(const TrainingModuleParameters& params)
{
	// Train, but using folds and testing afterwards
	// TODO: implement this
	std::cout << "Running training module in cross-validation mode." << std::endl;
	std::cout << "Loading data." << std::endl;
	VariableSizeMatrixType FeaturesOfAllSubjects;
	VariableLengthVectorType LabelsOfAllSubjects;
	std::vector<std::string> featureRowHeaders;
	std::vector<std::string> featureColHeaders;
	std::vector<std::string> labelRowHeaders;
	std::vector<std::string> labelColHeaders;
	typedef vnl_matrix<double> MatrixType;
	MatrixType dataMatrix;

	auto featureLoadingSucceeded = GetFeatureDataFromFile(params.inputFeaturesFile, FeaturesOfAllSubjects, featureRowHeaders, featureColHeaders);
	auto labelLoadingSucceeded = GetLabelsFromFile(params.inputLabelsFile, LabelsOfAllSubjects, labelRowHeaders, labelColHeaders);

	auto folds = params.folds;
	if (folds <= LabelsOfAllSubjects.Size() * 2 )
	{
		throw std::runtime_error("Fold size too great for the amount of samples -- ensure your data has enough samples, and contains examples of each class.");
	}
	std::vector<int> sampleIndices;
	for (int i = 0; i < LabelsOfAllSubjects.Size(); i++)
	{
		// Just a 1 to N range that can be shuffled for indexing into things randomly later
		sampleIndices.push_back(i);
	}
	std::random_shuffle(sampleIndices.begin(), sampleIndices.end());
	int testGroupStartIndex = 0;
	int testGroupSize = sampleIndices.size() / (double)folds;
	int groupRemainder = sampleIndices.size() % (int)folds;
	int trainGroupSize = sampleIndices.size() - testGroupSize;
	for (int currentFold = 0; currentFold < folds; currentFold++)
	{
		VariableSizeMatrixType trainingData;
		VariableSizeMatrixType testingData;
		VectorDouble trainingLabels;
		VectorDouble testingLabels;
		std::vector<std::string> testingSubjectHeaders;
		std::vector<std::string> trainingSubjectHeaders;

		trainingData.SetSize(trainGroupSize, FeaturesOfAllSubjects.Cols());
		testingData.SetSize(testGroupSize, FeaturesOfAllSubjects.Cols());

		for (int i = 0; i < testGroupSize; i++) // copy testing data samples
		{
			int testGroupIndex = i + testGroupStartIndex;
			testingLabels.push_back(LabelsOfAllSubjects[sampleIndices[testGroupIndex]]);
			testingSubjectHeaders.push_back(featureRowHeaders[sampleIndices[testGroupIndex]]);

			for (int j = 0; j < FeaturesOfAllSubjects.Cols(); j++)
			{
				testingData(i, j) = FeaturesOfAllSubjects(sampleIndices[testGroupIndex], j);
			}
		}
		for (int i = 0; i < testGroupStartIndex; i++) // copy first part of training data samples
		{
			trainingLabels.push_back(LabelsOfAllSubjects[i]);
			trainingSubjectHeaders.push_back(featureRowHeaders[sampleIndices[i]]);

			for (int j = 0; j < FeaturesOfAllSubjects.Cols(); j++)
			{
				trainingData(i, j) = FeaturesOfAllSubjects(sampleIndices[i], j);
			}
		}
		for (int i = testGroupStartIndex + testGroupSize; i < sampleIndices.size(); i++) // copy second part of training data samples
		{
			trainingLabels.push_back(LabelsOfAllSubjects[i]);
			trainingSubjectHeaders.push_back(featureRowHeaders[sampleIndices[i]]);

			for (int j = 0; j < FeaturesOfAllSubjects.Cols(); j++)
			{
				auto value = FeaturesOfAllSubjects(sampleIndices[i], j);
				trainingData(i - testGroupSize, j) = value;
			}
		}

		auto outputFoldDir = params.outputDirectory + "/Fold" + std::to_string(currentFold);
		if (!cbica::makeDirectory(outputFoldDir))
		{
			throw std::runtime_error("Couldn't create the fold subdirectory during cross-validation mode.");
		}

		WriteNumericFilesWithHeaders(outputFoldDir + "/raw-training-features.csv", trainingData, trainingSubjectHeaders, featureColHeaders);
		WriteLabelFilesWithHeaders(outputFoldDir + "/training-labels.csv", trainingLabels, trainingSubjectHeaders, labelColHeaders);
		WriteNumericFilesWithHeaders(outputFoldDir + "/raw-testing-features.csv", testingData, testingSubjectHeaders, featureColHeaders);
		WriteLabelFilesWithHeaders(outputFoldDir + "/ground-truth-labels.csv", testingLabels, testingSubjectHeaders, labelColHeaders);

		// Run training on this fold's training data
		auto thisFoldParams = params;
		thisFoldParams.outputDirectory = outputFoldDir;
		thisFoldParams.configurationType = CAPTK::ClassificationConfigurationType::CONF_TYPE_SPLIT_TRAIN;
		thisFoldParams.inputFeaturesFile = outputFoldDir + "/raw-training-features.csv";
		thisFoldParams.inputLabelsFile = outputFoldDir + "/training-labels.csv";
		RunSplitTraining(thisFoldParams);

		// Run testing/inference on this fold's test featureset using the trained model, get balanced accuracy vs testing labels
		thisFoldParams.modelDirectory = outputFoldDir;
		thisFoldParams.inputFeaturesFile = outputFoldDir + "/raw-testing-features.csv";
		thisFoldParams.configurationType = CAPTK::ClassificationConfigurationType::CONF_TYPE_SPLIT_TEST;
		thisFoldParams.testPredictionsAgainstProvidedLabels = true;
		thisFoldParams.inputLabelsFile = outputFoldDir + "/ground-truth-labels.csv";
		RunSplitTesting(thisFoldParams);

		std::cout << "Finished processing fold " + std::to_string(currentFold) << std::endl;
		testGroupStartIndex += testGroupSize; // Shift indexing up for next fold

	}

	std::cout << "Finished cross-validation. Please check the output directory for each fold for performance scores." << std::endl;

	return true;
}


std::tuple<bool, bool, bool> TrainingModule::GetFeatureDataFromFile(std::string featuresFilename, VariableSizeMatrixType& featuresMatrix, std::vector<std::string>& rowHeaders, std::vector<std::string>& colHeaders)
{

	// Return value :
	// first is true if successful overall
	// second is true if row headers were present (also populates if true)
	// third is true if column headers were present (also populates if true)
	std::tuple<bool, bool, bool> result = std::make_tuple(false, false, false);
	bool hasRowHeaders = false;
	bool hasColumnHeaders = false;

	MatrixType dataMatrix;
	try
	{
		CSVFileReaderType::Pointer readerMean = CSVFileReaderType::New();
		readerMean->SetFileName(featuresFilename);
		readerMean->SetFieldDelimiterCharacter(',');
		readerMean->HasColumnHeadersOff();
		readerMean->HasRowHeadersOff();
		readerMean->Parse();
		dataMatrix = readerMean->GetArray2DDataObject()->GetMatrix();

		// Check row and column header existence via heuristic: presence of all-NaNs

		for (int row = 0; row < dataMatrix.rows(); row++)
		{
			if (!std::isnan(dataMatrix[row][0]))
			{
				// Found a number in row headers, must assume this column meant to be numerical (i.e. no row headers).
				hasRowHeaders = false;
				break;
			}
			hasRowHeaders = true; // All entries in the first column are Not A Number
		}

		for (int col = 0; col < dataMatrix.cols(); col++)
		{
			if (!std::isnan(dataMatrix[0][col]))
			{
				// Found a number in column headers, must assume this row meant to be numerical (i.e. no column headers).
				hasColumnHeaders = false;
				break;
			}
			hasColumnHeaders = true;// All entries in the first row are Not A Number
		}

		if (hasRowHeaders)
		{
			readerMean->HasRowHeadersOn();
		}
		if (hasColumnHeaders)
		{
			readerMean->HasColumnHeadersOn();
		}
		// Re-parse with this new knowledge if headers are present at all
		if (hasRowHeaders || hasColumnHeaders)
		{
			readerMean->Parse();
			dataMatrix = readerMean->GetArray2DDataObject()->GetMatrix();
		}

		colHeaders = readerMean->GetArray2DDataObject()->GetColumnHeaders();
		rowHeaders = readerMean->GetArray2DDataObject()->GetRowHeaders();

		featuresMatrix.SetSize(dataMatrix.rows(), dataMatrix.columns());

		for (unsigned int i = 0; i < dataMatrix.rows(); i++)
			for (unsigned int j = 0; j < dataMatrix.cols(); j++)
				featuresMatrix(i, j) = dataMatrix(i, j);
	}
	catch (const std::exception& e1)
	{
		std::cerr << "Error reading the feature file in the input directory. Check permissions and existence. Error code : " + std::string(e1.what()) << std::endl;
		return std::make_tuple(false, false, false);
	}
	return std::make_tuple(true, hasRowHeaders, hasColumnHeaders);
}

std::tuple<bool, bool, bool> TrainingModule::GetLabelsFromFile(std::string labelsFilename, VariableLengthVectorType& labelsVector, std::vector<std::string>& rowHeaders, std::vector<std::string>& colHeaders)
{
	MatrixType dataMatrix;
	bool hasRowHeaders = false;
	bool hasColumnHeaders = false;
	try
	{
		CSVFileReaderType::Pointer readerMean = CSVFileReaderType::New();
		readerMean->SetFileName(labelsFilename);
		readerMean->SetFieldDelimiterCharacter(',');
		readerMean->HasColumnHeadersOff();
		readerMean->HasRowHeadersOff();
		readerMean->Parse();
		dataMatrix = readerMean->GetArray2DDataObject()->GetMatrix();
		
		// Check row and column header existence via heuristic: presence of all-NaNs

		for (int row = 0; row < dataMatrix.rows(); row++)
		{
			if (!std::isnan(dataMatrix[row][0]))
			{
				// Found a number in row headers, must assume this column meant to be numerical (i.e. no row headers).
				hasRowHeaders = false;
				break;
			}
			hasRowHeaders = true; // All entries in the first column are Not A Number
		}

		for (int col = 0; col < dataMatrix.cols(); col++)
		{
			if (!std::isnan(dataMatrix[0][col]))
			{
				// Found a number in column headers, must assume this row meant to be numerical (i.e. no column headers).
				hasColumnHeaders = false;
				break;
			}
			hasColumnHeaders = true;// All entries in the first row are Not A Number
		}

		if (hasRowHeaders)
		{
			readerMean->HasRowHeadersOn();
		}
		if (hasColumnHeaders)
		{
			readerMean->HasColumnHeadersOn();
		}
		// Re-parse with this new knowledge if headers are present at all
		if (hasRowHeaders || hasColumnHeaders)
		{
			readerMean->Parse();
			dataMatrix = readerMean->GetArray2DDataObject()->GetMatrix();
		}

		colHeaders = readerMean->GetArray2DDataObject()->GetColumnHeaders();
		rowHeaders = readerMean->GetArray2DDataObject()->GetRowHeaders();

		labelsVector.SetSize(dataMatrix.rows());

		for (unsigned int i = 0; i < dataMatrix.rows(); i++)
			labelsVector[i] = dataMatrix(i, 0);
	}
	catch (const std::exception& e1)
	{
		std::cerr << "Error reading the labels file in the input directory. Check permissions and existence. Error code : " + std::string(e1.what()) << std::endl;
		return std::make_tuple(false, false, false);
	}
	return std::make_tuple(true, hasRowHeaders, hasColumnHeaders);
}

double TrainingModule::GetBinaryClassificationBalancedAccuracy(VariableLengthVectorType& predictedLabels, VariableLengthVectorType& actualLabels)
{
	// calculate performance measures
	double TP = 0;
	double TN = 0;
	double FP = 0;
	double FN = 0;

	for (unsigned int index = 0; index < predictedLabels.GetSize(); index++)
	{
		if (predictedLabels[index] == 1 && actualLabels[index] == 1)
			TP++;
		else if (predictedLabels[index] == -1 && actualLabels[index] == -1)
			TN++;
		else if (predictedLabels[index] == 1 && actualLabels[index] == -1)
			FP++;
		else if (predictedLabels[index] == -1 && actualLabels[index] == 1)
			FN++;
		else
		{
			//If we get here we have some kind of bug in label loading/prediction
			//Or the initial labels were provided in the wrong form
			// TODO: Break out of this appropriately
			std::cerr << "Warning: labels do not match the binary classification problem. Please submit a bug report to CBICA Software." << std::endl;
		}
	}
	double overallAccuracy = (TP + TN) / (double)predictedLabels.GetSize();
	double sensitivity = TP / (TP + FN);
	double specificity = TN / (TN + FP);
	double balancedAccuracy = (sensitivity + specificity) / (double)2;

	return balancedAccuracy;
}

void TrainingModule::RemoveNans(VariableSizeMatrixType& mat)
{
	for (unsigned int index1 = 0; index1 < mat.Rows(); index1++)
	{
		for (unsigned int index2 = 0; index2 < mat.Cols(); index2++)
		{
			if (std::isnan(mat[index1][index2]))
				mat[index1][index2] = 0;
		}
	}
}
void TrainingModule::RemoveNans(VariableLengthVectorType& vec)
{
	for (unsigned int index1 = 0; index1 < vec.Size(); index1++)
	{
		if (std::isnan(vec[index1]))
			vec[index1] = 0;
	}
}

void TrainingModule::WriteNumericFilesWithHeaders(const std::string& filePath, const VariableSizeMatrixType& data, const std::vector<std::string>& rowHeaders, const std::vector<std::string>& colHeaders)
{
	std::ofstream myfile;
	myfile.open(filePath);

	if (rowHeaders.size() != data.Rows())
	{
		throw std::logic_error("Header/data row mismatch!");
	}
	if (colHeaders.size() != data.Cols())
	{
		throw std::logic_error("Header/data column mismatch!");
	}

	// First write out the whole set of column headers if needed
	if (!colHeaders.empty())
	{
		myfile << " ,"; // filler so we don't misalign
		for (unsigned int i = 0; i < colHeaders.size(); i++)
		{
			if (i == 0)
				myfile << colHeaders[i];
			else
				myfile << "," << colHeaders[i];
		}
		myfile << "\n";
	}

	bool useRowHeaders = true;
	if (rowHeaders.empty()) 
	{
		useRowHeaders = false;
	}

	// Now write data with row headers at start if needed
	for (unsigned int index1 = 0; index1 < data.Rows(); index1++)
	{
		if (useRowHeaders)
		{
			myfile << rowHeaders[index1] << ",";
		}

		for (unsigned int index2 = 0; index2 < data.Cols(); index2++)
		{
			if (index2 == 0)
				myfile << std::to_string(data[index1][index2]);
			else
				myfile << "," << std::to_string(data[index1][index2]);
		}
		myfile << "\n";
	}
	myfile.close();
}

void TrainingModule::WriteLabelFilesWithHeaders(const std::string& filePath, const VectorDouble& data, const std::vector<std::string>& rowHeaders, const std::vector<std::string>& colHeaders)
{
	std::ofstream myfile;
	myfile.open(filePath);

	if (!rowHeaders.empty() && (rowHeaders.size() != data.size()))
	{
		throw std::logic_error("Header/data row mismatch!");
	}

	if (colHeaders.size() > 1)
	{
		throw std::logic_error("Header/data column mismatch!");
	}

	// First write out the whole set of column headers if needed
	if (!colHeaders.empty())
	{
		myfile << " ,"; // filler so we don't misalign
		for (unsigned int i = 0; i < colHeaders.size(); i++)
		{
			if (i == 0)
				myfile << colHeaders[i];
			else
				myfile << "," << colHeaders[i];
		}
		myfile << "\n";
	}

	bool useRowHeaders = true;
	if (rowHeaders.empty())
	{
		useRowHeaders = false;
	}

	// Now write data with row headers at start if needed
	for (unsigned int index1 = 0; index1 < data.size(); index1++)
	{
		if (useRowHeaders)
		{
			myfile << rowHeaders[index1] << ",";
		}

		myfile << data[index1];
		myfile << "\n";
	}
	myfile.close();
}


/*void TrainingModule::RunThread(const TrainingModuleParameters& params)
{
	if (isCurrentlyRunning)
	{
		// The below is not great practice -- this doesn't allow the thread to clean up, etc.
		// However, QThread::quit() is not enough, nor is QThread::exit().
		// If the training module code were in a form where it could accept an abort signal and stop cleanly,
		// then we could safely use those.
		// But for now, the cost of waiting for the entire training module to finish when the user wants
		// to quit is too great.
		// TODO: Revisit this with QThread::currentThread->isInterruptionRequested in internal code
		// and in the other internal functions, break out if interruption is requested (esp in longer loops)

		// Need to stop any current processing on the worker thread before continuing.
		// The calling code should always check isCurrentlyRunning itself before getting here.
		m_workerThread->terminate();
	}
	isCurrentlyRunning = true;

	auto toBeExecuted = std::bind(&TrainingModule::Run, this, params);
	m_workerThread = QThread::create(toBeExecuted);

	connect(m_workerThread, &QThread::finished, m_workerThread, &QObject::deleteLater);
	connect(m_workerThread, SIGNAL(finished()), this, SLOT(onThreadFinished()));

	m_workerThread->start();
}*/

/*
void TrainingModule::onThreadFinished()
{
	isCurrentlyRunning = false;
	emit done(lastResult);
}*/
