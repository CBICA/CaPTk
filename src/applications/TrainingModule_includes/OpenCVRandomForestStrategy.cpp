/** 
\file OpenCVRandomForestStrategy.cpp

\brief The file containing the OpenCVRandomForestStrategy class. 
This class defines an implementation of the IClassifierStrategy interface to cover the linear SVM from OpenCV.



Author: Alexander Getka
https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2020 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html

*/
#pragma once


#include "OpenCVRandomForestStrategy.h"
#include "opencv2/core/core.hpp"
#include "opencv2/ml.hpp"
#include "itkCSVArray2DFileReader.h"
#include "CaPTkUtils.h"

typedef itk::CSVArray2DFileReader<double> ReaderType;
typedef itk::CSVArray2DFileReader<double> CSVFileReaderType;
typedef vnl_matrix<double> MatrixType;

VectorDouble OpenCVRandomForestStrategy::predict(const VariableSizeMatrixType& features, bool predictDistances, std::string modelDir)
{
	if (!isModelLoaded)
	{
		std::cerr << "Attempted to run prediction without a loaded model. Please submit a bug report to CBICA Software." << std::endl;
		throw std::logic_error("Attempted to run prediction without a model loaded.");
	}
	auto randomForest = _RF;

	// copy features into a cv mat for the random forest
	cv::Mat testingData = cv::Mat::zeros(features.Rows(), features.Cols(), CV_32F);
	for (int row_idx = 0; row_idx < testingData.rows; row_idx++)
	{
		for (int col_idx = 0; col_idx < testingData.cols; col_idx++)
		{
			float value = (float)features(row_idx, col_idx);
			testingData.at< float >(row_idx, col_idx) = value;
		}
	}
	
	VectorDouble predictionOutput;

	cv::Mat predicted(features.Rows(), 1, CV_32SC1);
	// Output distances (raw sums) if predictDistances is true, otherwise output class labels
	if (predictDistances)
	{
		predicted.convertTo(predicted, CV_32FC1);
		randomForest->predict(testingData, predicted, cv::ml::RTrees::RAW_OUTPUT);// prepare for float outputs from prediction
	}
	else
	{
		randomForest->predict(testingData, predicted);
		
		predicted.convertTo(predicted, CV_32FC1);
	}

	for (int i = 0; i < testingData.rows; i++)
	{
		float prediction = predicted.at<float>(i, 0);
		predictionOutput.push_back(prediction);
	}

	return predictionOutput;
}

bool OpenCVRandomForestStrategy::train(const VariableSizeMatrixType& features, const VectorDouble& labels, std::string outDir)
{

	// Copy features/labels into a cv::mat
	cv::Mat trainingFeatures = cv::Mat::zeros(features.Rows(), features.Cols(), CV_32FC1);
	cv::Mat trainingLabels = cv::Mat::zeros(labels.size(), 1, CV_32F);

	for (int copyDataCounter = 0; copyDataCounter < trainingFeatures.rows; copyDataCounter++)
	{
		trainingLabels.ptr< float >(copyDataCounter)[0] = labels[copyDataCounter];
		for (int copyDataCounter2 = 0; copyDataCounter2 < trainingFeatures.cols; copyDataCounter2++)
			trainingFeatures.ptr< float >(copyDataCounter)[copyDataCounter2] = features[copyDataCounter][copyDataCounter2];
	}
	trainingLabels.convertTo(trainingLabels, CV_32SC1);

	auto trainingData = cv::ml::TrainData::create(trainingFeatures, cv::ml::ROW_SAMPLE, trainingLabels);

	auto randomForest = cv::ml::RTrees::create();

	// Since FS changes variable count throughout the process, probably best to leave this alone, otherwise it'll scale really weirdly.
	// randomForest->setActiveVarCount(params.randomForestActiveVarCount);

	// // Construct termination criteria from params:
	// float/double params.randomForestEpsilon OOB error has to be lower than this to pass
	// int params.randomForestMaxIterations maximum amount of trees that can be made
	auto terminator = cv::TermCriteria(cv::TermCriteria::COUNT + cv::TermCriteria::EPS, params.randomForestMaxIterations, params.randomForestEpsilon);
	randomForest->setTermCriteria(terminator);
	randomForest->train(trainingData);
	

	_RF = randomForest; // cv::Ptr takes care of managing this 
	this->isModelLoaded = true;
	return true;
}


void OpenCVRandomForestStrategy::loadModel(std::string modelDir)
{
	auto modelDirToLoad = modelDir;
	if (modelDir.empty())
	{
		modelDirToLoad = params.modelDirectory;
	}

	_RF = cv::Algorithm::load<cv::ml::RTrees>(modelDirToLoad + "/OpenCV_RF_Model.xml");
	this->isModelLoaded = true;
}
void OpenCVRandomForestStrategy::saveModel(std::string outputDir)
{
	if (!isModelLoaded)
	{
		std::cerr << "Attempted to save an undefined model. Please submit a bug report to CBICA Software." << std::endl;
		throw std::logic_error("Attempted to save an undefined model.");
	}

	auto outputDirectory = outputDir;
	if (outputDir.empty())
	{
		outputDirectory = params.outputDirectory;
	}

	_RF->save(outputDirectory + "/OpenCV_RF_Model.xml");
}

std::string OpenCVRandomForestStrategy::getModelInfoAsString()
{

	// This should generally be called after some kind of optimized training.
	// So after running train(), this function will get info from the best model trained there.
	
	// We want to pull info from the latest model -- if unloaded we need to fail.
	if (!isModelLoaded)
	{
		std::cerr << "Attempted to read from an undefined model. Please submit a bug report to CBICA Software. " << std::endl;
		throw std::logic_error("Attempted to read from an undefined model.");
	}

	std::string result = "";

	// TODO: fix for RANDOM FOREST stuff
	// Doing so here allows us to query this from top-level TrainingModule,
	// instead of printing this info on every individual training call
	result += "No additional information about this model is available. \n";
	// Maybe print epsilon and max_iterations termination criteria here? OOB error? Var importance at final step?

	

	return result;
}