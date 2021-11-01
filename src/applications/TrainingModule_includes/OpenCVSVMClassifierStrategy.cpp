/** 
\file OpenCVSVMClassifierStrategy.cpp

\brief The file containing the OpenCVSVMClassifierStrategy class. 
This class defines an implementation of the IClassifierStrategy interface to cover the linear SVM from OpenCV.



Author: Alexander Getka
https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2020 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html

*/
#pragma once


#include "OpenCVSVMClassifierStrategy.h"
#include "opencv2/core/core.hpp"
#include "opencv2/ml.hpp"
#include "itkCSVArray2DFileReader.h"
#include "CaPTkUtils.h"

typedef itk::CSVArray2DFileReader<double> ReaderType;
typedef itk::CSVArray2DFileReader<double> CSVFileReaderType;
typedef vnl_matrix<double> MatrixType;

VectorDouble OpenCVSVMClassifierStrategy::predict(const VariableSizeMatrixType& features, bool predictDistances, std::string modelDir)
{
	if (!isModelLoaded)
	{
		std::cerr << "Attempted to run prediction without a loaded model. Please submit a bug report to CBICA Software." << std::endl;
		throw std::logic_error("Attempted to run prediction without a model loaded.");
	}
	auto svm = _SVM;

	// copy features into a cv mat for the SVM
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
		svm->predict(testingData, predicted, cv::ml::SVM::RAW_OUTPUT);// prepare for float outputs from prediction
	}
	else
	{
		svm->predict(testingData, predicted);
		
		predicted.convertTo(predicted, CV_32FC1);
	}

	for (int i = 0; i < testingData.rows; i++)
	{
		float prediction = predicted.at<float>(i, 0);
		predictionOutput.push_back(prediction);
	}

	return predictionOutput;
}

bool OpenCVSVMClassifierStrategy::train(const VariableSizeMatrixType& features, const VectorDouble& labels, std::string outDir)
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

	// convert to TrainData form for trainAuto compatibility
	auto trainingData = cv::ml::TrainData::create(trainingFeatures, cv::ml::ROW_SAMPLE, trainingLabels);

	auto svm = cv::ml::SVM::create();
	svm->setType(cv::ml::SVM::C_SVC);

	// Do training across the parameter grid

	// Check relevant parameters
	bool cSearchRelevant = true; // Used whenever C_SVC is
	// Now to the case-by-case ones
	bool gammaSearchRelevant = false, coefSearchRelevant = false, degreeSearchRelevant = false;
	bool nuSearchRelevant = false, pSearchRelevant = false; // Not used for C_SVC
	if (params.classifierType == CAPTK::ClassifierType::CLASS_TYPE_SVM_LINEAR)
	{
		svm->setKernel(cv::ml::SVM::LINEAR);
	}
	else if (params.classifierType == CAPTK::ClassifierType::CLASS_TYPE_SVM_RBF)
	{
		svm->setKernel(cv::ml::SVM::RBF);
		gammaSearchRelevant = true;
	}
	else if (params.classifierType == CAPTK::ClassifierType::CLASS_TYPE_SVM_POLYNOMIAL)
	{
		svm->setKernel(cv::ml::SVM::POLY);
		gammaSearchRelevant = true;
		coefSearchRelevant = true;
		degreeSearchRelevant = true;
		svm->setDegree(2); // Needs to be positive, OpenCV seems to complain easily
	}
	else if (params.classifierType == CAPTK::ClassifierType::CLASS_TYPE_SVM_SIGMOID)
	{
		svm->setKernel(cv::ml::SVM::SIGMOID);
		gammaSearchRelevant = true;
		coefSearchRelevant = true;
	}
	else if (params.classifierType == CAPTK::ClassifierType::CLASS_TYPE_SVM_CHISQUARED)
	{
		svm->setKernel(cv::ml::SVM::CHI2);
		gammaSearchRelevant = true;
	}
	else if (params.classifierType == CAPTK::ClassifierType::CLASS_TYPE_SVM_INTERSECTION)
	{
		svm->setKernel(cv::ml::SVM::INTER);
	}
	
	if (params.optimizationType == CAPTK::OptimizationType::OPT_TYPE_OFF)
	{
		// perform plain training with default parameters
		svm->train(trainingData);
	}
	else // Optimization ON
	{
		// TODO: replace TrainAuto with manual parameter-looping -- this needs to be controllable from params.
		// That way, FS method can determine when parameter optimization is done.

		// Construct grids from parameters
		// TODO: get LogStep for each from params/users
		auto cGrid = cv::ml::SVM::getDefaultGrid(cv::ml::SVM::C);
		cGrid.maxVal = std::pow(params.cLogBase, params.cMax);
		cGrid.logStep = params.cLogBase;
		cGrid.minVal = std::pow(params.cLogBase, params.cMin);
		auto gGrid = cv::ml::SVM::getDefaultGrid(cv::ml::SVM::GAMMA);
		gGrid.maxVal = std::pow(params.gLogBase, params.gMax);
		gGrid.logStep = params.gLogBase;
		gGrid.minVal = std::pow(params.gLogBase, params.gMin);

		// TODO: check the maximums here for log values as for c and g if trainAuto doesn't work out
		auto coefGrid = cv::ml::SVM::getDefaultGrid(cv::ml::SVM::COEF);
		coefGrid.maxVal = params.coefMax;
		coefGrid.logStep = params.coefLogBase;
		coefGrid.minVal = params.coefMin;
		auto degreeGrid = cv::ml::SVM::getDefaultGrid(cv::ml::SVM::DEGREE);

		//degreeGrid.maxVal = params.degreeMax;
		//degreeGrid.logStep = params.degreeLogBase;
		//degreeGrid.minVal = params.degreeMin;

		// TODO Replace?
		//svm->train(trainingData);
		svm->trainAuto(trainingData, 5, cGrid, gGrid,
			cv::ml::SVM::getDefaultGrid(cv::ml::SVM::P),
			cv::ml::SVM::getDefaultGrid(cv::ml::SVM::NU),
			coefGrid, degreeGrid, true);
	}

	/*
	// Handle cases where the user doesn't want to do a grid search on one or more params
	if (params.cMin == params.cMax)
	{
		svm->setC(std::pow(2, params.cMin));
		cSearchRelevant = false;
	}
	if (params.gMin == params.gMax)
	{
		svm->setGamma(std::pow(2, params.gMin));
		gammaSearchRelevant = false;
	}
	if (params.coefMin == params.coefMax)
	{
		svm->setCoef0(params.coefMin);
		coefSearchRelevant = false;
	}
	if (params.degreeMin == params.degreeMax)
	{
		svm->setDegree(params.degreeMin);
		degreeSearchRelevant = false;
	}
	*/

	/*
	// Handle cases where the user wants to grid-search
	// This is messy -- TODO: get a solution for this that is cleaner
	for (double c = params.cMin; c <= params.cMax; c++) // ;)
	{
		for (double g = params.gMin; g <= params.gMax; g++)
		{
			for (double coef = params.coefMin; coef <= params.coefMax; coef++)
			{
				for (double degree = params.degreeMin; degree <= params.degreeMax; degree++)
				{
					svm->setC(c);
					svm->setGamma(g);
					svm->setCoef0(coef);
					svm->setDegree(degree);
				}
			}
		}
	}*/

	
													 
	//svm->train(trainingData);

	

	_SVM = svm; // cv::Ptr takes care of managing this 
	this->isModelLoaded = true;
	return true;
}


void OpenCVSVMClassifierStrategy::loadModel(std::string modelDir)
{
	auto modelDirToLoad = modelDir;
	if (modelDir.empty())
	{
		modelDirToLoad = params.modelDirectory;
	}

	_SVM = cv::Algorithm::load<cv::ml::SVM>(modelDirToLoad + "/OpenCV_SVM_Model.xml");
	this->isModelLoaded = true;
}
void OpenCVSVMClassifierStrategy::saveModel(std::string outputDir)
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

	_SVM->save(outputDirectory + "/OpenCV_SVM_Model.xml");
}

std::string OpenCVSVMClassifierStrategy::getModelInfoAsString()
{

	// This should generally be called after some kind of optimized training.
	// for OpenCV SVMs, cv::ml::svm::trainAuto takes care of finding the optimal parameters (via cross-validation).
	// So after running train(), this function will get info from the best model trained there.
	
	// We want to pull info from the latest model -- if unloaded we need to fail.
	if (!isModelLoaded)
	{
		std::cerr << "Attempted to read from an undefined model. Please submit a bug report to CBICA Software. " << std::endl;
		throw std::logic_error("Attempted to read from an undefined model.");
	}

	std::string result = "";
	// Add some information on final-chosen hyperparameters here.
	// Doing so here allows us to query this from top-level TrainingModule,
	// instead of printing this info on every individual training call
	result += "Final optimal SVM hyperparameters found during training: \n";

	// Need to handle variations on OpenCV SVM -- "P" isn't relevant for Linear, etc.
	if (params.classifierType == CAPTK::ClassifierType::CLASS_TYPE_SVM_LINEAR)
	{
		result += "C = " + std::to_string(_SVM->getC()) + "\n";
	}
	else if (params.classifierType == CAPTK::ClassifierType::CLASS_TYPE_SVM_RBF)
	{
		result += "C = " + std::to_string(_SVM->getC()) + "\n";
		result += "Gamma = " + std::to_string(_SVM->getGamma()) + "\n";
	}
	

	return result;
}