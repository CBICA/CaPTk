/** 
\file OpenCVSGDSVMStrategy.h

\brief The header file containing the OpenCVSGDSVMStrategy class. 
This class defines an implementation of the IClassifierStrategy interface to cover the SGD-SVM (stochastic gradient descent-based SVM) from OpenCV.
This is specifically for binary classification.


Author: Alexander Getka
https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2021 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html

*/
#pragma once


#include "IClassifierStrategy.h"
#include "opencv2/ml.hpp"

class OpenCVSGDSVMStrategy : public IClassifierStrategy
{

	// 

public:

	// Run prediction using the current context.
	// Returns vector of prediction as type double. Caller is responsible for using this or writing to disk.
	// param modelFile Optionally specify a model filename, otherwise "SVM_Model.xml" will be used.
	// param predictDistances bool : true to output distances (raw scores), default (false) to just output classification.
	VectorDouble predict(const VariableSizeMatrixType& features, bool predictDistances = false, std::string modelFile = "");

	// Run training using the current context. 
	// Returns true if successful, false otherwise.
	// This method should pull up-to-date parameters pulled from the "params" variable on each run.
	// param outDir Optionally specify an output filename to use (e.g. for OpenCV XML model output).
	// All output files will be placed under the output directory specified in params.
	bool train(const VariableSizeMatrixType& features, const VectorDouble& labels, std::string outDir = "");

	void loadModel(std::string modelDir);
	void saveModel(std::string outputDir);
	std::string getModelInfoAsString();

protected:
	// This is the model saved internally.
	// predict() and getHyperparametersAsString() will refer to this.
	cv::Ptr<cv::ml::SVMSGD> _SGDSVM;

};