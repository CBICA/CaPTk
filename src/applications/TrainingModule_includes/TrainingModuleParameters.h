/** 
\file TrainingModuleParameters.h

\brief The header file containing the TrainingModuleParameters class. 
This minimal class describes a container of parameters/settings to be passed to the TrainingModule class. 
An internet search for the "Parameter Object" design pattern should explain adequately.
This pattern prevents function signatures from blowing up across the TrainingModule and keeps them local to here.


Author: Alexander Getka
https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2020 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html

*/
#pragma once

#include <string>

class TrainingModuleParameters
{
	// Addition of a new parameter to the training module needs a change to this class.
	// The GUI constructs an instance of this class from collected options. 
	// This instance is passed to the TrainingModule backend upon user confirmation.
	// The instance can be plumbed through the call stack as needed, acting as a "context". 
	// Then any code that needs an input parameter specifically gets it from the instance of this class.
	// This is, at most, 1 extra param to TrainingModule functions, and can sometimes eliminate redundant params.
	// As a rule, any calculated or loaded values (such as feature vectors, etc) should be passed separately.

	// Pro: Keeps function signatures small & (more) constant
	// Pro: Localized here rather than spread across Training-related and/or GUI code
	// Con: This is still a header change for each new parameter (recompiles).  

	// TBD: Consider changing to a dictionary/map-type data structure.
	// Could use strings (or [type, default_value, shortname, description] string tuples?) as an intermediate for conversion
	// Even just including default-getters might be a good idea for GUI 
	// Also: Additional checks for validity? Something like bool [var]_hasBeenSet
public:
	// These parameters must be set by the user interface generating this object.
	// TrainingModule code should only use these if they are valid for the given input.
	// (ex: don't use gamma [gmin/gmax] for any code that presupposes linear kernel SVM

	// Required
	CAPTK::ClassifierType classifierType;
	CAPTK::OptimizationType optimizationType;
	std::string inputFeaturesFile;
	std::string inputLabelsFile;
	std::string outputDirectory;
	std::string modelDirectory;
	CAPTK::ClassificationConfigurationType configurationType;
	CAPTK::FeatureSelectionType featureSelectionType;
	CAPTK::CrossValidationType crossValidationType;
	int folds;

	// if zero, select features until performance peaks. Else, only select up to this many.
	int maxNumberOfFeatures = 0; 

	// Options with valid defaults
	bool testPredictionsAgainstProvidedLabels = false;
	
	// SVM parameters for C_SVC with base/max/min for grid search
	double cMin = -5; double cLogBase = 2; double cMax = 5;
	double gMin = -5; double gLogBase = 2; double gMax = 5;
	double coefMin = 0.1; double coefLogBase = 14; double coefMax = 300;
	double degreeMin = 0.01; double degreeLogBase = 7; double degreeMax = 4;

	// Termination criteria for random forest
	double randomForestEpsilon = 0.1;
	int randomForestMaxIterations = 50;

};