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

	// TBD: Consider changing to a dictionary/map-type data structure to eliminate header changes
	// Could use strings (or type, value string tuples) as an intermediate for conversion

public:
	int classifierType;
	int optimizationType;
	std::string inputFeaturesFile;
	std::string inputLabelsFile;
	std::string outputDirectory;
	std::string modelDirectory;
	int configurationType;
	int featureSelectionType;
	int crossValidationType;
	int folds;
	int cMin, cMax;
	int gMin, gMax;

};