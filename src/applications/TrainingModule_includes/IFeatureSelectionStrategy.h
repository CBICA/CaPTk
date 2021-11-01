/** 
\file IFeatureSelectionStrategy.h

\brief This abstract class defines an interface that provides access to a cross-validation technique.
Classes that implement this represent an individual strategy.
This can be consumed by a FeatureSelectionStrategy (for FS techniques that perform cross-validation iteratively, for example).
This can also be used as a wrapper around the training process in general.


Author: Alexander Getka
https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2020 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html

*/
#pragma once
#include "CaPTkEnums.h" 
#include "CaPTkDefines.h"
#include <string>
#include <vector>
#include "TrainingModuleParameters.h"
#include "IClassifierStrategy.h"

// This class is responsible for performing feature selection, given a classifier strategy.
class IFeatureSelectionStrategy
{

// 

public:

TrainingModuleParameters params; // This can be accessed externally to change parameters
// e.g. EffectSizeFeatureSelectionStrategy a; a.SetParameters(userSpecifiedParams); a.params.someFSParameter = someValue;


// Copies parameters from a TrainingModuleParameters reference. 
// This sets the general context for execution of training and prediction. 
// If you need to over-ride some parameters in code (such as changing output for cross-validation) there are 2 options:
// 1. Change the params values here after-the-fact (a.params.outputDir = "/some/path/on/disk";)
// 2. Intercept the TrainingModuleParameters when you receive it, copy it, and modify it before passing.
void SetParameters(const TrainingModuleParameters& params) // Implementing here to keep this header-only
{
	this->params = params;
}
;

// Perform feature selection on the given data using the given classification method/strategy.
// This will pull updated parameters from params on each run. It should ONLY use params specifically related to feature selection!
// returns std::vector<int> of selected feature indices
virtual std::vector<int> PerformFeatureSelectionBasedTraining(
	const VariableSizeMatrixType& features, const VectorDouble& labels, IClassifierStrategy* classifier) = 0;

// virtual destructor ensures any derived destructors are called (avoids surprises down the line).
virtual ~IFeatureSelectionStrategy() {};

protected:


// This space is reserved for implementation-specific functions and variables.
// This won't be used here in the interface, but it can be in a class that implements this.

};