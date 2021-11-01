/** 
\file SelectOptimalFeatureSelectionStrategy.cpp

\brief The file containing the SelectOptimalFeatureSelectionStrategy class. 
This defines an implementation of the IFeatureSelectionStrategy interface that iterates over all other FS methods and chooses that which performs best.
This should be recommended for users who don't have an idea of the model they want to train, 
or who want to see the best possible performance that can be achieved with the tool.
It will take a multiplicatively long time to run, especially when also optimizing the ML strategy and/or other parameters.


Author: Alexander Getka
https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2020 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html

*/
#pragma once

#include "IFeatureSelectionStrategy.h"
#include "SelectOptimalFeatureSelectionStrategy.h"
#include "opencv2/core/core.hpp"
#include "opencv2/ml.hpp"



// Perform feature selection on the given data using the given classification method/strategy.
// This will pull updated parameters from params on each run. It should ONLY use params specifically related to feature selection!
std::vector<int> SelectOptimalFeatureSelectionStrategy::PerformFeatureSelectionBasedTraining(
	const VariableSizeMatrixType& features, const VectorDouble& labels, IClassifierStrategy* classifier)
{
	std::cout << "Entering feature selection optimization..." << std::endl;
	int numberOfSelectedFeatures = 0;
	std::vector<int> FinalSelectedFeatures;

	// TODO: Implement loop over other FS methods here, select best to populate


    std::cout << "No. of selected features: " << FinalSelectedFeatures.size() << std::endl;

	std::cout << "Completed feature selection optimization." << std::endl;
	return FinalSelectedFeatures;
}



