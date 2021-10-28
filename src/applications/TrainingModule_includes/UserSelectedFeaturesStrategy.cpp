/** 
\file SelectOptimalFeatureSelectionStrategy.cpp

\brief The  file containing the UserSelectedFeaturesStrategy class implementation.
.
This class defines a "feature selection" method that actually just directly passes user-specified features.
Users should be able to specify either feature names (if feature headers were provided)
or indices (usable in any case, although needs handling if we remove, say, subject IDs from the leftmost column).

Author: Alexander Getka
https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2020 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html

*/
#pragma once

#include "UserSelectedFeaturesStrategy.h"

std::vector<int> UserSelectedFeaturesStrategy::PerformFeatureSelectionBasedTraining(
	const VariableSizeMatrixType& features, const VectorDouble& labels, IClassifierStrategy* classifier)
{
	std::cout << "Entering feature selection optimization..." << std::endl;
	int numberOfSelectedFeatures = 0;

	// TODO: Implement letting users specify features via their name (e.g. T1_NCR_Bin5 or Age), not just index

    std::vector<int> FinalSelectedFeatures;

	// TODO: Implement selection with indices read from param
    // push back the feature set that produced the best performance
    //for (int index = 0; index < params.maxNumberOfFeatures; index++)
    //    FinalSelectedFeatures.push_back(indices[index]); 
    std::cout << "No. of selected features: " << FinalSelectedFeatures.size() << std::endl;

	std::cout << "Completed feature selection optimization." << std::endl;
	return FinalSelectedFeatures;
}



