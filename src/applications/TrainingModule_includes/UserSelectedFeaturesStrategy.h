/** 
\file UserSelectedFeaturesStrategy.h

\brief The header file containing the declaration of the UserSelectedFeaturesStrategy class.

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

#include "IFeatureSelectionStrategy.h"


class UserSelectedFeaturesStrategy : public IFeatureSelectionStrategy
{

// Perform feature selection on the given data using the given classification method/strategy.
// This will pull updated parameters from params on each run. It should ONLY use params specifically related to feature selection!
// returns std::vector<int> Vector of final selected feature indices
std::vector<int> PerformFeatureSelectionBasedTraining(
	const VariableSizeMatrixType& features, const VectorDouble& labels, IClassifierStrategy* classifier);

};