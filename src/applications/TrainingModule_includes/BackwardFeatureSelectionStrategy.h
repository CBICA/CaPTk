/** 
\file BackwardFeatureSelectionStrategy.h

\brief The header file containing the BackwardFeatureSelectionStrategy class. 
This defines an implementation of the IFeatureSelectionStrategy interface to cover sequential backwards based feature selection.
This is also referred to as "backward feature elimination".

Author: Alexander Getka
https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2020 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html

*/
#pragma once

#include "IFeatureSelectionStrategy.h"


class BackwardFeatureSelectionStrategy : public IFeatureSelectionStrategy
{

// Perform feature selection on the given data using the given classification method/strategy.
// This will pull updated parameters from params on each run. It should ONLY use params specifically related to feature selection!
// returns std::vector<int> Vector of final selected features indices
std::vector<int> PerformFeatureSelectionBasedTraining(
	const VariableSizeMatrixType& features, const VectorDouble& labels, IClassifierStrategy* classifier);

};