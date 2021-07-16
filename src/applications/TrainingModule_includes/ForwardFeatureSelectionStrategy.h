/** 
\file ForwardFeatureSelectionStrategy.h

\brief The header file containing the ForwardFeatureSelectionStrategy class. 
This defines an extension of the IFeatureSelectionStrategy interface to cover forward feature selection (FFS).
This is based on a refactoring of CaPTk TrainingModule code originally written by Saima Rathore.

Author: Alexander Getka
https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2020 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html

*/
#pragma once

#include "IFeatureSelectionStrategy.h"


class ForwardFeatureSelectionStrategy : public IFeatureSelectionStrategy
{

// Perform feature selection on the given data using the given classification method/strategy.
// This will pull updated parameters from params on each run. It should ONLY use params specifically related to feature selection!
// returns std::vector<int> Vector of final selected features indices
std::vector<int> PerformFeatureSelectionBasedTraining(
	const VariableSizeMatrixType& features, const VectorDouble& labels, IClassifierStrategy* classifier);

private:

std::vector<int> GetUnselectedFeatures(const std::vector<int>& selectedFeatures, int featureSize);




};