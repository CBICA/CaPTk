/** 
\file EffectSizeFeatureSelectionStrategy.h

\brief The header file containing the EffectSizeFeatureSelectionStrategy class. 
This defines an extension of the IFeatureSelectionStrategy interface to cover Effect-size based feature selection.
This is based on a refactoring of CaPTk TrainingModule code originally written by Saima Rathore.

Author: Alexander Getka
https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2020 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html

*/
#pragma once

#include "IFeatureSelectionStrategy.h"


class EffectSizeFeatureSelectionStrategy : public IFeatureSelectionStrategy
{

// Perform feature selection on the given data using the given classification method/strategy.
// This will pull updated parameters from params on each run. It should ONLY use params specifically related to feature selection!
// returns std::vector<int> Vector of final selected feature indices
std::vector<int> PerformFeatureSelectionBasedTraining(
	const VariableSizeMatrixType& features, const VectorDouble& labels, IClassifierStrategy* classifier);

// Some feature selection approaches integrate cross-validation

private:

VectorDouble CalculateEffectSizes(const VariableSizeMatrixType training_features, VectorDouble targets);


// template function to sort 
template <typename T>
std::vector<size_t> sortIndexesByValue(const std::vector<T>& v)
{
    // initialize original index locations
    std::vector<size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] > v[i2]; });

    return idx;
}

};