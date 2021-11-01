/** 
\file SelectOptimalFeatureSelectionStrategy.h

\brief The header file containing the SelectOptimalFeatureSelectionStrategy class. 
This defines an implementation of the IFeatureSelectionStrategy interface that iterates over all other FS methods and chooses that which performs best.
This should be recommended for users who don't have an idea of the model they want to train, 
or who want to see the best possible performance that can be achieved with the tool.
It will take a multiplicatively long time to run, especially when also optimizing the ML strategy and/or other parameters..

Author: Alexander Getka
https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2020 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html

*/
#pragma once

#include "IFeatureSelectionStrategy.h"


class SelectOptimalFeatureSelectionStrategy : public IFeatureSelectionStrategy
{

// Perform feature selection on the given data using the given classification method/strategy.
// This will pull updated parameters from params on each run. It should ONLY use params specifically related to feature selection!
// returns std::vector<int> Vector of final selected feature indices
std::vector<int> PerformFeatureSelectionBasedTraining(
	const VariableSizeMatrixType& features, const VectorDouble& labels, IClassifierStrategy* classifier);

// Some feature selection approaches integrate cross-validation

private:


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