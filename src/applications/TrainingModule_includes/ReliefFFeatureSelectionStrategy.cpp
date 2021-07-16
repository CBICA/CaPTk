/** 
\file ReliefFFeatureSelectionStrategy.cpp

\brief The file containing the ReliefFSelectionStrategy class. 
This defines an implementation of the IFeatureSelectionStrategy interface to cover  feature selection with the RELIEF-F algorithm.
An informal description from wikipedia: https://en.wikipedia.org/wiki/Relief_(feature_selection)#ReliefF_Algorithm
Or the paper for RELIEF-F at https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.56.4740&rep=rep1&type=pdf

Author: Alexander Getka
https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2020 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html

*/
#pragma once

#include "IFeatureSelectionStrategy.h"
#include "ReliefFFeatureSelectionStrategy.h"
#include "opencv2/core/core.hpp"
#include "opencv2/ml.hpp"



// Perform feature selection on the given data using the given classification method/strategy.
// This will pull updated parameters from params on each run. It should ONLY use params specifically related to feature selection!
std::vector<int> ReliefFFeatureSelectionStrategy::PerformFeatureSelectionBasedTraining(
	const VariableSizeMatrixType& features, const VectorDouble& labels, IClassifierStrategy* classifier)
{
	std::cout << "Entering RELIEF-F feature selection." << std::endl;
	int numberOfSelectedFeatures = 0;

    VariableSizeMatrixType scaledFeatures;
    scaledFeatures.SetSize(features.Rows(), features.Cols());
    
    // Algorithm notes that you should scale features to 0-1 first
    // TODO: DO THIS
    for (unsigned int i = 0; i < features.Cols(); i++)
    {
        double thisFeatureMaximum = 0;
        double thisFeatureMinimum = 0;
        for (unsigned int j = 0; j < features.Rows(); j++)
        {
            if (features(j, i) > thisFeatureMaximum)
            {
                thisFeatureMaximum = features(j, i);
            }
            else if (features(j, i) < thisFeatureMinimum)
            {
                thisFeatureMinimum = features(j, i);
            }
        }

        for (unsigned int j = 0; j < features.Rows(); j++)
        {
            scaledFeatures(j, i) = ( (features(j, i) - thisFeatureMinimum) / (thisFeatureMaximum - thisFeatureMinimum) );
        }



    }

    int numberOfIterations = labels.size(); // This will use all labels currently.
    // TODO: Allow users to potentially set a number of iterations (exhaustive or a random subset of selectable size)


    // create initial weight vector properly
    auto weights = std::vector<double>();
    weights.resize(scaledFeatures.Cols());
    std::fill(weights.begin(), weights.end(), 0.0);

    for (unsigned int i = 0; i < labels.size(); i++) // again, plain loop -- see note above for random subset
    {
        // Find the nearest neighbors via Manhattan (L1) norm aka "taxicab" norm (dX + dY)
        // NOT the Euclidean distance ( sqrt(dX^2 + dY^2) ) used by regular RELIEF!

        unsigned int nearestHitIndex;
        unsigned int nearestMissIndex;

        double currentLeastDistanceToSameClass = std::numeric_limits<double>::max(); // start high, shrink as we go
        double currentLeastDistanceToDifferentClass = std::numeric_limits<double>::max();
        for (unsigned int j = 0; j < labels.size(); j++) // Compare this case to each other case for calculating nearest
        {
            if (j == i) // Can't compare the case to itself
            {
                continue; 
            }

            double currentDistance = 0; // Loop over features to get distance
            for (unsigned int featureIndex = 0; featureIndex < scaledFeatures.Cols(); featureIndex++)
            {
                currentDistance += std::abs(scaledFeatures(j, featureIndex) - scaledFeatures(i, featureIndex));
            }

            if (labels[j] == labels[i]) // Same class -- "near hit"
            {
                if (currentDistance < currentLeastDistanceToSameClass)
                {
                    currentLeastDistanceToSameClass = currentDistance;
                    nearestHitIndex = j;
                }
            }
            else if (labels[j] != labels[i]) // different class -- "near miss"
            {
               if (currentDistance < currentLeastDistanceToDifferentClass)
               {
                   currentLeastDistanceToDifferentClass = currentDistance;
                   nearestMissIndex = j;
               }
            }
        }

        //  Update weight vector by iterating over each feature
        // nearHit = nearest same-class entry's Xi 
        // nearMiss = nearest different-class entry's Xi
        // Wi  = Wi  - abs(Xi - nearHit) + abs(Xi - nearMiss)

        for (unsigned int featureIndex = 0; featureIndex < features.Cols(); featureIndex++)
        {
            weights[featureIndex] = weights[featureIndex] 
                - std::abs(scaledFeatures(i, featureIndex) - scaledFeatures(nearestHitIndex, featureIndex))
                + std::abs(scaledFeatures(i, featureIndex) - scaledFeatures(nearestMissIndex, featureIndex));
        }

    }

    std::vector<double> relevances; // Weights divided by iteration count
    for (unsigned int i = 0; i < weights.size(); i++)
    {
        relevances.push_back(weights[i] / (double)numberOfIterations);
    }

    // According to the strict RELIEF-F algorithm, we should set a threshold, tau
    // But this is parametrization. For right now let's just take the top N.
    // TODO: expose this as an option/fix this
    std::vector<size_t> indices = sortIndexesByValue(relevances);

    std::vector<int> FinalSelectedFeatures;

    // push back the feature set that produced the best performance
    for (int index = 0; index < params.maxNumberOfFeatures; index++)
        FinalSelectedFeatures.push_back(indices[index]); 
    std::cout << "No. of selected features: " << FinalSelectedFeatures.size() << std::endl;

	std::cout << "Completed RELIEF-F feature selection." << std::endl;
	return FinalSelectedFeatures;
}



