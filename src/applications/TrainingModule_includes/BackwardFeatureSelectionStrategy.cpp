/** 
\file BackwardFeatureSelectionStrategy.h

\brief The implementation file containing the BackwardFeatureSelectionStrategy class. 
This defines an implementation of the IFeatureSelectionStrategy interface to cover sequential backwards feature selection.
This is also referred to as "recursive feature elimination".

Author: Alexander Getka
https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2020 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html

*/
#pragma once

#include "BackwardFeatureSelectionStrategy.h"
//TBD: Merge backward and forward FS into one SequentialFS strategy that accepts an evaluation function/strategy?
// e.g. pass in an F-test strategy to run univariate F-test as the elimination criterion

// Summary of this approach:
// User wants to keep a number of features (M)
// start with all features selected (N total features)
// For all N until N = M
// -> For each remaining feature:
// -> -> Train a model without that feature 
// -> -> Get Cross-validated Balanced Accuracy of that model
// -> Find the feature-elimination that produced the best Cross-validated Balanced Accuracy
// -> Eliminate that feature


// Perform feature selection on the given data using the given classification method/strategy.
// This will pull updated parameters from params on each run. It should ONLY use params specifically related to feature selection!
std::vector<int> BackwardFeatureSelectionStrategy::PerformFeatureSelectionBasedTraining(
	const VariableSizeMatrixType& features, const VectorDouble& labels, IClassifierStrategy* classifier)
{
    std::cout << "Entering recursive feature elimination." << std::endl;
    auto optimizationType = params.optimizationType;
    auto cvtype = params.crossValidationType;
    auto classifiertype = params.classifierType;

    std::vector<int> selectedFeatureIndices; // init to all features
    selectedFeatureIndices.resize(features.Cols());
    for (int i = 0; i < features.Cols(); i++)
    {
        selectedFeatureIndices[i] = i;
    }
    std::vector<double> CrossValidatedBalancedAccuraciesFinal;
    std::vector<int> eliminatedFeatureIndices;

    //feature selection mechanism
    //---------------------------
    // The idea is that we gradually select fewer features here, evaluating performance at each step
    // later, we will use this information to decide on a final number of features.
    while (selectedFeatureIndices.size() > 1)
    {
        std::vector<double> CrossValidatedBalancedAccuracies;
        // iterate over the remaining features to find the best performing elimination
        for (unsigned int featureNo = 0; featureNo < selectedFeatureIndices.size(); featureNo++)
        {
            VariableSizeMatrixType reducedFeatureSet;
            reducedFeatureSet.SetSize(labels.size(), selectedFeatureIndices.size());

            // copy all selected features except the current one under observation to the current feature set
            // basically we are testing elimination of the current feature
            for (unsigned int j = 0; j < reducedFeatureSet.Rows(); j++)
            {
                int copiedFeatureCount = 0;
                for (unsigned int k = 0; k < selectedFeatureIndices.size(); k++ )
                {
                    if (k != featureNo) // Non-eliminated features
                    {
                        double value = (int)(features(j, selectedFeatureIndices[k]) * 10000 + .5);
                        reducedFeatureSet(j, copiedFeatureCount) = (double)value / 10000;
                        copiedFeatureCount++; // use this to avoid leaving a gap in the matrix
                    }
                }
            }

            // Perform cross-validated training and push back crossvalidation performance result
            auto bestCVResultTuple = classifier->GetCrossValidationPerformanceMeasures(reducedFeatureSet, labels);
            double cvResultBalancedAccuracy = std::get<3>(bestCVResultTuple);

            CrossValidatedBalancedAccuracies.push_back(cvResultBalancedAccuracy);
        }
        //sort cross-validated balanced accuracies and pick the best elimination
        //TBD: This selection can potentially be broken out into a separate function call for flexibility
        int index = std::distance(CrossValidatedBalancedAccuracies.begin(), std::max_element(CrossValidatedBalancedAccuracies.begin(), CrossValidatedBalancedAccuracies.end()));
        CrossValidatedBalancedAccuraciesFinal.push_back(CrossValidatedBalancedAccuracies[index]);

        /*
        //check if we have three consecutive 1's 
        int sizeofcvaccuracies = CrossValidatedBalancedAccuraciesFinal.size();
        if (sizeofcvaccuracies >= 3)
        {
            if (CrossValidatedBalancedAccuraciesFinal[sizeofcvaccuracies - 1] == 1 && CrossValidatedBalancedAccuraciesFinal[sizeofcvaccuracies - 2] == 1 && CrossValidatedBalancedAccuraciesFinal[sizeofcvaccuracies - 3] == 1)
                break;
        }
        */

        std::cout << "RFE: CurrentSize=" << selectedFeatureIndices.size() << " Performance=" << CrossValidatedBalancedAccuracies[index] << std::endl;

        // eliminate the worst performing feature index from the selection
        eliminatedFeatureIndices.push_back(selectedFeatureIndices[index]);
        selectedFeatureIndices.erase(selectedFeatureIndices.begin() + index); // excise the eliminated element, shrink by 1 before repeating

    }

    // Find the final amount of features to keep based on best overall performance among feature set sizes
    //we are doing moving average to avoid local maxima. in pairs of 3, we average them and pick the middle one. 
    //TODO: consider selecting the first feature out of the three features involved in moving average
    VariableLengthVectorType MovingAverageOnCrossValidatedPerformance;
    MovingAverageOnCrossValidatedPerformance.SetSize(CrossValidatedBalancedAccuraciesFinal.size());
    for (unsigned int index = 0; index < MovingAverageOnCrossValidatedPerformance.Size(); index++)
        MovingAverageOnCrossValidatedPerformance[index] = 0;

    for (unsigned int index = 1; index < CrossValidatedBalancedAccuraciesFinal.size() - 1; index++)
        MovingAverageOnCrossValidatedPerformance[index] = (CrossValidatedBalancedAccuraciesFinal[index] + CrossValidatedBalancedAccuraciesFinal[index - 1] + CrossValidatedBalancedAccuraciesFinal[index + 1]) / 3;

    int max_performance_elimination_count = std::distance(CrossValidatedBalancedAccuraciesFinal.begin(), std::max_element(CrossValidatedBalancedAccuraciesFinal.begin(), CrossValidatedBalancedAccuraciesFinal.end()));
    int max_performance_feature_count = features.Cols() - max_performance_elimination_count;

    // according to the highest-performing elimination count, choose the final feature eliminations
    std::vector<int> finalEliminatedFeatureIndices;
    for (unsigned int index = 0; index < max_performance_elimination_count; index++)
    {
        finalEliminatedFeatureIndices.push_back(eliminatedFeatureIndices[index]);
    }

    // from the final feature elimination data, produce final remaining features
    std::vector<int> finalSelectedFeatureIndices;
    
    // For each feature...
    for (int i = 0; i < features.Cols(); i++)
    {
        bool iInEliminatedFeatures = false;
        // Consider the feature selected only if it isn't in the eliminated features.
        for (int j = 0; j < finalEliminatedFeatureIndices.size(); j++)
        {
            if (i != finalEliminatedFeatureIndices[j])
                iInEliminatedFeatures = true;
        }
        if (!iInEliminatedFeatures)
            finalSelectedFeatureIndices.push_back(i);

    }

    std::cout << "No. of selected features: " << finalSelectedFeatureIndices.size() << std::endl;
    std::cout << "Completed Backward Feature Elimination." << std::endl;
    return finalSelectedFeatureIndices;
}

