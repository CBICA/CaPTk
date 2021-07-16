/** 
\file ForwardFeatureSelectionStrategy.cpp

\brief The implementation file containing the ForwardFeatureSelectionStrategy class. 
This defines an implementation of the IFeatureSelectionStrategy interface to cover forward feature selection (FFS).
This is based on a refactoring of CaPTk TrainingModule code originally written by Saima Rathore.

Author: Alexander Getka
https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2020 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html

*/
#pragma once

#include "ForwardFeatureSelectionStrategy.h"



// Perform feature selection on the given data using the given classification method/strategy.
// This will pull updated parameters from params on each run. It should ONLY use params specifically related to feature selection!
std::vector<int> ForwardFeatureSelectionStrategy::PerformFeatureSelectionBasedTraining(
	const VariableSizeMatrixType& features, const VectorDouble& labels, IClassifierStrategy* classifier)
{
    std::cout << "Entering forward feature selection." << std::endl;

    // We need to tweak parameters to ignore optimization during FFS (we'll reset to these initial params afterward)
    auto initialParameters = classifier->params;
    auto paramsIgnoringOptimization = initialParameters;
    paramsIgnoringOptimization.optimizationType = CAPTK::OptimizationType::OPT_TYPE_OFF;
    paramsIgnoringOptimization.crossValidationType = CAPTK::CrossValidationType::CV_TYPE_RESUBSTITUTION;
    classifier->SetParameters(paramsIgnoringOptimization);

	std::vector<int> selectedFeatures;
	std::vector<int> unselectedFeatures = GetUnselectedFeatures(selectedFeatures, features.Cols());
    std::vector<double> CrossValidatedBalancedAccuraciesFinal;

    //feature selection mechanism
    //---------------------------
    // The idea is that we gradually select more features, one by one, evaluating performance every time we add another feature
    // Later we'll use that performance info to determine a final number of features.
    double bestCVAccuracySoFar = 0;
    while (unselectedFeatures.size() > 0)
    {
        std::vector<double> CrossValidatedBalancedAccuraciesAfterAddingUnselectedFeatures;
        
        // Iterate over each remaining feature we could add
        for (unsigned int featureNo = 0; featureNo < unselectedFeatures.size(); featureNo++)
        {
            VariableSizeMatrixType reducedFeatureSet;
            reducedFeatureSet.SetSize(labels.size(), selectedFeatures.size() + 1);

            //copy the already selected features to the current feature set
            for (unsigned int j = 0; j < reducedFeatureSet.Rows(); j++)
            {
                for (unsigned int k = 0; k < selectedFeatures.size(); k++)
                {
                    double value = (int)(features(j, selectedFeatures[k]) * 10000 + .5);
                    reducedFeatureSet(j, k) = (double)value / 10000;
                }
            }

            //copy the new feature to the current feature set
            for (unsigned int j = 0; j < reducedFeatureSet.Rows(); j++)
            {
                double value = (int)(features(j, unselectedFeatures[featureNo]) * 10000 + .5);
                reducedFeatureSet(j, reducedFeatureSet.Cols() - 1) = (double)value / 10000;
            }

            // Perform cross-validated training and push back crossvalidation performance result
            //classifier->train(reducedFeatureSet, labels); // Don't run optimization during this training
            auto bestCVResultTuple = classifier->GetCrossValidationPerformanceMeasures(reducedFeatureSet, labels);
            double cvResultBalancedAccuracy = std::get<3>(bestCVResultTuple);
            
            CrossValidatedBalancedAccuraciesAfterAddingUnselectedFeatures.push_back(cvResultBalancedAccuracy);
        }
        //sort cross-validated balanced accuracies and pick the best-performing selected feature
        int index = std::distance(CrossValidatedBalancedAccuraciesAfterAddingUnselectedFeatures.begin(), std::max_element(CrossValidatedBalancedAccuraciesAfterAddingUnselectedFeatures.begin(), CrossValidatedBalancedAccuraciesAfterAddingUnselectedFeatures.end()));

        CrossValidatedBalancedAccuraciesFinal.push_back(CrossValidatedBalancedAccuraciesAfterAddingUnselectedFeatures[index]);

        //check if we have three consecutive 1's 
        int sizeofcvaccuracies = CrossValidatedBalancedAccuraciesFinal.size();
        if (sizeofcvaccuracies >= 3)
        {
            if (CrossValidatedBalancedAccuraciesFinal[sizeofcvaccuracies - 1] == 1 && CrossValidatedBalancedAccuraciesFinal[sizeofcvaccuracies - 2] == 1 && CrossValidatedBalancedAccuraciesFinal[sizeofcvaccuracies - 3] == 1)
                break;
        }

        selectedFeatures.push_back(unselectedFeatures[index]);
        std::cout << "FFS: CurrentSize=" << selectedFeatures.size() << " Performance=" << CrossValidatedBalancedAccuraciesAfterAddingUnselectedFeatures[index] << std::endl;
        unselectedFeatures = GetUnselectedFeatures(selectedFeatures, features.Cols());

        // Only interrupt feature additions if the user actually specified a (non-zero) value 
        if ((params.maxNumberOfFeatures != 0) && (selectedFeatures.size() >= params.maxNumberOfFeatures))
        {
            // Don't waste time computing further feature selection measures
            break;
        }
    }
    // Now we evaluate the moving average of balanced accuracies and find the maximum point to determine final feature count
    //we are doing moving average to avoid local maxima. in pairs of 3, we average them and pick the middle one. 
    //TODO: consider selecting the first feature out of the three features involved in moving average

    VariableLengthVectorType MovingAverageOnCrossValidatedPerformance;
    MovingAverageOnCrossValidatedPerformance.SetSize(CrossValidatedBalancedAccuraciesFinal.size());
    for (unsigned int index = 0; index < MovingAverageOnCrossValidatedPerformance.Size(); index++)
        MovingAverageOnCrossValidatedPerformance[index] = 0;

    for (unsigned int index = 1; index < CrossValidatedBalancedAccuraciesFinal.size() - 1; index++)
        MovingAverageOnCrossValidatedPerformance[index] = (CrossValidatedBalancedAccuraciesFinal[index] + CrossValidatedBalancedAccuraciesFinal[index - 1] + CrossValidatedBalancedAccuraciesFinal[index + 1]) / 3;

    // We just find the best performing count among the features we've already selected.
    // NOTE: This section included from Saima's code -- doesn't seem to actually use the moving average.
    int max_performance_counter = std::distance(CrossValidatedBalancedAccuraciesFinal.begin(), std::max_element(CrossValidatedBalancedAccuraciesFinal.begin(), CrossValidatedBalancedAccuraciesFinal.end()));
    std::vector<int> finalSelectedFeatures;
    for (unsigned int index = 0; index <= max_performance_counter; index++)
        finalSelectedFeatures.push_back(selectedFeatures[index]);

    // Reset params
    classifier->SetParameters(initialParameters);

	std::cout << "No. of selected features: " << finalSelectedFeatures.size() << std::endl;
	std::cout << "Completed Forward Feature Selection." << std::endl;
	return finalSelectedFeatures;
}


// Return a vector of indices that are not in selectedFeatures, up to a number equal to featureSize.
std::vector<int> ForwardFeatureSelectionStrategy::GetUnselectedFeatures(const std::vector<int>& selectedFeatures, int featureSize)
{
    std::vector<int> unselectedFeatures;
    for (unsigned int featureCounter = 0; featureCounter < featureSize; featureCounter++)
    {
        bool found = false;
      
        for (unsigned int index2 = 0; index2 < selectedFeatures.size(); index2++)
        {
            if (selectedFeatures[index2] == featureCounter)
            {
                found = true;
                break;
            }
        }
        if (found == false)
            unselectedFeatures.push_back(featureCounter);
    }
    return unselectedFeatures;
}

