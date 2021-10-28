/** 
\file EffectSizeFeatureSelectionStrategy.cpp

\brief The file containing the EffectSizeFeatureSelectionStrategy class. 
This defines an implementation of the IFeatureSelectionStrategy interface to cover Effect-size based feature selection.
This is based on a refactoring of CaPTk TrainingModule code originally written by Saima Rathore.

Author: Alexander Getka
https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2020 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html

*/
#pragma once

#include "IFeatureSelectionStrategy.h"
#include "EffectSizeFeatureSelectionStrategy.h"



// Perform feature selection on the given data using the given classification method/strategy.
// This will pull updated parameters from params on each run. It should ONLY use params specifically related to feature selection!
std::vector<int> EffectSizeFeatureSelectionStrategy::PerformFeatureSelectionBasedTraining(
	const VariableSizeMatrixType& features, const VectorDouble& labels, IClassifierStrategy* classifier)
{
	std::cout << "Entering effect size-based feature selection." << std::endl;
	int numberOfSelectedFeatures = 0;
    std::vector<int> FinalSelectedFeatures;
    VectorDouble EffectSize = CalculateEffectSizes(features, labels);
    // Absolute value of effect size -- effect is important regardless of direction
    for (unsigned int eSizeCounter = 0; eSizeCounter < EffectSize.size(); eSizeCounter++)
    {
        if (EffectSize[eSizeCounter] < 0)
            EffectSize[eSizeCounter] = EffectSize[eSizeCounter] * -1;
    }

    // get a vector of indices into EffectSize.
    // First entry has the largest effect size, the last entry has the smallest
    std::vector<size_t> indices = sortIndexesByValue(EffectSize);

    std::vector<double> CrossValidatedBalancedAccuracies;
    std::vector<double> PerFeaturePerformance;
    
    for (unsigned int featureNumber = 0; featureNumber < EffectSize.size(); featureNumber++)
    {
        VariableSizeMatrixType currentFeatureSet;
        currentFeatureSet.SetSize(labels.size(), featureNumber + 1);

        //copy the current feature to the feature set
        for (unsigned int j = 0; j < currentFeatureSet.Rows(); j++)
            for (unsigned int k = 0; k < currentFeatureSet.Cols(); k++)
                currentFeatureSet(j, k) = features(j, indices[k]);

        // Perform cross-validated training and push back crossvalidation performance result
        classifier->TrainAndGetPerformanceMeasures(currentFeatureSet, labels);
        auto cvResultTuple = classifier->CalculatePerformanceMeasuresAgainstLabels(currentFeatureSet, labels);
        double cvResultSpecificity = std::get<3>(cvResultTuple);

        CrossValidatedBalancedAccuracies.push_back(cvResultSpecificity);
    }

    if (params.maxNumberOfFeatures > features.Cols())
    {
        params.maxNumberOfFeatures = features.Cols();
    }
    if (params.maxNumberOfFeatures > 0)
    {
        for (int index = 0; index <= params.maxNumberOfFeatures; index++)
            FinalSelectedFeatures.push_back(indices[index]);
    }
    else // Find best performance
    {
        // Find how many features to keep
        //we are doing moving average to avoid local maxima. in pairs of 3, we average them and pick the middle one. 
        //TODO: consider selecting the first feature out of the three features involved in moving average

        std::vector<double> MovingAverageOnCrossValidatedPerformance;
        MovingAverageOnCrossValidatedPerformance.push_back(0);  //can not use first index for averaging
        for (unsigned int index = 1; index < CrossValidatedBalancedAccuracies.size() - 1; index++)
            MovingAverageOnCrossValidatedPerformance.push_back((CrossValidatedBalancedAccuracies[index] + CrossValidatedBalancedAccuracies[index - 1] + CrossValidatedBalancedAccuracies[index + 1]) / 3);
        MovingAverageOnCrossValidatedPerformance.push_back(0);  //can not use last index for averaging

        int max_performance_counter = std::distance(MovingAverageOnCrossValidatedPerformance.begin(), std::max_element(MovingAverageOnCrossValidatedPerformance.begin(), MovingAverageOnCrossValidatedPerformance.end()));

        // push back the feature set that produced the best performance
        for (int index = 0; index <= max_performance_counter; index++)
            FinalSelectedFeatures.push_back(indices[index]);
    }

    std::cout << "No. of selected features: " << FinalSelectedFeatures.size() << std::endl;

	std::cout << "Completed effect size-based feature selection." << std::endl;
	return FinalSelectedFeatures;
}

// Calculates effect sizes for all features.
VectorDouble EffectSizeFeatureSelectionStrategy::CalculateEffectSizes(
	const VariableSizeMatrixType training_features, const VectorDouble target)
{
    // C1 corresponds to negatively classed samples, C2 to positively classed samples
    int NoOfSamplesC1 = 0;
    int NoOfSamplesC2 = 0;
    std::vector<double> indices_set1;
    std::vector<double> indices_set2;
    VariableSizeMatrixType features_set1;
    VariableSizeMatrixType features_set2;
    VariableLengthVectorType mean_set1;
    VariableLengthVectorType mean_set2;

    for (int index = 0; index < target.size(); index++)
    {
        if (target[index] == -1)
            NoOfSamplesC1++;
        else if (target[index] == 1)
            NoOfSamplesC2++;
    }
    features_set1.SetSize(NoOfSamplesC1, training_features.Cols());
    features_set2.SetSize(NoOfSamplesC2, training_features.Cols());
    mean_set1.SetSize(training_features.Cols());
    mean_set2.SetSize(training_features.Cols());

    NoOfSamplesC1 = 0;
    NoOfSamplesC2 = 0;
    for (unsigned int index = 0; index < target.size(); index++)
    {
        if (target[index] == -1)
        {
            for (unsigned int featureNo = 0; featureNo < training_features.Cols(); featureNo++)
                features_set1(NoOfSamplesC1, featureNo) = training_features(index, featureNo);
            NoOfSamplesC1++;
        }
        else if (target[index] == 1)
        {
            for (unsigned int featureNo = 0; featureNo < training_features.Cols(); featureNo++)
                features_set2(NoOfSamplesC2, featureNo) = training_features(index, featureNo);
            NoOfSamplesC2++;
        }
    }
    std::vector<double> EffectSize;
    for (unsigned int featureNo = 0; featureNo < training_features.Cols(); featureNo++)
    {
        // Calculate means
        double temp = 0.0;
        for (int sampleNo = 0; sampleNo < NoOfSamplesC1; sampleNo++)
            temp = temp + features_set1(sampleNo, featureNo);
        mean_set1[featureNo] = temp / NoOfSamplesC1;

        temp = 0.0;
        for (int sampleNo = 0; sampleNo < NoOfSamplesC2; sampleNo++)
            temp = temp + features_set2(sampleNo, featureNo);
        mean_set2[featureNo] = temp / NoOfSamplesC2;

        // Calculate sums of squared deviations for each class
        double sum1 = 0;
        double sum2 = 0;
        for (int sampleNo = 0; sampleNo < NoOfSamplesC1; sampleNo++)
            sum1 = sum1 + (features_set1(sampleNo, featureNo) - mean_set1[featureNo]) * (features_set1(sampleNo, featureNo) - mean_set1[featureNo]);

        for (int sampleNo = 0; sampleNo < NoOfSamplesC2; sampleNo++)
            sum2 = sum2 + (features_set2(sampleNo, featureNo) - mean_set2[featureNo]) * (features_set2(sampleNo, featureNo) - mean_set2[featureNo]);

        // Calculate effect size (needs clarification)
        double SC1 = sum1 / (NoOfSamplesC1 - 1);
        double SC2 = sum2 / (NoOfSamplesC2 - 1);
        double SP = ((NoOfSamplesC1 - 1) * SC1 + (NoOfSamplesC2 - 1) * SC2) / (NoOfSamplesC1 + NoOfSamplesC2 - 2); // Mean of squared deviations
        double currentvalue = (mean_set1[featureNo] - mean_set2[featureNo]) / sqrt(SP); // Difference in feature means between positive and negative cases, divided by the feature's RMSD
        if (std::isnan(currentvalue))
            EffectSize.push_back(0.0001); // avoid divide-by-zero while making the effect size minimal
        else
            EffectSize.push_back(currentvalue);
    }

    return EffectSize;
}


