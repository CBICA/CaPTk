/** 
\file RandomForestFeatureSelectionStrategy.cpp

\brief The file containing the RandomForestFeatureSelectionStrategy class. 
This defines an implementation of the IFeatureSelectionStrategy interface to cover OpenCV Random Forest- based feature selection.

Author: Alexander Getka
https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2020 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html

*/
#pragma once

#include "IFeatureSelectionStrategy.h"
#include "RandomForestFeatureSelectionStrategy.h"
#include "opencv2/core/core.hpp"
#include "opencv2/ml.hpp"



// Perform feature selection on the given data using the given classification method/strategy.
// This will pull updated parameters from params on each run. It should ONLY use params specifically related to feature selection!
std::vector<int> RandomForestFeatureSelectionStrategy::PerformFeatureSelectionBasedTraining(
	const VariableSizeMatrixType& features, const VectorDouble& labels, IClassifierStrategy* classifier)
{
	std::cout << "Entering random forest based feature selection." << std::endl;
	int numberOfSelectedFeatures = 0;

    // Copy features/labels into a cv::mat
    cv::Mat trainingFeatures = cv::Mat::zeros(features.Rows(), features.Cols(), CV_32FC1);
    cv::Mat trainingLabels = cv::Mat::zeros(labels.size(), 1, CV_32F);

    for (int copyDataCounter = 0; copyDataCounter < trainingFeatures.rows; copyDataCounter++)
    {
        trainingLabels.ptr< float >(copyDataCounter)[0] = labels[copyDataCounter];
        for (int copyDataCounter2 = 0; copyDataCounter2 < trainingFeatures.cols; copyDataCounter2++)
            trainingFeatures.ptr< float >(copyDataCounter)[copyDataCounter2] = features[copyDataCounter][copyDataCounter2];
    }
    trainingLabels.convertTo(trainingLabels, CV_32SC1);

    auto trainingData = cv::ml::TrainData::create(trainingFeatures, cv::ml::ROW_SAMPLE, trainingLabels);
    auto randomForest = cv::ml::RTrees::create();

    // Reads the same params as the RF ML strategy... consult on whether this is appropriate
    auto terminator = cv::TermCriteria(cv::TermCriteria::COUNT + cv::TermCriteria::EPS, params.randomForestMaxIterations, params.randomForestEpsilon);
    
    //check setActiveVarCount parametrization here... might be necessary only for FS
    randomForest->setCalculateVarImportance(true);
    randomForest->setTermCriteria(terminator);
    randomForest->train(trainingData);
    auto importances = randomForest->getVarImportance(); // double check this
    std::vector<double> importancesSortable;

    for (int i = 0; i < importances.rows; i++)
    {
        float prediction = importances.at<float>(i, 0);
        importancesSortable.push_back(prediction);
    }
    std::vector<size_t> indices = sortIndexesByValue(importancesSortable);

    std::vector<int> FinalSelectedFeatures;
    
    unsigned int numFeaturesToSelect;
    if (params.maxNumberOfFeatures > features.Cols())
    {
        params.maxNumberOfFeatures = features.Cols();
    }
    numFeaturesToSelect = (params.maxNumberOfFeatures > 0 ? params.maxNumberOfFeatures : features.Cols());
    // push back the feature set that produced the best performance
    for (int index = 0; index < numFeaturesToSelect; index++)
        FinalSelectedFeatures.push_back(indices[index]); 
    std::cout << "No. of selected features: " << FinalSelectedFeatures.size() << std::endl;

	std::cout << "Completed random forest-based feature selection." << std::endl;
	return FinalSelectedFeatures;
}



