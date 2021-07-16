/** 
\file IClassifierStrategy.h

\brief The header file containing the IClassifierStrategy class. 
This class defines an interface that provides access to a classifier.
When adding a new classifier to the TrainingModule (such as a new kernel of SVM or a new classifier like Random Forest), 
create a new class that implements this interface.



Author: Alexander Getka
https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2020 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html

*/
#pragma once
#include "CaPTkEnums.h" 
#include "CaPTkDefines.h"
#include <string>
#include "TrainingModuleParameters.h"


class IClassifierStrategy
{

// 

public:

TrainingModuleParameters params; // This can be accessed externally to change parameters
// e.g. OpenCVSVMClassifierStrategy a; a.SetParameters(userSpecifiedParams); a.params.outputDir = "/some/path/on/disk";

// Copies parameters from a TrainingModuleParameters reference. 
// This sets the general context for execution of training and prediction. 
// If you need to over-ride some parameters in code there are 2 options:
// 1. Change the params values after-the-fact (a.params.outputDir = "/some/path/on/disk";)
// 2. Intercept the TrainingModuleParameters when you receive it, copy it, and modify it before passing.
inline void SetParameters(const TrainingModuleParameters& params)
{
	this->params = params;
}
;

// bool: true if a model has been loaded from disk or trained, false otherwise.
// Classes implementing this interface should define a variable that contains the classifier in question
// e.g. Classifiers utilizing an OpenCV svm could have a "cv::Ptr<cv::ml::SVM> _SVM" field.
// Other methods such as predict and getHyperparametersAsString will use the loaded model.
bool isModelLoaded = false;

// Loads a model from a location on disk. 
// param modelDir Leave empty (default) to use the model directory set in params.
virtual void loadModel(std::string modelDir = "") = 0;

// Saves a model to a location on disk.
// param outputDir Leave empty (default) to use the output directory set in params.
virtual void saveModel(std::string outputDir = "") = 0;

// Run prediction using the current context.
// Returns vector of prediction as type double. Caller is responsible for using this or writing to disk.
virtual VectorDouble predict(const VariableSizeMatrixType& features, bool predictDistances = false, std::string modelDir = "") = 0;

// Run training using the current context. 
// Returns true if successful, false otherwise.
// This method should pull up-to-date parameters pulled from the "params" variable on each run.
virtual bool train(const VariableSizeMatrixType& features, const VectorDouble& labels, std::string outDir = "") = 0;

// Get model info (e.g. hyperparameters) of the loaded model into a string representation (can be displayed via GUI or CLI).
// returns std::string with model information (e.g. "Optimal C = 4,\n Optimal Gamma = 3 ...") from the most recent model.
virtual std::string getModelInfoAsString() = 0;

// virtual destructor ensures any derived destructors are called (avoids surprises down the line).
virtual ~IClassifierStrategy() {};

std::tuple<double, double, double, double> CalculatePerformanceMeasuresAgainstLabels(const VariableSizeMatrixType& features, const VectorDouble& givenLabels)
{
    VectorDouble predictedLabels = this->predict(features);

    // calculate performance measures
    double TP = 0;
    double TN = 0;
    double FP = 0;
    double FN = 0;
    // Any more variables to get from this and we should probably make a map instead
    std::tuple<double, double, double, double> result;

    for (unsigned int index = 0; index < predictedLabels.size(); index++)
    {
        if (predictedLabels[index] == 1 && givenLabels[index] == 1)
            TP++;
        else if (predictedLabels[index] == -1 && givenLabels[index] == -1)
            TN++;
        else if (predictedLabels[index] == 1 && givenLabels[index] == -1)
            FP++;
        else if (predictedLabels[index] == -1 && givenLabels[index] == 1)
            FN++;
        else
        {
            //If we get here we have some kind of bug in label loading/prediction
            //Or the initial labels were provided in the wrong form
            // TODO: Break out of this appropriately
            std::cerr << "Warning: labels do not match the binary classification problem. Please submit a bug report to CBICA Software." << std::endl;
        }
    }
    double overallAccuracy = (TP + TN) / (double)predictedLabels.size();
    double sensitivity = TP / (TP + FN);
    double specificity = TN / (TN + FP);
    double balancedAccuracy = (sensitivity + specificity) / (double)2;
    result = std::make_tuple(overallAccuracy, sensitivity, specificity, balancedAccuracy);

    return result;
}

std::tuple<double, double, double, double> GetCrossValidationPerformanceMeasures(const VariableSizeMatrixType& features, const VectorDouble& labels)
{
    if (params.crossValidationType == CAPTK::CrossValidationType::CV_TYPE_RESUBSTITUTION)
    {
        this->train(features, labels);
        auto cvResult = this->CalculatePerformanceMeasuresAgainstLabels(features, labels);

        return cvResult;
    }
    else // Internal 5-fold
    {
        double folds = 5.0; //double to ensure division is fine
        std::vector<int> sampleIndices;
        for (int i = 0; i < labels.size(); i++)
        {
            sampleIndices.push_back(i);
        }
        // TODO: Run performance measures across all folds and average?
        std::random_shuffle(sampleIndices.begin(), sampleIndices.end());
        int testGroupStartIndex = 0;
        int testGroupSize = sampleIndices.size() / folds;
        int groupRemainder = sampleIndices.size() % (int)folds;
        int trainGroupSize = sampleIndices.size() - testGroupSize;
        std::vector<std::tuple<double, double, double, double>> cvResultsByFold;
        for (int currentFold = 0; currentFold < folds; currentFold++)
        { 
            VariableSizeMatrixType trainingData;
            VariableSizeMatrixType testingData;
            VectorDouble trainingLabels;
            VectorDouble testingLabels;

            trainingData.SetSize(trainGroupSize, features.Cols());
            testingData.SetSize(testGroupSize, features.Cols());

            for (int i = 0; i < testGroupSize; i++) // copy testing data samples
            {
                int testGroupIndex = i + testGroupStartIndex;
                testingLabels.push_back(labels[sampleIndices[testGroupIndex]]);

                for (int j = 0; j < features.Cols(); j++)
                {
                    testingData(i, j) = features(sampleIndices[testGroupIndex], j);
                }
            }
            for (int i = 0; i < testGroupStartIndex; i++) // copy first part of training data samples
            {
                trainingLabels.push_back(labels[i]);

                for (int j = 0; j < features.Cols(); j++)
                {
                    trainingData(i, j) = features(sampleIndices[i], j);
                }
            }
            for (int i = testGroupStartIndex+testGroupSize; i < sampleIndices.size(); i++) // copy second part of training data samples
            {
                trainingLabels.push_back(labels[i]);

                for (int j = 0; j < features.Cols(); j++)
                {
                    auto value = features(sampleIndices[i], j);
                    trainingData(i-testGroupSize, j) = value;
                }
            }
            
            this->train(trainingData, trainingLabels);
            auto cvResult = this->CalculatePerformanceMeasuresAgainstLabels(testingData, testingLabels);
            cvResultsByFold.push_back(cvResult);
            testGroupStartIndex += testGroupSize;
        }
        std::tuple<double, double, double, double> averageCVResults;
        double overallAccuracySum = 0;
        double sensitivitySum = 0;
        double specificitySum = 0;
        double balancedAccuracySum = 0;
        for (int currentFold = 0; currentFold < folds; currentFold++)
        {
            overallAccuracySum += std::get<0>(cvResultsByFold[currentFold]);
            sensitivitySum += std::get<1>(cvResultsByFold[currentFold]);
            specificitySum += std::get<2>(cvResultsByFold[currentFold]);
            balancedAccuracySum += std::get<3>(cvResultsByFold[currentFold]);
        }
        averageCVResults = std::make_tuple(overallAccuracySum / folds, 
            sensitivitySum / folds,
            specificitySum / folds, 
            balancedAccuracySum / folds);
        return averageCVResults;

    }
}


private:

// This space is reserved for implementation-specific functions and details.
// This won't be used here in the interface, but it can be in a class that implements this.

};