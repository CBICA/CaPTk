/**
\file  FeatureExtractionClass.h

\brief Declaration of the FeatureExtractionClass

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

#pragma once

//#include "CAPTk.h"
#include "itkVariableSizeMatrix.h"
using VariableSizeMatrixType = itk::VariableSizeMatrix< double >;
using VectorDouble = std::vector < double >;

#define TRAINING_LABEL_NEAR 0 // [TBD] - convert to enum
#define TRAINING_LABEL_FAR 1

class  FeatureExtractionClass
{
public:
  //!Constructor
  FeatureExtractionClass();
  
  //!Destructor
  ~FeatureExtractionClass()
  {
    // do nothing
  }

  //!Returns the training data
  VariableSizeMatrixType GetTrainingData()
  {
    return mTrainingData;
  }

  //!Sets the value of training data
  void SetTrainingData(VariableSizeMatrixType &inputTrainingData)
  {
    mTrainingData = inputTrainingData;
  }

  //!Returns the test data
  VariableSizeMatrixType GetTestData()
  {
    return mTestData;
  }

  //!Sets the value of test data
  void SetTestData(VariableSizeMatrixType &inputTestData)
  {
    mTestData = inputTestData;
  }

  /**
  \brief Formulates the training data by using intensities of near and far regions, and adding corresponding labels
  \param nearintensitities Intensities of near voxels
  \param farintensitities Intensities of far voxels
  */
  void FormulateTrainingData(const std::vector< VectorDouble > &nearintensitities, const std::vector< VectorDouble > &farintensitities);
  void FormulateSurvivalTrainingData(const VariableSizeMatrixType &inputFeatures, std::vector<double> inputSurvival, VariableSizeMatrixType & SixModelFeatures, VariableSizeMatrixType & EighteenModelFeatures);
  void FormulatePseudoprogressionTrainingData(const VariableSizeMatrixType &inputFeatures, std::vector<double> inputSurvival, VariableSizeMatrixType & SixModelFeatures, VariableSizeMatrixType & EighteenModelFeatures);

  void FormulateEGFRTrainingData(const VariableSizeMatrixType &inputFeatures, std::vector<double> inputSurvival, VariableSizeMatrixType & SixModelFeatures);

  void FormulateMolecularTrainingData(const VariableSizeMatrixType &inputFeatures, std::vector<double> inputLabels, VariableSizeMatrixType & proneuralModelFeatures, VariableSizeMatrixType & neuralModelFeatures, VariableSizeMatrixType & messenchymalModelFeatures, VariableSizeMatrixType & classicalModelFeatures);

  /**
  \brief Formulates the test data by using intensities of near and far regions, and adding corresponding label
  \param testdata Intensities of test voxels
  */
  VariableSizeMatrixType FormulateTestData(const std::vector< VectorDouble > &testdata);

  /**
  \brief Resample training data to have balanced positive and negative samples
  \param trainingdata Training data
  \param nearsamples Number of near samples
  \param farsamples Number of far samples
  */
  VariableSizeMatrixType ResampleTrainingData(const VariableSizeMatrixType &trainingdata, const unsigned int nearsamples, const unsigned int farsamples);

  VariableSizeMatrixType mTrainingData;
  VariableSizeMatrixType mTestData;

};
