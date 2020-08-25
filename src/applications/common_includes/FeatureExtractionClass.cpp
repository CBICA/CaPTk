/**
\file  FeatureExtractionClass.cpp

\brief Implementation of the FeatureExtractionClass

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

#include "FeatureExtractionClass.h"
//#include "itkVariableSizeMatrix.h"
//#include "itkVariableLengthVector.h"
//#include "CAPTk.h"
#include "CaPTkEnums.h"

FeatureExtractionClass::FeatureExtractionClass()
{

}

void FeatureExtractionClass::FormulateSurvivalTrainingData(VectorDouble inputSurvival, VectorDouble & SixModelLabels, VectorDouble & EighteenModelLabels)
{
  for (unsigned int i = 0; i < inputSurvival.size(); i++)
  {
    if (inputSurvival[i] <= 6)
      SixModelLabels.push_back(1);
    else
      SixModelLabels.push_back(-1);

    if (inputSurvival[i] >= 18)
      EighteenModelLabels.push_back(1);
    else
      EighteenModelLabels.push_back(-1);
  }
}

VariableSizeMatrixType FeatureExtractionClass::ResampleTrainingData(const VariableSizeMatrixType &trainingdata, const unsigned int NumberOfNearSamples, const unsigned int NumberOfFarSamples)
{
  VariableSizeMatrixType sampledTrainingData;
  int NumberOfFeatures = trainingdata.Cols() - 1;
  unsigned int stepSize = NumberOfFarSamples / NumberOfNearSamples;
  int NumberOfSampledFarSamples = 0;
  unsigned int TotalSamples = NumberOfFarSamples + NumberOfNearSamples;

  for (unsigned int i = 0; i < NumberOfFarSamples; i = i + stepSize)
    NumberOfSampledFarSamples++;

  sampledTrainingData.SetSize(NumberOfNearSamples + NumberOfSampledFarSamples, NumberOfFeatures + 1);

  for (unsigned int i = 0; i < NumberOfNearSamples; i++)
  {
    for (unsigned int j = 0; j < trainingdata.Cols(); j++)
      sampledTrainingData(i, j) = trainingdata(i, j);
  }
  int counter = NumberOfNearSamples;

  for (unsigned int i = NumberOfNearSamples; i < TotalSamples; i = i + stepSize)
  {
    for (unsigned int j = 0; j < trainingdata.Cols(); j++)
      sampledTrainingData(counter, j) = trainingdata(i, j);
    
    counter++;
  }
  return sampledTrainingData;
}


void FeatureExtractionClass::FormulateTrainingData(const std::vector< VectorDouble > &nearintensitities, const std::vector< VectorDouble > &farintensitities)
{
  size_t NumberOfNearSamples = nearintensitities.size();
  size_t NumberOfFarSamples = farintensitities.size();
  size_t NumberOfFeatures = nearintensitities[0].size();

  mTrainingData.SetSize(NumberOfNearSamples + NumberOfFarSamples, NumberOfFeatures + 1);
  for (size_t i = 0; i < NumberOfNearSamples; i++)
  {
    for (size_t j = 0; j < NumberOfFeatures; j++)
      mTrainingData(i, j) = nearintensitities[i][j];

    mTrainingData(i, NumberOfFeatures) = TRAINING_LABEL_NEAR;
  }
  for (size_t i = 0; i < NumberOfFarSamples; i++)
  {
    for (size_t j = 0; j < NumberOfFeatures; j++)
      mTrainingData(i + NumberOfNearSamples, j) = farintensitities[i][j];

    mTrainingData(i + NumberOfNearSamples, NumberOfFeatures) = TRAINING_LABEL_FAR;
  }
}


VariableSizeMatrixType FeatureExtractionClass::FormulateTestData(const std::vector< VectorDouble > &testdata)
{
  size_t NumberOfFeatures = testdata[0].size();
  size_t NumberOfSamples = testdata.size();

  mTestData.SetSize(NumberOfSamples, NumberOfFeatures + 1);

  for (size_t i = 0; i < NumberOfSamples; i++)
  {
    for (size_t j = 0; j < NumberOfFeatures; j++)
      mTestData(i, j) = testdata[i][j];

    mTestData(i, NumberOfFeatures) = TRAINING_LABEL_FAR;
  }

  return mTestData;
}

void FeatureExtractionClass::FormulatePseudoprogressionTrainingData(std::vector<double> inputLabels, VectorDouble & PseudoModelLabels, VectorDouble & RecurrenceModelLabels)
{
  for (unsigned int i = 0; i < inputLabels.size(); i++)
  {
    if (inputLabels[i] == 1 || inputLabels[i]==2)
      PseudoModelLabels.push_back(1);
    else
      PseudoModelLabels.push_back(-1);

    if (inputLabels[i] == 5 || inputLabels[i] == 6)
      RecurrenceModelLabels.push_back(1);
    else
      RecurrenceModelLabels.push_back(-1);
  }
}




void FeatureExtractionClass::FormulateEGFRTrainingData(const VariableSizeMatrixType &inputFeatures, std::vector<double> inputSurvival, VariableSizeMatrixType & SixModelFeatures)
{
  std::vector<int> posIndices;
  std::vector<int> negIndices;

  for (unsigned int i = 0; i < inputSurvival.size(); i++)
  {
    if (inputSurvival[i] ==1)
      posIndices.push_back(i);
    else
      negIndices.push_back(i);
  }
  SixModelFeatures.SetSize(inputFeatures.Rows(), inputFeatures.Cols() + 1);

  for (unsigned int i = 0; i < posIndices.size(); i++)
  {
    unsigned int j = 0;
    for (j = 0; j < inputFeatures.Cols(); j++)
      SixModelFeatures(i, j) = inputFeatures(posIndices[i], j);
    SixModelFeatures(i, j) = 0;
  }
  for (unsigned int i = 0; i < negIndices.size(); i++)
  {
    unsigned int j = 0;
    for (j = 0; j < inputFeatures.Cols(); j++)
      SixModelFeatures(i + posIndices.size(), j) = inputFeatures(negIndices[i], j);
    SixModelFeatures(i + posIndices.size(), j) = 1;
  }
}

void FeatureExtractionClass::FormulateMolecularTrainingData(VectorDouble inputLabels, 
  VectorDouble & proneuralModelLabels, VectorDouble & neuralModelLabels,
  VectorDouble & messModelLabels, VectorDouble & classicalModelLabels)
{
  for (unsigned int i = 0; i < inputLabels.size(); i++)
  {
    if (inputLabels[i] == CAPTK::MOLECULAR_SUBTYPES::PRONEURAL)
      proneuralModelLabels.push_back(1);
    else
      proneuralModelLabels.push_back(-1);

    if (inputLabels[i] == CAPTK::MOLECULAR_SUBTYPES::NEURAL)
      neuralModelLabels.push_back(1);
    else
      neuralModelLabels.push_back(-1);

    if (inputLabels[i] == CAPTK::MOLECULAR_SUBTYPES::MESSENCHYMAL)
      messModelLabels.push_back(1);
    else
      messModelLabels.push_back(-1);

    if (inputLabels[i] == CAPTK::MOLECULAR_SUBTYPES::CLASSICAL)
      classicalModelLabels.push_back(1);
    else
      classicalModelLabels.push_back(-1);
  }
}



//
//
//void FeatureExtractionClass::FormulateMolecularTrainingData(const VariableSizeMatrixType &inputFeatures, std::vector<double> inputLabels, VariableSizeMatrixType & proneuralModelFeatures, VariableSizeMatrixType & neuralModelFeatures, VariableSizeMatrixType & messenchymalModelFeatures, VariableSizeMatrixType & classicalModelFeatures)
//{
//	std::vector<int> proneuralModelIndices;
//	std::vector<int> neuralModelIndices;
//	std::vector<int> messenchymalModelIndices;
//	std::vector<int> classicalModelIndices;
//
//	for (unsigned int i = 0; i < inputLabels.size(); i++)
//	{
//		if (inputLabels[i] == CAPTK::MOLECULAR_SUBTYPES::PRONEURAL)
//			proneuralModelIndices.push_back(i);
//		else if (inputLabels[i] == CAPTK::MOLECULAR_SUBTYPES::NEURAL)
//			neuralModelIndices.push_back(i);
//		else if (inputLabels[i] == CAPTK::MOLECULAR_SUBTYPES::MESSENCHYMAL)
//			messenchymalModelIndices.push_back(i);
//		else
//			classicalModelIndices.push_back(i);
//	}
//	classicalModelFeatures.SetSize(inputFeatures.Rows(), inputFeatures.Cols() + 1);
//	proneuralModelFeatures.SetSize(inputFeatures.Rows(), inputFeatures.Cols() + 1);
//	neuralModelFeatures.SetSize(inputFeatures.Rows(), inputFeatures.Cols() + 1);
//	messenchymalModelFeatures.SetSize(inputFeatures.Rows(), inputFeatures.Cols() + 1);
//
//
//	for (unsigned int i = 0; i < inputFeatures.Rows(); i++)
//	{
//		unsigned int j = 0;
//		for (j = 0; j < inputFeatures.Cols(); j++)
//		{
//			classicalModelFeatures(i, j) = inputFeatures(i, j);
//			neuralModelFeatures(i, j) = inputFeatures(i, j);
//			proneuralModelFeatures(i, j) = inputFeatures(i, j);
//			messenchymalModelFeatures(i, j) = inputFeatures(i, j);
//		}
//		if (inputLabels[i] == CAPTK::MOLECULAR_SUBTYPES::PRONEURAL)
//		{
//			proneuralModelFeatures(i, j) = 1;
//			neuralModelFeatures(i, j) = 0;
//			messenchymalModelFeatures(i, j) = 0;
//			classicalModelFeatures(i, j) = 0;
//		}
//		else if (inputLabels[i] == CAPTK::MOLECULAR_SUBTYPES::NEURAL)
//		{
//			proneuralModelFeatures(i, j) = 0;
//			neuralModelFeatures(i, j) = 1;
//			messenchymalModelFeatures(i, j) = 0;
//			classicalModelFeatures(i, j) = 0;
//		}
//		else if (inputLabels[i] == CAPTK::MOLECULAR_SUBTYPES::CLASSICAL)
//		{
//			proneuralModelFeatures(i, j) = 0;
//			neuralModelFeatures(i, j) = 0;
//			messenchymalModelFeatures(i, j) = 0;
//			classicalModelFeatures(i, j) = 1;
//		}
//		else
//		{
//			proneuralModelFeatures(i, j) = 0;
//			neuralModelFeatures(i, j) = 0;
//			messenchymalModelFeatures(i, j) = 1;
//			classicalModelFeatures(i, j) = 0;
//		}
//	}
//}
