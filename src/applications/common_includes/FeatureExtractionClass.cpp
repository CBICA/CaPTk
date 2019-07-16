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

void FeatureExtractionClass::FormulateSurvivalTrainingData(const VariableSizeMatrixType &inputFeatures, std::vector<double> inputSurvival, VariableSizeMatrixType & SixModelFeatures, VariableSizeMatrixType & EighteenModelFeatures)
{
  std::vector<int> SixModelLowerIndices;
  std::vector<int> SixModelHigherIndices;
  std::vector<int> EighteenModelLowerIndices;
  std::vector<int> EighteenModelHigherIndices;

  for (unsigned int i = 0; i < inputSurvival.size(); i++)
  {
    if (inputSurvival[i] <= 6)
      SixModelLowerIndices.push_back(i);
    else
      SixModelHigherIndices.push_back(i);

    if (inputSurvival[i] <= 18)
      EighteenModelLowerIndices.push_back(i);
    else
      EighteenModelHigherIndices.push_back(i);
  }
  SixModelFeatures.SetSize(inputFeatures.Rows(), inputFeatures.Cols() + 1);
  EighteenModelFeatures.SetSize(inputFeatures.Rows(), inputFeatures.Cols() + 1);

  for (unsigned int i = 0; i < SixModelLowerIndices.size(); i++)
  {
    unsigned int j = 0;
    for (j = 0; j < inputFeatures.Cols(); j++)
      SixModelFeatures(i, j) = inputFeatures(SixModelLowerIndices[i], j);
    SixModelFeatures(i, j) = 0;
  }
  for (unsigned int i = 0; i < SixModelHigherIndices.size(); i++)
  {
    unsigned int j = 0;
    for (j = 0; j < inputFeatures.Cols(); j++)
      SixModelFeatures(i + SixModelLowerIndices.size(), j) = inputFeatures(SixModelHigherIndices[i], j);
    SixModelFeatures(i + SixModelLowerIndices.size(), j) = 1;
  }

  for (unsigned int i = 0; i < EighteenModelLowerIndices.size(); i++)
  {
    unsigned int j = 0;
    for (j = 0; j < inputFeatures.Cols(); j++)
      EighteenModelFeatures(i, j) = inputFeatures(EighteenModelLowerIndices[i], j);
    EighteenModelFeatures(i, j) = 0;
  }
  for (unsigned int i = 0; i < EighteenModelHigherIndices.size(); i++)
  {
    unsigned int j = 0;
    for (j = 0; j < inputFeatures.Cols(); j++)
      EighteenModelFeatures(i + EighteenModelLowerIndices.size(), j) = inputFeatures(EighteenModelHigherIndices[i], j);
    EighteenModelFeatures(i + EighteenModelLowerIndices.size(), j) = 1;
  }
}


void FeatureExtractionClass::FormulatePseudoprogressionTrainingData(const VariableSizeMatrixType &inputFeatures, std::vector<double> inputSurvival, VariableSizeMatrixType & PseudoModelFeatures, VariableSizeMatrixType & RecurrenceModelFeatures)
{
  std::vector<int> PseudoModelLowerIndices;
  std::vector<int> PseudoModelHigherIndices;
  std::vector<int> RecurrenceModelLowerIndices;
  std::vector<int> RecurrenceModelHigherIndices;

  for (unsigned int i = 0; i < inputSurvival.size(); i++)
  {
    if (inputSurvival[i] ==1 || inputSurvival[i] == 2)
      PseudoModelLowerIndices.push_back(i);
    else
      PseudoModelHigherIndices.push_back(i);

    if (inputSurvival[i] == 5 || inputSurvival[i] == 6)
      RecurrenceModelLowerIndices.push_back(i);
    else
      RecurrenceModelHigherIndices.push_back(i);
  }
  PseudoModelFeatures.SetSize(inputFeatures.Rows(), inputFeatures.Cols() + 1);
  RecurrenceModelFeatures.SetSize(inputFeatures.Rows(), inputFeatures.Cols() + 1);

  for (unsigned int i = 0; i < PseudoModelLowerIndices.size(); i++)
  {
    unsigned int j = 0;
    for (j = 0; j < inputFeatures.Cols(); j++)
      PseudoModelFeatures(i, j) = inputFeatures(PseudoModelLowerIndices[i], j);
    PseudoModelFeatures(i, j) = 0;
  }
  for (unsigned int i = 0; i < PseudoModelHigherIndices.size(); i++)
  {
    unsigned int j = 0;
    for (j = 0; j < inputFeatures.Cols(); j++)
      PseudoModelFeatures(i + PseudoModelLowerIndices.size(), j) = inputFeatures(PseudoModelHigherIndices[i], j);
    PseudoModelFeatures(i + PseudoModelLowerIndices.size(), j) = 1;
  }

  for (unsigned int i = 0; i < RecurrenceModelLowerIndices.size(); i++)
  {
    unsigned int j = 0;
    for (j = 0; j < inputFeatures.Cols(); j++)
      RecurrenceModelFeatures(i, j) = inputFeatures(RecurrenceModelLowerIndices[i], j);
    RecurrenceModelFeatures(i, j) = 0;
  }
  for (unsigned int i = 0; i < RecurrenceModelHigherIndices.size(); i++)
  {
    unsigned int j = 0;
    for (j = 0; j < inputFeatures.Cols(); j++)
      RecurrenceModelFeatures(i + RecurrenceModelLowerIndices.size(), j) = inputFeatures(RecurrenceModelHigherIndices[i], j);
    RecurrenceModelFeatures(i + RecurrenceModelLowerIndices.size(), j) = 1;
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






void FeatureExtractionClass::FormulateMolecularTrainingData(const VariableSizeMatrixType &inputFeatures, std::vector<double> inputLabels, VariableSizeMatrixType & proneuralModelFeatures, VariableSizeMatrixType & neuralModelFeatures, VariableSizeMatrixType & messenchymalModelFeatures, VariableSizeMatrixType & classicalModelFeatures)
{
	std::vector<int> proneuralModelIndices;
	std::vector<int> neuralModelIndices;
	std::vector<int> messenchymalModelIndices;
	std::vector<int> classicalModelIndices;

	for (unsigned int i = 0; i < inputLabels.size(); i++)
	{
		if (inputLabels[i] == CAPTK::MOLECULAR_SUBTYPES::PRONEURAL)
			proneuralModelIndices.push_back(i);
		else if (inputLabels[i] == CAPTK::MOLECULAR_SUBTYPES::NEURAL)
			neuralModelIndices.push_back(i);
		else if (inputLabels[i] == CAPTK::MOLECULAR_SUBTYPES::MESSENCHYMAL)
			messenchymalModelIndices.push_back(i);
		else
			classicalModelIndices.push_back(i);
	}
	classicalModelFeatures.SetSize(inputFeatures.Rows(), inputFeatures.Cols() + 1);
	proneuralModelFeatures.SetSize(inputFeatures.Rows(), inputFeatures.Cols() + 1);
	neuralModelFeatures.SetSize(inputFeatures.Rows(), inputFeatures.Cols() + 1);
	messenchymalModelFeatures.SetSize(inputFeatures.Rows(), inputFeatures.Cols() + 1);


	for (unsigned int i = 0; i < inputFeatures.Rows(); i++)
	{
		unsigned int j = 0;
		for (j = 0; j < inputFeatures.Cols(); j++)
		{
			classicalModelFeatures(i, j) = inputFeatures(i, j);
			neuralModelFeatures(i, j) = inputFeatures(i, j);
			proneuralModelFeatures(i, j) = inputFeatures(i, j);
			messenchymalModelFeatures(i, j) = inputFeatures(i, j);
		}
		if (inputLabels[i] == CAPTK::MOLECULAR_SUBTYPES::PRONEURAL)
		{
			proneuralModelFeatures(i, j) = 1;
			neuralModelFeatures(i, j) = 0;
			messenchymalModelFeatures(i, j) = 0;
			classicalModelFeatures(i, j) = 0;
		}
		else if (inputLabels[i] == CAPTK::MOLECULAR_SUBTYPES::NEURAL)
		{
			proneuralModelFeatures(i, j) = 0;
			neuralModelFeatures(i, j) = 1;
			messenchymalModelFeatures(i, j) = 0;
			classicalModelFeatures(i, j) = 0;
		}
		else if (inputLabels[i] == CAPTK::MOLECULAR_SUBTYPES::CLASSICAL)
		{
			proneuralModelFeatures(i, j) = 0;
			neuralModelFeatures(i, j) = 0;
			messenchymalModelFeatures(i, j) = 0;
			classicalModelFeatures(i, j) = 1;
		}
		else
		{
			proneuralModelFeatures(i, j) = 0;
			neuralModelFeatures(i, j) = 0;
			messenchymalModelFeatures(i, j) = 1;
			classicalModelFeatures(i, j) = 0;
		}
	}
}
