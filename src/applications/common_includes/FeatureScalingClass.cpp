/**
\file  FeatureScalingClass.h

\brief Implementation of the FeatureScalingClass

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or http://www.med.upenn.edu/sbia/software/license.html

*/

#include <vtkVersion.h>
#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPCAStatistics.h"
#include "vtkStringArray.h"
#include "vtkTable.h"
#include "FeatureScalingClass.h"

FeatureScalingClass::FeatureScalingClass()
{

}

FeatureScalingClass::~FeatureScalingClass()
{
  mMeanVector.SetSize(0);
  mStdVector.SetSize(0);
}


VariableSizeMatrixType FeatureScalingClass::ScaleGivenTrainingFeatures(const VariableSizeMatrixType &inputdata)
{
  unsigned int NumberOfSamples = inputdata.Rows();
  unsigned int NumberOfFeatures = inputdata.Cols() - 1;

  //---------calculate mean and variance for each feature-----------------------------
  mMeanVector.SetSize(NumberOfFeatures);
  mStdVector.SetSize(NumberOfFeatures);
  for (unsigned int featureNo = 0; featureNo < NumberOfFeatures; featureNo++)
  {
    double temp = 0.0;
    for (unsigned int sampleNo = 0; sampleNo < NumberOfSamples; sampleNo++)
      temp += inputdata(sampleNo, featureNo);
    
    mMeanVector[featureNo] = temp / NumberOfSamples;
    double mean = mMeanVector[featureNo];
    temp = 0.0;
    for (unsigned int sampleNo = 0; sampleNo < NumberOfSamples; sampleNo++)
      temp += (inputdata(sampleNo, featureNo) - mean)*(inputdata(sampleNo, featureNo) - mean);
    
    mStdVector[featureNo] = std::sqrt(temp / (NumberOfSamples - 1));
  }
  //---------calculate z-score for each feature value-----------------------------
  VariableSizeMatrixType scaledFeatureSet;
  scaledFeatureSet.SetSize(NumberOfSamples, NumberOfFeatures + 1); //+1 to score the actual label

  for (unsigned int featureNo = 0; featureNo < NumberOfFeatures; featureNo++)
  {
    for (unsigned int sampleNo = 0; sampleNo < NumberOfSamples; sampleNo++)
      scaledFeatureSet(sampleNo, featureNo) = GetZScore(mMeanVector[featureNo], mStdVector[featureNo], inputdata(sampleNo, featureNo));
  }
  //---------------copy label to each sample----------------------------------------
  for (unsigned int sampleNo = 0; sampleNo < NumberOfSamples; sampleNo++)
    scaledFeatureSet(sampleNo, NumberOfFeatures) = inputdata(sampleNo, NumberOfFeatures);

  return scaledFeatureSet;
}

VariableSizeMatrixType FeatureScalingClass::ScaleGivenTestingFeatures(const VariableSizeMatrixType &inputdata)
{
  unsigned int  NumberOfSamples = inputdata.Rows();
  unsigned int  NumberOfFeatures = inputdata.Cols() - 1;

  //---------calculate z-score for each feature value-----------------------------
  VariableSizeMatrixType scaledFeatureSet;
  scaledFeatureSet.SetSize(NumberOfSamples, NumberOfFeatures + 1); //+1 to score the actual label

  for (unsigned int featureNo = 0; featureNo < NumberOfFeatures; featureNo++)
  {
    for (unsigned int sampleNo = 0; sampleNo < NumberOfSamples; sampleNo++)
      scaledFeatureSet(sampleNo, featureNo) = GetZScore(mMeanVector[featureNo], mStdVector[featureNo], inputdata(sampleNo, featureNo));
  }
  //---------------copy label to each sample----------------------------------------
  for (unsigned int sampleNo = 0; sampleNo < NumberOfSamples; sampleNo++)
    scaledFeatureSet(sampleNo, NumberOfFeatures) = inputdata(sampleNo, NumberOfFeatures);
  
  return scaledFeatureSet;
}
double FeatureScalingClass::GetZScore(const double &mean, const double &variance, const double &featureval)
{
  return (featureval - mean) / (variance + FLT_MIN);
}
void FeatureScalingClass::ScaleGivenTrainingFeatures(const VariableSizeMatrixType &inputdata, VariableSizeMatrixType &scaledFeatureSet, VariableLengthVectorType &meanVector, VariableLengthVectorType &stdVector)
{
  int NumberOfSamples = inputdata.Rows();
  int NumberOfFeatures = inputdata.Cols();

  //---------calculate mean and variance for each feature----------------
  meanVector.SetSize(NumberOfFeatures);
  stdVector.SetSize(NumberOfFeatures);
  for (int featureNo = 0; featureNo < NumberOfFeatures; featureNo++)
  {
    double temp = 0.0;
    for (int sampleNo = 0; sampleNo < NumberOfSamples; sampleNo++)
      temp = temp + inputdata(sampleNo, featureNo);
    meanVector[featureNo] = temp / NumberOfSamples;
    double mean = meanVector[featureNo];
    temp = 0.0;
    for (int sampleNo = 0; sampleNo < NumberOfSamples; sampleNo++)
      temp = temp + (inputdata(sampleNo, featureNo) - mean)*(inputdata(sampleNo, featureNo) - mean);
    double std = std::sqrt(temp / (NumberOfSamples - 1));
    stdVector[featureNo] = std;
  }
  //---------calculate z-score for each feature value-----------------------------
  scaledFeatureSet.SetSize(NumberOfSamples, NumberOfFeatures);
  for (int featureNo = 0; featureNo < NumberOfFeatures; featureNo++)
  {
    for (int sampleNo = 0; sampleNo < NumberOfSamples; sampleNo++)
      scaledFeatureSet(sampleNo, featureNo) = GetZScore(meanVector[featureNo], stdVector[featureNo], inputdata(sampleNo, featureNo));
  }
}
VariableSizeMatrixType FeatureScalingClass::ScaleGivenTestingFeatures(const VariableSizeMatrixType &inputdata, const VariableLengthVectorType &meandata, const VariableLengthVectorType &stddata)
{
  int NumberOfSamples = inputdata.Rows();
  int NumberOfFeatures = inputdata.Cols();

  //---------calculate z-score for each feature value-----------------------------
  VariableSizeMatrixType scaledFeatureSet;
  scaledFeatureSet.SetSize(NumberOfSamples, NumberOfFeatures); //+1 to score the actual label

  for (int featureNo = 0; featureNo < NumberOfFeatures; featureNo++)
  {
    for (int sampleNo = 0; sampleNo < NumberOfSamples; sampleNo++)
      scaledFeatureSet(sampleNo, featureNo) = GetZScore(meandata[featureNo], stddata[featureNo], inputdata(sampleNo, featureNo));
  }
  return scaledFeatureSet;
}