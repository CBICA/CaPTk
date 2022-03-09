/**
\file  TrainingModule.h

\brief The header file containing the TrainingModule class, used to build machine learning models
This class is a refactored framework based on the original TrainingModule code by Saima Rathore.
Author: Alexander Getka
Library Dependencies: ITK 4.7+ <br>

https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2020 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html

*/
#pragma once

//#include "cbicaUtilities.h"
//#include "FeatureReductionClass.h"
//#include "CAPTk.h"
#include "itkExtractImageFilter.h"
#include "itkCSVArray2DFileReader.h"
#include "FeatureScalingClass.h"
#include "CaPTkDefines.h"
#include "cbicaLogging.h"
#include "CaPTkEnums.h"
#include "TrainingModuleParameters.h"
#include "TrainingModuleResult.h"
#include <QObject>
#include <QThread>

#ifdef APP_BASE_CAPTK_H
#include "ApplicationBase.h"
#endif

typedef itk::Image< float, 3 > ImageType;
typedef std::tuple<VectorDouble, VectorDouble, VariableSizeMatrixType, VectorDouble, VectorDouble, VariableSizeMatrixType, VectorDouble,VectorDouble> FoldTupleType;
typedef std::map<int, FoldTupleType> MapType;


typedef std::tuple<VectorDouble, VectorDouble, VariableSizeMatrixType> TrainingDataTuple;
typedef std::map<int, TrainingDataTuple> TrainingMapType;

typedef std::tuple<VectorDouble, VectorDouble, VariableSizeMatrixType, VectorDouble> TestingDataTuple;
typedef std::map<int, TestingDataTuple> TestingMapType;

typedef itk::CSVArray2DFileReader<double> CSVFileReaderType;
typedef vnl_matrix<double> MatrixType;

class TrainingModule : public QObject
#ifdef APP_BASE_CAPTK_H
    : public ApplicationBase
#endif
{
   //Q_OBJECT

//private:
   //QThread* m_workerThread;

public:
  // atomic to ensure that race condition can't happen if this object is checked as it's being modified
  // isCurrentlyRunning should ALWAYS be checked before trying to access the result, otherwise it'll be very bad.
  //std::atomic_bool isCurrentlyRunning; // True if the thread is busy, false otherwise
  //TrainingModuleResult lastResult; // Automatically set when processing is finished (so we can run in another thread)

  cbica::Logging logger;



  TrainingModule();
  ~TrainingModule();

  // returns bool True if training completed without catching any errors, false otherwise
  const TrainingModuleResult Run(const TrainingModuleParameters& params);

  // Run the task in a separate thread -- useful for running in the background w/ progress bar, etc.
  // The caller takes responsibility for connecting any slots/signals necessary to make this useful. 
  //void RunThread(const TrainingModuleParameters& params);

  bool RunKFoldCrossValidation(const TrainingModuleParameters& params);
  bool RunSplitTraining(const TrainingModuleParameters& params);
  bool RunSplitTesting(const TrainingModuleParameters& params);

  void RemoveNans(VariableSizeMatrixType& mat);
  void RemoveNans(VariableLengthVectorType& vec);

  // given a set of predicted (binary-classifier) labels and actual labels, returns the balanced accuracy.
  static double GetBinaryClassificationBalancedAccuracy(VariableLengthVectorType& predicted, VariableLengthVectorType& actual);

  static void GetHeaderInformationFromFile(std::string filename);

  // param featuresMatrix : the matrix to load features into
  // returns True if successful, false otherwise
  static bool GetFeatureDataFromFile(std::string featuresFilename, VariableSizeMatrixType& featuresMatrix);

  // param labelsVector : the vector to load labels into
  // returns True if successful, false otherwise
  static bool GetLabelsFromFile(std::string labelsFilename, VariableLengthVectorType& labelsVector);

//public slots:
    //void onThreadFinished();
    //void onResultGenerated(TrainingModuleResult result);

//signals:
    //void updateProgress(int number);
    //void done(TrainingModuleResult result);

};
// Below are implementation details from the old Training Module, left only for reference purposes
/*
  VectorDouble TestData(const VariableSizeMatrixType inputFeatures, const TrainingModuleParameters& params);
  VectorDouble TrainData(const VariableSizeMatrixType inputFeatures, const VariableLengthVectorType inputLabels, 
      const TrainingModuleParameters& params);

  bool TrainData2(const VariableSizeMatrixType inputFeatures, const VariableLengthVectorType inputLabels,
    const TrainingModuleParameters& params);

  std::vector<int> UpdateUnselectedFeatures(std::vector<int> SelectedFeatures, int size);

  


  std::string mEighteenTrainedFile, mSixTrainedFile;

  VectorDouble CalculatePerformanceMeasures(VariableLengthVectorType predictedLabels, std::vector<double> GivenLabels);

  bool CheckPerformanceStatus(double ist, double second, double third, double fourth, double fifth, double sixth, double seventh, double eighth, double ninth, double tenth);
  VectorDouble CalculatePerformanceMeasures(VectorDouble predictedLabels, VectorDouble GivenLabels);

  VectorDouble CrossValidation(const VariableSizeMatrixType inputFeatures, const VariableLengthVectorType inputLabels, const TrainingModuleParameters& params);

  VectorDouble InternalCrossValidation(VariableSizeMatrixType inputFeatures, std::vector<double> inputLabels, double cValue, double gValue,int kerneltype);

  VectorDouble InternalCrossValidationResubstitution(VariableSizeMatrixType inputFeatures, std::vector<double> inputLabels, double cValue, double gValue, int kerneltype);

  VectorDouble SplitTrainTest(const VariableSizeMatrixType inputFeatures, const VariableLengthVectorType inputLabels, const TrainingModuleParameters& params, const int training_size);

  VectorDouble trainOpenCVSVM(const VariableSizeMatrixType &trainingDataAndLabels, const std::string &outputModelName, bool considerWeights, int ApplicationCallingSVM, double bestc, double bestg);

  VectorDouble testOpenCVSVM(const VariableSizeMatrixType &testingData, const std::string &inputModelName);

  VectorDouble CombineEstimates(const VariableLengthVectorType &estimates1, const VariableLengthVectorType &estimates2);

  VectorDouble EffectSizeFeatureSelection(const VariableSizeMatrixType training_features, std::vector<double> target);

  std::string CheckDataQuality(const VariableSizeMatrixType & FeaturesOfAllSubjects, const VariableLengthVectorType & LabelsOfAllSubjects);

  VectorDouble InternalCrossValidationSplitTrainTest(VariableSizeMatrixType inputFeatures, std::vector<double> inputLabels, double cValue, double gValue, int kerneltype, int counter, std::string outputfolder);
  
  std::vector<int> EffectSizeBasedFeatureSelection(const VariableSizeMatrixType inputdata, const VectorDouble labels, const TrainingModuleParameters& params, VectorDouble & crossvalidatedaccuracies);
  std::vector<int> SVMFFSBasedFeatureSelection(const VariableSizeMatrixType inputdata, const VectorDouble labels, const TrainingModuleParameters& params, 
    VectorDouble &crossvalidatedaccuracies);
//  std::vector<int> CorrelationBasedFeatureSelection(const VariableSizeMatrixType inputdata, const VectorDouble labels);


  template <typename T>
  std::vector<size_t> sort_indexes(const std::vector<T> &v);

  //void WriteCSVFiles(VariableSizeMatrixType inputdata, std::string filepath);
  //void WriteCSVFiles(VectorVectorDouble inputdata, std::string filepath);
  //void WriteCSVFiles(VariableLengthVectorType inputdata, std::string filepath);
*/

/*
template <typename T>
std::vector<size_t> TrainingModule::sort_indexes(const std::vector<T> &v)
{
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] > v[i2]; });

  return idx;
}

*/