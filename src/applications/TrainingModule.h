/**
\file  TrainingModule.h

\brief The header file containing the TrainingModule class, used to build machine learning models
Author: Saima Rathore
Library Dependecies: ITK 4.7+ <br>

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2016 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

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

class TrainingModule
#ifdef APP_BASE_CAPTK_H
  : public ApplicationBase
#endif
{
public:

  cbica::Logging logger;

  //! Default constructor
  TrainingModule() {};
  ~TrainingModule() {};

  VectorDouble TestData(const VariableSizeMatrixType inputFeatures, const std::string modelfolder, const int classifiertype, const std::string outputfolder);
  VectorDouble TrainData(const VariableSizeMatrixType inputFeatures, const VariableLengthVectorType inputLabels,
    const std::string outputfolder, const int classifiertype);

  bool TrainData2(const VariableSizeMatrixType inputFeatures, const VariableLengthVectorType inputLabels,
    const std::string outputfolder, const int classifiertype,const int featureselectiontype,
     const int optimizationType, const int crossvalidationType);

  std::vector<int> UpdateUnselectedFeatures(std::vector<int> SelectedFeatures, int size);

  bool Run(const std::string inputFeaturesFile,
    const std::string outputdirectory,
    const std::string inputLabelsFile,
    const std::string modeldirectory,
    const int classifiertype, const int foldtype,
    const int confType, const int featureselectiontype,
    const int optimizationType, const int crossvalidationType);


  std::string mEighteenTrainedFile, mSixTrainedFile;

  VectorDouble CalculatePerformanceMeasures(VariableLengthVectorType predictedLabels, std::vector<double> GivenLabels);

  bool CheckPerformanceStatus(double ist, double second, double third, double fourth, double fifth, double sixth, double seventh, double eighth, double ninth, double tenth);
  VectorDouble CalculatePerformanceMeasures(VectorDouble predictedLabels, VectorDouble GivenLabels);

  VectorDouble CrossValidation(const VariableSizeMatrixType inputFeatures, const VariableLengthVectorType inputLabels, const std::string outputfolder,
    const int classifiertype, const int foldtype,const int featureselectiontype,
    const int optimizationType, const int crossvalidationType);

  VectorDouble InternalCrossValidation(VariableSizeMatrixType inputFeatures, std::vector<double> inputLabels, double cValue, double gValue,int kerneltype);

  VectorDouble InternalCrossValidationResubstitution(VariableSizeMatrixType inputFeatures, std::vector<double> inputLabels, double cValue, double gValue, int kerneltype);

  VectorDouble SplitTrainTest(const VariableSizeMatrixType inputFeatures, const VariableLengthVectorType inputLabels, const std::string outputfolder, const int classifiertype, const int training_size);

  VectorDouble trainOpenCVSVM(const VariableSizeMatrixType &trainingDataAndLabels, const std::string &outputModelName, bool considerWeights, int ApplicationCallingSVM, double bestc, double bestg);

  VectorDouble testOpenCVSVM(const VariableSizeMatrixType &testingData, const std::string &inputModelName);

  VectorDouble CombineEstimates(const VariableLengthVectorType &estimates1, const VariableLengthVectorType &estimates2);

  VectorDouble EffectSizeFeatureSelection(const VariableSizeMatrixType training_features, std::vector<double> target);

  std::string CheckDataQuality(const VariableSizeMatrixType & FeaturesOfAllSubjects, const VariableLengthVectorType & LabelsOfAllSubjects);

  VectorDouble InternalCrossValidationSplitTrainTest(VariableSizeMatrixType inputFeatures, std::vector<double> inputLabels, double cValue, double gValue, int kerneltype, int counter, std::string outputfolder);
  
  std::vector<int> EffectSizeBasedFeatureSelection(const VariableSizeMatrixType inputdata, const VectorDouble labels, const int classifiertype, const int optimizationtype, const int cvtype, VectorDouble & crossvalidatedaccuracies);
  std::vector<int> SVMFFSBasedFeatureSelection(const VariableSizeMatrixType inputdata, const VectorDouble labels, const int classifiertype, 
    const int optimizationtype, const int cvtype, 
    VectorDouble &crossvalidatedaccuracies);
//  std::vector<int> CorrelationBasedFeatureSelection(const VariableSizeMatrixType inputdata, const VectorDouble labels);


  template <typename T>
  std::vector<size_t> sort_indexes(const std::vector<T> &v);

  //void WriteCSVFiles(VariableSizeMatrixType inputdata, std::string filepath);
  //void WriteCSVFiles(VectorVectorDouble inputdata, std::string filepath);
  //void WriteCSVFiles(VariableLengthVectorType inputdata, std::string filepath);

private:

};


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