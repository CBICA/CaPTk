/**
\file  TrainingModule.h

\brief The header file containing the SurvivaPredictor class, used to find Survival Prediction Index 
Author: Saima Rathore
Library Dependecies: ITK 4.7+ <br>

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2016 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/
#pragma once

#include "cbicaUtilities.h"
#include "FeatureReductionClass.h"
#include "CAPTk.h"
#include "itkExtractImageFilter.h"
#include "FeatureScalingClass.h"


#ifdef APP_BASE_CAPTK_H
#include "ApplicationBase.h"
#endif

typedef itk::Image< float, 3 > ImageType;
typedef std::tuple<VectorDouble, VectorDouble, VariableSizeMatrixType, VectorDouble, VectorDouble, VariableSizeMatrixType,VectorDouble> FoldTupleType;
typedef std::map<int, FoldTupleType> MapType;

typedef itk::CSVArray2DFileReader<double> CSVFileReaderType;
typedef vnl_matrix<double> MatrixType;

/**
\class TrainingModule

\brief Calculates Survival Prediction Index 

Reference:

@article{macyszyn2015imaging,
title={Imaging patterns predict patient survival and molecular subtype in glioblastoma via machine learning techniques},
author={Macyszyn, Luke and Akbari, Hamed and Pisapia, Jared M and Da, Xiao and Attiah, Mark and Pigrish, Vadim and Bi, Yingtao and Pal, Sharmistha and Davuluri, Ramana V and Roccograndi, Laura and others},
journal={Neuro-oncology},
volume={18},
number={3},
pages={417--425},
year={2015},
publisher={Society for Neuro-Oncology}
}

*/
class TrainingModule
#ifdef APP_BASE_CAPTK_H
  : public ApplicationBase
#endif
{
public:

  cbica::Logging logger;

  //! Default constructor
  TrainingModule(){};
 ~TrainingModule(){};

  bool Run(const std::string inputFeaturesFile, const std::string inputLabelsFile, const std::string outputdirectory,const int classifierType, const int foldtype);

  std::string mEighteenTrainedFile, mSixTrainedFile;

  VectorDouble CalculatePerformanceMeasures(VariableLengthVectorType predictedLabels, std::vector<double> GivenLabels);

  bool CheckPerformanceStatus(double ist, double second, double third, double fourth, double fifth, double sixth, double seventh, double eighth, double ninth, double tenth);
  VectorDouble CalculatePerformanceMeasures(VectorDouble predictedLabels, VectorDouble GivenLabels);
  
  VectorDouble CrossValidation(const VariableSizeMatrixType inputFeatures, const VariableLengthVectorType inputLabels, const std::string outputfolder,
    const int classifiertype, const int foldtype);

  VectorDouble InternalCrossValidation(VariableSizeMatrixType inputFeatures, std::vector<double> inputLabels, double cValue, double gValue,int kerneltype);


  VectorDouble SplitTrainTest(const VariableSizeMatrixType inputFeatures, const VariableLengthVectorType inputLabels, const std::string outputfolder, const int classifiertype, const int number_of_folds, const int traiing_size);

 VectorDouble trainOpenCVSVM(const VariableSizeMatrixType &trainingDataAndLabels, const std::string &outputModelName, bool considerWeights, int ApplicationCallingSVM, double bestc, double bestg);

  VectorDouble testOpenCVSVM(const VariableSizeMatrixType &testingData, const std::string &inputModelName);

  VectorDouble CombineEstimates(const VariableLengthVectorType &estimates1, const VariableLengthVectorType &estimates2);

  VectorDouble EffectSizeFeatureSelection(const VariableSizeMatrixType training_features, std::vector<double> target);

  VectorDouble InternalCrossValidationSplitTrainTest(VariableSizeMatrixType inputFeatures, std::vector<double> inputLabels, double cValue, double gValue, int kerneltype, int counter, std::string outputfolder);
  template <typename T>
  std::vector<size_t> sort_indexes(const std::vector<T> &v);

private:

};





