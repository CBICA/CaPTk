#pragma once

#include "CaPTkDefines.h"
#include "CaPTkEnums.h"
#include "opencv2/core/core.hpp"
#include "opencv2/ml.hpp"
/**
\brief Train and save the SVM classifier

\param trainingDataAndLabels Input training data with last column as training labels
\param outputModelName File to save the model
*/
inline VectorDouble trainOpenCVSVM(const VariableSizeMatrixType &trainingDataAndLabels, const std::string &outputModelName, bool considerWeights, int ApplicationCallingSVM)
{
  cv::Mat trainingData = cv::Mat::zeros(trainingDataAndLabels.Rows(), trainingDataAndLabels.Cols() - 1, CV_32FC1),
    trainingLabels = cv::Mat::zeros(trainingDataAndLabels.Rows(), 1, CV_32FC1);

  //// debugging purposes
  //std::ofstream file;
  //std::string base, path, ext;
  //cbica::splitFileName(outputModelName, path, base, ext);
  //file.open((path + base + "_trainingAndLabels.csv").c_str());
  //for (size_t i = 0; i < trainingDataAndLabels.Rows(); i++)
  //{
  //  for (size_t j = 0; j < trainingDataAndLabels.Cols(); j++)
  //  {
  //    file << trainingDataAndLabels(i, j) << ",";
  //  }
  //  file << "\n";
  //}
  //file.close();

  cv::Mat trainingWeights = cv::Mat::zeros(2, 1, CV_32FC1);
  size_t label1Counter = 0, label2Counter = 0;
  // fast cv::Mat access 
  for (unsigned int i = 0; i < trainingDataAndLabels.Rows(); ++i)
  {
    for (unsigned int j = 0; j < trainingDataAndLabels.Cols() - 1; ++j)
    {
      trainingData.ptr< float >(i)[j] = trainingDataAndLabels(i, j);
    }

    // last column consists of labels
    trainingLabels.ptr< float >(i)[0] = trainingDataAndLabels(i, trainingDataAndLabels.Cols() - 1);
    if (considerWeights)
    {
      // start counter to assign weights
      if (trainingLabels.ptr< float >(i)[0] == 0)
      {
        label1Counter++;
      }
      else
      {
        label2Counter++;
      }
    }
  }

  // slow cv::Mat access 
  //for (size_t i = 0; i < trainingDataAndLabels.Rows(); i++)
  //{
  //  for (size_t j = 0; j < trainingDataAndLabels.Cols() - 1; j++)
  //  {
  //    trainingData.at< float >(i, j) = trainingDataAndLabels(i, j);
  //  }
  //
  //  // last column consists of labels
  //  trainingLabels.at< float >(i,0) = trainingDataAndLabels(i, trainingDataAndLabels.Cols() - 1);
  //
  //  if (considerWeights)
  //  {
  //    // start counter to assign weights
  //    if (trainingDataAndLabels(i, trainingDataAndLabels.Cols() - 1) == 0)
  //    {
  //      label1Counter++;
  //    }
  //    else
  //    {
  //      label2Counter++;
  //    }
  //  }
  //}
  //cv::Mat diff;
  //cv::compare(trainingData, trainingData_temp, diff, cv::CMP_NE);
  //int nz = cv::countNonZero(diff);

  //cv::compare(trainingLabels, trainingLabels_temp, diff, cv::CMP_NE);
  //int nz2 = cv::countNonZero(diff);

  trainingLabels.convertTo(trainingLabels, CV_32SC1);

  auto svm = cv::ml::SVM::create();
  svm->setType(cv::ml::SVM::C_SVC);
  //svm->setC(1);
  //svm->setGamma(0.01);
  // parameters for auto_train
  cv::ml::ParamGrid grid_C(-5, 5, 2); // Parameter C of a SVM optimization problem (C_SVC / EPS_SVR / NU_SVR)
  cv::ml::ParamGrid grid_Gamma(-5, 5, 2); // Parameter \gamma of a kernel function (POLY / RBF / SIGMOID / CHI2)
  cv::ml::ParamGrid grid_P(-5, 5, 2); // Parameter \epsilon of a SVM optimization problem (EPS_SVR)
  cv::ml::ParamGrid grid_Nu(-5, 5, 2); // Parameter \nu of a SVM optimization problem (NU_SVC / ONE_CLASS / NU_SVR)
  cv::ml::ParamGrid grid_Degree(0, 5, 1); // Parameter degree of a kernel function (POLY)
  cv::ml::ParamGrid grid_Coeff0(-5, 5, 2); // this is Parameter coef0 of a kernel function (POLY / SIGMOID)
  if (ApplicationCallingSVM == CAPTK::ApplicationCallingSVM::Recurrence)
  {
    svm->setKernel(cv::ml::SVM::RBF); // using this produces terrible results for recurrence
  }
  else
  {
    svm->setKernel(cv::ml::SVM::LINEAR);
  }
  svm->setC(1);
  svm->setGamma(0.01);
  //svm->setKernel(cv::ml::SVM::POLY); // this crashes

  // assign penalties
  //if (considerWeights)
  //{
  //  trainingWeights.ptr< float >(0)[0] = label1Counter;
  //  trainingWeights.ptr< float >(1)[0] = label2Counter;
  //  //trainingWeights.at< float >(0, 0) = label1Counter;
  //  //trainingWeights.at< float >(1, 0) = label2Counter;
  //  svm->setClassWeights(trainingWeights);
  //}
  bool res = true;
  std::string msg;

  try
  {
    res = svm->train(trainingData, cv::ml::ROW_SAMPLE, trainingLabels);
    //res = svm->trainAuto(cv::ml::TrainData::create(trainingData, cv::ml::ROW_SAMPLE, trainingLabels), 
    // 10, grid_C, grid_Gamma, grid_P, grid_Nu, grid_Coeff0, grid_Degree, true);
  }
  catch (cv::Exception ex)
  {
    msg = ex.what();
  }
  //---------------------------------find distances on training adat----------------------------------
  VectorDouble returnVec;
  returnVec.resize(trainingData.rows);
  cv::Mat predicted(1, 1, CV_32F); // not sure if this is required if we do train_auto for cross-validation
  if (ApplicationCallingSVM == CAPTK::ApplicationCallingSVM::Recurrence) // only recurrence needs the distance map
  {
    for (int i = 0; i < trainingData.rows; i++)
    {
      cv::Mat sample = trainingData.row(i);
      //float value = svm->predict(sample, cv::Mat(), cv::ml::StatModel::RAW_OUTPUT);`
      /*float value = */svm->predict(sample, predicted, true /*cv::ml::StatModel::RAW_OUTPUT*/);
      returnVec[i] = /*predicted.at<float>(0, 0)*/predicted.ptr< float >(0)[0];
    }
  }
  else
  {
    returnVec.push_back(1.0); // just so survival application doesn't thrown an error
  }
  //--------------------------------------------------------------------------------------------------
  //svm->setTermCriteria(cv::TermCriteria(cv::TermCriteria::MAX_ITER, (int)1e7, 1e-6));
  if (res)
  {
    //cv::Mat sv = svm->getUncompressedSupportVectors();
    svm->save(outputModelName);
  }
  else
  {
    //return false;
  }

  return returnVec;
}

/**
\brief Load the SVM classifier and get the distances from the hyperplane

\param testingData Input training data with last column as training labels
\param inputModelName File to save the model
\return Distances of classification
*/
inline VectorDouble testOpenCVSVM(const VariableSizeMatrixType &testingData, const std::string &inputModelName)
{
  auto svm = cv::Algorithm::load<cv::ml::SVM>(inputModelName);

  //std::ofstream file;
  //file.open("Z:/Projects/testingData.csv");
  //for (size_t i = 0; i < testingData.Rows(); i++)
  //{
  //  for (size_t j = 0; j < testingData.Cols(); j++)
  //  {
  //    file << testingData(i, j) << ",";
  //  }
  //  file << "\n";
  //}
  //file.close();

  VectorDouble returnVecScore;
  VectorDouble returnVecLabel;
  cv::Mat testingDataMat = cv::Mat::zeros(testingData.Rows(), testingData.Cols(), CV_32FC1), outputProbs;

  // fast cv::Mat access
  for (unsigned int i = 0; i < testingData.Rows(); i++)
  {
    for (unsigned int j = 0; j < testingData.Cols(); j++)
    {
      testingDataMat.ptr< float >(i)[j] = testingData(i, j);
    }
  }

  // slow cv::Mat access
  //for (size_t i = 0; i < testingData.Rows(); i++)
  //{
  //  for (size_t j = 0; j < testingData.Cols() - 1; j++)
  //  {
  //    testingDataMat.at< float >(i, j) = testingData(i, j);
  //  }
  //}

  // see http://docs.opencv.org/trunk/db/d7d/classcv_1_1ml_1_1StatModel.html#af1ea864e1c19796e6264ebb3950c0b9a for details regarding why '1'
  //svm->predict(testingDataMat, returnVec, 1);

  returnVecScore.resize(testingDataMat.rows);
  returnVecLabel.resize(testingDataMat.rows);
  //this segment of code iterates through all the test samples and assigns predicted scores
  cv::Mat predicted(1, 1, CV_32F);
  for (int i = 0; i < testingDataMat.rows; i++)
  {
    cv::Mat sample = testingDataMat.row(i);
    svm->predict(sample, predicted, true/*cv::ml::StatModel::RAW_OUTPUT*/);
    returnVecScore[i] = predicted.ptr< float >(0)[0];
  }
  //this segment of code iterates through all the test samples and assigns predicted labels
  for (int i = 0; i < testingDataMat.rows; i++)
  {
    cv::Mat sample = testingDataMat.row(i);
    svm->predict(sample, predicted, false/*cv::ml::StatModel::RAW_OUTPUT*/);
    returnVecLabel[i] = predicted.ptr< float >(0)[0];
  }
  for (int i = 0; i < returnVecLabel.size(); i++)
    returnVecScore[i]= std::abs(returnVecScore[i])*returnVecLabel[i];

  return returnVecScore;
}

/**
\brief Load the SVM classifier and predict the probabilities

\param testingData Input training data with last column as training labels
\param inputModelName File to save the model
\return Probabilities of classification
*/
//inline VectorDouble testOpenCVSVM_Probs(const VariableSizeMatrixType &testingData, const std::string &inputModelName)
//{
//  auto svm = cv::Algorithm::load<cv::ml::SVM>(inputModelName);
//
//  VectorDouble returnVec;
//  returnVec.resize(testingData.Rows());
//  cv::Mat testingDataMat = cv::Mat::zeros(testingData.Rows(), testingData.Cols(), CV_32FC1), outputProbs;
//
//  //// fast cv::Mat access -- needs to be checked
//  //float *p;
//  //for (int i = 0; i < testingData.Rows(); ++i)
//  //{
//  //  p = testingDataMat.ptr< float >(i);
//  //  for (int j = 0; j < testingData.Cols(); ++j)
//  //  {
//  //    p[j] = testingData(i, j);
//  //  }
//  //}
//
//  for (size_t i = 0; i < testingData.Rows(); i++)
//  {
//    for (size_t j = 0; j < testingData.Cols(); j++)
//    {
//      testingDataMat.at< float >(i, j) = testingData(i, j);
//    }
//  }
//
//  // see http://docs.opencv.org/trunk/db/d7d/classcv_1_1ml_1_1StatModel.html#af1ea864e1c19796e6264ebb3950c0b9a for details regarding why '1'
//  svm->predict(testingDataMat, returnVec, 1);
//
//  float max_1, max_2;
//  for (size_t i = 0; i < testingDataMat.rows; i++)
//  {
//    float temp;
//    //svm->
//  }
//
//  for (size_t i = 0; i < returnVec.size(); i++)
//  {
//    //returnVec[i] = 1.0 / (1.0 + exp( returnVec[i])); // gives probs for class 1
//    returnVec[i] = 1.0 / (1.0 + exp(-returnVec[i])); // gives probs for class 2
//  }
//  
//  return returnVec;
//}

/**
\brief Guess Image Type

\param str String to guess
\return deduced type
*/