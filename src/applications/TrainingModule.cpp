/**
\file  TrainingSimulator.cpp

\brief The source file containing the SurvivaPredictor class, used to find Survival Prediction Index
Author: Saima Rathore
Library Dependecies: ITK 4.7+ <br>

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edutri

Copyright (c) 2016 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software/license.html

*/

#include "TrainingModule.h"
#include "opencv2/core/core.hpp"
#include "opencv2/ml.hpp"
#include <cmath>
#include "CaPTkUtils.h"

bool TrainingModule::CheckPerformanceStatus(double ist, double second, double third, double fourth, double fifth, double sixth, double seventh, double eighth, double ninth, double tenth)
{
  if (ist<tenth && second<tenth && third<tenth && fourth<tenth && fifth < tenth && sixth <tenth && seventh<tenth && eighth<tenth && ninth < tenth)
    return false;
  else
    return true;
}

VectorDouble TrainingModule::CalculatePerformanceMeasures(VectorDouble predictedLabels, VectorDouble GivenLabels)
{
  //calcualte performance measures
  double TP = 0;
  double TN = 0;
  double FP = 0;
  double FN = 0;
  VectorDouble result;

  for (int index = 0; index< predictedLabels.size(); index++)
  {
    if (predictedLabels[index] == 1 && GivenLabels[index] == 1)
      TP++;
    else if (predictedLabels[index] == -1 && GivenLabels[index] == -1)
      TN++;
    else if (predictedLabels[index] == 1 && GivenLabels[index] == -1)
      FP++;
    else if (predictedLabels[index] == -1 && GivenLabels[index] == 1)
      FN++;
    else
    {
    }
  }
  result.push_back((TP + TN) / predictedLabels.size());
  double sen = TP / (TP + FN);
  double spe = TN / (TN + FP);
  result.push_back(sen);
  result.push_back(spe);

  result.push_back((sen + spe) / 2);
  return result;
}

VectorDouble TrainingModule::CalculatePerformanceMeasures(VariableLengthVectorType predictedLabels, std::vector<double> GivenLabels)
{
  //calcualte performance measures
  double TP = 0;
  double TN = 0;
  double FP = 0;
  double FN = 0;
  VectorDouble result;

  for (unsigned int index = 0; index< predictedLabels.Size(); index++)
  {
    if (predictedLabels[index] == 1 && GivenLabels[index] == 1)
      TP++;
    else if (predictedLabels[index] == -1 && GivenLabels[index] == -1)
      TN++;
    else if (predictedLabels[index] == 1 && GivenLabels[index] == -1)
      FP++;
    else if (predictedLabels[index] == -1 && GivenLabels[index] == 1)
      FN++;
    else
    {
    }
  }
  result.push_back((TP + TN) / predictedLabels.Size());
  double sen = TP / (TP + FN);
  double spe = TN / (TN + FP);
  result.push_back(sen);
  result.push_back(spe);

  result.push_back((sen + spe) / 2);
  return result;
}

VectorDouble TrainingModule::InternalCrossValidation(VariableSizeMatrixType inputFeatures, std::vector<double> inputLabels, double cValue, double gValue, int kerneltype)
{
  VariableLengthVectorType predictedLabels;
  predictedLabels.SetSize(inputLabels.size());


  int fold_size = inputFeatures.Rows() / 5;
  //make loops for training and validation
  for (int index = 0; index < 5; index++)
  {
    VariableSizeMatrixType trainingdata;
    VariableSizeMatrixType testingdata;
    std::vector<int> trainingindices;
    std::vector<int> testingindices;
    std::vector<int> traininglabels;
    std::vector<int> testinglabels;

    for (int index2 = 0; index2 < fold_size; index2++)
    {
      int currentindex = index * fold_size + index2;
      testingindices.push_back(currentindex);
    }
    for (int index3 = 0; index3 < inputLabels.size(); index3++)
    {
      int found = 0;
      for (int index4 = 0; index4 < testingindices.size(); index4++)
      {
        if (index3 == testingindices[index4])
        {
          found = 1;
          break;
        }
      }
      if (found == 0)
        trainingindices.push_back(index3);
    }


    //copy training and test data from given dataset
    cv::Mat trainingData = cv::Mat::zeros(trainingindices.size(), inputFeatures.Cols(), CV_32FC1), trainingLabels = cv::Mat::zeros(trainingindices.size(), 1, CV_32FC1);
    cv::Mat testingData = cv::Mat::zeros(testingindices.size(), inputFeatures.Cols(), CV_32FC1), testingLabels = cv::Mat::zeros(testingindices.size(), 1, CV_32FC1);
    int trainingCounter = -1;
    int testingCounter = -1;
    for (int index = 0; index < trainingData.rows; index++)
    {
      trainingCounter++;
      trainingLabels.ptr< float >(trainingCounter)[0] = inputLabels[trainingindices[index]];
      for (int featureNo = 0; featureNo < trainingData.cols; featureNo++)
        trainingData.ptr< float >(trainingCounter)[featureNo] = inputFeatures(trainingindices[index], featureNo);
    }
    for (int index = 0; index < testingData.rows; index++)
    {
      testingCounter++;
      testingLabels.ptr< float >(testingCounter)[0] = inputLabels[testingindices[index]];
      for (int featureNo = 0; featureNo < testingData.cols; featureNo++)
        testingData.ptr< float >(testingCounter)[featureNo] = inputFeatures(testingindices[index], featureNo);
    }
    trainingLabels.convertTo(trainingLabels, CV_32SC1);
    //make an SVM model
    auto svm = cv::ml::SVM::create();
    svm->setType(cv::ml::SVM::C_SVC);
    svm->setC(cValue);

    if (kerneltype == 2)
    {
      svm->setGamma(gValue);
      svm->setKernel(cv::ml::SVM::RBF);
    }
    else
      svm->setKernel(cv::ml::SVM::LINEAR);
    bool res = true;
    std::string msg;

    try
    {
      res = svm->train(trainingData, cv::ml::ROW_SAMPLE, trainingLabels);
    }
    catch (cv::Exception ex)
    {
      msg = ex.what();
    }
    //apply SVM model on test data
    cv::Mat predicted(1, 1, CV_32F);
    for (int i = 0; i < testingData.rows; i++)
    {
      cv::Mat sample = testingData.row(i);
      svm->predict(sample, predicted, false);
      predictedLabels[testingindices[i]] = predicted.ptr< float >(0)[0];
    }
  }
  VectorDouble results = CalculatePerformanceMeasures(predictedLabels, inputLabels);
  return results;
}

VectorDouble TrainingModule::InternalCrossValidationResubstitution(VariableSizeMatrixType inputFeatures, std::vector<double> inputLabels, double cValue, double gValue, int kerneltype)
{
  VariableLengthVectorType predictedLabels;
  predictedLabels.SetSize(inputLabels.size());
  VariableLengthVectorType predictedDistances;
  predictedDistances.SetSize(inputLabels.size());


    //copy training and test data from given dataset
    cv::Mat trainingData = cv::Mat::zeros(inputFeatures.Rows(), inputFeatures.Cols(), CV_32FC1), trainingLabels = cv::Mat::zeros(inputFeatures.Rows(), 1, CV_32FC1);
    for (int index = 0; index < trainingData.rows; index++)
    {
      trainingLabels.ptr< float >(index)[0] = inputLabels[index];
      for (int featureNo = 0; featureNo < trainingData.cols; featureNo++)
        trainingData.ptr< float >(index)[featureNo] = inputFeatures(index, featureNo);
    }

    trainingLabels.convertTo(trainingLabels, CV_32SC1);
    //make an SVM model
    auto svm = cv::ml::SVM::create();
    svm->setType(cv::ml::SVM::C_SVC);
    //svm->setTermCriteria(cv::TermCriteria(cv::TermCriteria::MAX_ITER, 100, 1e-6));
    svm->setC(cValue);

    if (kerneltype == 2)
    {
      svm->setGamma(gValue);
      svm->setKernel(cv::ml::SVM::RBF);
    }
    else
      svm->setKernel(cv::ml::SVM::LINEAR);
    bool res = true;
    std::string msg;

    try
    {
      res = svm->train(trainingData, cv::ml::ROW_SAMPLE, trainingLabels);
    }
    catch (cv::Exception ex)
    {
      msg = ex.what();
    }
    //apply SVM model on test data
    cv::Mat predicted(1, 1, CV_32F);
    for (int i = 0; i < trainingData.rows; i++)
    {
      cv::Mat sample = trainingData.row(i);
      svm->predict(sample, predicted, false);
      predictedLabels[i] = predicted.ptr< float >(0)[0];
    }
    for (int i = 0; i < trainingData.rows; i++)
    {
      cv::Mat sample = trainingData.row(i);
      svm->predict(sample, predicted, true);
      predictedDistances[i] = predicted.ptr< float >(0)[0];
    }

  VectorDouble results = CalculatePerformanceMeasures(predictedLabels, inputLabels);

  return results;
}

VectorDouble TrainingModule::InternalCrossValidationSplitTrainTest(VariableSizeMatrixType inputFeatures, std::vector<double> inputLabels, double cValue, double gValue, int kerneltype, int counter, std::string outputfolder)
{
  VariableLengthVectorType predictedLabels;
  predictedLabels.SetSize(inputLabels.size());

  //copy training and test data from given dataset
  cv::Mat trainingData = cv::Mat::zeros(inputFeatures.Rows(), inputFeatures.Cols(), CV_32FC1), trainingLabels = cv::Mat::zeros(inputLabels.size(), 1, CV_32FC1);
  cv::Mat testingData = cv::Mat::zeros(inputFeatures.Rows(), inputFeatures.Cols(), CV_32FC1), testingLabels = cv::Mat::zeros(inputLabels.size(), 1, CV_32FC1);

  for (int index = 0; index < trainingData.rows; index++)
  {
    trainingLabels.ptr< float >(index)[0] = inputLabels[index];
    testingLabels.ptr< float >(index)[0] = inputLabels[index];
    for (int featureNo = 0; featureNo < trainingData.cols; featureNo++)
    {
      double x = inputFeatures(index, featureNo);
      x = floor(x * 100 + 0.5) / 100;
      trainingData.ptr< float >(index)[featureNo] = x;
      testingData.ptr< float >(index)[featureNo] = x;
    }
  }

  std::ofstream myfile;
  myfile.open(outputfolder + "/" + std::to_string(counter) + ".csv");
  for (unsigned int index1 = 0; index1 < inputFeatures.Rows(); index1++)
  {
    for (unsigned int index2 = 0; index2 < inputFeatures.Cols(); index2++)
    {
      if (index2 == 0)
        myfile << std::to_string(inputFeatures[index1][index2]);
      else
        myfile << "," << std::to_string(inputFeatures[index1][index2]);
    }
    myfile << "\n";
  }
  myfile.close();

  myfile.open(outputfolder + "/" + std::to_string(counter) + "_labels.csv");
  for (unsigned int index1 = 0; index1 < inputFeatures.Rows(); index1++)
  {
    if (index1 == 0)
      myfile << std::to_string(inputLabels[index1]);
    else
      myfile << "," << std::to_string(inputLabels[index1]);
  }
  myfile << "\n";
  myfile.close();

  trainingLabels.convertTo(trainingLabels, CV_32SC1);
  //make an SVM model
  auto svm = cv::ml::SVM::create();
  svm->setType(cv::ml::SVM::C_SVC);
  svm->setC(cValue);

  if (kerneltype == 2)
  {
    svm->setGamma(gValue);
    svm->setKernel(cv::ml::SVM::RBF);
  }
  else
    svm->setKernel(cv::ml::SVM::LINEAR);
  bool res = true;
  std::string msg;

  try
  {
    res = svm->train(trainingData, cv::ml::ROW_SAMPLE, trainingLabels);
  }
  catch (cv::Exception ex)
  {
    msg = ex.what();
  }
  //apply SVM model on test data
  cv::Mat predicted(1, 1, CV_32F);
  for (int i = 0; i < testingData.rows; i++)
  {
    cv::Mat sample = testingData.row(i);
    svm->predict(sample, predicted, false);
    predictedLabels[i] = predicted.ptr< float >(0)[0];
  }

  myfile.open(outputfolder + "/" + std::to_string(counter) + "_predictions.csv");
  for (unsigned int index1 = 0; index1 < predictedLabels.Size(); index1++)
  {
    if (index1 == 0)
      myfile << std::to_string(predictedLabels[index1]);
    else
      myfile << "," << std::to_string(predictedLabels[index1]);
  }
  myfile << "\n";
  myfile.close();

  VectorDouble results = CalculatePerformanceMeasures(predictedLabels, inputLabels);
  return results;
}

VectorDouble TrainingModule::CrossValidation(const VariableSizeMatrixType inputFeatures, const VariableLengthVectorType inputLabels, const std::string outputfolder, const int classifiertype, const int number_of_folds, const int featureselectiontype,const int optimizationType, const int crossvalidationType)
{
  MapType FoldingDataMap;
  int fold_size = inputLabels.Size() / number_of_folds;
  int remainder = inputLabels.Size() % number_of_folds;

  //make loops for stratifing training and validation sets, CVO partition like structure
  for (int index = 0; index < number_of_folds; index++)
  {
    //find training and testing indices
    std::vector<double> trainingindices;
    std::vector<double> traininglabels;
    VariableSizeMatrixType trainingfeatures;

    std::vector<double> testingindices;
    std::vector<double> testinglabels;
    std::vector<double> predictedlabels;
    std::vector<double> predicteddistances;
    VariableSizeMatrixType testingfeatures;


    for (int index2 = 0; index2 < fold_size; index2++)
    {
      int currentindex = index * fold_size + index2;
      testingindices.push_back(currentindex);
    }
    //copy the remaining subjects in the last fold
    if (index == number_of_folds - 1)
    {
      for (int remainingLoop = 0; remainingLoop < remainder; remainingLoop++)
        testingindices.push_back(inputLabels.Size() - 1 - remainingLoop);
    }

    for (unsigned int index3 = 0; index3 < inputLabels.Size(); index3++)
    {
      int found = 0;
      for (unsigned int index4 = 0; index4 < testingindices.size(); index4++)
      {
        if (index3 == testingindices[index4])
        {
          found = 1;
          break;
        }
      }
      if (found == 0)
        trainingindices.push_back(index3);
    }
    //find training and testing labels and features corresponding to training and testing indices
    trainingfeatures.SetSize(trainingindices.size(), inputFeatures.Cols());
    for (unsigned int i = 0; i < trainingindices.size(); i++)
    {
      traininglabels.push_back(inputLabels[trainingindices[i]]);
      for (unsigned int j = 0; j < trainingfeatures.Cols(); j++)
        trainingfeatures(i, j) = inputFeatures(trainingindices[i], j);
    }
    testingfeatures.SetSize(testingindices.size(), inputFeatures.Cols());
    for (unsigned int i = 0; i < testingindices.size(); i++)
    {
      testinglabels.push_back(inputLabels[testingindices[i]]);
      predictedlabels.push_back(-1);
      predicteddistances.push_back(0);
      for (unsigned int j = 0; j < testingfeatures.Cols(); j++)
        testingfeatures(i, j) = inputFeatures(testingindices[i], j);
    }
    std::tuple<VectorDouble, VectorDouble, VariableSizeMatrixType, VectorDouble, VectorDouble, VariableSizeMatrixType, VectorDouble,VectorDouble> new_tuple(trainingindices, traininglabels, trainingfeatures, testingindices, testinglabels, testingfeatures, predictedlabels,predicteddistances);
    FoldingDataMap[index] = new_tuple;
  }
  for (int index = 0; index < number_of_folds; index++)
  {
     std::cout << "***************************** index=" << index << "*********************************" << std::endl << std::endl;
    cbica::Logging(loggerFile, "***************************** index=" + std::to_string(index) + "*********************************\n\n");
    VariableSizeMatrixType scaledFeatureSet = std::get<2>(FoldingDataMap[index]);
    std::vector<double> traininglabels = std::get<1>(FoldingDataMap[index]);


    std::vector<int> FinalSelectedFeatures;
    VectorDouble crossvalidatedAccuracies;
    if (featureselectiontype == 1)
      FinalSelectedFeatures = EffectSizeBasedFeatureSelection(scaledFeatureSet, traininglabels, classifiertype, optimizationType,crossvalidationType, crossvalidatedAccuracies);
    //else if (featureselectiontype == 2)
    //  FinalSelectedFeatures = CorrelationBasedFeatureSelection(scaledFeatureSet,traininglabels,classifiertype);
    else if (featureselectiontype == 3)
      FinalSelectedFeatures = SVMFFSBasedFeatureSelection(scaledFeatureSet, traininglabels, classifiertype, optimizationType, crossvalidationType, crossvalidatedAccuracies);
    //else if (featureselectiontype == 4)
    //  FinalSelectedFeatures = SVMRFEBasedFeatureSelection(scaledFeatureSet,traininglabels,classifiertype);
    //else
    //{
    //  std::cout << "The selected feature selection method is not avaialble inside CaPTk. Please select another method..." << std::endl;
    //}


    //Write the selected features
     WriteCSVFiles(FinalSelectedFeatures, outputfolder + "/SelectedFeatures_FoldNo_" + std::to_string(index + 1) + ".csv");
     WriteCSVFiles(crossvalidatedAccuracies, outputfolder + "/CrossValidatedAccuracies_FoldNo_" + std::to_string(index + 1) + ".csv");

    //optimize classifiers parameters on the selected feature set
    VariableSizeMatrixType FinalSelectedFeatureSet;
    FinalSelectedFeatureSet.SetSize(scaledFeatureSet.Rows(), FinalSelectedFeatures.size());

    //copy the already selected features to the current feature set
    for (unsigned int j = 0; j < FinalSelectedFeatureSet.Rows(); j++)
      for (unsigned int k = 0; k < FinalSelectedFeatures.size(); k++)
        FinalSelectedFeatureSet(j, k) = scaledFeatureSet(j, FinalSelectedFeatures[k]);

    double bestCV = 0;
    double bestC = 1;
    double bestG = 1 / FinalSelectedFeatures.size();
    if (classifiertype == 2)
    {
      for (double cValue = -5; cValue <= 5; cValue = cValue + 1)
        for (double gValue = -5; gValue <= 5; gValue = gValue + 1)
        {
          //this cross-validation mechanism trains and tests on the same dataset. 
          // we have another function avialable called 'InternalCrossValidation' that does cross-validation via 5-fold cross-validation. 
          //TODO: Consider adding an option to ask the user what type of internal cross-validation they want. 

          VectorDouble result = InternalCrossValidationResubstitution(FinalSelectedFeatureSet, traininglabels, pow(2, cValue), pow(2, gValue), classifiertype);
          if (result[3] > bestCV)
          {
            bestC = pow(2, cValue);
            bestG = pow(2, gValue);
            bestCV = result[3];
          }
        }
    }
    else
    {
      for (double cValue = -5; cValue <= 5; cValue = cValue + 1)
      {
        VectorDouble result = InternalCrossValidationResubstitution(FinalSelectedFeatureSet, traininglabels, pow(2, cValue), 0.1, classifiertype);
        if (result[3] > bestCV)
        {
          bestC = pow(2, cValue);
          bestCV = result[3];
        }
      }
    }
    std::cout << "Optimal C=" << bestC << std::endl;
    std::cout << "Optimal gamma=" << bestG << std::endl;






    //train a new model on selcted features and optimal values of classifiers' parameters
    //copy training and test data from given dataset
    cv::Mat trainingData = cv::Mat::zeros(FinalSelectedFeatureSet.Rows(), FinalSelectedFeatureSet.Cols(), CV_32FC1);
    cv::Mat trainingLabels = cv::Mat::zeros(FinalSelectedFeatureSet.Rows(), 1, CV_32FC1);
    cv::Mat testingData = cv::Mat::zeros(std::get<5>(FoldingDataMap[index]).Rows(), FinalSelectedFeatureSet.Cols(), CV_32FC1);
    cv::Mat testingLabels = cv::Mat::zeros(std::get<4>(FoldingDataMap[index]).size(), 1, CV_32FC1);

    int trainingCounter = 0;
    int testingCounter = 0;
    for (int copyDataCounter = 0; copyDataCounter < trainingData.rows; copyDataCounter++)
    {
      trainingLabels.ptr< float >(copyDataCounter)[0] = std::get<1>(FoldingDataMap[index])[copyDataCounter];
      for (int copyDataCounter2 = 0; copyDataCounter2 <trainingData.cols; copyDataCounter2++)
        trainingData.ptr< float >(copyDataCounter)[copyDataCounter2] = FinalSelectedFeatureSet(copyDataCounter, copyDataCounter2);
    }
    for (int copyDataCounter = 0; copyDataCounter < testingData.rows; copyDataCounter++)
    {
      testingLabels.ptr< float >(copyDataCounter)[0] = std::get<4>(FoldingDataMap[index])[copyDataCounter];
      for (int copyDataCounter2 = 0; copyDataCounter2 <testingData.cols; copyDataCounter2++)
        testingData.ptr< float >(copyDataCounter)[copyDataCounter2] = std::get<5>(FoldingDataMap[index])(copyDataCounter, FinalSelectedFeatures[copyDataCounter2]);
    }
    if (classifiertype == 2)
      // std::cout << "index=" << index << ", Best C = " << bestC << ", Best G=" << bestG << ", Best CV=" << bestCV << ", Features=" << numberOfSelectedFeatures << std::endl;
      cbica::Logging(loggerFile, "index=" + std::to_string(index) + ", Best C = " + std::to_string(bestC) + ", Best G = " + std::to_string(bestG) + ", Best CV = " + std::to_string(bestCV) + ", Features=" + std::to_string(FinalSelectedFeatures.size()) + "\n");
    else
      // std::cout << "index=" << index << ", Best C = " << bestC << ", Best CV=" << bestCV << ", Features=" << numberOfSelectedFeatures << std::endl;
      cbica::Logging(loggerFile, "index=" + std::to_string(index) + ", Best C = " + std::to_string(bestC) + ", Features=" + std::to_string(FinalSelectedFeatures.size()) + "\n");

    trainingLabels.convertTo(trainingLabels, CV_32SC1);
    auto svm = cv::ml::SVM::create();
    svm->setType(cv::ml::SVM::C_SVC);
    svm->setC(bestC);

    if (classifiertype == 2)
    {
      svm->setGamma(bestG);
      svm->setKernel(cv::ml::SVM::RBF);
    }
    else
      svm->setKernel(cv::ml::SVM::LINEAR);

    bool res = true;
    std::string msg;
    VectorDouble predictedLabels;
    VectorDouble predictedDistances;
    try
    {
      res = svm->train(trainingData, cv::ml::ROW_SAMPLE, trainingLabels);
      svm->save(outputfolder + "/Model_FoldNo_" + std::to_string(index + 1) + ".xml");
    }
    catch (cv::Exception ex)
    {
      msg = ex.what();
    }
    //apply SVM model on test data to get labels
    cv::Mat predicted(1, 1, CV_32F);
    for (int i = 0; i < testingData.rows; i++)
    {
      cv::Mat sample = testingData.row(i);
      svm->predict(sample, predicted, false);
      predictedLabels.push_back(predicted.ptr< float >(0)[0]);
    }
    //apply SVM model on test data to get score
    for (int i = 0; i < testingData.rows; i++)
    {
      cv::Mat sample = testingData.row(i);
      svm->predict(sample, predicted, true);
      predictedDistances.push_back(predicted.ptr< float >(0)[0]);
    }
    for (int i = 0; i < predictedDistances.size(); i++)
      predictedDistances[i] = std::abs(predictedDistances[i])*predictedLabels[i];
    
    std::tuple<VectorDouble, VectorDouble, VariableSizeMatrixType, VectorDouble, VectorDouble, VariableSizeMatrixType, VectorDouble,VectorDouble> new_tuple_duplicate = FoldingDataMap[index];
    std::get<6>(new_tuple_duplicate) = predictedLabels;
    std::get<7>(new_tuple_duplicate) = predictedDistances;
    FoldingDataMap[index] = new_tuple_duplicate;
  }

  //combining the predicted results from all the map entries
  VectorDouble FinalTargetLabels;
  VectorDouble FinalPredictedLabels;
  VectorDouble FinalPredictedDistances;
  for (auto const &mapiterator : FoldingDataMap)
  {
    VectorDouble TargetLabels = std::get<4>(mapiterator.second);
    VectorDouble PredictedLabels = std::get<6>(mapiterator.second);
    VectorDouble PredictedDistances = std::get<7>(mapiterator.second);
    for (int index = 0; index < TargetLabels.size(); index++)
    {
      FinalTargetLabels.push_back(TargetLabels[index]);
      FinalPredictedLabels.push_back(PredictedLabels[index]);
      FinalPredictedDistances.push_back(PredictedDistances[index]);
    }
  }

  //calcualte final performance and write in the csv file
  VectorDouble FinalPerformance = CalculatePerformanceMeasures(FinalPredictedLabels, FinalTargetLabels);
  WriteCSVFiles(FinalPredictedLabels, outputfolder + "/predicted_labels.csv");
  WriteCSVFiles(FinalPredictedDistances, outputfolder + "/predicted_distances.csv");

  std::ofstream myfile;
  myfile.open(outputfolder + "/performance.csv");
  myfile << "Accuracy,Sensitivity,Specificity,BalancedAccuracy \n";
  myfile << std::to_string(FinalPerformance[0]) + "," + std::to_string(FinalPerformance[1]) + "," + std::to_string(FinalPerformance[2]) + "," + std::to_string(FinalPerformance[3]) + "\n";
  myfile.close();


  return FinalPerformance;
}



VectorDouble TrainingModule::EffectSizeFeatureSelection(const VariableSizeMatrixType training_features, std::vector<double> target)
{
  //make set 1and set2
  int NoOfSamplesC1 = 0;
  int NoOfSamplesC2 = 0;
  std::vector<double> indices_set1;
  std::vector<double> indices_set2;
  VariableSizeMatrixType features_set1;
  VariableSizeMatrixType features_set2;
  VariableLengthVectorType mean_set1;
  VariableLengthVectorType mean_set2;

  for (int index = 0; index < target.size(); index++)
  {
    if (target[index] == -1)
      NoOfSamplesC1++;
    else if (target[index] == 1)
      NoOfSamplesC2++;
  }
  features_set1.SetSize(NoOfSamplesC1, training_features.Cols());
  features_set2.SetSize(NoOfSamplesC2, training_features.Cols());
  mean_set1.SetSize(training_features.Cols());
  mean_set2.SetSize(training_features.Cols());

  NoOfSamplesC1 = 0;
  NoOfSamplesC2 = 0;
  for (unsigned int index = 0; index < target.size(); index++)
  {
    if (target[index] == -1)
    {
      for (unsigned int featureNo = 0; featureNo < training_features.Cols(); featureNo++)
        features_set1(NoOfSamplesC1, featureNo) = training_features(index, featureNo);
      NoOfSamplesC1++;
    }
    else if (target[index] == 1)
    {
      for (unsigned int featureNo = 0; featureNo < training_features.Cols(); featureNo++)
        features_set2(NoOfSamplesC2, featureNo) = training_features(index, featureNo);
      NoOfSamplesC2++;
    }
  }
  std::vector<double> EffectSize;
  for (unsigned int featureNo = 0; featureNo < training_features.Cols(); featureNo++)
  {
    double temp = 0.0;
    for (int sampleNo = 0; sampleNo < NoOfSamplesC1; sampleNo++)
      temp = temp + features_set1(sampleNo, featureNo);
    mean_set1[featureNo] = temp / NoOfSamplesC1;

    temp = 0.0;
    for (int sampleNo = 0; sampleNo < NoOfSamplesC2; sampleNo++)
      temp = temp + features_set2(sampleNo, featureNo);
    mean_set2[featureNo] = temp / NoOfSamplesC2;


    double sum1 = 0;
    double sum2 = 0;
    for (int sampleNo = 0; sampleNo < NoOfSamplesC1; sampleNo++)
      sum1 = sum1 + (features_set1(sampleNo, featureNo) - mean_set1[featureNo])*(features_set1(sampleNo, featureNo) - mean_set1[featureNo]);

    for (int sampleNo = 0; sampleNo < NoOfSamplesC2; sampleNo++)
      sum2 = sum2 + (features_set2(sampleNo, featureNo) - mean_set2[featureNo])*(features_set2(sampleNo, featureNo) - mean_set2[featureNo]);

    double SC1 = sum1 / (NoOfSamplesC1 - 1);
    double SC2 = sum2 / (NoOfSamplesC2 - 1);
    double SP = ((NoOfSamplesC1 - 1)*SC1 + (NoOfSamplesC2 - 1)*SC2) / (NoOfSamplesC1 + NoOfSamplesC2 - 2);
    double currentvalue = (mean_set1[featureNo] - mean_set2[featureNo]) / sqrt(SP);
    if(std::isnan(currentvalue))
      EffectSize.push_back(0.0001);
    else
      EffectSize.push_back(currentvalue);
  }
  //std::vector<size_t> indices = sort_indexes(EffectSize);
  //VariableSizeMatrixType selected_feature_set;

  //for (int index1 = 0; index1 < training_features.Rows(); index1++)
  //	for (int index = 0; index < no_of_features; index++)
  //		selected_feature_set(index1, index) = training_features(index1, indices[index]);

  return EffectSize;
}
//
//VectorDouble TrainingModule::CorrelationBasedFeatureSelection(const VariableSizeMatrixType training_features, std::vector<double> target)
//{
//  //make set 1and set2
//  int NoOfSamplesC1 = 0;
//  int NoOfSamplesC2 = 0;
//  std::vector<double> indices_set1;
//  std::vector<double> indices_set2;
//  VariableSizeMatrixType features_set1;
//  VariableSizeMatrixType features_set2;
//  VariableLengthVectorType mean_set1;
//  VariableLengthVectorType mean_set2;
//
//  for (int index = 0; index < target.size(); index++)
//  {
//    if (target[index] == -1)
//      NoOfSamplesC1++;
//    else if (target[index] == 1)
//      NoOfSamplesC2++;
//  }
//  features_set1.SetSize(NoOfSamplesC1, training_features.Cols());
//  features_set2.SetSize(NoOfSamplesC2, training_features.Cols());
//  mean_set1.SetSize(training_features.Cols());
//  mean_set2.SetSize(training_features.Cols());
//
//  NoOfSamplesC1 = 0;
//  NoOfSamplesC2 = 0;
//  for (unsigned int index = 0; index < target.size(); index++)
//  {
//    if (target[index] == -1)
//    {
//      for (unsigned int featureNo = 0; featureNo < training_features.Cols(); featureNo++)
//        features_set1(NoOfSamplesC1, featureNo) = training_features(index, featureNo);
//      NoOfSamplesC1++;
//    }
//    else if (target[index] == 1)
//    {
//      for (unsigned int featureNo = 0; featureNo < training_features.Cols(); featureNo++)
//        features_set2(NoOfSamplesC2, featureNo) = training_features(index, featureNo);
//      NoOfSamplesC2++;
//    }
//  }
//  std::vector<double> EffectSize;
//  for (unsigned int featureNo = 0; featureNo < training_features.Cols(); featureNo++)
//  {
//    double temp = 0.0;
//    for (int sampleNo = 0; sampleNo < NoOfSamplesC1; sampleNo++)
//      temp = temp + features_set1(sampleNo, featureNo);
//    mean_set1[featureNo] = temp / NoOfSamplesC1;
//
//    temp = 0.0;
//    for (int sampleNo = 0; sampleNo < NoOfSamplesC2; sampleNo++)
//      temp = temp + features_set2(sampleNo, featureNo);
//    mean_set2[featureNo] = temp / NoOfSamplesC2;
//
//
//    double sum1 = 0;
//    double sum2 = 0;
//    for (int sampleNo = 0; sampleNo < NoOfSamplesC1; sampleNo++)
//      sum1 = sum1 + (features_set1(sampleNo, featureNo) - mean_set1[featureNo])*(features_set1(sampleNo, featureNo) - mean_set1[featureNo]);
//
//    for (int sampleNo = 0; sampleNo < NoOfSamplesC2; sampleNo++)
//      sum2 = sum2 + (features_set2(sampleNo, featureNo) - mean_set2[featureNo])*(features_set2(sampleNo, featureNo) - mean_set2[featureNo]);
//
//    double SC1 = sum1 / (NoOfSamplesC1 - 1);
//    double SC2 = sum2 / (NoOfSamplesC2 - 1);
//    double SP = ((NoOfSamplesC1 - 1)*SC1 + (NoOfSamplesC2 - 1)*SC2) / (NoOfSamplesC1 + NoOfSamplesC2 - 2);
//    double currentvalue = (mean_set1[featureNo] - mean_set2[featureNo]) / sqrt(SP);
//    if (std::isnan(currentvalue))
//      EffectSize.push_back(0.0001);
//    else
//      EffectSize.push_back(currentvalue);
//  }
//  //std::vector<size_t> indices = sort_indexes(EffectSize);
//  //VariableSizeMatrixType selected_feature_set;
//
//  //for (int index1 = 0; index1 < training_features.Rows(); index1++)
//  //	for (int index = 0; index < no_of_features; index++)
//  //		selected_feature_set(index1, index) = training_features(index1, indices[index]);
//
//  return EffectSize;
//}



//----------------------------------------------------------------------------------------------
VectorDouble TrainingModule::CombineEstimates(const VariableLengthVectorType &estimates1, const VariableLengthVectorType &estimates2)
{
  VectorDouble returnVec;
  returnVec.resize(estimates1.Size());
  for (size_t i = 0; i < estimates1.Size(); i++)
  {
    float temp_abs, temp_pos1, temp_neg1, temp_1, temp_2;
    // estimate for 1st vector
    if (std::abs(estimates1[i]) < 2)
    {
      temp_abs = estimates1[i];
    }
    else
    {
      temp_abs = 0;
    }

    if (estimates1[i] > 1)
    {
      temp_pos1 = 1;
    }
    else
    {
      temp_pos1 = 0;
    }

    if (estimates1[i] < -1)
    {
      temp_neg1 = 1;
    }
    else
    {
      temp_neg1 = 0;
    }
    temp_1 = temp_abs + (temp_pos1 - temp_neg1);

    // estimate for 2nd vector, all temp values are getting overwritten
    if (std::abs(estimates2[i]) < 2)
    {
      temp_abs = estimates2[i];
    }
    else
    {
      temp_abs = 0;
    }

    if (estimates2[i] > 1)
    {
      temp_pos1 = 1;
    }
    else
    {
      temp_pos1 = 0;
    }

    if (estimates2[i] < -1)
    {
      temp_neg1 = 1;
    }
    else
    {
      temp_neg1 = 0;
    }
    temp_2 = temp_abs + (temp_pos1 - temp_neg1);

    // combine the two
    returnVec[i] = temp_1 + temp_2;
  }

  return returnVec;
}

VectorDouble TrainingModule::testOpenCVSVM(const VariableSizeMatrixType &testingData, const std::string &inputModelName)
{
  //auto svm = cv::Algorithm::load<cv::ml::SVM>(inputModelName);

  ////std::ofstream file;
  ////file.open("Z:/Projects/testingData.csv");
  ////for (size_t i = 0; i < testingData.Rows(); i++)
  ////{
  ////  for (size_t j = 0; j < testingData.Cols(); j++)
  ////  {
  ////    file << testingData(i, j) << ",";
  ////  }
  ////  file << "\n";
  ////}
  ////file.close();

  VectorDouble returnVec;
  ////VariableSizeMatrixType returnMat;
  ////returnMat.SetSize(testingData.Rows(), 1);
  ////returnVec.resize(testingData.Rows());
  //cv::Mat testingDataMat = cv::Mat::zeros(testingData.Rows(), testingData.Cols() - 1, CV_32FC1), outputProbs;

  //// fast cv::Mat access
  //for (unsigned int i = 0; i < testingData.Rows(); ++i)
  //{
  //  for (unsigned int j = 0; j < testingData.Cols(); ++j)
  //  {
  //    testingDataMat.ptr< float >(i)[j] = testingData(i, j);
  //  }
  //}

  //// slow cv::Mat access
  ////for (size_t i = 0; i < testingData.Rows(); i++)
  ////{
  ////  for (size_t j = 0; j < testingData.Cols() - 1; j++)
  ////  {
  ////    testingDataMat.at< float >(i, j) = testingData(i, j);
  ////  }
  ////}

  //// see http://docs.opencv.org/trunk/db/d7d/classcv_1_1ml_1_1StatModel.html#af1ea864e1c19796e6264ebb3950c0b9a for details regarding why '1'
  ////svm->predict(testingDataMat, returnVec, 1);

  //returnVec.resize(testingDataMat.rows);
  //cv::Mat predicted(1, 1, CV_32F);
  //for (int i = 0; i < testingDataMat.rows; i++)
  //{
  //  cv::Mat sample = testingDataMat.row(i);
  //  /*float value = */svm->predict(sample, predicted, true/*cv::ml::StatModel::RAW_OUTPUT*/);
  //  //float p = /*predicted.at<float>(0, 0)*/predicted.ptr< float >(0)[0];
  //  returnVec[i] = predicted.ptr< float >(0)[0];
  //}

  ////for (size_t i = 0; i < outputProbs.rows; i++)
  ////{
  ////  returnMat[i] = outputProbs.at<float>(i, 0);
  ////}

  return returnVec;
}




bool TrainingModule::Run(const std::string inputFeaturesFile,
  const std::string outputdirectory,
  const std::string inputLabelsFile,
  const std::string modeldirectory,
  const int classifiertype, const int foldtype,
  const int confType, const int featureselectiontype,
  const int optimizationType, const int crossvalidationType)
{
  std::cout << "Training module." << std::endl;
  //reading features and labels from the input data
  VariableSizeMatrixType FeaturesOfAllSubjects;
  VariableLengthVectorType LabelsOfAllSubjects;
  typedef vnl_matrix<double> MatrixType;
  MatrixType dataMatrix;

  try
  {
    CSVFileReaderType::Pointer readerMean = CSVFileReaderType::New();
    readerMean->SetFileName(inputFeaturesFile);
    readerMean->SetFieldDelimiterCharacter(',');
    readerMean->HasColumnHeadersOff();
    readerMean->HasRowHeadersOff();
    readerMean->Parse();
    dataMatrix = readerMean->GetArray2DDataObject()->GetMatrix();
    FeaturesOfAllSubjects.SetSize(dataMatrix.rows() - 1, dataMatrix.columns() - 1);

    for (unsigned int i = 1; i < dataMatrix.rows(); i++)
      for (unsigned int j = 1; j < dataMatrix.cols(); j++)
        FeaturesOfAllSubjects(i - 1, j - 1) = dataMatrix(i, j);
  }
  catch (const std::exception& e1)
  {
    std::cerr << "Cannot find the feature file in the input directory. Error code : " + std::string(e1.what()) << "\n";
    return false;
  }

  //This segment of the code removes static features from the feature set. However, due to the removal of static features, the size of the feature set becomes smaller and introduces discreency between the training and the test dataset (that user may use later) in terms of the size and the selected features. 
  //therefore, this code is commented until we find a better solution

  //std::cout << "Remove static features from the feature set." << std::endl;
  //unsigned int NumberOfSamples = FeaturesOfAllSubjects.Rows();
  //unsigned int NumberOfFeatures = FeaturesOfAllSubjects.Cols();
  //std::vector<int> mStdVector;
  //for (unsigned int featureNo = 0; featureNo < NumberOfFeatures; featureNo++)
  //{
  //  double temp = 0.0;
  //  for (unsigned int sampleNo = 0; sampleNo < NumberOfSamples; sampleNo++)
  //    temp += FeaturesOfAllSubjects(sampleNo, featureNo);
  //  double mean = temp / NumberOfSamples;
  //  temp = 0.0;
  //  for (unsigned int sampleNo = 0; sampleNo < NumberOfSamples; sampleNo++)
  //    temp += (FeaturesOfAllSubjects(sampleNo, featureNo) - mean)*(FeaturesOfAllSubjects(sampleNo, featureNo) - mean);
  //  if (std::sqrt(temp / (NumberOfSamples - 1)) == 0)
  //    mStdVector.push_back(featureNo);
  //}

  //FeaturesOfAllSubjectsAfterRemovingStaticFeatures.SetSize(FeaturesOfAllSubjects.Rows(), FeaturesOfAllSubjects.Cols() - mStdVector.size());

  //int featureCounter = 0;
  //for (unsigned int featureNo = 0; featureNo < NumberOfFeatures; featureNo++)
  //{
  //  //check for the presence of this feature vector in the list of features to exclude
  //  std::vector<int>::iterator it;
  //  it = std::find(mStdVector.begin(), mStdVector.end(), featureNo);
  //  if (it != mStdVector.end())
  //    continue;

  //  for (unsigned int sampleNo = 0; sampleNo < NumberOfSamples; sampleNo++)
  //    FeaturesOfAllSubjectsAfterRemovingStaticFeatures(sampleNo, featureCounter) = FeaturesOfAllSubjects(sampleNo, featureNo);
  //  featureCounter++;
  //}
  if (confType != CAPTK::ClassificationConfigurationType::CONF_TYPE_SPLIT_TEST)
  {
    try
    {
      CSVFileReaderType::Pointer readerMean = CSVFileReaderType::New();
      readerMean->SetFileName(inputLabelsFile);
      readerMean->SetFieldDelimiterCharacter(',');
      readerMean->HasColumnHeadersOff();
      readerMean->HasRowHeadersOff();
      readerMean->Parse();
      dataMatrix = readerMean->GetArray2DDataObject()->GetMatrix();
      LabelsOfAllSubjects.SetSize(dataMatrix.rows() - 1);

      for (unsigned int i = 1; i < dataMatrix.rows(); i++)
        LabelsOfAllSubjects[i - 1] = dataMatrix(i, 0);
    }
    catch (const std::exception& e1)
    {
      std::cerr << "Cannot find the labels file in the input directory. Error code : " + std::string(e1.what()) << "\n";
      return false;
    }
  }

  std::cout << "Data loaded." << std::endl;

  TrainingModule mTrainingSimulator;
  VectorDouble FinalResult;

  if (confType == CAPTK::ClassificationConfigurationType::CONF_TYPE_KFOLD_CV)
  {  
    //z-scoring of input features and saving corresponding mean and standard deviation in the output directory
    FeatureScalingClass mFeaturesScaling;
    VariableSizeMatrixType scaledFeatureSet;
    VariableLengthVectorType meanVector;
    VariableLengthVectorType stdVector;
    mFeaturesScaling.ScaleGivenTrainingFeatures(FeaturesOfAllSubjects, scaledFeatureSet, meanVector, stdVector);

    //remove the nan values
    for (unsigned int index1 = 0; index1 < scaledFeatureSet.Rows(); index1++)
    {
      for (unsigned int index2 = 0; index2 < scaledFeatureSet.Cols(); index2++)
      {
        if (std::isnan(scaledFeatureSet[index1][index2]))
          scaledFeatureSet[index1][index2] = 0;
      }
    }
    for (unsigned int index1 = 0; index1 < meanVector.Size(); index1++)
    {
      if (std::isnan(meanVector[index1]))
        meanVector[index1] = 0;
      if (std::isnan(stdVector[index1]))
        stdVector[index1] = 0;
    }
    WriteCSVFiles(scaledFeatureSet, outputdirectory + "/scaled_feature_set.csv");
    WriteCSVFiles(meanVector,outputdirectory + "/zscore_mean.csv");
    WriteCSVFiles(stdVector, outputdirectory + "/zscore_std.csv");

    std::cout << "Scaling parameters written." << std::endl;
    FinalResult = mTrainingSimulator.CrossValidation(scaledFeatureSet, LabelsOfAllSubjects, outputdirectory, classifiertype, foldtype,featureselectiontype,optimizationType,crossvalidationType);
  }
  else if (confType == CAPTK::ClassificationConfigurationType::CONF_TYPE_DOUBLE)
  {
    //scaling of input features and saving corresponding mean and standard deviation in the output directory
    FeatureScalingClass mFeaturesScaling;
    VariableSizeMatrixType scaledFeatureSet;
    VariableLengthVectorType meanVector;
    VariableLengthVectorType stdVector;
    mFeaturesScaling.ScaleGivenTrainingFeatures(FeaturesOfAllSubjects, scaledFeatureSet, meanVector, stdVector);


    //remove the nan values
    for (unsigned int index1 = 0; index1 < scaledFeatureSet.Rows(); index1++)
    {
      for (unsigned int index2 = 0; index2 < scaledFeatureSet.Cols(); index2++)
      {
        if (std::isnan(scaledFeatureSet[index1][index2]))
          scaledFeatureSet[index1][index2] = 0;
      }
    }
    for (unsigned int index1 = 0; index1 < meanVector.Size(); index1++)
    {
      if (std::isnan(meanVector[index1]))
        meanVector[index1] = 0;
      if (std::isnan(stdVector[index1]))
        stdVector[index1] = 0;
    }



    std::ofstream myfile;
    myfile.open(outputdirectory + "/scaled_feature_set.csv");
    for (unsigned int index1 = 0; index1 < scaledFeatureSet.Rows(); index1++)
    {
      std::string onerow = "";
      for (unsigned int index2 = 0; index2 < scaledFeatureSet.Cols(); index2++)
      {
        if (index2 == 0)
          onerow = std::to_string(scaledFeatureSet[index1][index2]);
        else
          onerow = onerow + "," + std::to_string(scaledFeatureSet[index1][index2]);
      }
      onerow = onerow + "\n";
      myfile << onerow;
    }
    myfile.close();

    myfile.open(outputdirectory + "/zscore_mean.csv");
    for (unsigned int index1 = 0; index1 < meanVector.Size(); index1++)
      myfile << std::to_string(meanVector[index1]) + "\n";
    myfile.close();
    myfile.open(outputdirectory + "/zscore_std.csv");
    for (unsigned int index1 = 0; index1 < stdVector.Size(); index1++)
      myfile << std::to_string(stdVector[index1]) + "\n";
    myfile.close();
    std::cout << "Scaling parameters written." << std::endl;
    FinalResult = mTrainingSimulator.SplitTrainTest(scaledFeatureSet, LabelsOfAllSubjects, outputdirectory, classifiertype, foldtype);
  }
  else if (confType == CAPTK::ClassificationConfigurationType::CONF_TYPE_SPLIT_TRAIN)   //train
  {
    if (mTrainingSimulator.TrainData2(FeaturesOfAllSubjects, LabelsOfAllSubjects, outputdirectory, classifiertype,featureselectiontype,optimizationType, crossvalidationType) == true)
      std::cout << "Training finished successfully: Trained model saved in the specified output directory." << std::endl;
  }
  else   //test
  {
    std::cout << "testing code calling" << std::endl;
    FinalResult = mTrainingSimulator.TestData(FeaturesOfAllSubjects, modeldirectory, classifiertype, outputdirectory);
  }
  //std::cout << "Accuray=" << FinalResult[0] << std::endl;
  //std::cout << "Sensitivity=" << Finalresult[3] << std::endl;
  //std::cout << "Specificity=" << FinalResult[2] << std::endl;
  //std::cout << "Balanced Accuracy=" << FinalResult[3] << std::endl;

  //cbica::Logging(loggerFile, "Accuracy=" + std::to_string(FinalResult[0]) + "\n");
  //cbica::Logging(loggerFile, "Sensitivity=" + std::to_string(Finalresult[3]) + "\n");
  //cbica::Logging(loggerFile, "Specificity=" + std::to_string(FinalResult[2]) + "\n");
  //cbica::Logging(loggerFile, "Balanced Accuracy=" + std::to_string(FinalResult[3]) + "\n");

  return true;
}









VectorDouble TrainingModule::SplitTrainTest(const VariableSizeMatrixType inputFeatures, const VariableLengthVectorType inputLabels,
  const std::string outputfolder, const int classifiertype, const int training_size)
{
  MapType FoldingDataMap;
  std::vector<double> trainingindices;
  std::vector<double> traininglabels;
  std::vector<double> testingindices;
  std::vector<double> testinglabels;

  std::vector<double> predictedlabels;
  std::vector<double> predicteddistances;

  //make loops for training and validation, CVO partition like structure
  for (int index = 0; index < training_size; index++)
    trainingindices.push_back(index);

  for (unsigned int index = training_size; index < inputLabels.Size(); index++)
    testingindices.push_back(index);

  //std::cout << "testigindices" << testingindices.size() << std::endl;

  VariableSizeMatrixType trainingfeatures;
  VariableSizeMatrixType testingfeatures;

  trainingfeatures.SetSize(trainingindices.size(), inputFeatures.Cols());
  for (unsigned int i = 0; i < trainingindices.size(); i++)
  {
    traininglabels.push_back(inputLabels[trainingindices[i]]);
    for (unsigned int j = 0; j < trainingfeatures.Cols(); j++)
      trainingfeatures(i, j) = inputFeatures(trainingindices[i], j);
  }

  testingfeatures.SetSize(testingindices.size(), inputFeatures.Cols());
  for (unsigned int i = 0; i < testingindices.size(); i++)
  {
    testinglabels.push_back(inputLabels[testingindices[i]]);
    predictedlabels.push_back(-1);
    for (unsigned int j = 0; j < testingfeatures.Cols(); j++)
      testingfeatures(i, j) = inputFeatures(testingindices[i], j);
  }
  std::tuple<VectorDouble, VectorDouble, VariableSizeMatrixType, VectorDouble, VectorDouble, VariableSizeMatrixType, VectorDouble,VectorDouble> new_tuple(trainingindices, traininglabels, trainingfeatures, testingindices, testinglabels, testingfeatures, predictedlabels,predicteddistances);
  FoldingDataMap[0] = new_tuple;

  std::cout << "Folding map populated" << std::endl;

  //feature selection mechanism
  //---------------------------
  VectorDouble EffectSize = EffectSizeFeatureSelection(std::get<2>(FoldingDataMap[0]), std::get<1>(FoldingDataMap[0]));
  for (unsigned int eSizeCounter = 0; eSizeCounter < EffectSize.size(); eSizeCounter++)
  {
    if (EffectSize[eSizeCounter] < 0)
      EffectSize[eSizeCounter] = EffectSize[eSizeCounter] * -1;
  }
  std::vector<size_t> indices = sort_indexes(EffectSize);

  std::ofstream myfile;
  myfile.open(outputfolder + "/effctsizes.csv");
  for (unsigned int index1 = 0; index1 < EffectSize.size(); index1++)
    myfile << std::to_string(EffectSize[indices[index1]]) + "," + std::to_string(indices[index1]) + "\n";
  myfile.close();

  VariableSizeMatrixType CrossValidatedPerformances;
  CrossValidatedPerformances.SetSize(std::get<2>(FoldingDataMap[0]).Cols(), 8);
  for (unsigned int index1 = 0; index1 < CrossValidatedPerformances.Rows(); index1++)
    for (unsigned int index2 = 0; index2 < CrossValidatedPerformances.Cols(); index2++)
      CrossValidatedPerformances[index1][index2] = 0;

  for (unsigned int featureNo = 15; featureNo < std::get<2>(FoldingDataMap[0]).Cols(); featureNo++)
  {
    //std::cout << featureNo << std::endl;

    VariableSizeMatrixType reducedFeatureSet;
    reducedFeatureSet.SetSize(std::get<1>(FoldingDataMap[0]).size(), featureNo + 1);

    //copy the already selected features to the reduced feature set
    for (unsigned int j = 0; j < reducedFeatureSet.Rows(); j++)
      for (unsigned int k = 0; k < reducedFeatureSet.Cols(); k++)
        reducedFeatureSet(j, k) = std::get<2>(FoldingDataMap[0])(j, indices[k]);

    //check performance using cross-validation
    VectorDouble performance = InternalCrossValidation(reducedFeatureSet, traininglabels, 1, 0.01, 1);
    CrossValidatedPerformances(featureNo, 3) = performance[3];
  }
  std::cout << "Feature Selection Done!!!" << std::endl;

  for (unsigned int performanceNo = 17; performanceNo < CrossValidatedPerformances.Rows() - 2; performanceNo++)
  {
    CrossValidatedPerformances[performanceNo][2] = (CrossValidatedPerformances[performanceNo][3] +
      CrossValidatedPerformances[performanceNo - 1][3] +
      CrossValidatedPerformances[performanceNo - 2][3] +
      CrossValidatedPerformances[performanceNo + 1][3] +
      CrossValidatedPerformances[performanceNo + 2][3]) / 5;
  }
  double maxAveagePerformance = CrossValidatedPerformances[0][2];
  int maxFeatureNumber = -1;
  for (unsigned int performanceNo = 17; performanceNo < CrossValidatedPerformances.Rows(); performanceNo++)
  {
    if (CrossValidatedPerformances[performanceNo][2] >= maxAveagePerformance)
    {
      maxAveagePerformance = CrossValidatedPerformances[performanceNo][2];
      maxFeatureNumber = performanceNo;
    }
  }
  //std::cout << "maxFeatureNumber" << maxFeatureNumber << std::endl;
  //std::cout << "maxAveagePerformance:" << maxAveagePerformance << std::endl;

  //std::ofstream myfile;
  //myfile.open("E:/Projects/PSU/CrossValidationTrainingData.csv");
  //for (unsigned int index1 = 0; index1 < CrossValidatedPerformances.Rows(); index1++)
  //{
  //  for (unsigned int index2 = 0; index2 < CrossValidatedPerformances.Cols(); index2++)
  //  {
  //    if (index2 == 0)
  //      myfile << std::to_string(CrossValidatedPerformances[index1][index2]);
  //    else
  //      myfile << "," << std::to_string(CrossValidatedPerformances[index1][index2]);
  //  }
  //  myfile << "\n";
  //}
  //myfile.close();


  //Write the selected features
  myfile.open(outputfolder + "/selectedfeatures.csv");
  for (int index1 = 0; index1 < maxFeatureNumber; index1++)
    myfile << std::to_string(EffectSize[indices[index1]]) + "," + std::to_string(indices[index1]) + "\n";
  myfile.close();

  //Testing of the replicaiton cohort

  int featureNo = maxFeatureNumber;
  //std::cout << "No. of selected features:" << featureNo << std::endl;
  VariableSizeMatrixType reducedFeatureSet;
  reducedFeatureSet.SetSize(std::get<1>(FoldingDataMap[0]).size(), featureNo + 1);

  //copy the already selected features to the reduced feature set
  for (unsigned int j = 0; j < reducedFeatureSet.Rows(); j++)
    for (unsigned int k = 0; k < reducedFeatureSet.Cols(); k++)
    {
      reducedFeatureSet(j, k) = std::get<2>(FoldingDataMap[0])(j, indices[k]);
    }

  double bestCV = 0;
  double bestC = 1;
  for (double cValue = 1; cValue <= 5; cValue = cValue + 2)
  {
    VectorDouble result = InternalCrossValidation(reducedFeatureSet, traininglabels, pow(2, cValue), 0.01, classifiertype);
    float value = (int)(result[0] * 10 + .5);
    double currentcv = (float)value / 100;

    if (currentcv > bestCV)
    {
      bestC = pow(2, cValue);
      bestCV = currentcv;
      //std::cout << "best c: " << bestC << std::endl;
      //std::cout << "best cv: " << currentcv << std::endl;
    }
  }

  //model training and testing
  cv::Mat trainingData = cv::Mat::zeros(std::get<2>(FoldingDataMap[0]).Rows(), reducedFeatureSet.Cols(), CV_32FC1);
  cv::Mat trainingLabels = cv::Mat::zeros(std::get<2>(FoldingDataMap[0]).Rows(), 1, CV_32FC1);
  cv::Mat testingData = cv::Mat::zeros(std::get<5>(FoldingDataMap[0]).Rows(), reducedFeatureSet.Cols(), CV_32FC1);
  cv::Mat testingLabels = cv::Mat::zeros(std::get<5>(FoldingDataMap[0]).Rows(), 1, CV_32FC1);

  for (int copyDataCounter = 0; copyDataCounter < trainingData.rows; copyDataCounter++)
  {
    trainingLabels.ptr< float >(copyDataCounter)[0] = std::get<1>(FoldingDataMap[0])[copyDataCounter];
    for (int copyDataCounter2 = 0; copyDataCounter2 < trainingData.cols; copyDataCounter2++)
      trainingData.ptr< float >(copyDataCounter)[copyDataCounter2] = reducedFeatureSet[copyDataCounter][copyDataCounter2];
  }
  for (int copyDataCounter = 0; copyDataCounter < testingData.rows; copyDataCounter++)
  {
    testingLabels.ptr< float >(copyDataCounter)[0] = std::get<4>(FoldingDataMap[0])[copyDataCounter];
    for (int copyDataCounter2 = 0; copyDataCounter2 < testingData.cols; copyDataCounter2++)
      testingData.ptr< float >(copyDataCounter)[copyDataCounter2] = std::get<5>(FoldingDataMap[0])(copyDataCounter, indices[copyDataCounter2]);
  }
  ////--------------------just to write in files
  //VariableSizeMatrixType datatowrite;
  //datatowrite.SetSize(testingData.rows, testingData.cols);
  //for (int copyDataCounter = 0; copyDataCounter < testingData.rows; copyDataCounter++)
  //{
  //  for (int copyDataCounter2 = 0; copyDataCounter2 < testingData.cols; copyDataCounter2++)
  //  {
  //    datatowrite[copyDataCounter][copyDataCounter2] = std::get<5>(FoldingDataMap[0])(copyDataCounter, indices[copyDataCounter2]);
  //    //std::cout<<"Selected feature index: "<<indices[copyDataCounter2]<<std::endl;
  //  }
  //}
  //myfile.open("E:/Projects/PSU/ScaledTestingData.csv");
  //for(int index1=0;index1<datatowrite.Rows();index1++)
  //{
  //  for (int index2=0;index2<datatowrite.Cols();index2++)
  //  {
  //    if (index2 == 0)
  //      myfile << std::to_string(datatowrite[index1][index2]);
  //    else
  //      myfile << "," << std::to_string(datatowrite[index1][index2]);
  //  }
  //  myfile << "\n";
  //}
  //myfile.close();
  //------------------------------------------

  trainingLabels.convertTo(trainingLabels, CV_32SC1);
  auto svm = cv::ml::SVM::create();
  svm->setType(cv::ml::SVM::C_SVC);
  svm->setC(bestC);
  svm->setKernel(cv::ml::SVM::LINEAR);

  bool res = true;
  std::string msg;
  VectorDouble predictedLabels;
  VectorDouble predictedDistances;
  try
  {
    res = svm->train(trainingData, cv::ml::ROW_SAMPLE, trainingLabels);
    svm->save(outputfolder + "/SVM_Model.xml");
  }
  catch (cv::Exception ex)
  {
    msg = ex.what();
  }

  //apply SVM model on test data
  cv::Mat predicted(1, 1, CV_32F);
  for (int i = 0; i < testingData.rows; i++)
  {
    cv::Mat sample = testingData.row(i);
    svm->predict(sample, predicted, false);
    predictedLabels.push_back(predicted.ptr< float >(0)[0]);
    svm->predict(sample, predicted, true);
    predictedDistances.push_back(predicted.ptr< float >(0)[0]);
  }
  //calculate final performance and write in the csv file
  VectorDouble FinalPerformance = CalculatePerformanceMeasures(predictedLabels, std::get<4>(FoldingDataMap[0]));


  CrossValidatedPerformances(featureNo, 4) = FinalPerformance[0];
  CrossValidatedPerformances(featureNo, 5) = FinalPerformance[1];
  CrossValidatedPerformances(featureNo, 6) = FinalPerformance[2];
  CrossValidatedPerformances(featureNo, 7) = FinalPerformance[3];
  //std::cout << "predicted balanced" << FinalPerformance[3] << std::endl;



  myfile.open(outputfolder + "/predicted_distances.csv");
  for (unsigned int index1 = 0; index1 < predictedDistances.size(); index1++)
    myfile << std::to_string(std::abs(predictedDistances[index1])*predictedLabels[index1]) + "," + std::to_string(std::get<4>(FoldingDataMap[0])[index1]) + "\n";
  myfile.close();

  return FinalPerformance;
}



VectorDouble TrainingModule::TrainData(const VariableSizeMatrixType inputFeatures, const VariableLengthVectorType inputLabels,
  const std::string outputdirectory, const int classifiertype)
{

  //feature scaling mechanism
  //---------------------------
  std::cout << "Scaling started" << std::endl;
  FeatureScalingClass mFeaturesScaling;
  VariableSizeMatrixType scaledFeatureSet;
  VariableLengthVectorType meanVector;
  VariableLengthVectorType stdVector;
  mFeaturesScaling.ScaleGivenTrainingFeatures(inputFeatures, scaledFeatureSet, meanVector, stdVector);


  //remove the nan values
  for (unsigned int index1 = 0; index1 < scaledFeatureSet.Rows(); index1++)
  {
    for (unsigned int index2 = 0; index2 < scaledFeatureSet.Cols(); index2++)
    {
      if (std::isnan(scaledFeatureSet[index1][index2]))
        scaledFeatureSet[index1][index2] = 0;
    }
  }
  for (unsigned int index1 = 0; index1 < meanVector.Size(); index1++)
  {
    if (std::isnan(meanVector[index1]))
      meanVector[index1] = 0;
    if (std::isnan(stdVector[index1]))
      stdVector[index1] = 0;
  }
  std::cout << "Starting writing scaled parameters." << std::endl;
  std::ofstream myfile;
  myfile.open(outputdirectory + "/scaled_feature_set.csv");
  for (unsigned int index1 = 0; index1 < scaledFeatureSet.Rows(); index1++)
  {
    std::string onerow = "";
    for (unsigned int index2 = 0; index2 < scaledFeatureSet.Cols(); index2++)
    {
      if (index2 == 0)
        onerow = std::to_string(scaledFeatureSet[index1][index2]);
      else
        onerow = onerow + "," + std::to_string(scaledFeatureSet[index1][index2]);
    }
    onerow = onerow + "\n";
    myfile << onerow;
  }
  myfile.close();

  myfile.open(outputdirectory + "/zscore_mean.csv");
  for (unsigned int index1 = 0; index1 < meanVector.Size(); index1++)
    myfile << std::to_string(meanVector[index1]) + "\n";
  myfile.close();
  myfile.open(outputdirectory + "/zscore_std.csv");
  for (unsigned int index1 = 0; index1 < stdVector.Size(); index1++)
    myfile << std::to_string(stdVector[index1]) + "\n";
  myfile.close();
  std::cout << "Scaling parameters written." << std::endl;



  //copy the training labels
  //---------------------------
  std::vector<double> trainingindices;
  std::vector<double> traininglabels;
  for (unsigned int index = 0; index < inputLabels.Size(); index++)
    traininglabels.push_back(inputLabels[index]);
    

  //sorting based on effect sizes
  //-----------------------------
  VectorDouble EffectSize = EffectSizeFeatureSelection(inputFeatures, traininglabels);  //to convert to vector double
  for (unsigned int eSizeCounter = 0; eSizeCounter < EffectSize.size(); eSizeCounter++)
  {
    if (EffectSize[eSizeCounter] < 0)
      EffectSize[eSizeCounter] = EffectSize[eSizeCounter] * -1;
  }
  std::vector<size_t> indices = sort_indexes(EffectSize);

  myfile.open(outputdirectory + "/effctsizes.csv");
  for (unsigned int index1 = 0; index1 < EffectSize.size(); index1++)
    myfile << std::to_string(EffectSize[indices[index1]]) + "," + std::to_string(indices[index1]) + "\n";
  myfile.close();

  VariableSizeMatrixType CrossValidatedPerformances;
  CrossValidatedPerformances.SetSize(inputFeatures.Cols(), 8);
  for (unsigned int index1 = 0; index1 < CrossValidatedPerformances.Rows(); index1++)
    for (unsigned int index2 = 0; index2 < CrossValidatedPerformances.Cols(); index2++)
      CrossValidatedPerformances[index1][index2] = 0;

  //feature selection mechanism
  //---------------------------
  
  for (unsigned int featureNo = 0; featureNo < inputFeatures.Cols(); featureNo++)
  {
    //std::cout << featureNo << std::endl;

    VariableSizeMatrixType reducedFeatureSet;
    reducedFeatureSet.SetSize(traininglabels.size(), featureNo + 1);

    //copy the already selected features to the reduced feature set
    for (unsigned int j = 0; j < reducedFeatureSet.Rows(); j++)
      for (unsigned int k = 0; k < reducedFeatureSet.Cols(); k++)
        reducedFeatureSet(j, k) = inputFeatures(j, indices[k]);

    //check performance using cross-validation
    VectorDouble performance = InternalCrossValidation(reducedFeatureSet, traininglabels, 1, 0.01, 1);
    CrossValidatedPerformances(featureNo, 3) = performance[3];
  }
  std::cout << "Feature Selection Done!!!" << std::endl;

  for (unsigned int performanceNo = 17; performanceNo < CrossValidatedPerformances.Rows() - 2; performanceNo++)
  {
    CrossValidatedPerformances[performanceNo][2] = (CrossValidatedPerformances[performanceNo][3] +
      CrossValidatedPerformances[performanceNo - 1][3] +
      CrossValidatedPerformances[performanceNo - 2][3] +
      CrossValidatedPerformances[performanceNo + 1][3] +
      CrossValidatedPerformances[performanceNo + 2][3]) / 5;
  }
  double maxAveagePerformance = CrossValidatedPerformances[0][2];
  int maxFeatureNumber = -1;
  for (unsigned int performanceNo = 17; performanceNo < CrossValidatedPerformances.Rows(); performanceNo++)
  {
    if (CrossValidatedPerformances[performanceNo][2] >= maxAveagePerformance)
    {
      maxAveagePerformance = CrossValidatedPerformances[performanceNo][2];
      maxFeatureNumber = performanceNo;
    }
  }
  //std::cout << "maxFeatureNumber" << maxFeatureNumber << std::endl;
  //std::cout << "maxAveagePerformance:" << maxAveagePerformance << std::endl;

  //std::ofstream myfile;
  //myfile.open("E:/Projects/PSU/CrossValidationTrainingData.csv");
  //for (unsigned int index1 = 0; index1 < CrossValidatedPerformances.Rows(); index1++)
  //{
  //  for (unsigned int index2 = 0; index2 < CrossValidatedPerformances.Cols(); index2++)
  //  {
  //    if (index2 == 0)
  //      myfile << std::to_string(CrossValidatedPerformances[index1][index2]);
  //    else
  //      myfile << "," << std::to_string(CrossValidatedPerformances[index1][index2]);
  //  }
  //  myfile << "\n";
  //}
  //myfile.close();


  //Write the selected features
  myfile.open(outputdirectory + "/selectedfeatures.csv");
  for (int index1 = 0; index1 <= maxFeatureNumber; index1++)
    myfile << std::to_string(EffectSize[indices[index1]]) + "," + std::to_string(indices[index1]) + "\n";
  myfile.close();

  //Testing of the replicaiton cohort

  int featureNo = maxFeatureNumber;
  //std::cout << "No. of selected features:" << featureNo << std::endl;
  VariableSizeMatrixType reducedFeatureSet;
  reducedFeatureSet.SetSize(traininglabels.size(), featureNo + 1);

  //copy the already selected features to the reduced feature set
  for (unsigned int j = 0; j < reducedFeatureSet.Rows(); j++)
    for (unsigned int k = 0; k < reducedFeatureSet.Cols(); k++)
    {
      reducedFeatureSet(j, k) = inputFeatures(j, indices[k]);
    }

  double bestCV = 0;
  double bestC = 1;
  for (double cValue = 1; cValue <= 5; cValue = cValue + 2)
  {
    VectorDouble result = InternalCrossValidation(reducedFeatureSet, traininglabels, pow(2, cValue), 0.01, classifiertype);
    //float value = (int)(result[0] * 10 + .5);
    //double currentcv = (float)value / 100;

    if (result[0] > bestCV)
    {
      bestC = pow(2, cValue);
      bestCV = result[0];
      //std::cout << "best c: " << bestC << std::endl;
      //std::cout << "best cv: " << currentcv << std::endl;
    }
  }

  //model training and testing
  cv::Mat trainingData = cv::Mat::zeros(reducedFeatureSet.Rows(), reducedFeatureSet.Cols(), CV_32FC1);
  cv::Mat trainingLabels = cv::Mat::zeros(reducedFeatureSet.Rows(), 1, CV_32FC1);

  for (int copyDataCounter = 0; copyDataCounter < trainingData.rows; copyDataCounter++)
  {
    trainingLabels.ptr< float >(copyDataCounter)[0] =traininglabels[copyDataCounter];
    for (int copyDataCounter2 = 0; copyDataCounter2 < trainingData.cols; copyDataCounter2++)
      trainingData.ptr< float >(copyDataCounter)[copyDataCounter2] = reducedFeatureSet[copyDataCounter][copyDataCounter2];
  }
  ////--------------------just to write in files
  //VariableSizeMatrixType datatowrite;
  //datatowrite.SetSize(testingData.rows, testingData.cols);
  //for (int copyDataCounter = 0; copyDataCounter < testingData.rows; copyDataCounter++)
  //{
  //  for (int copyDataCounter2 = 0; copyDataCounter2 < testingData.cols; copyDataCounter2++)
  //  {
  //    datatowrite[copyDataCounter][copyDataCounter2] = std::get<5>(FoldingDataMap[0])(copyDataCounter, indices[copyDataCounter2]);
  //    //std::cout<<"Selected feature index: "<<indices[copyDataCounter2]<<std::endl;
  //  }
  //}
  //myfile.open("E:/Projects/PSU/ScaledTestingData.csv");
  //for(int index1=0;index1<datatowrite.Rows();index1++)
  //{
  //  for (int index2=0;index2<datatowrite.Cols();index2++)
  //  {
  //    if (index2 == 0)
  //      myfile << std::to_string(datatowrite[index1][index2]);
  //    else
  //      myfile << "," << std::to_string(datatowrite[index1][index2]);
  //  }
  //  myfile << "\n";
  //}
  //myfile.close();
  //------------------------------------------

  trainingLabels.convertTo(trainingLabels, CV_32SC1);
  auto svm = cv::ml::SVM::create();
  svm->setType(cv::ml::SVM::C_SVC);
  svm->setC(bestC);
  svm->setKernel(cv::ml::SVM::LINEAR);

  bool res = true;
  std::string msg;
  VectorDouble predictedLabels;
  VectorDouble predictedDistances;
  try
  {
    res = svm->train(trainingData, cv::ml::ROW_SAMPLE, trainingLabels);
    svm->save(outputdirectory + "/SVM_Model.xml");
  }
  catch (cv::Exception ex)
  {
    msg = ex.what();
  }

  ////apply SVM model on test data
  //cv::Mat predicted(1, 1, CV_32F);
  //for (int i = 0; i < testingData.rows; i++)
  //{
  //  cv::Mat sample = testingData.row(i);
  //  svm->predict(sample, predicted, false);
  //  predictedLabels.push_back(predicted.ptr< float >(0)[0]);
  //  svm->predict(sample, predicted, true);
  //  predictedDistances.push_back(predicted.ptr< float >(0)[0]);
  //}
  ////calculate final performance and write in the csv file
  VectorDouble FinalPerformance;
  //= CalculatePerformanceMeasures(predictedLabels, std::get<4>(FoldingDataMap[0]));


  //CrossValidatedPerformances(featureNo, 4) = FinalPerformance[0];
  //CrossValidatedPerformances(featureNo, 5) = FinalPerformance[1];
  //CrossValidatedPerformances(featureNo, 6) = FinalPerformance[2];
  //CrossValidatedPerformances(featureNo, 7) = FinalPerformance[3];
  ////std::cout << "predicted balanced" << FinalPerformance[3] << std::endl;



  //myfile.open(outputfolder + "/predicted_distances.csv");
  //for (unsigned int index1 = 0; index1 < predictedDistances.size(); index1++)
  //  myfile << std::to_string(std::abs(predictedDistances[index1])*predictedLabels[index1]) + "," + std::to_string(std::get<4>(FoldingDataMap[0])[index1]) + "\n";
  //myfile.close();

  return FinalPerformance;
}

bool TrainingModule::TrainData2(const VariableSizeMatrixType inputFeatures, const VariableLengthVectorType inputLabels,
  const std::string outputdirectory, const int classifiertype,const int featureselectiontype,
  const int optimizationType, const int crossvalidationType)
{
  //scaling of the given features
  //---------------------------
  std::cout << "Scaling started" << std::endl;
  FeatureScalingClass mFeaturesScaling;
  VariableSizeMatrixType scaledFeatureSet;
  VariableLengthVectorType meanVector;
  VariableLengthVectorType stdVector;
  mFeaturesScaling.ScaleGivenTrainingFeatures(inputFeatures, scaledFeatureSet, meanVector, stdVector);

  //remove the nan values introduced after the scaling process
  for (unsigned int index1 = 0; index1 < scaledFeatureSet.Rows(); index1++)
  {
    for (unsigned int index2 = 0; index2 < scaledFeatureSet.Cols(); index2++)
    {
      if (std::isnan(scaledFeatureSet[index1][index2]))
        scaledFeatureSet[index1][index2] = 0;
    }
  }
  for (unsigned int index1 = 0; index1 < meanVector.Size(); index1++)
  {
    if (std::isnan(meanVector[index1]))
      meanVector[index1] = 0;
    if (std::isnan(stdVector[index1]))
      stdVector[index1] = 0;
  }
  WriteCSVFiles(scaledFeatureSet,outputdirectory + "/scaled_feature_set.csv");
  WriteCSVFiles(meanVector,outputdirectory + "/zscore_mean.csv");
  WriteCSVFiles(stdVector, outputdirectory + "/zscore_std.csv");
  std::cout << "Scaling parameters written." << std::endl;
  
  //copy the training labels
  //---------------------------
  std::vector<double> trainingindices;
  std::vector<double> traininglabels;
  for (unsigned int index = 0; index < inputLabels.Size(); index++)
    traininglabels.push_back(inputLabels[index]);

  std::vector<int> FinalSelectedFeatures;
  VectorDouble crossvalidatedAccuracies; 
  if (featureselectiontype == 1)
    FinalSelectedFeatures = EffectSizeBasedFeatureSelection(scaledFeatureSet,traininglabels,classifiertype, optimizationType,crossvalidationType, crossvalidatedAccuracies);
  //else if (featureselectiontype == 2)
  //  FinalSelectedFeatures = CorrelationBasedFeatureSelection(scaledFeatureSet,traininglabels,classifiertype); 
  else if (featureselectiontype == 3)
  FinalSelectedFeatures = SVMFFSBasedFeatureSelection(scaledFeatureSet, traininglabels, classifiertype, optimizationType,crossvalidationType, crossvalidatedAccuracies);
  //else if (featureselectiontype == 4)
  //  FinalSelectedFeatures = SVMRFEBasedFeatureSelection(scaledFeatureSet,traininglabels,classifiertype);
  else
  {
    std::cout << "The selected feature selection method is not avaialble inside CaPTk. Please select another method..." << std::endl;
    return false;
  }


  //Write the selected features and crossvalidated accuracies
  WriteCSVFiles(FinalSelectedFeatures,outputdirectory + "/SelectedFeatures.csv");
  WriteCSVFiles(crossvalidatedAccuracies, outputdirectory + "/CrossValidatedAccuracies.csv");

  //optimize classifiers parameters on the selected feature set
  VariableSizeMatrixType FinalSelectedFeatureSet;
  FinalSelectedFeatureSet.SetSize(traininglabels.size(), FinalSelectedFeatures.size());

  //copy the already selected features to the current feature set
  for (unsigned int j = 0; j < FinalSelectedFeatureSet.Rows(); j++)
    for (unsigned int k = 0; k < FinalSelectedFeatures.size(); k++)
      FinalSelectedFeatureSet(j, k) = scaledFeatureSet(j, FinalSelectedFeatures[k]);

  double bestCV = 0;
  double bestC = 1;
  double bestG = 1 / FinalSelectedFeatures.size();
  if (classifiertype == 2)
  {
    VectorDouble result;
    for (double cValue = -5; cValue <= 5; cValue = cValue + 1)
      for (double gValue = -5; gValue <= 5; gValue = gValue + 1)
      {
        VectorDouble result = InternalCrossValidation(FinalSelectedFeatureSet, traininglabels, pow(2, cValue), pow(2, gValue), classifiertype);
        if (result[3] > bestCV)
        {
          bestC = pow(2, cValue);
          bestG = pow(2, gValue);
          bestCV = result[3];
        }
      }
  }
  else
  {
    for (double cValue = -5; cValue <= 5; cValue = cValue + 1)
    {
      VectorDouble result = InternalCrossValidation(FinalSelectedFeatureSet, traininglabels, pow(2, cValue), 0.01, classifiertype);
      if (result[3] > bestCV)
      {
        bestC = pow(2, cValue);
        bestCV = result[3];
      }
    }
  }
  std::cout << "Optimal C=" << bestC << std::endl;
  std::cout << "Optimal gamma=" << bestG << std::endl;

  //model training and testing
  cv::Mat trainingData = cv::Mat::zeros(FinalSelectedFeatureSet.Rows(), FinalSelectedFeatureSet.Cols(), CV_32FC1);
  cv::Mat trainingLabels = cv::Mat::zeros(FinalSelectedFeatureSet.Rows(), 1, CV_32FC1);

  for (int copyDataCounter = 0; copyDataCounter < trainingData.rows; copyDataCounter++)
  {
    trainingLabels.ptr< float >(copyDataCounter)[0] = traininglabels[copyDataCounter];
    for (int copyDataCounter2 = 0; copyDataCounter2 < trainingData.cols; copyDataCounter2++)
      trainingData.ptr< float >(copyDataCounter)[copyDataCounter2] = FinalSelectedFeatureSet[copyDataCounter][copyDataCounter2];
  }

  trainingLabels.convertTo(trainingLabels, CV_32SC1);

  //make an SVM model
  auto svm = cv::ml::SVM::create();
  svm->setType(cv::ml::SVM::C_SVC);
  svm->setC(bestC);

  if (classifiertype == 2)
  {
    svm->setGamma(bestG);
    svm->setKernel(cv::ml::SVM::RBF);
  }
  else
    svm->setKernel(cv::ml::SVM::LINEAR);

  bool res;
  std::string msg;
  VectorDouble predictedLabels;
  VectorDouble predictedDistances;
  try
  {
    res = svm->train(trainingData, cv::ml::ROW_SAMPLE, trainingLabels);
    svm->save(outputdirectory + "/SVM_Model.xml");
  }
  catch (cv::Exception ex)
  {
    msg = ex.what();
  }
  //apply SVM model on training data
  cv::Mat predicted(1, 1, CV_32F);
  for (int i = 0; i < trainingData.rows; i++)
  {
    cv::Mat sample = trainingData.row(i);
    svm->predict(sample, predicted, false);
    predictedLabels.push_back(predicted.ptr< float >(0)[0]);
    svm->predict(sample, predicted, true);
    predictedDistances.push_back(predicted.ptr< float >(0)[0]);
  }
  for (int i = 0; i < predictedLabels.size(); i++)
    predictedDistances[i] = std::abs(predictedDistances[i])*predictedLabels[i];
  //Write the selected features and crossvalidated accuracies
  WriteCSVFiles(predictedLabels, outputdirectory + "/predicted_labels_training.csv");
  WriteCSVFiles(predictedDistances, outputdirectory + "/predicted_distances_training.csv");

  return res;
}


VectorDouble TrainingModule::TestData(const VariableSizeMatrixType inputFeatures,  const std::string modeldirectory, const int classifiertype, const std::string outputfolder)
{
  typedef itk::CSVArray2DFileReader<double> ReaderType;
  VectorDouble results;
  CSVFileReaderType::Pointer reader = CSVFileReaderType::New();

  MatrixType meanMatrix;
  VariableLengthVectorType mean;
  VariableLengthVectorType stddevition;
  try
  {
    reader->SetFileName(modeldirectory + "/zscore_mean.csv");
    reader->SetFieldDelimiterCharacter(',');
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();
    reader->Parse();
    meanMatrix = reader->GetArray2DDataObject()->GetMatrix();

    mean.SetSize(meanMatrix.size());
    for (unsigned int i = 0; i < meanMatrix.size(); i++)
      mean[i] = meanMatrix(0,i);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Error in reading the file: " + modeldirectory + "/zscore_mean.csv. Error code : " + std::string(e1.what()));
  }

  MatrixType stdMatrix;
  try
  {
    reader->SetFileName(modeldirectory + "/zscore_std.csv");
    reader->SetFieldDelimiterCharacter(',');
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();
    reader->Parse();
    stdMatrix = reader->GetArray2DDataObject()->GetMatrix();

    stddevition.SetSize(stdMatrix.size());
    for (unsigned int i = 0; i < stdMatrix.size(); i++)
      stddevition[i] = stdMatrix(0,i);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Error in reading the file: " + modeldirectory + "/zscore_std.csv. Error code : " + std::string(e1.what()));
  }

  //feature selection process for test data
  std::vector<double> selectedfeaturesindices;
  MatrixType dataMatrix;
  try
  {
    reader->SetFileName(modeldirectory + "/SelectedFeatures.csv");
    reader->SetFieldDelimiterCharacter(',');
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();
    reader->Parse();
    dataMatrix = reader->GetArray2DDataObject()->GetMatrix();

    for (unsigned int i = 0; i < dataMatrix.cols(); i++)
      selectedfeaturesindices.push_back(dataMatrix(0,i));
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Error in reading the file: " + modeldirectory + "/selectedfeatures.csv. Error code : " + std::string(e1.what()));
  }


  FeatureScalingClass mFeatureScalingLocalPtr;
  VariableSizeMatrixType ScaledTestingData = mFeatureScalingLocalPtr.ScaleGivenTestingFeatures(inputFeatures, mean, stddevition);

  std::cout << "Testing data scaled." << std::endl;
  VariableSizeMatrixType FinalSelectedFeatureSet;
  FinalSelectedFeatureSet.SetSize(ScaledTestingData.Rows(), selectedfeaturesindices.size());

  //copy the already selected features to the current feature set
  for (unsigned int j = 0; j < ScaledTestingData.Rows(); j++)
    for (unsigned int k = 0; k < FinalSelectedFeatureSet.Cols(); k++)
      FinalSelectedFeatureSet(j, k) = ScaledTestingData(j, selectedfeaturesindices[k]);
  std::cout << "Testing data features selected." << std::endl;

  MapType FoldingDataMap;
  //read from csv the selected features
  std::vector<double> testingindices;
  std::vector<double> testinglabels;

  std::vector<double> predictedlabels;

  auto svm = cv::Algorithm::load<cv::ml::SVM>(modeldirectory +"/SVM_Model.xml" );
  cv::Mat testingData = cv::Mat::zeros(FinalSelectedFeatureSet.Rows(), FinalSelectedFeatureSet.Cols(), CV_32FC1);
  cv::Mat testingLabels = cv::Mat::zeros(FinalSelectedFeatureSet.Rows(), 1, CV_32FC1);

  for (int copyDataCounter = 0; copyDataCounter < testingData.rows; copyDataCounter++)
  {
    testingLabels.ptr< float >(copyDataCounter)[0] = 1;
    for (int copyDataCounter2 = 0; copyDataCounter2 < testingData.cols; copyDataCounter2++)
      testingData.ptr< float >(copyDataCounter)[copyDataCounter2] = FinalSelectedFeatureSet(copyDataCounter, copyDataCounter2);
  }
  //////--------------------just to write in files
  std::ofstream myfile;
  VariableSizeMatrixType datatowrite;
  datatowrite.SetSize(testingData.rows, testingData.cols);
  for (int copyDataCounter = 0; copyDataCounter < testingData.rows; copyDataCounter++)
  {
    for (int copyDataCounter2 = 0; copyDataCounter2 < testingData.cols; copyDataCounter2++)
    {
      datatowrite[copyDataCounter][copyDataCounter2] = FinalSelectedFeatureSet(copyDataCounter, copyDataCounter2);
      //std::cout<<"Selected feature index: "<<indices[copyDataCounter2]<<std::endl;
    }
  }
  myfile.open(outputfolder + "/ScaledTestingData.csv");
  for(unsigned int index1=0;index1<datatowrite.Rows();index1++)
  {
    for (unsigned int index2=0;index2<datatowrite.Cols();index2++)
    {
      if (index2 == 0)
        myfile << std::to_string(datatowrite[index1][index2]);
      else
        myfile << "," << std::to_string(datatowrite[index1][index2]);
    }
    myfile << "\n";
  }
  myfile.close();
  ////------------------------------------------

  //trainingLabels.convertTo(trainingLabels, CV_32SC1);
  //auto svm = cv::ml::SVM::create();
  //svm->setType(cv::ml::SVM::C_SVC);
  //svm->setC(bestC);
  //svm->setKernel(cv::ml::SVM::LINEAR);

  //bool res = true;
  //std::string msg;
  VectorDouble predictedLabels;
  VectorDouble predictedDistances;
  //try
  //{
  //  res = svm->train(trainingData, cv::ml::ROW_SAMPLE, trainingLabels);
  //  svm->save(outputfolder + "/SVM_Model.xml");
  //}
  //catch (cv::Exception ex)
  //{
  //  msg = ex.what();
  //}
  

  //apply SVM model on test data
  cv::Mat predicted(1, 1, CV_32F);
  for (int i = 0; i < testingData.rows; i++)
  {
    cv::Mat sample = testingData.row(i);
    svm->predict(sample, predicted, false);
    predictedLabels.push_back(predicted.ptr< float >(0)[0]);
    svm->predict(sample, predicted, true);
    predictedDistances.push_back(predicted.ptr< float >(0)[0]);
  }
  //calculate final performance and write in the csv file
  //VectorDouble FinalPerformance = CalculatePerformanceMeasures(predictedLabels, std::get<4>(FoldingDataMap[0]));


  //CrossValidatedPerformances(featureNo, 4) = FinalPerformance[0];
  //CrossValidatedPerformances(featureNo, 5) = FinalPerformance[1];
  //CrossValidatedPerformances(featureNo, 6) = FinalPerformance[2];
  //CrossValidatedPerformances(featureNo, 7) = FinalPerformance[3];
  //std::cout << "predicted balanced" << FinalPerformance[3] << std::endl;



  myfile.open(outputfolder + "/predicted_distances.csv");
  for (unsigned int index1 = 0; index1 < predictedDistances.size(); index1++)
    myfile << std::to_string(std::abs(predictedDistances[index1])*predictedLabels[index1]) + "\n";
  myfile.close();
  myfile.open(outputfolder + "/predicted_labels.csv");
  for (unsigned int index1 = 0; index1 < predictedDistances.size(); index1++)
    myfile << std::to_string(predictedLabels[index1]) + "\n";
  myfile.close();

  return results;
}
//----------------------------------------------------------------------------------------------------------
std::vector<int> TrainingModule::UpdateUnselectedFeatures(std::vector<int> SelectedFeatures,int featuresize)
{
  std::vector<int> UnselectedFeatures;
  for (unsigned int featureCounter = 0; featureCounter < featuresize; featureCounter++)
  {
    int found = 0;
    for (unsigned int index2 = 0; index2 < SelectedFeatures.size(); index2++)
    {
      if (SelectedFeatures[index2] == featureCounter)
      {
        found = 1;
        break;
      }
    }
    if (found == 0)
      UnselectedFeatures.push_back(featureCounter);
  }
  return UnselectedFeatures;
}
//-----------------------------------------------------------------------------------------------------------
std::vector<int> TrainingModule::SVMFFSBasedFeatureSelection(const VariableSizeMatrixType inputdata, const VectorDouble labels, const int classifiertype,const int optimizationtype,const int cvtype,VectorDouble &crossvalidatedaccuracies)
{
  std::vector<int> SelectedFeatures;
  std::vector<double> CrossValidatedBalancedAccuraciesFinal;
  std::vector<int> UnselectedFeatures = UpdateUnselectedFeatures(SelectedFeatures, inputdata.Cols());

  //feature selection mechanism
  //---------------------------
  while (UnselectedFeatures.size() > 0)
  {
    std::vector<double> CrossValidatedBalancedAccuraciesAfterAddingUnselectedFeatures;
    for (unsigned int featureNo = 0; featureNo < UnselectedFeatures.size(); featureNo++)
    {
      VariableSizeMatrixType reducedFeatureSet;
      reducedFeatureSet.SetSize(labels.size(), SelectedFeatures.size() + 1);

      //copy the already selected features to the current feature set
      for (unsigned int j = 0; j < reducedFeatureSet.Rows(); j++)
      {
        for (unsigned int k = 0; k < SelectedFeatures.size(); k++)
        {
          double value = (int)(inputdata(j, SelectedFeatures[k]) * 10000 + .5);
          reducedFeatureSet(j, k) = (double)value / 10000;
        }
      }

      //copy the new feature to the current feature set
      for (unsigned int j = 0; j < reducedFeatureSet.Rows(); j++)
      {
        double value = (int)(inputdata(j, UnselectedFeatures[featureNo]) * 10000 + .5);
        reducedFeatureSet(j, reducedFeatureSet.Cols() - 1) = (double)value / 10000;
      }

      //check crossvalidated performance after adding the current feature
      double bestCV = 0;
      double bestC = 1;
      double bestG = 1 / reducedFeatureSet.Cols();
      if (optimizationtype == 1)
      {
        if (classifiertype == 2)
        {
          for (double cValue = -5; cValue <= 5; cValue = cValue + 1)
            for (double gValue = -5; gValue <= 5; gValue = gValue + 1)
            {
              VectorDouble result;
              if (cvtype == 1)
                result = InternalCrossValidationResubstitution(reducedFeatureSet, labels, pow(2, cValue), pow(2, gValue), classifiertype);
              else
                result = InternalCrossValidation(reducedFeatureSet, labels, pow(2, cValue), pow(2, gValue), classifiertype);

              if (result[3] > bestCV)
                bestCV = result[3];
            }
        }
        else
        {
          for (double cValue = -5; cValue <=5; cValue = cValue + 1)
          {
            VectorDouble result;
            if (cvtype == 1)
              result = InternalCrossValidationResubstitution(reducedFeatureSet, labels, pow(2, cValue), 0.01, classifiertype);
            else
              result = InternalCrossValidation(reducedFeatureSet, labels, pow(2, cValue), 0.01, classifiertype);

            if (result[3] > bestCV)
              bestCV = result[3];
          }
        }
      }
      else
      {
        VectorDouble result;
        if (cvtype == 1)
          result = InternalCrossValidationResubstitution(reducedFeatureSet, labels, bestC, bestG, classifiertype);
        else
          result = InternalCrossValidation(reducedFeatureSet, labels, bestC, bestG, classifiertype);
        bestCV = result[3];
      }
      CrossValidatedBalancedAccuraciesAfterAddingUnselectedFeatures.push_back(bestCV);
    }
    //sort cross-validated balanced accuracies and pick the best selected feature
    int index = std::distance(CrossValidatedBalancedAccuraciesAfterAddingUnselectedFeatures.begin(), std::max_element(CrossValidatedBalancedAccuraciesAfterAddingUnselectedFeatures.begin(), CrossValidatedBalancedAccuraciesAfterAddingUnselectedFeatures.end()));
    CrossValidatedBalancedAccuraciesFinal.push_back(CrossValidatedBalancedAccuraciesAfterAddingUnselectedFeatures[index]);
    
    //check if we have three consecutive 1's 
    int sizeofcvaccuracies = CrossValidatedBalancedAccuraciesFinal.size();
    if (sizeofcvaccuracies >= 3)
    {
      if (CrossValidatedBalancedAccuraciesFinal[sizeofcvaccuracies - 1] == 1 && CrossValidatedBalancedAccuraciesFinal[sizeofcvaccuracies - 2] == 1 && CrossValidatedBalancedAccuraciesFinal[sizeofcvaccuracies - 3] == 1)
        break;
    }
    //if (sizeofcvaccuracies >= 2)
    //{
    //  if (CrossValidatedBalancedAccuraciesFinal[sizeofcvaccuracies - 1]< CrossValidatedBalancedAccuraciesFinal[sizeofcvaccuracies - 2] )
    //    break;
    //}

    SelectedFeatures.push_back(UnselectedFeatures[index]);
    std::cout << "CurrentSize=" << SelectedFeatures.size() << " Performance=" << CrossValidatedBalancedAccuraciesAfterAddingUnselectedFeatures[index] << std::endl;
    //for (int counter = 0; counter < SelectedFeatures.size(); counter++)
    //  std::cout << SelectedFeatures[counter] << std::endl;
    UnselectedFeatures = UpdateUnselectedFeatures(SelectedFeatures, inputdata.Cols());
  }

  std::cout << "Feature Selection Done!!!" << std::endl;
  crossvalidatedaccuracies = CrossValidatedBalancedAccuraciesFinal;

  //we are doing moving average to avoid local maxima. in pairs of 3, we average them and pick the middle one. 
  //TODO: consider selecting the first feature out of the three features involved in moving average
  VariableLengthVectorType MovingAverageOnCrossValidatedPerformance;
  MovingAverageOnCrossValidatedPerformance.SetSize(CrossValidatedBalancedAccuraciesFinal.size());
  for (int index = 0; index < MovingAverageOnCrossValidatedPerformance.Size(); index++)
    MovingAverageOnCrossValidatedPerformance[index] = 0;

  for (unsigned int index = 1; index < CrossValidatedBalancedAccuraciesFinal.size() - 1; index++)
    MovingAverageOnCrossValidatedPerformance[index] = (CrossValidatedBalancedAccuraciesFinal[index] + CrossValidatedBalancedAccuraciesFinal[index - 1] + CrossValidatedBalancedAccuraciesFinal[index + 1]) / 3;

  int max_performance_counter = std::distance(CrossValidatedBalancedAccuraciesFinal.begin(), std::max_element(CrossValidatedBalancedAccuraciesFinal.begin(), CrossValidatedBalancedAccuraciesFinal.end()));
  std::vector<int> FinalSelectedFeatures;
  for (int index = 0; index <= max_performance_counter; index++)
    FinalSelectedFeatures.push_back(SelectedFeatures[index]);
  return FinalSelectedFeatures;
}
std::vector<int> TrainingModule::EffectSizeBasedFeatureSelection(const VariableSizeMatrixType inputdata, const VectorDouble labels,const int classifiertype,const int optimizationtype, const int cvtype,VectorDouble & crossvalidatedaccuracies)
{
  double numberOfSelectedFeatures = 0;
  VectorDouble EffectSize = EffectSizeFeatureSelection(inputdata, labels);
  for (unsigned int eSizeCounter = 0; eSizeCounter < EffectSize.size(); eSizeCounter++)
  {
    if (EffectSize[eSizeCounter] < 0)
      EffectSize[eSizeCounter] = EffectSize[eSizeCounter] * -1;
  }
  std::vector<size_t> indices = sort_indexes(EffectSize);

  //keep on copy the  features
  std::vector<double> CrossValidatedBalancedAccuracies;
  std::vector<double> PerFeaturePerformance;
  for (unsigned int featureNo = 0; featureNo < EffectSize.size(); featureNo++)
  {
    VariableSizeMatrixType currentFeatureSet;
    currentFeatureSet.SetSize(labels.size(), featureNo + 1);

    //copy the current feature to the feature set
    for (unsigned int j = 0; j < currentFeatureSet.Rows(); j++)
      for (unsigned int k = 0; k < currentFeatureSet.Cols(); k++)
        currentFeatureSet(j, k) = inputdata(j, indices[k]);

    //check crossvalidated performance after adding the current feature
    double bestCV = 0;
    double bestC = 1;
    double numerator = 1;
    double bestG = numerator / ((double)currentFeatureSet.Cols());
    if (optimizationtype == 1)
    {
      if (classifiertype == 2)
      {
        //for now, we have constant ranges to search optimal values of C and Gamma parameters
        //should be changed in future to get these ranges from users
        for (double cValue = -5; cValue <= 5; cValue = cValue + 1)
          for (double gValue = -5; gValue <= 5; gValue = gValue + 1)
          {
            VectorDouble result;
            if (cvtype == 1)
              result = InternalCrossValidationResubstitution(currentFeatureSet, labels, pow(2, cValue), pow(2, gValue), classifiertype);
            else
              result = InternalCrossValidation(currentFeatureSet, labels, pow(2, cValue), pow(2, gValue), classifiertype);

            if (result[3] > bestCV)
              bestCV = result[3];
          }
      }
      else
      {
        for (double cValue = -5; cValue <= 5; cValue = cValue + 1)
        {
          VectorDouble result;
          if (cvtype == 1)
            result = InternalCrossValidationResubstitution(currentFeatureSet, labels, pow(2, cValue), 0.01, classifiertype);
          else
            result = InternalCrossValidation(currentFeatureSet, labels, pow(2, cValue), 0.01, classifiertype);

          if (result[3] > bestCV)
            bestCV = result[3];
        }
      }
    }
    else
    {
      VectorDouble result;
      if (cvtype == 1)
        result = InternalCrossValidationResubstitution(currentFeatureSet, labels, bestC, bestG, classifiertype);
      else
        result = InternalCrossValidation(currentFeatureSet, labels, bestC, bestG, classifiertype);
      bestCV = result[3];
    }
    CrossValidatedBalancedAccuracies.push_back(bestCV);
  }
  crossvalidatedaccuracies = CrossValidatedBalancedAccuracies;
  std::cout << "Feature Selection Done!" << std::endl;
  //we are doing moving average to avoid local maxima. in pairs of 3, we average them and pick the middle one. 
  //TODO: consider selecting the first feature out of the three features involved in moving average

  std::vector<double> MovingAverageOnCrossValidatedPerformance;
  MovingAverageOnCrossValidatedPerformance.push_back(0);  //can not use first index for averaging
  for (unsigned int index = 1; index < CrossValidatedBalancedAccuracies.size() - 1; index++)
    MovingAverageOnCrossValidatedPerformance.push_back((CrossValidatedBalancedAccuracies[index] + CrossValidatedBalancedAccuracies[index - 1] + CrossValidatedBalancedAccuracies[index + 1]) / 3);
  MovingAverageOnCrossValidatedPerformance.push_back(0);  //can not use last index for averaging

  int max_performance_counter = std::distance(MovingAverageOnCrossValidatedPerformance.begin(), std::max_element(MovingAverageOnCrossValidatedPerformance.begin(), MovingAverageOnCrossValidatedPerformance.end()));
  std::vector<int> FinalSelectedFeatures;
  for (int index = 0; index <= max_performance_counter; index++)
    FinalSelectedFeatures.push_back(indices[index]);
  std::cout << "No. of selected features!" << FinalSelectedFeatures.size() << std::endl;
  return FinalSelectedFeatures;
}
