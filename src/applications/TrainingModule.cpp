/**
\file  TrainingSimulator.cpp

\brief The source file containing the SurvivaPredictor class, used to find Survival Prediction Index
Author: Saima Rathore
Library Dependecies: ITK 4.7+ <br>

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2016 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

#include "TrainingModule.h"
#include "opencv2/core/core.hpp"
#include "opencv2/ml.hpp"


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

VectorDouble TrainingModule::InternalCrossValidation(VariableSizeMatrixType inputFeatures, std::vector<double> inputLabels, double cValue, double gValue,int kerneltype)
{
  VariableLengthVectorType predictedLabels;
  predictedLabels.SetSize(inputLabels.size());


  int fold_size = inputFeatures.Rows() / 10;
  //make loops for training and validation
  for (int index = 0; index < 10; index++)
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
      x= floor(x * 100 + 0.5) / 100;
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

VectorDouble TrainingModule::CrossValidation(const VariableSizeMatrixType inputFeatures, const VariableLengthVectorType inputLabels,const std::string outputfolder,const int classifiertype,const int number_of_folds)
{
  MapType FoldingDataMap;
  int fold_size = inputLabels.Size() / number_of_folds;
  int remainder = inputLabels.Size() % number_of_folds;

  //make loops for training and validation, CVO partition like structure
  for (int index = 0; index < number_of_folds; index++)
  {
    //find training and testing indices
    std::vector<double> trainingindices;
    std::vector<double> traininglabels;
    VariableSizeMatrixType trainingfeatures;

    std::vector<double> testingindices;
    std::vector<double> testinglabels;
    std::vector<double> predictedlabels;
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
    //find training and testing labels and features
    trainingfeatures.SetSize(trainingindices.size(), inputFeatures.Cols());
    for (int i = 0; i < trainingindices.size(); i++)
    {
      traininglabels.push_back(inputLabels[trainingindices[i]]);
      for (unsigned int j = 0; j < trainingfeatures.Cols(); j++)
        trainingfeatures(i, j) = inputFeatures(trainingindices[i], j);
    }
    testingfeatures.SetSize(testingindices.size(), inputFeatures.Cols());
    for (int i = 0; i < testingindices.size(); i++)
    {
      testinglabels.push_back(inputLabels[testingindices[i]]);
      predictedlabels.push_back(-1);
      for (unsigned int j = 0; j < testingfeatures.Cols(); j++)
        testingfeatures(i, j) = inputFeatures(testingindices[i], j);
    }
    std::tuple<VectorDouble, VectorDouble, VariableSizeMatrixType, VectorDouble, VectorDouble, VariableSizeMatrixType, VectorDouble> new_tuple(trainingindices, traininglabels, trainingfeatures, testingindices, testinglabels, testingfeatures, predictedlabels);
    FoldingDataMap[index] = new_tuple;
  }
  for (int index = 0; index < number_of_folds; index++)
  {
    std::cout << "***************************** index=" << index << "*********************************" << std::endl << std::endl;
    // cbica::Logging(loggerFile, "***************************** index=" + std::to_string(index) + "*********************************\n\n");

    //feature selection mechanism
    //---------------------------
    double numberOfSelectedFeatures = 0;
    VectorDouble EffectSize = EffectSizeFeatureSelection(std::get<2>(FoldingDataMap[index]), std::get<1>(FoldingDataMap[index]));
    for (int eSizeCounter = 0; eSizeCounter < EffectSize.size(); eSizeCounter++)
    {
      if (EffectSize[eSizeCounter] < 0)
        EffectSize[eSizeCounter] = EffectSize[eSizeCounter] * -1;
    }

    std::vector<size_t> indices = sort_indexes(EffectSize);

    //keep on copy the  features
    std::vector<double> PerFeaturePerformance;
    for (int featureNo = 0; featureNo < EffectSize.size(); featureNo++)
    {
      VariableSizeMatrixType reducedFeatureSet;
      reducedFeatureSet.SetSize(std::get<1>(FoldingDataMap[index]).size(), featureNo + 1);

      //copy the reduced features to the reduced feature set
      for (unsigned int j = 0; j < reducedFeatureSet.Rows(); j++)
        for (unsigned int k = 0; k < reducedFeatureSet.Cols(); k++)
          reducedFeatureSet(j, k) = std::get<2>(FoldingDataMap[index])(j, indices[k]);

      //check performance using cross-validation
      VectorDouble performance = InternalCrossValidation(reducedFeatureSet, std::get<1>(FoldingDataMap[index]), 1, 0.01,classifiertype);
      PerFeaturePerformance.push_back(performance[3]);
      std::cout << "feature #: " << featureNo << "performance=" <<performance[3]<< std::endl;
      // cbica::Logging(loggerFile, "feature #: " + std::to_string(featureNo) + " performance=" + std::to_string(performance[3]) + "\n");

      //check whether to contiue the loop or break due to no increase in perfrmance
      if (featureNo >= 10)
      {
        if (CheckPerformanceStatus(PerFeaturePerformance[featureNo], PerFeaturePerformance[featureNo - 1], PerFeaturePerformance[featureNo - 2], PerFeaturePerformance[featureNo - 3], PerFeaturePerformance[featureNo - 4], PerFeaturePerformance[featureNo - 5], PerFeaturePerformance[featureNo - 6], PerFeaturePerformance[featureNo - 7], PerFeaturePerformance[featureNo - 8], PerFeaturePerformance[featureNo - 9]) == false)
        {
          numberOfSelectedFeatures = featureNo - 9;
          break;
        }
        else
          numberOfSelectedFeatures = featureNo - 9;
      }
    }
    VariableSizeMatrixType selectedFeatureSet;
    selectedFeatureSet.SetSize(std::get<1>(FoldingDataMap[index]).size(), numberOfSelectedFeatures);
    for (unsigned int j = 0; j < selectedFeatureSet.Rows(); j++)
      for (int k = 0; k < numberOfSelectedFeatures; k++)
        selectedFeatureSet(j, k) = std::get<2>(FoldingDataMap[index])(j, k);



    //selection of the optimal values of classifiers
    double bestC = 0;
    double bestG = 0;
    double bestCV = 0;
    if (classifiertype == 2)
    {
      for (double cValue = 1; cValue <= 10; cValue++)
      {
        for (double gValue = 0.01; gValue <= 5; gValue = gValue + 0.05)
        {
          VectorDouble result = InternalCrossValidation(selectedFeatureSet, std::get<1>(FoldingDataMap[index]), cValue, gValue,classifiertype);
          if (result[3] > bestCV)
          {
            bestC = cValue;
            bestG = gValue;
            bestCV = result[3];
          }
        }
      }
    }
    else
    {
      for (double cValue = 1; cValue <= 10; cValue++)
      {
        VectorDouble result = InternalCrossValidation(selectedFeatureSet, std::get<1>(FoldingDataMap[index]), cValue, 0.01,classifiertype);
        if (result[3] > bestCV)
        {
          bestC = cValue;
          bestCV = result[3];
        }
      }
    }

    
    
    //train a new model on selcted features and optimal values of classifiers' parameters
    //copy training and test data from given dataset
    cv::Mat trainingData = cv::Mat::zeros(selectedFeatureSet.Rows(), selectedFeatureSet.Cols(), CV_32FC1);
    cv::Mat trainingLabels = cv::Mat::zeros(std::get<1>(FoldingDataMap[index]).size(), 1, CV_32FC1);
    cv::Mat testingData = cv::Mat::zeros(std::get<5>(FoldingDataMap[index]).Rows(), selectedFeatureSet.Cols(), CV_32FC1);
    cv::Mat testingLabels = cv::Mat::zeros(std::get<4>(FoldingDataMap[index]).size(), 1, CV_32FC1);

    int trainingCounter = 0;
    int testingCounter = 0;
    for (int copyDataCounter = 0; copyDataCounter < trainingData.rows; copyDataCounter++)
    {
      trainingLabels.ptr< float >(copyDataCounter)[0] = std::get<1>(FoldingDataMap[index])[copyDataCounter];
      for (int copyDataCounter2 = 0; copyDataCounter2 <trainingData.cols; copyDataCounter2++)
        trainingData.ptr< float >(copyDataCounter)[copyDataCounter2] = selectedFeatureSet(copyDataCounter, copyDataCounter2);
    }
    for (int copyDataCounter = 0; copyDataCounter < testingData.rows; copyDataCounter++)
    {
      testingLabels.ptr< float >(copyDataCounter)[0] = std::get<4>(FoldingDataMap[index])[copyDataCounter];
      for (int copyDataCounter2 = 0; copyDataCounter2 <testingData.cols; copyDataCounter2++)
        testingData.ptr< float >(copyDataCounter)[copyDataCounter2] = std::get<5>(FoldingDataMap[index])(copyDataCounter, indices[copyDataCounter2]);
    }
   if(classifiertype==2)
      std::cout << "index=" << index << ", Best C = " << bestC << ", Best G=" << bestG << ", Best CV=" << bestCV << ", Features=" << numberOfSelectedFeatures << std::endl;
      // cbica::Logging(loggerFile, "index=" + std::to_string(index) + ", Best C = " + std::to_string(bestC) + ", Best G = " + std::to_string(bestG) + ", Best CV = " + std::to_string(bestCV) + ", Features=" + std::to_string(numberOfSelectedFeatures) + "\n" );
   else
    std::cout << "index=" << index << ", Best C = " << bestC << ", Best CV=" << bestCV << ", Features=" << numberOfSelectedFeatures << std::endl;
    // cbica::Logging(loggerFile, "index=" + std::to_string(index) + ", Best C = " + std::to_string(bestC) + ", Features=" + std::to_string(numberOfSelectedFeatures) + "\n" );

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
    try
    {
      res = svm->train(trainingData, cv::ml::ROW_SAMPLE, trainingLabels);
      svm->save(outputfolder + "/Model_FoldNo_" +std::to_string(index+1)+".xml");
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
    }
    std::tuple<VectorDouble, VectorDouble, VariableSizeMatrixType, VectorDouble, VectorDouble, VariableSizeMatrixType, VectorDouble> new_tuple_duplicate = FoldingDataMap[index];
    std::get<6>(new_tuple_duplicate) = predictedLabels;
    FoldingDataMap[index] = new_tuple_duplicate;
  }

  //combining the predicted results from all the map entries
  VectorDouble FinalTargetLabels;
  VectorDouble FinalPredictedLabels;
  for (auto const &mapiterator : FoldingDataMap)
  {
    VectorDouble TargetLabels = std::get<4>(mapiterator.second);
    VectorDouble PredictedLabels = std::get<6>(mapiterator.second);
    for (int index = 0; index < TargetLabels.size(); index++)
    {
      FinalTargetLabels.push_back(TargetLabels[index]);
      FinalPredictedLabels.push_back(PredictedLabels[index]);
    }
  }

  //calcualte final performance and write in the csv file
  VectorDouble FinalPerformance = CalculatePerformanceMeasures(FinalPredictedLabels, FinalTargetLabels);
  std::ofstream myfile;
  myfile.open(outputfolder + "/predicted_labels.csv");
  for (unsigned int index1 = 0; index1 < FinalPredictedLabels.size(); index1++)
    myfile << std::to_string(FinalPredictedLabels[index1]) + "\n";
  myfile.close();

  myfile.open(outputfolder + "/performance.csv");
  myfile << "Accuracy,Sensitivity,Specificity,BalancedAccuracy \n";
  myfile << std::to_string(FinalPerformance[0])+","+ std::to_string(FinalPerformance[1]) + ","+std::to_string(FinalPerformance[2]) + "," + std::to_string(FinalPerformance[3]) +"\n";
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
  for (int index = 0; index < target.size(); index++)
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
    EffectSize.push_back((mean_set1[featureNo] - mean_set2[featureNo]) / sqrt(SP));
  }
  //std::vector<size_t> indices = sort_indexes(EffectSize);
  //VariableSizeMatrixType selected_feature_set;

  //for (int index1 = 0; index1 < training_features.Rows(); index1++)
  //	for (int index = 0; index < no_of_features; index++)
  //		selected_feature_set(index1, index) = training_features(index1, indices[index]);

  //return selected_feature_set;
  //EffectSize(find(isnan(EffectSize))) = 0.0001;
  return EffectSize;
}




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




bool TrainingModule::Run(const std::string inputFeaturesFile, const std::string inputLabelsFile, const std::string outputdirectory,const int classifiertype, const int foldtype)
{
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
    FeaturesOfAllSubjects.SetSize(dataMatrix.rows()-1, dataMatrix.columns()-1);

    for (unsigned int i = 1; i < dataMatrix.rows(); i++)
      for (unsigned int j = 1; j < dataMatrix.cols(); j++)
        FeaturesOfAllSubjects(i-1, j-1) = dataMatrix(i, j);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Cannot find the file 'features.csv' in the input directory. Error code : " + std::string(e1.what()));
    return false;
  }

  try
  {
    CSVFileReaderType::Pointer readerMean = CSVFileReaderType::New();
    readerMean->SetFileName(inputLabelsFile);
    readerMean->SetFieldDelimiterCharacter(',');
    readerMean->HasColumnHeadersOff();
    readerMean->HasRowHeadersOff();
    readerMean->Parse();
    dataMatrix = readerMean->GetArray2DDataObject()->GetMatrix();
    LabelsOfAllSubjects.SetSize(dataMatrix.rows()-1);

    for (unsigned int i = 1; i < dataMatrix.rows(); i++)
      LabelsOfAllSubjects[i-1] = dataMatrix(i, 0);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Cannot find the file 'features.csv' in the input directory. Error code : " + std::string(e1.what()));
    return false;
  }


  //scaling of input features and saving corresponding mean and standard deviation in the output directory
  FeatureScalingClass mFeaturesScaling;
  VariableSizeMatrixType scaledFeatureSet;
  VariableLengthVectorType meanVector;
  VariableLengthVectorType stdVector;
  mFeaturesScaling.ScaleGivenTrainingFeatures(FeaturesOfAllSubjects, scaledFeatureSet, meanVector, stdVector);


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

  TrainingModule mTrainingSimulator;
  VectorDouble FinalResult = mTrainingSimulator.CrossValidation(scaledFeatureSet, LabelsOfAllSubjects, outputdirectory,classifiertype,foldtype);
  //VectorDouble FinalResult = mTrainingSimulator.SplitTrainTest(scaledFeatureSet, LabelsOfAllSubjects, outputdirectory, classifiertype, foldtype,40);

   std::cout << "Accuray=" << FinalResult[0] << std::endl;
   std::cout << "Sensitivity=" << FinalResult[1] << std::endl;
   std::cout << "Specificity=" << FinalResult[2] << std::endl;
   std::cout << "Balanced Accuracy=" << FinalResult[3] << std::endl;

  //cbica::Logging(loggerFile, "Accuracy=" + std::to_string(FinalResult[0]) + "\n");
  //cbica::Logging(loggerFile, "Sensitivity=" + std::to_string(FinalResult[1]) + "\n");
  //cbica::Logging(loggerFile, "Specificity=" + std::to_string(FinalResult[2]) + "\n");
  //cbica::Logging(loggerFile, "Balanced Accuracy=" + std::to_string(FinalResult[3]) + "\n");

  return true;
}

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








VectorDouble TrainingModule::SplitTrainTest(const VariableSizeMatrixType inputFeatures, const VariableLengthVectorType inputLabels, 
                                            const std::string outputfolder, const int classifiertype, const int number_of_folds,const int training_size)
{
  MapType FoldingDataMap;
  std::vector<double> trainingindices;
  std::vector<double> traininglabels;
  std::vector<double> testingindices;
  std::vector<double> testinglabels;

  std::vector<double> predictedlabels;

  //make loops for training and validation, CVO partition like structure
  for (int index = 0; index < training_size; index++)
    trainingindices.push_back(index);

   for (unsigned int index = training_size; index < inputLabels.Size(); index++)
     testingindices.push_back(index);

   VariableSizeMatrixType trainingfeatures;
   VariableSizeMatrixType testingfeatures;

    //find training and testing labels and features
    trainingfeatures.SetSize(trainingindices.size(), inputFeatures.Cols());
    for (int i = 0; i < trainingindices.size(); i++)
    {
      traininglabels.push_back(inputLabels[trainingindices[i]]);
      for (unsigned int j = 0; j < trainingfeatures.Cols(); j++)
        trainingfeatures(i, j) = inputFeatures(trainingindices[i], j);
    }
    testingfeatures.SetSize(testingindices.size(), inputFeatures.Cols());
    for (int i = 0; i < testingindices.size(); i++)      
    {
      testinglabels.push_back(inputLabels[testingindices[i]]);
      predictedlabels.push_back(-1);
      for (unsigned int j = 0; j < testingfeatures.Cols(); j++)
        testingfeatures(i, j) = inputFeatures(testingindices[i], j);
    }
    std::tuple<VectorDouble, VectorDouble, VariableSizeMatrixType, VectorDouble, VectorDouble, VariableSizeMatrixType, VectorDouble> new_tuple(trainingindices, traininglabels, trainingfeatures, testingindices, testinglabels, testingfeatures, predictedlabels);
    FoldingDataMap[0] = new_tuple;

    
  for (unsigned int index = 0; index < FoldingDataMap.size(); index++)
  {
    std::cout << "***************************** index=" << index << "*********************************" << std::endl << std::endl;

    //feature selection mechanism
    //---------------------------
    //double numberOfSelectedFeatures = 44;
    //VectorDouble EffectSize = EffectSizeFeatureSelection(std::get<2>(FoldingDataMap[index]), std::get<1>(FoldingDataMap[index]));
    //for (int eSizeCounter = 0; eSizeCounter < EffectSize.size(); eSizeCounter++)
    //{
    //  if (EffectSize[eSizeCounter] < 0)
    //    EffectSize[eSizeCounter] = EffectSize[eSizeCounter] * -1;
    //}
    //std::vector<size_t> indices = sort_indexes(EffectSize);

    std::ofstream myfile;
    myfile.open(outputfolder + "/test_data.csv");
    for (unsigned int index1 = 0; index1 < std::get<2>(FoldingDataMap[index]).Rows(); index1++)
    {
      for (unsigned int index2 = 0; index2 < std::get<2>(FoldingDataMap[index]).Cols(); index2++)
      {
        if (index2 == 0)
          myfile << std::to_string(std::get<2>(FoldingDataMap[index])[index1][index2]);
        else
          myfile << "," << std::to_string(std::get<2>(FoldingDataMap[index])[index1][index2]);
      }
      myfile << "\n";
    }
    myfile.close();


    ////keep on copy the  features
    std::vector<double> selectedfeatures;
    std::vector<double> selectedperformances;
    int counter = 1;
    bool improvement = true;

    while (improvement == true)
    {
      std::vector<double> indices;
      std::vector<double> performances;

      //iterate over the remaining features
      for (int featureNo = 93; featureNo < 96; featureNo++)
        //for (int featureNo = 0; featureNo < std::get<2>(FoldingDataMap[index]).Cols(); featureNo++)
      {

        if (std::find(selectedfeatures.begin(), selectedfeatures.end(), featureNo) != selectedfeatures.end() == true)
          continue;

          VariableSizeMatrixType reducedFeatureSet;
          reducedFeatureSet.SetSize(std::get<1>(FoldingDataMap[index]).size(), selectedfeatures.size() + 1);

          //copy the already selected features to the reduced feature set
          for (unsigned int j = 0; j < reducedFeatureSet.Rows(); j++)
            for (int k = 0; k < selectedfeatures.size(); k++)
              reducedFeatureSet(j, k) = std::get<2>(FoldingDataMap[index])(j, selectedfeatures[k]);

          //copy the new feature
          for (unsigned int j = 0; j < reducedFeatureSet.Rows(); j++)
            reducedFeatureSet(j, selectedfeatures.size()) = std::get<2>(FoldingDataMap[index])(j, featureNo);
          
          //check performance using cross-validation
          VectorDouble performance = InternalCrossValidationSplitTrainTest(reducedFeatureSet, std::get<1>(FoldingDataMap[index]), 1, 0.1, classifiertype,counter, outputfolder);
          counter++;
          performances.push_back(performance[0]);
          indices.push_back(featureNo);
      }
      //find maximum performance
      double maxPerf = performances[0];
      int maxLoc = -1;
      for (int featureNo = 0; featureNo < performances.size(); featureNo++)
      {
        if (performances[featureNo] >= maxPerf)
        {
          maxPerf = performances[featureNo];
          maxLoc = featureNo;
        }
      }
      //check for performance improvement
      //if (selectedfeatures.size() >= 6)
      //{
      //  double averageper = (maxPerf + selectedperformances[selectedperformances.size() - 1] + selectedperformances[selectedperformances.size() - 2] + selectedperformances[selectedperformances.size() - 3] + selectedperformances[selectedperformances.size() - 4])/5;
      //  if (averageper <= selectedperformances[selectedperformances.size() - 5])
      //  {
      //    improvement = false;
      //    break;
      //  }
      //}
      
      if (selectedfeatures.size() >= 44)
          improvement = false;

      //add if there is improvement
      if (improvement == true)
      {
        selectedfeatures.push_back(indices[maxLoc]);
        selectedperformances.push_back(maxPerf);
      }
      std::cout << selectedfeatures.size() << " : " << indices[maxLoc] <<" : "<< maxPerf << std::endl;
    }
    
      
      int numberOfSelectedFeatures = 10;

      //std::cout << "feature #: " << featureNo << "performance=" <<performance[3]<< std::endl;

      //check whether to contiue the loop or break due to no increase in perfrmance
    //  if (featureNo >= 10)
    //  {
    //    if (CheckPerformanceStatus(PerFeaturePerformance[featureNo], PerFeaturePerformance[featureNo - 1], PerFeaturePerformance[featureNo - 2], PerFeaturePerformance[featureNo - 3], PerFeaturePerformance[featureNo - 4], PerFeaturePerformance[featureNo - 5], PerFeaturePerformance[featureNo - 6], PerFeaturePerformance[featureNo - 7], PerFeaturePerformance[featureNo - 8], PerFeaturePerformance[featureNo - 9]) == false)
    //    {
    //      numberOfSelectedFeatures = featureNo - 9;
    //      break;
    //    }
    //    else
    //      numberOfSelectedFeatures = featureNo - 9;
    //  }
    //}
    VariableSizeMatrixType selectedFeatureSet;
    selectedFeatureSet.SetSize(std::get<1>(FoldingDataMap[index]).size(), selectedfeatures.size());
    for (unsigned int j = 0; j < selectedFeatureSet.Rows(); j++)
      for (unsigned int k = 0; k < selectedFeatureSet.Cols(); k++)
        selectedFeatureSet(j, k) = std::get<2>(FoldingDataMap[index])(j, selectedfeatures[k]);



    //selection of the optimal values of classifiers
    double bestC = 0;
    double bestG = 0;
    double bestCV = 0;
    if (classifiertype == 2)
    {
      for (double cValue = 1; cValue <= 5; cValue++)
      {
        for (double gValue = 0.01; gValue <= 5; gValue = gValue + 0.05)
        {
          VectorDouble result = InternalCrossValidation(selectedFeatureSet, std::get<1>(FoldingDataMap[index]), cValue, gValue, classifiertype);
          if (result[3] > bestCV)
          {
            bestC = cValue;
            bestG = gValue;
            bestCV = result[3];
          }
        }
      }
    }
    else
    {
      for (double cValue = 1; cValue <= 5; cValue++)
      {
        VectorDouble result = InternalCrossValidation(selectedFeatureSet, std::get<1>(FoldingDataMap[index]), cValue, 0.1, classifiertype);
        if (result[0] > bestCV)
        {
          bestC = pow(2, cValue);
          bestCV = result[0];
        }
      }
    }



    //train a new model on selcted features and optimal values of classifiers' parameters
    //copy training and test data from given dataset
    cv::Mat trainingData = cv::Mat::zeros(selectedFeatureSet.Rows(), selectedFeatureSet.Cols(), CV_32FC1);
    cv::Mat trainingLabels = cv::Mat::zeros(std::get<1>(FoldingDataMap[index]).size(), 1, CV_32FC1);
    cv::Mat testingData = cv::Mat::zeros(std::get<5>(FoldingDataMap[index]).Rows(), selectedFeatureSet.Cols(), CV_32FC1);
    cv::Mat testingLabels = cv::Mat::zeros(std::get<4>(FoldingDataMap[index]).size(), 1, CV_32FC1);

    int trainingCounter = 0;
    int testingCounter = 0;
    for (int copyDataCounter = 0; copyDataCounter < trainingData.rows; copyDataCounter++)
    {
      trainingLabels.ptr< float >(copyDataCounter)[0] = std::get<1>(FoldingDataMap[index])[copyDataCounter];
      for (int copyDataCounter2 = 0; copyDataCounter2 <trainingData.cols; copyDataCounter2++)
        trainingData.ptr< float >(copyDataCounter)[copyDataCounter2] = selectedFeatureSet(copyDataCounter, copyDataCounter2);
    }
    //for (int copyDataCounter = 0; copyDataCounter < testingData.rows; copyDataCounter++)
    //{
    //  testingLabels.ptr< float >(copyDataCounter)[0] = std::get<4>(FoldingDataMap[index])[copyDataCounter];
    //  for (int copyDataCounter2 = 0; copyDataCounter2 <testingData.cols; copyDataCounter2++)
    //    testingData.ptr< float >(copyDataCounter)[copyDataCounter2] = std::get<5>(FoldingDataMap[index])(copyDataCounter, indices[copyDataCounter2]);
    //}
    //if (classifiertype == 2)
    //  std::cout << "index=" << index << ", Best C = " << bestC << ", Best G=" << bestG << ", Best CV=" << bestCV << ", Features=" << numberOfSelectedFeatures << std::endl;
    //else
    //  std::cout << "index=" << index << ", Best C = " << bestC << ", Best CV=" << bestCV << ", Features=" << numberOfSelectedFeatures << std::endl;

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
    try
    {
      res = svm->train(trainingData, cv::ml::ROW_SAMPLE, trainingLabels);
      svm->save(outputfolder + "/Model_FoldNo_" + std::to_string(index + 1) + ".xml");
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
    }
    std::tuple<VectorDouble, VectorDouble, VariableSizeMatrixType, VectorDouble, VectorDouble, VariableSizeMatrixType, VectorDouble> new_tuple_duplicate = FoldingDataMap[index];
    std::get<6>(new_tuple_duplicate) = predictedLabels;
    FoldingDataMap[index] = new_tuple_duplicate;
  }

  //combining the predicted results from all the map entries
  VectorDouble FinalTargetLabels;
  VectorDouble FinalPredictedLabels;
  for (auto const &mapiterator : FoldingDataMap)
  {
    VectorDouble TargetLabels = std::get<4>(mapiterator.second);
    VectorDouble PredictedLabels = std::get<6>(mapiterator.second);
    for (int index = 0; index < TargetLabels.size(); index++)
    {
      FinalTargetLabels.push_back(TargetLabels[index]);
      FinalPredictedLabels.push_back(PredictedLabels[index]);
    }
  }

  //calculate final performance and write in the csv file
  VectorDouble FinalPerformance = CalculatePerformanceMeasures(FinalPredictedLabels, FinalTargetLabels);
  std::ofstream myfile;
  myfile.open(outputfolder + "/predicted_labels.csv");
  for (unsigned int index1 = 0; index1 < FinalPredictedLabels.size(); index1++)
    myfile << std::to_string(FinalPredictedLabels[index1]) + "\n";
  myfile.close();

  myfile.open(outputfolder + "/performance.csv");
  myfile << "Accuracy,Sensitivity,Specificity,BalancedAccuracy \n";
  myfile << std::to_string(FinalPerformance[0]) + "," + std::to_string(FinalPerformance[1]) + "," + std::to_string(FinalPerformance[2]) + "," + std::to_string(FinalPerformance[3]) + "\n";
  myfile.close();


  return FinalPerformance;
}
