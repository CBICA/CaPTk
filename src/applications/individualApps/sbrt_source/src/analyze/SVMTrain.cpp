/**
\file  SVMTrain.cpp

\brief Implementation of SVMTrain class

https://www.cbica.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2016 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.cbica.upenn.edu/sbia/software/license.html

*/

#include "itkCSVArray2DFileReader.h"
#include "SVMTrain.h"

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))


SVMTrain::SVMTrain()
{
  mNumberOfFolds = 10;
  param.svm_type = C_SVC;
  param.kernel_type = RBF;
  param.degree = 3;
  //param.gamma = 0;
  param.coef0 = 0;
  param.nu = 0.5;
  param.cache_size = 100;
  //param.C = 1;
  param.eps = 1e-3;
  param.p = 0.1;
  param.shrinking = 1;
  param.probability = 1;
  param.nr_weight = 0;
  param.weight_label = NULL;
  param.weight = NULL;
  mCrossValidation = 1;
  bestc = 0;
  bestg = 0;
  bestaccuracy = 0;
}


int SVMTrain::main_train(VariableSizeMatrixType &trainingdata)
{
  std::string model_file_name;
  FILE* o;
#ifdef _WIN32
  fopen_s(&o, "iterations.txt", "w");
#else
  o = fopen("iterations.txt", "w");
#endif
  
  for (double c = 10; c <= 100; c = c + 10)
  {
    //for (double g= -5; g<=5; g=g+2)
    //{
    param.gamma = 0.0769231;
    //param.gamma = std::pow(2, g);
    param.C = c;// std::pow(2, c);

    const char *error_msg;
    // double overall_accuracy = 0;
    read_problem_from_matrix(trainingdata);
    error_msg = svm_check_parameter(&prob, &param);
    if (error_msg)
    {
      fprintf(stderr, "ERROR: %s\n", error_msg);
      exit(1);
    }
    if (mCrossValidation)
    {
      double accuracy = do_cross_validation();
      model = svm_train(&prob, &param);
      std::string filename = "ModelFile_" + std::to_string(c) + "_" + std::to_string(param.gamma) + ".model";
      model_file_name = filename;      
      if (svm_save_model(model_file_name.c_str(), model))
      {
        std::cerr << "Can't save model to file '" << model_file_name << "'\n";
        exit(EXIT_FAILURE);
      }
      fprintf(o, "%f %f %f \n", param.C, param.gamma, accuracy);
      if (accuracy > bestaccuracy)
      {
        bestaccuracy = accuracy;
        bestg = param.gamma;
        bestc = c;
      }
    }
    else
    {
      model = svm_train(&prob, &param);
    }
    svm_destroy_param(&param);
    free(prob.y);
    free(prob.x);
    free(x_space);
    svm_free_and_destroy_model(&model);
    //}
  }
  if (mCrossValidation)
  {
    param.gamma = std::pow(2, bestg);
    param.C = std::pow(2, bestc);

    const char *error_msg;
    // double overall_accuracy = 0;
    read_problem_from_matrix(trainingdata);
    error_msg = svm_check_parameter(&prob, &param);
    if (error_msg)
    {
      fprintf(stderr, "ERROR: %s\n", error_msg);
      exit(1);
    }

    model = svm_train(&prob, &param);
    std::string filename = "FinalModelFile_" + std::to_string(bestc) + "_" + std::to_string(bestg) + ".model";
    model_file_name = filename;
    if (svm_save_model(model_file_name.c_str(), model))
    {
      std::cerr << "Can't save model to file '" << model_file_name << "'\n";
      exit(EXIT_FAILURE);
    }
    svm_destroy_param(&param);
    free(prob.y);
    free(prob.x);
    free(x_space);
    svm_free_and_destroy_model(&model);
  }
  fclose(o);
  return 1;
}


double SVMTrain::do_cross_validation()
{
  int i;
  int total_correct = 0;
  double total_error = 0;
  double sumv = 0, sumy = 0, sumvv = 0, sumyy = 0, sumvy = 0;
  double *target = Malloc(double, prob.l);
  double accuracy = 0;

  svm_cross_validation(&prob, &param, mNumberOfFolds, target);
  if (param.svm_type == EPSILON_SVR ||
    param.svm_type == NU_SVR)
  {
    for (i = 0; i<prob.l; i++)
    {
      double y = prob.y[i];
      double v = target[i];
      total_error += (v - y)*(v - y);
      sumv += v;
      sumy += y;
      sumvv += v*v;
      sumyy += y*y;
      sumvy += v*y;
    }
    printf("Cross Validation Mean squared error = %g\n", total_error / prob.l);
    printf("Cross Validation Squared correlation coefficient = %g\n",
      ((prob.l*sumvy - sumv*sumy)*(prob.l*sumvy - sumv*sumy)) /
      ((prob.l*sumvv - sumv*sumv)*(prob.l*sumyy - sumy*sumy))
      );
  }
  else
  {
    for (i = 0; i<prob.l; i++)
      if (target[i] == prob.y[i])
        ++total_correct;
    printf("Cross Validation Accuracy = %g%%\n", 100.0*total_correct / prob.l);
    accuracy = 100.0*total_correct / prob.l;
  }
  free(target);
  return accuracy;
}


void SVMTrain::read_problem_from_matrix(VariableSizeMatrixType &trainingdata)
{
  size_t elements;
  prob.l = trainingdata.Rows();
  elements = trainingdata.Rows() * trainingdata.Cols();

  prob.y = Malloc(double, prob.l);
  prob.x = Malloc(struct svm_node *, prob.l);
  x_space = Malloc(struct svm_node, elements);

  int j = 0;
  for (int i = 0; i < prob.l; i++)
  {
    prob.x[i] = &x_space[j];
    prob.y[i] = trainingdata(i, trainingdata.Cols() - 1);
    unsigned int inputEntries = 0;
    while (inputEntries < trainingdata.Cols() - 1)
    {
      x_space[j].index = inputEntries + 1;
      x_space[j].value = trainingdata(i, inputEntries);
      inputEntries++;
      ++j;
    }
    x_space[j++].index = -1;

  }
  //param.gamma = 1.0 / trainingdata.Cols();
}
