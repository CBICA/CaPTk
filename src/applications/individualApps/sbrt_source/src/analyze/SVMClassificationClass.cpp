/**
\file  SVMClassificationClass.cpp

\brief Implementation of SVMClassificationClass

https://www.cbica.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2016 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.cbica.upenn.edu/sbia/software/license.html

*/

#include <algorithm>
#include <stdio.h>

#include "SVMClassificationClass.h"
#include "svm.h"

#if (_WIN32)
static const char  cSeparator = '\\';
//  static const char* cSeparators = "\\/";
#else
static const char  cSeparator = '/';
//  static const char* cSeparators = "/";
#endif

struct svm_parameter param;       // set by parse_command_line
struct svm_problem prob;          // set by read_problem
struct svm_model *model;
struct svm_model *best_model;
struct svm_node *x_space;
int mCrossValidation;
int mNumberOfFolds = 10;
int bestc;
int bestg;
int bestaccuracy;
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
// read in a problem (in svmlight format)

static char *line = NULL;
static int max_line_len;

static char* readline(FILE *input)
{
  int len;

  if (fgets(line, max_line_len, input) == NULL)
    return NULL;

  while (strrchr(line, '\n') == NULL)
  {
    max_line_len *= 2;
    line = (char *)realloc(line, max_line_len);
    len = (int)strlen(line);
    if (fgets(line + len, max_line_len - len, input) == NULL)
      break;
  }
  return line;
}

void exit_input_error(int line_num)
{
  //fprintf(stderr,"Wrong input format at line %d\n", line_num);
  exit(1);
}

void read_problem_txt(const char *filename)
{
  int max_index, inst_max_index, i;
  size_t elements, j;
  FILE *fp;
#ifdef _WIN32
  fopen_s(&fp, filename, "r");
#else
  fp = fopen(filename, "r");
#endif

  char *endptr;
  char *idx, *val, *label;

  if (fp == NULL)
  {
    fprintf(stderr, "can't open input file %s\n", filename);
    exit(1);
  }

  prob.l = 0;
  elements = 0;

  max_line_len = 1024;
  line = Malloc(char, max_line_len);
  while (readline(fp) != NULL)
  {
    char* p = NULL; 
#ifdef _WIN32
    p = strtok_s(line, " \t", NULL);
#else
    p = strtok(line, " \t"); // label
#endif

    // features
    while (1)
    {
#ifdef _WIN32
      p = strtok_s(NULL, " \t", NULL);
#else
      p = strtok(NULL, " \t"); // label
#endif

      if (p == NULL || *p == '\n') // check '\n' as ' ' may be after the last feature
        break;
      ++elements;
    }
    ++elements;
    ++prob.l;
  }
  rewind(fp);

  prob.y = Malloc(double, prob.l);
  prob.x = Malloc(struct svm_node *, prob.l);
  x_space = Malloc(struct svm_node, elements);

  max_index = 0;
  j = 0;
  for (i = 0; i<prob.l; i++)
  {
    inst_max_index = -1; // strtol gives 0 if wrong format, and precomputed kernel has <index> start from 0
    readline(fp);
    prob.x[i] = &x_space[j];
#ifdef _WIN32
    label = strtok_s(line, " \t\n", NULL);
#else
    label = strtok(line, " \t\n");
#endif

    if (label == NULL) // empty line
      exit_input_error(i + 1);

    prob.y[i] = strtod(label, &endptr);
    if (endptr == label || *endptr != '\0')
      exit_input_error(i + 1);

    while (1)
    {
#ifdef _WIN32
      idx = strtok_s(NULL, ":", NULL);
#else
      idx = strtok(NULL, ":");
#endif

#ifdef _WIN32
      val = strtok_s(NULL, " \t", NULL);
#else
      val = strtok(NULL, " \t");
#endif


      if (val == NULL)
        break;

      errno = 0;
      x_space[j].index = (int)strtol(idx, &endptr, 10);
      if (endptr == idx || errno != 0 || *endptr != '\0' || x_space[j].index <= inst_max_index)
        exit_input_error(i + 1);
      else
        inst_max_index = x_space[j].index;

      errno = 0;
      x_space[j].value = strtod(val, &endptr);
      if (endptr == val || errno != 0 || (*endptr != '\0' && !isspace(*endptr)))
        exit_input_error(i + 1);

      ++j;
    }

    if (inst_max_index > max_index)
      max_index = inst_max_index;
    x_space[j++].index = -1;
  }

  if (param.gamma == 0 && max_index > 0)
    param.gamma = 1.0 / max_index;

  if (param.kernel_type == PRECOMPUTED)
    for (i = 0; i<prob.l; i++)
    {
      if (prob.x[i][0].index != 0)
      {
        fprintf(stderr, "Wrong input format: first column must be 0:sample_serial_number\n");
        exit(1);
      }
      if ((int)prob.x[i][0].value <= 0 || (int)prob.x[i][0].value > max_index)
      {
        fprintf(stderr, "Wrong input format: sample_serial_number out of range\n");
        exit(1);
      }
    }

  fclose(fp);
}

double do_cross_validation()
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

void read_problem_from_matrix(VariableSizeMatrixType &trainingdata)
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
      std::cout << x_space[j].value;
      inputEntries++;
      ++j;
    }
    x_space[j++].index = -1;

  }
  //param.gamma = 1.0 / trainingdata.Cols();
}

double main_train(VariableSizeMatrixType &trainingdata, double cVal, double gVal, const std::string &outputDirectory)
{
  param.svm_type = C_SVC;
  param.kernel_type = LINEAR;
  param.degree = 3;
  param.coef0 = 0;
  param.nu = 0.5;
  param.cache_size = 100;
  param.eps = 1e-3;
  param.p = 0.1;
  param.shrinking = 1;
  param.probability = 1;
  param.nr_weight = 0;
  param.weight_label = NULL;
  param.weight = NULL;
  param.C = cVal;
  param.gamma = gVal;

  //FILE* o;
  //o = fopen("iterations.txt", "w");
  // int mCrossValidation = 1;
  // double bestc = 0;
  // double bestg = 0;
  // double bestaccuracy = 0;
  int cross_validation = 1;
  double accuracy = 0;

  std::string input_file_name = "TrainingData.txt";
  std::string finalOutputDirectory = outputDirectory + std::string(&cSeparator) + "FinalModelFile.model";
  const char * model_file_name = finalOutputDirectory.c_str();
  const char *error_msg;
  read_problem_from_matrix(trainingdata);
  error_msg = svm_check_parameter(&prob, &param);

  if (error_msg)
  {
    fprintf(stderr, "ERROR: %s\n", error_msg);
    exit(1);
  }

  if (cross_validation)
  {
    accuracy = do_cross_validation();
    //=----------additional lines to remove---------------
    model = svm_train(&prob, &param);
    if (svm_save_model(model_file_name, model))
    {
      fprintf(stderr, "can't save model to file %s\n", model_file_name);
      exit(1);
    }
    svm_free_and_destroy_model(&model);
    //----------------------------------------------------
  }
  else
  {
    model = svm_train(&prob, &param);
    if (svm_save_model(model_file_name, model))
    {
      fprintf(stderr, "can't save model to file %s\n", model_file_name);
      exit(1);
    }
    svm_free_and_destroy_model(&model);
  }
  return accuracy;


  //for (double c = 10; c <= 100; c = c + 10)
  //{
  //     //for (double g= -5; g<=5; g=g+2)
  //     //{
  //     param.svm_type = C_SVC;
  //     param.kernel_type = RBF;
  //     param.degree = 3;
  //     param.coef0 = 0;
  //     param.nu = 0.5;
  //     param.cache_size = 100;
  //     param.eps = 1e-3;
  //     param.p = 0.1;
  //     param.shrinking = 1;
  //     param.probability = 1;
  //     param.nr_weight = 0;
  //     param.weight_label = NULL;
  //     param.weight = NULL;

  //     param.gamma = 0.0769231;
  //     //param.gamma = std::pow(2, g);
  //     param.C = c;// std::pow(2, c);

  //     const char *error_msg;
  //     double overall_accuracy = 0;
  //     read_problem_from_matrix(trainingdata);
  //     error_msg = svm_check_parameter(&prob, &param);
  //     if (error_msg)
  //     {
  //            fprintf(stderr, "ERROR: %s\n", error_msg);
  //            exit(1);
  //     }
  //     if (mCrossValidation)
  //     {
  //            double accuracy = do_cross_validation();
  //            model = svm_train(&prob, &param);
  //            std::string filename = "ModelFile_" + cbica::toString<int>(c) + "_" + cbica::toString<int>(param.gamma) + ".model";
  //            strcpy(model_file_name, filename.c_str());
  //            if (svm_save_model(model_file_name, model))
  //            {
  //                   fprintf(stderr, "can't save model to file %s\n", model_file_name);
  //                   exit(1);
  //            }
  //            fprintf(o, "%f %f %f \n", param.C, param.gamma, accuracy);
  //            if (accuracy > bestaccuracy)
  //            {
  //                   bestaccuracy = accuracy;
  //                   bestg = param.gamma;
  //                   bestc = c;
  //            }
  //     }
  //     else
  //     {
  //            model = svm_train(&prob, &param);
  //     }
  //     svm_destroy_param(&param);
  //     free(prob.y);
  //     free(prob.x);
  //     free(x_space);
  //     svm_free_and_destroy_model(&model);
  //     //}
  //}
  //if (mCrossValidation)
  //{
  //     param.gamma = std::pow(2, bestg);
  //     param.C = std::pow(2, bestc);

  //     const char *error_msg;
  //     double overall_accuracy = 0;
  //     read_problem_from_matrix(trainingdata);
  //     error_msg = svm_check_parameter(&prob, &param);
  //     if (error_msg)
  //     {
  //            fprintf(stderr, "ERROR: %s\n", error_msg);
  //            exit(1);
  //     }

  //     model = svm_train(&prob, &param);
  //     std::string filename = "FinalModelFile_" + cbica::toString<int>(bestc) + "_" + cbica::toString<int>(bestg) + ".model";
  //     strcpy(model_file_name, filename.c_str());
  //     if (svm_save_model(model_file_name, model))
  //     {
  //            fprintf(stderr, "can't save model to file %s\n", model_file_name);
  //            exit(1);
  //     }
  //     svm_destroy_param(&param);
  //     free(prob.y);
  //     free(prob.x);
  //     free(x_space);
  //     svm_free_and_destroy_model(&model);
  //}
  //fclose(o);
  //return 1;
}





SVMClassificationClass::SVMClassificationClass()
{
  mTrainingClassObject = new SVMTrain();
  mTestingClassObject = new SVMTest();
  mTrainedModel = NULL;
  mModelFile = "";
  mLastEncounteredError = "";
}
SVMClassificationClass::~SVMClassificationClass()
{
}

int SVMClassificationClass::Training(VariableSizeMatrixType &trainingdata, const std::string &outputDirectory)
{
  //double accuracy[10];
  //int index = 0;
  //for (int i = 30; i <=100; i = i + 10)
  //{
  //     accuracy[index] = main_train(trainingdata, i);
  //     index++;
  //}

  ////mTrainingClassObject->main_train(trainingdata);
  //return 1;

  ////if (!system(NULL))
  ////{
  ////   return;
  ////}
  FILE* o;
#ifdef _WIN32
  fopen_s(&o, "iterations.txt", "w");
#else
  o = fopen("iterations.txt", "w");
#endif

  double gBest = 0;
  double cBest = 0;
  double best_accuracy = 0;
  for (double c = -5; c <= 5; c = c + 2)
  {
    for (double g = -5; g <= 5; g = g + 2)
    {
      double gValue = std::pow(2, g);
      double cValue = std::pow(2, c);
      double accuracy = main_train(trainingdata, cValue, gValue, outputDirectory);
      ////mTrainingClassObject->main_train(trainingdata);
      //=======
      //     //double gBest = 0;
      //     //double cBest = 0;
      //     //double best_accuracy = 0;
      //     //for (double c =10; c <= 100; c = c + 10)
      //     //{
      //     ///*   for (double g = -5; g <= 5; g = g + 2)
      //     //     {*/
      //     //            double gValue = 0.0769231;// std::pow(2, g);
      //     //            double cValue = c;// std::pow(2, c);
      //>>>>>>> .r98
      //
      //<<<<<<< .mine
      //int accuracy = system(systemCall.c_str());
      //if (accuracy == 0)
      //{
      //     mLastEncounteredError = "Error in the training process.";
      //     return 0;
      //}
      //else
      if (accuracy >= best_accuracy)
      {
        best_accuracy = accuracy;
        gBest = gValue;
        cBest = cValue;
      }
      fprintf(o, "%f %f %f \n", cValue, gValue, accuracy);
      //=======
      //            //mTrainingClassObject->main_train(trainingdata);

      //            int accuracy = system(systemCall.c_str());
      //            if (accuracy == 0)
      //            {
      //                   mLastEncounteredError = "Error in the training process.";
      //                   return 0;
      //            }
      //            else if (accuracy >= best_accuracy)
      //            {
      //                   best_accuracy = accuracy;
      //                   gBest = gValue;
      //                   cBest = cValue;
      //            }
      //            fprintf(o, "%f %f %d \n", cValue, gValue, accuracy);
      //>>>>>>> .r98
      //
      //<<<<<<< .mine
    }
  }
  main_train(trainingdata, cBest, gBest, outputDirectory);
  //int status = 0;
  //status  = system(systemCall.c_str());
  //if (status == 0)
  //     mLastEncounteredError = "Error in the training process.";

  //=======
  //     //     }
  //     ////}
  //     //int status = 0;
  //     //status  = system(systemCall.c_str());
  //     //if (status == 0)
  //     //     mLastEncounteredError = "Error in the training process.";
  //     //
  //>>>>>>> .r98

  //fclose(o);
  return 1;
}

VectorVectorDouble SVMClassificationClass::Testing(VariableSizeMatrixType &testdata, bool existingclassification, std::string filename)
{
  VectorVectorDouble result;
  if (existingclassification)
    mTestingClassObject->SetModelFileName(mModelFile);
  result = mTestingClassObject->main_test(testdata, existingclassification, filename);

  return result;
}
