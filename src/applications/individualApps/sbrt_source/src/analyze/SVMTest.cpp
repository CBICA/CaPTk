/**
\file  SVMTest.cpp

\brief Implementation of SVMTest class

https://www.cbica.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2016 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.cbica.upenn.edu/sbia/software/license.html

*/

#include "SVMTest.h"
#include "itkCSVArray2DFileReader.h"

VectorVectorDouble SVMTest::predict(FILE *output, VariableSizeMatrixType &testdata)
{
  int correct = 0;
  int total = 0;
  double error = 0;
  double sump = 0, sumt = 0, sumpp = 0, sumtt = 0, sumpt = 0;

  int svm_type = svm_get_svm_type(model);
  int nr_class = svm_get_nr_class(model);
  double *prob_estimates = NULL;
  int j;



  if (predict_probability)
  {
    if (svm_type == NU_SVR || svm_type == EPSILON_SVR)
#if defined (__MSC_VER__)
      info
#else
      printf
#endif
      ("Prob. model for test data: target value = predicted value + z,\nz: Laplace distribution e^(-|z|/sigma)/(2sigma),sigma=%g\n", svm_get_svr_probability(model));
    else
    {
      int *labels = (int *)malloc(nr_class*sizeof(int));
      svm_get_labels(model, labels);
      prob_estimates = (double *)malloc(nr_class*sizeof(double));
      fprintf(output, "labels");
      for (j = 0; j<nr_class; j++)
        fprintf(output, " %d", labels[j]);
      fprintf(output, "\n");
      free(labels);
    }
  }

  VectorVectorDouble results;
  // can this be made parallel?
  for (unsigned int NumberOfSamples = 0; NumberOfSamples < testdata.Rows(); NumberOfSamples++)
  {
    VectorDouble resultForOneSample;
    int i = 0;
    double target_label, predict_label;
    target_label = testdata(NumberOfSamples, testdata.Cols() - 1);

    unsigned int inputEntries = 0;
    while (inputEntries < testdata.Cols() - 1)
    {
      x[i].index = inputEntries + 1;
      x[i].value = testdata(NumberOfSamples, inputEntries);
      inputEntries++;
      ++i;
    }
    x[i].index = -1;

    if (predict_probability && (svm_type == C_SVC || svm_type == NU_SVC))
    {
      predict_label = svm_predict_probability(model, x, prob_estimates);
      fprintf(output, "%g", predict_label);
      for (j = 0; j<nr_class; j++)
        fprintf(output, " %g", prob_estimates[j]);

      resultForOneSample.push_back(predict_label);
      for (j = 0; j<nr_class; j++)
        resultForOneSample.push_back(prob_estimates[j]);

      results.push_back(resultForOneSample);
      fprintf(output, "\n");
    }
    else
    {
      predict_label = svm_predict(model, x);
      //fprintf(output, "%g\n", predict_label);
    }

    if (predict_label == target_label)
      ++correct;
    error += (predict_label - target_label)*(predict_label - target_label);
    sump += predict_label;
    sumt += target_label;
    sumpp += predict_label*predict_label;
    sumtt += target_label*target_label;
    sumpt += predict_label*target_label;
    ++total;
  }
  if (svm_type == NU_SVR || svm_type == EPSILON_SVR)
  {
#if defined (__MSC_VER__)
    info
#else
    printf
#endif
      ("Mean squared error = %g (regression)\n", error / total);
#if defined (__MSC_VER__)
    info
#else
    printf
#endif
      ("Squared correlation coefficient = %g (regression)\n",
      ((total*sumpt - sump*sumt)*(total*sumpt - sump*sumt)) /
      ((total*sumpp - sump*sump)*(total*sumtt - sumt*sumt))
      );
  }
  else
#if defined (__MSC_VER__)
    info
#else
    printf
#endif
    ("Accuracy = %g%% (%d/%d) (classification)\n",
    (double)correct / total * 100, correct, total);
  if (predict_probability)
    free(prob_estimates);

  return results;
}

VectorVectorDouble SVMTest::main_test(VariableSizeMatrixType &testdata, bool existing, std::string &filename)
{
  FILE /**input,*/ *outputfile;
  predict_probability = 1;

#ifdef _WIN32
  fopen_s(&outputfile, "output.txt", "w");
#else
  outputfile = fopen("output.txt", "w");
#endif

  if (outputfile == NULL)
  {
    fprintf(stderr, "can't open output file");
    exit(1);
  }
  if (existing == false)
    mModelFileName = filename;
  
  model = svm_load_model(filename.c_str());
  
  if (model == NULL)
  {
    fprintf(stderr, "can't open model file");
    exit(1);
  }

  x = (struct svm_node *) malloc(max_nr_attr*sizeof(struct svm_node));
  if (predict_probability)
  {
    if (svm_check_probability_model(model) == 0)
    {
      fprintf(stderr, "Model does not support probabiliy estimates\n");
      exit(1);
    }
  }
  else
  {
    if (svm_check_probability_model(model) != 0)
#if defined (__MSC_VER__)
      info
#else
      printf
#endif
      ("Model supports probability estimates, but disabled in prediction.\n");
  }

  VectorVectorDouble result = predict(outputfile, testdata);
  svm_free_and_destroy_model(&model);
 // free(x);
//  fclose(outputfile);
  return result;
}
