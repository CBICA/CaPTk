/**
\file  SVMTrain.h

\brief Declaration of SVMTrain class

https://www.cbica.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2016 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.cbica.upenn.edu/sbia/software/license.html

*/

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include "svm.h"
//#include "CAPTk.h"

class SVMTrain
{
public:
  struct svm_parameter param;		// set by parse_command_line
  struct svm_problem prob;		// set by read_problem
  struct svm_model *model;
  struct svm_model *best_model;
  struct svm_node *x_space;
  int mCrossValidation;
  int mNumberOfFolds;
  char *line;
  int max_line_len;
  int bestc;
  int bestg;
  int bestaccuracy;



  int GetBestGValue()
  {
    return bestg;
  }
  int GetBestCValue()
  {
    return bestc;
  }
  int GetBestAccuray()
  {
    return bestaccuracy;
  }


  SVMTrain();
  double do_cross_validation();
  int main_train(VariableSizeMatrixType &trainingdata);
  void parse_command_line();
  void read_problem_from_matrix(VariableSizeMatrixType &trainingdata);
};
