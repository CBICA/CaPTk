/**
\file  SVMTest.h

\brief Declaration of SVMTest class

https://www.cbica.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2016 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.cbica.upenn.edu/sbia/software/license.html

*/

//#include "stdafx.h"
#include "svm.h"
//#include "CAPTk.h"


class SVMTest
{
public:
  struct svm_node *x;
  struct svm_model* model;
  int max_nr_attr;
  int predict_probability;
  char *line;
  int max_line_len;
  std::string mModelFileName;


  SVMTest()
  {
    max_nr_attr = 200;
    predict_probability = 1;
  }

#if defined(__MSC_VER__)
  int(*info)(const char *fmt, ...) = &printf; // GCC handles c++11 tag in weird way
#endif

  void exit_input_error(int line_num);
  VectorVectorDouble predict(FILE *output, VariableSizeMatrixType &testdata);
  VectorVectorDouble main_test(VariableSizeMatrixType &testdata, bool existing, std::string &modelFileName);

  void SetModelFileName(std::string &filename)
  {
    mModelFileName = filename;
  }
};
