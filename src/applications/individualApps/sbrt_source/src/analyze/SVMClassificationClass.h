/**
\file  SVMClassificationClass.cpp

\brief Declaration of SVMClassificationClass

https://www.cbica.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2016 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.cbica.upenn.edu/sbia/software/license.html

*/

#pragma once

#include "stdio.h"
#include "SVMTrain.h"
#include "SVMTest.h"
//#include "CAPTk.h"
#include <vector>
#include <itkVariableLengthVector.h>
#include <itkVector.h>
#include <itkVariableSizeMatrix.h>
#include "svm.h"

//typedef std::vector< std::vector < double > > VectorVectorDouble;
//typedef itk::VariableSizeMatrix< double > VariableSizeMatrixType;

class SVMClassificationClass
{
public:
  SVMTrain * mTrainingClassObject;
  SVMTest *  mTestingClassObject;
  svm_model * mTrainedModel;
  std::string mModelFile;
  std::string mLastEncounteredError;

  std::string GetLastEncounteredError()
  {
    return mLastEncounteredError;
  }
  std::string GetModelFileName()
  {
    return mModelFile;
  }
  void SetModelFileName(std::string filename)
  {
    mModelFile = filename;
  }
  SVMClassificationClass();
  ~SVMClassificationClass();

  int Training(VariableSizeMatrixType &trainingdata, const std::string &outputDirectory);
  VectorVectorDouble Testing(VariableSizeMatrixType &testdata, bool classmethod, std::string modelFileName);
};
