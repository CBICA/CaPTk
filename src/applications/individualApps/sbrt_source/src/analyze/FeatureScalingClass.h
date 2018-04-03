/**
\file  FeatureScalingClass.h

\brief Declaration of the FeatureScalingClass

https://www.cbica.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2016 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.cbica.upenn.edu/sbia/software/license.html

*/

#pragma once

#include "iostream"
//#include "CAPTk.h"
#include <itkVariableLengthVector.h>
#include <itkVector.h>
#include <itkVariableSizeMatrix.h>


using namespace std;
typedef itk::VariableSizeMatrix< double > VariableSizeMatrixType;
typedef itk::VariableLengthVector< double > VariableLengthVectorType;

class  FeatureScalingClass
{
public:
  //!Constructor
  FeatureScalingClass(int mnum_features);
  
  //!Destructor
  ~FeatureScalingClass();

  /**
  \brief Scales given training features
  \param inputdata Training data
  */
  VariableSizeMatrixType ScaleGivenTrainingFeatures(const VariableSizeMatrixType &inputdata);


  /**
  \brief Scales given test features
  \param inputdata Test data
  */
  VariableSizeMatrixType ScaleGivenTestingFeatures(const VariableSizeMatrixType &inputdata);

  //!Returns mean of current feature set
  VariableLengthVectorType GetMeanVector()
  {
    return mMeanVector;
  }
  
  //!Returns standard deviation of current feature set
  VariableLengthVectorType GetStdVector()
  {
    return mStdVector;
  }

  /**
  \brief Sets the mean and standard deviation vectors
  \param meanvector Vector of means of given features
  \param stdvector Vector of standard deviation of given features
  */
  void SetParameters(VariableLengthVectorType &meanvector, VariableLengthVectorType &stdvector)
  {
    mMeanVector = meanvector;
    mStdVector = stdvector;
  }

  /**
  \brief Returns the z-score of a feature value
  \param mean Mean of the feature vector
  \param variance Variance of the feature vector
  \param featureval Value of the feature
  */
  double GetZScore(const double &mean, const double &variance, const double &featureval);

  /**
  \brief Returns mean of a feature specified by given index
  \param index Index of the feature in the feature vector
  */
  double GetMeanForGivenFeatureIndex(const int &index)
  {
    return mMeanVector[index];
  }

  /**
  \brief Returns standard deviation of a feature specified by given index
  \param index Index of the feature in the feature vector
  */
  double GetStdForGivenFeatureIndex(const int &index)
  {
    return mStdVector[index];
  }

  //!Resets the parameters of the class
  void ResetParameters()
  {
    mMeanVector.SetSize(0);
    mStdVector.SetSize(0);
  }

private:
  VariableLengthVectorType mMeanVector;
  VariableLengthVectorType mStdVector;
  int mnum_features;
};
