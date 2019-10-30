/**
\file  FeatureReductionClass.h

\brief Declaration of the FeatureReductionClass

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

#pragma once

#include "iostream"
#include "vtkSmartPointer.h"
#include "vnl/vnl_matrix.h"
//#include "vnl_matrix.h"
#include "itkVariableSizeMatrix.h"
#include "itkVariableLengthVector.h"
//#include "vtkTable.h"
//#include "CAPTk.h"
using VectorVectorDouble = std::vector< std::vector < double > >;
using VariableSizeMatrixType = itk::VariableSizeMatrix< double >;
using VariableLengthVectorType = itk::VariableLengthVector< double >;

class vtkTable;

#define NO_OF_PCA_FEATURES 45 // TODO: this selects the number of PCA features that are generated - this needs to be dynamic and picked up from the size[3] of the loaded perfusion image; search for this variable and "45" everywhere and replace with perfusionImage_size[3]

class  FeatureReductionClass
{
public:
  //!Constructor
  FeatureReductionClass();

  //!Destructor
  ~FeatureReductionClass();

  /**
  \brief Finds the first few discerning principal components in perfusion data
  \param intensities Input perfusion data
  */
  vtkSmartPointer< vtkTable >  GetDiscerningPerfusionTimePointsFullPCA(VectorVectorDouble &intensities, VariableSizeMatrixType &TransformationMatrix, VariableLengthVectorType &MeanVector);

  vtkSmartPointer<vtkTable> GetDiscerningPerfusionTimePoints(VectorVectorDouble &intensities);
  vtkSmartPointer<vtkTable> GetDiscerningPerfusionTimePointsDynamic(VectorVectorDouble &intensities);
  vtkSmartPointer< vtkTable >  GetDiscerningPerfusionTimePoints(vnl_matrix<double> &intensities);
  vtkSmartPointer<vtkTable> GetDiscerningPerfusionTimePoints(VectorVectorDouble &intensities, VariableSizeMatrixType &TransformationMatrix, VariableLengthVectorType &MeanVector);
  vtkSmartPointer<vtkTable> GetDiscerningPerfusionTimePointsForPSU(VectorVectorDouble &intensities, VariableSizeMatrixType &TransformationMatrix, VariableLengthVectorType &MeanVector);
  /**
  \brief Applies the exsiting PCA model (developed on training data) on test data
  \param intensities Test data
  */
  VectorVectorDouble ApplyPCAOnTestData(VectorVectorDouble &intensities);
  VectorVectorDouble  ApplyPCAOnTestDataWithGivenTransformations(VectorVectorDouble &intensities, VariableSizeMatrixType & TransformationMatrix, VariableLengthVectorType & MeanVector);
  /**
  \brief Calculates the transpose of an input matrix
  \param inputmatrix Input matrix
  */
  VariableSizeMatrixType MatrixTranspose(VariableSizeMatrixType &inputmatrix);

  /**
  \brief Calculates average of each feature present in input data
  \param inputdata Input matrix (in vtkTable format)
  */
  VariableLengthVectorType ComputeMeanOfGivenFeatureVectors(vtkSmartPointer< vtkTable > &inputdata);

  /**
  \brief Calculates average of each feature present in input data
  \param inputdata Input matrix (in vector format)
  */
  VariableLengthVectorType ComputeMeanOfGivenFeatureVectors(VectorVectorDouble &inputdata);

  VariableLengthVectorType ComputeMeanOfGivenFeatureVectors(vnl_matrix<double> & inputdata);

  //!Returns the PCA transformation matrix
  VariableSizeMatrixType GetPCATransformationMatrix()
  {
    return PCATransformationMatrix;
  }

  //!Returns the avergae perfusion signal
  VariableLengthVectorType GetPerfusionMeanVector()
  {
    return mPMeanvector;
  }

  /**
  \brief Sets the values of PCA transformation matrix and average perfusion signal
  \param pcaMatrix PCA transformation matrix
  \param pcaMean Avergae perfusion signal
  */
  void SetParameters(VariableSizeMatrixType &pcaMatrix, VariableLengthVectorType &pcaMean)
  {
    PCATransformationMatrix = pcaMatrix;
    mPMeanvector = pcaMean;
  }

  /**
  \brief Resets the parameters
  */
  void ResetParameters()
  {
    PCATransformationMatrix.SetSize(0, 0);
    mPMeanvector.SetSize(0);
  }

private:
  VariableSizeMatrixType PCATransformationMatrix;
  VariableLengthVectorType mPMeanvector;

};
