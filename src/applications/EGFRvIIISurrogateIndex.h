/**
\file  EGFRvIIISurrogateIndex.h

\brief The header file containing EGFR status estimation functionality

Library Dependecies: ITK 4.7+ <br>
Header Dependencies: CAPTk.h, FeatureReductionClass.h, ApplicationBase.h

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
see SBIA_noncommercial_license.txt file.

*/
#ifndef _EGFRvIIISurrogateIndex_h_
#define _EGFRvIIISurrogateIndex_h_

#include "CaPTkEnums.h"
#include "CaPTkDefines.h"
#include "FeatureReductionClass.h"
#include "vtkVariant.h"
#include "vtkTable.h"
#ifdef APP_BASE_CAPTK_H
#include "ApplicationBase.h"
#endif

#define EGFR_PCS 3 // number of principal components defined for PHI Estimation

/**
\class EGFRStatusPredictor

\brief A small description of the application
This application evaluates the Epidermal Growth Factor Receptor splice variant III (EGFRvIII) 
status in primary glioblastoma, by quantitative pattern analysis of the spatial heterogeneity 
of peritumoral perfusion dynamics from Dynamic Susceptibility Contrast (DSC) MRI scans, 
through the Peritumoral Heterogeneity Index.

Reference:
@article{Bakas.2017,
title={In vivo detection of EGFRvIII in glioblastoma via perfusion magnetic resonance imaging signature consistent with deep peritumoral infiltration: the phi-index},
author={S.Bakas, H.Akbari, J.Pisapia, M.Martinez-Lage, M.Rozycki, S.Rathore, N.Dahmane, D.M.O'Rourke, C.Davatzikos},
journal={Clinical Cancer Research},
pages={4724-4734},
year={2017}
}
*/
class EGFRStatusPredictor
#ifdef APP_BASE_CAPTK_H
  : public ApplicationBase
#endif
{
public:
  //! Default constructor
  EGFRStatusPredictor();

  //! Default destructor
  ~EGFRStatusPredictor();

  /**
  \brief Calculates the EGFR status based on given input image and near/far points

  \param perfImagePointer Perfusion image in NIfTI format
  \param nearIndices Indices of the near region
  \param farIndices Indices of the far region
  \param imagetype Image type whether NIfTI or DICOM
  */
  template< class ImageType = ImageTypeFloat3D, class PerfusionImageType = ImageTypeFloat4D >
  VectorDouble PredictEGFRStatus(typename PerfusionImageType::Pointer perfImagePointer, 
    std::vector<typename ImageType::IndexType> &nearIndices, 
    std::vector<typename ImageType::IndexType> &farIndices);

  /**
  \brief Loads perfusion data for the near and far indices from given NIfTI/DICOM image

  \param perfImagePointer Perfusion image in NIfTI format
  \param nearIndices Indices of the near region
  \param farIndices Indices of the far region
  \param pNearIntensities A vector of perfusion data of near indices
  \param pFarIntensities A vector of perfusion data of far indices
  \param imagetype Image type whether NIfTI or DICOM
  */
  template< class ImageType = ImageTypeFloat3D, class PerfusionImageType = ImageTypeFloat4D >
  void LoadPerfusionData(typename PerfusionImageType::Pointer perfImagePointer, 
    std::vector<typename ImageType::IndexType> &nearIndices, 
    std::vector<typename ImageType::IndexType> &farIndices,
    VectorVectorDouble & pNearIntensities, 
    VectorVectorDouble & pFarIntensities);

  /**
  \brief Removes the indices (near and far) from EGFR calculation if they do not have corresponding perfusion data
  \param rpNearIntensities A vector of perfusion data of near indices
  \param rpFarIntensities A vector of perfusion data of far indices
  \param percentageNear Percentage of near indices which have corresponding perfusion data
  \param percentageFar  Percentage of far indices which have corresponding perfusion data
  */
  void CalculateQualifiedIndices(VectorVectorDouble &rpNearIntensities, VectorVectorDouble &rpFarIntensities, 
    double & percentageNear, double & percentageFar);

  /**
  \brief Calculates an average perfusion signal
  \param PerfusionIntensities A vector of perfusion data
  \param avgSignal Average perfusion signal
  */
  void CalculateAveragePerfusionSignal(const VectorVectorDouble &PerfusionIntensities, VectorDouble &avgSignal);

  /**
  \brief Calculates maximum drop in the average perfusion signal
  \param avgSignal Average perfusion signal
  */
  double CalculateMaximumDrop(const VectorDouble &avgSignal);

  /**
  \brief Calculates BhattaCharya distance between the perfusion signals of near and far regions
  \param rpNearIntensities Perfusion data of near indices
  \param rpFarIntensities Perfusion data of far indices
  */
  double CalculateBhattacharyaCoefficient(const VectorVectorDouble &rpNearIntensities, const VectorVectorDouble &rpFarIntensities);


  /**
  \brief Calculates transpose of given input matrix
  \param  inputmatrix Input data
  */
  VariableSizeMatrixType MatrixTranspose(const VariableSizeMatrixType &inputmatrix);

  /**
  \brief Calculates covariance matrix of given input matrix
  \param  inputData Input data
  */
  VariableSizeMatrixType GetCovarianceMatrix(const VectorVectorDouble &inputData);

  /**
  \brief Calculates Cholesky factorization of given input data
  \param inputData Input data
  */
  VariableSizeMatrixType GetCholeskyFactorization(VariableSizeMatrixType &inputData);

  /**
  \brief Calculates sum of two input matrices
  \param matrix1 Input matrix1
  \param matrix2 Input matrix2
  */
  VariableSizeMatrixType GetSumOfTwoMatrice(const VariableSizeMatrixType &matrix1, const VariableSizeMatrixType &matrix2);
};

template<class ImageType, class PerfusionImageType>
VectorDouble EGFRStatusPredictor::PredictEGFRStatus(typename PerfusionImageType::Pointer perfImagePointer, 
  std::vector<typename ImageType::IndexType> &nearIndices, std::vector<typename ImageType::IndexType> &farIndices)
{
  VectorDouble EGFRStatusParams;
  EGFRStatusParams.push_back(0);
  EGFRStatusParams.push_back(0);
  EGFRStatusParams.push_back(0);
  EGFRStatusParams.push_back(0);
  EGFRStatusParams.push_back(0);

  VectorVectorDouble pNearIntensities;
  VectorVectorDouble pFarIntensities;
  VectorDouble avgNearSignal;
  VectorDouble avgFarSignal;

#ifdef APP_BASE_CAPTK_H
  messageUpdate("EGFRvIII Status Prediction");
  progressUpdate(0);
#endif

  LoadPerfusionData<ImageType, PerfusionImageType>(perfImagePointer, nearIndices, farIndices, pNearIntensities, pFarIntensities);
  double nearPercentage = 0;
  double farPercentage = 0;
  CalculateQualifiedIndices(pNearIntensities, pFarIntensities, nearPercentage, farPercentage);

  if (pNearIntensities.size() == 0 || pFarIntensities.size() == 0)
    return EGFRStatusParams;

  CalculateAveragePerfusionSignal(pNearIntensities, avgNearSignal);
  CalculateAveragePerfusionSignal(pFarIntensities, avgFarSignal);

  double maxNearDrop = CalculateMaximumDrop(avgNearSignal);
  double maxFarDrop = CalculateMaximumDrop(avgFarSignal);

#ifdef APP_BASE_CAPTK_H
  progressUpdate(30);
#endif
  FeatureReductionClass obj;
  vtkSmartPointer<vtkTable> rpNear = obj.GetDiscerningPerfusionTimePointsDynamic(pNearIntensities);
  vtkSmartPointer<vtkTable> rpFar = obj.GetDiscerningPerfusionTimePointsDynamic(pFarIntensities);


#ifdef APP_BASE_CAPTK_H
  progressUpdate(40);
#endif

  VectorVectorDouble rpNearIntensities;
  VectorVectorDouble rpFarIntensities;

  for (unsigned int i = 0; i < rpNear.GetPointer()->GetNumberOfRows(); i++)
  {
    VectorDouble oneSubjectPCs;
    for (int j = 0; j < EGFR_PCS; j++)
      oneSubjectPCs.push_back(rpNear->GetValue(i, j).ToDouble());
    rpNearIntensities.push_back(oneSubjectPCs);
  }

#ifdef APP_BASE_CAPTK_H
  progressUpdate(70);
#endif

  for (unsigned int i = 0; i < rpFar.GetPointer()->GetNumberOfRows(); i++)
  {
    VectorDouble oneSubjectPCs;
    for (int j = 0; j < EGFR_PCS; j++)
      oneSubjectPCs.push_back(rpFar->GetValue(i, j).ToDouble());
    rpFarIntensities.push_back(oneSubjectPCs);
  }
#ifdef APP_BASE_CAPTK_H
  progressUpdate(80);
#endif

  double bDistance = CalculateBhattacharyaCoefficient(rpNearIntensities, rpFarIntensities);

#ifdef APP_BASE_CAPTK_H
  progressUpdate(100);
#endif
  EGFRStatusParams[0] = bDistance;
  EGFRStatusParams[1] = maxNearDrop;
  EGFRStatusParams[2] = maxFarDrop;
  EGFRStatusParams[3] = pNearIntensities.size();
  EGFRStatusParams[4] = pFarIntensities.size();

  return EGFRStatusParams;
}

template<class ImageType, class PerfusionImageType>
void EGFRStatusPredictor::LoadPerfusionData(typename PerfusionImageType::Pointer perfImagePointer, 
  std::vector<typename ImageType::IndexType> &nearIndices, std::vector<typename ImageType::IndexType> &farIndices,
  VectorVectorDouble & pNearIntensities, VectorVectorDouble & pFarIntensities)
{
  auto timeStamps = perfImagePointer->GetLargestPossibleRegion().GetSize()[3];

  for (unsigned int i = 0; i < nearIndices.size(); i++)
  {
    VectorDouble perfusionIntensitiesPerVoxel;
   typename PerfusionImageType::IndexType perfVoxelIndex;
      perfVoxelIndex[0] = nearIndices[i][0];
      perfVoxelIndex[1] = nearIndices[i][1];
      perfVoxelIndex[2] = nearIndices[i][2];
      for (unsigned int j = 0; j < timeStamps; j++)
      {
        perfVoxelIndex[3] = j;
        perfusionIntensitiesPerVoxel.push_back(std::round(static_cast<double>(perfImagePointer.GetPointer()->GetPixel(perfVoxelIndex))));
      }
    pNearIntensities.push_back(perfusionIntensitiesPerVoxel);
  }
  for (unsigned int i = 0; i < farIndices.size(); i++)
  {
    VectorDouble perfusionIntensitiesPerVoxel;
      typename PerfusionImageType::IndexType perfVoxelIndex;
      perfVoxelIndex[0] = farIndices[i][0];
      perfVoxelIndex[1] = farIndices[i][1];
      perfVoxelIndex[2] = farIndices[i][2];

      for (unsigned int j = 0; j < timeStamps; j++)
      {
        perfVoxelIndex[3] = j;
        perfusionIntensitiesPerVoxel.push_back(std::round(static_cast<double>(perfImagePointer.GetPointer()->GetPixel(perfVoxelIndex))));
      }
    pFarIntensities.push_back(perfusionIntensitiesPerVoxel);
  }
}

#endif
