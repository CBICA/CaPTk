/**
\file PerfusionPCA.h

This file holds the declaration of the class PerfusionPCA.

https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/captk/license.html
*/

#ifndef _PerfusionPCA_h_
#define _PerfusionPCA_h_

#include "iostream"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIterator.h"
#include "cbicaUtilities.h"
#include "FeatureReductionClass.h"
#include "vtkTable.h"
#include "vtkVariant.h"
#include <itkExtractImageFilter.h>
#include "cbicaLogging.h"
#include "CaPTkDefines.h"
#include "CaPTkEnums.h"
#include "CaPTkUtils.h"
#include "itkCSVArray2DFileReader.h"
#include "cbicaITKSafeImageIO.h"

#ifdef APP_BASE_CaPTk_H
#include "ApplicationBase.h"
#endif

typedef itk::Image< float, 3 > ImageType;
typedef itk::Image< float, 4 > PerfusionImageType;
typedef itk::CSVArray2DFileReader<double> CSVFileReaderType;
typedef vnl_matrix<double> MatrixType;

typedef std::tuple< std::vector<ImageType::IndexType>, VariableSizeMatrixType> PerfusionTupleType;
typedef std::map<int, PerfusionTupleType> PerfusionMapType;


class PerfusionPCA
#ifdef APP_BASE_CaPTk_H
	: public ApplicationBase
#endif
{
public:

	enum ErrorCode
	{
		IncorrectInputStructure = 0,
		DifferentTimePoints,
		NoError
	};

  cbica::Logging logger;
  //! Default constructor
  PerfusionPCA()
  {
    logger.UseNewFile(loggerFile);
  }
  ~PerfusionPCA() {};
  
  /*template<class ImageTypeFloat4D, class ImageTypeFloat3D>
  std::vector<typename ImageTypeFloat3D::Pointer> Run(typename ImageTypeFloat3D::Pointer maskImagePointerNifti, typename ImageTypeFloat4D::Pointer perfImagePointerNifti);
*/
  ErrorCode LoadData(std::string &inValidSubject);
  bool TrainNewPerfusionModel(const int number, const std::string inputdirectory, const std::string outputdirectory);
  ErrorCode ApplyExistingPCAModel(const int number, const std::string inputdirectory, const std::string outputdirectory,const std::string ModelDirectoryName);
  
  template<class PerfusionImageType, class ImageType>
  VariableSizeMatrixType LoadPerfusionData(typename ImageType::Pointer maskImagePointerNifti, typename PerfusionImageType::Pointer perfImagePointerNifti, std::vector< typename ImageType::IndexType> &indices);

  PerfusionMapType CombineAndCalculatePerfusionPCA(PerfusionMapType PerfusionDataMap, VariableSizeMatrixType &TransformationMatrix, VariableLengthVectorType &MeanVector, vtkSmartPointer<vtkTable> &ReducedPCAs, vtkSmartPointer<vtkDoubleArray> &variance);
  VariableSizeMatrixType ColumnWiseScaling(VariableSizeMatrixType inputdata);

  PerfusionMapType CombineAndCalculatePerfusionPCAForTestData(PerfusionMapType PerfusionDataMap, VariableSizeMatrixType &TransformationMatrix, VariableLengthVectorType &MeanVector);

  //write vtk array
  void WritevtkArray(vtkDoubleArray* inputdata, std::string filepath);

  //determine # PCs based on user provided variance threshold
  int DetermineNumberOfPCsFromVariance(vtkSmartPointer<vtkDoubleArray> variance);

  void SetVarianceThreshold(float threshold);

  void SetNumberOfPCs(int pcs);

  void RequestPerfusionDataWholePopulation(bool request);

  void WriteNumberOfPCs(int n, std::string filepath);

  int ReadNumberOfPCsFromModel(std::string filepath);

  void LoadQualifiedSubjectsFromGivenDirectoryForPCA(const std::string directoryname);

  bool HasValidSubjects();
private:

	int m_TotalTimePoints = 0;

	//map of subject id to tuple of vector of mask indices and matrix of perfusion pixel data 
	//size of perfusion data map is equal to the number of subjects
	PerfusionMapType m_PerfusionDataMap;

	float m_VarianceThreshold = 0.0;
	int m_NumberOfPCs = 0;
	bool m_VarianceThresholdDefined = false;
	bool m_NumberOfPCsDefined = false;
	bool m_PerfusionDataForWholePopulationRequested = false;

	std::vector<std::map<CAPTK::ImageModalityType, std::string>> m_ValidSubjectList;
};

template<class PerfusionImageType, class ImageType>
VariableSizeMatrixType PerfusionPCA::LoadPerfusionData(typename ImageType::Pointer maskImagePointerNifti, typename PerfusionImageType::Pointer perfImagePointerNifti, std::vector<typename ImageType::IndexType> &qualifiedIndices)
{
  VectorVectorDouble pIntensities;
  //----------------------find the voxels of the mask------------------------------
  typedef itk::ImageRegionIteratorWithIndex <ImageTypeFloat3D> IteratorType;
  IteratorType maskIt(maskImagePointerNifti, maskImagePointerNifti->GetLargestPossibleRegion());

  maskIt.GoToBegin();
  int mask_counter = 0;
  while (!maskIt.IsAtEnd())
  {
    if (maskIt.Get() == 1)
    {
      mask_counter++;
      qualifiedIndices.push_back(maskIt.GetIndex());
    }
    ++maskIt;
  }
  //--------------------------populate the covariance matrix --------------------------------
  typename ImageTypeFloat4D::RegionType region = perfImagePointerNifti->GetLargestPossibleRegion();
  VariableSizeMatrixType revisedcovariancematrix;
  revisedcovariancematrix.SetSize(mask_counter, region.GetSize()[3]);
  for (size_t i = 0; i < qualifiedIndices.size(); i++)
  {
    typename ImageTypeFloat4D::IndexType regionIndex;
    regionIndex[0] = qualifiedIndices[i][0];
    regionIndex[1] = qualifiedIndices[i][1];
    regionIndex[2] = qualifiedIndices[i][2];
    for (size_t j = 0; j < region.GetSize()[3]; j++)
    {
      regionIndex[3] = j;
      revisedcovariancematrix(i, j) = perfImagePointerNifti->GetPixel(regionIndex);
    }

  }
  return revisedcovariancematrix;

}



#endif
