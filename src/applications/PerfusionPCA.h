/**
\file PerfusionPCA.h

This file holds the declaration of the class PerfusionPCA.

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or http://www.med.upenn.edu/sbia/software/license.html
*/
#ifndef _PerfusionPCA_h_
#define _PerfusionPCA_h_

#include "iostream"
#include "vtkImageData.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIterator.h"
#include "cbicaUtilities.h"
#include "FeatureReductionClass.h"
#include "vtkTable.h"
#include "vtkVariant.h"
//#include "CaPTk.h"
#include <itkExtractImageFilter.h>
#include "cbicaLogging.h"
#include "CaPTkDefines.h"
#include "CaPTkEnums.h"
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

  void PrepareNewPCAModel(const int number, const std::string inputdirectory, const std::string outputdirectory,std::vector<std::map<CAPTK::ImageModalityType, std::string>> QualifiedSubjects);
  void ApplyExistingPCAModel(const int number, const std::string inputdirectory, const std::string outputdirectory, std::vector<std::map<CAPTK::ImageModalityType, std::string>> QualifiedSubjects,const std::string ModelDirectoryName);
  
  template<class PerfusionImageType, class ImageType>
  VariableSizeMatrixType LoadPerfusionData(typename ImageType::Pointer maskImagePointerNifti, typename PerfusionImageType::Pointer perfImagePointerNifti, std::vector< typename ImageType::IndexType> &indices);

  PerfusionMapType CombineAndCalculatePerfusionPCA(PerfusionMapType PerfusionDataMap, VariableSizeMatrixType &TransformationMatrix, VariableLengthVectorType &MeanVector);
  VariableSizeMatrixType ColumnWiseScaling(VariableSizeMatrixType inputdata);

private:

};

PerfusionMapType PerfusionPCA::CombineAndCalculatePerfusionPCA(PerfusionMapType PerfusionDataMap, VariableSizeMatrixType &TransformationMatrix, VariableLengthVectorType &MeanVector)
{
  PerfusionMapType RevisedPerfusionMap;

  std::vector<int> sizes;
  VectorVectorDouble CombinedPerfusionFeaturesMap;
  for (auto const &mapiterator : PerfusionDataMap)
  {
    VariableSizeMatrixType Features = std::get<1>(mapiterator.second);
    sizes.push_back(Features.Rows());
    for (unsigned int i = 0; i < Features.Rows(); i++)
    {
      VectorDouble oneVector;
      for (unsigned int j = 0; j < 45; j++)
        oneVector.push_back(Features(i, j));
      CombinedPerfusionFeaturesMap.push_back(oneVector);
    }
  }
  FeatureReductionClass m_featureReduction;
  vtkSmartPointer<vtkTable> ReducedPCAs = m_featureReduction.GetDiscerningPerfusionTimePoints(CombinedPerfusionFeaturesMap, TransformationMatrix, MeanVector);

  int start = 0;
  for (unsigned int index = 0; index<sizes.size(); index++)// for (auto const &mapiterator : PerfusionDataMap) 
  {
    VariableSizeMatrixType OnePatietnPerfusionData;
    OnePatietnPerfusionData.SetSize(sizes[index], 10);

    if (index != 0)
      start = start + sizes[index - 1];

    for (int i = start; i < start + sizes[index]; i++)
      for (unsigned int j = 0; j < 10; j++)
        OnePatietnPerfusionData(i - start, j) = ReducedPCAs->GetValue(i, j).ToDouble();

    OnePatietnPerfusionData = ColumnWiseScaling(OnePatietnPerfusionData);
    PerfusionTupleType new_tuple(std::get<0>(PerfusionDataMap[index]), OnePatietnPerfusionData);
    RevisedPerfusionMap[index] = new_tuple;
  }
  return RevisedPerfusionMap;
}

VariableSizeMatrixType PerfusionPCA::ColumnWiseScaling(VariableSizeMatrixType inputdata)
{
  //data(:, i) = (data(:, i) - min(data(:, i))). / (max(data(:, i) - min(data(:, i))));

  int NumberOfSamples = inputdata.Rows();
  int NumberOfFeatures = inputdata.Cols();
  VariableSizeMatrixType outputdata;
  outputdata.SetSize(NumberOfSamples, NumberOfFeatures);

  //---------calculate mean and variance for each feature----------------
  VariableLengthVectorType minVector;
  VariableLengthVectorType maxVector;
  minVector.SetSize(NumberOfFeatures);
  maxVector.SetSize(NumberOfFeatures);

  for (int featureNo = 0; featureNo < NumberOfFeatures; featureNo++)
  {
    double min = inputdata(0, featureNo);
    double max = inputdata(0, featureNo);
    for (int sampleNo = 0; sampleNo < NumberOfSamples; sampleNo++)
    {
      if (inputdata(sampleNo, featureNo) < min)
        min = inputdata(sampleNo, featureNo);
      if (inputdata(sampleNo, featureNo) > max)
        max = inputdata(sampleNo, featureNo);
    }
    for (int sampleNo = 0; sampleNo < NumberOfSamples; sampleNo++)
      outputdata(sampleNo, featureNo) = ((inputdata(sampleNo, featureNo) - min) * 255) / (max - min);
  }

  return outputdata;
}

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

void PerfusionPCA::ApplyExistingPCAModel(const int number, const std::string inputdirectory, const std::string outputdirectory, std::vector<std::map<CAPTK::ImageModalityType, std::string>> trainingsubjects,const std::string modelDirectoryName)
{
}

void PerfusionPCA::PrepareNewPCAModel(const int number,const std::string inputdirectory,const std::string outputdirectory, std::vector<std::map<CAPTK::ImageModalityType, std::string>> trainingsubjects)
{
  PerfusionMapType PerfusionDataMap;


  //Extracting perfusion data of all the patients and putting in PerfusionDataMap
  for (unsigned int sid = 0; sid < trainingsubjects.size(); sid++)
  {
    std::cout << "Loading Perfusion Image: " << sid << std::endl;
    std::map<CAPTK::ImageModalityType, std::string> currentsubject = trainingsubjects[sid];
    ImageType::Pointer LabelImagePointer = cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG]));
    auto perfImagePointerNifti = cbica::ReadImage<PerfusionImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION]));
    std::vector<ImageType::IndexType> indices;

    VariableSizeMatrixType perfusionData = LoadPerfusionData<PerfusionImageType, ImageType>(LabelImagePointer, perfImagePointerNifti, indices);
    PerfusionTupleType new_tuple(indices, perfusionData);
    PerfusionDataMap[sid] = new_tuple;
  }

  //combining perfusion data, calcualting PCA
  //-----------------------------------------
  VariableSizeMatrixType TransformationMatrix;
  VariableLengthVectorType MeanVector;
  PerfusionMapType perfFeatures = CombineAndCalculatePerfusionPCA(PerfusionDataMap, TransformationMatrix, MeanVector);

  std::ofstream myfile;
  myfile.open(outputdirectory + "/PCA_PERF.csv");
  for (unsigned int index1 = 0; index1 < TransformationMatrix.Rows(); index1++)
  {
    for (unsigned int index2 = 0; index2 < TransformationMatrix.Cols(); index2++)
    {
      if (index2 == 0)
        myfile << std::to_string(TransformationMatrix[index1][index2]);
      else
        myfile << "," << std::to_string(TransformationMatrix[index1][index2]);
    }
    myfile << "\n";
  }
  myfile.close();
  myfile.open(outputdirectory + "/Mean_PERF.csv");
  for (unsigned int index1 = 0; index1 < MeanVector.Size(); index1++)
    myfile << std::to_string(MeanVector[index1]) << ",";
  myfile << "\n";
  myfile.close();

  //Putting back in images of respective patients
  //---------------------------------------------

  std::vector<std::vector<ImageType::Pointer>> RevisedPerfusionImagesOfAllPatients;

  for (unsigned int sid = 0; sid < trainingsubjects.size(); sid++)
  {
    std::cout << "Revising Perfusion Image: " << sid << std::endl;
    std::map<CAPTK::ImageModalityType, std::string> currentsubject = trainingsubjects[sid];
    auto perfImagePointerNifti = cbica::ReadImage<PerfusionImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION]));
    ImageTypeFloat4D::RegionType region = perfImagePointerNifti.GetPointer()->GetLargestPossibleRegion();
    ImageTypeFloat4D::IndexType regionIndex;
    ImageTypeFloat4D::SizeType regionSize;
    regionSize[0] = region.GetSize()[0];
    regionSize[1] = region.GetSize()[1];
    regionSize[2] = region.GetSize()[2];
    regionSize[3] = 0;
    regionIndex[0] = 0;
    regionIndex[1] = 0;
    regionIndex[2] = 0;
    regionIndex[3] = 0;

    std::vector<ImageType::Pointer> OnePatientperfusionImages;
    for (int i = 0; i < number; i++)
    {
      regionIndex[3] = i;
      ImageTypeFloat4D::RegionType desiredRegion(regionIndex, regionSize);
      auto filter = itk::ExtractImageFilter< ImageTypeFloat4D, ImageTypeFloat3D >::New();
      filter->SetExtractionRegion(desiredRegion);
      filter->SetInput(perfImagePointerNifti);
      filter->SetDirectionCollapseToIdentity();
      filter->Update();
      ImageType::Pointer CurrentTimePoint = filter->GetOutput();

      itk::ImageRegionIteratorWithIndex <ImageType> imageIt(CurrentTimePoint, CurrentTimePoint->GetLargestPossibleRegion());
      imageIt.GoToBegin();
      while (!imageIt.IsAtEnd())
      {
        imageIt.Set(0);
        ++imageIt;
      }
      std::vector<ImageType::IndexType> indices = std::get<0>(perfFeatures[sid]);
      VariableSizeMatrixType revisedPerfData = std::get<1>(perfFeatures[sid]);
      for (int j = 0; j < indices.size(); j++)
        CurrentTimePoint.GetPointer()->SetPixel(indices[j], revisedPerfData(j, i));

      OnePatientperfusionImages.push_back(CurrentTimePoint);
      //cbica::WriteImage<ImageType>(CurrentTimePoint, outputdirectory + std::to_string(sid) + "_" + std::to_string(i) + ".nii.gz");
    }
    RevisedPerfusionImagesOfAllPatients.push_back(OnePatientperfusionImages);
  }
}


#endif