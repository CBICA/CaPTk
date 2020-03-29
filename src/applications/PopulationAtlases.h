/**
\file  PopulationAtlases.h

\brief The header file containing the SurvivaPredictor class, used to find Survival Prediction Index
Library Dependencies: ITK 4.7+ <br>

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/


#ifndef _PopulationAtlases_h_
#define _PopulationAtlases_h_

//#include "CAPTk.h"
#include "NiftiDataManager.h"
#include "FeatureReductionClass.h"
#include "FeatureScalingClass.h"
#include "FeatureExtractionClass.h"
#include "itkCSVArray2DFileReader.h"
#include "itkConnectedComponentImageFilter.h"
#include "cbicaLogging.h"
#include "CaPTkEnums.h"

#ifdef APP_BASE_CAPTK_H
#include "ApplicationBase.h"
#endif

typedef itk::Image< float, 3 > ImageType;
typedef itk::CSVArray2DFileReader<double> CSVFileReaderType;
typedef vnl_matrix<double> MatrixType;


/**
\class PopulationAtlases

\brief Calculates Survival Prediction Index

Reference:

@article{macyszyn2015imaging,
title={Imaging patterns predict patient survival and molecular subtype in glioblastoma via machine learning techniques},
author={Macyszyn, Luke and Akbari, Hamed and Pisapia, Jared M and Da, Xiao and Attiah, Mark and Pigrish, Vadim and Bi, Yingtao and Pal, Sharmistha and Davuluri, Ramana V and Roccograndi, Laura and others},
journal={Neuro-oncology},
volume={18},
number={3},
pages={417--425},
year={2015},
publisher={Society for Neuro-Oncology}
}

*/
class PopulationAtlases
#ifdef APP_BASE_CAPTK_H
	: public ApplicationBase
#endif
{
public:
	//! Default constructor
	PopulationAtlases()
	{
		mLastErrorMessage = "";
	}

	//! Default destructor
	~PopulationAtlases() {};

	cbica::Logging logger;
	std::string mLastErrorMessage;

	/**
	\brief Calculates the spatial location features for given tumor images

	\param labelImagePointer		Pointer to image having segmentation labels in patient space
	\param atlasImagePointer		Pointer to image having segmentation labels in atlas space
	*/
	
	std::vector<typename ImageType::Pointer> GeneratePopualtionAtlas(const std::vector<std::string>, const std::vector<int>, const std::string, const int no_of_atlases, const std::string);
	
  bool CalculateSpatialLocationFeatures(const std::vector<std::string> image_paths,
    const std::vector<int> atlas_labels,
    const std::string inputatlas,
    const int number_of_regions,
    VariableSizeMatrixType & AllLocationFeatures,
    const std::string outputdirectory);
  
  template <class ImageType>
	typename ImageType::Pointer ReadNiftiImage(const std::string &filename);
	

  template<class ImageType>
  std::vector<typename ImageType::Pointer> RevisedTumorArea(const typename ImageType::Pointer &labelImagePointer, const int areathreshold);

  template<class ImageType>
  VectorDouble GetSpatialLocationFeaturesForFixedNumberOfRegions(const typename ImageType::Pointer &labelImagePointer,
    const typename ImageType::Pointer &jacobtemplateImagePointer,
    const int numberofregions);

  void Run()
	{

	}

private:

};

//template<class ImageType>
//VectorDouble PopulationAtlases::GetSpatialLocationFeatures(const typename ImageType::Pointer &labelImagePointer, const typename ImageType::Pointer &jacobtemplateImagePointer)
//{
//	std::vector<typename ImageType::Pointer> RevisedImages = RevisedTumorArea<ImageType>(labelImagePointer);
//	typename ImageType::Pointer tumorImage = RevisedImages[0];
//	typename ImageType::Pointer etumorImage = RevisedImages[1];
//	typename ImageType::Pointer ncrImage = RevisedImages[2];
//
//
//	typename ImageType::Pointer localizeImage = ImageType::New();
//	localizeImage->CopyInformation(labelImagePointer);
//	localizeImage->SetRequestedRegion(labelImagePointer->GetLargestPossibleRegion());
//	localizeImage->SetBufferedRegion(labelImagePointer->GetBufferedRegion());
//	localizeImage->Allocate();
//	localizeImage->FillBuffer(0);
//
//	//mulitply revised tumor image with the atlas ROI image
//	typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
//	IteratorType tumorIt(tumorImage, tumorImage->GetLargestPossibleRegion());
//	IteratorType localizeIt(localizeImage, localizeImage->GetLargestPossibleRegion());
//	IteratorType atlasIt(jacobtemplateImagePointer, jacobtemplateImagePointer->GetLargestPossibleRegion());
//
//	atlasIt.GoToBegin();
//	tumorIt.GoToBegin();
//	localizeIt.GoToBegin();
//
//	while (!tumorIt.IsAtEnd())
//	{
//		localizeIt.Set(tumorIt.Get()*atlasIt.Get());
//		++localizeIt;
//		++tumorIt;
//		++atlasIt;
//	}
//	//find number of voxels in 9 ROIs
//	VectorDouble location;
//	int tumorSize = 0;
//	for (unsigned int i = 0; i < 9; i++)
//	{
//		int counter = 0;
//		localizeIt.GoToBegin();
//		while (!localizeIt.IsAtEnd())
//		{
//			if (localizeIt.Get() == i + 1)
//			{
//				counter++;
//				tumorSize++;
//			}
//			++localizeIt;
//		}
//		location.push_back(counter);
//	}
//	//find percentages in 9 ROIs
//	for (unsigned int i = 0; i < 9; i++)
//		location[i] = (location[i] * 100) / tumorSize;
//
//	return location;
//}


template <class ImageType>
typename ImageType::Pointer PopulationAtlases::ReadNiftiImage(const std::string &filename)
{
	typedef itk::ImageFileReader<ImageType> ImageReaderType;
	typename ImageReaderType::Pointer reader = ImageReaderType::New();
	reader->SetFileName(filename);

	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject& e)
	{
		std::cerr << "Error caught: " << e.what() << "\n";
		//cbica::Logging(loggerFile, "Error caught during testing: " + std::string(e.GetDescription()));
		exit(EXIT_FAILURE);
	}
return reader->GetOutput();
}
template<class ImageType>
VectorDouble PopulationAtlases::GetSpatialLocationFeaturesForFixedNumberOfRegions(const typename ImageType::Pointer &labelImagePointer, 
  const typename ImageType::Pointer &jacobtemplateImagePointer, 
  const int numberofregions)
{
  std::vector<typename ImageType::Pointer> RevisedImages = RevisedTumorArea<ImageType>(labelImagePointer, 100);
  typename ImageType::Pointer tumorImage = RevisedImages[0];
  typename ImageType::Pointer etumorImage = RevisedImages[1];
  typename ImageType::Pointer ncrImage = RevisedImages[2];

  typename ImageType::Pointer localizeImage = ImageType::New();
  localizeImage->CopyInformation(labelImagePointer);
  localizeImage->SetRequestedRegion(labelImagePointer->GetLargestPossibleRegion());
  localizeImage->SetBufferedRegion(labelImagePointer->GetBufferedRegion());
  localizeImage->Allocate();
  localizeImage->FillBuffer(0);

  //mulitply revised tumor image with the atlas ROI image
  typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
  IteratorType tumorIt(tumorImage, tumorImage->GetLargestPossibleRegion());
  IteratorType localizeIt(localizeImage, localizeImage->GetLargestPossibleRegion());
  IteratorType atlasIt(jacobtemplateImagePointer, jacobtemplateImagePointer->GetLargestPossibleRegion());

  atlasIt.GoToBegin();
  tumorIt.GoToBegin();
  localizeIt.GoToBegin();
  int counter = 0;

  while (!tumorIt.IsAtEnd())
  {
    if (tumorIt.Get() > 0)
      counter++;
    localizeIt.Set(tumorIt.Get()*atlasIt.Get());
    ++localizeIt;
    ++tumorIt;
    ++atlasIt;
    
  }
  std::cout << "Number of tumor voxels=" << counter << std::endl;
  //find number of voxels in pre-defined number of ROIs
  VectorDouble location;
  int tumorSize = 0;
  for (int i = 0; i < numberofregions; i++)
  {
    int counter = 0;
    localizeIt.GoToBegin();
    while (!localizeIt.IsAtEnd())
    {
      if (localizeIt.Get() == i + 1)
      {
        counter++;
        tumorSize++;
      }
      ++localizeIt;
    }
    location.push_back(counter);
  }
  for (int i = 0; i < numberofregions; i++)
    std::cout << "Raw voxels in region # " << i+1 << " = " << location[i] << std::endl;
    //find percentages in predefined number of ROIs
  for (int i = 0; i < numberofregions; i++)
    location[i] = (location[i] * 100) / tumorSize;

  return location;
}

template<class ImageType>
std::vector<typename ImageType::Pointer> PopulationAtlases::RevisedTumorArea(const typename ImageType::Pointer &labelImagePointer, const int areathreshold)
{
  std::vector<typename ImageType::Pointer> RevisedImages;

  typename ImageType::Pointer tumorImage = ImageType::New();
  tumorImage->CopyInformation(labelImagePointer);
  tumorImage->SetRequestedRegion(labelImagePointer->GetLargestPossibleRegion());
  tumorImage->SetBufferedRegion(labelImagePointer->GetBufferedRegion());
  tumorImage->Allocate();
  tumorImage->FillBuffer(0);

  typename ImageType::Pointer etumorImage = ImageType::New();
  etumorImage->CopyInformation(labelImagePointer);
  etumorImage->SetRequestedRegion(labelImagePointer->GetLargestPossibleRegion());
  etumorImage->SetBufferedRegion(labelImagePointer->GetBufferedRegion());
  etumorImage->Allocate();
  etumorImage->FillBuffer(0);

  typename ImageType::Pointer ncrImage = ImageType::New();
  ncrImage->CopyInformation(labelImagePointer);
  ncrImage->SetRequestedRegion(labelImagePointer->GetLargestPossibleRegion());
  ncrImage->SetBufferedRegion(labelImagePointer->GetBufferedRegion());
  ncrImage->Allocate();
  ncrImage->FillBuffer(0);


  typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
  IteratorType imIt(labelImagePointer, labelImagePointer->GetLargestPossibleRegion());
  IteratorType ncrIt(ncrImage, ncrImage->GetLargestPossibleRegion());
  IteratorType etIt(etumorImage, etumorImage->GetLargestPossibleRegion());
  IteratorType tumorIt(tumorImage, tumorImage->GetLargestPossibleRegion());
  imIt.GoToBegin();
  tumorIt.GoToBegin();
  etIt.GoToBegin();
  ncrIt.GoToBegin();

  while (!tumorIt.IsAtEnd())
  {
    if (imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::TUMOR || imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::NONENHANCING)
      tumorIt.Set(CAPTK::VOXEL_STATUS::ON);
    else
      tumorIt.Set(CAPTK::VOXEL_STATUS::OFF);

    if (imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::TUMOR)
      etIt.Set(CAPTK::VOXEL_STATUS::ON);
    else
      etIt.Set(CAPTK::VOXEL_STATUS::OFF);

    if (imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::NONENHANCING)
      ncrIt.Set(CAPTK::VOXEL_STATUS::ON);
    else
      ncrIt.Set(CAPTK::VOXEL_STATUS::OFF);

    ++tumorIt;
    ++etIt;
    ++imIt;
    ++ncrIt;
  }

  typedef itk::Image< unsigned short, 3 > OutputImageType;

  typedef itk::ConnectedComponentImageFilter <ImageType, OutputImageType> ConnectedComponentImageFilterType;
  typename ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
  connected->FullyConnectedOn();
  connected->SetInput(tumorImage);
  connected->Update();
  OutputImageType::Pointer labeledImage = connected->GetOutput();

  connected->GetObjectCount();
  std::vector<int> sizes;
  typedef itk::ImageRegionIteratorWithIndex <OutputImageType> OutputIteratorType;
  OutputIteratorType lbimIt(labeledImage, labeledImage->GetLargestPossibleRegion());

  for (unsigned int i = 0; i < connected->GetObjectCount(); i++)
  {
    int counter = 0;
    lbimIt.GoToBegin();
    while (!lbimIt.IsAtEnd())
    {
      if (lbimIt.Get() == i + 1)
        counter++;
      ++lbimIt;
    }
    sizes.push_back(counter);
  }
  for (unsigned int i = 0; i < connected->GetObjectCount(); i++)
  {
    if (sizes[i] <= areathreshold)
    {
      lbimIt.GoToBegin();
      tumorIt.GoToBegin();
      while (!lbimIt.IsAtEnd())
      {
        if (lbimIt.Get() == i + 1)
          tumorIt.Set(CAPTK::VOXEL_STATUS::OFF);
        ++lbimIt;
        ++tumorIt;
      }
    }
  }
  tumorIt.GoToBegin();
  etIt.GoToBegin();
  ncrIt.GoToBegin();

  while (!tumorIt.IsAtEnd())
  {
    etIt.Set(etIt.Get()*tumorIt.Get());
    ncrIt.Set(ncrIt.Get()*tumorIt.Get());

    ++tumorIt;
    ++etIt;
    ++ncrIt;
  }
  RevisedImages.push_back(tumorImage);
  RevisedImages.push_back(etumorImage);
  RevisedImages.push_back(ncrImage);

  return RevisedImages;
}

#endif







