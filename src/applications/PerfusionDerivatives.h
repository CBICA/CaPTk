/**
\file  PerfusionDerivatives.h

\brief The header file containing the PerfusionDerivatives class, used to calculate PSR, RCBV, and PH
Author: Saima Rathore
Library Dependecies: ITK 4.7+ <br>

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software/license.html

*/


#ifndef _PerfusionDerivatives_h_
#define _PerfusionDerivatives_h_

#include "cbicaUtilities.h"
#include "FeatureReductionClass.h"
#include "itkExtractImageFilter.h"
#include "cbicaLogging.h"
#include "CaPTkDefines.h"
#include "CaPTkUtils.h"
#include "cbicaITKSafeImageIO.h"


#ifdef APP_BASE_CaPTk_H
#include "ApplicationBase.h"
#endif

#define NO_OF_PCS 5 // total number of principal components used 


/**
\class PerfusionDerivatives

\brief Calculates Perfusion Derivatives

Reference:

@article{paulson2008comparison,
title={Comparison of dynamic susceptibility-weighted contrast-enhanced MR methods: recommendations for measuring relative cerebral blood volume in brain tumors},
author={Paulson, Eric S and Schmainda, Kathleen M},
journal={Radiology},
volume={249},
number={2},
pages={601--613},
year={2008},
publisher={Radiological Society of North America}
}

}

*/

class PerfusionDerivatives
#ifdef APP_BASE_CaPTk_H
  : public ApplicationBase
#endif
{
public:

	cbica::Logging logger;
  double baseline_start_threshold=0;
  double baseline_end_threshold = 20;
  double recovery_start_threshold=66;
  double recovery_end_threshold=88;
  //! Default constructor
  PerfusionDerivatives()
  {
	  logger.UseNewFile(loggerFile);
  }

  //! Default destructor
  ~PerfusionDerivatives() {};


  template< class ImageType = ImageTypeFloat3D, class PerfusionImageType = ImageTypeFloat4D >
  std::vector<typename ImageType::Pointer> Run(std::string perfImagePointerNifti, bool rcbv, bool psr, bool ph, const std::string outputdirectory);

  template< class ImageType = ImageTypeFloat3D, class PerfusionImageType = ImageTypeFloat4D >
  typename ImageType::Pointer CalculatePerfusionVolumeMean(typename PerfusionImageType::Pointer perfImagePointerNifti, int start, int end);

  template< class ImageType = ImageTypeFloat3D, class PerfusionImageType = ImageTypeFloat4D >
  typename ImageType::Pointer CalculatePH(typename PerfusionImageType::Pointer perfImagePointerNifti);

  template< class ImageType = ImageTypeFloat3D, class PerfusionImageType = ImageTypeFloat4D >
  typename ImageType::Pointer CalculateRCBV(typename PerfusionImageType::Pointer perfImagePointerNifti);

  template< class ImageType, class PerfusionImageType>
  typename ImageType::Pointer CalculateSignalRecovery(typename PerfusionImageType::Pointer perfImagePointerNifti);

  template< class ImageType, class PerfusionImageType>
  typename ImageType::Pointer GetOneImageVolume(typename PerfusionImageType::Pointer perfImagePointerNifti, int index);

  template< class ImageType, class PerfusionImageType >
  std::vector<double> CalculateAveragePerfusionCurve(typename PerfusionImageType::Pointer perfImagePointerNifti);

  template<class ImageType, class PerfusionImageType >
  bool IsPerfusionQualityGood(typename PerfusionImageType::Pointer perfImagePointerNifti, const std::string outputdirectory);
};

template< class ImageType, class PerfusionImageType >
std::vector<double> PerfusionDerivatives::CalculateAveragePerfusionCurve(typename PerfusionImageType::Pointer perfImagePointerNifti)
{
  typename PerfusionImageType::RegionType region = perfImagePointerNifti->GetLargestPossibleRegion();
  std::vector<double> AverageCurve(region.GetSize()[3],0);
  typedef itk::ImageRegionIteratorWithIndex <PerfusionImageType> IteratorType;
  IteratorType perfIt(perfImagePointerNifti, perfImagePointerNifti->GetLargestPossibleRegion());
  perfIt.GoToBegin();
  while (!perfIt.IsAtEnd())
  {
    AverageCurve[perfIt.GetIndex()[3]]= AverageCurve[perfIt.GetIndex()[3]]+ perfIt.Get();
    ++perfIt;
  }
  int SizeOfOneVolume = region.GetSize()[0] * region.GetSize()[1] * region.GetSize()[2];
  for(unsigned int index=0;index<AverageCurve.size();index++)
    AverageCurve[index] = std::round(AverageCurve[index]/SizeOfOneVolume);
  return AverageCurve;
}

template<class ImageType, class PerfusionImageType >
bool PerfusionDerivatives::IsPerfusionQualityGood(typename PerfusionImageType::Pointer perfImagePointerNifti,const std::string outputdirectory)
{
  std::vector<double> averagecurve = CalculateAveragePerfusionCurve<ImageType, PerfusionImageType>(perfImagePointerNifti);
  WriteCSVFiles(averagecurve, outputdirectory + "/AeragePerfusionCurve.csv");
  
  int baseline_start  = perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[3] * baseline_start_threshold/100;
  int baseline_end    = perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[3] * baseline_end_threshold/100;
  int recovery_start  = perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[3] * recovery_start_threshold/100;
  int recovery_end    = perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[3] * recovery_end_threshold/100;
  int index_min       = std::distance(averagecurve.begin(), std::min_element(averagecurve.begin(), averagecurve.end()));
  std::cout << "Number of time-points: " << perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[3] << std::endl;
  std::cout << "Baseline start: " << baseline_start << std::endl;
  std::cout << "Baseline end: " << baseline_end << std::endl;
  std::cout << "Recovery start: " << recovery_start << std::endl;
  std::cout << "Recovery end: " << recovery_end << std::endl;
  std::cout << "Drop location: " << index_min << std::endl;

  if(index_min<=baseline_end || index_min>=recovery_start)
  {
    std::cout << "Drop of the curve does not lie in between the baseline and the recovery signal." << std::endl;
    return false;
  }
  return true;
}

template< class ImageType, class PerfusionImageType >
std::vector<typename ImageType::Pointer> PerfusionDerivatives::Run(std::string perfusionFile, bool rcbv, bool psr, bool ph, const std::string outputdirectory)
{
  std::vector<typename ImageType::Pointer> perfusionDerivatives;
  typename PerfusionImageType::Pointer perfImagePointerNifti;
  try
  {
    perfImagePointerNifti = cbica::ReadImage<PerfusionImageType>(perfusionFile);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Unable to open the given DSC-MRI file. Error code : " + std::string(e1.what()));
    return perfusionDerivatives;
  }
  if (IsPerfusionQualityGood<ImageType, PerfusionImageType>(perfImagePointerNifti,outputdirectory) == false)
    return perfusionDerivatives;
  std::cout << "Perfusion curve validated!!!" << std::endl;

  perfusionDerivatives.push_back(NULL);
  perfusionDerivatives.push_back(NULL);
  perfusionDerivatives.push_back(NULL);

  typedef itk::ImageDuplicator< PerfusionImageType > DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(perfImagePointerNifti);
  duplicator->Update();
  typename PerfusionImageType::Pointer FirstCopyImage = duplicator->GetOutput();

  typename DuplicatorType::Pointer duplicator1 = DuplicatorType::New();
  duplicator1->SetInputImage(perfImagePointerNifti);
  duplicator1->Update();
  typename PerfusionImageType::Pointer SecondCopyImage = duplicator1->GetOutput();

  try
  {
    if (psr == true)
      perfusionDerivatives[0] = this->CalculateSignalRecovery<ImageType, PerfusionImageType>(FirstCopyImage);

    if (ph == true)
      perfusionDerivatives[1] = this->CalculatePH<ImageType, PerfusionImageType>(SecondCopyImage);

    if (rcbv == true)
      perfusionDerivatives[2] = this->CalculateRCBV<ImageType, PerfusionImageType>(perfImagePointerNifti);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Unable to calculate perfusion descriptors. Error code : " + std::string(e1.what()));
    return perfusionDerivatives;
  }
  return perfusionDerivatives;
}
template< class ImageType, class PerfusionImageType>
typename ImageType::Pointer PerfusionDerivatives::GetOneImageVolume(typename PerfusionImageType::Pointer perfImagePointerNifti, int index)
{
  typename PerfusionImageType::RegionType region = perfImagePointerNifti->GetLargestPossibleRegion();
  typename PerfusionImageType::IndexType regionIndex;
  typename PerfusionImageType::SizeType regionSize;
  regionSize[0] = region.GetSize()[0];
  regionSize[1] = region.GetSize()[1];
  regionSize[2] = region.GetSize()[2];
  regionSize[3] = 0;
  regionIndex[0] = 0;
  regionIndex[1] = 0;
  regionIndex[2] = 0;
  regionIndex[3] = index;
  typename PerfusionImageType::RegionType desiredRegion(regionIndex, regionSize);

  typedef itk::ExtractImageFilter< PerfusionImageType, ImageType > FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetExtractionRegion(desiredRegion);
  filter->SetInput(perfImagePointerNifti);

  filter->SetDirectionCollapseToSubmatrix(); // This is required.
  filter->Update();
  return filter->GetOutput();
}

template< class ImageType, class PerfusionImageType >
typename ImageType::Pointer PerfusionDerivatives::CalculateSignalRecovery(typename PerfusionImageType::Pointer perfImagePointerNifti)
{

  typename ImageType::Pointer A = cbica::GetExtractedImages<PerfusionImageType,ImageType>(perfImagePointerNifti)[0];
  typename ImageType::Pointer B = cbica::GetExtractedImages<PerfusionImageType, ImageType>(perfImagePointerNifti)[0];
  typename ImageType::Pointer C = cbica::GetExtractedImages<PerfusionImageType, ImageType>(perfImagePointerNifti)[0];

  for (unsigned int x = 0; x < perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[0]; x++)
    for (unsigned int y = 0; y < perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[1]; y++)
      for (unsigned int z = 0; z < perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[2]; z++)
      {
        typename PerfusionImageType::IndexType index4D;
        index4D[0] = x;
        index4D[1] = y;
        index4D[2] = z;

        typename ImageType::IndexType index3D;
        index3D[0] = x;
        index3D[1] = y;
        index3D[2] = z;
        //mean of baseline of the signal
        int baseline_start = perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[3] * baseline_start_threshold/100;
        int baseline_end = perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[3] * baseline_end_threshold/100;

        double sum = 0;
        for (unsigned int k = baseline_start; k <= baseline_end; k++)
        {
          index4D[3] = k;
          sum = sum + perfImagePointerNifti->GetPixel(index4D);
        }
        //total number of volumes is equal to the difference between the baseline_end and baseline_start, plus 1 
        A.GetPointer()->SetPixel(index3D, sum / (baseline_end-baseline_start+1));
        //---------------------------------------minimum vector------------------------------------
        std::vector<double> local_measures;
        for (unsigned int k = 0; k < perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[3]; k++)
        {
          index4D[3] = k;
          local_measures.push_back(perfImagePointerNifti->GetPixel(index4D));
        }
        double min_value = *std::min_element(std::begin(local_measures), std::end(local_measures));
        B.GetPointer()->SetPixel(index3D, min_value);

        //taking 66%- 88% of the signal
        int recovery_start = perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[3] * recovery_start_threshold/100;
        int recovery_end   = perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[3] * recovery_end_threshold/100;
        sum = 0;
        for (unsigned int k = recovery_start; k <= recovery_end; k++)
        {
          index4D[3] = k;
          sum = sum + perfImagePointerNifti->GetPixel(index4D);
        }
        C.GetPointer()->SetPixel(index3D, sum / (recovery_end-recovery_start+1));
      }

  typename ImageType::Pointer PSR = cbica::GetExtractedImages<PerfusionImageType, ImageType>(perfImagePointerNifti)[0];

  typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
  IteratorType aIt(A, A->GetLargestPossibleRegion());
  IteratorType bIt(B, B->GetLargestPossibleRegion());
  IteratorType cIt(C, C->GetLargestPossibleRegion());
  IteratorType psrIt(PSR, PSR->GetLargestPossibleRegion());

  aIt.GoToBegin();
  bIt.GoToBegin();
  cIt.GoToBegin();
  psrIt.GoToBegin();

  while (!aIt.IsAtEnd())
  {
    double val = 0;
    if ((aIt.Get() - bIt.Get()) != 0)
      val = ((cIt.Get() - bIt.Get()) / (aIt.Get() - bIt.Get())) * 255;

    psrIt.Set(val);
    ++aIt;
    ++bIt;
    ++cIt;
    ++psrIt;
  }

  //scaling the image between 0-255
  psrIt.GoToBegin();
  double min_val = psrIt.Get();
  double max_val = psrIt.Get();
  while (!psrIt.IsAtEnd())
  {
    if (psrIt.Get() > max_val)
      max_val = psrIt.Get();
    if (psrIt.Get() < min_val)
      min_val = psrIt.Get();
    ++psrIt;
  }
  std::cout << "Min:" << min_val << std::endl;
  std::cout << "Max:" << max_val << std::endl;

  typename ImageType::Pointer ScaledPSR = ImageType::New();
  ScaledPSR->CopyInformation(PSR);
  ScaledPSR->SetRequestedRegion(PSR->GetLargestPossibleRegion());
  ScaledPSR->SetBufferedRegion(PSR->GetBufferedRegion());
  ScaledPSR->Allocate();
  ScaledPSR->FillBuffer(0);
  IteratorType scaledpsrIt(ScaledPSR, ScaledPSR->GetLargestPossibleRegion());
  psrIt.GoToBegin();
  scaledpsrIt.GoToBegin();

  while (!psrIt.IsAtEnd())
  {
    if (psrIt.Get() <= 255)
      scaledpsrIt.Set(std::round(psrIt.Get()));
    else
      scaledpsrIt.Set(255);

    ++psrIt;
    ++scaledpsrIt;
  }
  return ScaledPSR;
}



template< class ImageType, class PerfusionImageType >
typename ImageType::Pointer PerfusionDerivatives::CalculatePH(typename PerfusionImageType::Pointer perfImagePointerNifti)
{
  //taking baseline of the signal
  int baseline_start  = perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[3] * baseline_start_threshold/100;
  int baseline_end    = perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[3] * baseline_end_threshold/100;

  std::vector<typename ImageType::Pointer> perfusionVolumesVectorAll = cbica::GetExtractedImages<PerfusionImageType, ImageType>(perfImagePointerNifti);
  std::vector<typename ImageType::Pointer> perfusionVolumesVectorBaselines;
  for (int x = baseline_start; x <= baseline_end; x++)
    perfusionVolumesVectorBaselines.push_back(perfusionVolumesVectorAll[x]);

  typename ImageType::Pointer A = cbica::GetExtractedImages<PerfusionImageType, ImageType>(perfImagePointerNifti)[0];
  for (unsigned int x = 0; x < perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[0]; x++)
    for (unsigned int y = 0; y < perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[1]; y++)
      for (unsigned int z = 0; z < perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[2]; z++)
      {
        typename ImageType::IndexType index;
        index[0] = x;
        index[1] = y;
        index[2] = z;
        double sum = 0;
        for (unsigned int i = 0; i < perfusionVolumesVectorBaselines.size(); i++)
          sum = sum + perfusionVolumesVectorBaselines[i]->GetPixel(index);
        A.GetPointer()->SetPixel(index, sum / (baseline_end-baseline_start+1));
      }
  //---------------------------------------minimum vector------------------------------------
  typename ImageType::Pointer B = cbica::GetExtractedImages<PerfusionImageType, ImageType>(perfImagePointerNifti)[0];
  for (unsigned int x = 0; x < perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[0]; x++)
    for (unsigned int y = 0; y < perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[1]; y++)
      for (unsigned int z = 0; z < perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[2]; z++)
      {
        typename PerfusionImageType::IndexType perfIndex;
        typename ImageType::IndexType index;
        index[0] = x;
        index[1] = y;
        index[2] = z;

        perfIndex[0] = x;
        perfIndex[1] = y;
        perfIndex[2] = z;

        std::vector<double> local_measures;
        for (unsigned int k = 0; k < perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[3]; k++)
        {
          perfIndex[3] = k;
          local_measures.push_back(perfImagePointerNifti->GetPixel(perfIndex));
        }
        double min_value = *std::min_element(std::begin(local_measures), std::end(local_measures));
        B.GetPointer()->SetPixel(index, min_value);
      }
  //---------------------------------------------------------------------------------------------------
  typename ImageType::Pointer PH = cbica::GetExtractedImages<PerfusionImageType, ImageType>(perfImagePointerNifti)[0];
  typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
  IteratorType aIt(A, A->GetLargestPossibleRegion());
  IteratorType bIt(B, B->GetLargestPossibleRegion());
  IteratorType phIt(PH, PH->GetLargestPossibleRegion());

  aIt.GoToBegin();
  bIt.GoToBegin();
  phIt.GoToBegin();

  while (!aIt.IsAtEnd())
  {
    phIt.Set(std::abs(aIt.Get() - bIt.Get()));
    ++aIt;
    ++bIt;
    ++phIt;
  }
  return PH;
}

template< class ImageType, class PerfusionImageType >
typename ImageType::Pointer PerfusionDerivatives::CalculateRCBV(typename PerfusionImageType::Pointer perfImagePointerNifti)
{
	typename PerfusionImageType::RegionType region = perfImagePointerNifti->GetLargestPossibleRegion();
  //6th time-point is used here considering the fact that perfusion signal becomes stable at this point. might replace with percentage in future
  typename ImageType::Pointer MASK = cbica::GetExtractedImages<PerfusionImageType, ImageType>(perfImagePointerNifti)[6];
	//-------------------------------
	//step 1
	typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
	IteratorType maskIt(MASK, MASK->GetLargestPossibleRegion());
	maskIt.GoToBegin();
	while (!maskIt.IsAtEnd())
	{
		if (maskIt.Get() < 30)
			maskIt.Set(0);
		++maskIt;
	}
	typedef itk::ImageFileWriter< ImageType > WriterType;
	typename WriterType::Pointer writer1 = WriterType::New();
	//writer1->SetFileName("MASK.nii.gz");
	//writer1->SetInput(MASK);
	//writer1->Update();

  //taking baseline of the signal
  int baseline_start = perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[3] * baseline_start_threshold/100;
  int baseline_end = perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[3] * baseline_end_threshold/100;

	typename ImageType::Pointer A = CalculatePerfusionVolumeMean<ImageType, PerfusionImageType>(perfImagePointerNifti, baseline_start, baseline_end);
	double eps = 2.2204e-16;
	//double eps = 0;
	typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
	IteratorType aIt(A, A->GetLargestPossibleRegion());
	aIt.GoToBegin();
	while (!aIt.IsAtEnd())
	{
		aIt.Set(aIt.Get() + eps);
		++aIt;
	}
	//writer1->SetFileName("A.nii.gz");
	//writer1->SetInput(A);
	//writer1->Update();





	////-------------------------------------
	//step 2
	double TE = 0.05;
	//for (unsigned int i = 0; i < 1; i++)
	// {
	//   typename ImageType::Pointer currentVolume = GetOneImageVolume<ImageType, PerfusionImageType>(perfImagePointerNifti, i);

	//writer1->SetFileName("FirstVolumeBeforeAnything.nii.gz");
	//writer1->SetInput(currentVolume);
	//writer1->Update();

	//   IteratorType volIt(currentVolume, currentVolume->GetLargestPossibleRegion());
	//   IteratorType aIt(A, A->GetLargestPossibleRegion());
	//   volIt.GoToBegin();
	//   aIt.GoToBegin();
	//   while (!volIt.IsAtEnd())
	//   {
	//	double currentValue;
	//	if (aIt.Get() == 0)
	//		currentValue = 0;
	//	else
	//		currentValue= volIt.Get() / aIt.Get();

	//     if (currentValue > 1)
	//       currentValue = 1;

	//     volIt.Set((-1 * std::log(currentValue)) / TE);
	//     ++volIt;
	//     ++aIt;
	//   }

	//writer1->SetFileName("FirstVolumeBeforeAddingTo4D.nii.gz");
	//writer1->SetInput(currentVolume);
	//writer1->Update();

	for (unsigned int x = 0; x < region.GetSize()[0]; x++)
	{
		for (unsigned int y = 0; y < region.GetSize()[1]; y++)
		{
			for (unsigned int z = 0; z < region.GetSize()[2]; z++)
			{
				for (unsigned int l = 0; l < region.GetSize()[3]; l++)
				{
					typename PerfusionImageType::IndexType Index4D;
					Index4D[0] = x;
					Index4D[1] = y;
					Index4D[2] = z;
					Index4D[3] = l;

					typename ImageType::IndexType Index3D;
					Index3D[0] = x;
					Index3D[1] = y;
					Index3D[2] = z;

					double currentValue = perfImagePointerNifti.GetPointer()->GetPixel(Index4D) / A.GetPointer()->GetPixel(Index3D);
					perfImagePointerNifti.GetPointer()->SetPixel(Index4D, currentValue);
				}
			}
		}
	}
	//typename ImageType::Pointer currentVolume = GetOneImageVolume<ImageType, PerfusionImageType>(perfImagePointerNifti, 0);
	//writer1->SetFileName("FirstVolumeAfterDivision.nii.gz");
	//writer1->SetInput(currentVolume);
	//writer1->Update();

	for (unsigned int x = 0; x < region.GetSize()[0]; x++)
	{
		for (unsigned int y = 0; y < region.GetSize()[1]; y++)
		{
			for (unsigned int z = 0; z < region.GetSize()[2]; z++)
			{
				for (unsigned int l = 0; l < region.GetSize()[3]; l++)
				{
					typename PerfusionImageType::IndexType Index4D;
					Index4D[0] = x;
					Index4D[1] = y;
					Index4D[2] = z;
					Index4D[3] = l;

					if (perfImagePointerNifti.GetPointer()->GetPixel(Index4D)>1)
						perfImagePointerNifti.GetPointer()->SetPixel(Index4D, 1);
				}
			}
		}
	}
	//currentVolume = GetOneImageVolume<ImageType, PerfusionImageType>(perfImagePointerNifti, 0);
	//writer1->SetFileName("FirstVolumeAfterScaling.nii.gz");
	//writer1->SetInput(currentVolume);
	//writer1->Update();

	for (unsigned int x = 0; x < region.GetSize()[0]; x++)
	{
		for (unsigned int y = 0; y < region.GetSize()[1]; y++)
		{
			for (unsigned int z = 0; z < region.GetSize()[2]; z++)
			{
				for (unsigned int l = 0; l < region.GetSize()[3]; l++)
				{
					typename PerfusionImageType::IndexType Index4D;
					Index4D[0] = x;
					Index4D[1] = y;
					Index4D[2] = z;
					Index4D[3] = l;

					double value = 0;
					if (perfImagePointerNifti.GetPointer()->GetPixel(Index4D) == 0)
						value = 0;
					else
						value = -1 * std::log(perfImagePointerNifti.GetPointer()->GetPixel(Index4D)) / TE;

					perfImagePointerNifti.GetPointer()->SetPixel(Index4D, value);
				}
			}
		}
	}

	//currentVolume = GetOneImageVolume<ImageType, PerfusionImageType>(perfImagePointerNifti, 0);
	//writer1->SetFileName("FirstVolumeAfterLogging.nii.gz");
	//writer1->SetInput(currentVolume);
	//writer1->Update();

	////----------------------convert to matrix type instead of image type-----------------------------------------------
	VariableSizeMatrixType perfusionImage;
	VariableSizeMatrixType perfusionImageIndices;
	VariableLengthVectorType maskImage;
	int maskIndicesGreaterThan30Counter = 0;

	perfusionImage.SetSize(region.GetSize()[0] * region.GetSize()[1] * region.GetSize()[2], region.GetSize()[3]);
	perfusionImageIndices.SetSize(region.GetSize()[0] * region.GetSize()[1] * region.GetSize()[2], 3);
	maskImage.SetSize(region.GetSize()[0] * region.GetSize()[1] * region.GetSize()[2], 1);
	for (unsigned int x = 0; x < region.GetSize()[0]; x++)
	{
		for (unsigned int y = 0; y < region.GetSize()[1]; y++)
		{
			for (unsigned int z = 0; z < region.GetSize()[2]; z++)
			{
				typename ImageType::IndexType index3D;
				index3D[0] = x;
				index3D[1] = y;
				index3D[2] = z;
				if (MASK.GetPointer()->GetPixel(index3D) <= 30)
					continue;

				typename PerfusionImageType::IndexType index4D;
				index4D[0] = x;
				index4D[1] = y;
				index4D[2] = z;
				for (unsigned int k = 0; k < region.GetSize()[3]; k++)
				{
					index4D[3] = k;
					perfusionImage(maskIndicesGreaterThan30Counter, k) = perfImagePointerNifti.GetPointer()->GetPixel(index4D);
				}
				perfusionImageIndices(maskIndicesGreaterThan30Counter, 0) = x;
				perfusionImageIndices(maskIndicesGreaterThan30Counter, 1) = y;
				perfusionImageIndices(maskIndicesGreaterThan30Counter, 2) = z;

				maskIndicesGreaterThan30Counter++;
			}
		}
	}
	//////-----------------------------------------------------------------------------------------------------------------
	VariableLengthVectorType meanPerfusionImage;
	meanPerfusionImage.SetSize(region.GetSize()[3], 1);
	for (unsigned int x = 0; x < perfusionImage.Cols(); x++)
	{
		double local_sum = 0;
		for (int y = 0; y < maskIndicesGreaterThan30Counter; y++)
			local_sum = local_sum + perfusionImage(y, x);
		meanPerfusionImage[x] = local_sum / maskIndicesGreaterThan30Counter;
	}
	//////----------------------------------------------------------------------------------------------------------------
	for (unsigned int x = 0; x < meanPerfusionImage.Size(); x++)
	{
		//if (isinf(meanPerfusionImage[x]) == 1)
		// meanPerfusionImage[x] = (meanPerfusionImage[x - 1] + meanPerfusionImage[x + 1]) / 2;
		if (meanPerfusionImage[x] < 0)
			meanPerfusionImage[x] = 0;
	}
	//////-----------------------------------------------------------------------------------------------------------------
	//////find maximum element
	int max_index = 0;
	double max_value = meanPerfusionImage[0];
	for (unsigned int x = 0; x < meanPerfusionImage.Size(); x++)
		if (meanPerfusionImage[x]>max_value)
		{
			max_value = meanPerfusionImage[x];
			max_index = x;
		}
	////////-----------------------------------------------------------------------------------------------------------------
	////////mean perfusion vector until maximum index
	std::vector<double> meanPerfusionUntilMaximum;
	for (int x = 0; x <= max_index; x++)
		meanPerfusionUntilMaximum.push_back(meanPerfusionImage[x]);

	//mean perfusion vector after maximum index
	std::vector<double> meanPerfusionAfterMaximum;
	for (unsigned int x = max_index; x <meanPerfusionImage.Size(); x++)
		meanPerfusionAfterMaximum.push_back(meanPerfusionImage[x]);
	//////-----------------------------------------------------------------------------------------------------------------
	//find the minimum until the maximum value
	int point1 = 0;
	double min_value = meanPerfusionUntilMaximum[0];
	for (unsigned int x = 0; x < meanPerfusionUntilMaximum.size(); x++)
		if (meanPerfusionUntilMaximum[x] < min_value)
		{
			min_value = meanPerfusionUntilMaximum[x];
			point1 = x;
		}
	//-----------------------------------------------------------------------------------------------------------------
	//find the minimum after the maximum value
	int point2 = 0;
	min_value = meanPerfusionAfterMaximum[0];
	for (unsigned int x = 0; x < meanPerfusionAfterMaximum.size(); x++)
		if (meanPerfusionAfterMaximum[x] < min_value)
		{
			min_value = meanPerfusionAfterMaximum[x];
			point2 = x;
		}
	double min_value_increased = min_value + 0.5;
	//-----------------------------------------------------------------------------------------------------------------
	//find the value between minumim and min+0.5 after the maximum value
	std::vector<int> min_index;
	for (unsigned int x = 0; x < meanPerfusionAfterMaximum.size(); x++)
		if (meanPerfusionAfterMaximum[x] < min_value_increased)
		{
			min_index.push_back(x);
			break;
		}
	point2 = max_index + min_index[0];
	//////-----------------------------------------------------------------------------------------------------------------
  typename ImageType::Pointer rCBV = cbica::GetExtractedImages<PerfusionImageType, ImageType>(perfImagePointerNifti)[0];
	std::vector<double> rCBVImage;
	for (unsigned int x = 0; x < region.GetSize()[0]; x++)
		for (unsigned int y = 0; y < region.GetSize()[1]; y++)
			for (unsigned int z = 0; z < region.GetSize()[2]; z++)
			{
				typename ImageType::IndexType index3D;
				index3D[0] = x;
				index3D[1] = y;
				index3D[2] = z;

				typename PerfusionImageType::IndexType index4D;
				index4D[0] = x;
				index4D[1] = y;
				index4D[2] = z;
				double sum = 0;
				//////Set rCBV=0 where image =0;
				if (MASK.GetPointer()->GetPixel(index3D) == 0)
				{
					rCBV->SetPixel(index3D, 0);
					rCBVImage.push_back(0);
				}
				else
				{
					for (int i = point1; i <= point2; i++)
					{
						index4D[3] = i;
						sum = sum + perfImagePointerNifti.GetPointer()->GetPixel(index4D);
					}
					rCBV->SetPixel(index3D, sum);
					rCBVImage.push_back(sum);
				}
			}
	//writer1->SetFileName("RCBV.nii.gz");
	//writer1->SetInput(rCBV);
	//writer1->Update();
	//////---------------------------------------------------------------------------------------------------------------
	std::vector<double> rCBV_copy = rCBVImage;
	std::sort(rCBV_copy.begin(), rCBV_copy.end());

	int index = std::distance(rCBV_copy.begin(), std::max_element(rCBV_copy.begin(), rCBV_copy.end()));
	double ww = std::round(rCBV_copy[std::round(index - 0.001*index)]);

	////----------------------------------------------------------------------------------------------------------------
	////Multiply rCBV with 255 and divide by ww
	IteratorType rcbvIt(rCBV, rCBV->GetLargestPossibleRegion());
	rcbvIt.GoToBegin();
	while (!rcbvIt.IsAtEnd())
	{
		rcbvIt.Set((rcbvIt.Get() * 255) / ww);
		++rcbvIt;
	}
	return rCBV;
}

template< class ImageType, class PerfusionImageType >
typename ImageType::Pointer PerfusionDerivatives::CalculatePerfusionVolumeMean(typename PerfusionImageType::Pointer perfImagePointerNifti, int start, int end)
{
  typename ImageType::Pointer outputImage = cbica::GetExtractedImages<PerfusionImageType, ImageType>(perfImagePointerNifti)[0];

	int no_of_slices = end - start + 1;
	ImageTypeFloat4D::RegionType region = perfImagePointerNifti->GetLargestPossibleRegion();
	for (unsigned int i = 0; i < region.GetSize()[0]; i++)
		for (unsigned int j = 0; j < region.GetSize()[1]; j++)
			for (unsigned int k = 0; k < region.GetSize()[2]; k++)
			{
				typename ImageType::IndexType index3D;
				index3D[0] = i;
				index3D[1] = j;
				index3D[2] = k;

				double local_sum = 0;

				for (int l = start; l <= end; l++)
				{
					typename PerfusionImageType::IndexType index4D;
					index4D[0] = i;
					index4D[1] = j;
					index4D[2] = k;
					index4D[3] = l;
					local_sum = local_sum + perfImagePointerNifti.GetPointer()->GetPixel(index4D);
				}

				double value = std::round(local_sum / no_of_slices);
				outputImage->SetPixel(index3D, value);
			}
	return outputImage;
}

#endif