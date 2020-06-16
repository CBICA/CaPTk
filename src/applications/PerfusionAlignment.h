/**
\file  PerfusionAlignment.h

\brief The header file containing the PerfusionAlignment class, used to align DSC-MRI curves. 
Author: Saima Rathore
Library Dependecies: ITK 4.7+ <br>

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software/license.html

*/


#ifndef _PerfusionAlignment_h_
#define _PerfusionAlignment_h_

#include "cbicaUtilities.h"
#include "FeatureReductionClass.h"
//#include "CaPTk.h"
#include "cbicaLogging.h"
#include "CaPTkDefines.h"
#include "itkExtractImageFilter.h"
#include "NiftiDataManager.h"
#include "DicomMetadataReader.h"

#ifdef APP_BASE_CaPTk_H
#include "ApplicationBase.h"
#endif


#include <cstdio>
#include <cstdlib>
#include <vector>
#include "spline.h"


/**
\class PerfusionAlignment

\brief Calculates Aligned Perfusion Curve

*/

class PerfusionAlignment
#ifdef APP_BASE_CaPTk_H
  : public ApplicationBase
#endif
{
public:

	cbica::Logging logger;
  //! Default constructor
  PerfusionAlignment()
  {
	  logger.UseNewFile(loggerFile);
  }

  //! Default destructor
  ~PerfusionAlignment() {};
  
  NiftiDataManager mNiftiLocalPtr;

  //This function calculates characteristics of a given DSC-MRI curve
  void GetParametersFromTheCurve(std::vector<double> curve, double &base, double&drop, double &maxval, double &minval);

  //This function retrieves the value of a specific tag from the header of given dicom file
  std::string ReadMetaDataTags(std::string filepath,std::string tagstrng);
  
  //This function interpolates the given DSC-MRI curve based on th given parameters
  std::vector<double> GetInterpolatedCurve(std::vector<double> averagecurve, double timeinseconds, double totaltimeduration);

  //This is the main function of PerfusionAlignment class that is called from .cxx file with all the required parameters.  
  template< class ImageType = ImageTypeFloat3D, class PerfusionImageType = ImageTypeFloat4D >
  std::vector<typename ImageType::Pointer> Run(std::string perfImagePointerNifti, std::string dicomFile, std::string t1ceFile, int pointsbeforedrop, int pointsafterdrop, std::vector<double> & OriginalCurve, std::vector<double> & RevisedCurve,const double echotime);

//This function calculates average 3D image of the time-points of 4D DSC-MRI image specified by the start and end parameters
  template< class ImageType, class PerfusionImageType >
  std::vector<double> CalculatePerfusionVolumeMean(typename PerfusionImageType::Pointer perfImagePointerNifti, typename ImageType::Pointer ImagePointerMask, int start, int end);

  //This function calculates 3D standard deviation image of the time-points of 4D DSC-MRI image specified by the start and end parameters
  template< class ImageType = ImageTypeFloat3D, class PerfusionImageType = ImageTypeFloat4D >
  typename ImageType::Pointer CalculatePerfusionVolumeStd(typename PerfusionImageType::Pointer perfImagePointerNifti, int start, int end);

  
  template< class ImageType = ImageTypeFloat3D, class PerfusionImageType = ImageTypeFloat4D >
  typename ImageType::Pointer CalculateRCBV(typename PerfusionImageType::Pointer perfImagePointerNifti, const double TE);

  //This function returns 3D image volume of 4D DSC-MRI image specified by the index parameter
  template< class ImageType, class PerfusionImageType>
  typename ImageType::Pointer GetOneImageVolume(typename PerfusionImageType::Pointer perfImagePointerNifti, int index);
};

template< class ImageType, class PerfusionImageType >
std::vector<typename ImageType::Pointer> PerfusionAlignment::Run(std::string perfusionFile, std::string dicomFile, std::string t1ceFile, int pointsbeforedrop,int pointsafterdrop,std::vector<double> & OriginalCurve, std::vector<double> & RevisedCurve,const double echotime)
{
  std::vector<typename ImageType::Pointer> PerfusionAlignment;
  typename PerfusionImageType::Pointer perfImagePointerNifti;
  typename ImageType::Pointer t1ceImagePointer;
  try
  {
    perfImagePointerNifti = cbica::ReadImage<PerfusionImageType>(perfusionFile);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Unable to open the given DSC-MRI file. Error code : " + std::string(e1.what()));
    return PerfusionAlignment;
  }
  try
  {
    t1ceImagePointer = cbica::ReadImage< ImageType >(t1ceFile);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Unable to open the given T1 post-contrast enhanced image file. Error code : " + std::string(e1.what()));
    return PerfusionAlignment;
  }


  try
  {
    //get original curve
    typename ImageType::Pointer MASK = CalculatePerfusionVolumeStd<ImageType, PerfusionImageType>(perfImagePointerNifti, 0, 9); //values do not matter here
    std::vector<double> averagecurve = CalculatePerfusionVolumeMean<ImageType, PerfusionImageType>(perfImagePointerNifti, MASK, 0, 9); //values do not matter here
    OriginalCurve = averagecurve;
    //read dicom to get values of tags 
    std::string timeinseconds = ReadMetaDataTags(dicomFile, "0018|0080");
    typename PerfusionImageType::RegionType region = perfImagePointerNifti->GetLargestPossibleRegion();

    ////get interpolated curve
    //std::vector<double> interpolatedcurve = GetInterpolatedCurve(averagecurve,std::stof(timeinseconds),totaltimeduration);

    //std::ofstream myfile;
    //myfile.open("E:/original_curve.csv");
    //for (unsigned int index1 = 0; index1 < averagecurve.size(); index1++)
    //      myfile << std::to_string(averagecurve[index1])<< "\n";
    //myfile.close();

    // Resize
    typename PerfusionImageType::SizeType inputSize = perfImagePointerNifti->GetLargestPossibleRegion().GetSize();
    typename PerfusionImageType::SizeType outputSize;
    outputSize[0] = inputSize[0];
    outputSize[1] = inputSize[1];
    outputSize[2] = inputSize[2];
    //outputSize[3] = (std::stof(timeinseconds) * 2 / 1000) * region.GetSize()[3];
    outputSize[3] = echotime * region.GetSize()[3];

    typename PerfusionImageType::SpacingType outputSpacing;
    outputSpacing[0] = perfImagePointerNifti->GetSpacing()[0];
    outputSpacing[1] = perfImagePointerNifti->GetSpacing()[1];
    outputSpacing[2] = perfImagePointerNifti->GetSpacing()[2];
    outputSpacing[3] = perfImagePointerNifti->GetSpacing()[3] * (static_cast<double>(inputSize[3]) / static_cast<double>(outputSize[3]));
    

    typedef itk::IdentityTransform<double, 4> TransformType;
    TransformType::Pointer _pTransform = TransformType::New();
    _pTransform->SetIdentity();


    typedef itk::ResampleImageFilter<PerfusionImageType, PerfusionImageType> ResampleImageFilterType;
    typename ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();
    resample->SetInput(perfImagePointerNifti);
    resample->SetSize(outputSize);
    resample->SetOutputSpacing(outputSpacing);
    resample->SetOutputOrigin(perfImagePointerNifti->GetOrigin());
    resample->SetOutputDirection(perfImagePointerNifti->GetDirection());
    resample->SetTransform(_pTransform);
    resample->UpdateLargestPossibleRegion();

    averagecurve = CalculatePerfusionVolumeMean<ImageType, PerfusionImageType>(resample->GetOutput(), MASK, 0, 9); //values do not matter here
    RevisedCurve = averagecurve;
    double base, drop, maxcurve, mincurve;
    GetParametersFromTheCurve(averagecurve, base, drop, maxcurve, mincurve);

    std::cout << "base: " << base << " drop: " << drop << " min: " << mincurve << " max: " << maxcurve << std::endl;
    //write the corresponding perfusion 3D images
    typename PerfusionImageType::RegionType region1 = resample->GetOutput()->GetLargestPossibleRegion();
    typename PerfusionImageType::IndexType regionIndex;
    typename PerfusionImageType::SizeType regionSize;
    regionSize[0] = region1.GetSize()[0];
    regionSize[1] = region1.GetSize()[1];
    regionSize[2] = region1.GetSize()[2];
    regionSize[3] = 0;
    regionIndex[0] = 0;
    regionIndex[1] = 0;
    regionIndex[2] = 0;


    for (int index = drop - pointsbeforedrop; index < drop + pointsafterdrop; index++)
    {
      typename ImageType::Pointer NewImage = ImageType::New();
      NewImage->CopyInformation(t1ceImagePointer);
      NewImage->SetRequestedRegion(t1ceImagePointer->GetLargestPossibleRegion());
      NewImage->SetBufferedRegion(t1ceImagePointer->GetBufferedRegion());
      NewImage->Allocate();
      NewImage->FillBuffer(0);
    
      regionIndex[3] = index;
      typename PerfusionImageType::RegionType desiredRegion(regionIndex, regionSize);
      auto filter = itk::ExtractImageFilter< PerfusionImageType, ImageType >::New();
      filter->SetExtractionRegion(desiredRegion);
      filter->SetInput(resample->GetOutput());
      filter->SetDirectionCollapseToIdentity();
      filter->Update();
      typename ImageType::Pointer CurrentTimePoint = filter->GetOutput();

      itk::ImageRegionIteratorWithIndex <ImageType> imageIt(CurrentTimePoint, CurrentTimePoint->GetLargestPossibleRegion());
      itk::ImageRegionIteratorWithIndex <ImageType> newIt(NewImage, NewImage->GetLargestPossibleRegion());
      imageIt.GoToBegin(); newIt.GoToBegin();
      while (!imageIt.IsAtEnd())
      {
        newIt.Set(imageIt.Get());
        ++imageIt;
        ++newIt;
      }
      PerfusionAlignment.push_back(NewImage);
    }
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Unable to perform perfusion alignment. Error code : " + std::string(e1.what()));
    return PerfusionAlignment;
  }
  return PerfusionAlignment;
}

template< class ImageType, class PerfusionImageType>
typename ImageType::Pointer PerfusionAlignment::GetOneImageVolume(typename PerfusionImageType::Pointer perfImagePointerNifti, int index)
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

  filter->SetDirectionCollapseToIdentity(); // This is required.
  filter->Update();
  return filter->GetOutput();
}

template< class ImageType, class PerfusionImageType >
typename ImageType::Pointer PerfusionAlignment::CalculateRCBV(typename PerfusionImageType::Pointer perfImagePointerNifti,double EchoTime)
{
  typename PerfusionImageType::RegionType region = perfImagePointerNifti->GetLargestPossibleRegion();
	typename ImageType::Pointer MASK = GetOneImageVolume<ImageType, PerfusionImageType>(perfImagePointerNifti, 6);
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

	typename ImageType::Pointer A = CalculatePerfusionVolumeMean<ImageType, PerfusionImageType>(perfImagePointerNifti, 0, 9);
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
	//typedef vnl_matrix<double> MatrixType;
	//MatrixType data;
	//data.set_size(885791 , 3);
	//for (unsigned int i = 0; i < 885791; i++)
	//{
	// data(i, 0) = perfusionImageIndices(i, 0);
	// data(i, 1) = perfusionImageIndices(i, 1);
	// data(i, 2) = perfusionImageIndices(i, 2);
	//}
	//typedef itk::CSVNumericObjectFileWriter<double, 885791,3> WriterTypeVector;
	//WriterTypeVector::Pointer writerv = WriterTypeVector::New();
	//writerv->SetFileName("AllPerfusionIndices.csv");
	//writerv->SetInput(&data);
	//writerv->Write();


	//data.set_size(885791, 45);
	//for (unsigned int i = 0; i < 885791; i++)
	// for (unsigned int j = 0; j < 45; j++)
	//  data(i, j) = perfusionImage(i, j);
	//typedef itk::CSVNumericObjectFileWriter<double, 885791, 45> WriterTypeVectorData;
	//WriterTypeVectorData::Pointer writerdata = WriterTypeVectorData::New();
	//writerdata->SetFileName("AllPerfusion.csv");
	//writerdata->SetInput(&data);
	//writerdata->Write();
	//////-----------------------------------------------------------------------------------------------------------------
	VariableLengthVectorType meanPerfusionImage;
	meanPerfusionImage.SetSize(45, 1);
	for (unsigned int x = 0; x < perfusionImage.Cols(); x++)
	{
		double local_sum = 0;
		for (int y = 0; y < maskIndicesGreaterThan30Counter; y++)
			local_sum = local_sum + perfusionImage(y, x);
		meanPerfusionImage[x] = local_sum / maskIndicesGreaterThan30Counter;
	}
	//data.set_size(45,1);
	//for (unsigned int i = 0; i < meanPerfusionImage.Size(); i++)
	// data(i, 0) = meanPerfusionImage[i];
	//typedef itk::CSVNumericObjectFileWriter<double, 45,1> WriterTypeVector1;
	//WriterTypeVector1::Pointer writerv1 = WriterTypeVector1::New();
	//writerv1->SetFileName("AveragePerfusion.csv");
	//writerv1->SetInput(&data);
	//writerv1->Write();







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
	typename ImageType::Pointer rCBV = GetOneImageVolume<ImageType, PerfusionImageType>(perfImagePointerNifti, 0);
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
	////Multiply rCBV with 255 adn divide by ww
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
std::vector<double> PerfusionAlignment::CalculatePerfusionVolumeMean(typename PerfusionImageType::Pointer perfImagePointerNifti, typename ImageType::Pointer ImagePointerMask, int start, int end)
{
  std::vector<double> AverageCurve;
  typename ImageType::Pointer outputImage = cbica::GetExtractedImages<PerfusionImageType, ImageType>(perfImagePointerNifti)[0];
	ImageTypeFloat4D::RegionType region = perfImagePointerNifti->GetLargestPossibleRegion();
  for (unsigned int volumes = 0; volumes < region.GetSize()[3]; volumes++)
  {
    double volume_mean = 0;
    for (unsigned int i = 0; i < region.GetSize()[0]; i++)
      for (unsigned int j = 0; j < region.GetSize()[1]; j++)
        for (unsigned int k = 0; k < region.GetSize()[2]; k++)
        {
          typename ImageType::IndexType index3D;
          index3D[0] = i;
          index3D[1] = j;
          index3D[2] = k;

          if (ImagePointerMask.GetPointer()->GetPixel(index3D) == 0)
            continue;

          typename PerfusionImageType::IndexType index4D;
          index4D[0] = i;
          index4D[1] = j;
          index4D[2] = k;
          index4D[3] = volumes;
          volume_mean = volume_mean + perfImagePointerNifti.GetPointer()->GetPixel(index4D);
        }
    AverageCurve.push_back(std::round(volume_mean / region.GetSize()[3]));
  }
  return AverageCurve;
}

template< class ImageType, class PerfusionImageType >
typename ImageType::Pointer PerfusionAlignment::CalculatePerfusionVolumeStd(typename PerfusionImageType::Pointer perfImagePointerNifti, int start, int end)
{
  typename ImageType::Pointer outputImage = cbica::GetExtractedImages<PerfusionImageType, ImageType>(perfImagePointerNifti)[0];

  int no_of_slices = end - start + 1;
  ImageTypeFloat4D::RegionType region = perfImagePointerNifti->GetLargestPossibleRegion();
  for (unsigned int i = 0; i < region.GetSize()[0]; i++)
  {
    for (unsigned int j = 0; j < region.GetSize()[1]; j++)
    {
      for (unsigned int k = 0; k < region.GetSize()[2]; k++)
      {
        typename ImageType::IndexType index3D;
        index3D[0] = i;
        index3D[1] = j;
        index3D[2] = k;

        //calculate mean values
        double local_sum = 0;
        for (int l = 0; l < region.GetSize()[3]; l++)
        {
          typename PerfusionImageType::IndexType index4D;
          index4D[0] = i;
          index4D[1] = j;
          index4D[2] = k;
          index4D[3] = l;
          local_sum = local_sum + perfImagePointerNifti.GetPointer()->GetPixel(index4D);
        }
        double meanvalue = std::round(local_sum / no_of_slices);

        //calculate standard deviation
        double temp = 0.0;
        for (int l = 0; l < region.GetSize()[3]; l++)
        {
          typename PerfusionImageType::IndexType index4D;
          index4D[0] = i;
          index4D[1] = j;
          index4D[2] = k;
          index4D[3] = l;
          temp += (perfImagePointerNifti.GetPointer()->GetPixel(index4D) - meanvalue)*(perfImagePointerNifti.GetPointer()->GetPixel(index4D) - meanvalue);
        }
        double stddeve = std::sqrt(temp / (region.GetSize()[3] - 1));
        if(stddeve>0)
          outputImage->SetPixel(index3D, 1);
        else
          outputImage->SetPixel(index3D, 0);
      }
    }
  }
  return outputImage;
}

std::string PerfusionAlignment::ReadMetaDataTags(std::string filepaths,std::string tag)
{
  DicomMetadataReader dcmMetaReader;
  dcmMetaReader.SetFilePath(filepaths);
  bool readstatus = dcmMetaReader.ReadMetaData();
  std::string label;
  std::string value;

  if (readstatus)
  {
    bool tagFound = dcmMetaReader.GetTagValue(tag, label, value);
    if (tagFound)
      return value;
    else
      return "";
  }
  else
  {
    return "";
  }
}

std::vector<double> PerfusionAlignment::GetInterpolatedCurve(std::vector<double> averagecurve, double timeinseconds, double totaltimeduration)
{
  return averagecurve;
}
void PerfusionAlignment::GetParametersFromTheCurve(std::vector<double> curve, double &base, double&drop, double &maxval, double &minval)
{
  std::vector<double> CER;
  //curve = vq1;
  //CER = curve * 100 / (mean(curve(4:8)) - min(curve));
  //MAX(i) = mean(curve(4:8));
  //MIN(i) = min(curve);
  //Base1(i) = mean(curve(4:8));

  //Base2(i) = mean(CER(4:8));
  //CER = CER - (mean(CER(4:8))) + 300;

  //Drop = find(CER == min(CER));
  //SP_Drop(i) = Drop;
  //CER = CER(Drop - 17:Drop + 36);
  //plot(CER);      hold on;
  //CERData(i, :) = CER;
  float max = *max_element(curve.begin(), curve.end());
  float min = *min_element(curve.begin(), curve.end());


  for (int index = 0; index < curve.size(); index++)
    CER.push_back((curve[index] * 100) / (((curve[3] + curve[4] + curve[5] + curve[6] + curve[7]) / 5) - min));

  maxval = (curve[3] + curve[4] + curve[5] + curve[6] + curve[7]) / 5;
  minval = min;
  base = (curve[3] + curve[4] + curve[5] + curve[6] + curve[7]) / 5;

  double average_new_cer = (CER[3] + CER[4] + CER[5] + CER[6] + CER[7]) / 5;
  for (int index = 0; index < CER.size(); index++)
    CER[index] = CER[index] - average_new_cer + 300;

  drop = std::min_element(CER.begin(), CER.end()) - CER.begin();
}
#endif

