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
#include "cbicaLogging.h"
#include "CaPTkDefines.h"
#include "CaPTkUtils.h"
#include "itkExtractImageFilter.h"
#include "DicomMetadataReader.h"
#include "cbicaITKSafeImageIO.h"
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
  
  //This function calculates characteristics of a given DSC-MRI curve
  void GetParametersFromTheCurve(std::vector<double> curve, double &base, double&drop, double &maxval, double &minval);

  //This function retrieves the value of a specific tag from the header of given dicom file
  std::string ReadMetaDataTags(std::string filepath,std::string tagstrng);
  
  //This function interpolates the given DSC-MRI curve based on th given parameters
  std::vector<double> GetInterpolatedCurve(std::vector<double> averagecurve, double timeinseconds, double totaltimeduration);

  //This is the main function of PerfusionAlignment class that is called from .cxx file with all the required parameters.  
  template< class ImageType = ImageTypeFloat3D, class PerfusionImageType = ImageTypeFloat4D >
  std::vector<typename ImageType::Pointer> Run(std::string perfImagePointerNifti, std::string t1ceFile, int pointsbeforedrop, int pointsafterdrop, std::vector<double> & OriginalCurve, std::vector<double> & RevisedCurve,const double echotime);

//This function calculates average 3D image of the time-points of 4D DSC-MRI image specified by the start and end parameters
  template< class ImageType, class PerfusionImageType >
  std::vector<double> CalculatePerfusionVolumeMean(typename PerfusionImageType::Pointer perfImagePointerNifti, typename ImageType::Pointer ImagePointerMask, int start, int end);

  //This function calculates 3D standard deviation image of the time-points of 4D DSC-MRI image specified by the start and end parameters
  template< class ImageType = ImageTypeFloat3D, class PerfusionImageType = ImageTypeFloat4D >
  typename ImageType::Pointer CalculatePerfusionVolumeStd(typename PerfusionImageType::Pointer perfImagePointerNifti, int start, int end);

};

template< class ImageType, class PerfusionImageType >
std::vector<typename ImageType::Pointer> PerfusionAlignment::Run(std::string perfusionFile, std::string t1ceFile, int pointsbeforedrop,int pointsafterdrop,std::vector<double> & OriginalCurve, std::vector<double> & RevisedCurve,const double echotime)
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
    //std::string timeinseconds = ReadMetaDataTags(dicomFile, "0018|0080");
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

