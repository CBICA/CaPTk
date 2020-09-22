/**
\file  PerfusionAlignment.h

\brief The header file containing the PerfusionAlignment class, used to align DSC-MRI curves. 
Author: Saima Rathore
Library Dependecies: ITK 4.7+ <br>

https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/captk/license.html

*/

#ifndef _PerfusionAlignment_h_
#define _PerfusionAlignment_h_

#include "cbicaUtilities.h"
#include "cbicaLogging.h"
#include "CaPTkDefines.h"
#include "CaPTkUtils.h"
#include "itkExtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "DicomMetadataReader.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaITKUtilities.h"
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
  template< class ImageType, class PerfusionImageType >
  std::vector<typename ImageType::Pointer> Run(std::string perfusionFile,
    int pointsbeforedrop, int pointsafterdrop,
    std::vector<double> & OriginalCurve,
    std::vector<double> & InterpolatedCurve,
    std::vector<double> & RevisedCurve,
    std::vector<double> & TruncatedCurve,
    const double timeresolution);

  //This function calculates average 3D image of the time-points of 4D DSC-MRI image specified by the start and end parameters
  template< class ImageType, class PerfusionImageType >
  std::vector<double> CalculatePerfusionVolumeMean(typename PerfusionImageType::Pointer perfImagePointerNifti, typename ImageType::Pointer ImagePointerMask, int start, int end);

  //This function normalizes the baseline value of DSC-MRI curves
  template< class PerfusionImageType >
  typename PerfusionImageType::Pointer NormalizeBaselineValueForNonZeroVoxels(typename PerfusionImageType::Pointer perfImagePointerNifti, typename PerfusionImageType::Pointer ImagePointerMask, double max_value);

  //This function calculates 3D standard deviation image of the time-points of 4D DSC-MRI image specified by the start and end parameters
  template< class ImageType = ImageTypeFloat3D, class PerfusionImageType = ImageTypeFloat4D >
  typename ImageType::Pointer CalculatePerfusionVolumeStd(typename PerfusionImageType::Pointer perfImagePointerNifti, typename ImageType::Pointer firstVolume, int start, int end);

};

template< class ImageType, class PerfusionImageType >
std::vector<typename ImageType::Pointer> PerfusionAlignment::Run(std::string perfusionFile, 
  int pointsbeforedrop,int pointsafterdrop,
  std::vector<double> & OriginalCurve, 
  std::vector<double> & InterpolatedCurve,
  std::vector<double> & RevisedCurve,
  std::vector<double> & TruncatedCurve,
  const double timeresolution)
{
  std::vector<typename ImageType::Pointer> PerfusionAlignment;
  typename PerfusionImageType::Pointer perfImagePointerNifti;

  try
  {
    perfImagePointerNifti = cbica::ReadImage<PerfusionImageType>(perfusionFile);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Unable to open the given DSC-MRI file. Error code : " + std::string(e1.what()));
    return PerfusionAlignment;
  }

  auto perfusionImageVolumes = cbica::GetExtractedImages< PerfusionImageType, ImageType >(perfImagePointerNifti);
  auto t1ceImagePointer = perfusionImageVolumes[0];

  try
  {
    //get original curve
    std::cout << "Calculating mean and std-dev from perfusion image.\n";
    typename ImageType::Pointer MASK = CalculatePerfusionVolumeStd<ImageType, PerfusionImageType>(perfImagePointerNifti, t1ceImagePointer, 0, 9); //values do not matter here
    OriginalCurve = CalculatePerfusionVolumeMean<ImageType, PerfusionImageType>(perfImagePointerNifti, MASK, 0, 9); //values do not matter here
    
    std::cout << "Started resampling.\n";
    // Resize
    auto inputSize = perfImagePointerNifti->GetLargestPossibleRegion().GetSize();
    auto outputSize = inputSize;
    outputSize[3] = timeresolution * inputSize[3];

    auto inputSpacing = perfImagePointerNifti->GetSpacing();
    auto outputSpacing = inputSpacing;
    outputSpacing[3] = inputSpacing[3] * (static_cast<double>(inputSize[3]) / static_cast<double>(outputSize[3]));

    auto resampledPerfusion = cbica::ResampleImage< PerfusionImageType >(perfImagePointerNifti, outputSpacing, outputSize);
    auto resampledPerfusion_volumes = cbica::GetExtractedImages< PerfusionImageType, ImageType >(resampledPerfusion);
    std::vector< typename ImageType::Pointer > mask_volumes(resampledPerfusion_volumes.size());
    for (size_t i = 0; i < mask_volumes.size(); i++)
    {
      mask_volumes[i] = MASK;
      mask_volumes[i]->DisconnectPipeline();
    }
    auto mask_4d = cbica::GetJoinedImage< ImageType, PerfusionImageType >(mask_volumes, resampledPerfusion->GetSpacing()[3]);

    auto masker = itk::MaskImageFilter< PerfusionImageType, PerfusionImageType >::New();
    masker->SetInput(resampledPerfusion);
    masker->SetMaskImage(mask_4d);
    masker->Update();
    resampledPerfusion = masker->GetOutput();
    resampledPerfusion_volumes = cbica::GetExtractedImages< PerfusionImageType, ImageType >(resampledPerfusion);

    InterpolatedCurve = CalculatePerfusionVolumeMean<ImageType, PerfusionImageType>(resampledPerfusion, MASK, 0, 9); //values do not matter here

    double base, drop, maxcurve, mincurve;
    GetParametersFromTheCurve(InterpolatedCurve, base, drop, maxcurve, mincurve);
    std::cout << "Curve characteristics after interpolation::: base = " << base << "; drop = " << drop << "; min = " << mincurve << "; max = " << maxcurve << std::endl;

    auto rescaledImage = NormalizeBaselineValueForNonZeroVoxels< PerfusionImageType >(resampledPerfusion, mask_4d, maxcurve);

    RevisedCurve = CalculatePerfusionVolumeMean<ImageType, PerfusionImageType>(rescaledImage, MASK, 0, 9); //values do not matter here
    GetParametersFromTheCurve(RevisedCurve, base, drop, maxcurve, mincurve);
    std::cout << "Curve characteristics after base normalization::: base = " << base << "; drop = " << drop << "; min = " << mincurve << "; max = " << maxcurve << std::endl;

    auto resample_normalized_volumes = cbica::GetExtractedImages< PerfusionImageType, ImageType >(rescaledImage);

    if ((drop - pointsbeforedrop) <= 0)
    {
      std::cerr << "Drop has been estimated at '" << drop << "' but time-points before drop is given as '" << pointsbeforedrop << "', which is not possible.\n";
      return PerfusionAlignment;
    }
    if ((drop + pointsafterdrop) >= resample_normalized_volumes.size())
    {
      std::cerr << "Drop has been estimated at '" << drop << "' and total number of time-points after resampling are '" << resample_normalized_volumes.size() << "' but time-points after drop is given as '" << pointsafterdrop << "', which is not possible.\n";
      return PerfusionAlignment;
    }

    for (unsigned int index = drop - pointsbeforedrop; index <= drop + pointsafterdrop; index++)
    {
      auto NewImage = resample_normalized_volumes[index];
      NewImage->DisconnectPipeline();

      PerfusionAlignment.push_back(NewImage);
      TruncatedCurve.push_back(RevisedCurve[index]);
    }
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Unable to perform perfusion alignment. Error code : " + std::string(e1.what()));
    return PerfusionAlignment;
  }

  auto joinedImage = cbica::GetJoinedImage< ImageType, PerfusionImageType >(PerfusionAlignment);

  return PerfusionAlignment;
}

template< class ImageType, class PerfusionImageType >
std::vector<double> PerfusionAlignment::CalculatePerfusionVolumeMean(typename PerfusionImageType::Pointer perfImagePointerNifti, typename ImageType::Pointer ImagePointerMask, int start, int end)
{
  std::vector<double> AverageCurve;
  auto perfusionImageSize = perfImagePointerNifti->GetLargestPossibleRegion().GetSize();
  for (unsigned int volumes = 0; volumes < perfusionImageSize[3]; volumes++)
  {
    double volume_mean = 0;
    size_t nonzero_counter = 0;
    for (unsigned int i = 0; i < perfusionImageSize[0]; i++)
      for (unsigned int j = 0; j < perfusionImageSize[1]; j++)
        for (unsigned int k = 0; k < perfusionImageSize[2]; k++)
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
          nonzero_counter++;
        }
    std::cout << volumes << " : " << std::round(volume_mean / nonzero_counter) << std::endl; 
    // std::cout << volumes << " : " << std::round(volume_mean / perfusionImageSize[3]) << std::endl; // original implementation
    AverageCurve.push_back(std::round(volume_mean / nonzero_counter));
  }
  return AverageCurve;
}

template< class ImageType, class PerfusionImageType >
typename ImageType::Pointer PerfusionAlignment::CalculatePerfusionVolumeStd(typename PerfusionImageType::Pointer perfImagePointerNifti, typename ImageType::Pointer firstVolume, int start, int end)
{
  typename ImageType::Pointer outputImage = firstVolume;
  auto perfusionImageSize = perfImagePointerNifti->GetLargestPossibleRegion().GetSize();
  for (unsigned int i = 0; i < perfusionImageSize[0]; i++)
  {
    for (unsigned int j = 0; j < perfusionImageSize[1]; j++)
    {
      for (unsigned int k = 0; k < perfusionImageSize[2]; k++)
      {
        typename ImageType::IndexType index3D;
        index3D[0] = i;
        index3D[1] = j;
        index3D[2] = k;

        // calculate mean values along time-scale
        double local_sum = 0;
        for (int l = 0; l < perfusionImageSize[3]; l++)
        {
          typename PerfusionImageType::IndexType index4D;
          index4D[0] = i;
          index4D[1] = j;
          index4D[2] = k;
          index4D[3] = l;
          local_sum = local_sum + perfImagePointerNifti.GetPointer()->GetPixel(index4D);
        }
        auto meanOnTimeScale = std::round(local_sum / perfusionImageSize[3]);

        // calculate standard deviation
        double temp = 0.0;
        for (int l = 0; l < perfusionImageSize[3]; l++)
        {
          typename PerfusionImageType::IndexType index4D;
          index4D[0] = i;
          index4D[1] = j;
          index4D[2] = k;
          index4D[3] = l;
          auto currentVoxel = perfImagePointerNifti.GetPointer()->GetPixel(index4D);
          temp += (currentVoxel - meanOnTimeScale)*(currentVoxel - meanOnTimeScale);
        }
        double stddeve = std::sqrt(temp / (perfusionImageSize[3] - 1));
        if (stddeve > 10)
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

  base = (curve[3] + curve[4] + curve[5] + curve[6] + curve[7]) / 5;
  maxval = *max_element(curve.begin(), curve.end());
  for (unsigned int index = 0; index < curve.size(); index++)
    if (curve[index] == 0)
      curve[index] = maxval;
  minval = *min_element(curve.begin(), curve.end());
  drop = std::min_element(curve.begin(), curve.end()) - curve.begin();

  //for (int index = 0; index < curve.size(); index++)
  //  CER.push_back((curve[index] * 100) / (((curve[3] + curve[4] + curve[5] + curve[6] + curve[7]) / 5) - min));

  //maxval = (curve[3] + curve[4] + curve[5] + curve[6] + curve[7]) / 5;
  //minval = min;
  //base = (curve[3] + curve[4] + curve[5] + curve[6] + curve[7]) / 5;

  //double average_new_cer = (CER[3] + CER[4] + CER[5] + CER[6] + CER[7]) / 5;
  //for (int index = 0; index < CER.size(); index++)
  //  CER[index] = CER[index] - average_new_cer + 300;

  //drop = std::min_element(CER.begin(), CER.end()) - CER.begin();
}
template< class PerfusionImageType >
typename PerfusionImageType::Pointer PerfusionAlignment::NormalizeBaselineValueForNonZeroVoxels(typename PerfusionImageType::Pointer perfImagePointerNifti, typename PerfusionImageType::Pointer ImagePointerMask, double max_value)
{
  double min_value = 0;
  double base_value = 300;
  auto outputImage = cbica::CreateImage< PerfusionImageType >(perfImagePointerNifti);
  itk::ImageRegionConstIterator< PerfusionImageType > inputImageIterator(perfImagePointerNifti, perfImagePointerNifti->GetLargestPossibleRegion()),
    inputMaskIterator(ImagePointerMask, ImagePointerMask->GetLargestPossibleRegion());
  itk::ImageRegionIterator< PerfusionImageType > outputImageIterator(outputImage, outputImage->GetLargestPossibleRegion());

  for (inputMaskIterator.GoToBegin(); !inputMaskIterator.IsAtEnd(); ++inputMaskIterator)
  {
    if (inputMaskIterator.Get() != 0)
    {
      auto currentIndex = inputMaskIterator.GetIndex();
      inputImageIterator.SetIndex(currentIndex);
      outputImageIterator.SetIndex(currentIndex);

      auto currentImageValue = inputImageIterator.Get();
      auto current_value = (inputImageIterator.Get() - min_value) * base_value / (max_value - min_value);
      outputImageIterator.Set(current_value);
    }
  }

  return outputImage;
}
#endif

