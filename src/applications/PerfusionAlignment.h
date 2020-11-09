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
#include "itkStatisticsImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkChangeInformationImageFilter.h"

#include "DicomMetadataReader.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaITKUtilities.h"
#include "cbicaStatistics.h"
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
class vtkChartXY;

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
  std::pair< std::vector<typename ImageType::Pointer>, typename ImageType::Pointer > Run(std::string perfusionFile, std::string maskFile,
    int pointsbeforedrop, int pointsafterdrop,
    std::vector<double> & OriginalCurve,
    std::vector<double> & InterpolatedCurve,
    std::vector<double> & RevisedCurve,
    std::vector<double> & TruncatedCurve,
    const double timeresolution,
    const double time_outputPerfTime,
    const float scale_maxIntensityBeforeDrop,
    const float scale_intensityDropInMeanCurve);

  //This function calculates average 3D image of the time-points of 4D DSC-MRI image specified by the start and end parameters
  template< class ImageType, class PerfusionImageType >
  std::vector<double> CalculatePerfusionVolumeMean(typename PerfusionImageType::Pointer perfImagePointerNifti, typename ImageType::Pointer ImagePointerMask);

  //This function shifts the baseline value of DSC-MRI curves
  template< class PerfusionImageType >
  typename PerfusionImageType::Pointer ShiftBaselineValueForNonZeroVoxels(typename PerfusionImageType::Pointer perfImagePointerNifti, typename PerfusionImageType::Pointer ImagePointerMask, double max_value, float baseLine);

  //This function scales the drop value of DSC-MRI curves
  template< class PerfusionImageType >
  typename PerfusionImageType::Pointer ScaleDropValue(typename PerfusionImageType::Pointer perfImagePointerNifti, typename PerfusionImageType::Pointer ImagePointerMask, double max_value, double min_curve, double multiplier);

  //This function calculates 3D standard deviation image of the time-points of 4D DSC-MRI image specified by the start and end parameters
  template< class ImageType = ImageTypeFloat3D, class PerfusionImageType = ImageTypeFloat4D >
  typename ImageType::Pointer CalculatePerfusionVolumeStd(typename PerfusionImageType::Pointer perfImagePointerNifti, typename ImageType::Pointer firstVolume, int start, int end, float stdDev_threshold);

  //! This function calculates the std-dev of the perfusion signal across the 4D timepoints
  template< class ImageType = ImageTypeFloat3D >
  std::pair<typename ImageType::Pointer, std::vector< typename ImageType::PixelType > > GetStdDevFrom4DImage(std::vector< typename ImageType::Pointer > &perfusionImageVolumes, typename ImageType::Pointer maskImage)
  {
    auto stdDevImage = cbica::CreateImage< ImageType >(perfusionImageVolumes[0]);

    auto size = maskImage->GetLargestPossibleRegion().GetSize();
    auto totalSize = 1;
    for (size_t d = 0; d < ImageType::ImageDimension; d++)
    {
      totalSize *= size[d];
    }
    std::vector< typename ImageType::PixelType > stdDevVector(totalSize);

    using IteratorType = itk::ImageRegionConstIteratorWithIndex< ImageType >;
    itk::ImageRegionIteratorWithIndex< ImageType > stdDevIterator(stdDevImage, stdDevImage->GetLargestPossibleRegion());

    // initialize the image iterators
    IteratorType maskIterator(maskImage, maskImage->GetLargestPossibleRegion());
    std::vector< IteratorType > volumeIterators(perfusionImageVolumes.size());
    for (size_t i = 0; i < perfusionImageVolumes.size(); i++)
    {
      volumeIterators[i] = IteratorType(perfusionImageVolumes[i], perfusionImageVolumes[i]->GetLargestPossibleRegion());
    }

    // iterate through each voxel and get 
    size_t sizeCounter = 0;
    for (maskIterator.GoToBegin(); !maskIterator.IsAtEnd(); ++maskIterator)
    {
      auto currentIndex = maskIterator.GetIndex();
      if (maskIterator.Get() > 0)
      {
        std::vector< typename ImageType::PixelType > currentVoxels(volumeIterators.size());
        for (size_t i = 0; i < volumeIterators.size(); i++)
        {
          volumeIterators[i].SetIndex(currentIndex);
          currentVoxels[i] = volumeIterators[i].Get();
        }

        cbica::Statistics< float >  statsCalculator(currentVoxels);
        auto currentStdDev = statsCalculator.GetStandardDeviation();
        stdDevIterator.SetIndex(currentIndex);
        stdDevIterator.Set(currentStdDev);
        stdDevVector[sizeCounter] = currentStdDev;
        m_nonZeroVoxelsPerVolume++;
      }
      else
      {
        stdDevVector[sizeCounter] = 0;
      }
      sizeCounter++;
    }

    return std::make_pair(stdDevImage, stdDevVector);
  }

  /**
  \brief fill the chart with supplied perfusion alignment data
  \param chart pointer to vtkChartXY that will be filled
  \param map pointer to map of data name to pointer of vector of data that needs to be plotted
  */
  void Createchart(vtkChartXY* chart, std::map<std::string, std::vector<float>*>* map);

private:
  bool m_negativesDetected = false;
  size_t m_nonZeroVoxelsPerVolume = 0;
};

template< class ImageType, class PerfusionImageType >
std::pair< std::vector<typename ImageType::Pointer>, typename ImageType::Pointer > PerfusionAlignment::Run(std::string perfusionFile, std::string maskFile,
  int pointsbeforedrop,int pointsafterdrop,
  std::vector<double> & OriginalCurve, 
  std::vector<double> & InterpolatedCurve,
  std::vector<double> & RevisedCurve,
  std::vector<double> & TruncatedCurve,
  const double timeresolution,
  const double time_outputPerfTime,
  const float scale_maxIntensityBeforeDrop,
  const float scale_intensityDropInMeanCurve)
{
  std::vector<typename ImageType::Pointer> PerfusionAlignment;
  typename PerfusionImageType::Pointer perfImagePointerNifti;

  typename ImageType::Pointer maskImage = ImageType::New();

  try
  {
    perfImagePointerNifti = cbica::ReadImage<PerfusionImageType>(perfusionFile);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Unable to open the given DSC-MRI file. Error code : " + std::string(e1.what()));
    return std::make_pair(PerfusionAlignment, maskImage);
  }

  auto inputSpacing = perfImagePointerNifti->GetSpacing();
  inputSpacing[3] = timeresolution;

  auto infoChanger = itk::ChangeInformationImageFilter< PerfusionImageType >::New();
  infoChanger->SetInput(perfImagePointerNifti);
  infoChanger->SetOutputSpacing(inputSpacing);
  infoChanger->ChangeSpacingOn();
  infoChanger->Update();
  perfImagePointerNifti = infoChanger->GetOutput();

  auto perfusionImageVolumes = cbica::GetExtractedImages< PerfusionImageType, ImageType >(perfImagePointerNifti);
  auto t1ceImagePointer = perfusionImageVolumes[0];

  try
  {
    //get original curve
    std::cout << "Calculating mean and std-dev from perfusion image.\n";
    if (cbica::isFile(maskFile))
    {
      maskImage = cbica::ReadImage< ImageType >(maskFile);
      if (!cbica::ImageSanityCheck< ImageType >(t1ceImagePointer, maskImage, 0.1, 0))
      {
        std::cerr << "Mask and first perfusion volume are not in the same space, please check again or choose one of the default masking options.\n";
        return std::make_pair(PerfusionAlignment, maskImage);
      }
    }
    else
    {
      auto thresholder = itk::OtsuThresholdImageFilter< ImageType, ImageType >::New();
      thresholder->SetInput(t1ceImagePointer);
      thresholder->SetOutsideValue(1);
      thresholder->SetInsideValue(0);
      thresholder->Update();
      maskImage = thresholder->GetOutput();
    }

    if (maskImage.IsNull())
    {
      std::cerr << "Mask image needs to be defined at this point.\n";
      return std::make_pair(PerfusionAlignment, maskImage);
    }

    //maskImage = CalculatePerfusionVolumeStd<ImageType, PerfusionImageType>(perfImagePointerNifti, t1ceImagePointer, 0, 9, stdDev);
    // put an error check here
    // if numberOfNonZeroVoxelsInMask > 0.5 * totalNumberOfVoxels
    // print_error << "Warning: the mask is larger than expected volume of brain, please perform quality-check after process completion.\n"
    OriginalCurve = CalculatePerfusionVolumeMean<ImageType, PerfusionImageType>(perfImagePointerNifti, maskImage); 
    std::cout << "Started resampling.\n";
    // Resize
    auto outputSpacing = inputSpacing;
    outputSpacing[3] = time_outputPerfTime;

    auto resampledPerfusion = cbica::ResampleImage< PerfusionImageType >(perfImagePointerNifti, outputSpacing, "linear", true);
    auto resampledPerfusion_volumes = cbica::GetExtractedImages< PerfusionImageType, ImageType >(resampledPerfusion);
    
    InterpolatedCurve = CalculatePerfusionVolumeMean<ImageType, PerfusionImageType>(resampledPerfusion, maskImage);

    double base, drop, maxcurve, mincurve;
    GetParametersFromTheCurve(InterpolatedCurve, base, drop, maxcurve, mincurve);
    std::cout << "Curve characteristics after interpolation::: base = " << base << "; drop = " << drop << "; min = " << mincurve << "; max = " << maxcurve << std::endl;

    resampledPerfusion = ScaleDropValue< PerfusionImageType >(resampledPerfusion, cbica::CreateImage< PerfusionImageType >(resampledPerfusion, 1), base, mincurve, scale_intensityDropInMeanCurve);
    RevisedCurve = CalculatePerfusionVolumeMean<ImageType, PerfusionImageType>(resampledPerfusion, maskImage);
    GetParametersFromTheCurve(RevisedCurve, base, drop, maxcurve, mincurve);
    std::cout << "Curve characteristics after drop scaling ::: base = " << base << "; drop = " << drop << "; min = " << mincurve << "; max = " << maxcurve << std::endl;

    auto shiftedImage = ShiftBaselineValueForNonZeroVoxels< PerfusionImageType >(resampledPerfusion, cbica::CreateImage< PerfusionImageType >(resampledPerfusion, 1), base, scale_maxIntensityBeforeDrop);
    auto shifted_normalized_volumes = cbica::GetExtractedImages< PerfusionImageType, ImageType >(shiftedImage);

    RevisedCurve = CalculatePerfusionVolumeMean<ImageType, PerfusionImageType>(shiftedImage, maskImage); 
    GetParametersFromTheCurve(RevisedCurve, base, drop, maxcurve, mincurve);
    std::cout << "Curve characteristics after base shifting::: base = " << base << "; drop = " << drop << "; min = " << mincurve << "; max = " << maxcurve << std::endl;

    if ((drop - pointsbeforedrop) <= 0)
    {
      std::cerr << "Drop has been estimated at '" << drop << "' but time-points before drop is given as '" << pointsbeforedrop << "', which is not possible.\n";
      return std::make_pair(PerfusionAlignment, maskImage);
    }
    if ((drop + pointsafterdrop) >= shifted_normalized_volumes.size())
    {
      std::cerr << "Drop has been estimated at '" << drop << "' and total number of time-points after resampling are '" << shifted_normalized_volumes.size() << "' but time-points after drop is given as '" << pointsafterdrop << "', which is not possible.\n";
      return std::make_pair(PerfusionAlignment, maskImage);
    }

    if (m_negativesDetected)
    {
      std::cerr << "WARNING: There are negative values in the shifted image, which might lead to erroneous results during downstream analyses; please check the baseline value.\n";
    }

    for (unsigned int index = drop - pointsbeforedrop; index <= drop + pointsafterdrop; index++)
    {
      auto NewImage = shifted_normalized_volumes[index];
      NewImage->DisconnectPipeline();

      PerfusionAlignment.push_back(NewImage);
      TruncatedCurve.push_back(RevisedCurve[index]);
    }
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Unable to perform perfusion alignment. Error code : " + std::string(e1.what()));
    return std::make_pair(PerfusionAlignment, maskImage);
  }

  return std::make_pair(PerfusionAlignment, maskImage);
}

template< class ImageType, class PerfusionImageType >
std::vector<double> PerfusionAlignment::CalculatePerfusionVolumeMean(typename PerfusionImageType::Pointer perfImagePointerNifti, typename ImageType::Pointer ImagePointerMask)
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
    AverageCurve.push_back(std::round(volume_mean / nonzero_counter));
  }
  return AverageCurve;
}

template< class ImageType, class PerfusionImageType >
typename ImageType::Pointer PerfusionAlignment::CalculatePerfusionVolumeStd(typename PerfusionImageType::Pointer perfImagePointerNifti, typename ImageType::Pointer firstVolume, int start, int end, float stdDev_threshold)
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
        if (stddeve > stdDev_threshold)
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
typename PerfusionImageType::Pointer PerfusionAlignment::ShiftBaselineValueForNonZeroVoxels(typename PerfusionImageType::Pointer perfImagePointerNifti, typename PerfusionImageType::Pointer ImagePointerMask, double base_average, float baseLine)
{
  auto outputImage = perfImagePointerNifti;
  outputImage->DisconnectPipeline();
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

      auto output_value = (inputImageIterator.Get() - base_average) + baseLine;
      if (!m_negativesDetected)
      {
        if (output_value < 0)
        {
          m_negativesDetected = true;
        }
      }
      outputImageIterator.Set(output_value);
    }
  }
  return outputImage;
}
template< class PerfusionImageType >
typename PerfusionImageType::Pointer PerfusionAlignment::ScaleDropValue(typename PerfusionImageType::Pointer perfImagePointerNifti, typename PerfusionImageType::Pointer ImagePointerMask, double base_average, double min_curve, double multiplier)
{
  auto outputImage = perfImagePointerNifti;
  outputImage->DisconnectPipeline();
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

      auto output_value = (inputImageIterator.Get() * multiplier) / (base_average - min_curve);
      outputImageIterator.Set(output_value);
    }
  }
  return outputImage;
}

#endif

