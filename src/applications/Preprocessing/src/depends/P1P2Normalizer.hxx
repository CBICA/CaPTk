#pragma once

#include <cmath>

#include "P1P2Normalizer.h"

#include "itkDivideImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkStatisticsImageFilter.h"

template< class TImageType >
void P1P2Normalizer< TImageType >::SetInputImage(typename TImageType::Pointer image)
{
  m_inputImage = image;
  m_algorithmDone = false;
}


template< class TImageType >
void P1P2Normalizer< TImageType >::Update()
{
  if (!m_algorithmDone)
  {
    // estimate statistics
    auto stats_originalImage = GetStatisticsForImage(m_inputImage, false);

    auto thresholder = itk::ThresholdImageFilter< TImageType >::New();
    thresholder->SetInput(m_inputImage);
    thresholder->ThresholdBelow(stats_originalImage["Mean"]);
    thresholder->SetOutsideValue(0);
    thresholder->Update();
    auto m_inputImage_thresh = thresholder->GetOutput();
    auto stats_thresholdedImage = GetStatisticsForImage(thresholder->GetOutput(), false);

    auto maskUpdater = itk::BinaryThresholdImageFilter< TImageType, TImageType >::New();
    maskUpdater->SetInput(m_inputImage);
    maskUpdater->SetLowerThreshold(stats_originalImage["Mean"]);
    maskUpdater->SetInsideValue(1);
    maskUpdater->SetOutsideValue(0);
    maskUpdater->Update();
    m_mask = maskUpdater->GetOutput();

    // initialize the histogram    
    using ImageToHistogramFilterType = itk::Statistics::MaskedImageToHistogramFilter< TImageType, TImageType >;

    const unsigned int numberOfComponents = 1; // we are always assuming a greyscale image
    typename ImageToHistogramFilterType::HistogramType::SizeType size(numberOfComponents);
    size.Fill(1000); // maximum calculated before will be the max in the histogram

    auto filter = ImageToHistogramFilterType::New();
    typename ImageToHistogramFilterType::HistogramType::MeasurementVectorType min(numberOfComponents), max(numberOfComponents);
    min.fill(stats_thresholdedImage["Min"]);
    max.fill(stats_thresholdedImage["Max"]);

    filter->SetInput(m_inputImage_thresh);
    filter->SetMaskImage(m_mask);
    filter->SetMaskValue(1);
    filter->SetHistogramSize(size);
    filter->SetHistogramBinMinimum(min);
    filter->SetHistogramBinMaximum(max);
    filter->SetMarginalScale(1000); // this is for accuracy
    filter->Update();

    auto histogram = filter->GetOutput();

    auto lower = histogram->Quantile(0, m_quantLower);
    auto upper = histogram->Quantile(0, m_quantUpper);

    auto thresholder2 = itk::ThresholdImageFilter< TImageType >::New();
    thresholder2->SetInput(m_inputImage);
    thresholder2->ThresholdBelow(lower);
    thresholder2->SetOutsideValue(lower);
    thresholder2->Update();

    auto thresholder3 = itk::ThresholdImageFilter< TImageType >::New();
    thresholder3->SetInput(thresholder2->GetOutput());
    thresholder3->ThresholdAbove(upper);
    thresholder3->SetOutsideValue(upper);
    thresholder3->Update();

    auto subtractor = itk::SubtractImageFilter< TImageType >::New();
    subtractor->SetInput1(thresholder3->GetOutput());
    subtractor->SetConstant2(lower);
    subtractor->Update();

    auto divider = itk::DivideImageFilter< TImageType, TImageType, TImageType >::New();
    divider->SetInput1(subtractor->GetOutput());
    divider->SetConstant2(upper);
    divider->Update();

    m_output = divider->GetOutput();

    m_algorithmDone = true;
  }
}

template< class TImageType >
typename TImageType::Pointer P1P2Normalizer< TImageType >::GetOutput()
{
  if (!m_algorithmDone)
  {
    Update();
  }
  return m_output;
}

template< class TImageType >
std::map< std::string, double > P1P2Normalizer< TImageType >::GetStatisticsForImage(const typename TImageType::Pointer m_inputImage, bool considerMask)
{
  std::map< std::string, double > results;

  auto statsCalculator = itk::StatisticsImageFilter< TImageType >::New();
  statsCalculator->SetInput(m_inputImage);
  try
  {
    statsCalculator->Update();
  }
  catch (const std::exception&e)
  {
    std::cerr << "Error caught during stats calculation: " << e.what() << "\n";
    return results;
  }

  results["Max"] = statsCalculator->GetMaximum();
  results["Min"] = statsCalculator->GetMinimum();
  results["Std"] = sqrtf(statsCalculator->GetVariance());
  results["Mean"] = statsCalculator->GetMean();
  
  return results;
}