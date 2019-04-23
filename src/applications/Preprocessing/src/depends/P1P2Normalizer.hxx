#pragma once

#include "P1P2Normalizer.h"

#include "itkDivideImageFilter.h"

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
    thresholder->Update();

    auto m_inputImage_meanThresh = thresholder->GetOutput();
    auto stats_thresholdedImage = GetStatisticsForImage(m_inputImage_meanThresh, false);

    m_mask = cbica::CreateImage< TImageType >(m_inputImage_meanThresh, 1);

    // initialize the histogram    
    using ImageToHistogramFilterType = itk::Statistics::MaskedImageToHistogramFilter< TImageType, TImageType >;

    const unsigned int numberOfComponents = 1; // we are always assuming a greyscale image
    typename ImageToHistogramFilterType::HistogramType::SizeType size(numberOfComponents);
    size.Fill(stats_thresholdedImage["Max"]); // maximum calculated before will be the max in the histogram

    auto filter = ImageToHistogramFilterType::New();
    typename ImageToHistogramFilterType::HistogramType::MeasurementVectorType min(numberOfComponents), max(numberOfComponents);
    min.fill(stats_thresholdedImage["Min"]);
    max.fill(stats_thresholdedImage["Max"]);

    filter->SetInput(m_inputImage_meanThresh);
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

    thresholder->SetInput(m_inputImage);
    thresholder->ThresholdBelow(lower);
    thresholder->SetOutsideValue(lower);
    thresholder->Update();

    thresholder->SetInput(thresholder->GetOutput());
    thresholder->ThresholdBelow(upper);
    thresholder->SetOutsideValue(upper);
    thresholder->Update();

    auto subtractor = itk::SubtractImageFilter< TImageType >::New();
    subtractor->SetInput1(thresholder->GetOutput());
    subtractor->SetConstant2(lower);
    subtractor->Update();

    auto divider = itk::DivideImageFilter< TImageType, TImageType, TImageType >::New();
    divider->SetInput(subtractor->GetOutput());
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
  return m_output;
}

template< class TImageType >
std::map< std::string, double > P1P2Normalizer< TImageType >::GetStatisticsForImage(const typename TImageType::Pointer m_inputImage, bool considerMask)
{
  std::map< std::string, double > results;
  std::vector< typename TImageType::PixelType > nonZeroPixels;
  // mean, stdDev, max

  TConstIteratorType  imageIterator(m_inputImage, m_inputImage->GetLargestPossibleRegion());
  imageIterator.GoToBegin();
  while (!imageIterator.IsAtEnd())
  {
    auto currentPixel = imageIterator.Get();
    if (considerMask)
    {
      if (currentPixel > 0)
      {
        nonZeroPixels.push_back(currentPixel);
      }
    }
    else
    {
      nonZeroPixels.push_back(currentPixel);
    }

    ++imageIterator;
  }
  cbica::Statistics< typename TImageType::PixelType > calculator;
  calculator.SetInput(nonZeroPixels);

  results["Max"] = calculator.GetMaximum();
  results["Min"] = calculator.GetMinimum();
  results["Std"] = calculator.GetStandardDeviation();
  results["Mean"] = calculator.GetMean();

  return results;
}