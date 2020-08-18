/**
\file FeatureExtraction.hxx
\brief Contains the implementations of class FeatureExtraction
*/

#pragma once

//#include "FeatureExtraction.h"
//#include "itkDOMNodeXMLReader.h"
//#include "itkDOMNodeXMLWriter.h"
//#include "itkLabelStatisticsImageFilter.h"
//#include "itkLabelGeometryImageFilter.h"
//#include "itkBinaryImageToLabelMapFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkMaskedImageToHistogramFilter.h"
//#include "itkRegionOfInterestImageFilter.h"
#include "itkRoundImageFilter.h"

#include "itkEnhancedHistogramToRunLengthFeaturesFilter.h"
#include "itkEnhancedScalarImageToRunLengthFeaturesFilter.h"

#include "itkHistogramToRunLengthFeaturesFilter.h"
#include "itkScalarImageToRunLengthFeaturesFilter.h"
#include "itkBoundingBox.h"


//#include "itkOpenCVImageBridge.h"

#include "NGLDMFeatures.h"
#include "NGTDMFeatures.h"
#include "MorphologicFeatures.h"
#include "HistogramFeatures.h"
#include "FractalBoxCount.h"
#include "GaborWavelets.h"
#include "LawsMeasures.h"
#include "EdgeEnhancement.h"
#include "PowerSpectrum.h"
#include "LBPMeasures.h"
#include "GLSZMFeatures.h"
#include "GLCMFeatures.h"
#include "GLRLMFeatures.h"
//#include "FractalBoxCount_template.h"

#include "cbicaProgressBar.h"

//TBD
#include <math.h> //for debugging
//TBD

#ifdef _OPENMP
#include <omp.h>
#endif

template< class TImage >
float FeatureExtraction< TImage >::GetMaximumDistanceWithinTheDefinedROI(const typename TImage::Pointer image, const typename TImage::Pointer mask)
{
  // ref: https://github.com/MachadoLF/TextureAnalysisExtension/blob/603aac5c58d6e95e510300c7c67a6b87c143c5d1/TextureProcessing/RunLengthFeat/RunLengthFeat.cxx#L73
  using BoundingBoxType = itk::BoundingBox<unsigned long, TImage::ImageDimension >;
  auto bbox = BoundingBoxType::New();
  auto points = BoundingBoxType::PointsContainer::New();
  itk::Point< float, TImage::ImageDimension > point;

  unsigned int idx = 0;

  itk::ImageRegionIteratorWithIndex< TImage > ItI(image, image->GetLargestPossibleRegion());

  for (ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI)
  {
    if (mask->GetPixel(ItI.GetIndex()) == 1)
    {
      image->TransformIndexToPhysicalPoint(ItI.GetIndex(), point);
      points->InsertElement(idx++, point);
    }
  }
  bbox->SetPoints(points);
  bbox->ComputeBoundingBox();
  auto pointMin = bbox->GetMinimum();
  auto pointMax = bbox->GetMaximum();
  return pointMin.EuclideanDistanceTo(pointMax);
}

template< class TImage >
template< class TVolumeImage >
void FeatureExtraction< TImage >::CalculateVolumetric(const typename TVolumeImage::Pointer mask, std::map<std::string, double>& featurevec)
{
  int count = 0;
  itk::ImageRegionIteratorWithIndex< TVolumeImage > interIt(mask, mask->GetLargestPossibleRegion());
  for (interIt.GoToBegin(); !interIt.IsAtEnd(); ++interIt)
  {
    if (interIt.Get() > 0)
      count++;
  }
  auto spacing = mask->GetSpacing();
  double voxvol = 1;
  for (size_t i = 0; i < TVolumeImage::ImageDimension; i++)
  {
    voxvol *= spacing[i];
  }
  double volume = count * voxvol;
  featurevec["Pixels"] = count;

  if (TVolumeImage::ImageDimension == 3)
  {
    featurevec["Volume"] = volume;
  }
  else if (TVolumeImage::ImageDimension == 2)
  {
    featurevec["Area"] = volume; // this is actually the area since the computation has happened in 2D
  }
  else
  {
    // do nothing
  }
}


template< class TImage >
void FeatureExtraction< TImage >::CalculateFractalDimensions(const typename TImage::Pointer itkImage, std::map< std::string, double >& featurevec, bool latticePatch)
{
  FractalBoxCount< TImage > fractalDimensionCalculator;
  fractalDimensionCalculator.SetInputImage(itkImage);
  fractalDimensionCalculator.SetRadius(m_Radius);
  fractalDimensionCalculator.SetLatticePointStatus(latticePatch);
  fractalDimensionCalculator.SetStartingIndex(m_currentLatticeStart);
  if (m_debug)
  {
    fractalDimensionCalculator.EnableDebugMode();
  }
  fractalDimensionCalculator.Update();
  auto temp = fractalDimensionCalculator.GetOutput();
  for (auto const& f : temp)
  {
    featurevec[f.first] = f.second;
  }
}


template< class TImage >
void FeatureExtraction< TImage >::CalculateLawsMeasures(const typename TImage::Pointer itkImage, std::map< std::string, double >& featurevec)
{
  LawsMeasures< TImage > lawsMeasuresCalculator;
  lawsMeasuresCalculator.SetInputImage(itkImage);
  lawsMeasuresCalculator.SetStartingIndex(m_currentLatticeStart);
  if (m_debug)
  {
    lawsMeasuresCalculator.EnableDebugMode();
  }
  lawsMeasuresCalculator.Update();
  auto temp = lawsMeasuresCalculator.GetOutput();
  for (auto const& f : temp)
  {
    featurevec[f.first] = f.second;
  }
}


template< class TImage >
void FeatureExtraction< TImage >::CalculateEdgeEnhancement(const typename TImage::Pointer itkImage, std::map< std::string, double >& featurevec)
{
  EdgeEnhancement< TImage > edgeEnhancementCalculator;
  edgeEnhancementCalculator.SetInputImage(itkImage);
  if (m_Radius == -1)
  {
    edgeEnhancementCalculator.SetRadius(m_Radius_float);
  }
  else
  {
    edgeEnhancementCalculator.SetRadius(m_Radius);
  }
  edgeEnhancementCalculator.SetStartingIndex(m_currentLatticeStart);
  edgeEnhancementCalculator.SetETA(m_edgesETA);
  edgeEnhancementCalculator.SetEpsilon(m_edgesEpsilon);
  if (m_debug)
  {
    edgeEnhancementCalculator.EnableDebugMode();
  }
  edgeEnhancementCalculator.Update();
  auto temp = edgeEnhancementCalculator.GetOutput();
  for (auto const& f : temp)
  {
    featurevec[f.first] = f.second;
  }
}


template< class TImage >
void FeatureExtraction< TImage >::CalculateLBP(const typename TImage::Pointer itkImage, const typename TImage::Pointer mask, std::map< std::string, double >& featurevec)
{
  LBPMeasures< TImage > lbpCalculator;
  if (m_Radius == -1)
  {
    lbpCalculator.SetRadius(m_Radius_float);
  }
  else
  {
    lbpCalculator.SetRadius(m_Radius);
  }
  lbpCalculator.SetInputImage(itkImage);
  lbpCalculator.SetInputMask(mask);
  lbpCalculator.SetNeighbors(m_neighborhood);
  lbpCalculator.SetLBPStyle(m_LBPStyle);
  if (m_debug)
  {
    lbpCalculator.EnableDebugMode();
  }
  lbpCalculator.Update();
  auto temp = lbpCalculator.GetOutput();
  for (auto const& f : temp)
  {
    featurevec[f.first] = f.second;
  }
}


template< class TImage >
void FeatureExtraction< TImage >::CalculatePowerSpectrum(const typename TImage::Pointer itkImage, std::map< std::string, double >& featurevec)
{
  PowerSpectrum< TImage > powerSpectrumCalculator;
  powerSpectrumCalculator.SetInputImage(itkImage);
  powerSpectrumCalculator.SetCenter(m_centerIndexString);
  powerSpectrumCalculator.SetStartingIndex(m_currentLatticeStart);
  if (m_debug)
  {
    powerSpectrumCalculator.EnableDebugMode();
  }
  powerSpectrumCalculator.Update();
  auto temp = powerSpectrumCalculator.GetOutput();
  for (auto const& f : temp)
  {
    featurevec[f.first] = f.second;
  }
}


template< class TImage >
void FeatureExtraction< TImage >::CalculateGaborWavelets(const typename TImage::Pointer itkImage, std::map< std::string, double >& featurevec, bool latticePatch)
{
  GaborWavelets< TImage > gaborWaveletCalculator;
  gaborWaveletCalculator.SetInputImage(itkImage);
  if (m_Radius == -1)
  {
    gaborWaveletCalculator.SetRadius(m_Radius_float);
  }
  else
  {
    gaborWaveletCalculator.SetRadius(m_Radius);
  }
  gaborWaveletCalculator.SetStartingIndex(m_currentLatticeStart);
  gaborWaveletCalculator.SetDirections(m_Direction);
  gaborWaveletCalculator.SetLevel(m_gaborLevel);
  gaborWaveletCalculator.SetFMax(m_gaborFMax);
  gaborWaveletCalculator.SetGamma(m_gaborGamma);
  if (m_debug)
  {
    gaborWaveletCalculator.EnableDebugMode();
  }
  gaborWaveletCalculator.Update();
  auto temp = gaborWaveletCalculator.GetOutput();
  for (auto const& f : temp)
  {
    featurevec[f.first] = f.second;
  }
}

template< class TImage >
template< class TShapeImage >
void FeatureExtraction< TImage >::CalculateMorphologic(const typename TImage::Pointer image, const typename TShapeImage::Pointer mask1, const typename TImage::Pointer mask2, std::map<std::string, double>& featurevec)
{
  MorphologicFeatures< TImage, TShapeImage > morphologicCalculator;
  morphologicCalculator.SetInputImage(image);
  morphologicCalculator.SetInputMask(mask2);
  morphologicCalculator.SetMaskShape(mask1);
  morphologicCalculator.SetStartingIndex(m_currentLatticeStart);
  morphologicCalculator.SetRange(m_Range);
  if (m_debug)
  {
    morphologicCalculator.EnableDebugMode();
  }
  if (m_morphologicCalculateFeret)
  {
    morphologicCalculator.EnableCalculateFeretDiameter();
  }
  morphologicCalculator.Update();
  auto temp = morphologicCalculator.GetOutput();
  if (temp.empty())
  {
    WriteErrorFile("");
    return;
  }
  for (auto const& f : temp)
  {
    featurevec[f.first] = f.second;
  }
}

template< class TImage >
void FeatureExtraction< TImage >::CalculateNGLDM(const typename TImage::Pointer itkImage,
  const typename TImage::Pointer maskImage, OffsetVectorPointer offset, std::map<std::string, double>& featurevec)
{
  NGLDMFeatures< TImage > ngldmCalculator;
  ngldmCalculator.SetInputImage(itkImage);
  ngldmCalculator.SetInputMask(maskImage);
  ngldmCalculator.SetNumBins(m_Bins);
  //ngldmCalculator.SetRange(m_Radius); //chebyshev distance delta
  //if (m_histogramBinningType == "Uniform")
  //{
    ngldmCalculator.SetMinimum(m_minimumToConsider);
    ngldmCalculator.SetMaximum(m_maximumToConsider);
  //}
  //else
  //{
  //  ngldmCalculator.SetHistogramTypeEqual();
  //}
  if (m_debug)
  {
    ngldmCalculator.EnableDebugMode();
  }
  ngldmCalculator.SetDistanceMax(GetMaximumDistanceWithinTheDefinedROI(itkImage, maskImage));
  //ngldmCalculator.Update();
  //std::cout << "[DEBUG] FeatureExtraction.hxx::NGLDM::calculator.GetRange() = " << ngldmCalculator.GetRange() << std::endl;

  auto temp = ngldmCalculator.GetOutput();
  for (auto const& f : temp)
  {
    featurevec[f.first] = f.second;
  }
}


template< class TImage >
void FeatureExtraction< TImage >::CalculateNGTDM(const typename TImage::Pointer itkImage,
  const typename TImage::Pointer maskImage, OffsetVectorPointer offset, std::map<std::string, double>& featurevec)
{
  //std::cout << "[DEBUG] FeatureExtraction.hxx::NGTDM" << std::endl;

  NGTDMFeatures< TImage > ngtdmCalculator;
  ngtdmCalculator.SetInputImage(itkImage);
  ngtdmCalculator.SetInputMask(maskImage);
  ngtdmCalculator.SetNumBins(m_Bins);
  ngtdmCalculator.SetRange(m_Radius);
  //if (m_histogramBinningType == "Uniform")
  //{
    ngtdmCalculator.SetMinimum(m_minimumToConsider);
    ngtdmCalculator.SetMaximum(m_maximumToConsider);
  //}
  //else
  //{
  //  ngtdmCalculator.SetHistogramTypeEqual();
  //}
  ngtdmCalculator.SetStartingIndex(m_currentLatticeStart);
  ngtdmCalculator.SetRange(m_Range);
  ngtdmCalculator.Update();
  if (m_debug)
  {
    ngtdmCalculator.EnableDebugMode();
  }

  auto temp = ngtdmCalculator.GetOutput();
  for (auto const& f : temp)
  {
    featurevec[f.first] = f.second;
  }
}


template< class TImage >
void FeatureExtraction< TImage >::CalculateGLSZM(const typename TImage::Pointer itkImage, const typename TImage::Pointer maskImage, OffsetVectorPointer offset, std::map<std::string, double>& featurevec)
{
  GLSZMFeatures< TImage > glszmCalculator;
  glszmCalculator.SetInputImage(itkImage);
  glszmCalculator.SetInputMask(maskImage);
  glszmCalculator.SetNumBins(m_Bins);
  glszmCalculator.SetMaxSize(m_Range);
  //if (m_histogramBinningType == "Uniform")
  //{
    glszmCalculator.SetMinimum(m_minimumToConsider);
    glszmCalculator.SetMaximum(m_maximumToConsider);
  //}
  //else
  //{
  //  glszmCalculator.SetHistogramTypeEqual();
  //}
  glszmCalculator.SetStartingIndex(m_currentLatticeStart);
  glszmCalculator.SetOffsets(offset);
  if (m_debug)
  {
    glszmCalculator.EnableDebugMode();
  }
  glszmCalculator.Update();
  auto temp = glszmCalculator.GetOutput();
  for (auto const& f : temp)
  {
    featurevec[f.first] = f.second;
  }
}


template< class TImage >
void FeatureExtraction< TImage >::CalculateHistogram(const typename TImage::Pointer image, const typename TImage::Pointer mask, std::map< std::string, double >& featurevec, bool latticePatch)
{
  /// histogram calculation from ITK -- for texture feature pipeline 
  typename TImage::PixelType min, max;

  if (m_QuantizationExtent == "Image")
  {
    min = m_minimumToConsider;
    max = m_maximumToConsider;
  }
  else if (m_QuantizationExtent == "ROI")
  {
    min = m_statistics_local.GetMinimum();
    max = m_statistics_local.GetMaximum();
  }


  using HistogramFilterType = itk::Statistics::MaskedImageToHistogramFilter< TImage, TImage >;
  using HistogramMeasurementType = typename HistogramFilterType::HistogramType::MeasurementVectorType;

  HistogramMeasurementType lowerBound(1), upperBound(1);
  upperBound.Fill(max);

  typename HistogramFilterType::HistogramType::SizeType size(1); // this is always grayscale

  auto histogramCalculator = HistogramFilterType::New();
  histogramCalculator->SetInput(image);
  histogramCalculator->SetMaskImage(mask);
  histogramCalculator->SetMaskValue(1);

  float minToConsider;
  if (m_Bins_min == std::numeric_limits<float>::max())
  {
    minToConsider = min;
  }
  else
  {
    minToConsider = m_Bins_min;
  }
  lowerBound.Fill(minToConsider);
  switch (m_histogramBinningType)
  {
  case HistogramBinningType::FixedBinNumber:
  {
    histogramCalculator->SetHistogramBinMinimum(lowerBound);
    histogramCalculator->SetHistogramBinMaximum(upperBound);
    size.Fill(m_Bins);
    break;
  }
  case HistogramBinningType::FixedBinSize:
  {
    histogramCalculator->SetHistogramBinMinimum(lowerBound);
    histogramCalculator->SetHistogramBinMaximum(upperBound);
    float actualBins = static_cast<float>(upperBound[0] - lowerBound[0]) / static_cast<float>(m_Bins); // here, the 'm_Bins' holds the bin size, not the total number of bins
    size.Fill(actualBins);
    break;
  }
  case HistogramBinningType::Equal:
  {
    size.Fill(m_Bins); // no need to set the minim and maximum in this case
    break;
  }
  default:
    break;
  }
  histogramCalculator->SetHistogramSize(size);

  try
  {
    histogramCalculator->Update();
  }
  catch (itk::ExceptionObject & error)
  {
    std::cerr << "Error: " << error << "\n";
  }

  auto histogram = histogramCalculator->GetOutput();

  auto totalFrequency = histogram->GetTotalFrequency();
  std::vector< typename TImage::PixelType > histogramVector/*(histogram->GetTotalFrequency())*/;
  double entropy = 0, uniformity = 0;

  for (size_t bin = 0; bin < histogram->Size(); ++bin)
  {
    auto currentFrequency = histogram->GetFrequency(bin, 0);
    for (size_t f = 0; f < currentFrequency; ++f)
    {
      histogramVector.push_back(bin);
    }
    auto prob = static_cast<double>(currentFrequency) / static_cast<double>(totalFrequency);
    featurevec["Bin-" + std::to_string(bin) + "_Probability"] = prob;
    featurevec["Bin-" + std::to_string(bin) + "_Frequency"] = currentFrequency;
    //featurevec["Bin-" + std::to_string(bin) + "_Max"] = histogram->GetBinMax(0, bin);
    //featurevec["Bin-" + std::to_string(bin) + "_Min"] = histogram->GetBinMin(0, bin);

    if (prob != 0)
    {
      entropy += prob * std::log2(prob);
    }
    uniformity += std::pow(prob, 2); // IBSI 3.4.19
  }

  entropy = -entropy; // IBSI 3.4.18

  cbica::Statistics< typename TImage::PixelType > histogramStatsCalculator;
  histogramStatsCalculator.SetInput(histogramVector);

  TConstIteratorType imageIterator(image, image->GetBufferedRegion()), maskIterator(mask, mask->GetBufferedRegion());
  auto fifthPercentileValue = histogram->Quantile(0, 0.05)/*histogramStatsCalculator.GetNthPercentileElement(5)*/;
  auto ninetyFifthPercentileValue = histogram->Quantile(0, 0.95)/*histogramStatsCalculator.GetNthPercentileElement(95)*/;

  auto fifthPercentileMean = 0.0;
  auto ninetyFifthPercentileMean = 0.0;
  double fifthN = 0.0;
  double ninetyFifthN = 0.0;
  double N = 0.0;

  for (maskIterator.GoToBegin(); !maskIterator.IsAtEnd(); ++maskIterator)
  {
    if (maskIterator.Get() > 0)
    {
      imageIterator.SetIndex(maskIterator.GetIndex());
      auto currentValue = imageIterator.Get();

      N++;

      if (currentValue <= fifthPercentileValue)
      {
        fifthPercentileMean += currentValue;
        fifthN++;
      }
      else if (currentValue >= ninetyFifthPercentileValue)
      {
        ninetyFifthPercentileMean += currentValue;
        ninetyFifthN++;
      }
    }
  }

  // calculate final statistics
  fifthPercentileMean /= fifthN;
  ninetyFifthPercentileMean /= ninetyFifthN;

  featurevec["Minimum"] = histogramStatsCalculator.GetMinimum();
  featurevec["Maximum"] = histogramStatsCalculator.GetMaximum();
  featurevec["Mean"] = histogramStatsCalculator.GetMean();
  featurevec["Sum"] = histogramStatsCalculator.GetSum();
  featurevec["Mode"] = histogramStatsCalculator.GetMode();
  featurevec["Median"] = histogramStatsCalculator.GetMedian();
  featurevec["Variance"] = histogramStatsCalculator.GetVariance();
  featurevec["StandardDeviation"] = histogramStatsCalculator.GetStandardDeviation();
  featurevec["Skewness"] = histogramStatsCalculator.GetSkewness();
  featurevec["Kurtosis"] = histogramStatsCalculator.GetKurtosis();
  featurevec["Range"] = histogramStatsCalculator.GetRange();
  featurevec["Energy"] = histogramStatsCalculator.GetEnergy();
  featurevec["Entropy"] = entropy;
  featurevec["Uniformity"] = uniformity;
  featurevec["RootMeanSquare"] = histogramStatsCalculator.GetRootMeanSquare();
  featurevec["TenthPercentile"] = histogramStatsCalculator.GetNthPercentileElement(10);
  featurevec["NinetiethPercentile"] = histogramStatsCalculator.GetNthPercentileElement(90);
  featurevec["FifthPercentile"] = fifthPercentileValue;
  featurevec["NinetyFifthPercentile"] = ninetyFifthPercentileValue;
  featurevec["FifthPercentileMean"] = fifthPercentileMean;
  featurevec["NinetyFifthPercentileMean"] = ninetyFifthPercentileMean;
  featurevec["TwentyFifthPercentile"] = histogramStatsCalculator.GetNthPercentileElement(25);
  featurevec["SeventyFifthPercentile"] = histogramStatsCalculator.GetNthPercentileElement(75);
  featurevec["InterQuartileRange"] = histogramStatsCalculator.GetInterQuartileRange();
  featurevec["MeanAbsoluteDeviation"] = histogramStatsCalculator.GetMeanAbsoluteDeviation();
  featurevec["RobustMeanAbsoluteDeviation1090"] = histogramStatsCalculator.GetRobustMeanAbsoluteDeviation(10, 90);
  featurevec["MedianAbsoluteDeviation"] = histogramStatsCalculator.GetMedianAbsoluteDeviation();
  featurevec["CoefficientOfVariation"] = histogramStatsCalculator.GetCoefficientOfVariation();
  featurevec["QuartileCoefficientOfVariation"] = histogramStatsCalculator.GetQuartileCoefficientOfDispersion();
}


template< class TImage >
void FeatureExtraction< TImage >::CalculateGLRLM(const typename TImage::Pointer image, const typename TImage::Pointer mask, OffsetVectorPointer offset, std::map<std::string, double>& featurevec, bool latticePatch)
{
  GLRLMFeatures< TImage > glrlmCalculator;
  glrlmCalculator.SetInputImage(image);
  glrlmCalculator.SetInputMask(mask);
  if (m_Bins_min != std::numeric_limits<float>::max())
  {
    m_minimumToConsider = m_Bins_min;
  }
  glrlmCalculator.SetHistogramBinningType(m_histogramBinningType);
  switch (m_histogramBinningType)
  {
  case HistogramBinningType::FixedBinNumber:
  {
    glrlmCalculator.SetMinimum(m_minimumToConsider);
    glrlmCalculator.SetMaximum(m_maximumToConsider);
    glrlmCalculator.SetNumBins(m_Bins);
    break;
  }
  case HistogramBinningType::FixedBinSize:
  {
    glrlmCalculator.SetMinimum(m_minimumToConsider);
    glrlmCalculator.SetMaximum(m_maximumToConsider);
    float actualBins = static_cast<float>(m_maximumToConsider - m_minimumToConsider) / static_cast<float>(m_Bins); // here, the 'm_Bins' holds the bin size, not the total number of bins
    glrlmCalculator.SetNumBins(actualBins);
    break;
  }
  case HistogramBinningType::Equal:
  {
    glrlmCalculator.SetNumBins(m_Bins);
    break;
  }
  default:
    break;
  }
  glrlmCalculator.SetOffsets(offset);
  glrlmCalculator.SetOffsetSelectorType(m_offsetSelect);
  if (m_debug)
  {
    glrlmCalculator.EnableDebugMode();
  }
  if (latticePatch)
  {
    auto maxStep = m_latticeSizeImage[0];
    for (size_t d = 1; d < TImageType::ImageDimension; d++)
    {
      if (maxStep < m_latticeSizeImage[d])
      {
        maxStep = m_latticeSizeImage[d];
      }
    }
    glrlmCalculator.SetDistanceMax(std::sqrt(2) * (maxStep - 1));
  }
  else // get the maximum possible distance within the defied ROI bounding box
  {
    glrlmCalculator.SetDistanceMax(GetMaximumDistanceWithinTheDefinedROI(image, mask));
  }
  glrlmCalculator.Update();
  auto temp = glrlmCalculator.GetOutput();
  for (auto const& f : temp)
  {
    featurevec[f.first] = f.second;
  }
}


template< class TImage >
void FeatureExtraction< TImage >::CalculateGLCM(const typename TImage::Pointer image, const typename TImage::Pointer mask, OffsetVectorPointer offset, std::map<std::string, double>& featurevec, bool latticePatch)
{
  GLCMFeatures< TImage > glcmCalculator;
  glcmCalculator.SetInputImage(image);
  glcmCalculator.SetInputMask(mask);
  if (m_Bins_min != std::numeric_limits<float>::max())
  {
    m_minimumToConsider = m_Bins_min;
  }
  glcmCalculator.SetHistogramBinningType(m_histogramBinningType);
  switch (m_histogramBinningType)
  {
  case HistogramBinningType::FixedBinNumber:
  {
    glcmCalculator.SetMinimum(m_minimumToConsider);
    glcmCalculator.SetMaximum(m_maximumToConsider);
    glcmCalculator.SetNumBins(m_Bins);
    break;
  }
  case HistogramBinningType::FixedBinSize:
  {
    glcmCalculator.SetMinimum(m_minimumToConsider);
    glcmCalculator.SetMaximum(m_maximumToConsider);
    float actualBins = static_cast<float>(m_maximumToConsider - m_minimumToConsider) / static_cast<float>(m_Bins); // here, the 'm_Bins' holds the bin size, not the total number of bins
    glcmCalculator.SetNumBins(actualBins);
    break;
  }
  case HistogramBinningType::Equal:
  {
    glcmCalculator.SetNumBins(m_Bins);
    break;
  }
  default:
    break;
  }
  glcmCalculator.SetOffsets(offset);
  glcmCalculator.SetOffsetSelectorType(m_offsetSelect);
  if (m_debug)
  {
    glcmCalculator.EnableDebugMode();
  }
  glcmCalculator.Update();
  auto temp = glcmCalculator.GetOutput();
  for (auto const& f : temp)
  {
    featurevec[f.first] = f.second;
  }
}


template< class TImage >
void FeatureExtraction< TImage >::CalculateIntensity(std::vector< typename TImage::PixelType >& nonZeroVoxels, std::map< std::string, double >& featurevec, bool latticePatch)
{
  cbica::Statistics< typename TImage::PixelType > statisticsCalculatorToUse;
  if (m_QuantizationExtent == "Image")
  {
    statisticsCalculatorToUse = m_statistics_global[m_currentROIValue];
  }
  else if (m_QuantizationExtent == "ROI")
  {
    statisticsCalculatorToUse = m_statistics_local;
  }

  m_minimumToConsider = statisticsCalculatorToUse.GetMinimum();
  m_maximumToConsider = statisticsCalculatorToUse.GetMaximum();

  featurevec["Minimum"] = statisticsCalculatorToUse.GetMinimum();
  featurevec["Maximum"] = statisticsCalculatorToUse.GetMaximum();
  featurevec["Mean"] = statisticsCalculatorToUse.GetMean();
  featurevec["Sum"] = statisticsCalculatorToUse.GetSum();
  featurevec["Mode"] = statisticsCalculatorToUse.GetMode();
  featurevec["Median"] = statisticsCalculatorToUse.GetMedian();
  featurevec["Variance"] = statisticsCalculatorToUse.GetVariance();
  featurevec["StandardDeviation"] = statisticsCalculatorToUse.GetStandardDeviation();
  featurevec["Skewness"] = statisticsCalculatorToUse.GetSkewness();
  featurevec["Kurtosis"] = statisticsCalculatorToUse.GetKurtosis();
  featurevec["Range"] = statisticsCalculatorToUse.GetRange();
  featurevec["Energy"] = statisticsCalculatorToUse.GetEnergy();
  featurevec["RootMeanSquare"] = statisticsCalculatorToUse.GetRootMeanSquare();
  featurevec["TenthPercentile"] = statisticsCalculatorToUse.GetNthPercentileElement(10);
  featurevec["NinetiethPercentile"] = statisticsCalculatorToUse.GetNthPercentileElement(90);
  featurevec["InterQuartileRange"] = statisticsCalculatorToUse.GetInterQuartileRange();
  featurevec["MeanAbsoluteDeviation"] = statisticsCalculatorToUse.GetMeanAbsoluteDeviation();
  //featurevec["RobustMeanAbsoluteDeviation1090"] = statisticsCalculatorToUse.GetRobustMeanAbsoluteDeviation(10, 90); // commented because of IBSI
  featurevec["MedianAbsoluteDeviation"] = statisticsCalculatorToUse.GetMedianAbsoluteDeviation();
  featurevec["CoefficientOfVariation"] = statisticsCalculatorToUse.GetCoefficientOfVariation();
  featurevec["QuartileCoefficientOfVariation"] = statisticsCalculatorToUse.GetQuartileCoefficientOfDispersion();

  return;
}


template< class TImage >
typename TImage::Pointer FeatureExtraction< TImage >::GetPatchedImage(const typename TImage::Pointer inputImage)
{
  typename TImage::RegionType region(m_currentLatticeStart, m_latticeSizeImage);

  auto extractor = itk::ExtractImageFilter< TImage, TImage >::New();
  extractor->SetInput(inputImage);
  extractor->SetExtractionRegion(region);
  extractor->Update();

  return extractor->GetOutput();
}


template< class TImage >
int FeatureExtraction< TImage >::GetRadiusInImageCoordinates(float radiusInWorldCoordinates)
{
  auto spacing = m_Mask->GetSpacing();
  int minRad = 0;
  for (size_t d = 0; d < TImage::ImageDimension; d++)
  {
    auto temp = radiusInWorldCoordinates / spacing[d];
    auto tempRad = minRad;
    // this is a contingency in cases where the radius has been initialized to be less than the pixel spacing
    // or if it has been initialized as a negative number
    if (temp < 1)
    {
      tempRad = 1;
    }
    else
    {
      tempRad = std::round(temp);
    }

    // we shall keep minimum calculated radius
    if (minRad < tempRad)
    {
      minRad = tempRad;
    }
  }
  return minRad;
 }

template< class TImage >
void FeatureExtraction< TImage >::SetFeatureParam(std::string featureFamily)
{
  auto temp = m_Features.find(featureFamily);
  if (temp != m_Features.end())
  {
    auto featurefamily = temp->second;
    auto parameters = std::get<1>(featurefamily);

    for (auto const& ent1 : parameters) // loop through all parameters for the featureFamily
    {
      auto const& outer_key = ent1.first;
      auto inner_map = ent1.second;

      auto currentValue = inner_map["Value"]; // get the current value of the parameter; other fields are Commends,Default,Range,Type

      if (!outer_key.empty())
      {
        if (outer_key == ParamsString[Dimension])
        {
          m_Dimension = std::atoi(currentValue.c_str());
          if (TImage::ImageDimension == 2) // if 2D image is detected, ensure that the dimension also becomes 2
          {
            m_Dimension = 2;
          }
        }
        else if (outer_key == ParamsString[Axis])
        {
          m_Axis = currentValue;
        }
        else if (outer_key == ParamsString[Radius])
        {
          auto temp = cbica::stringSplit(currentValue, ":");
          if (temp.size() == 1) // single value calculation
          {
            if (currentValue.find(".") != std::string::npos) // this means that the distance is float
            {
              m_Radius_range.push_back(
                GetRadiusInImageCoordinates(
                  std::atof(currentValue.c_str())));
            }
            else
            {
              m_Radius_range.clear();
              m_Radius_range.push_back(std::atoi(currentValue.c_str()));
            }
          }
          else
          {
            // sanity check
            if (temp.size() != 3)
            {
              std::cerr << "Range needs to be in the format 'Min:Step:Max'.\n";
              exit(EXIT_FAILURE);
            }

            std::vector< int > tempRange;

            // check for world coordinates in full set
            bool worldRadDetected = false;
            for (size_t i = 0; i < temp.size(); i++)
            {
              // sanity checks
              if (temp[i].empty())
              {
                std::cerr << "Cannot pass an empty argument to the radius range.\n";
                exit(EXIT_FAILURE);
              }
              else if (temp[i].find("-") != std::string::npos)
              {
                std::cerr << "Radius cannot be negative.\n";
                exit(EXIT_FAILURE);
              }

              if (temp[i].find(".") != std::string::npos) // this means that the distance is float
              {
                worldRadDetected = true;
                break;
              }
            }

            // if a single value is detected in world coordinates, process the entire set the same way
            for (size_t i = 0; i < temp.size(); i++)
            {
              if (worldRadDetected)
              {
                tempRange.push_back(
                  GetRadiusInImageCoordinates(
                    std::atof(temp[i].c_str())));
              }
              else
              {
                tempRange.push_back(std::atoi(temp[i].c_str()));
              }
            }

            int min = tempRange[0],
              max = tempRange[2],
              range = tempRange[1];

            if (min > max) // fail-safe in case someone passes 'Max:Step:Min'
            {
              std::swap(min, max);
            }
            m_Radius_range.clear();
            // populate the full range
            for (int rad = min; rad <= max; rad += range)
            {
              m_Radius_range.push_back(rad);
            }
          } // else-loop end
        }
        else if (outer_key == ParamsString[Neighborhood])
        {
          m_neighborhood = std::atoi(currentValue.c_str());
        }
        else if (outer_key == ParamsString[Bins])
        {
          auto temp = cbica::stringSplit(currentValue, ":");
          if (temp.size() == 1) // single value calculation
          {
            m_Bins_range.clear();
            m_Bins_range.push_back(std::atoi(currentValue.c_str()));
          }
          else
          {
            // sanity check
            if (temp.size() != 3)
            {
              std::cerr << "Range needs to be in the format 'Min:Step:Max'.\n";
              exit(EXIT_FAILURE);
            }
            for (size_t bin = 0; bin < temp.size(); bin++)
            {
              // sanity checks
              if (temp[bin].empty())
              {
                std::cerr << "Cannot pass an empty argument to the bin range.\n";
                exit(EXIT_FAILURE);
              }
              else if (temp[bin].find("-") != std::string::npos)
              {
                std::cerr << "Bins cannot be negative.\n";
                exit(EXIT_FAILURE);
              }
              else if (temp[bin].find(".") != std::string::npos)
              {
                std::cerr << "Bins need to be integer values.\n";
                exit(EXIT_FAILURE);
              }
            }
            int min = std::atoi(temp[0].c_str()),
              max = std::atoi(temp[2].c_str()),
              range = std::atoi(temp[1].c_str());

            if (min > max) // fail-safe in case someone passes 'Max:Step:Min'
            {
              std::swap(min, max);
            }
            m_Bins_range.clear();
            // populate the full range
            for (int bin = min; bin <= max; bin += range)
            {
              m_Bins_range.push_back(bin);
            }
          } // else-loop end
        }
        else if (outer_key == ParamsString[Bins_Min])
        {
          m_Bins_min = std::atof(currentValue.c_str());
          if (m_Bins_min == -666) // the default value coming from the parameter file
          {
            m_Bins_min = std::numeric_limits<float>::max();
          }
        }
        else if (outer_key == ParamsString[Directions])
        {
          if (currentValue.find("|") != std::string::npos)
          {
            m_offsetString = cbica::stringSplit(currentValue, "|");
          }
          else
          {
            m_Direction = std::atoi(currentValue.c_str());
          }
        }
        else if (outer_key == ParamsString[Offset])
        {
          m_offsetSelect = currentValue;
        }
        else if (outer_key == ParamsString[Range])
        {
          m_Range = std::atof(currentValue.c_str());
        }
        else if (outer_key == ParamsString[LatticeWindow])
        {
          m_latticeWindow = std::atof(currentValue.c_str());
        }
        else if (outer_key == ParamsString[LatticeStep])
        {
          m_latticeStep = std::atof(currentValue.c_str());
        }
        else if (outer_key == ParamsString[LatticeBoundary])
        {
          if (currentValue == "FluxNeumann")
          {
            m_fluxNeumannEnabled = true;
          }
          else if (currentValue == "ZeroPadding")
          {
            m_zeroPaddingEnabled = true;
          }
        }
        else if (outer_key == ParamsString[LatticePatchBoundary])
        {
          if (currentValue == "ROI")
          {
            m_patchOnRoiEnabled = true;
          }
          if (currentValue == "None")
          {
            m_patchOnRoiEnabled = true;
            m_patchBoundaryDisregarded = true;
          }
        }

        else if (outer_key == ParamsString[LatticeFullImage])
        {
          if (currentValue == "1")
          {
            m_patchFullImageComputation = true;
          }
        }
        else if (outer_key == ParamsString[GaborFMax])
        {
          m_gaborFMax = std::atof(currentValue.c_str());
        }
        else if (outer_key == ParamsString[GaborGamma])
        {
          m_gaborGamma = std::atof(currentValue.c_str());
        }
        else if (outer_key == ParamsString[GaborLevel])
        {
          m_gaborLevel = std::atoi(currentValue.c_str());
        }
        else if (outer_key == ParamsString[EdgesETA])
        {
          m_edgesETA = std::atof(currentValue.c_str());
        }
        else if (outer_key == ParamsString[EdgesEpsilon])
        {
          m_edgesEpsilon = std::atof(currentValue.c_str());
        }
        else if (outer_key == ParamsString[QuantizationExtent])
        {
          m_QuantizationExtent = currentValue;
        }
        else if (outer_key == ParamsString[QuantizationType])
        {
          auto temp = currentValue;
          std::transform(temp.begin(), temp.end(), temp.begin(), ::tolower);
          if (temp == "fixedbinnumber")
          {
            m_histogramBinningType = HistogramBinningType::FixedBinNumber;
          }
          else if (temp == "fixedbinsize")
          {
            m_histogramBinningType = HistogramBinningType::FixedBinSize;
          }
          else if (temp == "equal")
          {
            m_histogramBinningType = HistogramBinningType::Equal;
          }
          else
          {
            std::cerr << "Unsupported binning type selected; defaulting to FixedBinNumber.\n";
          }
        }
        else if (outer_key == ParamsString[Resampling])
        {
          m_resamplingResolution = std::atof(currentValue.c_str());
        }
        else if (outer_key == ParamsString[ResamplingInterpolator_Image])
        {
          m_resamplingInterpolator_Image = currentValue;
        }
        else if (outer_key == ParamsString[ResamplingInterpolator_Mask])
        {
          m_resamplingInterpolator_Mask = currentValue;
        }
        else if (outer_key == ParamsString[SliceComputation])
        {
          if (currentValue == "1")
          {
            m_SliceComputation = true;
          }
        }
        else if (outer_key == ParamsString[NaNHandling])
        {
          if (currentValue == "0")
          {
            m_keepNaNs = false;
          }
        }
        else if (outer_key == ParamsString[LBPStyle])
        {
          m_LBPStyle = std::atoi(currentValue.c_str());
        }
        else if (outer_key == ParamsString[MorphologicFeret])
        {
          auto temp = std::atoi(currentValue.c_str());
          if (temp == 1)
          {
            m_morphologicCalculateFeret = true;
          }
        }
      }
    } // end of feature parameter map iterator 
  } // sanity check
}


template< class TImage >
void FeatureExtraction< TImage >::SetSelectedROIsAndLabels(std::string roi, std::string roi_labels)
{
  if (!roi.empty() && roi != "all")
  {
    auto tempstr = cbica::stringSplit(roi, "|");
    auto tempstr2 = cbica::stringSplit(roi, ",");
    if (tempstr2.size() > tempstr.size())
    {
      tempstr = tempstr2;
    }
    for (size_t i = 0; i < tempstr.size(); i++)
    {
      m_roi.push_back(std::stoi(tempstr[i]));
      if (roi_labels.empty())
      {
        m_roiLabels.push_back(tempstr[i]);
      }
    }
  }

  if (!roi_labels.empty() && roi_labels != "all")
  {
    m_roiLabels = cbica::stringSplit(roi_labels, "|");
    auto tempstr2 = cbica::stringSplit(roi_labels, ",");
    if (tempstr2.size() > m_roiLabels.size())
    {
      m_roiLabels = tempstr2;
    }
  }

  if (m_debug)
  {
    m_logger.Write("ROI Values: " + roi);
    m_logger.Write("ROI Values Size: " + std::to_string(m_roi.size()));
    m_logger.Write("ROI Labels: " + roi_labels);
    m_logger.Write("ROI Labels Size: " + std::to_string(m_roiLabels.size()));
  }

  if (m_roiLabels.size() != m_roi.size())
  {
    std::string errorString = "The selected roi and the provided roi labels doesn't match";

    auto exeName = cbica::getExecutableName();
    std::transform(exeName.begin(), exeName.end(), exeName.begin(), ::tolower);

    if (exeName.find("captk") != std::string::npos) // TBD this needs a better check than simply "captk", preferably related to qt
    {
      //ShowErrorMessage(errorString);
      //        return featurevec;
    }
    else
    {
      m_logger.WriteError(errorString);
      WriteErrorFile(errorString);
      //exit(EXIT_FAILURE);
      return;
    }
  }
  m_algorithmDone = false;
}


template< class TImage >
void FeatureExtraction< TImage >::SetSelectedROIsAndLabels(std::vector< std::string > roi, std::vector< std::string > roi_labels)
{
  if (m_roiLabels.size() != m_roi.size())
  {
    std::string errorString = "The selected roi and the provided roi labels doesn't match";

    auto exeName = cbica::getExecutableName();
    std::transform(exeName.begin(), exeName.end(), exeName.begin(), ::tolower);

    if (exeName.find("captk") != std::string::npos) // TBD this needs a better check than simply "captk", preferably related to qt
    {
      //ShowErrorMessage(errorString);
      //        return featurevec;
    }
    else
    {
      m_logger.WriteError(errorString);
      WriteErrorFile(errorString);
      //exit(EXIT_FAILURE);
      return;
    }
  }
  else
  {
    for (size_t i = 0; i < roi.size(); i++)
    {
      m_roi.push_back(std::atoi(roi[i].c_str()));
    }
    if (roi.size() == roi_labels.size()) // use roi names, when defined
    {
      for (size_t i = 0; i < roi.size(); i++)
      {
        m_roiLabels.push_back(roi_labels[i]);
      }
    }
    else // otherwise, populate with default string values
    {
      for (size_t i = 0; i < roi.size(); i++)
      {
        m_roiLabels.push_back(roi[i]);
      }
    }
  }
  m_algorithmDone = false;
}


template< class TImage >
void FeatureExtraction< TImage >::SetRequestedFeatures(std::string filename, std::string selected_features)
{
  //parse through feature param file and create a map  of features and params in file
  // make all features available in param file as true 

  std::string parsed;

  //Now check for user input if any feature has been de-selected and make it false
  if (!selected_features.empty())
  {

    std::stringstream feature_stringstream(selected_features);
    while (getline(feature_stringstream, parsed, '|'))
    {
      m_featureParams = FeatureMap(filename, parsed).getFeatureMap();
      m_Features[parsed] = std::make_tuple(true, m_featureParams, parsed, parsed, std::map < std::string, double >());
    }
  }
  else
  {
    for (size_t i = 0; i < FeatureMax; i++)
    {
      m_featureParams = FeatureMap(filename, FeatureFamilyString[i]).getFeatureMap();

      m_Features[FeatureFamilyString[i]] = std::make_tuple(!m_featureParams.empty(), m_featureParams, FeatureFamilyString[i], FeatureFamilyString[i], std::map < std::string, double >());
    }
  }
  m_algorithmDone = false;
}


template< class TImage >
void FeatureExtraction< TImage >::SetRequestedFeatures(std::string filename, std::map<std::string, bool> selected_features)
{
  //parse through feature param file and create a map  of features and params in file
  // make all features available in param file as true 

  std::string parsed;

  //Now check for user input if any feature has been de-selected and make it false
  if (!selected_features.empty())
  {
    for (const auto& f : selected_features)
    {
      m_featureParams = FeatureMap(filename, f.first).getFeatureMap();
      m_Features[f.first] = std::make_tuple(f.second, m_featureParams, f.first, f.first, std::map < std::string, double >());
    }
  }
  m_algorithmDone = false;
}


template< class TImage >
void FeatureExtraction< TImage >::SetRequestedFeatures(std::map< std::string, std::vector< std::map<std::string, std::string> > >  featuresFromUI, std::map<std::string, bool> selected_features)
{
  //parse through feature param file and create a map  of features and params in file
  // make all features available in param file as true 

  std::string parsed;

  for (auto& currentFeature : featuresFromUI) // iterating over the feature families, i.e., GLCM, Intensity, ...
  {
    auto selectedFeatureFlagStruct = selected_features.find(currentFeature.first); // this is to check for user input in UI (if selected or not)

    std::map<std::string, std::map<std::string, std::string>> temp;
    for (size_t j = 0; j < currentFeature.second.size(); j++)
    {

      std::map< std::string, std::string > currentFeature_ParamsAndVals;
      std::string paramName;
      for (auto& currentFeature_Parameter : currentFeature.second[j]) // each parameter within the feature family
      {
        if (currentFeature_Parameter.first != "ParamName")
        {
          currentFeature_ParamsAndVals[currentFeature_Parameter.first] = currentFeature_Parameter.second;
        }
        else
        {
          paramName = currentFeature_Parameter.second;
        }
      }
      temp[paramName] = currentFeature_ParamsAndVals;
    }

    // error check for something weird happening on UI
    if (!currentFeature.second.empty())
    {
      m_Features[currentFeature.first] = std::make_tuple(selectedFeatureFlagStruct->second, // whether the feature is to be extracted or not
        temp, // parameters and respective values
        currentFeature.first, currentFeature.first, // these are the modality and roi label names, which get overwritten with the correct values in the "Update" function
        std::map < std::string, double >());
    }
  }
  m_algorithmDone = false;
}


template< class TImage >
void FeatureExtraction< TImage >::WriteFeatures(const std::string & modality, const std::string & label, const std::string & featureFamily,
  const std::map< std::string, double > & featureList, const std::string & parameters, typename TImage::IndexType centerIndex, bool featureMapWriteForLattice, float weight)
{
  // a precaution
  if (featureList.empty())
  {
    return;
  }
  if (m_outputFile.empty())
  {
    m_outputFile = cbica::createTmpDir() + "/" + m_patientID + "_FEOutput.csv";
    if (m_debug)
    {
      m_logger.WriteError("Output file has not been initialized; saving in '" + m_outputFile + "'");
    }
    SetOutputFilename(m_outputFile);
  }
  //std::ofstream myfile;

  for (auto const& f : featureList)
  {
    auto roiLabelFeatureFamilyFeature = modality + "_" + label + "_" + featureFamily + "_" + f.first;
    bool nanDetected = false;
    if (std::isnan(f.second) || (f.second != f.second))
    {
      if (m_debug)
      {
        m_logger.Write("NAN DETECTED: " + m_patientID + "_" + roiLabelFeatureFamilyFeature);
      }
      nanDetected = true;
      //std::cerr << "NAN DETECTED: " << m_patientID + "_" + roiLabelFeatureFamilyFeature + "_" + "CenterIdx_" + m_centerIndexString << "\n";
    }
    if ((std::isinf(f.second)))
    {
      if (m_debug)
      {
        m_logger.Write("INF DETECTED: " + m_patientID + "_" + roiLabelFeatureFamilyFeature);
      }
      //std::cerr << "INF DETECTED: " << m_patientID + "_" + roiLabelFeatureFamilyFeature + "_" + "CenterIdx_" + m_centerIndexString << "\n";
    }
    if (featureMapWriteForLattice) // if lattice computation has been request AND current ROI has a defined grid node
    {
      auto weightedFeature = f.second * weight;

      if (m_patchBoundaryDisregarded)
      {
        if (weight == 1)
        {
          m_LatticeFeatures[roiLabelFeatureFamilyFeature].push_back(weightedFeature); // calculate the averages
        }
      }
      else
      {
        m_LatticeFeatures[roiLabelFeatureFamilyFeature].push_back(weightedFeature); // calculate the averages
      }

      if (m_writeFeatureMaps)
      {
        auto currentModalityFeatureFamilyFeature = modality + "_" + roiLabelFeatureFamilyFeature;
        auto centerIndexToPopulate = centerIndex;
        for (size_t d = 0; d < TImage::ImageDimension; d++)
        {
          centerIndexToPopulate[d] = std::round(centerIndexToPopulate[d] / m_latticeStepImage[d]);
        }
        // if the corresponding feature map is null, initialize it
        if (m_downscaledFeatureMaps[currentModalityFeatureFamilyFeature].IsNull() || (m_downscaledFeatureMaps[currentModalityFeatureFamilyFeature]->GetBufferedRegion().GetSize()[0] == 0))
        {
          m_downscaledFeatureMaps[currentModalityFeatureFamilyFeature] = cbica::CreateImage< TImage >(m_featureMapBaseImage);
        }
        if (weightedFeature != 0)
        {
          m_downscaledFeatureMaps[currentModalityFeatureFamilyFeature]->SetPixel(centerIndexToPopulate, weightedFeature);
        }
        // write the weighted mask as well
        auto latticeWeightedMaskName = m_patientID + "_" + label + "_Lattice_Weighted_Mask";
        if (m_downscaledFeatureMaps[latticeWeightedMaskName].IsNull() || (m_downscaledFeatureMaps[latticeWeightedMaskName]->GetBufferedRegion().GetSize()[0] == 0))
        {
          m_downscaledFeatureMaps[latticeWeightedMaskName] = cbica::CreateImage< TImage >(m_featureMapBaseImage);
        }
        if (weight != 0)
        {
          m_downscaledFeatureMaps[latticeWeightedMaskName]->SetPixel(centerIndexToPopulate, weight);
        }
      }
    }
    else
    {
      // if NaN has been detec and the user has asked to keep them, write the feature
      // or if the feature has non-NaN value
      if ((nanDetected && m_keepNaNs) || !nanDetected)
      {
        if (m_outputVerticallyConcatenated)
        {
          auto tempParams = parameters;
          if ((featureFamily == "GLCM") || (featureFamily == "GLRLM"))
          {
            if ((f.first.find("Correlation") != std::string::npos) ||
              (f.first.find("LongRunEmphasis") != std::string::npos) ||
              (f.first.find("GreyLevelNonuniformity") != std::string::npos) ||
              (f.first.find("RunLengthNonuniformity") != std::string::npos) ||
              (f.first.find("LongRunLowGreyLevelEmphasis") != std::string::npos) ||
              (f.first.find("LongRunHighGreyLevelEmphasis") != std::string::npos) ||
              (f.first.find("TotalRuns") != std::string::npos) ||
              (f.first.find("RunPercentage") != std::string::npos))
            {
              tempParams += ": IBSI non-compliant";
            }
          }
          m_finalOutputToWrite += m_patientID + "," + modality + "," + label + "," + featureFamily + "," + f.first +
            "," + cbica::to_string_precision(f.second) + "," + tempParams + "\n";
        }

        // for training file, populate these 2 member variables
        m_trainingFile_featureNames += roiLabelFeatureFamilyFeature + ",";
        m_trainingFile_features += cbica::to_string_precision(f.second) + ",";
      }
    }
  }

  //  if (m_outputVerticallyConcatenated)
  //  {
  //#ifndef WIN32
  //    myfile.flush();
  //#endif
  //    myfile.close();
}


template< class TImage >
void FeatureExtraction< TImage >::SetInputImages(std::vector< typename TImage::Pointer > images, std::string modality)
{
  m_inputImages = images;

  if (!modality.empty())
  {
    m_modality = cbica::stringSplit(modality, "|");
  }
  m_algorithmDone = false;
}


template< class TImage >
void FeatureExtraction< TImage >::SetNewLogFile(const std::string & logFile)
{
  m_logger.UseNewFile(logFile);
}


template< class TImage >
void FeatureExtraction< TImage >::SetInputImages(std::vector< typename TImage::Pointer > images, std::vector< std::string > & modality)
{
  if (images.size() != modality.size())
  {
    m_logger.Write("Number of Images and number of modalities are not same; SubjectID: " + m_patientID);
    WriteErrorFile("Number of Images and number of modalities are not same");
    //exit(EXIT_FAILURE);
    return;
  }
  m_inputImages = images;
  m_modality = modality;
  m_algorithmDone = false;
}


template< class TImage >
typename TImage::Pointer FeatureExtraction< TImage >::GetSelectedSlice(typename TImage::Pointer mask, std::string axis)
{
  auto allSlices = GetSelectedSlice(mask);

  auto axis_wrap = axis; // to ensure different cases are handled
  std::transform(axis_wrap.begin(), axis_wrap.end(), axis_wrap.begin(), ::tolower);
  auto axis_wrap_ctr = axis_wrap.c_str();
  if (std::strcmp(axis_wrap_ctr, "x") == 0)
  {
    return allSlices[0];
  }
  else if (std::strcmp(axis_wrap_ctr, "y") == 0)
  {
    return allSlices[1];
  }
  else
  {
    return allSlices[2];
  }
}

template< class TImage >
std::vector< typename TImage::Pointer > FeatureExtraction< TImage >::GetSelectedSlice(typename TImage::Pointer mask)
{
  std::vector< typename TImage::Pointer > maxImageSlices;

  maxImageSlices.resize(TImage::ImageDimension);
  for (size_t dim = 0; dim < TImage::ImageDimension; dim++)
  {
    maxImageSlices[dim] = cbica::CreateImage< TImage >(mask);
  }
  typename TImage::SizeType originalSize = mask->GetLargestPossibleRegion().GetSize();
  auto maxVoxels = originalSize;
  maxVoxels.Fill(0);
  typename TImage::IndexType desiredIndexFinal;

  //std::cout << "[DEBUG] FeatureExtraction.hxx::GetSelectedSlice - started" << std::endl;
  if (originalSize.Dimension == 3)
  {
    for (int dim = 0; dim < TImage::ImageDimension; dim++) // dimension-loop
    {
      maxVoxels[dim] = 0;
      for (size_t i = 0; i < originalSize[dim]; i++)
      {
        //std::cout << "[DEBUG] FeatureExtraction.hxx::GetSelectedSlice - Checking axis [" << dim << "] - slice [" << i << "]" << std::endl;
        typename TImage::RegionType desiredRegion;
        typename TImage::SizeType desiredSize = originalSize;
        desiredSize[dim] = 0;
        typename TImage::IndexType desiredIndex;
        desiredIndex.Fill(0);
        desiredIndex[dim] = i;

        desiredRegion.SetIndex(desiredIndex);
        desiredRegion.SetSize(desiredSize);
        auto extractor = itk::ExtractImageFilter< TImage, ImageType2D >::New();
        extractor->SetInput(mask);
        extractor->SetDirectionCollapseToIdentity();
        extractor->SetExtractionRegion(desiredRegion);
        extractor->Update();

        itk::ImageRegionConstIterator< ImageType2D > iterator(extractor->GetOutput(), extractor->GetOutput()->GetLargestPossibleRegion());
        size_t currentNonZero = 0;
        for (iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator)
        {
          if (iterator.Get() != 0)
          {
            currentNonZero++;
          }
        }

        if (currentNonZero > maxVoxels[dim])
        {
          maxVoxels[dim] = currentNonZero;
          desiredIndexFinal[dim] = i;
        }
      }
    }
  }

  if (originalSize.Dimension == 3)
  {
    itk::ImageRegionIteratorWithIndex< TImage > iteratorExtractor(mask, mask->GetLargestPossibleRegion());

    for (iteratorExtractor.GoToBegin(); !iteratorExtractor.IsAtEnd(); ++iteratorExtractor)
    {
      auto currentMaskValue = iteratorExtractor.Get();
      if (currentMaskValue > 0) // only proceed if non-zero index
      {
        auto idx = iteratorExtractor.GetIndex();
        for (int dim = 0; dim < TImage::ImageDimension; dim++) // loop through all the axes
        {
          if (idx[dim] == desiredIndexFinal[dim]) // if there is a match with the final index along each axis
          {
            maxImageSlices[dim]->SetPixel(idx, currentMaskValue);
          }
        }
      }
    }

  }

  return maxImageSlices;
}


template< class TImage >
void FeatureExtraction< TImage >::Update()
{
  if (!m_algorithmDone)
  {
    auto t1 = std::chrono::high_resolution_clock::now();
    
    m_initializedTimestamp = cbica::getCurrentLocalTimestamp();

    if (m_debug)
    {
      m_logger.Write("Checking mask validity (whether it is empty or not)");
    }

    if (!m_maskValidated)
    {
      if (!m_roi.empty())
      {
        TConstIteratorType maskIt(m_Mask, m_Mask->GetBufferedRegion());
        for (size_t x = 0; x < m_roi.size(); x++)
        {
          maskIt.GoToBegin();
          while (!maskIt.IsAtEnd())
          {
            if (maskIt.Get() == m_roi[x])
              break;
            ++maskIt;
            if (maskIt.IsAtEnd())
            {
              std::cerr << "The ROI for calculation, '" << std::to_string(m_roi[x]) << "' does not exist in the mask; SubjectID: " << m_patientID << "\n";
              WriteErrorFile("The ROI for calculation, '" + std::to_string(m_roi[x]) + "' does not exist in the mask.");
              //exit(EXIT_FAILURE);
              return;
            }
          }
        }
      }
      m_maskValidated = true;
    }

    bool imagesAreOkay = true;
    if (!m_inputImages.empty())
    {
      auto tempSize = m_inputImages[0]->GetBufferedRegion().GetSize();
      for (size_t i = 1; i < m_inputImages.size(); i++)
      {
        auto currentSize = m_inputImages[i]->GetBufferedRegion().GetSize();
        for (size_t d = 0; d < TImage::ImageDimension; d++)
        {
          if (tempSize[d] != currentSize[d])
          {
            m_logger.WriteError("Size Mismatch with images, cannot process.");
            imagesAreOkay = false;
          }
        }
      }
    }

    if (imagesAreOkay)
    {
      // Check if input mask is not null
      auto minMaxCal = itk::MinimumMaximumImageCalculator< TImage >::New();
      minMaxCal->SetImage(m_Mask);
      minMaxCal->ComputeMaximum();
      if (minMaxCal->GetMaximum() == 0)
      {
        std::string errorString = "Mask hasn't been initialized";

        auto exeName = cbica::getExecutableName();
        std::transform(exeName.begin(), exeName.end(), exeName.begin(), ::tolower);
        ////ShowErrorMessage("exeName = " + exeName);

        if (exeName.find("captk") != std::string::npos) // TBD this needs a better check than simply "captk", preferably related to qt
        {
          //  ShowErrorMessage(errorString);
          //  return m_outputFeatureVector;
        }
        else
        {
          m_logger.WriteError(errorString);
          WriteErrorFile(errorString);
          //exit(EXIT_FAILURE);
          return;
        }
      }

      // get the lattice properties, if any
      {
        auto temp = m_Features.find(FeatureFamilyString[Lattice]);
        if (temp != m_Features.end())
        {
          if (std::get<0>(temp->second)) // if the feature family has been selected in the GUI
          {
            m_LatticeComputation = true;
            SetFeatureParam(FeatureFamilyString[Lattice]);
            // all the computation is happening in m_roiConstructor
          }
        }
      }

      // get the quantization properties, if any
      {
        auto temp = m_Features.find(FeatureFamilyString[Generic]);
        if (temp != m_Features.end())
        {
          if (std::get<0>(temp->second)) // if the feature family has been selected in the GUI
          {
            SetFeatureParam(FeatureFamilyString[Generic]);
          }
        }
      }

      if (m_resamplingResolution > 0)
      {
        if (m_debug)
        {
          std::cout << "[DEBUG] Performing resampling of image(s) and mask.\n";
        }
        if (m_resamplingResolution >= 2)
        {
          std::cout << "Feature Extraction is happening on very coarsely sampled inputs, please consider lowering the sampling rate for increased accuracy.\n";
        }
        for (size_t i = 0; i < m_inputImages.size(); i++)
        {
          m_inputImages[i] = cbica::ResampleImage< TImage >(m_inputImages[i], m_resamplingResolution, m_resamplingInterpolator_Image);
          if (m_debug || m_writeIntermediateFiles)
          {
            if (m_debug)
            {
              std::cout << "[DEBUG] Writing resampled image(s) to the output directory.\n";
            }
            cbica::WriteImage< TImage >(m_inputImages[i], cbica::normPath(m_outputPath + "/image_" + m_modality[i] +
              "_resampled_" + std::to_string(m_resamplingResolution) + "-" + m_resamplingInterpolator_Image +
              "_" + m_initializedTimestamp + ".nii.gz"));
          }
        }
        m_Mask = cbica::ResampleImage< TImage >(m_Mask, m_resamplingResolution, m_resamplingInterpolator_Mask);
        if (m_resamplingInterpolator_Mask.find("Nearest") == std::string::npos)
        {
          auto roundingFilter = itk::RoundImageFilter< TImage, TImage >::New();
          roundingFilter->SetInput(m_Mask);
          roundingFilter->Update();
          m_Mask = roundingFilter->GetOutput();
        }
        if (m_debug || m_writeIntermediateFiles)
        {
          if (m_debug)
          {
            std::cout << "[DEBUG] Writing resampled mask to the output directory.\n";
          }
          cbica::WriteImage< TImage >(m_Mask, cbica::normPath(m_outputPath +
            "/mask_resampled_" + std::to_string(m_resamplingResolution) + "-" + m_resamplingInterpolator_Mask + 
            "_" + m_initializedTimestamp + ".nii.gz"));
        }
      }

      // sanity check for the generated ROIs and image to ensure they are in same physical space
      if (!cbica::ImageSanityCheck< TImage >(m_inputImages[0], m_Mask))
      {
        std::cerr << "ERROR: the image and mask are not in the same physical space; please check properties, resample the images are re-try.\n";
        exit(EXIT_FAILURE);
      }

      if (m_debug)
      {
        m_logger.Write("Started Construction of ROIs");
      }

      // set the ROIConstructor up
      m_roiConstructor.SetInputMask(m_Mask);
      m_roiConstructor.SetNewLogFile(m_logger.getLoggingFileName());
      m_roiConstructor.SetSelectedROIsAndLabels(m_roi, m_roiLabels);
      m_roiConstructor.SetLatticeGridStep(m_latticeStep); // if lattice features have not been requested, this gets initialized as zero and no patches are computed
      m_roiConstructor.SetLatticeWindowSize(m_latticeWindow); // if lattice features have not been requested, this gets initialized as zero and no patches are computed
      m_roiConstructor.SetBoundaryCondition(m_fluxNeumannEnabled);
      m_roiConstructor.SetPatchConstructionConditionROI(m_patchOnRoiEnabled);
      m_roiConstructor.SetPatchConstructionConditionNone(m_patchBoundaryDisregarded);
      m_roiConstructor.Update();
      auto allROIs = m_roiConstructor.GetOutput();
      m_LatticeComputation = m_roiConstructor.IsLatticeEnabled(); // checking whether lattice has been enabled or not
      m_latticeStepImage = m_roiConstructor.GetLatticeStepImage();
      auto temp = m_roiConstructor.GetLatticeRadius();

      for (size_t d = 0; d < TImage::ImageDimension; d++)
      {
        m_latticeSizeImage[d] = temp[d] * 2 + 1;
      }

      auto featureMapImageSize = m_Mask->GetBufferedRegion().GetSize(); // size of the (down-sampled) feature map
      auto inputImageSize = m_Mask->GetBufferedRegion().GetSize(); // size of the (down-sampled) feature map
      auto featureMapImageSpacing = m_Mask->GetSpacing(); // spacing of the (down-sampled) feature map
                                //auto featureMapImageSize_world = cbica::GetDistances< TImage >(m_Mask); // size of the (down-sampled) feature map in world coordinates
      
      if (m_LatticeComputation && m_writeFeatureMaps) // check if writing of feature maps has been requested or not
      {
        if (m_debug)
        {
          m_logger.Write("Initializing output base image for feature maps");
        }

        for (size_t i = 0; i < TImage::ImageDimension; i++)
        {
          auto temp = static_cast<float>(featureMapImageSize[i]) / static_cast<float>(m_latticeStepImage[i]);
          if (fmodf(temp, 1) == 0)
          {
            featureMapImageSize[i] = temp + 1;
          }
          else
          {
            featureMapImageSize[i] = std::floor(temp);
          }
          featureMapImageSpacing[i] = inputImageSize[i] / featureMapImageSize[i] * featureMapImageSpacing[i];
        }

        // initialize the feature map output -- this is only used for the lattice feature maps
        m_featureMapBaseImage = cbica::CreateImage< TImage >(m_Mask);
        auto tempDir1 = m_featureMapBaseImage->GetDirection();

        auto resampler = itk::ResampleImageFilter< TImage, TImage >::New();
        resampler->SetInput(m_featureMapBaseImage);
        resampler->SetSize(featureMapImageSize);
        resampler->SetTransform(itk::IdentityTransform< double, TImage::ImageDimension >::New());
        resampler->SetOutputDirection(m_featureMapBaseImage->GetDirection());
        resampler->SetOutputSpacing(featureMapImageSpacing);
        resampler->UpdateLargestPossibleRegion();
        m_featureMapBaseImage = resampler->GetOutput();
      }

      if (m_debug)
      {
        m_logger.Write("Starting feature extraction for every image and every ROI");
      }

      if (m_threads == -1)
      {
        m_threads = omp_get_max_threads();
      }

      size_t j = 0;
      if (!m_LatticeComputation)
      {
        m_threads = 1; // no need for multi-threading if lattice is disabled
      }
      else
      {
        if (!m_patchFullImageComputation && (allROIs.size() > m_roi.size()))
        {
          j += m_roi.size();
        }
      }

      cbica::ProgressBar progressBar(allROIs.size() * m_inputImages.size());
      if (!m_debug)
      {
        std::cout << "Starting computation of selected features.\n";
        progressBar.display();
      }
      //#pragma omp parallel for num_threads(m_threads)
      for (/*j has been initialized earlier*/; j < allROIs.size(); j++)
      {
        bool volumetricFeaturesExtracted = false, morphologicFeaturesExtracted = false;
        for (size_t i = 0; i < m_inputImages.size(); i++)
        {
          auto writeFeatureMapsAndLattice = m_LatticeComputation && allROIs[j].latticeGridPoint;
          // construct the mask and the non-zero image values for each iteration - saves a *lot* of memory
          auto currentMask = cbica::CreateImage< TImage >(m_Mask), currentMask_patch = cbica::CreateImage< TImage >(m_Mask);
          auto currentInputImage = m_inputImages[i], currentInputImage_patch = m_inputImages[i];
          m_currentLatticeCenter = allROIs[j].centerIndex;
          m_currentLatticeStart = m_currentLatticeCenter;
          if (allROIs[j].latticeGridPoint)
          {
            for (size_t d = 0; d < TImage::ImageDimension; d++)
            {
              m_currentLatticeStart[d] = m_currentLatticeStart[d] - std::floor(m_latticeSizeImage[d] / 2); // floor is done because m_latticeSizeImage has a '1' which has been added
            }
          }
          m_currentROIValue = allROIs[j].value;
          m_centerIndexString = "(" + std::to_string(m_currentLatticeCenter[0]);
          for (size_t d = 1; d < TImage::ImageDimension; d++)
          {
            m_centerIndexString += "|" + std::to_string(m_currentLatticeCenter[d]);
          }
          m_centerIndexString += ")";

          if (allROIs[j].latticeGridPoint)
          {
            currentInputImage_patch = GetPatchedImage(m_inputImages[i]);
            currentMask_patch = cbica::CreateImage< TImage >(currentInputImage_patch, 1);
          }
          TIteratorType currentMaskIterator(currentMask_patch, currentMask_patch->GetBufferedRegion()); // LargestRegion is apparently not defined for 2D images
          TConstIteratorType currentImageIterator(m_inputImages[i], m_inputImages[i]->GetBufferedRegion());

          if (allROIs[j].latticeGridPoint)
          {
            auto temp = static_cast<float>(j + i) / static_cast<float>(allROIs.size() + m_inputImages.size());
            if (m_debug)
            {
              m_logger.Write("Percentage done: " + std::to_string(temp * 100));
            }
            //
            //auto totalMemory = static_cast< float >(cbica::getTotalMemory()) / 10e8;
            //auto usedMemory = static_cast< float >(cbica::getCurrentlyUsedMemory());
            //auto freeMem = totalMemory - usedMemory;
            //m_logger.Write("Approximate free memory on machine: '" + std::to_string(freeMem / 10e8) + "' Gb.'");
          }
          else
          {
            if (m_debug)
            {
              m_logger.Write("Calculating Features for modality '" + m_modality[i] + "' and ROI '" + allROIs[j].label + "'");
            }
          }

          m_currentNonZeroImageValues.clear();
          // initialize the mask of the current ROI and get the non-zero pixel/voxel values from the current image
          for (size_t m = 0; m < allROIs[j].nonZeroIndeces.size(); m++)
          {
            currentMaskIterator.SetIndex(allROIs[j].nonZeroIndeces[m]);
            currentMaskIterator.Set(1);

            currentImageIterator.SetIndex(allROIs[j].nonZeroIndeces[m]);
            m_currentNonZeroImageValues.push_back(currentImageIterator.Get());
          }

          // calculate intensity features are always calculated 
          {
            auto tempT1 = std::chrono::high_resolution_clock::now();

            if (m_QuantizationExtent == "Image")
            {
              if (!allROIs[j].latticeGridPoint)
              {
                m_statistics_global[m_currentROIValue].SetInput(m_currentNonZeroImageValues);
              }
            }
            else
            {
              m_statistics_local.SetInput(m_currentNonZeroImageValues);
            }

            auto temp = m_Features.find(FeatureFamilyString[Intensity]);
            SetFeatureParam("Intensity");
            std::get<2>(temp->second) = m_modality[i];
            std::get<3>(temp->second) = allROIs[j].label;
            CalculateIntensity(m_currentNonZeroImageValues, std::get<4>(temp->second), allROIs[j].latticeGridPoint);

            std::string currentFeatureFamily = "Intensity";
            if (!m_Bins_range.empty() && !m_Radius_range.empty())
            {
              currentFeatureFamily += "_Bins-" +
                std::to_string(m_Bins_range[0]) + "_Radius-" + std::to_string(m_Radius_range[0]);
            }

            WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily, std::get<4>(temp->second), "N.A.", m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

            if (m_debug)
            {
              auto tempT2 = std::chrono::high_resolution_clock::now();
              m_logger.Write("Intensity Features for modality '" + m_modality[i] + "' and ROI '" + allROIs[j].label + "' calculated in " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(tempT2 - tempT1).count()) + " milliseconds");
            }
          }

          // get the maximum sliced mask for 3 axes - not needed for 2D
          std::vector< typename TImage::Pointer > currentMask_patch_axisImages;
          if (!writeFeatureMapsAndLattice) // if not lattice image, we do the slice computation
          {
            if ((TImage::ImageDimension == 3) && m_SliceComputation)
            {
              currentMask_patch_axisImages = GetSelectedSlice(currentMask_patch);
            }
          }
          else
          {
            currentMask_patch_axisImages.push_back(currentMask_patch);
          }

          // iterate over the entire feature family enum
          for (size_t f = 1/*Intensity features already calculated*/; f < FeatureMax; ++f)
          {
            SetFeatureParam(FeatureFamilyString[f]);
            switch (f)
            {
              if (m_debug)
              {
                std::cout << "[DEBUG] FeatureExtraction.hxx::SetFeatureParam::FeatureFamilyString[" << f << "]" << std::endl;
              }
            // case Intensity is not needed since it always calculated
            case Histogram:
            {
              auto temp = m_Features.find(FeatureFamilyString[f]);
              if (temp != m_Features.end())
              {
                if (std::get<0>(temp->second))
                {
                  auto tempT1 = std::chrono::high_resolution_clock::now();

                  std::string currentFeatureFamily = std::string(FeatureFamilyString[f]) + "_Bins-" +
                    std::to_string(m_Bins_range[0]);
                  if (!m_Radius_range.empty())
                  {
                    currentFeatureFamily += "_Radius-" + std::to_string(m_Radius_range[0]);
                  }

                  //auto local_map = std::get<1>(temp->second);
                  std::get<2>(temp->second) = m_modality[i];
                  std::get<3>(temp->second) = allROIs[j].label;
                  for (size_t b = 0; b < m_Bins_range.size(); b++)
                  {
                    m_Bins = m_Bins_range[b];
                    auto m_Bins_string = std::to_string(m_Bins);
                    CalculateHistogram(currentInputImage_patch, currentMask_patch, std::get<4>(temp->second), allROIs[j].latticeGridPoint);

                    WriteFeatures(m_modality[i], allROIs[j].label, std::string(currentFeatureFamily) + "_Bins-" + m_Bins_string, std::get<4>(temp->second),
                      "Bins=" + m_Bins_string, m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);
                  }
                  if (m_debug)
                  {
                    auto tempT2 = std::chrono::high_resolution_clock::now();
                    m_logger.Write("Histogram Features for modality '" + m_modality[i] + "' and ROI '" + allROIs[j].label + "' calculated in " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(tempT2 - tempT1).count()) + " milliseconds");
                  }
                }
              }
              break;
            }
            case Morphologic:
            {
              auto temp = m_Features.find(FeatureFamilyString[f]);
              if (temp != m_Features.end())
              {
                if (std::get<0>(temp->second))
                {
                  if (!morphologicFeaturesExtracted) // because this feature is to be extracted per ROI and not per ROI & modality
                  {
                    auto tempT1 = std::chrono::high_resolution_clock::now();

                    std::get<2>(temp->second) = "ALL";
                    std::get<3>(temp->second) = allROIs[j].label;

                    std::string currentFeatureFamily = std::string(FeatureFamilyString[f]);

                    if (!m_Bins_range.empty() && !m_Radius_range.empty())
                    {
                      currentFeatureFamily += "_Bins-" +
                        std::to_string(m_Bins_range[0]) + "_Radius-" + std::to_string(m_Radius_range[0]);
                    }

                    if (TImage::ImageDimension == 3)
                    {
                      CalculateMorphologic<TImage>(currentInputImage_patch, currentMask_patch, currentMask_patch, std::get<4>(temp->second));
                      WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily, std::get<4>(temp->second),
                        "Axis=3D;Dimension=3D", m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

                      if (!writeFeatureMapsAndLattice && m_SliceComputation)
                      {
                        //std::string currentFeatureFamily = FeatureFamilyString[f];
                        CalculateMorphologic<TImage>(currentInputImage_patch, currentMask_patch_axisImages[0], currentMask_patch, std::get<4>(temp->second));
                        WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily + "_X", std::get<4>(temp->second),
                          "Axis=X;Dimension=2D", m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

                        CalculateMorphologic<TImage>(currentInputImage_patch, currentMask_patch_axisImages[1], currentMask_patch, std::get<4>(temp->second));
                        WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily + "_Y", std::get<4>(temp->second),
                          "Axis=Y;Dimension=2D", m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

                        CalculateMorphologic<TImage>(currentInputImage_patch, currentMask_patch_axisImages[2], currentMask_patch, std::get<4>(temp->second));
                        WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily + "_Z", std::get<4>(temp->second),
                          "Axis=Z;Dimension=2D", m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);
                      }
                    }
                    else
                    {
                      CalculateMorphologic<TImage>(currentInputImage_patch, currentMask_patch, currentMask_patch, std::get<4>(temp->second));
                      WriteFeatures(m_modality[i], allROIs[j].label, FeatureFamilyString[f], std::get<4>(temp->second),
                        "Axis=2D;Dimension=2D", m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);
                    }
                    
                    if (m_debug)
                    {
                      auto tempT2 = std::chrono::high_resolution_clock::now();
                      m_logger.Write("Morphologic Features for modality '" + m_modality[i] + "' and ROI '" + allROIs[j].label + "' calculated in " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(tempT2 - tempT1).count()) + " milliseconds");
                    }
                    morphologicFeaturesExtracted = true;
                  }
                }
              }
              break;
            }
            case Volumetric:
            {
              auto temp = m_Features.find(FeatureFamilyString[f]);
              if (temp != m_Features.end())
              {
                if (std::get<0>(temp->second))
                {
                  if (!volumetricFeaturesExtracted) // because this feature is to be extracted per ROI and not per ROI & modality
                  {
                    auto tempT1 = std::chrono::high_resolution_clock::now();

                    std::get<2>(temp->second) = "ALL";
                    std::get<3>(temp->second) = allROIs[j].label;
                    
                    std::string currentFeatureFamily = std::string(FeatureFamilyString[f]);

                    if (!m_Bins_range.empty() && !m_Radius_range.empty())
                    {
                      currentFeatureFamily += "_Bins-" +
                        std::to_string(m_Bins_range[0]) + "_Radius-" + std::to_string(m_Radius_range[0]);
                    }

                    if (TImage::ImageDimension == 3)
                    {
                      CalculateVolumetric<TImage>(currentMask_patch, std::get<4>(temp->second));
                      WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily, std::get<4>(temp->second),
                        "Axis=3D;Dimension=3D", m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

                      if (!writeFeatureMapsAndLattice && m_SliceComputation)
                      {
                        //std::string currentFeatureFamily = FeatureFamilyString[f];
                        CalculateVolumetric<TImage>(currentMask_patch_axisImages[0], std::get<4>(temp->second));
                        WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily + "_X", std::get<4>(temp->second),
                          "Axis=X;Dimension=2D", m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

                        CalculateVolumetric<TImage>(currentMask_patch_axisImages[1], std::get<4>(temp->second));
                        WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily + "_Y", std::get<4>(temp->second),
                          "Axis=Y;Dimension=2D", m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

                        CalculateVolumetric<TImage>(currentMask_patch_axisImages[2], std::get<4>(temp->second));
                        WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily + "_Z", std::get<4>(temp->second),
                          "Axis=Z;Dimension=2D", m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);
                      }
                    }
                    else
                    {
                      CalculateVolumetric<TImage>(currentMask_patch, std::get<4>(temp->second));
                      WriteFeatures(m_modality[i], allROIs[j].label, FeatureFamilyString[f], std::get<4>(temp->second),
                        "Axis=2D;Dimension=2D", m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);
                    }

                    if (m_debug)
                    {
                      auto tempT2 = std::chrono::high_resolution_clock::now();
                      m_logger.Write("Volumetric Features for modality '" + m_modality[i] + "' and ROI '" + allROIs[j].label + "' calculated in " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(tempT2 - tempT1).count()) + " milliseconds");
                    }
                    volumetricFeaturesExtracted = true;
                  }
                }
              }
              break;
            }
            case GLCM:
            {
              auto temp = m_Features.find(FeatureFamilyString[f]);
              if (temp != m_Features.end())
              {
                if (std::get<0>(temp->second))
                {
                  auto tempT1 = std::chrono::high_resolution_clock::now();

                  std::get<2>(temp->second) = m_modality[i];
                  std::get<3>(temp->second) = allROIs[j].label;

                  for (size_t r = 0; r < m_Radius_range.size(); r++)
                  {
                    m_Radius = m_Radius_range[r];
                    auto m_Radius_string = std::to_string(m_Radius);

                    auto offsets = GetOffsetVector(m_Radius, /*m_Direction*/27);
                    //auto offsets_2D = GetOffsetVector(m_Radius, /*m_Direction*/8);

                    for (size_t b = 0; b < m_Bins_range.size(); b++)
                    {
                      m_Bins = m_Bins_range[b];
                      auto m_Bins_string = std::to_string(m_Bins);
                      std::string currentFeatureFamily = std::string(FeatureFamilyString[f]) + "_Bins-" + m_Bins_string + "_Radius-" + m_Radius_string;
                      if (TImage::ImageDimension == 3)
                      {
                        CalculateGLCM(currentInputImage_patch, currentMask_patch, offsets[0], std::get<4>(temp->second), allROIs[j].latticeGridPoint);
                        WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily, std::get<4>(temp->second),
                          "Axis=3D;Dimension=3D;Bins=" + m_Bins_string + ";Directions=" + std::to_string(m_Direction) +
                          ";Radius=" + m_Radius_string + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

                        if (!writeFeatureMapsAndLattice && m_SliceComputation)
                        {
                          CalculateGLCM(currentInputImage_patch, currentMask_patch_axisImages[0], offsets[1], std::get<4>(temp->second), allROIs[j].latticeGridPoint);
                          WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily + "_X", std::get<4>(temp->second),
                            "Axis=X;Dimension=2D;Bins=" + m_Bins_string + ";Directions=" + std::to_string(m_Direction) +
                            ";Radius=" + m_Radius_string + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

                          CalculateGLCM(currentInputImage_patch, currentMask_patch_axisImages[1], offsets[2], std::get<4>(temp->second), allROIs[j].latticeGridPoint);
                          WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily + "_Y", std::get<4>(temp->second),
                            "Axis=Y;Dimension=2D;Bins=" + m_Bins_string + ";Directions=" + std::to_string(m_Direction) +
                            ";Radius=" + m_Radius_string + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

                          CalculateGLCM(currentInputImage_patch, currentMask_patch_axisImages[2], offsets[3], std::get<4>(temp->second), allROIs[j].latticeGridPoint);
                          WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily + "_Z", std::get<4>(temp->second),
                            "Axis=2;Dimension=2D;Bins=" + m_Bins_string + ";Directions=" + std::to_string(m_Direction) +
                            ";Radius=" + m_Radius_string + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);
                        }
                      }
                      else
                      {
                        CalculateGLCM(currentInputImage_patch, currentMask_patch, offsets[0], std::get<4>(temp->second), allROIs[j].latticeGridPoint);
                        WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily + m_Bins_string, std::get<4>(temp->second),
                          "Axis=" + m_Axis + ";Dimension=" + std::to_string(m_Dimension) + ";Bins=" + m_Bins_string + ";Directions=" + std::to_string(m_Direction) +
                          ";Radius=" + m_Radius_string + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice);
                      }
                    } // end bin-loop
                  } // end radius-loop

                  if (m_debug)
                  {
                    auto tempT2 = std::chrono::high_resolution_clock::now();
                    m_logger.Write("GLCM Features for modality '" + m_modality[i] + "' and ROI '" + allROIs[j].label + "' calculated in " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(tempT2 - tempT1).count()) + " milliseconds");
                  }
                }
              }
              break;
            }
            case GLRLM:
            {
              auto temp = m_Features.find(FeatureFamilyString[f]);
              if (temp != m_Features.end())
              {
                if (std::get<0>(temp->second))
                {
                  auto tempT1 = std::chrono::high_resolution_clock::now();

                  std::get<2>(temp->second) = m_modality[i];
                  std::get<3>(temp->second) = allROIs[j].label;

                  for (size_t r = 0; r < m_Radius_range.size(); r++)
                  {
                    m_Radius = m_Radius_range[r];
                    auto m_Radius_string = std::to_string(m_Radius);

                    auto offsets = GetOffsetVector(m_Radius, /*m_Direction*/27);
                    //auto offsets_2D = GetOffsetVector(m_Radius, /*m_Direction*/8);

                    for (size_t b = 0; b < m_Bins_range.size(); b++)
                    {
                      m_Bins = m_Bins_range[b];
                      auto m_Bins_string = std::to_string(m_Bins);
                      std::string currentFeatureFamily = std::string(FeatureFamilyString[f]) + "_Bins-" + m_Bins_string + "_Radius-" + m_Radius_string;
                      if (TImage::ImageDimension == 3)
                      {
                        CalculateGLRLM(currentInputImage_patch, currentMask_patch, offsets[0], std::get<4>(temp->second), allROIs[j].latticeGridPoint);
                        WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily, std::get<4>(temp->second),
                          "Axis=3D;Dimension=3D;Bins=" + std::to_string(m_Bins) + ";Directions=" + std::to_string(m_Direction) +
                          ";Radius=" + m_Radius_string + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

                        if (!writeFeatureMapsAndLattice && m_SliceComputation)
                        {
                          CalculateGLRLM(currentInputImage_patch, currentMask_patch_axisImages[0], offsets[1], std::get<4>(temp->second), allROIs[j].latticeGridPoint);
                          WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily + "_X", std::get<4>(temp->second),
                            "Axis=X;Dimension=2D;Bins=" + std::to_string(m_Bins) + ";Directions=" + std::to_string(m_Direction) +
                            ";Radius=" + m_Radius_string + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

                          CalculateGLRLM(currentInputImage_patch, currentMask_patch_axisImages[1], offsets[2], std::get<4>(temp->second), allROIs[j].latticeGridPoint);
                          WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily + "_Y", std::get<4>(temp->second),
                            "Axis=Y;Dimension=2D;Bins=" + std::to_string(m_Bins) + ";Directions=" + std::to_string(m_Direction) +
                            ";Radius=" + m_Radius_string + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

                          CalculateGLRLM(currentInputImage_patch, currentMask_patch_axisImages[2], offsets[3], std::get<4>(temp->second), allROIs[j].latticeGridPoint);
                          WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily + "_Z", std::get<4>(temp->second),
                            "Axis=2;Dimension=2D;Bins=" + std::to_string(m_Bins) + ";Directions=" + std::to_string(m_Direction) +
                            ";Radius=" + m_Radius_string + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);
                        }
                      }
                      else
                      {
                        CalculateGLRLM(currentInputImage_patch, currentMask_patch, offsets[0], std::get<4>(temp->second), allROIs[j].latticeGridPoint);
                        WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily, std::get<4>(temp->second),
                          "Axis=" + m_Axis + ";Dimension=" + std::to_string(m_Dimension) + ";Bins=" + std::to_string(m_Bins) + ";Directions=" + std::to_string(m_Direction) +
                          ";Radius=" + m_Radius_string + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice);
                      }
                    } // end bin-loop
                  } // end radius-loop

                  if (m_debug)
                  {
                    auto tempT2 = std::chrono::high_resolution_clock::now();
                    m_logger.Write("GLRLM Features for modality '" + m_modality[i] + "' and ROI '" + allROIs[j].label + "' calculated in " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(tempT2 - tempT1).count()) + " milliseconds");
                  }
                }
              }
              break;
            }
            case GLSZM:
            {
              auto temp = m_Features.find(FeatureFamilyString[f]);
              if (temp != m_Features.end())
              {
                if (std::get<0>(temp->second))
                {
                  auto tempT1 = std::chrono::high_resolution_clock::now();

                  std::get<2>(temp->second) = m_modality[i];
                  std::get<3>(temp->second) = allROIs[j].label;

                  for (size_t r = 0; r < m_Radius_range.size(); r++)
                  {
                    m_Radius = m_Radius_range[r];
                    auto m_Radius_string = std::to_string(m_Radius);

                    auto offsets = GetOffsetVector(m_Radius, /*m_Direction*/27);
                    //auto offsets_2D = GetOffsetVector(m_Radius, /*m_Direction*/8);

                    for (size_t b = 0; b < m_Bins_range.size(); b++)
                    {
                      m_Bins = m_Bins_range[b];
                      auto m_Bins_string = std::to_string(m_Bins);
                      std::string currentFeatureFamily = std::string(FeatureFamilyString[f]) + "_Bins-" + m_Bins_string + "_Radius-" + m_Radius_string;
                      if (TImage::ImageDimension == 3)
                      {
                        CalculateGLSZM(currentInputImage_patch, currentMask_patch, offsets[0], std::get<4>(temp->second));
                        WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily, std::get<4>(temp->second),
                          "Axis=3D;Dimension=3D;Bins=" + m_Bins_string + ";Directions=" + std::to_string(m_Direction) +
                          ";Radius=" + m_Radius_string + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

                        if (!writeFeatureMapsAndLattice && m_SliceComputation)
                        {
                          CalculateGLSZM(currentInputImage_patch, currentMask_patch_axisImages[0], offsets[1], std::get<4>(temp->second));
                          WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily + "_X", std::get<4>(temp->second),
                            "Axis=X;Dimension=2D;Bins=" + m_Bins_string + ";Directions=" + std::to_string(m_Direction) +
                            ";Radius=" + m_Radius_string + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

                          CalculateGLSZM(currentInputImage_patch, currentMask_patch_axisImages[1], offsets[2], std::get<4>(temp->second));
                          WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily + "_Y", std::get<4>(temp->second),
                            "Axis=Y;Dimension=2D;Bins=" + m_Bins_string + ";Directions=" + std::to_string(m_Direction) +
                            ";Radius=" + m_Radius_string + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

                          CalculateGLSZM(currentInputImage_patch, currentMask_patch_axisImages[2], offsets[3], std::get<4>(temp->second));
                          WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily + "_Z", std::get<4>(temp->second),
                            "Axis=2;Dimension=2D;Bins=" + m_Bins_string + ";Directions=" + std::to_string(m_Direction) +
                            ";Radius=" + m_Radius_string + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);
                        }
                      }
                      else
                      {
                        CalculateGLSZM(currentInputImage_patch, currentMask_patch, offsets[0], std::get<4>(temp->second));
                        WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily, std::get<4>(temp->second),
                          "Axis=" + m_Axis + ";Dimension=" + std::to_string(m_Dimension) + ";Bins=" + m_Bins_string + ";Directions=" + std::to_string(m_Direction) +
                          ";Radius=" + m_Radius_string + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice);
                      }
                    } // end bin-loop
                  } // end radius-loop

                  if (m_debug)
                  {
                    auto tempT2 = std::chrono::high_resolution_clock::now();
                    m_logger.Write("GLSZM Features for modality '" + m_modality[i] + "' and ROI '" + allROIs[j].label + "' calculated in " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(tempT2 - tempT1).count()) + " milliseconds");
                  }
                }
              }
              break;
            }
            case NGTDM:
            {
              //std::cout << "[DEBUG] FeatureExtraction.hxx::case NGTDM" << std::endl;
              auto temp = m_Features.find(FeatureFamilyString[f]);
              if (temp != m_Features.end())
              {
                if (std::get<0>(temp->second))
                {
                  auto tempT1 = std::chrono::high_resolution_clock::now();

                  std::get<2>(temp->second) = m_modality[i];
                  std::get<3>(temp->second) = allROIs[j].label;

                  for (size_t r = 0; r < m_Radius_range.size(); r++)
                  {
                    m_Radius = m_Radius_range[r];
                    auto m_Radius_string = std::to_string(m_Radius);

                    auto offsets = GetOffsetVector(m_Radius, /*m_Direction*/27);
                    //auto offsets_2D = GetOffsetVector(m_Radius, /*m_Direction*/8);

                    for (size_t b = 0; b < m_Bins_range.size(); b++)
                    {
                      m_Bins = m_Bins_range[b];
                      auto m_Bins_string = std::to_string(m_Bins);
                      std::string currentFeatureFamily = std::string(FeatureFamilyString[f]) + "_Bins-" + m_Bins_string + "_Radius-" + m_Radius_string;
                      if (TImage::ImageDimension == 3)
                      {
                        std::string currentFeatureFamily = FeatureFamilyString[f];
                        CalculateNGTDM(currentInputImage_patch, currentMask_patch, offsets[0], std::get<4>(temp->second));
                        WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily, std::get<4>(temp->second),
                          "Axis=3D;Dimension=3D;Bins=" + m_Bins_string + ";Directions=" + std::to_string(m_Direction) +
                          ";Radius=" + m_Radius_string + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

                        if (!writeFeatureMapsAndLattice && m_SliceComputation)
                        {
                          CalculateNGTDM(currentInputImage_patch, currentMask_patch_axisImages[0], offsets[1], std::get<4>(temp->second));
                          WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily + "_X", std::get<4>(temp->second),
                            "Axis=X;Dimension=2D;Bins=" + m_Bins_string + ";Directions=" + std::to_string(m_Direction) +
                            ";Radius=" + m_Radius_string + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

                          CalculateNGTDM(currentInputImage_patch, currentMask_patch_axisImages[1], offsets[2], std::get<4>(temp->second));
                          WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily + "_Y", std::get<4>(temp->second),
                            "Axis=Y;Dimension=2D;Bins=" + m_Bins_string + ";Directions=" + std::to_string(m_Direction) +
                            ";Radius=" + m_Radius_string + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

                          CalculateNGTDM(currentInputImage_patch, currentMask_patch_axisImages[2], offsets[3], std::get<4>(temp->second));
                          WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily + "_Z", std::get<4>(temp->second),
                            "Axis=2;Dimension=2D;Bins=" + m_Bins_string + ";Directions=" + std::to_string(m_Direction) +
                            ";Radius=" + m_Radius_string + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);
                        }
                      }
                      else
                      {
                        CalculateNGTDM(currentInputImage_patch, currentMask_patch, offsets[0], std::get<4>(temp->second));
                        WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily, std::get<4>(temp->second),
                          "Axis=" + m_Axis + ";Dimension=" + std::to_string(m_Dimension) + ";Bins=" + m_Bins_string + ";Directions=" + std::to_string(m_Direction) +
                          ";Radius=" + m_Radius_string + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice);
                      }
                    } // end bin-loop
                  } // end radius-loop

                  if (m_debug)
                  {
                    auto tempT2 = std::chrono::high_resolution_clock::now();
                    m_logger.Write("NGTDM Features for modality '" + m_modality[i] + "' and ROI '" + allROIs[j].label + "' calculated in " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(tempT2 - tempT1).count()) + " milliseconds");
                  }
                }
              }
              break;
            }
            case NGLDM:
            {
              //std::cout << "[DEBUG] FeatureExtraction.hxx::case NGLDM" << std::endl;
              auto temp = m_Features.find(FeatureFamilyString[f]);
              if (temp != m_Features.end())
              {
                if (std::get<0>(temp->second))
                {
                  auto tempT1 = std::chrono::high_resolution_clock::now();

                  std::get<2>(temp->second) = m_modality[i];
                  std::get<3>(temp->second) = allROIs[j].label;

                  for (size_t r = 0; r < m_Radius_range.size(); r++)
                  {
                    m_Radius = m_Radius_range[r];
                    auto m_Radius_string = std::to_string(m_Radius);

                    auto offsets = GetOffsetVector(m_Radius, /*m_Direction*/27);
                    //auto offsets_2D = GetOffsetVector(m_Radius, /*m_Direction*/8);

                    for (size_t b = 0; b < m_Bins_range.size(); b++)
                    {
                      m_Bins = m_Bins_range[b];
                      auto m_Bins_string = std::to_string(m_Bins);
                      std::string currentFeatureFamily = std::string(FeatureFamilyString[f]) + "_Bins-" + m_Bins_string + "_Radius-" + m_Radius_string;

                      if (TImage::ImageDimension == 3)
                      {
                        std::string currentFeatureFamily = FeatureFamilyString[f];
                        //CalculateNGLDM(currentInputImage_patch, currentMask_patch, offsets[0], std::get<4>(temp->second));
                        //WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily, std::get<4>(temp->second),
                        //  "Axis=3D;Dimension=3D;Bins=" + std::to_string(m_Bins) + ";Directions=" + std::to_string(m_Direction) +
                        //  ";Radius=" + std::to_string(m_Radius) + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

                        //if (!writeFeatureMapsAndLattice && m_SliceComputation)
                        //{
                          // CalculateNGLDM(currentInputImage_patch, currentMask_patch_axisImages[0], offsets_2D, std::get<4>(temp->second));
                          // WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily + "_X", std::get<4>(temp->second),
                          //   "Axis=X;Dimension=2D;Bins=" + std::to_string(m_Bins) + ";Directions=" + std::to_string(m_Direction) +
                          //   ";Radius=" + std::to_string(m_Radius) + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);
                          // 
                          // CalculateNGLDM(currentInputImage_patch, currentMask_patch_axisImages[1], offsets_2D, std::get<4>(temp->second));
                          // WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily + "_Y", std::get<4>(temp->second),
                          //   "Axis=Y;Dimension=2D;Bins=" + std::to_string(m_Bins) + ";Directions=" + std::to_string(m_Direction) +
                          //   ";Radius=" + std::to_string(m_Radius) + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);
                          // 
                          // CalculateNGLDM(currentInputImage_patch, currentMask_patch_axisImages[2], offsets_2D, std::get<4>(temp->second));
                          // WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily + "_Z", std::get<4>(temp->second),
                          //   "Axis=2;Dimension=2D;Bins=" + std::to_string(m_Bins) + ";Directions=" + std::to_string(m_Direction) +
                          //   ";Radius=" + std::to_string(m_Radius) + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);
                        //}
                      }
                      else
                      {
                        if (m_debug)
                        {
                          std::cout << "[DEBUG] NGLDM - Not yet implemented for non-3D" << std::endl;
                        }

                        // CalculateNGTDM(currentInputImage_patch, currentMask_patch, offsets, std::get<4>(temp->second));
                        // WriteFeatures(m_modality[i], allROIs[j].label, FeatureFamilyString[f], std::get<4>(temp->second),
                        //   "Axis=" + m_Axis + ";Dimension=" + std::to_string(m_Dimension) + ";Bins=" + std::to_string(m_Bins) + ";Directions=" + std::to_string(m_Direction) +
                        //   ";Radius=" + std::to_string(m_Radius) + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice);
                      }
                    } // end bin-loop
                  } // end radius-loop

                  if (m_debug)
                  {
                    auto tempT2 = std::chrono::high_resolution_clock::now();
                    m_logger.Write("NGLDM Features for modality '" + m_modality[i] + "' and ROI '" + allROIs[j].label + "' calculated in " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(tempT2 - tempT1).count()) + " milliseconds");
                  }
                }
              }
              break;
            }
            case FractalDimension:
            {
              if (TImage::ImageDimension == 2)
              {
                auto temp = m_Features.find(FeatureFamilyString[f]);
                if (temp != m_Features.end())
                {
                  if (std::get<0>(temp->second))
                  {
                    /// TBD: this is only for lattice points because the computation class has not been optimized (memory usage is going into the Gb range)
                    if (allROIs[j].latticeGridPoint)
                    {
                      auto tempT1 = std::chrono::high_resolution_clock::now();

                      std::get<2>(temp->second) = m_modality[i];
                      std::get<3>(temp->second) = allROIs[j].label;

                      CalculateFractalDimensions(currentInputImage_patch, std::get<4>(temp->second), allROIs[j].latticeGridPoint);

                      WriteFeatures(m_modality[i], allROIs[j].label, FeatureFamilyString[f], std::get<4>(temp->second), "",
                        m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

                      if (m_debug)
                      {
                        auto tempT2 = std::chrono::high_resolution_clock::now();
                        m_logger.Write("Fractal Dimension Features for modality '" + m_modality[i] + "' and ROI '" + allROIs[j].label + "' calculated in " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(tempT2 - tempT1).count()) + " milliseconds");
                      }
                    }
                  }
                }
              }
              break;
            }
            case Gabor:
            {
              if (TImage::ImageDimension == 2)
              {
                auto temp = m_Features.find(FeatureFamilyString[f]);
                if (temp != m_Features.end())
                {
                  if (std::get<0>(temp->second))
                  {
                    /// TBD: this is only for lattice points because the computation class has not been optimized (memory usage is going into the Gb range)
                    if (allROIs[j].latticeGridPoint)
                    {
                      auto tempT1 = std::chrono::high_resolution_clock::now();

                      std::get<2>(temp->second) = m_modality[i];
                      std::get<3>(temp->second) = allROIs[j].label;

                      for (size_t r = 0; r < m_Radius_range.size(); r++)
                      {
                        m_Radius = m_Radius_range[r];
                        auto m_Radius_string = std::to_string(m_Radius);

                        CalculateGaborWavelets(currentInputImage_patch, std::get<4>(temp->second), allROIs[j].latticeGridPoint);

                        WriteFeatures(m_modality[i], allROIs[j].label, std::string(FeatureFamilyString[f]) + "_Radius-" + m_Radius_string, std::get<4>(temp->second),
                          "Radius=" + std::to_string(m_Radius) + ";FMax=" + std::to_string(m_gaborFMax) + ";Gamma=" + std::to_string(m_gaborGamma) +
                          ";Directions=" + m_Radius_string + ";Level=" + std::to_string(m_gaborLevel), m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);
                      } // end radius-loop

                      if (m_debug)
                      {
                        auto tempT2 = std::chrono::high_resolution_clock::now();
                        m_logger.Write("Gabor Wavelet Features for modality '" + m_modality[i] + "' and ROI '" + allROIs[j].label + "' calculated in " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(tempT2 - tempT1).count()) + " milliseconds");
                      }
                    }
                  }
                }
              }
              break;
            }
            case Laws:
            {
              if (TImage::ImageDimension == 2)
              {
                auto temp = m_Features.find(FeatureFamilyString[f]);
                if (temp != m_Features.end())
                {
                  if (std::get<0>(temp->second))
                  {
                    /// TBD: this is only for lattice points because the computation class has not been optimized (memory usage is going into the Gb range)
                    if (allROIs[j].latticeGridPoint)
                    {
                      auto tempT1 = std::chrono::high_resolution_clock::now();

                      std::get<2>(temp->second) = m_modality[i];
                      std::get<3>(temp->second) = allROIs[j].label;

                      CalculateLawsMeasures(currentInputImage_patch, std::get<4>(temp->second));

                      WriteFeatures(m_modality[i], allROIs[j].label, FeatureFamilyString[f], std::get<4>(temp->second), "",
                        m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

                      if (m_debug)
                      {
                        auto tempT2 = std::chrono::high_resolution_clock::now();
                        m_logger.Write("Law Measurement Features for modality '" + m_modality[i] + "' and ROI '" + allROIs[j].label + "' calculated in " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(tempT2 - tempT1).count()) + " milliseconds");
                      }
                    }
                  }
                }
              }
              break;
            }
            case Edges:
            {
              if (TImage::ImageDimension == 2)
              {
                auto temp = m_Features.find(FeatureFamilyString[f]);
                if (temp != m_Features.end())
                {
                  if (std::get<0>(temp->second))
                  {
                    /// TBD: this is only for lattice points because the computation class has not been optimized (memory usage is going into the Gb range)
                    if (allROIs[j].latticeGridPoint)
                    {
                      auto tempT1 = std::chrono::high_resolution_clock::now();

                      std::get<2>(temp->second) = m_modality[i];
                      std::get<3>(temp->second) = allROIs[j].label;

                      for (size_t r = 0; r < m_Radius_range.size(); r++)
                      {
                        m_Radius = m_Radius_range[r];
                        auto m_Radius_string = std::to_string(m_Radius);

                        CalculateEdgeEnhancement(currentInputImage_patch, std::get<4>(temp->second));

                        WriteFeatures(m_modality[i], allROIs[j].label, std::string(FeatureFamilyString[f]) + "_Radius-" + m_Radius_string, std::get<4>(temp->second),
                          "ETA=" + std::to_string(m_edgesETA) + ";Epsilon=" + std::to_string(m_edgesEpsilon) + ";Radius=" + m_Radius_string,
                          m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);
                      } // end radius-loop

                      if (m_debug)
                      {
                        auto tempT2 = std::chrono::high_resolution_clock::now();
                        m_logger.Write("Edge EnhancementFeatures for modality '" + m_modality[i] + "' and ROI '" + allROIs[j].label + "' calculated in " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(tempT2 - tempT1).count()) + " milliseconds");
                      }
                    }
                  }
                }
              }
              break;
            }
            case Power:
            {
              if (TImage::ImageDimension == 2)
              {
                auto temp = m_Features.find(FeatureFamilyString[f]);
                if (temp != m_Features.end())
                {
                  if (std::get<0>(temp->second))
                  {
                    /// TBD: this is only for lattice points because the computation class has not been optimized (memory usage is going into the Gb range)
                    if (allROIs[j].latticeGridPoint)
                    {
                      auto tempT1 = std::chrono::high_resolution_clock::now();

                      std::get<2>(temp->second) = m_modality[i];
                      std::get<3>(temp->second) = allROIs[j].label;

                      CalculatePowerSpectrum(currentInputImage_patch, std::get<4>(temp->second));

                      WriteFeatures(m_modality[i], allROIs[j].label, FeatureFamilyString[f], std::get<4>(temp->second), "",
                        m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

                      if (m_debug)
                      {
                        auto tempT2 = std::chrono::high_resolution_clock::now();
                        m_logger.Write("Power Spectrum Features for modality '" + m_modality[i] + "' and ROI '" + allROIs[j].label + "' calculated in " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(tempT2 - tempT1).count()) + " milliseconds");
                      }
                    }
                  }
                }
              }
              break;
            }
            case LBP:
            {
              auto temp = m_Features.find(FeatureFamilyString[f]);
              if (temp != m_Features.end())
              {
                if (std::get<0>(temp->second))
                {
                  auto tempT1 = std::chrono::high_resolution_clock::now();

                  std::get<2>(temp->second) = m_modality[i];
                  std::get<3>(temp->second) = allROIs[j].label;

                  for (size_t r = 0; r < m_Radius_range.size(); r++)
                  {
                    m_Radius = m_Radius_range[r];
                    auto m_Radius_string = std::to_string(m_Radius);

                    std::string currentFeatureFamily = std::string(FeatureFamilyString[f]) + "_Radius-" 
                      + std::to_string(m_Radius_range[r]);

                    if (!m_Bins_range.empty())
                    {
                      currentFeatureFamily += "_Bins-" +
                        std::to_string(m_Bins_range[0]);
                    }

                    CalculateLBP(currentInputImage_patch, currentMask_patch, std::get<4>(temp->second));

                    WriteFeatures(m_modality[i], allROIs[j].label, currentFeatureFamily, std::get<4>(temp->second),
                      "Neighborhood=" + std::to_string(m_neighborhood) + ";Radius=" + m_Radius_string + ";Style=" + std::to_string(m_LBPStyle), m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);
                  } // end radius-loop

                  if (m_debug)
                  {
                    auto tempT2 = std::chrono::high_resolution_clock::now();
                    m_logger.Write("LBP Features for modality '" + m_modality[i] + "' and ROI '" + allROIs[j].label + "' calculated in " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(tempT2 - tempT1).count()) + " milliseconds");
                  }
                }
              }
              break;
            }
            default: // undefined Feature
              break;
            }
          } // end of feature iteration   

          if (!m_debug)
          {
            ++progressBar;
            progressBar.display();
          }
        } // End of image iteration loop

      } // end of allROIs iteration

      if (!m_debug)
      {
        ++progressBar;
        progressBar.display();
        progressBar.done();
        std::cout << "Finished calculating features, writing the output.\n";
      }

      // calculate 1st order statistics of the different lattice features
      for (auto const& entry : m_LatticeFeatures)
      {
        auto currentPatientModalityROIFeatureFamilyFeature = m_patientID + "_" + entry.first + "_";
        auto currentFeatureValues = entry.second;

        cbica::Statistics< double > firstOrderCalculator(currentFeatureValues);
        auto currentMax = cbica::to_string_precision(firstOrderCalculator.GetMaximum());
        auto currentMin = cbica::to_string_precision(firstOrderCalculator.GetMinimum());
        auto currentVar = cbica::to_string_precision(firstOrderCalculator.GetVariance());
        auto currentStdDev = cbica::to_string_precision(firstOrderCalculator.GetStandardDeviation());
        auto currentSkew = cbica::to_string_precision(firstOrderCalculator.GetSkewness());
        auto currentKurt = cbica::to_string_precision(firstOrderCalculator.GetKurtosis());
        auto currentMean = cbica::to_string_precision(firstOrderCalculator.GetMean());
        auto currentMedian = cbica::to_string_precision(firstOrderCalculator.GetMedian());

        // write the above features into m_output
        if (m_outputVerticallyConcatenated)
        {
          currentPatientModalityROIFeatureFamilyFeature = cbica::stringReplace(currentPatientModalityROIFeatureFamilyFeature, "_", ",");

          m_finalOutputToWrite += currentPatientModalityROIFeatureFamilyFeature + "Max," + currentMax + ",\n";
          m_finalOutputToWrite += currentPatientModalityROIFeatureFamilyFeature + "Min," + currentMin + ",\n";
          m_finalOutputToWrite += currentPatientModalityROIFeatureFamilyFeature + "Variance," + currentVar + ",\n";
          m_finalOutputToWrite += currentPatientModalityROIFeatureFamilyFeature + "StdDev," + currentStdDev + ",\n";
          m_finalOutputToWrite += currentPatientModalityROIFeatureFamilyFeature + "Skewness," + currentSkew + ",\n";
          m_finalOutputToWrite += currentPatientModalityROIFeatureFamilyFeature + "Kurtosis," + currentKurt + ",\n";
          m_finalOutputToWrite += currentPatientModalityROIFeatureFamilyFeature + "Mean," + currentMean + ",\n";
          m_finalOutputToWrite += currentPatientModalityROIFeatureFamilyFeature + "Median," + currentMedian + ",\n";
        }
        else
        {
          currentPatientModalityROIFeatureFamilyFeature = cbica::stringReplace(currentPatientModalityROIFeatureFamilyFeature, m_patientID + "_", "");

          m_trainingFile_featureNames += currentPatientModalityROIFeatureFamilyFeature + "Max" + ",";
          m_trainingFile_features += currentMax + ",";
          m_trainingFile_featureNames += currentPatientModalityROIFeatureFamilyFeature + "Min" + ",";
          m_trainingFile_features += currentMin + ",";
          m_trainingFile_featureNames += currentPatientModalityROIFeatureFamilyFeature + "Variance" + ",";
          m_trainingFile_features += currentVar + ",";
          m_trainingFile_featureNames += currentPatientModalityROIFeatureFamilyFeature + "StdDev" + ",";
          m_trainingFile_features += currentStdDev + ",";
          m_trainingFile_featureNames += currentPatientModalityROIFeatureFamilyFeature + "Skewness" + ",";
          m_trainingFile_features += currentSkew + ",";
          m_trainingFile_featureNames += currentPatientModalityROIFeatureFamilyFeature + "Kurtosis" + ",";
          m_trainingFile_features += currentKurt + ",";
          m_trainingFile_featureNames += currentPatientModalityROIFeatureFamilyFeature + "Mean" + ",";
          m_trainingFile_features += currentMean + ",";
          m_trainingFile_featureNames += currentPatientModalityROIFeatureFamilyFeature + "Median" + ",";
          m_trainingFile_features += currentMedian + ",";
        }
      } // end of lattice features loop

      // write the features for training
      if (m_outputVerticallyConcatenated)
      {
        if (!cbica::isFile(m_outputFile)) // if file is not present & no lattice features are extracted, write the CSV headers 
        {
          m_finalOutputToWrite = "SubjectID,Modality,ROILabel,FeatureFamily,Feature,Value,Parameters\n" + m_finalOutputToWrite;
        }
        std::ofstream myfile;
        myfile.open(m_outputFile, std::ios_base::app);
        // check for locks in a cluster environment
        while (!myfile.is_open())
        {
          cbica::sleep(100);
          myfile.open(m_outputFile, std::ios_base::app);
        }
        myfile << m_finalOutputToWrite;
#ifndef WIN32
        myfile.flush();
#endif
        myfile.close();
      }
      else
      {
        m_trainingFile_featureNames.pop_back(); // since the last character is always a ","
        m_trainingFile_features.pop_back(); // since the last character is always a ","
        auto featureNamesVec = cbica::stringSplit(m_trainingFile_featureNames, ",");
        auto featureVec = cbica::stringSplit(m_trainingFile_features, ",");
        if (featureNamesVec.size() != featureVec.size())
        {
          m_logger.WriteError("Something went wrong and the featureNames (" + std::to_string(featureNamesVec.size()) +
            ") and featureVec (" + std::to_string(featureVec.size()) + ") are not of same size.");
          //exit(EXIT_FAILURE);
          WriteErrorFile("featureNames and featureVector are not of the same size");
          return;
        }

        bool firstRun = true;
        if (cbica::isFile(m_outputFile))
        {
          firstRun = false;

          //m_logger.Write("Training File detected from previous run, writing a new one.");
          //auto temp = cbica::replaceString(cbica::getCurrentLocalDateAndTime(), ",", "");
          //temp = cbica::replaceString(temp, ":", "");
          //m_outputFile = m_outputPath + cbica::getFilenameBase(m_outputFile) + temp + "_forTraining.csv";
        }

        std::ofstream myfile;
        myfile.open(m_outputFile, std::ios_base::app);
        // check for locks in a cluster environment
        while (!myfile.is_open())
        {
          cbica::sleep(100);
          myfile.open(m_outputFile, std::ios_base::app);
        }

        if (firstRun) // write the feature names if this is the first run of the file
        {
          myfile << "SubjectID," + m_trainingFile_featureNames + "\n"; // the assumption is that the filename for m_outputFile is used in the same study, i.e., has same # of features
        }
        myfile.close();
        myfile.open(m_outputFile, std::ios_base::app);
        myfile << m_patientID + "," + m_trainingFile_features + "\n";
#ifndef WIN32
        myfile.flush();
#endif
        myfile.close();
      }

      if (m_writeFeatureMaps && !m_downscaledFeatureMaps.empty())
      {
        m_logger.Write("Writing Feature Maps");

        auto featureMapsPath = m_outputPath + "/featureMaps_" + m_initializedTimestamp;
        cbica::createDir(featureMapsPath);
        for (auto const& entry : m_downscaledFeatureMaps)
        {
          auto currentDownscaledFileName = entry.first;
          cbica::WriteImage< TImage >(entry.second, cbica::normPath(featureMapsPath + "/" + currentDownscaledFileName + ".nii.gz"));
        }
      }

      m_algorithmDone = true;

      auto t2 = std::chrono::high_resolution_clock::now();
      std::cout << "Total computation time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " milliseconds\n";
    } // end imagesAreOkay check
  } // end algorithmDone check
} // end update function
