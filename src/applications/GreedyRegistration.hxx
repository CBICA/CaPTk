/**
\file FeatureExtraction.hxx

\brief Contains the implementations of class FeatureExtraction

*/

#pragma once

#include "FeatureExtraction.h"
#include "itkDOMNodeXMLReader.h"
#include "itkDOMNodeXMLWriter.h"

template< class TImage >
template< class TVolumeImage >
void FeatureExtraction< TImage >::CalculateVolumetric(typename TVolumeImage::Pointer mask, std::map<std::string, double>& featurevec)
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
template< class TShapeImage >
void FeatureExtraction< TImage >::CalculateShape(typename TShapeImage::Pointer mask, std::map<std::string, double>& featurevec)
{
  if (m_Dimension == 3)
  {
    typedef short LabelType;
    typedef itk::Image< LabelType, TShapeImage::ImageDimension > OutputImageType;
    typedef itk::ShapeLabelObject< LabelType, TShapeImage::ImageDimension > ShapeLabelObjectType;
    typedef itk::LabelMap< ShapeLabelObjectType > LabelMapType;

    typedef itk::ConnectedComponentImageFilter < TShapeImage, OutputImageType > ConnectedComponentImageFilterType;
    typedef itk::LabelImageToShapeLabelMapFilter < OutputImageType, LabelMapType > I2LType;
    typename ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
    connected->SetInput(mask);
    connected->Update();

    typename I2LType::Pointer i2l = I2LType::New();
    i2l->SetInput(connected->GetOutput());
    i2l->SetComputePerimeter(true);
    i2l->Update();
    typename LabelMapType::Pointer labelMap = i2l->GetOutput();
    //std::cout << " has " << labelMap->GetNumberOfLabelObjects() << " labels." << std::endl;
    //typename ShapeLabelObjectType::Pointer labelObject = labelMap->GetNthLabelObject(0);

    //std::cout << labelObject->GetPrincipalMoments() << labelObject->GetElongation() <<
    //labelObject->GetPerimeter() << labelObject->GetRoundness() << labelObject->GetFlatness();
    //auto spacing = mask->GetSpacing();
    //double voxvol = spacing[0] * spacing[1] * spacing[2];
    //double volume = labelObject->GetNumberOfPixels() * voxvol;

    int numbers = labelMap->GetNumberOfLabelObjects();
    std::vector <double> eccentricity1;
    std::vector <double> eccentricity2;
    std::vector <double> roundness;
    std::vector <double> flatness;
    std::vector <double> elongation;
    std::vector <double> perimeter;

    for (int i = 1; i < numbers; i++)
    {
      typename ShapeLabelObjectType::Pointer labelObject = labelMap->GetNthLabelObject(i);

      if (labelObject->GetNumberOfPixels() < 20)
        continue;
      auto princComps = labelObject->GetPrincipalMoments();
      eccentricity1.push_back(sqrt(1 - std::pow((princComps[0] / princComps[2]), 2)));
      eccentricity2.push_back(sqrt(1 - std::pow((princComps[0] / princComps[1]), 2)));
      roundness.push_back(labelObject->GetRoundness());
      flatness.push_back(labelObject->GetFlatness());
      elongation.push_back(labelObject->GetElongation());
      perimeter.push_back(labelObject->GetPerimeter());

    }
    double mean_ecc1 = 0;
    double mean_ecc2 = 0;
    double mean_round = 0;
    double mean_flat = 0;
    double mean_perim = 0;
    double mean_elong = 0;
    for (unsigned int i = 0; i < eccentricity1.size(); i++)
    {
      mean_ecc1 = mean_ecc1 + eccentricity1[i];
      mean_ecc2 = mean_ecc2 + eccentricity2[i];
      mean_round = mean_round + roundness[i];
      mean_flat = mean_flat + flatness[i];
      mean_perim = mean_perim + perimeter[i];
      mean_elong = mean_elong + elongation[i];
    }

    featurevec["Eccentricity"] = mean_ecc1 / eccentricity1.size();
    featurevec["Elongation"] = mean_elong / elongation.size();
    featurevec["Perimeter"] = mean_perim;
    featurevec["Roundness"] = mean_round / roundness.size();
    featurevec["Flatness"] = mean_flat / flatness.size();
  }
  if (m_Dimension == 2)
  {
    typedef  short LabelType;
    typedef itk::Image< LabelType, TShapeImage::ImageDimension > OutputImageType;
    typedef itk::ShapeLabelObject< LabelType, TShapeImage::ImageDimension > ShapeLabelObjectType;
    typedef itk::LabelMap< ShapeLabelObjectType > LabelMapType;

    typedef itk::ConnectedComponentImageFilter < TShapeImage, OutputImageType > ConnectedComponentImageFilterType;
    typedef itk::LabelImageToShapeLabelMapFilter < OutputImageType, LabelMapType > I2LType;
    typename ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
    connected->SetInput(mask);
    connected->Update();

    typename I2LType::Pointer i2l = I2LType::New();
    i2l->SetInput(connected->GetOutput());
    i2l->SetComputePerimeter(true);
    i2l->Update();
    typename LabelMapType::Pointer labelMap = i2l->GetOutput();
    typename ShapeLabelObjectType::Pointer labelObject = labelMap->GetNthLabelObject(0);
    //auto spacing = mask->GetSpacing();
    //double pixel_spacing = spacing[0] * spacing[1];
    //double area = labelObject->GetNumberOfPixels() * pixel_spacing;

    int numbers = labelMap->GetNumberOfLabelObjects();
    std::vector <double> eccentricity1;
    std::vector <double> eccentricity2;
    std::vector <double> roundness;
    std::vector <double> flatness;
    std::vector <double> elongation;
    std::vector <double> perimeter;

    for (int i = 1; i < numbers; i++)
    {
      typename ShapeLabelObjectType::Pointer labelObject = labelMap->GetNthLabelObject(i);

      if (labelObject->GetNumberOfPixels() < 20) // this is to ensure really small portions of an ROI do not get picked 
        continue;

      auto princComps = labelObject->GetPrincipalMoments();
      eccentricity1.push_back(sqrt(1 - std::pow((princComps[0] / princComps[2]), 2)));
      eccentricity2.push_back(sqrt(1 - std::pow((princComps[0] / princComps[1]), 2)));
      roundness.push_back(labelObject->GetRoundness());
      flatness.push_back(labelObject->GetFlatness());
      elongation.push_back(labelObject->GetElongation());
      perimeter.push_back(labelObject->GetPerimeter());
    }

    double mean_ecc1 = 0;
    double mean_ecc2 = 0;
    double mean_round = 0;
    double mean_flat = 0;
    double mean_perim = 0;
    double mean_elong = 0;
    for (unsigned int i = 0; i < eccentricity1.size(); i++)
    {
      mean_ecc1 = mean_ecc1 + eccentricity1[i];
      mean_ecc2 = mean_ecc2 + eccentricity2[i];
      mean_round = mean_round + roundness[i];
      mean_flat = mean_flat + flatness[i];
      mean_perim = mean_perim + perimeter[i];
      mean_elong = mean_elong + elongation[i];
    }

    featurevec["Eccentricity"] = mean_ecc1 / eccentricity1.size();
    featurevec["Elongation"] = mean_elong / elongation.size();
    featurevec["Perimeter"] = mean_perim;
    featurevec["Roundness"] = mean_round / roundness.size();
    featurevec["Flatness"] = mean_flat / flatness.size();
  }
}


template< class TImage >
void FeatureExtraction< TImage >::CalculateNGTDM(typename TImage::Pointer itkImage,
  typename TImage::Pointer maskImage, OffsetVector *offset, std::map<std::string, double>& featurevec)
{
  //typedef TImage MaskType;
  typedef itk::Statistics::EnhancedScalarImageToNeighbourhoodGreyLevelDifferenceFeaturesFilter< TImage > FilterType;
  typedef itk::MinimumMaximumImageCalculator< TImage > MinMaxComputerType;
  typedef typename FilterType::NeighbourhoodGreyLevelDifferenceFeaturesFilterType TextureFilterType;

  typename FilterType::Pointer filter = FilterType::New();

  typename FilterType::OffsetVector::Pointer newOffset = FilterType::OffsetVector::New();
  auto oldOffsets = filter->GetOffsets();
  auto oldOffsetsIterator = oldOffsets->Begin();
  while (oldOffsetsIterator != oldOffsets->End())
  {
    bool continueOuterLoop = false;
    typename FilterType::OffsetType offset = oldOffsetsIterator->Value();
    for (size_t i = 0; i < TImage::ImageDimension; ++i)
    {
      if (/*params.m_Direction == i + 2 &&*/ offset[i] != 0)
      {
        continueOuterLoop = true;
      }
    }
    oldOffsetsIterator++;
    if (continueOuterLoop)
      newOffset->push_back(offset);
  }
  filter->SetOffsets(newOffset);


  // All features are required
  typename FilterType::FeatureNameVectorPointer requestedFeatures = FilterType::FeatureNameVector::New();
  requestedFeatures->push_back(TextureFilterType::Coarseness);
  requestedFeatures->push_back(TextureFilterType::Contrast);
  requestedFeatures->push_back(TextureFilterType::Busyness);
  requestedFeatures->push_back(TextureFilterType::Complexity);
  requestedFeatures->push_back(TextureFilterType::Strength);

  typename MinMaxComputerType::Pointer minMaxComputer = MinMaxComputerType::New();
  minMaxComputer->SetImage(itkImage);
  minMaxComputer->Compute();

  auto duplicator_image = itk::ImageDuplicator< TImage >::New();
  auto duplicator_mask = itk::ImageDuplicator< TImage >::New();
  duplicator_image->SetInputImage(itkImage);
  duplicator_image->Update();
  duplicator_mask->SetInputImage(maskImage);
  duplicator_mask->Update();

  auto testimage = duplicator_image->GetOutput();
  auto testmask = duplicator_mask->GetOutput();
  filter->SetInput(testimage);
  filter->SetMaskImage(testmask);
  filter->SetRequestedFeatures(requestedFeatures);
  //int rangeOfPixels = m_Range;
  filter->SetPixelValueMinMax(m_statistics.GetMinimum(), m_statistics.GetMaximum());
  filter->SetNumberOfBinsPerAxis(m_Bins);

  //filter->SetDistanceValueMinMax(0, m_Range);

  filter->Update();


  auto featureMeans = filter->GetFeatureMeans();
  auto featureStd = filter->GetFeatureStandardDeviations();

  //std::ostringstream  ss;
  //ss << rangeOfPixels;
  //std::string strRange = ss.str();
  for (std::size_t i = 0; i < featureMeans->size(); ++i)
  {
    switch (i)
    {
    case TextureFilterType::Coarseness:
      featurevec["Coarsness"] = featureMeans->ElementAt(i);
      break;
    case TextureFilterType::Contrast:
      featurevec["NeighbourContrast"] = featureMeans->ElementAt(i);
      break;
    case TextureFilterType::Busyness:
      featurevec["Busyness"] = featureMeans->ElementAt(i);
      break;
    case TextureFilterType::Complexity:
      featurevec["Complexity"] = featureMeans->ElementAt(i);
      break;
    case TextureFilterType::Strength:
      featurevec["Strength"] = featureMeans->ElementAt(i);
      break;
    default:
      break;
    }
  }
}


template< class TImage >
void FeatureExtraction< TImage >::CalculateGLSZM(typename TImage::Pointer itkImage, typename TImage::Pointer maskImage, std::map<std::string, double>& featurevec)
{
  typedef itk::Statistics::EnhancedScalarImageToSizeZoneFeaturesFilter< TImage > FilterType;
  typedef typename FilterType::SizeZoneFeaturesFilterType TextureFilterType;

  typename FilterType::Pointer filter = FilterType::New();

  typename FilterType::OffsetVector::Pointer newOffset = FilterType::OffsetVector::New();
  auto oldOffsets = filter->GetOffsets();
  auto oldOffsetsIterator = oldOffsets->Begin();
  while (oldOffsetsIterator != oldOffsets->End())
  {
    bool continueOuterLoop = false;
    typename FilterType::OffsetType offset = oldOffsetsIterator->Value();
    for (size_t i = 0; i < TImage::ImageDimension; ++i)
    {
      if (/*params.m_Direction == i + 2 &&*/ offset[i] != 0)
      {
        continueOuterLoop = true;
      }
    }
    /*if (params.m_Direction == 1)
    {
    offset[0] = 0;
    offset[1] = 0;
    offset[2] = 1;
    newOffset->push_back(offset);
    break;
    }*/

    oldOffsetsIterator++;
    if (continueOuterLoop)
      newOffset->push_back(offset);
  }
  filter->SetOffsets(newOffset);


  // All features are required
  typename FilterType::FeatureNameVectorPointer requestedFeatures = FilterType::FeatureNameVector::New();
  requestedFeatures->push_back(TextureFilterType::SmallZoneEmphasis);
  requestedFeatures->push_back(TextureFilterType::LargeZoneEmphasis);
  requestedFeatures->push_back(TextureFilterType::GreyLevelNonuniformity);
  requestedFeatures->push_back(TextureFilterType::SizeZoneNonuniformity);
  requestedFeatures->push_back(TextureFilterType::LowGreyLevelZoneEmphasis);
  requestedFeatures->push_back(TextureFilterType::HighGreyLevelZoneEmphasis);
  requestedFeatures->push_back(TextureFilterType::SmallZoneLowGreyLevelEmphasis);
  requestedFeatures->push_back(TextureFilterType::SmallZoneHighGreyLevelEmphasis);
  requestedFeatures->push_back(TextureFilterType::LargeZoneLowGreyLevelEmphasis);
  requestedFeatures->push_back(TextureFilterType::LargeZoneHighGreyLevelEmphasis);
  requestedFeatures->push_back(TextureFilterType::ZonePercentage);
  requestedFeatures->push_back(TextureFilterType::GreyLevelVariance);
  requestedFeatures->push_back(TextureFilterType::SizeZoneVariance);
  requestedFeatures->push_back(TextureFilterType::ZoneEntropy);

  filter->SetInput(itkImage);
  filter->SetMaskImage(maskImage);
  filter->SetRequestedFeatures(requestedFeatures);
  int rangeOfPixels = m_Range;
  /*if (rangeOfPixels < 2)
  rangeOfPixels = 256;

  if (params.m_UseCtRange)
  {
  filter->SetPixelValueMinMax((TPixel)(-1024.5),(TPixel)(3096.5));
  filter->SetNumberOfBinsPerAxis(3096.5+1024.5);
  } else*/
  {
    filter->SetPixelValueMinMax(m_statistics.GetMinimum(), m_statistics.GetMaximum());
    filter->SetNumberOfBinsPerAxis(m_Bins);
  }

  //filter->SetDistanceValueMinMax(0, rangeOfPixels);

  filter->Update();

  auto featureMeans = filter->GetFeatureMeans();
  auto featureStd = filter->GetFeatureStandardDeviations();

  //filter->Delete();

  std::ostringstream  ss;
  ss << rangeOfPixels;
  std::string strRange = ss.str();
  for (std::size_t i = 0; i < featureMeans->size(); ++i)
  {
    switch (i)
    {
    case TextureFilterType::SmallZoneEmphasis:
      featurevec["ZE"] = featureMeans->ElementAt(i);
      break;
    case TextureFilterType::LargeZoneEmphasis:
      featurevec["LZE"] = featureMeans->ElementAt(i);
      break;
    case TextureFilterType::GreyLevelNonuniformity:
      featurevec["GLN"] = featureMeans->ElementAt(i);
      break;
    case TextureFilterType::SizeZoneNonuniformity:
      featurevec["ZSN"] = featureMeans->ElementAt(i);
      break;
    case TextureFilterType::LowGreyLevelZoneEmphasis:
      featurevec["LGZE"] = featureMeans->ElementAt(i);
      break;
    case TextureFilterType::HighGreyLevelZoneEmphasis:
      featurevec["HGZE"] = featureMeans->ElementAt(i);
      break;
    case TextureFilterType::SmallZoneLowGreyLevelEmphasis:
      featurevec["SZLGE"] = featureMeans->ElementAt(i);
      break;
    case TextureFilterType::SmallZoneHighGreyLevelEmphasis:
      featurevec["SZHGE"] = featureMeans->ElementAt(i);
      break;
    case TextureFilterType::LargeZoneLowGreyLevelEmphasis:
      featurevec["LZLGE"] = featureMeans->ElementAt(i);
      break;
    case TextureFilterType::LargeZoneHighGreyLevelEmphasis:
      featurevec["LZHGE"] = featureMeans->ElementAt(i);
      break;
    case TextureFilterType::ZonePercentage:
      featurevec["ZP"] = featureMeans->ElementAt(i);
      break;
    case TextureFilterType::GreyLevelVariance:
      featurevec["GLV"] = featureMeans->ElementAt(i);
      break;
    case TextureFilterType::SizeZoneVariance:
      featurevec["SZV"] = featureMeans->ElementAt(i);
      break;
    case TextureFilterType::ZoneEntropy:
      featurevec["ZE"] = featureMeans->ElementAt(i);
      break;
    default:
      break;
    }
  }
}


template< class TImage >
void FeatureExtraction< TImage >::CalculateHistogram(typename TImage::Pointer image, typename TImage::Pointer mask, std::map< std::string, double > &featurevec)
{
  const float maxRescaleVal = 1000;
  typedef itk::RescaleIntensityImageFilter< TImage, TImage > RescaleFilterType;
  typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput(image);
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(maxRescaleVal);
  rescaleFilter->Update();

  std::vector<double> intensities;
  itk::ImageRegionConstIterator <TImage> imageIt(rescaleFilter->GetOutput(), rescaleFilter->GetOutput()->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator <TImage> maskIt(mask, mask->GetLargestPossibleRegion());
  imageIt.GoToBegin();
  maskIt.GoToBegin();

  while (!imageIt.IsAtEnd())
  {
    if (maskIt.Get() > 0)
      intensities.push_back(std::round(imageIt.Get()));
    ++imageIt;
    ++maskIt;
  }

  //-------------------------------
  double interval = maxRescaleVal / m_Bins;
  double final_interval = (int)(interval * 100);
  final_interval = (double)final_interval / 100;


  std::vector<double> finalBins;
  std::vector<std::vector<double>> Ranges;
  double current_index = 0;
  for (int i = 0; i < m_Bins; i++)
  {
    std::vector<double> onerange;
    onerange.push_back(current_index);
    current_index = current_index + final_interval;
    onerange.push_back(current_index);
    Ranges.push_back(onerange);

    if (static_cast<int>(Ranges.size()) == m_Bins)
      Ranges[Ranges.size() - 1][1] = maxRescaleVal;
  }
  //toadd the last one
  for (unsigned int j = 0; j < Ranges.size(); j++)
  {
    std::vector<double> onerange = Ranges[j];
    int counter = 0;
    for (unsigned int i = 0; i < intensities.size(); i++)
    {
      if (onerange[0] == 0)
      {
        if (intensities[i] >= onerange[0] && intensities[i] <= onerange[1])
          counter = counter + 1;
      }
      else
      {
        if (intensities[i] > onerange[0] && intensities[i] <= onerange[1])
          counter = counter + 1;
      }
    }
    finalBins.push_back(counter);
  }
  for (unsigned int j = 0; j < finalBins.size(); j++)
  {
    finalBins[j] = (finalBins[j] * 100) / intensities.size();
    featurevec["Bin_" + std::to_string(j)] = finalBins[j];
    featurevec["BinEndIntensity_" + std::to_string(j)] = Ranges[j][1];
  }
}


template< class TImage >
void FeatureExtraction< TImage >::CalculateGLRLM(typename TImage::Pointer image, typename TImage::Pointer mask, OffsetVector *offset, std::map<std::string, double> &featurevec)
{
  using HistogramFrequencyContainerType = itk::Statistics::DenseFrequencyContainer2;
  using RunLengthFilterType = itk::Statistics::ScalarImageToRunLengthFeaturesFilter< TImage, HistogramFrequencyContainerType >;
  using RunLengthMatrixGenerator = typename RunLengthFilterType::RunLengthMatrixFilterType;
  using RunLengthFeatures = typename RunLengthFilterType::RunLengthFeaturesFilterType;
  using InternalRunLengthFeatureName = typename RunLengthFeatures::RunLengthFeatureName;

  typename  RunLengthMatrixGenerator::Pointer matrix_generator = RunLengthMatrixGenerator::New();
  matrix_generator->SetMaskImage(mask);
  matrix_generator->SetInput(image);
  matrix_generator->SetInsidePixelValue(1);
  matrix_generator->SetPixelValueMinMax(m_statistics.GetMinimum(), m_statistics.GetMaximum());
  //removed SetDistanceValueMinMax
  //matrix_generator->SetDistanceValueMinMax(0, m_Range); // TOCHECK - why is this only between 0-4? P
  matrix_generator->SetNumberOfBinsPerAxis(m_Bins); // TOCHECK - needs to be statistically significant

  typename  RunLengthFeatures::Pointer runLengthMatrixCalculator = RunLengthFeatures::New();

  typename  RunLengthFilterType::FeatureNameVectorPointer requestedFeatures = RunLengthFilterType::FeatureNameVector::New();
  requestedFeatures->push_back(RunLengthFeatures::ShortRunEmphasis);
  requestedFeatures->push_back(RunLengthFeatures::LongRunEmphasis);
  requestedFeatures->push_back(RunLengthFeatures::GreyLevelNonuniformity);
  requestedFeatures->push_back(RunLengthFeatures::RunLengthNonuniformity);
  requestedFeatures->push_back(RunLengthFeatures::LowGreyLevelRunEmphasis);
  requestedFeatures->push_back(RunLengthFeatures::HighGreyLevelRunEmphasis);
  requestedFeatures->push_back(RunLengthFeatures::ShortRunLowGreyLevelEmphasis);
  requestedFeatures->push_back(RunLengthFeatures::ShortRunHighGreyLevelEmphasis);
  requestedFeatures->push_back(RunLengthFeatures::LongRunLowGreyLevelEmphasis);
  requestedFeatures->push_back(RunLengthFeatures::LongRunHighGreyLevelEmphasis);

  std::vector <std::string> featurename;
  featurename.push_back("SRE");
  featurename.push_back("LRE");
  featurename.push_back("GLN");
  featurename.push_back("RLN");
  featurename.push_back("LGRE");
  featurename.push_back("HGRE");
  featurename.push_back("SRLGE");
  featurename.push_back("SRHGE");
  featurename.push_back("LRLGE");
  featurename.push_back("LRHGE");

  typename  OffsetVector::ConstIterator offsetIt;
  size_t offsetNum = 0, featureNum = 0;

  if (m_offsetSelect == "Average")
  {
    std::vector<double> tempfeatures/*(requestedFeatures->size(), 0)*/;
    tempfeatures.resize(requestedFeatures->size());
    //matrix_generator->SetOffsets(offset);
    for (offsetIt = offset->Begin(); offsetIt != offset->End(); offsetIt++, offsetNum++)
    {
      matrix_generator->SetOffset(offsetIt.Value());
      matrix_generator->Update();

      runLengthMatrixCalculator->SetInput(matrix_generator->GetOutput());
      runLengthMatrixCalculator->Update();
      typename RunLengthFilterType::FeatureNameVector::ConstIterator fnameIt;
      //std::ostringstream ss;
      featureNum = 0;
      for (fnameIt = requestedFeatures->Begin(); fnameIt != requestedFeatures->End(); ++fnameIt, featureNum++)
      {
        //ss << offsetNum;
        tempfeatures[featureNum] = tempfeatures[featureNum] + runLengthMatrixCalculator->GetFeature((InternalRunLengthFeatureName)fnameIt.Value());
        //featurevec[featurename[featureNum]] = runLengthMatrixCalculator->GetFeature((InternalRunLengthFeatureName)fnameIt.Value());
      }
    }
    for (size_t i = 0; i < tempfeatures.size(); i++)
    {
      tempfeatures[i] = tempfeatures[i] / offset->size();
      featurevec[featurename[i]] = tempfeatures[i];
    }
  }
  else
  {
    for (offsetIt = offset->Begin(); offsetIt != offset->End(); offsetIt++, offsetNum++)
    {
      matrix_generator->SetOffset(offsetIt.Value());
      matrix_generator->Update();

      runLengthMatrixCalculator->SetInput(matrix_generator->GetOutput());
      runLengthMatrixCalculator->Update();
      typename RunLengthFilterType::FeatureNameVector::ConstIterator fnameIt;
      featureNum = 0;
      for (fnameIt = requestedFeatures->Begin(); fnameIt != requestedFeatures->End(); fnameIt++, featureNum++)
      {
        featurevec[featurename[featureNum] + "_Offset_" + std::to_string(offsetNum)] = runLengthMatrixCalculator->GetFeature((InternalRunLengthFeatureName)fnameIt.Value());
      }
    }
  }
}


template< class TImage >
void FeatureExtraction< TImage >::CalculateGLCM(typename TImage::Pointer image, typename TImage::Pointer mask, OffsetVector *offset, std::map<std::string, double> &featurevec)
{
  using Image2CoOccuranceType = itk::Statistics::ScalarImageToCooccurrenceMatrixFilter < TImage >;
  using HistogramType = typename Image2CoOccuranceType::HistogramType;
  using Hist2FeaturesType = itk::Statistics::HistogramToTextureFeaturesFilter< HistogramType >;

  double contrast = 0, correl = 0, ener = 0, entro = 0, homo = 0, clustershade = 0, clusterprominance = 0, autocorr = 0;

  if (m_offsetSelect == "Average")
  {

    for (size_t i = 0; i < offset->size(); i++)
    {
      typename  Image2CoOccuranceType::Pointer glcmGenerator = Image2CoOccuranceType::New();
      glcmGenerator->SetNumberOfBinsPerAxis(m_Bins); //reasonable number of bins
      glcmGenerator->SetPixelValueMinMax(m_statistics.GetMinimum(), m_statistics.GetMaximum());
      glcmGenerator->SetMaskImage(mask);
      glcmGenerator->SetInput(image);
      auto featureCalc = Hist2FeaturesType::New();

      glcmGenerator->SetOffset(offset->at(i));
      glcmGenerator->Update();
      featureCalc->SetInput(glcmGenerator->GetOutput());
      featureCalc->Update();

      contrast += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::Inertia));
      correl += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::Correlation));
      ener += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::Energy));
      homo += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::InverseDifferenceMoment));
      entro += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::Entropy));
      clustershade += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::ClusterShade));
      clusterprominance += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::ClusterProminence));
      autocorr += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::HaralickCorrelation));
    }

    contrast = contrast / offset->size();
    correl = correl / offset->size();
    ener = ener / offset->size();
    homo = homo / offset->size();
    entro = entro / offset->size();
    clusterprominance = clusterprominance / offset->size();
    clustershade = clustershade / offset->size();
    autocorr = autocorr / offset->size();

    featurevec["Correlation"] = correl;
    featurevec["Contrast"] = contrast;
    featurevec["Entropy"] = entro;
    featurevec["Homogenety"] = homo;
    featurevec["Clustershade"] = clustershade;
    featurevec["Clusterprominence"] = clusterprominance;
    featurevec["Autocorrelation"] = autocorr;
    featurevec["Energy"] = ener;
  }
  else
  {
    typename  Image2CoOccuranceType::Pointer glcmGenerator = Image2CoOccuranceType::New();
    glcmGenerator->SetNumberOfBinsPerAxis(m_Bins); //reasonable number of bins
    glcmGenerator->SetPixelValueMinMax(m_statistics.GetMinimum(), m_statistics.GetMaximum());
    glcmGenerator->SetMaskImage(mask);
    glcmGenerator->SetInput(image);
    auto featureCalc = Hist2FeaturesType::New();

    for (size_t i = 0; i < offset->size(); i++)
    {
      glcmGenerator->SetOffset(offset->at(i));
      glcmGenerator->Update();
      featureCalc->SetInput(glcmGenerator->GetOutput());
      featureCalc->Update();

      auto tempStr = std::to_string(i);
      featurevec[std::string("Correlation") + "_Offset_" + tempStr] = featureCalc->GetFeature(Hist2FeaturesType::Correlation);
      featurevec[std::string("Energy") + "_Offset_" + tempStr] = featureCalc->GetFeature(Hist2FeaturesType::Energy);
      featurevec[std::string("Contrast") + "_Offset_" + tempStr] = featureCalc->GetFeature(Hist2FeaturesType::Inertia);
      featurevec[std::string("Entropy") + "_Offset_" + tempStr] = featureCalc->GetFeature(Hist2FeaturesType::Entropy);
      featurevec[std::string("Homogenety") + "_Offset_" + tempStr] = featureCalc->GetFeature(Hist2FeaturesType::InverseDifferenceMoment);
      featurevec[std::string("Clustershade") + "_Offset_" + tempStr] = featureCalc->GetFeature(Hist2FeaturesType::ClusterShade);
      featurevec[std::string("Clusterprominence") + "_Offset_" + tempStr] = featureCalc->GetFeature(Hist2FeaturesType::ClusterProminence);
      featurevec[std::string("Autocorrelation") + "_Offset_" + tempStr] = featureCalc->GetFeature(Hist2FeaturesType::HaralickCorrelation);
    }
  }

  // TODO: Sung to add his GLCM extraction code here
  //featurevec[std::string("Correlation]) + "_Sung"] = 0;
}


template< class TImage >
void FeatureExtraction< TImage >::CalculateIntensity(std::vector< typename TImage::PixelType >& nonZeroVoxels, std::map< std::string, double > &featurevec)
{
  m_statistics.SetInput(nonZeroVoxels);

  featurevec["Minimum"] = m_statistics.GetMinimum();
  featurevec["Maximum"] = m_statistics.GetMaximum();
  featurevec["SumAverage"] = m_statistics.GetMean();
  featurevec["Variance"] = m_statistics.GetVariance();
  featurevec["STD"] = m_statistics.GetStandardDeviation();
  featurevec["Skewness"] = m_statistics.GetSkewness();
  featurevec["Kurtosis"] = m_statistics.GetKurtosis();

  return;
}


template< class TImage >
void FeatureExtraction< TImage >::SetFeatureParam(std::string selected_feature)
{
  auto featurefamily = m_Features.find(selected_feature)->second;
  auto p = std::get<1>(featurefamily);

  if (p.find(ParamsString[Dimension]) != p.end())
  {
    auto  pp = p[ParamsString[Dimension]];
    m_Dimension = stoi(pp["Value"]);
  }
  if (p.find(ParamsString[Axis]) != p.end())
  {
    auto  pp = p[ParamsString[Axis]];
    m_Axis = pp["Value"];
  }
  if (p.find(ParamsString[Bins]) != p.end())
  {
    auto  pp = p[ParamsString[Bins]];
    m_Bins = stoi(pp["Value"]);
  }
  if (p.find(ParamsString[Radius]) != p.end())
  {
    auto  pp = p[ParamsString[Radius]];
    m_Radius = stoi(pp["Value"]);
  }
  if (p.find(ParamsString[Neighborhood]) != p.end())
  {
    auto  pp = p[ParamsString[Neighborhood]];
    m_neighborhood = stoi(pp["Value"]);
  }
  if (p.find(ParamsString[Offset]) != p.end())
  {
    auto  pp = p[ParamsString[Offset]];
    m_offsetSelect = pp["Value"];
  }
  if (p.find(ParamsString[Range]) != p.end())
  {
    auto  pp = p[ParamsString[Range]];
    m_Range = std::stoi(pp["Value"]);
  }
  if (p.find(ParamsString[Directions]) != p.end())
  {
    auto  pp = p[ParamsString[Directions]];
    m_Direction = std::stoi(pp["Value"]);
  }
}


template< class TImage >
void FeatureExtraction< TImage >::SetSelectedLabels(std::string roi, std::string roi_labels)
{
  if (!roi.empty() && roi != "all")
  {
    std::vector<std::string> tempstr = cbica::stringSplit(roi, "|");
    for (size_t i = 0; i< tempstr.size(); i++)
    {
      m_roi.push_back(std::stoi(tempstr[i]));
      if (roi_labels.empty())
      {
        m_labelname.push_back(tempstr[i]);
      }
    }
  }

  if (!roi_labels.empty() && roi_labels != "all")
  {
    m_labelname = cbica::stringSplit(roi_labels, "|");
  }

  if (m_debug)
  {
    std::cout << "[DEBUG] ROI Values: " << roi << "\n";
    std::cout << "[DEBUG] ROI Values Size: " << std::to_string(m_roi.size()) << "\n";
    std::cout << "[DEBUG] ROI Labels: " << roi_labels << "\n";
    std::cout << "[DEBUG] ROI Labels Size: " << std::to_string(m_labelname.size()) << "\n";
  }

  if (m_labelname.size() != m_roi.size())
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
      std::cerr << errorString << "\n";
      exit(EXIT_FAILURE);
    }
  }

}


template< class TImage >
void FeatureExtraction< TImage >::SetSelectedLabels(std::vector< std::string > roi, std::vector< std::string > roi_labels)
{
  m_roi = roi;
  m_labelname = roi_labels;

  if (m_labelname.size() != m_roi.size())
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
      std::cerr << errorString << "\n";
      exit(EXIT_FAILURE);
    }
  }

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
    for (const auto &f : selected_features)
    {
      m_featureParams = FeatureMap(filename, f.first).getFeatureMap();
      m_Features[f.first] = std::make_tuple(f.second, m_featureParams, f.first, f.first, std::map < std::string, double >());
    }
  }
}


template< class TImage >
void FeatureExtraction< TImage >::SetRequestedFeatures(std::map< std::string, std::vector< std::map<std::string, std::string> > >  featuresFromUI, std::map<std::string, bool> selected_features)
{
  //parse through feature param file and create a map  of features and params in file
  // make all features available in param file as true 

  std::string parsed;

  for (auto &currentFeature : featuresFromUI) // iterating over the feature families, i.e., GLCM, Intensity, ...
  {
    auto selectedFeatureFlagStruct = selected_features.find(currentFeature.first); // this is to check for user input in UI (if selected or not)

    std::map<std::string, std::map<std::string, std::string>> temp;
    for (size_t j = 0; j < currentFeature.second.size(); j++)
    {
      std::map< std::string, std::string > currentFeature_ParamsAndVals;
      std::string paramName;
      for (auto &currentFeature_Parameter : currentFeature.second[j]) // each parameter within the feature family
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

    m_Features[currentFeature.first] = std::make_tuple(selectedFeatureFlagStruct->second, // whether the feature is to be extracted or not
      temp, // parameters and respective values
      currentFeature.first, currentFeature.first, // these are the modality and roi label names, which get overwritten with the correct values in the "Update" function
      std::map < std::string, double >());
  }
}


template< class TImage >
void FeatureExtraction< TImage >::WriteFeatureList(std::string Filename, std::vector< FeatureType > featurevec, std::string outputfiletype)
{
  std::string path, base, extension;
  cbica::splitFileName(Filename, path, base, extension);

  // use m_Features to write the parameters and its values here

  if (extension == ".csv")
  {
    std::ofstream myfile;
    myfile.open(Filename, std::ofstream::out | std::ofstream::app);

    for (size_t j = 0; j < featurevec.size(); j++)
    {
      FeatureType tempfeature = featurevec[j];
      for (auto const &ent1 : tempfeature)
      {
        auto temptuple = ent1.second;
        auto features = std::get<4>(temptuple);
        for (auto const &f : features)
        {
          auto str_f_second = std::to_string(f.second);

          std::string toWrite;
          if (!m_patientID.empty())
          {
            toWrite = m_patientID + "_";
          }
          myfile << toWrite + std::get<2>(temptuple) + "_" + ent1.first + "_" + std::get<3>(temptuple) + "_" + f.first + "," + str_f_second + "\n";
#ifndef WIN32
          myfile.flush();
#endif
        }
        //myfile << "\n";
      }
    }
    myfile.close();
  }
  /// TBD: do something for XML
}


template< class TImage >
void FeatureExtraction< TImage >::SetInputImages(std::vector< typename TImage::Pointer > images, std::string modality)
{
  m_inputImages = images;

  if (!modality.empty())
  {
    m_modality = cbica::stringSplit(modality, "|");
  }
}


template< class TImage >
void FeatureExtraction< TImage >::SetInputImages(std::vector< typename TImage::Pointer > images, std::vector< std::string >& modality)
{
  if (images.size() != modality.size())
  {
    std::cout << "Number of Images and number of modalities are not same.\n";
    exit(EXIT_FAILURE);
  }
  m_inputImages = images;
  m_modality = modality;
}


template< class TImage >
ImageType2D::Pointer FeatureExtraction< TImage >::GetSelectedSlice(typename TImage::Pointer mask, std::string axis)
{
  std::vector< typename ImageType2D::Pointer > maxImageSlices;

  maxImageSlices.resize(3);
  typename TImage::SizeType originalSize = mask->GetLargestPossibleRegion().GetSize();
  auto maxVoxels = originalSize;
  maxVoxels.Fill(0);

  if (originalSize.Dimension == 3)
  {
    for (int dim = 0; dim < 3; dim++) // dimension-loop
    {
      maxVoxels[dim] = 0;
      for (size_t i = 0; i < originalSize[dim]; i++)
      {
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

          auto duplicator = itk::ImageDuplicator< ImageType2D >::New();
          duplicator->SetInputImage(extractor->GetOutput());
          duplicator->Update();
          maxImageSlices[dim] = duplicator->GetOutput();
        }
      }
    }
  }
  if (axis == "x")
  {
    return maxImageSlices[0];
  }
  else if (axis == "y")
  {
    return maxImageSlices[1];
  }
  else
  {
    return maxImageSlices[2];
  }
}


template< class TImage >
std::vector< typename FeatureExtraction< TImage >::ROIProperties > FeatureExtraction< TImage >::GetAllROIs(const typename TImage::Pointer inputImage, const typename TImage::Pointer total_mask)
{
  std::vector< ROIProperties > m_roiProperties; // this needs to be defined as a member variable

  TConstIteratorType totalMaskIterator(total_mask, total_mask->GetLargestPossibleRegion());
  TConstIteratorType imageIterator(inputImage, inputImage->GetLargestPossibleRegion());
  totalMaskIterator.GoToBegin();

  if (!m_roi.empty())
  {
    m_roiProperties.resize(m_roi.size());

    while (!totalMaskIterator.IsAtEnd())
    {
      auto currentTotalMaskImageValue = totalMaskIterator.Get();
      if (currentTotalMaskImageValue > 0)
      {
        auto currentIndex = totalMaskIterator.GetIndex();

        for (size_t i = 0; i < m_roi.size(); i++)
        {
          if (currentTotalMaskImageValue == static_cast< typename TImage::PixelType >(m_roi[i]))
          {
            if (m_roiProperties[i].maskImage.IsNull() || (m_roiProperties[i].maskImage->GetLargestPossibleRegion().GetSize()[0] == 0)) // check if this image has been initialized before
            {
              m_roiProperties[i].maskImage = cbica::CreateImage< TImage >(total_mask); // create a blank image taking the properties from total_mask for the current ROI
            }
            TIteratorType maskIterator(m_roiProperties[i].maskImage, m_roiProperties[i].maskImage->GetLargestPossibleRegion());

            maskIterator.SetIndex(currentIndex);
            imageIterator.SetIndex(currentIndex);

            m_roiProperties[i].value = m_roi[i]; // set the value of the current ROI 

            maskIterator.Set(1);
            m_roiProperties[i].nonZeroPixelValues.push_back(imageIterator.Get()); // store the non-zero pixel/voxel values of the image for the current ROI
            m_roiProperties[i].m_nonZeroIndeces.push_back(currentIndex); // store the non-zero pixel/voxel indeces of the image for the current ROI
          }
        }
      }
      ++totalMaskIterator;
    }
  }
  else
  {
    while (!totalMaskIterator.IsAtEnd())
    {
      auto currentTotalMaskImageValue = totalMaskIterator.Get();
      if (currentTotalMaskImageValue > 0)
      {
        if (m_roi.empty())
        {
          m_roi.push_back(currentTotalMaskImageValue);
          m_labelname.push_back(std::to_string(currentTotalMaskImageValue));
        }
        else if (std::find(m_roi.begin(), m_roi.end(), currentTotalMaskImageValue) == m_roi.end()) // if the current mask value is not found in m_roi, add it
        {
          m_roi.push_back(currentTotalMaskImageValue);
          m_labelname.push_back(std::to_string(currentTotalMaskImageValue));
        }
      }
      ++totalMaskIterator;
    }
    return GetAllROIs(inputImage, total_mask);
  }

  return m_roiProperties;
}

template< class TImage >
void FeatureExtraction< TImage >::Update()
{
  //    Check if input mask is not null
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
      std::cerr << errorString << "\n";
      exit(EXIT_FAILURE);
    }
  }


  for (size_t i = 0; i < m_inputImages.size(); i++)
  {
    auto currentImageROIs = GetAllROIs(m_inputImages[i], m_Mask);

    /* Iteration over requested number of roi*/
    bool intensityFeaturesCalculated = false;
    for (size_t j = 0; j < m_labelname.size(); j++)
    {
      /* Intensity features are calculated regardless of selection in the param file*/
      if (!intensityFeaturesCalculated)
      {
        auto temp = m_Features.find("Intensity");
        std::get<2>(temp->second) = m_modality[i];
        std::get<3>(temp->second) = m_labelname[j];
        CalculateIntensity(currentImageROIs[j].nonZeroPixelValues, std::get<4>(temp->second));
        intensityFeaturesCalculated = true;
        std::cout << "Intensity Features for modality '" << m_modality[i] << "' and ROI '" << m_labelname[j] << "' calculated.\n";
      }

      /*Iterate over m_feture data structure */
      for (auto const &mapiterator : m_Features)
      {
        SetFeatureParam(mapiterator.first);
        if (mapiterator.first == "Morphologic")
        {
          auto temp = m_Features.find(mapiterator.first);
          if (std::get<0>(temp->second))
          {
            std::get<2>(temp->second) = m_modality[i];
            std::get<3>(temp->second) = m_labelname[j];

            /* this dimensionality reduction applies only to shape and Volumetric features */
            if (TImage::ImageDimension == 3)
            {
              if (m_Dimension == 2) // extracts slice with maximum area along the specified axis
              {
                auto selected_axis_image = GetSelectedSlice(currentImageROIs[j].maskImage, m_Axis);
                CalculateShape<ImageType2D>(selected_axis_image, std::get<4>(temp->second));
              }
              else
              {
                CalculateShape<TImage>(currentImageROIs[j].maskImage, std::get<4>(temp->second));
              }
            }
            else
            {
              CalculateShape<TImage>(currentImageROIs[j].maskImage, std::get<4>(temp->second));
            }
            std::cout << "Morphologic Features for modality '" << m_modality[i] << "' and ROI '" << m_labelname[j] << "' calculated.\n";
          }
        }

        if (mapiterator.first == "Histogram")
        {
          auto temp = m_Features.find(mapiterator.first);
          if (std::get<0>(temp->second))
          {
            //auto local_map = std::get<1>(temp->second);
            std::get<2>(temp->second) = m_modality[i];
            std::get<3>(temp->second) = m_labelname[j];
            CalculateHistogram(m_inputImages[i], currentImageROIs[j].maskImage, std::get<4>(temp->second));
            std::cout << "Histogram Features for modality '" << m_modality[i] << "' and ROI '" << m_labelname[j] << "' calculated.\n";
          }
        }

        if (mapiterator.first == "Volumetric")
        {
          auto temp = m_Features.find(mapiterator.first);
          if (std::get<0>(temp->second))
          {
            std::get<2>(temp->second) = m_modality[i];
            std::get<3>(temp->second) = m_labelname[j];

            /* this dimensionality reduction applies only to shape and Volumetric features */
            if (TImage::ImageDimension == 3)
            {
              if (m_Dimension == 2)
              {
                ImageType2D::Pointer selected_axis_image = GetSelectedSlice(currentImageROIs[j].maskImage, m_Axis);
                CalculateVolumetric<ImageType2D>(selected_axis_image, std::get<4>(temp->second));
              }
              else
              {
                CalculateVolumetric<TImage>(currentImageROIs[j].maskImage, std::get<4>(temp->second));
              }
            }
            else
            {
              CalculateVolumetric<TImage>(currentImageROIs[j].maskImage, std::get<4>(temp->second));
            }
            std::cout << "Volumetric Features for modality '" << m_modality[i] << "' and ROI '" << m_labelname[j] << "' calculated.\n";
          }
        }

        if (mapiterator.first == "GLRLM")
        {
          auto temp = m_Features.find(mapiterator.first);
          if (std::get<0>(temp->second))
          {
            std::get<2>(temp->second) = m_modality[i];
            std::get<3>(temp->second) = m_labelname[j];

            auto offsets = GetOffsetVector(m_Radius, m_Direction);

            CalculateGLRLM(m_inputImages[i], currentImageROIs[j].maskImage, offsets, std::get<4>(temp->second));
            std::cout << "GLRLM Features for modality '" << m_modality[i] << "' and ROI '" << m_labelname[j] << "' calculated.\n";
          }
        }

        if (mapiterator.first == "GLCM")
        {
          auto temp = m_Features.find(mapiterator.first);
          if (std::get<0>(temp->second))
          {
            std::get<2>(temp->second) = m_modality[i];
            std::get<3>(temp->second) = m_labelname[j];

            auto offsets = GetOffsetVector(m_Radius, m_Direction);

            CalculateGLCM(m_inputImages[i], currentImageROIs[j].maskImage, offsets, std::get<4>(temp->second));
            std::cout << "GLCM Features for modality '" << m_modality[i] << "' and ROI '" << m_labelname[j] << "' calculated.\n";
          }
        }

        if (mapiterator.first == "GLSZM")
        {
          auto temp = m_Features.find(mapiterator.first);
          if (std::get<0>(temp->second))
          {
            std::get<2>(temp->second) = m_modality[i];
            std::get<3>(temp->second) = m_labelname[j];

            auto offsets = GetOffsetVector(m_Radius, m_Direction);

            CalculateGLSZM(m_inputImages[i], currentImageROIs[j].maskImage, std::get<4>(temp->second));
            std::cout << "GLSZM Features for modality '" << m_modality[i] << "' and ROI '" << m_labelname[j] << "' calculated.\n";
          }
        }

        if (mapiterator.first == "NGTDM")
        {
          auto temp = m_Features.find(mapiterator.first);
          if (std::get<0>(temp->second))
          {
            std::get<2>(temp->second) = m_modality[i];
            std::get<3>(temp->second) = m_labelname[j];

            auto offsets = GetOffsetVector(m_Radius, m_Direction);

            CalculateNGTDM(m_inputImages[i], currentImageROIs[j].maskImage, offsets, std::get<4>(temp->second));
            std::cout << "NGTDM Features for modality '" << m_modality[i] << "' and ROI '" << m_labelname[j] << "' calculated.\n";
          }
        }

        if (TImage::ImageDimension == 3)
        {
          if (mapiterator.first == "LBP")
          {
            //auto temp = m_Features.find(mapiterator.first);

            std::set<int> x_value, y_value, z_value;
            for (unsigned int k = 0; k < currentImageROIs[j].m_nonZeroIndeces.size(); k++)
            {
              x_value.insert(currentImageROIs[j].m_nonZeroIndeces[k][0]);
              y_value.insert(currentImageROIs[j].m_nonZeroIndeces[k][1]);
              z_value.insert(currentImageROIs[j].m_nonZeroIndeces[k][2]);
            }

            const ImageTypeFloat3D::SizeType  size = { { x_value.size(), y_value.size(), z_value.size() } };
            const ImageTypeFloat3D::IndexType start = { { *x_value.begin(), *y_value.begin(), *z_value.begin() } };
            //Pad image before calculating LBP

            typename TImage::SizeType lowerExtendRegion;
            lowerExtendRegion.Fill(m_Radius);

            typename TImage::SizeType upperExtendRegion;
            upperExtendRegion.Fill(m_Radius);

            typename TImage::PixelType constantPixel = 0;

            auto padFilter = /*typename*/ itk::ConstantPadImageFilter < TImage, TImage >::New();
            padFilter->SetInput(currentImageROIs[j].maskImage);
            //padFilter->SetPadBound(outputRegion); // Calls SetPadLowerBound(region) and SetPadUpperBound(region)
            padFilter->SetPadLowerBound(lowerExtendRegion);
            padFilter->SetPadUpperBound(upperExtendRegion);
            padFilter->SetConstant(constantPixel);
            padFilter->Update();
            typename TImage::Pointer  lbproi = padFilter->GetOutput();
            //const typename  TImage::SizeType image_size = currentImageROIs[j].maskImage->GetLargestPossibleRegion().GetSize();
            int x1 = -1;  int y1 = -1; int z1 = -1;
            for (unsigned int z = start[2] - m_Radius; z < start[2] + size[2] + m_Radius; z1++, z++)
            {
              y1 = -1;
              for (unsigned int y = start[1] - m_Radius; y < start[1] + m_Radius + size[1]; y1++, y++)
              {
                x1 = -1;
                for (unsigned int x = start[0] - m_Radius; x < start[0] + m_Radius + size[0]; x1++, x++)
                {
                  //TImage::IndexType ind1 = { x, y, z };
                  //TImage::IndexType ind2 = { x1, y1, z1 };
                  //float  pixelValue;
                  //if (x < image_size[0] && y < image_size[1] && z < image_size[2])
                  //{
                  //  pixelValue = currentImageROIs[j].maskImage->GetPixel(ind1);
                  //  lbproi->SetPixel(ind2, pixelValue);
                  //}
                  //else
                  //{
                  //  pixelValue = 0;  // x component
                  //  lbproi->SetPixel(ind2, pixelValue);
                  //}
                }
              }
            }

            //LBPFeatures lbpfeatures;

            //lbpfeatures.calculateLBP<TImage>(m_inputImages[i], currentImageROIs[j].maskImage, lbproi, m_Radius, m_neighborhood, m_modality[i], std::to_string(m_roi[j]), std::get<4>(temp->second));
            std::cout << "LBP Features for modality '" << m_modality[i] << "' and ROI '" << m_labelname[j] << "' calculated.\n";
          }
        }

      } // end of feature iteration      
      m_outputFeatureVector.push_back(m_Features);
    }/* End of iteration over roi*/


  } /*End of image iteration loop*/

  if (m_debug)
  {
    std::cout << "[DEBUG] Writing Features.\n";
  }
  WriteFeatureList(m_File, m_outputFeatureVector/*, "csv"*/);
}