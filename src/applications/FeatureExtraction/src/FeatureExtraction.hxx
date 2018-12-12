/**
\file FeatureExtraction.hxx

\brief Contains the implementations of class FeatureExtraction

*/

#pragma once

//#include "FeatureExtraction.h"
#include "itkDOMNodeXMLReader.h"
#include "itkDOMNodeXMLWriter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkBinaryImageToLabelMapFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkMaskedImageToHistogramFilter.h"
#include "itkRegionOfInterestImageFilter.h"

#include "itkEnhancedHistogramToRunLengthFeaturesFilter.h"
#include "itkEnhancedScalarImageToRunLengthFeaturesFilter.h"

//#include "itkOpenCVImageBridge.h"

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
//#include "FractalBoxCount_template.h"

#ifdef _OPENMP
#include <omp.h>
#endif

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
void FeatureExtraction< TImage >::CalculateFractalDimensions(const typename TImage::Pointer itkImage, std::map< std::string, double > &featurevec, bool latticePatch)
{
  FractalBoxCount< TImage > fractalDimensionCalculator;
  fractalDimensionCalculator.SetInputImage(itkImage);
  fractalDimensionCalculator.SetRadius(m_Radius);
  fractalDimensionCalculator.SetLatticePointStatus(latticePatch);
  fractalDimensionCalculator.SetStartingIndex(m_currentLatticeStart);
  fractalDimensionCalculator.Update();
  auto temp = fractalDimensionCalculator.GetOutput();
  for (auto const &f : temp)
  {
    featurevec[f.first] = f.second;
  }
}


template< class TImage >
void FeatureExtraction< TImage >::CalculateLawsMeasures(const typename TImage::Pointer itkImage, std::map< std::string, double > &featurevec)
{
  LawsMeasures< TImage > lawsMeasuresCalculator;
  lawsMeasuresCalculator.SetInputImage(itkImage);
  lawsMeasuresCalculator.SetStartingIndex(m_currentLatticeStart);
  lawsMeasuresCalculator.Update();
  auto temp = lawsMeasuresCalculator.GetOutput();
  for (auto const &f : temp)
  {
    featurevec[f.first] = f.second;
  }
}


template< class TImage >
void FeatureExtraction< TImage >::CalculateEdgeEnhancement(const typename TImage::Pointer itkImage, std::map< std::string, double > &featurevec)
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
  edgeEnhancementCalculator.Update();
  auto temp = edgeEnhancementCalculator.GetOutput();
  for (auto const &f : temp)
  {
    featurevec[f.first] = f.second;
  }
}


template< class TImage >
void FeatureExtraction< TImage >::CalculateLBP(const typename TImage::Pointer itkImage, const typename TImage::Pointer mask, std::map< std::string, double > &featurevec)
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
  lbpCalculator.SetNeighbors(m_neighborhood);
  lbpCalculator.SetLBPStyle(m_LBPStyle);
  lbpCalculator.Update();
  auto temp = lbpCalculator.GetOutput();
  for (auto const &f : temp)
  {
    featurevec[f.first] = f.second;
  }
}


template< class TImage >
void FeatureExtraction< TImage >::CalculatePowerSpectrum(const typename TImage::Pointer itkImage, std::map< std::string, double > &featurevec)
{
  PowerSpectrum< TImage > powerSpectrumCalculator;
  powerSpectrumCalculator.SetInputImage(itkImage);
  powerSpectrumCalculator.SetCenter(m_centerIndexString);
  powerSpectrumCalculator.SetStartingIndex(m_currentLatticeStart);
  powerSpectrumCalculator.Update();
  auto temp = powerSpectrumCalculator.GetOutput();
  for (auto const &f : temp)
  {
    featurevec[f.first] = f.second;
  }
}


template< class TImage >
void FeatureExtraction< TImage >::CalculateGaborWavelets(const typename TImage::Pointer itkImage, std::map< std::string, double > &featurevec, bool latticePatch)
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
  gaborWaveletCalculator.Update();
  auto temp = gaborWaveletCalculator.GetOutput();
  for (auto const &f : temp)
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
  morphologicCalculator.Update();
  auto temp = morphologicCalculator.GetOutput();
  for (auto const &f : temp)
  {
    featurevec[f.first] = f.second;
  }
}

template< class TImage >
void FeatureExtraction< TImage >::CalculateNGTDM(const typename TImage::Pointer itkImage,
  const typename TImage::Pointer maskImage, OffsetVector *offset, std::map<std::string, double>& featurevec)
{
  NGTDMFeatures< TImage > calculator;
  calculator.SetInputImage(itkImage);
  calculator.SetInputMask(maskImage);
  calculator.SetNumBins(m_Bins);
  calculator.SetRange(m_Radius);
  calculator.SetMinimum(m_minimumToConsider);
  calculator.SetMaximum(m_maximumToConsider);
  calculator.SetStartingIndex(m_currentLatticeStart);
  calculator.Update();

  //auto temp = calculator.GetOutput();
  //double double_Strength = calculator.GetStrength();
  //double double_Complexity = calculator.GetComplexity();
  //double double_Coarsness = calculator.GetCoarsness();
  //double double_Contrast = calculator.GetContrast();
  //double double_Busyness = calculator.GetBusyness();
  //std::cout << "\n" << std::endl;
  //std::cout << "Strength =" << double_Strength << std::endl;
  //std::cout << "Complexity =" << double_Complexity << std::endl;
  //std::cout << "Coarsness =" << double_Coarsness << std::endl;
  //std::cout << "Contrast =" << double_Contrast << std::endl;
  //std::cout << "Busyness =" << double_Busyness << std::endl;

  featurevec["Strength"] = calculator.GetStrength();
  featurevec["Complexity"] = calculator.GetComplexity();
  featurevec["Coarseness"] = calculator.GetCoarsness();
  featurevec["Constrast"] = calculator.GetContrast();
  featurevec["Busyness"] = calculator.GetBusyness();
  /* commenting out old codes calling NGLDM
  typedef itk::Statistics::EnhancedScalarImageToNeighbourhoodGreyLevelDifferenceFeaturesFilter< TImage > FilterType;
  typedef typename FilterType::NeighbourhoodGreyLevelDifferenceFeaturesFilterType TextureFilterType;
  //typedef itk::MinimumMaximumImageCalculator< TImage > MinMaxComputerType;

  using MatrixGenerator = itk::Statistics::EnhancedScalarImageToNeighbourhoodGreyLevelDifferenceMatrixFilter< TImage >;
  using FeatureCalculator = itk::Statistics::EnhancedHistogramToNeighbourhoodGreyLevelDifferenceFeaturesFilter< typename MatrixGenerator::HistogramType >;

  typename MatrixGenerator::Pointer matrixFilter = MatrixGenerator::New();
  typename FeatureCalculator::Pointer featureFilter = FeatureCalculator::New();
  matrixFilter->SetPixelValueMinMax(m_minimumToConsider, m_maximumToConsider);
  matrixFilter->SetNumberOfBinsPerAxis(m_Bins);
  matrixFilter->SetInput(itkImage);
  matrixFilter->SetMaskImage(maskImage);

  //voxel count
  //unsigned long numberOfVoxels = 0;
  //itk::ImageRegionConstIterator<TImage> voxelCountIter(maskImage, maskImage->GetLargestPossibleRegion());
  //while (!voxelCountIter.IsAtEnd())
  //{
   // if (voxelCountIter.Get() > 0)
    //  ++numberOfVoxels;
   // ++voxelCountIter;
  //}

  //add offsets except self voxel
  typename OffsetVector::ConstIterator offsetIt;
  std::vector <OffsetType> tVector;
  for (offsetIt = offset->Begin(); offsetIt != offset->End(); offsetIt++)
  {

    tVector.push_back(offsetIt.Value());
  }
  matrixFilter->AddOffsets(tVector);
  matrixFilter->Update();

  //typename FilterType::FeatureNameVectorPointer requestedFeatures = FilterType::FeatureNameVector::New();

  //requestedFeatures->push_back(TextureFilterType::Coarseness);
  //requestedFeatures->push_back(TextureFilterType::Contrast);
  //requestedFeatures->push_back(TextureFilterType::Busyness);
  //requestedFeatures->push_back(TextureFilterType::Complexity);
  //requestedFeatures->push_back(TextureFilterType::Strength);
  //requestedFeatures->push_back(20);


  featureFilter->SetNumberOfVoxels(m_currentNonZeroImageValues.size());
  featureFilter->SetInput(matrixFilter->GetOutput());
  featureFilter->SetSiMatrix(matrixFilter->GetSiMatrix());
  featureFilter->Update();

  auto Coarseness = featureFilter->GetCoarseness();
  auto Contrast = featureFilter->GetContrast();
  auto Business = featureFilter->GetBusyness();
  auto Complexity = featureFilter->GetComplexity();
  auto Strength = featureFilter->GetStrength();
  std::cout << "\n Coarseness = " << Coarseness << std::endl;
  std::cout << "\n Contrast = " << Contrast << std::endl;
  std::cout << "\n Business = " << Business << std::endl;
  std::cout << "\n Complexity = " << Complexity << std::endl;
  std::cout << "\n Strength = " << Strength << std::endl;
  */

  //typename MatrixGenerator::OffsetVector::Pointer newOffset = MatrixGenerator::OffsetVector::New();
  //auto oldOffsets = matrixFilter->GetOffsets();
  //auto oldOffsetsIterator = oldOffsets->Begin();
  //while (oldOffsetsIterator != oldOffsets->End())
  //{
  //  bool continueOuterLoop = false;
  //  typename MatrixGenerator::OffsetType offset = oldOffsetsIterator->Value();
  //  for (size_t i = 0; i < TImage::ImageDimension; ++i)
  //  {
  //    //if (/*params.m_Direction == i + 2 &&*/)
  //    if (offset[i] != 0)
  //    {
  //      continueOuterLoop = true;
  //    }
  //  }
  //  oldOffsetsIterator++;
  //  if (continueOuterLoop)
  //    newOffset->push_back(offset);
  //}
  //matrixFilter->SetOffsets(newOffset);
  //matrixFilter->SetInput(itkImage);
  //matrixFilter->SetMaskImage(maskImage);
  //matrixFilter->SetPixelValueMinMax(m_minimumToConsider, m_maximumToConsider);
  //matrixFilter->SetNumberOfBinsPerAxis(m_Bins);

  //matrixFilter->Update();

  //unsigned long numberOfVoxels = 0;
  //TConstIteratorType voxelCountIter(maskImage, maskImage->GetLargestPossibleRegion());
  //while (!voxelCountIter.IsAtEnd())
  //{
  //  if (voxelCountIter.Get() > 0)
  //    ++numberOfVoxels;
  //  ++voxelCountIter;
  //}

  //featureFilter->SetInput(matrixFilter->GetOutput());
  //featureFilter->SetSiMatrix(matrixFilter->GetSiMatrix());
  //featureFilter->SetNumberOfVoxels(numberOfVoxels);
  //featureFilter->Update();

  //typedef short NeighbourhoodGreyLevelDifferenceFeatureName;
  //typedef itk::VectorContainer<unsigned char, NeighbourhoodGreyLevelDifferenceFeatureName> FeatureNameVector;
  //typedef typename FeatureNameVector::Pointer      FeatureNameVectorPointer;
  //FeatureNameVectorPointer requestedFeatures = FeatureNameVector::New();

  //requestedFeatures->push_back(FeatureCalculator::Coarseness);
  //requestedFeatures->push_back(FeatureCalculator::Contrast);
  //requestedFeatures->push_back(FeatureCalculator::Busyness);
  //requestedFeatures->push_back(FeatureCalculator::Complexity);
  //requestedFeatures->push_back(FeatureCalculator::Strength);

  //typedef typename FeatureCalculator::NeighbourhoodGreyLevelDifferenceFeatureName
  //  InternalNeighbourhoodGreyLevelDifferenceFeatureName;

  //int numFeatures = requestedFeatures->size();
  //double *features = new double[numFeatures];

  //itk::VectorContainer< unsigned char, double >::Pointer featureMeans, featureStd;

  //int /*offsetNum, */featureNum;
  //typename FeatureNameVector::ConstIterator fnameIt;
  //for (fnameIt = requestedFeatures->Begin(), featureNum = 0;
  //  fnameIt != requestedFeatures->End(); fnameIt++, featureNum++)
  //{
  //  InternalNeighbourhoodGreyLevelDifferenceFeatureName tn = (InternalNeighbourhoodGreyLevelDifferenceFeatureName)fnameIt.Value();
  //  double xx = featureFilter->GetFeature(tn);

  //  features[featureNum] = xx;
  //  featureMeans->push_back(xx);
  //}


  // can't directly set this->m_RequestedFeatures since it is const!

  //typename FilterType::Pointer filter = FilterType::New();

  //typename FilterType::OffsetVector::Pointer newOffset = FilterType::OffsetVector::New();
  //auto oldOffsets = filter->GetOffsets();
  //auto oldOffsetsIterator = oldOffsets->Begin();
  //while (oldOffsetsIterator != oldOffsets->End())
  //{
  //  bool continueOuterLoop = false;
  //  typename FilterType::OffsetType offset = oldOffsetsIterator->Value();
  //  for (size_t i = 0; i < TImage::ImageDimension; ++i)
  //  {
  //    if (/*params.m_Direction == i + 2 &&*/ offset[i] != 0)
  //    {
  //      continueOuterLoop = true;
  //    }
  //  }
  //  oldOffsetsIterator++;
  //  if (continueOuterLoop)
  //    newOffset->push_back(offset);
  //}
  //filter->SetOffsets(newOffset);


  //// All features are required
  //typename FilterType::FeatureNameVectorPointer requestedFeatures = FilterType::FeatureNameVector::New();
  //requestedFeatures->push_back(TextureFilterType::Coarseness);
  //requestedFeatures->push_back(TextureFilterType::Contrast);
  //requestedFeatures->push_back(TextureFilterType::Busyness);
  //requestedFeatures->push_back(TextureFilterType::Complexity);
  //requestedFeatures->push_back(TextureFilterType::Strength);

  //typename MinMaxComputerType::Pointer minMaxComputer = MinMaxComputerType::New();
  //minMaxComputer->SetImage(itkImage);
  //minMaxComputer->Compute();

  //auto duplicator_image = itk::ImageDuplicator< TImage >::New();
  //auto duplicator_mask = itk::ImageDuplicator< TImage >::New();
  //duplicator_image->SetInputImage(itkImage);
  //duplicator_image->Update();
  //duplicator_mask->SetInputImage(maskImage);
  //duplicator_mask->Update();

  //auto testimage = duplicator_image->GetOutput();
  //auto testmask = duplicator_mask->GetOutput();
  //filter->SetInput(testimage);
  //filter->SetMaskImage(testmask);
  //filter->SetRequestedFeatures(requestedFeatures);
  ////int rangeOfPixels = m_Range;
  //filter->SetPixelValueMinMax(m_minimumToConsider, m_maximumToConsider);
  //filter->SetNumberOfBinsPerAxis(m_Bins);

  ////filter->SetDistanceValueMinMax(0, m_Range);

  //filter->Update();

  //auto featureMeans = featureFilter->GetFeatureMeans();
  //auto featureStd = featureFilter->GetFeatureStandardDeviations();

  //std::ostringstream  ss;
  ////ss << rangeOfPixels;
  //std::string strRange = ss.str();
  //for (std::size_t i = 0; i < featureMeans->size(); ++i)
  //{
  //  switch (i)
  //  {
  //  case TextureFilterType::Coarseness:
  //    featurevec["Coarsness"] = featureMeans->ElementAt(i);
  //    break;
  //  case TextureFilterType::Contrast:
  //    featurevec["NeighbourContrast"] = featureMeans->ElementAt(i);
  //    break;
  //  case TextureFilterType::Busyness:
  //    featurevec["Busyness"] = featureMeans->ElementAt(i);
  //    break;
  //  case TextureFilterType::Complexity:
  //    featurevec["Complexity"] = featureMeans->ElementAt(i);
  //    break;
  //  case TextureFilterType::Strength:
  //    featurevec["Strength"] = featureMeans->ElementAt(i);
  //    break;
  //  default:
  //    break;
  //  }
  //}

}


template< class TImage >
void FeatureExtraction< TImage >::CalculateGLSZM(const typename TImage::Pointer itkImage, const typename TImage::Pointer maskImage, OffsetVector *offset, std::map<std::string, double>& featurevec)
{
  GLSZMFeatures< TImage > glszmCalculator;
  glszmCalculator.SetInputImage(itkImage);
  glszmCalculator.SetInputMask(maskImage);
  glszmCalculator.SetNumBins(m_Bins);
  glszmCalculator.SetRadius(m_Radius);
  glszmCalculator.SetRadius(m_Radius_float);
  glszmCalculator.SetMaxSize(m_Range);
  glszmCalculator.SetMinimum(m_minimumToConsider);
  glszmCalculator.SetMaximum(m_maximumToConsider);
  glszmCalculator.SetStartingIndex(m_currentLatticeStart);
  glszmCalculator.SetOffsets(offset);
  glszmCalculator.Update();
  auto temp = glszmCalculator.GetOutput();
  for (auto const &f : temp)
  {
    featurevec[f.first] = f.second;
    //TBD - for debugging GLSZM features
    //std::cout << "[DEBUG] GLSZM Feature '" << f.first << "' = " << f.second << "\n";
    //TBD - for debugging GLSZM features
  }
  //using FilterType = itk::Statistics::EnhancedScalarImageToSizeZoneFeaturesFilter< TImage >;
  //typedef typename FilterType::SizeZoneFeaturesFilterType TextureFilterType;

  //auto filter = FilterType::New();
  //filter->SetInput(itkImage);
  //filter->SetMaskImage(maskImage);
  //filter->SetInsidePixelValue(1); //maskImage thrown into CalculateGLSZM filter should have already been converted to binary (0 for outside and 1 for inside) mask image, so set inside pixel value of maskImage to 1
  //filter->SetNumberOfBinsPerAxis(m_Bins); //for quantization of grey levels
  //filter->SetPixelValueMinMax(m_minimumToConsider, m_maximumToConsider); //Set the min and max(inclusive) pixel value that will be used for feature calculations.Optional; for default value see above.
  //// All features are required
  //typename FilterType::FeatureNameVectorPointer requestedFeatures = FilterType::FeatureNameVector::New();
  //requestedFeatures->push_back(TextureFilterType::SmallZoneEmphasis);
  //requestedFeatures->push_back(TextureFilterType::LargeZoneEmphasis);
  //requestedFeatures->push_back(TextureFilterType::GreyLevelNonuniformity);
  //requestedFeatures->push_back(TextureFilterType::GreyLevelNonuniformityNormalized);
  //requestedFeatures->push_back(TextureFilterType::SizeZoneNonuniformity);
  //requestedFeatures->push_back(TextureFilterType::SizeZoneNonuniformityNormalized);
  //requestedFeatures->push_back(TextureFilterType::LowGreyLevelZoneEmphasis);
  //requestedFeatures->push_back(TextureFilterType::HighGreyLevelZoneEmphasis);
  //requestedFeatures->push_back(TextureFilterType::SmallZoneLowGreyLevelEmphasis);
  //requestedFeatures->push_back(TextureFilterType::SmallZoneHighGreyLevelEmphasis);
  //requestedFeatures->push_back(TextureFilterType::LargeZoneLowGreyLevelEmphasis);
  //requestedFeatures->push_back(TextureFilterType::LargeZoneHighGreyLevelEmphasis);
  //requestedFeatures->push_back(TextureFilterType::ZonePercentage);
  //requestedFeatures->push_back(TextureFilterType::GreyLevelVariance);
  //requestedFeatures->push_back(TextureFilterType::SizeZoneVariance);
  //requestedFeatures->push_back(TextureFilterType::ZoneEntropy);

  //filter->SetRequestedFeatures(requestedFeatures);
  //filter->Update();
  ///*int rangeOfPixels = m_Range;
  //if (rangeOfPixels < 2)
  //rangeOfPixels = 256;

  //if (params.m_UseCtRange)
  //{
  //filter->SetPixelValueMinMax((TPixel)(-1024.5),(TPixel)(3096.5));
  //filter->SetNumberOfBinsPerAxis(3096.5+1024.5);
  //} else*/
  ////{
  ////  filter->SetPixelValueMinMax(m_minimumToConsider, m_maximumToConsider);
  ////  filter->SetNumberOfBinsPerAxis(m_Bins);
  ////}

  ////filter->SetDistanceValueMinMax(0, rangeOfPixels);

  ////filter->Update();

  //auto featureMeans = filter->GetFeatureMeans();
  //auto featureStd = filter->GetFeatureStandardDeviations();

  //std::cout << "\n" << std::endl;
  //for (std::size_t i = 0; i < featureMeans->size(); ++i)
  //{
  //  switch (i)
  //  {
  //  case TextureFilterType::SmallZoneEmphasis:
  //    featurevec["ZE_Mean"] = featureMeans->ElementAt(i);
  //    featurevec["ZE_STD"] = featureStd->ElementAt(i);
  //    //std::cout << "[" << i << "]|SmallZoneEmphasis:" << featureMeans->ElementAt(i) << std::endl;
  //    //std::cout << "[" << i << "]|SmallZoneEmphasis_STD:" << featureStd->ElementAt(i) << std::endl;
  //    break;
  //  case TextureFilterType::LargeZoneEmphasis:
  //    featurevec["LZE_Mean"] = featureMeans->ElementAt(i);
  //    featurevec["LZE_STD"] = featureStd->ElementAt(i);
  //    //std::cout << "[" << i << "]|LargeZoneEmphasis:" << featureMeans->ElementAt(i) << std::endl;
  //    //std::cout << "[" << i << "]|LargeZoneEmphasis_STD:" << featureStd->ElementAt(i) << std::endl;
  //    break;
  //  case TextureFilterType::GreyLevelNonuniformity:
  //    featurevec["GLN_Mean"] = featureMeans->ElementAt(i);
  //    featurevec["GLN_STD"] = featureStd->ElementAt(i);
  //    //std::cout << "[" << i << "]|GreyLevelNonuniformity:" << featureMeans->ElementAt(i) << std::endl;
  //    //std::cout << "[" << i << "]|GreyLevelNonuniformity_STD:" << featureStd->ElementAt(i) << std::endl;
  //    break;
  //  case TextureFilterType::GreyLevelNonuniformityNormalized:
  //    featurevec["GLNNorm_Mean"] = featureMeans->ElementAt(i);
  //    featurevec["GLNNorm_STD"] = featureStd->ElementAt(i);
  //    //std::cout << "[" << i << "]|GreyLevelNonuniformityNormalized:" << featureMeans->ElementAt(i) << std::endl;
  //    //std::cout << "[" << i << "]|GreyLevelNonuniformityNormalized_STD:" << featureStd->ElementAt(i) << std::endl;
  //    break;
  //  case TextureFilterType::SizeZoneNonuniformity:
  //    featurevec["ZSN_Mean"] = featureMeans->ElementAt(i);
  //    featurevec["ZSN_STD"] = featureStd->ElementAt(i);
  //    //std::cout << "[" << i << "]|SizeZoneNonuniformity:" << featureMeans->ElementAt(i) << std::endl;
  //    //std::cout << "[" << i << "]|SizeZoneNonuniformity_STD:" << featureStd->ElementAt(i) << std::endl;
  //    break;
  //  case TextureFilterType::SizeZoneNonuniformityNormalized:
  //    featurevec["ZSNNorm_Mean"] = featureMeans->ElementAt(i);
  //    featurevec["ZSNNorm_STD"] = featureStd->ElementAt(i);
  //    //std::cout << "[" << i << "]|SizeZoneNonuniformityNormalized:" << featureMeans->ElementAt(i) << std::endl;
  //    //std::cout << "[" << i << "]|SizeZoneNonuniformityNormalized_STD:" << featureStd->ElementAt(i) << std::endl;
  //    break;
  //  case TextureFilterType::LowGreyLevelZoneEmphasis:
  //    featurevec["LGZE_Mean"] = featureMeans->ElementAt(i);
  //    featurevec["LGZE_STD"] = featureStd->ElementAt(i);
  //    //std::cout << "[" << i << "]|LowGreyLevelZoneEmphasis:" << featureMeans->ElementAt(i) << std::endl;
  //    //std::cout << "[" << i << "]|LowGreyLevelZoneEmphasis_STD:" << featureStd->ElementAt(i) << std::endl;
  //    break;
  //  case TextureFilterType::HighGreyLevelZoneEmphasis:
  //    featurevec["HGZE_Mean"] = featureMeans->ElementAt(i);
  //    featurevec["HGZE_STD"] = featureStd->ElementAt(i);
  //    //std::cout << "[" << i << "]|HighGreyLevelZoneEmphasis:" << featureMeans->ElementAt(i) << std::endl;
  //    //std::cout << "[" << i << "]|HighGreyLevelZoneEmphasis_STD:" << featureStd->ElementAt(i) << std::endl;
  //    break;
  //  case TextureFilterType::SmallZoneLowGreyLevelEmphasis:
  //    featurevec["SZLGE_Mean"] = featureMeans->ElementAt(i);
  //    featurevec["SZLGE_STD"] = featureStd->ElementAt(i);
  //    //std::cout << "[" << i << "]|SmallZoneLowGreyLevelEmphasis:" << featureMeans->ElementAt(i) << std::endl;
  //    //std::cout << "[" << i << "]|SmallZoneLowGreyLevelEmphasis_STD:" << featureStd->ElementAt(i) << std::endl;
  //    break;
  //  case TextureFilterType::SmallZoneHighGreyLevelEmphasis:
  //    featurevec["SZHGE_Mean"] = featureMeans->ElementAt(i);
  //    featurevec["SZHGE_STD"] = featureStd->ElementAt(i);
  //    //std::cout << "[" << i << "]|SmallZoneHighGreyLevelEmphasis:" << featureMeans->ElementAt(i) << std::endl;
  //    //std::cout << "[" << i << "]|SmallZoneHighGreyLevelEmphasis_STD:" << featureStd->ElementAt(i) << std::endl;
  //    break;
  //  case TextureFilterType::LargeZoneLowGreyLevelEmphasis:
  //    featurevec["LZLGE_Mean"] = featureMeans->ElementAt(i);
  //    featurevec["LZLGE_STD"] = featureStd->ElementAt(i);
  //    //std::cout << "[" << i << "]|LargeZoneLowGreyLevelEmphasis:" << featureMeans->ElementAt(i) << std::endl;
  //    //std::cout << "[" << i << "]|LargeZoneLowGreyLevelEmphasis_STD:" << featureStd->ElementAt(i) << std::endl;
  //    break;
  //  case TextureFilterType::LargeZoneHighGreyLevelEmphasis:
  //    featurevec["LZHGE_Mean"] = featureMeans->ElementAt(i);
  //    featurevec["LZHGE_STD"] = featureStd->ElementAt(i);
  //    //std::cout << "[" << i << "]|LargeZoneHighGreyLevelEmphasis:" << featureMeans->ElementAt(i) << std::endl;
  //    //std::cout << "[" << i << "]|LargeZoneHighGreyLevelEmphasis_STD:" << featureStd->ElementAt(i) << std::endl;
  //    break;
  //  case TextureFilterType::ZonePercentage:
  //    featurevec["ZP_Mean"] = featureMeans->ElementAt(i);
  //    featurevec["ZP_STD"] = featureStd->ElementAt(i);
  //    //std::cout << "[" << i << "]|ZonePercentage:" << featureMeans->ElementAt(i) << std::endl;
  //    //std::cout << "[" << i << "]|ZonePercentage_STD:" << featureStd->ElementAt(i) << std::endl;
  //    break;
  //  case TextureFilterType::GreyLevelVariance:
  //    featurevec["GLV_Mean"] = featureMeans->ElementAt(i);
  //    featurevec["GLV_STD"] = featureStd->ElementAt(i);
  //    //std::cout << "[" << i << "]|GreyLevelVariance:" << featureMeans->ElementAt(i) << std::endl;
  //    //std::cout << "[" << i << "]|GreyLevelVariance_STD:" << featureStd->ElementAt(i) << std::endl;
  //    break;
  //  case TextureFilterType::SizeZoneVariance:
  //    featurevec["SZV_Mean"] = featureMeans->ElementAt(i);
  //    featurevec["SZV_STD"] = featureStd->ElementAt(i);
  //    //std::cout << "[" << i << "]|SizeZoneVariance:" << featureMeans->ElementAt(i) << std::endl;
  //    //std::cout << "[" << i << "]|SizeZoneVariance_STD:" << featureStd->ElementAt(i) << std::endl;
  //    break;
  //  case TextureFilterType::ZoneEntropy:
  //    featurevec["ZE_Mean"] = featureMeans->ElementAt(i);
  //    featurevec["ZE_STD"] = featureStd->ElementAt(i);
  //    //std::cout << "[" << i << "]|ZoneEntropy:" << featureMeans->ElementAt(i) << std::endl;
  //    //std::cout << "[" << i << "]|ZoneEntropy_STD:" << featureStd->ElementAt(i) << std::endl;
  //    break;
  //  default:
  //    break;
  //  }
  //}
}


template< class TImage >
void FeatureExtraction< TImage >::CalculateHistogram(const typename TImage::Pointer image, const typename TImage::Pointer mask, std::map< std::string, double > &featurevec, bool latticePatch)
{
  /// histogram calculation from ITK -- for texture feature pipeline
  typename TImage::PixelType min, max;

  if (m_QuantizationType == "Image")
  {
    min = m_minimumToConsider;
    max = m_maximumToConsider;
  }
  if (m_QuantizationType == "ROI")
  {
    min = m_statistics_local.GetMinimum();
    max = m_statistics_local.GetMaximum();
  }


  using HistogramFilterType = itk::Statistics::MaskedImageToHistogramFilter< TImage, TImage >;
  using HistogramMeasurementType = typename HistogramFilterType::HistogramType::MeasurementVectorType;

  HistogramMeasurementType lowerBound(m_Bins), upperBound(m_Bins);
  lowerBound.Fill(min);
  upperBound.Fill(max);

  typename HistogramFilterType::HistogramType::SizeType size(1); // this is always grayscale
  size.Fill(m_Bins);

  auto histogramCalculator = HistogramFilterType::New();
  histogramCalculator->SetInput(image);
  histogramCalculator->SetMaskImage(mask);
  histogramCalculator->SetMaskValue(1);
  histogramCalculator->SetHistogramBinMinimum(lowerBound);
  histogramCalculator->SetHistogramBinMaximum(upperBound);
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
  /// histogram calculation from ITK -- for texture feature pipeline ends

  //for (size_t i = 0; i < histogram->GetSize()[0]; i++)
  //{
  //  featurevec["Bin_" + std::to_string(i) + "_Frequency"] = histogram->GetFrequency(i);
  //  featurevec["Bin_" + std::to_string(i) + "_Max"] = histogram->GetBinMax(0, i);
  //  featurevec["Bin_" + std::to_string(i) + "_Min"] = histogram->GetBinMin(0, i);
  //}

  //const float maxRescaleVal = 1000;
  //typedef itk::RescaleIntensityImageFilter< TImage, TImage > RescaleFilterType;
  //typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  //rescaleFilter->SetInput(image);
  //rescaleFilter->SetOutputMinimum(0);
  //rescaleFilter->SetOutputMaximum(maxRescaleVal);
  //rescaleFilter->Update();

  //std::vector<double> intensities;
  //itk::ImageRegionConstIterator <TImage> imageIt(rescaleFilter->GetOutput(), rescaleFilter->GetOutput()->GetLargestPossibleRegion());
  //itk::ImageRegionConstIterator <TImage> maskIt(mask, mask->GetLargestPossibleRegion());
  //imageIt.GoToBegin();
  //maskIt.GoToBegin();

  //while (!imageIt.IsAtEnd())
  //{
  //  if (maskIt.Get() > 0)
  //    intensities.push_back(std::round(imageIt.Get()));
  //  ++imageIt;
  //  ++maskIt;
  //}

  ////-------------------------------
  //double interval = maxRescaleVal / m_Bins;
  //double final_interval = (int)(interval * 100);
  //final_interval = (double)final_interval / 100;


  //std::vector<double> finalBins;
  //std::vector<std::vector<double>> Ranges;
  //double current_index = 0;
  //for (int i = 0; i < m_Bins; i++)
  //{
  //  std::vector<double> onerange;
  //  onerange.push_back(current_index);
  //  current_index = current_index + final_interval;
  //  onerange.push_back(current_index);
  //  Ranges.push_back(onerange);

  //  if (static_cast<int>(Ranges.size()) == m_Bins)
  //    Ranges[Ranges.size() - 1][1] = maxRescaleVal;
  //}
  ////toadd the last one
  //for (unsigned int j = 0; j < Ranges.size(); j++)
  //{
  //  std::vector<double> onerange = Ranges[j];
  //  int counter = 0;
  //  for (unsigned int i = 0; i < intensities.size(); i++)
  //  {
  //    if (onerange[0] == 0)
  //    {
  //      if (intensities[i] >= onerange[0] && intensities[i] <= onerange[1])
  //        counter = counter + 1;
  //    }
  //    else
  //    {
  //      if (intensities[i] > onerange[0] && intensities[i] <= onerange[1])
  //        counter = counter + 1;
  //    }
  //  }
  //  finalBins.push_back(counter);
  //}
  //for (unsigned int j = 0; j < finalBins.size(); j++)
  //{
  //  finalBins[j] = (finalBins[j] * 100) / intensities.size();
  //  featurevec["Bin_" + std::to_string(j)] = finalBins[j];
  //  featurevec["BinEndIntensity_" + std::to_string(j)] = Ranges[j][1];
  //}

  //HistogramFeatures< TImage > histogramCalculator;
  //histogramCalculator.SetInputImage(image);
  //histogramCalculator.SetInputMask(mask);
  //histogramCalculator.SetStartingIndex(m_currentLatticeStart);
  //histogramCalculator.SetMinimum(min);
  //histogramCalculator.SetMaximum(min);
  //histogramCalculator.Update();
  //auto temp = histogramCalculator.GetOutput();
  //for (auto const &f : temp)
  //{
  //  featurevec[f.first] = f.second;
  //}
}


template< class TImage >
void FeatureExtraction< TImage >::CalculateGLRLM(const typename TImage::Pointer image, const typename TImage::Pointer mask, OffsetVector *offset, std::map<std::string, double> &featurevec, bool latticePatch)
{
  using HistogramFrequencyContainerType = itk::Statistics::DenseFrequencyContainer2;
  using RunLengthFilterType = itk::Statistics::EnhancedScalarImageToRunLengthFeaturesFilter< TImage, HistogramFrequencyContainerType >;
  //using RunLengthFilterType = itk::Statistics::ScalarImageToRunLengthFeaturesFilter< TImage, HistogramFrequencyContainerType >;
  using RunLengthMatrixGenerator = typename RunLengthFilterType::RunLengthMatrixFilterType;
  using RunLengthFeatures = typename RunLengthFilterType::RunLengthFeaturesFilterType;

  typename  RunLengthMatrixGenerator::Pointer matrix_generator = RunLengthMatrixGenerator::New();
  matrix_generator->SetInput(image);
  matrix_generator->SetMaskImage(mask);
  matrix_generator->SetInsidePixelValue(1);
  matrix_generator->SetPixelValueMinMax(m_minimumToConsider, m_maximumToConsider);

  if (latticePatch)
  {
    auto maxStep = m_latticeSizeImage[0];
    for (size_t d = 1; d < TImage::ImageDimension; d++)
    {
      if (maxStep < m_latticeSizeImage[d])
      {
        maxStep = m_latticeSizeImage[d];
      }
    }
    matrix_generator->SetDistanceValueMinMax(0, std::sqrt(2) * (maxStep - 1));
  }
  //removed SetDistanceValueMinMax
  //matrix_generator->SetDistanceValueMinMax(0, m_Range); // TOCHECK - why is this only between 0-4? P
  matrix_generator->SetNumberOfBinsPerAxis(m_Bins); // TOCHECK - needs to be statistically significant

  typename  RunLengthFeatures::Pointer runLengthMatrixCalculator = RunLengthFeatures::New();

  typename  OffsetVector::ConstIterator offsetIt;
  size_t offsetNum = 0;

  auto size = image->GetBufferedRegion().GetSize();
  double size_total = size[0];
  for (size_t d = 1; d < TImage::ImageDimension; d++)
  {
    size_total *= size[d];
  }

  if ((m_offsetSelect == "Average") || (m_offsetSelect == "Individual"))
  {
    double sre = 0, lre = 0, gln = 0, glnn = 0, rln = 0, rlnn = 0, rp = 0, lglre = 0, hglre = 0, srlgle = 0, srhgle = 0, lrlgle = 0, lrhgle = 0,
      runs = 0, glv = 0, rlv = 0, re = 0;

    for (offsetIt = offset->Begin(); offsetIt != offset->End(); offsetIt++, offsetNum++)
    {
      matrix_generator->SetOffset(offsetIt.Value());
      matrix_generator->Update();

      runLengthMatrixCalculator->SetInput(matrix_generator->GetOutput());
      runLengthMatrixCalculator->Update();

      if (m_offsetSelect == "Average")
      {
        sre += runLengthMatrixCalculator->GetShortRunEmphasis();
        lre += runLengthMatrixCalculator->GetLongRunEmphasis();
        gln += runLengthMatrixCalculator->GetGreyLevelNonuniformity();
        rln += runLengthMatrixCalculator->GetRunLengthNonuniformity();
        lglre += runLengthMatrixCalculator->GetLowGreyLevelRunEmphasis();
        hglre += runLengthMatrixCalculator->GetHighGreyLevelRunEmphasis();
        srlgle += runLengthMatrixCalculator->GetShortRunLowGreyLevelEmphasis();
        srhgle += runLengthMatrixCalculator->GetShortRunHighGreyLevelEmphasis();
        lrlgle += runLengthMatrixCalculator->GetLongRunLowGreyLevelEmphasis();
        lrhgle += runLengthMatrixCalculator->GetLongRunHighGreyLevelEmphasis();
        runs += runLengthMatrixCalculator->GetTotalNumberOfRuns();
        rp += static_cast<double>(runLengthMatrixCalculator->GetTotalNumberOfRuns()) / static_cast<double>(m_currentNonZeroImageValues.size());
        rlnn += runLengthMatrixCalculator->GetRunLengthNonuniformityNormalized();
        glnn += runLengthMatrixCalculator->GetGreyLevelNonuniformityNormalized();
        glv += runLengthMatrixCalculator->GetGreyLevelVariance();
        rlv += runLengthMatrixCalculator->GetRunLengthVariance();
        re += runLengthMatrixCalculator->GetRunEntropy();
      }
      else // individual
      {
        featurevec["ShortRunEmphasis_Offset_" + std::to_string(offsetNum)] = runLengthMatrixCalculator->GetShortRunEmphasis();
        featurevec["LongRunEmphasis_Offset_" + std::to_string(offsetNum)] = runLengthMatrixCalculator->GetLongRunEmphasis();
        featurevec["GreyLevelNonuniformity_Offset_" + std::to_string(offsetNum)] = runLengthMatrixCalculator->GetGreyLevelNonuniformity();
        featurevec["RunLengthNonuniformity_Offset_" + std::to_string(offsetNum)] = runLengthMatrixCalculator->GetRunLengthNonuniformity();
        featurevec["LowGreyLevelRunEmphasis_Offset_" + std::to_string(offsetNum)] = runLengthMatrixCalculator->GetLowGreyLevelRunEmphasis();
        featurevec["HighGreyLevelRunEmphasis_Offset_" + std::to_string(offsetNum)] = runLengthMatrixCalculator->GetHighGreyLevelRunEmphasis();
        featurevec["ShortRunLowGreyLevelEmphasis_Offset_" + std::to_string(offsetNum)] = runLengthMatrixCalculator->GetShortRunLowGreyLevelEmphasis();
        featurevec["ShortRunHighGreyLevelEmphasis_Offset_" + std::to_string(offsetNum)] = runLengthMatrixCalculator->GetShortRunHighGreyLevelEmphasis();
        featurevec["LongRunLowGreyLevelEmphasis_Offset_" + std::to_string(offsetNum)] = runLengthMatrixCalculator->GetLongRunLowGreyLevelEmphasis();
        featurevec["LongRunHighGreyLevelEmphasis_Offset_" + std::to_string(offsetNum)] = runLengthMatrixCalculator->GetLongRunHighGreyLevelEmphasis();
        featurevec["TotalRuns_Offset_" + std::to_string(offsetNum)] = runLengthMatrixCalculator->GetTotalNumberOfRuns();
        featurevec["RunPercentage_Offset_" + std::to_string(offsetNum)] = featurevec["TotalRuns_Offset_" + std::to_string(offsetNum)] / static_cast<double>(m_currentNonZeroImageValues.size());
        featurevec["RunLengthNonuniformityNormalized_Offset_" + std::to_string(offsetNum)] = runLengthMatrixCalculator->GetRunLengthNonuniformityNormalized();
        featurevec["GreyLevelNonuniformityNormalized_Offset_" + std::to_string(offsetNum)] = runLengthMatrixCalculator->GetGreyLevelNonuniformityNormalized();
        featurevec["GreyLevelVariance_Offset_" + std::to_string(offsetNum)] = runLengthMatrixCalculator->GetGreyLevelVariance();
        featurevec["RunLengthVariance_Offset_" + std::to_string(offsetNum)] = runLengthMatrixCalculator->GetRunLengthVariance();
        featurevec["RunEntropy_Offset_" + std::to_string(offsetNum)] = runLengthMatrixCalculator->GetRunEntropy();
      }
    }
    if (m_offsetSelect == "Average")
    {
      sre /= offset->size();
      lre /= offset->size();
      gln /= offset->size();
      rln /= offset->size();
      lglre /= offset->size();
      hglre /= offset->size();
      srlgle /= offset->size();
      srhgle /= offset->size();
      lrlgle /= offset->size();
      lrhgle /= offset->size();
      rp /= offset->size();
      runs /= offset->size();
      rlnn /= offset->size();
      glnn /= offset->size();
      glv /= offset->size();
      rlv /= offset->size();
      re /= offset->size();

      featurevec["ShortRunEmphasis"] = sre;
      featurevec["LongRunEmphasis"] = lre;
      featurevec["GreyLevelNonuniformity"] = gln;
      featurevec["RunLengthNonuniformity"] = rln;
      featurevec["RunPercentage"] = rp;
      featurevec["LowGreyLevelRunEmphasis"] = lglre;
      featurevec["HighGreyLevelRunEmphasis"] = hglre;
      featurevec["ShortRunLowGreyLevelEmphasis"] = srlgle;
      featurevec["ShortRunHighGreyLevelEmphasis"] = srhgle;
      featurevec["LongRunLowGreyLevelEmphasis"] = lrlgle;
      featurevec["LongRunHighGreyLevelEmphasis"] = lrhgle;
      featurevec["TotalRuns"] = runs;
      featurevec["RunLengthNonuniformityNormalized"] = rlnn;
      featurevec["GreyLevelNonuniformityNormalized"] = glnn;
      featurevec["GreyLevelVariance_Offset"] = glv;
      featurevec["RunLengthVariance_Offset"] = rlv;
      featurevec["RunEntropy"] = re;
    }
  }
  else if ((m_offsetSelect == "ITKDefault") || (m_offsetSelect == "Combined"))
  {
    matrix_generator->SetOffsets(offset);
    matrix_generator->Update();

    auto temp = matrix_generator->GetOutput();

    runLengthMatrixCalculator->SetInput(matrix_generator->GetOutput());
    runLengthMatrixCalculator->Update();

    featurevec["ShortRunEmphasis"] = runLengthMatrixCalculator->GetShortRunEmphasis();
    featurevec["LongRunEmphasis"] = runLengthMatrixCalculator->GetLongRunEmphasis();
    featurevec["GreyLevelNonuniformity"] = runLengthMatrixCalculator->GetGreyLevelNonuniformity();
    featurevec["RunLengthNonuniformity"] = runLengthMatrixCalculator->GetRunLengthNonuniformity();
    featurevec["LowGreyLevelRunEmphasis"] = runLengthMatrixCalculator->GetLowGreyLevelRunEmphasis();
    featurevec["HighGreyLevelRunEmphasis"] = runLengthMatrixCalculator->GetHighGreyLevelRunEmphasis();
    featurevec["ShortRunLowGreyLevelEmphasis"] = runLengthMatrixCalculator->GetShortRunLowGreyLevelEmphasis();
    featurevec["ShortRunHighGreyLevelEmphasis"] = runLengthMatrixCalculator->GetShortRunHighGreyLevelEmphasis();
    featurevec["LongRunLowGreyLevelEmphasis"] = runLengthMatrixCalculator->GetLongRunLowGreyLevelEmphasis();
    featurevec["LongRunHighGreyLevelEmphasis"] = runLengthMatrixCalculator->GetLongRunHighGreyLevelEmphasis();
    featurevec["TotalRuns"] = runLengthMatrixCalculator->GetTotalNumberOfRuns();
    featurevec["RunPercentage"] = featurevec["TotalRuns"] / static_cast<double>(offset->size() * m_currentNonZeroImageValues.size());
    featurevec["RunLengthNonuniformityNormalized"] = runLengthMatrixCalculator->GetRunLengthNonuniformityNormalized();
    featurevec["GreyLevelNonuniformityNormalized"] = runLengthMatrixCalculator->GetGreyLevelNonuniformityNormalized();
    featurevec["GreyLevelVariance"] = runLengthMatrixCalculator->GetGreyLevelVariance();
    featurevec["RunLengthVariance"] = runLengthMatrixCalculator->GetRunLengthVariance();
    featurevec["RunEntropy"] = runLengthMatrixCalculator->GetRunEntropy();
  }
  else
  {
    // not defined, so don't do anything to featurevec
  }

}


template< class TImage >
void FeatureExtraction< TImage >::CalculateGLCM(const typename TImage::Pointer image, const typename TImage::Pointer mask, OffsetVector *offset, std::map<std::string, double> &featurevec, bool latticePatch)
{
  using Image2CoOccuranceType = itk::Statistics::ScalarImageToCooccurrenceMatrixFilter < TImage >;
  using HistogramType = typename Image2CoOccuranceType::HistogramType;
  using Hist2FeaturesType = itk::Statistics::HistogramToTextureFeaturesFilter< HistogramType >;

  double contrast = 0, correl = 0, ener = 0, entro = 0, homo = 0, clustershade = 0, clusterprominance = 0, autocorr = 0;

  auto image_wrap = image;
  auto mask_wrap = mask;

  //if (latticePatch)
  //{
  //  image_wrap = GetPatchedImage(image);
  //  mask_wrap = cbica::CreateImage< TImage >(image_wrap, 1);
  //}

  if (m_offsetSelect == "Average")
  {
    for (size_t i = 0; i < offset->size(); i++)
    {
      typename  Image2CoOccuranceType::Pointer glcmGenerator = Image2CoOccuranceType::New();
      glcmGenerator->SetNumberOfBinsPerAxis(m_Bins); //reasonable number of bins
      glcmGenerator->SetPixelValueMinMax(m_minimumToConsider, m_maximumToConsider);
      glcmGenerator->SetMaskImage(mask_wrap);
      glcmGenerator->SetInput(image_wrap);
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

    featurevec["Energy"] = ener;
    featurevec["Entropy"] = entro;
    featurevec["Correlation"] = correl;
    featurevec["Homogeneity"] = homo; // also called "inverse difference moment"
    featurevec["Contrast"] = contrast; // also called "inertia"
    featurevec["ClusterShade"] = clustershade;
    featurevec["ClusterProminence"] = clusterprominance;
    featurevec["AutoCorrelation"] = autocorr; // called "haralick"
  }
  else if ((m_offsetSelect == "ITKDefault") || (m_offsetSelect == "Combined"))
  {
    typename  Image2CoOccuranceType::Pointer glcmGenerator = Image2CoOccuranceType::New();
    glcmGenerator->SetNumberOfBinsPerAxis(m_Bins); //reasonable number of bins
    glcmGenerator->SetPixelValueMinMax(m_minimumToConsider, m_maximumToConsider);
    glcmGenerator->SetMaskImage(mask_wrap);
    glcmGenerator->SetInput(image_wrap);
    auto featureCalc = Hist2FeaturesType::New();

    glcmGenerator->SetOffsets(offset);
    glcmGenerator->Update();
    featureCalc->SetInput(glcmGenerator->GetOutput());
    featureCalc->Update();

    featurevec["Energy"] = static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::Energy));
    featurevec["Entropy"] = static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::Entropy));
    featurevec["Correlation"] = static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::Correlation));
    featurevec["Homogeneity"] = static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::InverseDifferenceMoment)); // also called "difference moment"
    featurevec["Contrast"] = static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::Inertia)); // also called "inertia"
    featurevec["ClusterShade"] = static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::ClusterShade));
    featurevec["ClusterProminence"] = static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::ClusterProminence));
    featurevec["AutoCorrelation"] = static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::HaralickCorrelation)); // called "haralick"
  }
  else
  {
    typename  Image2CoOccuranceType::Pointer glcmGenerator = Image2CoOccuranceType::New();
    glcmGenerator->SetNumberOfBinsPerAxis(m_Bins); //reasonable number of bins
    glcmGenerator->SetPixelValueMinMax(m_minimumToConsider, m_maximumToConsider);
    glcmGenerator->SetMaskImage(mask_wrap);
    glcmGenerator->SetInput(image_wrap);
    auto featureCalc = Hist2FeaturesType::New();

    for (size_t i = 0; i < offset->size(); i++)
    {
      glcmGenerator->SetOffset(offset->at(i));
      glcmGenerator->Update();
      featureCalc->SetInput(glcmGenerator->GetOutput());
      featureCalc->Update();

      auto tempStr = "_Offset_" + std::to_string(i);
      featurevec[std::string("Energy") + tempStr] = featureCalc->GetFeature(Hist2FeaturesType::Energy);
      featurevec[std::string("Entropy") + tempStr] = featureCalc->GetFeature(Hist2FeaturesType::Entropy);
      featurevec[std::string("Correlation") + tempStr] = featureCalc->GetFeature(Hist2FeaturesType::Correlation);
      featurevec[std::string("Homogeneity") + tempStr] = featureCalc->GetFeature(Hist2FeaturesType::InverseDifferenceMoment); // also called "difference moment"
      featurevec[std::string("Contrast") + tempStr] = featureCalc->GetFeature(Hist2FeaturesType::Inertia); // also called "inertia"
      featurevec[std::string("ClusterShade") + tempStr] = featureCalc->GetFeature(Hist2FeaturesType::ClusterShade);
      featurevec[std::string("ClusterProminence") + tempStr] = featureCalc->GetFeature(Hist2FeaturesType::ClusterProminence);
      featurevec[std::string("AutoCorrelation") + tempStr] = featureCalc->GetFeature(Hist2FeaturesType::HaralickCorrelation); // called "haralick"
    }
  }
  // TODO: Sung to add his GLCM extraction code here
  //featurevec[std::string("Correlation]) + "_Sung"] = 0;
}


template< class TImage >
void FeatureExtraction< TImage >::CalculateIntensity(std::vector< typename TImage::PixelType >& nonZeroVoxels, std::map< std::string, double > &featurevec, bool latticePatch)
{
  cbica::Statistics< typename TImage::PixelType > statisticsCalculatorToUse;
  if (m_QuantizationType == "Image")
  {
    statisticsCalculatorToUse = m_statistics_global[m_currentROIValue];
  }
  if (m_QuantizationType == "ROI")
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
  //featurevec["RobustMeanAbsoluteDeviation1090"] = statisticsCalculatorToUse.GetRobustMeanAbsoluteDeviation(10, 90);
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

  //auto temp1 = extractor->GetOutput();

  //auto roiExtractor = itk::RegionOfInterestImageFilter< TImage, TImage >::New();
  //roiExtractor->SetRegionOfInterest(region);
  //roiExtractor->SetInput(inputImage);
  //roiExtractor->Update();

  //auto temp2 = extractor->GetOutput();

  return extractor->GetOutput();
}

template< class TImage >
void FeatureExtraction< TImage >::SetFeatureParam(std::string featureFamily)
{
  auto temp = m_Features.find(featureFamily);
  if (temp != m_Features.end())
  {
    auto featurefamily = temp->second;
    auto parameters = std::get<1>(featurefamily);

    for (auto const &ent1 : parameters) // loop through all parameters for the featureFamily
    {
      auto const &outer_key = ent1.first;
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
          if (currentValue.find(".") != std::string::npos) // this means that the distance is float
          {
            m_Radius_float = std::atof(currentValue.c_str());
            m_Radius = -1;
          }
          else
          {
            m_Radius = std::atoi(currentValue.c_str());
            m_Radius_float = -1;
          }
        }
        else if (outer_key == ParamsString[Neighborhood])
        {
          m_neighborhood = std::atoi(currentValue.c_str());
        }
        else if (outer_key == ParamsString[Bins])
        {
          m_Bins = std::atoi(currentValue.c_str());
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
          m_Range = std::atoi(currentValue.c_str());
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
        else if (outer_key == ParamsString[QuantizationType])
        {
          m_QuantizationType = currentValue;
        }
        else if (outer_key == ParamsString[LBPStyle])
        {
          m_LBPStyle = std::atoi(currentValue.c_str());
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
      exit(EXIT_FAILURE);
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
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    for (size_t i = 0; i < roi.size(); i++)
    {
      m_roi.push_back(std::atoi(roi[i].c_str()));
      m_roiLabels.push_back(roi_labels[i]);
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
    for (const auto &f : selected_features)
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
  m_algorithmDone = false;
}


template< class TImage >
void FeatureExtraction< TImage >::WriteFeatures(const std::string &modality, const std::string &label, const std::string &featureFamily,
  const std::map< std::string, double > &featureList, const std::string &parameters, typename TImage::IndexType centerIndex, bool featureMapWriteForLattice, float weight)
{
  if (m_outputFile.empty())
  {
    m_outputFile = cbica::createTmpDir() + "featureExtractionOutput.csv";
    m_logger.WriteError("Output file has not been initialized; saving in '" + m_outputFile + "'");
    SetOutputFilename(m_outputFile);
  }
  std::ofstream myfile;

  for (auto const &f : featureList)
  {
    auto roiLabelFeatureFamilyFeature = modality + "_" + label + "_" + featureFamily + "_" + f.first;
    if (std::isnan(f.second) || (f.second != f.second))
    {
      m_logger.Write("NAN DETECTED: " + m_patientID + "_" + roiLabelFeatureFamilyFeature);
      std::cerr << "NAN DETECTED: " << m_patientID + "_" + roiLabelFeatureFamilyFeature + "_" + "CenterIdx_" + m_centerIndexString << "\n";
    }
    if ((std::isinf(f.second)))
    {
      m_logger.Write("INF DETECTED: " + m_patientID + "_" + roiLabelFeatureFamilyFeature);
      std::cerr << "INF DETECTED: " << m_patientID + "_" + roiLabelFeatureFamilyFeature + "_" + "CenterIdx_" + m_centerIndexString << "\n";
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
      if (m_outputVerticallyConcatenated)
      {
        if (!cbica::isFile(m_outputFile)) // if file is not present, write the CSV headers 
        {
          myfile.open(m_outputFile, std::ios_base::app);
          // check for locks in a cluster environment
          while (!myfile.is_open())
          {
            cbica::sleep(100);
            myfile.open(m_outputFile, std::ios_base::out | std::ios_base::app);
          }
          myfile << "SubjectID_Modality_ROILabel_FeatureFamily_Feature,Value,Parameters\n";
#ifndef WIN32
          myfile.flush();
#endif
          myfile.close();
        }
        //else // otherwise, append
        {
          myfile.open(m_outputFile, std::ofstream::out | std::ofstream::app);
          // check for locks in a cluster environment
          while (!myfile.is_open())
          {
            cbica::sleep(100);
            myfile.open(m_outputFile, std::ios_base::out | std::ios_base::app);
          }
          myfile << m_patientID + "_" + modality + "_" + label + "_" + featureFamily + "_" + f.first +
            "," + cbica::to_string_precision(f.second) + "," + parameters + "\n";
        }
#ifndef WIN32
        myfile.flush();
#endif
        myfile.close();
      }

      // for training file, populate these 2 member variables
      m_trainingFile_featureNames += roiLabelFeatureFamilyFeature + ",";
      m_trainingFile_features += cbica::to_string_precision(f.second) + ",";
    }
  }
  if (m_outputVerticallyConcatenated)
  {
#ifndef WIN32
    myfile.flush();
#endif
    myfile.close();
  }
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
void FeatureExtraction< TImage >::SetNewLogFile(const std::string &logFile)
{
  m_logger.UseNewFile(logFile);
}


template< class TImage >
void FeatureExtraction< TImage >::SetInputImages(std::vector< typename TImage::Pointer > images, std::vector< std::string >& modality)
{
  if (images.size() != modality.size())
  {
    m_logger.Write("Number of Images and number of modalities are not same");
    exit(EXIT_FAILURE);
  }
  m_inputImages = images;
  m_modality = modality;
  m_algorithmDone = false;
}


template< class TImage >
typename TImage::Pointer FeatureExtraction< TImage >::GetSelectedSlice(typename TImage::Pointer mask, std::string axis)
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
          //std::cout << "[DEBUG] FeatureExtraction.hxx::GetSelectedSlice - slice [" << i << "] -> maxVoxels[" << dim << "] = " << currentNonZero << std::endl;
          //std::cout << "[DEBUG] FeatureExtraction.hxx::GetSelectedSlice - maxVoxels[" << dim << "] = " << maxVoxels[dim] << std::endl;
          //std::cout << "[DEBUG] FeatureExtraction.hxx::GetSelectedSlice - desiredIndexFinal[" << dim << "] = " << desiredIndexFinal[dim] << std::endl;
          // auto duplicator = itk::ImageDuplicator< ImageType2D >::New();
          // duplicator->SetInputImage(extractor->GetOutput());
          // duplicator->Update();
          // maxImageSlices[dim] = duplicator->GetOutput(); // TBD: duplicator->GetOutput() changed to fix compilation issue using vcpkg-cbica
        }
      }
    }
  }

  if (originalSize.Dimension == 3)
  {
    for (int dim = 0; dim < TImage::ImageDimension; dim++) // dimension-loop
    {
      //typename TImage::RegionType desiredRegion;
      //typename TImage::SizeType desiredSize = originalSize;
      ////desiredSize[dim] = 0;
      //typename TImage::IndexType desiredIndex;
      ////desiredIndex.Fill(0);
      ////desiredIndex[dim] = desiredIndexFinal[dim];
      ////desiredIndex = desiredIndexFinal;
      //desiredIndex.Fill(0);
      //desiredIndex[dim] = desiredIndexFinal[dim];

      //std::cout << "[DEBUG] FeatureExtraction.hxx::GetSelectedSlice - maxVoxels[" << dim << "] = " << maxVoxels[dim] << std::endl;
      //std::cout << "[DEBUG] FeatureExtraction.hxx::GetSelectedSlice - desiredIndex[" << dim << "] = " << desiredIndex[dim] << std::endl;
      //std::cout << "[DEBUG] FeatureExtraction.hxx::GetSelectedSlice - axis [" << dim << "] - searching for max region slice. - slice = " << desiredIndex[dim] << std::endl;
      //desiredRegion.SetIndex(desiredIndex);
      //std::cout << "[DEBUG] FeatureExtraction.hxx::GetSelectedSlice - SetIndex(" << desiredIndex << ")" << std::endl;
      //desiredRegion.SetSize(desiredSize);
      //std::cout << "[DEBUG] FeatureExtraction.hxx::GetSelectedSlice - SetSize(" << desiredSize << ")" << std::endl;
      ////auto extractor = itk::ExtractImageFilter< TImage, ImageType2D >::New();
      //auto extractor = itk::ExtractImageFilter< TImage, TImage >::New();
      //extractor->SetInput(mask);
      //std::cout << "[DEBUG] FeatureExtraction.hxx::GetSelectedSlice - SetInput Finished." << std::endl;
      //extractor->SetDirectionCollapseToIdentity();
      //std::cout << "[DEBUG] FeatureExtraction.hxx::GetSelectedSlice - SetDirectionCollapseToIdentity Finished." << std::endl;
      //extractor->SetExtractionRegion(desiredRegion);
      //std::cout << "[DEBUG] FeatureExtraction.hxx::GetSelectedSlice - SetExtractionRegion Finished." << std::endl;
      //extractor->Update();
      //std::cout << "[DEBUG] FeatureExtraction.hxx::GetSelectedSlice - Update Finished." << std::endl;

      //maxImageSlices[dim] = extractor->GetOutput();

      itk::ImageRegionIteratorWithIndex< TImage > iteratorExtractor(mask, mask->GetLargestPossibleRegion());
      itk::ImageRegionIteratorWithIndex< TImage > iteratormaxImageSlices(maxImageSlices[dim], maxImageSlices[dim]->GetLargestPossibleRegion());
      // itk::ImageRegionConstIterator< ImageType2D > iterator(extractor->GetOutput(), extractor->GetOutput()->GetLargestPossibleRegion());
      for (iteratormaxImageSlices.GoToBegin(); !iteratormaxImageSlices.IsAtEnd(); ++iteratormaxImageSlices)
      {
        auto idx = iteratormaxImageSlices.GetIndex();
        if (idx[dim] == desiredIndexFinal[dim]) {
          iteratorExtractor.SetIndex(idx);
          //std::cout << "[DEBUG] FeatureExtraction.hxx::GetSelectedSlice - index[" << idx << "] = " << iteratorExtractor.Get() << std::endl;
          iteratormaxImageSlices.Set(iteratorExtractor.Get());
        }
        //if (iteratorExtractor.Get() > 0) {
        //  std::cout << "[DEBUG] FeatureExtraction.hxx::GetSelectedSlice - index[" << idx << "] = " << iteratorExtractor.Get() << std::endl;
        //  iteratormaxImageSlices.Set(iteratorExtractor.Get());
        //}
        //std::cout << " - iteratormaxImageSlices - index[" << idx << "] = " << iteratormaxImageSlices.Get() << std::endl;
      }
    }
  }

  //return the correct element of the 3 element vector depending on axis defined by input
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
void FeatureExtraction< TImage >::Update()
{
  if (!m_algorithmDone)
  {
    auto t1 = std::chrono::high_resolution_clock::now();

    if (m_debug)
    {
      m_logger.Write("Checking mask validity (whether it is empty or not)");
    }

    if (!m_maskValidated)
    {
      TConstIteratorType maskIt(m_Mask, m_Mask->GetBufferedRegion());
      if (m_roi.size() != 0)
      {
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
              std::cerr << "The ROI for calculation, '" << std::to_string(m_roi[x]) << "' does not exist in the mask.\n";
              exit(EXIT_FAILURE);
            }
          }
        }
      }
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
          exit(EXIT_FAILURE);
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
        auto temp = m_Features.find(FeatureFamilyString[Quantization]);
        if (temp != m_Features.end())
        {
          if (std::get<0>(temp->second)) // if the feature family has been selected in the GUI
          {
            m_LatticeComputation = true;
            SetFeatureParam(FeatureFamilyString[Quantization]);
          }
        }
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
        if (!m_patchFullImageComputation)
        {
          j += m_roi.size();
        }
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
            m_logger.Write("Percentage done: " + std::to_string(temp * 100));
            //auto totalMemory = static_cast< float >(cbica::getTotalMemory()) / 10e8;
            //auto usedMemory = static_cast< float >(cbica::getCurrentlyUsedMemory());
            //auto freeMem = totalMemory - usedMemory;
            //m_logger.Write("Approximate free memory on machine: '" + std::to_string(freeMem / 10e8) + "' Gb.'");
          }
          else
          {
            m_logger.Write("Calculating Features for modality '" + m_modality[i] + "' and ROI '" + allROIs[j].label + "'");
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

            if (m_QuantizationType == "Image")
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
            std::get<2>(temp->second) = m_modality[i];
            std::get<3>(temp->second) = allROIs[j].label;
            CalculateIntensity(m_currentNonZeroImageValues, std::get<4>(temp->second), allROIs[j].latticeGridPoint);
            WriteFeatures(m_modality[i], allROIs[j].label, FeatureFamilyString[Intensity], std::get<4>(temp->second), "N.A.", m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

            if (m_debug)
            {
              auto tempT2 = std::chrono::high_resolution_clock::now();
              m_logger.Write("Intensity Features for modality '" + m_modality[i] + "' and ROI '" + allROIs[j].label + "' calculated in " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(tempT2 - tempT1).count()) + " milliseconds");
            }
          }

          // iterate over the entire feature family enum
          for (size_t f = 1/*Intensity features already calculated*/; f < FeatureMax; ++f)
          {
            SetFeatureParam(FeatureFamilyString[f]);
            switch (f)
            {
              // case Intensity is not needed since it always calculated
            case Histogram:
            {
              auto temp = m_Features.find(FeatureFamilyString[f]);
              if (temp != m_Features.end())
              {
                if (std::get<0>(temp->second))
                {
                  auto tempT1 = std::chrono::high_resolution_clock::now();

                  //auto local_map = std::get<1>(temp->second);
                  std::get<2>(temp->second) = m_modality[i];
                  std::get<3>(temp->second) = allROIs[j].label;
                  CalculateHistogram(currentInputImage_patch, currentMask_patch, std::get<4>(temp->second), allROIs[j].latticeGridPoint);
                  WriteFeatures(m_modality[i], allROIs[j].label, FeatureFamilyString[f], std::get<4>(temp->second),
                    "Bins=" + std::to_string(m_Bins), m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

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

                    std::get<2>(temp->second) = m_modality[i];
                    std::get<3>(temp->second) = allROIs[j].label;

                    /* this dimensionality reduction applies only to shape and Volumetric features */
                    if (TImage::ImageDimension == 3)
                    {
                      if (m_Dimension == 2) // extracts slice with maximum area along the specified axis
                      {
                        //std::cout << "[DEBUG] FeatureExtraction.hxx - calling GetSelectedSlice" << std::endl;
                        auto selected_axis_image = GetSelectedSlice(currentMask_patch, m_Axis);
                        //std::cout << "[DEBUG] FeatureExtraction.hxx - called GetSelectedSlice" << std::endl;
                        //CalculateMorphologic<ImageType2D>(currentInputImage_patch, selected_axis_image, currentMask_patch, std::get<4>(temp->second)); //old with 2D
                        CalculateMorphologic<TImage>(currentInputImage_patch, selected_axis_image, currentMask_patch, std::get<4>(temp->second));
                      }
                      else
                      {
                        CalculateMorphologic<TImage>(currentInputImage_patch, currentMask_patch, currentMask_patch, std::get<4>(temp->second));
                      }
                    }
                    else
                    {
                      CalculateMorphologic<TImage>(currentInputImage_patch, currentMask_patch, currentMask_patch, std::get<4>(temp->second));
                    }
                    WriteFeatures(m_modality[i], allROIs[j].label, FeatureFamilyString[f], std::get<4>(temp->second),
                      "Axis=" + m_Axis + ";Dimension=" + std::to_string(m_Dimension), m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

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

                    std::get<2>(temp->second) = m_modality[i];
                    std::get<3>(temp->second) = allROIs[j].label;

                    /* this dimensionality reduction applies only to shape and Volumetric features */
                    if (TImage::ImageDimension == 3)
                    {
                      if (m_Dimension == 2)
                      {
                        //ImageType2D::Pointer selected_axis_image = GetSelectedSlice(currentMask_patch, m_Axis);
                        //CalculateVolumetric<ImageType2D>(selected_axis_image, std::get<4>(temp->second));
                        auto selected_axis_image = GetSelectedSlice(currentMask_patch, m_Axis);
                        CalculateVolumetric<TImage>(selected_axis_image, std::get<4>(temp->second));
                      }
                      else
                      {
                        CalculateVolumetric<TImage>(currentMask_patch, std::get<4>(temp->second));
                      }
                    }
                    else
                    {
                      CalculateVolumetric<TImage>(currentMask_patch, std::get<4>(temp->second));
                    }
                    WriteFeatures(m_modality[i], allROIs[j].label, FeatureFamilyString[f], std::get<4>(temp->second),
                      "Axis=" + m_Axis + ";Dimension=" + std::to_string(m_Dimension), m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

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

                  auto offsets = GetOffsetVector(m_Radius, m_Direction);
                  /* this dimensionality reduction applies only to shape and Volumetric features */
                  if (TImage::ImageDimension == 3)
                  {
                    if (m_Dimension == 2) // extracts slice with maximum area along the specified axis
                    {
                      //std::cout << "[DEBUG] FeatureExtraction.hxx - calling GetSelectedSlice" << std::endl;
                      auto selected_axis_image = GetSelectedSlice(currentMask_patch, m_Axis);
                      offsets = GetOffsetVector(m_Radius, 26); // because anything other than 26 doesn't work properly for GLSZM computation
                      //cbica::WriteImage< TImage >(selected_axis_image, "tmp.nii.gz");
                      //std::cout << "[DEBUG] FeatureExtraction.hxx - called GetSelectedSlice" << std::endl;
                      //CalculateMorphologic<ImageType2D>(currentInputImage_patch, selected_axis_image, currentMask_patch, std::get<4>(temp->second)); //old with 2D
                      //CalculateMorphologic<TImage>(currentInputImage_patch, selected_axis_image, currentMask_patch, std::get<4>(temp->second));
                      CalculateGLCM(currentInputImage_patch, selected_axis_image, offsets, std::get<4>(temp->second), allROIs[j].latticeGridPoint);
                    }
                    else
                    {
                      CalculateGLCM(currentInputImage_patch, currentMask_patch, offsets, std::get<4>(temp->second), allROIs[j].latticeGridPoint);
                    }
                  }
                  else
                  {
                    CalculateGLCM(currentInputImage_patch, currentMask_patch, offsets, std::get<4>(temp->second), allROIs[j].latticeGridPoint);
                  }

                  WriteFeatures(m_modality[i], allROIs[j].label, FeatureFamilyString[f], std::get<4>(temp->second),
                    "Axis=" + m_Axis + ";Dimension=" + std::to_string(m_Dimension) + ";Bins=" + std::to_string(m_Bins) + ";Directions=" + std::to_string(m_Direction) +
                    ";Radius=" + std::to_string(m_Radius) + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice);

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

                  auto offsets = GetOffsetVector(m_Radius, m_Direction);
                  /* this dimensionality reduction applies only to shape and Volumetric features */
                  if (TImage::ImageDimension == 3)
                  {
                    if (m_Dimension == 2) // extracts slice with maximum area along the specified axis
                    {
                      //std::cout << "[DEBUG] FeatureExtraction.hxx - calling GetSelectedSlice" << std::endl;
                      auto selected_axis_image = GetSelectedSlice(currentMask_patch, m_Axis);
                      offsets = GetOffsetVector(m_Radius, 26); // because anything other than 26 doesn't work properly for GLSZM computation
                      //cbica::WriteImage< TImage >(selected_axis_image, "tmp.nii.gz");
                      //std::cout << "[DEBUG] FeatureExtraction.hxx - called GetSelectedSlice" << std::endl;
                      //CalculateMorphologic<ImageType2D>(currentInputImage_patch, selected_axis_image, currentMask_patch, std::get<4>(temp->second)); //old with 2D
                      //CalculateMorphologic<TImage>(currentInputImage_patch, selected_axis_image, currentMask_patch, std::get<4>(temp->second));
                      CalculateGLRLM(currentInputImage_patch, selected_axis_image, offsets, std::get<4>(temp->second), allROIs[j].latticeGridPoint);
                    }
                    else
                    {
                      CalculateGLRLM(currentInputImage_patch, currentMask_patch, offsets, std::get<4>(temp->second), allROIs[j].latticeGridPoint);
                    }
                  }
                  else
                  {
                    CalculateGLRLM(currentInputImage_patch, currentMask_patch, offsets, std::get<4>(temp->second), allROIs[j].latticeGridPoint);
                  }

                  WriteFeatures(m_modality[i], allROIs[j].label, FeatureFamilyString[f], std::get<4>(temp->second),
                    "Axis=" + m_Axis + ";Dimension=" + std::to_string(m_Dimension) + ";Bins=" + std::to_string(m_Bins) + ";Directions=" + std::to_string(m_Direction) +
                    ";Radius=" + std::to_string(m_Radius) + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

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

                  auto offsets = GetOffsetVector(m_Radius, m_Direction);
                  /* this dimensionality reduction applies only to shape and Volumetric features */
                  if (TImage::ImageDimension == 3)
                  {
                    if (m_Dimension == 2) // extracts slice with maximum area along the specified axis
                    {
                      //std::cout << "[DEBUG] FeatureExtraction.hxx - calling GetSelectedSlice" << std::endl;
                      offsets = GetOffsetVector(m_Radius, 26);
                      auto selected_axis_image = GetSelectedSlice(currentMask_patch, m_Axis);
                      offsets = GetOffsetVector(m_Radius, 26); // because anything other than 26 doesn't work properly for GLSZM computation
                      //cbica::WriteImage< TImage >(selected_axis_image, "tmp.nii.gz");
                      //std::cout << "[DEBUG] FeatureExtraction.hxx - called GetSelectedSlice" << std::endl;
                      //CalculateMorphologic<ImageType2D>(currentInputImage_patch, selected_axis_image, currentMask_patch, std::get<4>(temp->second)); //old with 2D
                      //CalculateMorphologic<TImage>(currentInputImage_patch, selected_axis_image, currentMask_patch, std::get<4>(temp->second));
                      CalculateGLSZM(currentInputImage_patch, selected_axis_image, offsets, std::get<4>(temp->second));
                    }
                    else
                    {
                      CalculateGLSZM(currentInputImage_patch, currentMask_patch, offsets, std::get<4>(temp->second));
                    }
                  }
                  else
                  {
                    CalculateGLSZM(currentInputImage_patch, currentMask_patch, offsets, std::get<4>(temp->second));
                  }

                  WriteFeatures(m_modality[i], allROIs[j].label, FeatureFamilyString[f], std::get<4>(temp->second),
                    "Axis=" + m_Axis + ";Dimension=" + std::to_string(m_Dimension) + ";Bins=" + std::to_string(m_Bins) + ";Directions=" + std::to_string(m_Direction) +
                    ";Radius=" + std::to_string(m_Radius) + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

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
              auto temp = m_Features.find(FeatureFamilyString[f]);
              if (temp != m_Features.end())
              {
                if (std::get<0>(temp->second))
                {
                  auto tempT1 = std::chrono::high_resolution_clock::now();

                  std::get<2>(temp->second) = m_modality[i];
                  std::get<3>(temp->second) = allROIs[j].label;

                  auto offsets = GetOffsetVector(m_Radius, m_Direction);
                  /* this dimensionality reduction applies only to shape and Volumetric features */
                  if (TImage::ImageDimension == 3)
                  {
                    if (m_Dimension == 2) // extracts slice with maximum area along the specified axis
                    {
                      //std::cout << "[DEBUG] FeatureExtraction.hxx - calling GetSelectedSlice" << std::endl;
                      auto selected_axis_image = GetSelectedSlice(currentMask_patch, m_Axis);
                      offsets = GetOffsetVector(m_Radius, 26); // because anything other than 26 doesn't work properly for GLSZM computation
                      //cbica::WriteImage< TImage >(selected_axis_image, "tmp.nii.gz");
                      //std::cout << "[DEBUG] FeatureExtraction.hxx - called GetSelectedSlice" << std::endl;
                      //CalculateMorphologic<ImageType2D>(currentInputImage_patch, selected_axis_image, currentMask_patch, std::get<4>(temp->second)); //old with 2D
                      //CalculateMorphologic<TImage>(currentInputImage_patch, selected_axis_image, currentMask_patch, std::get<4>(temp->second));
                      CalculateNGTDM(currentInputImage_patch, selected_axis_image, offsets, std::get<4>(temp->second));
                    }
                    else
                    {
                      CalculateNGTDM(currentInputImage_patch, currentMask_patch, offsets, std::get<4>(temp->second));
                    }
                  }
                  else
                  {
                    CalculateNGTDM(currentInputImage_patch, currentMask_patch, offsets, std::get<4>(temp->second));
                  }

                  WriteFeatures(m_modality[i], allROIs[j].label, FeatureFamilyString[f], std::get<4>(temp->second),
                    "Axis=" + m_Axis + ";Dimension=" + std::to_string(m_Dimension) + ";Bins=" + std::to_string(m_Bins) + ";Directions=" + std::to_string(m_Direction) +
                    ";Radius=" + std::to_string(m_Radius) + ";OffsetType=" + m_offsetSelect, m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

                  if (m_debug)
                  {
                    auto tempT2 = std::chrono::high_resolution_clock::now();
                    m_logger.Write("NGTDM Features for modality '" + m_modality[i] + "' and ROI '" + allROIs[j].label + "' calculated in " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(tempT2 - tempT1).count()) + " milliseconds");
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

                      CalculateGaborWavelets(currentInputImage_patch, std::get<4>(temp->second), allROIs[j].latticeGridPoint);
                      WriteFeatures(m_modality[i], allROIs[j].label, FeatureFamilyString[f], std::get<4>(temp->second),
                        "Radius=" + std::to_string(m_Radius) + ";FMax=" + std::to_string(m_gaborFMax) + ";Gamma=" + std::to_string(m_gaborGamma) +
                        ";Directions=" + std::to_string(m_Direction) + ";Level=" + std::to_string(m_gaborLevel), m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

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

                      CalculateEdgeEnhancement(currentInputImage_patch, std::get<4>(temp->second));
                      WriteFeatures(m_modality[i], allROIs[j].label, FeatureFamilyString[f], std::get<4>(temp->second),
                        "ETA=" + std::to_string(m_edgesETA) + ";Epsilon=" + std::to_string(m_edgesEpsilon) + ";Radius=" + std::to_string(m_Radius),
                        m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

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

                  CalculateLBP(currentInputImage_patch, currentMask_patch, std::get<4>(temp->second));
                  WriteFeatures(m_modality[i], allROIs[j].label, FeatureFamilyString[f], std::get<4>(temp->second),
                    "Neighborhood=" + std::to_string(m_neighborhood) + ";Radius=" + std::to_string(m_Radius) + ";Style=" + std::to_string(m_LBPStyle), m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);

                  if (m_debug)
                  {
                    auto tempT2 = std::chrono::high_resolution_clock::now();
                    m_logger.Write("GLRLM Features for modality '" + m_modality[i] + "' and ROI '" + allROIs[j].label + "' calculated in " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(tempT2 - tempT1).count()) + " milliseconds");
                  }


                  /////// OLD CODE

                  // do LBP calculation -- which is NOT what is being done below

                  /* TBD this portion needs to be changed by code from Saima
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

                  auto padFilter = itk::ConstantPadImageFilter < TImage, TImage >::New();
                  padFilter->SetInput(currentMask);
                  //padFilter->SetPadBound(outputRegion); // Calls SetPadLowerBound(region) and SetPadUpperBound(region)
                  padFilter->SetPadLowerBound(lowerExtendRegion);
                  padFilter->SetPadUpperBound(upperExtendRegion);
                  padFilter->SetConstant(constantPixel);
                  padFilter->Update();
                  typename TImage::Pointer  lbproi = padFilter->GetOutput();
                  //const typename  TImage::SizeType image_size = currentMask->GetLargestPossibleRegion().GetSize();
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
                  //  pixelValue = currentMask->GetPixel(ind1);
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

                  //lbpfeatures.calculateLBP<TImage>(currentInputImage_patch, currentMask_patch, lbproi, m_Radius, m_neighborhood, m_modality[i], std::to_string(m_roi[j]), std::get<4>(temp->second));
                  //WriteFeatures(m_modality[i],allROIs[j].label, FeatureFamilyString[f], std::get<4>(temp->second), m_currentLatticeCenter, writeFeatureMapsAndLattice, allROIs[j].weight);
                  if (m_debug)
                  {
                  m_logger.Write("LBP Features for modality '" + m_modality[i] + "' and ROI '" +allROIs[j].label + "' calculated");
                  }
                  */
                }
              }
              break;
            }
            default: // undefined Feature
              break;
            }
          } // end of feature iteration   
        } // End of image iteration loop

      } // end of allROIs iteration

      // calculate 1st order statistics of the different lattice features
      for (auto const &entry : m_LatticeFeatures)
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
          std::ofstream myfile;
          myfile.open(m_outputFile, std::ofstream::out | std::ofstream::app);

          myfile << currentPatientModalityROIFeatureFamilyFeature + "Max" +
            "," + currentMax + "," + "\n";
          myfile << currentPatientModalityROIFeatureFamilyFeature + "Min" +
            "," + currentMin + "," + "\n";
          myfile << currentPatientModalityROIFeatureFamilyFeature + "Variance" +
            "," + currentVar + "," + "\n";
          myfile << currentPatientModalityROIFeatureFamilyFeature + "StdDev" +
            "," + currentStdDev + "," + "\n";
          myfile << currentPatientModalityROIFeatureFamilyFeature + "Skewness" +
            "," + currentSkew + "," + "\n";
          myfile << currentPatientModalityROIFeatureFamilyFeature + "Kurtosis" +
            "," + currentKurt + "," + "\n";
          myfile << currentPatientModalityROIFeatureFamilyFeature + "Mean" +
            "," + currentMean + "," + "\n";
          myfile << currentPatientModalityROIFeatureFamilyFeature + "Median" +
            "," + currentMedian + "," + "\n";
        }
        else
        {
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
      if (!m_outputVerticallyConcatenated)
      {
        m_trainingFile_featureNames.pop_back(); // since the last character is always a ","
        m_trainingFile_features.pop_back(); // since the last character is always a ","
        auto featureNamesVec = cbica::stringSplit(m_trainingFile_featureNames, ",");
        auto featureVec = cbica::stringSplit(m_trainingFile_features, ",");
        if (featureNamesVec.size() != featureVec.size())
        {
          m_logger.WriteError("Something went wrong and the featureNames (" + std::to_string(featureNamesVec.size()) +
            ") and featureVec (" + std::to_string(featureVec.size()) + ") is not of same size.");
          exit(EXIT_FAILURE);
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

        for (auto const &entry : m_downscaledFeatureMaps)
        {
          auto currentDownscaledFileName = entry.first;
          cbica::WriteImage< TImage >(entry.second, m_outputPath + "/" + currentDownscaledFileName + ".nii.gz");
        }
      }

      m_algorithmDone = true;

      auto t2 = std::chrono::high_resolution_clock::now();
      std::cout << "FE took " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " milliseconds\n";
    }
  }
}
