#pragma once
#include <iostream>
#include <iterator>
#include <stdlib.h>
#include <string>
#include <algorithm>
#include <exception>
#include <typeinfo>
#include <stdexcept>
#include <stdio.h>
#include <numeric>
#include <functional>
#include <vector> 




#include "cbicaCmdParser.h"

#include "cbicaUtilities.h"
#include "itkImageFileReader.h"
#include "itkScalarImageToRunLengthFeaturesFilter.h"
#include "itkHistogramToTextureFeaturesFilter.h"
#include <itkVectorImage.h>
#include "itkScalarImageToCooccurrenceMatrixFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVectorContainer.h"
#include "itkStatisticsImageFilter.h"
#include "itkImageToHistogramFilter.h"
#include <itkSliceBySliceImageFilter.h>
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"

const unsigned int  Dimension = 3;
typedef  float           PixelType;
typedef itk::Image< float, 3 > ImageTypeFloat3D;
typedef ImageTypeFloat3D::OffsetType OffsetType;
typedef itk::VectorContainer< unsigned char, OffsetType > OffsetVector;
typedef OffsetVector::Pointer  OffsetVectorPointer;

typedef itk::Statistics::ScalarImageToCooccurrenceMatrixFilter<ImageTypeFloat3D>Image2CoOccuranceType;
typedef Image2CoOccuranceType::HistogramType HistogramType;
typedef itk::Statistics::HistogramToTextureFeaturesFilter<HistogramType> Hist2FeaturesType1;

typedef itk::Statistics::DenseFrequencyContainer2 HistogramFrequencyContainerType;
typedef itk::Statistics::ScalarImageToRunLengthFeaturesFilter
<ImageTypeFloat3D, HistogramFrequencyContainerType> RunLengthFilterType;


typedef short                                    RunLengthFeatureName;
typedef itk::VectorContainer<unsigned char, RunLengthFeatureName>  FeatureNameVector;

typedef  short                                LabelType;
typedef itk::Image< LabelType, Dimension >            OutputImageType;
typedef itk::ShapeLabelObject< LabelType, 3 > ShapeLabelObjectType;
typedef itk::LabelMap< ShapeLabelObjectType >         LabelMapType;

class generateTextureFeatures
{
public:
  typedef enum
  {
    Brain,
    Breast,
    Lung,
    Generic
  }Organ_type;

  
  generateTextureFeatures();
  ~generateTextureFeatures();

  template <class TImageType = ImageTypeFloat3D>
  void setInputs(int type, std::vector<typename TImageType::Pointer>* Input_image, typename TImageType::Pointer mask);


  template <class TImageType = ImageTypeFloat3D, typename Offsets = OffsetVectorPointer>
  void calculateTextureFeatures(typename TImageType::Pointer image, typename TImageType::Pointer mask, int min, int max, Offsets offset, std::vector<std::string> &glcm_featurename, std::vector<double>& glcm_feature);

  template <class TImageType = ImageTypeFloat3D, typename Offsets = OffsetVector >
  void calculateRunLength(typename TImageType::Pointer image, typename TImageType::Pointer mask, int min, int max,
    Offsets offset, std::vector<std::string>& runLength_featurename, std::vector<double>& runLength_feature);

  template <class TImageType = ImageTypeFloat3D >
  void ShapeFeatures(typename TImageType::Pointer mask, std::vector<std::string>& runLength_featurename, std::vector<double>& runLength_feature);

  template <class TImageType = ImageTypeFloat3D>
  void Intensity_featurs(int type, typename TImageType::Pointer image, typename TImageType::Pointer mask);


  template <class TImageType = ImageTypeFloat3D>
  void set_organ_type(int type, std::vector<typename TImageType::Pointer>& Input_image, typename TImageType::Pointer mask);


private:
  double minimum_intensity, maximum_intensity;
  OffsetVectorPointer offsets = OffsetVector::New();

  //** Add all known modality of images that will be used by different applications 
  ImageTypeFloat3D::Pointer image;  ImageTypeFloat3D::Pointer mask;
  ImageTypeFloat3D::Pointer Pet_image;

  std::vector<std::string> *Featurename;
  std::vector<std::string> *Featurevec;

  //these could be user changeable parameters
  int hist_bin, num_neighbours;


};


/**
\calculating all 8 Glcm features from given image and mask
*/
template <class TImageType, typename Offsets>
void generateTextureFeatures::calculateTextureFeatures(typename TImageType::Pointer inputImage, typename TImageType::Pointer mask, int min, int max, Offsets offset, std::vector<std::string> &glcm_featurename,
  std::vector<double>& glcm_feature)
{

  typename TImageType::Pointer inertia = TImageType::New();
  typename TImageType::Pointer correlation = TImageType::New();
  typename TImageType::Pointer energy = TImageType::New();
  typename TImageType::Pointer entropy = TImageType::New();
  inertia->SetRegions(inputImage->GetLargestPossibleRegion());
  correlation->SetRegions(inputImage->GetLargestPossibleRegion());
  energy->SetRegions(inputImage->GetLargestPossibleRegion());
  entropy->SetRegions(inputImage->GetLargestPossibleRegion());
  //	typedef itk::VectorContainer< unsigned char, OffsetType > OffsetVector;
  OffsetVectorPointer offsets ;
  offsets = offset;



  for (int i = 0; i < offsets->Size(); i++)
  {
    Image2CoOccuranceType::Pointer glcmGenerator = Image2CoOccuranceType::New();
    glcmGenerator->SetOffset(offsets->at(i));
    glcmGenerator->SetNumberOfBinsPerAxis(16); //reasonable number of bins
    glcmGenerator->SetPixelValueMinMax(min, max); //for input UCHAR pixel type
    Hist2FeaturesType1::Pointer featureCalc = Hist2FeaturesType1::New();
    glcmGenerator->SetInput(mask);
    glcmGenerator->Update();
    featureCalc->SetInput(glcmGenerator->GetOutput());
    featureCalc->Update();


    std::stringstream ss;
    ss << offsets->GetElement(i);
    std::string current_offset = ss.str();
   
    glcm_featurename.push_back("correlation" + current_offset);
    glcm_featurename.push_back("energy" + current_offset);
    glcm_featurename.push_back("contrast" + current_offset);
    glcm_featurename.push_back("entropy" + current_offset);
    //glcm_featurename.push_back("inversediffmom" + current_offset);
    //glcm_featurename.push_back("Clustershade" + current_offset);
    //glcm_featurename.push_back("ClusterProminence" + current_offset);
    //glcm_featurename.push_back("haarlickCorrelation" + current_offset);

    glcm_feature.push_back(featureCalc->GetFeature(Hist2FeaturesType1::Correlation));
    glcm_feature.push_back(featureCalc->GetFeature(Hist2FeaturesType1::Energy));
    glcm_feature.push_back(featureCalc->GetFeature(Hist2FeaturesType1::Inertia));
    glcm_feature.push_back(featureCalc->GetFeature(Hist2FeaturesType1::Entropy));
  //  glcm_feature.push_back(featureCalc->GetFeature(Hist2FeaturesType1::InverseDifferenceMoment));
  //  glcm_feature.push_back(featureCalc->GetFeature(Hist2FeaturesType1::ClusterShade));
 //   glcm_feature.push_back(featureCalc->GetFeature(Hist2FeaturesType1::ClusterProminence));
   // glcm_feature.push_back(featureCalc->GetFeature(Hist2FeaturesType1::HaralickCorrelation));

  }

}


/**
\calculating Runlength features from given image and mask
*/
template <class TImageType, typename Offsets >
void generateTextureFeatures::calculateRunLength(typename TImageType::Pointer image, typename  TImageType::Pointer mask, int min, int max,
  Offsets offset, std::vector<std::string>& runLength_featurename, std::vector<double>& runLength_feature)
{
  RunLengthFilterType::Pointer texFilter = RunLengthFilterType::New();
  OffsetVectorPointer offsets;
  offsets = offset;

  typedef  RunLengthFilterType::RunLengthMatrixFilterType RunLengthMatrixGenerator;
  typedef  RunLengthFilterType::RunLengthFeaturesFilterType RunLengthFeatures;
  typename  RunLengthMatrixGenerator::Pointer matrix_generator = RunLengthMatrixGenerator::New();
  typedef RunLengthFilterType::RunLengthFeaturesFilterType    RunLengthFeatureName;
  typedef RunLengthFeatureName::RunLengthFeatureName InternalRunLengthFeatureName;
  RunLengthFilterType::FeatureNameVectorPointer requestedFeatures = RunLengthFilterType::FeatureNameVector::New();

  requestedFeatures->push_back(RunLengthFeatureName::ShortRunEmphasis);
  requestedFeatures->push_back(RunLengthFeatureName::LongRunEmphasis);
  requestedFeatures->push_back(
    RunLengthFeatureName::GreyLevelNonuniformity);
  requestedFeatures->push_back(
    RunLengthFeatureName::RunLengthNonuniformity);
  requestedFeatures->push_back(
    RunLengthFeatureName::LowGreyLevelRunEmphasis);
  requestedFeatures->push_back(
    RunLengthFeatureName::HighGreyLevelRunEmphasis);
  requestedFeatures->push_back(
    RunLengthFeatureName::ShortRunLowGreyLevelEmphasis);
  requestedFeatures->push_back(
    RunLengthFeatureName::ShortRunHighGreyLevelEmphasis);
  requestedFeatures->push_back(
    RunLengthFeatureName::LongRunLowGreyLevelEmphasis);
  requestedFeatures->push_back(
    RunLengthFeatureName::LongRunHighGreyLevelEmphasis);
  std::vector <std::string> featurename;
  featurename.push_back("SRE"); featurename.push_back("LRE");
  featurename.push_back("GLN"); featurename.push_back("RLN");
  featurename.push_back("LGRE"); featurename.push_back("HGRE");
  featurename.push_back("SRLGE");  featurename.push_back("SRHGE");
  featurename.push_back("LRLGE");  featurename.push_back("LRHGE");


  OffsetVector::ConstIterator offsetIt;
  size_t offsetNum, featureNum;
  size_t numOffsets = offsets->size();
  size_t numFeatures = requestedFeatures->size();
  double **features;
  features = new double *[numOffsets];

  for (size_t i = 0; i < numOffsets; i++)
  {
    features[i] = new double[numFeatures];
  }

  for (offsetIt = offsets->Begin(), offsetNum = 0;
    offsetIt != offsets->End(); offsetIt++, offsetNum++){

    matrix_generator->SetOffset(offsetIt.Value());
    matrix_generator->SetInput(mask);
    matrix_generator->SetInsidePixelValue(1);
    matrix_generator->SetPixelValueMinMax(min, max);
    matrix_generator->SetDistanceValueMinMax(0, 4);
    matrix_generator->SetNumberOfBinsPerAxis(10);
    matrix_generator->Update();

    std::stringstream ss;

    RunLengthFeatures::Pointer runLengthMatrixCalculator = RunLengthFeatures::New();
    runLengthMatrixCalculator->SetInput(matrix_generator->GetOutput());
    runLengthMatrixCalculator->Update();
    FeatureNameVector::ConstIterator fnameIt;
    for (fnameIt = requestedFeatures->Begin(), featureNum = 0;
      fnameIt != requestedFeatures->End(); fnameIt++, featureNum++)
    {

      ss << offsetIt.Value();
      std::string current_offset = ss.str();
      runLength_featurename.push_back(featurename.at(fnameIt.Value()) + current_offset);
      runLength_feature.push_back(runLengthMatrixCalculator->GetFeature((InternalRunLengthFeatureName)fnameIt.Value()));
      features[offsetNum][featureNum] = runLengthMatrixCalculator->GetFeature((InternalRunLengthFeatureName)fnameIt.Value());
    }
  }
}



/**
\calculating shape features of the given tumoer mask
*/
template <class TImageType  >
void generateTextureFeatures::ShapeFeatures(typename TImageType::Pointer mask, std::vector<std::string>& shape_featurename, std::vector<double>& shape_feature)
{

  typedef itk::ConnectedComponentImageFilter < ImageTypeFloat3D, OutputImageType >
    ConnectedComponentImageFilterType;
  typedef itk::LabelImageToShapeLabelMapFilter < OutputImageType, LabelMapType >
    I2LType;
  ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
  connected->SetInput(mask);
  connected->Update();


  I2LType::Pointer i2l = I2LType::New();
  i2l->SetInput(connected->GetOutput());
  i2l->SetComputePerimeter(true);
  i2l->Update();
  LabelMapType *labelMap = i2l->GetOutput();

  ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(0);

  shape_featurename.push_back("no_of_pix");
  shape_featurename.push_back("principleMomemt_1");
  shape_featurename.push_back("principleMomemt_2");
  shape_featurename.push_back("principleMomemt_3");
  shape_featurename.push_back("Elongation");
  shape_featurename.push_back("Perimeter");
  shape_featurename.push_back("Roundness_3d");
  shape_featurename.push_back("Flatness");
  // shape_featurename.push_back("Roundness_2d");

  shape_feature.push_back(labelObject->GetNumberOfPixels());
  shape_feature.push_back(labelObject->GetPrincipalMoments().operator[](0));
  shape_feature.push_back(labelObject->GetPrincipalMoments().operator[](1));
  shape_feature.push_back(labelObject->GetPrincipalMoments().operator[](2));
  shape_feature.push_back(labelObject->GetElongation());
  shape_feature.push_back(labelObject->GetPerimeter());
  shape_feature.push_back(labelObject->GetRoundness());
  shape_feature.push_back(labelObject->GetFlatness());


  /**********************************/
  // * for calculating the roundness crossing through each slice. Must be modified *******/

  //typedef itk::SliceBySliceImageFilter<ImageType, 2DImageType> SliceImageFilter;
  //SliceImageFilter slice_image = SliceImageFilter::New();
  //slice_image->SetInputFilter(mask);


}