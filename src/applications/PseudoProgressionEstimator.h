/**
\file PseudoProgressionEstimator.h

This file holds the declaration of the class PseudoProgressionEstimator.

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html
*/

#ifndef _PseudoProgressionEstimator_h_
#define _PseudoProgressionEstimator_h_

//#include "CAPTk.h"
#include "NiftiDataManager.h"
#include "OutputWritingManager.h"
#include "FeatureReductionClass.h"
#include "FeatureScalingClass.h"
#include "FeatureExtractionClass.h"
#include "fProgressDialog.h"
#include "itkMeanImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkShapeLabelObject.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkScalarImageToCooccurrenceMatrixFilter.h"
#include "itkHistogramToTextureFeaturesFilter.h"
#include "itkScalarImageToRunLengthFeaturesFilter.h"
#include "itkCSVArray2DFileReader.h"
#include "itkConnectedComponentImageFilter.h"
#include "CaPTkEnums.h"
#include "CaPTkUtils.h"
#include "CaPTkClassifierUtils.h"
#include "cbicaLogging.h"

#ifdef APP_BASE_CAPTK_H
#include "ApplicationBase.h"
#endif
#define NO_OF_PCS 5 // total number of principal components used 

typedef itk::Image< float, 3 > ImageType;
typedef itk::Image< float, 4 > PerfusionImageType;
typedef itk::CSVArray2DFileReader<double> CSVFileReaderType;
typedef vnl_matrix<double> MatrixType;

typedef std::tuple< std::vector<ImageType::IndexType>, VariableSizeMatrixType> PerfusionTupleType;
typedef std::map<int, PerfusionTupleType> PerfusionMapType;

#define PSP_NO_OF_FEATURES 1040
#define TXT_NO_OF_FEATURES 810

/**
\class PseudoProgressionEstimator
\brief A small description of PseudoProgressionEstimator
Give more details about the class (like its usage, default parameters, etc.) along with the correct reference paper.

References:

@inproceedings{,
title={Imaging Surrogates of Infiltration Obtained Via Multiparametric Imaging Pattern Analysis Predict Subsequent Location of Recurrence of Glioblastoma},
author={Akbari, Hamed MD, PhD; Macyszyn, Luke MD, MA; Da, Xiao MSc; Bilello, Michel MD, PhD; Wolf, Ronald L. MD, PhD; Martinez-Lage, Maria MD; Biros, George PhD; Alonso-Basanta, Michelle MD, PhD; O'Rourke, Donald M. MD; Davatzikos, Christos PhD},
pages={572�580},
year={2016},
organization={Neurosurgery}
}

@inproceedings{,
title={Pattern Analysis of Dynamic Susceptibility Contrast-enhanced MR Imaging Demonstrates Peritumoral Tissue Heterogeneity},
author={Hamed Akbari, MD, PhD Luke Macyszyn, MD, MA Xiao Da, MS Ronald L. Wolf, MD, PhD Michel Bilello, MD, PhD Ragini Verma, PhD Donald M. O�Rourke, MD Christos Davatzikos, PhD},
pages={502-510},
year={2014},
organization={RSNA Radiology}
}
*/
class PseudoProgressionEstimator
#ifdef APP_BASE_CAPTK_H
  public ApplicationBase
#endif
{
public:
  /**
  \brief Default Constructor
  */
  PseudoProgressionEstimator()
  {
    mPseudoTrainedFile = "PSU_SVM_Model.xml";
    mRecurrenceTrainedFile = "REC_SVM_Model.xml";
    logger.UseNewFile(loggerFile);
  };

  //! Default Destructor
  ~PseudoProgressionEstimator();

  //! Gets last error
  std::string mLastEncounteredError;
  std::string mPseudoTrainedFile;
  std::string mRecurrenceTrainedFile;
  cbica::Logging logger;


  VariableSizeMatrixType GetModelSelectedFeatures(VariableSizeMatrixType & ScaledFeatureSetAfterAddingLabel, VariableLengthVectorType & psuSelectedFeatures);
  VectorVectorDouble CombineAllThePerfusionFeaures(VectorVectorDouble T1IntensityHistogram,
    VectorVectorDouble TCIntensityHistogram,
    VectorVectorDouble T1TCIntensityHistogram,
    VectorVectorDouble T2IntensityHistogram,
    VectorVectorDouble FLIntensityHistogram,
    VectorVectorDouble T2FLIntensityHistogram,
    VectorVectorDouble AXIntensityHistogram,
    VectorVectorDouble FAIntensityHistogram,
    VectorVectorDouble RDIntensityHistogram,
    VectorVectorDouble TRIntensityHistogram,
    VectorVectorDouble PHIntensityHistogram,
    VectorVectorDouble PSIntensityHistogram,
    VectorVectorDouble RCIntensityHistogram,
    VectorVectorDouble PCA1IntensityHistogram,
    VectorVectorDouble PCA2IntensityHistogram,
    VectorVectorDouble PCA3IntensityHistogram,
    VectorVectorDouble PCA4IntensityHistogram,
    VectorVectorDouble PCA5IntensityHistogram,
    VectorVectorDouble PCA6IntensityHistogram,
    VectorVectorDouble PCA7IntensityHistogram,
    VectorVectorDouble PCA8IntensityHistogram,
    VectorVectorDouble PCA9IntensityHistogram,
    VectorVectorDouble PCA10IntensityHistogram,
    std::string output);

  vtkSmartPointer<vtkTable> MakePCAMatrix(int NumberOfFeatures, int NumberOfSamples);

  VariableSizeMatrixType ColumnWiseScaling(VariableSizeMatrixType PerfusionData);

  PerfusionMapType CombineAndCalculatePerfusionPCA(PerfusionMapType PerfusionDataMap, VariableSizeMatrixType &TransformationMatrix, VariableLengthVectorType &MeanVector);
  PerfusionMapType CombineAndCalculatePerfusionPCAForTestData(PerfusionMapType PerfusionDataMap, VariableSizeMatrixType &TransformationMatrix, VariableLengthVectorType &MeanVector);

  template<class PerfusionImageType, class ImageType>
  VariableSizeMatrixType LoadPerfusionData(typename ImageType::Pointer maskImagePointerNifti, typename PerfusionImageType::Pointer perfImagePointerNifti, std::vector< typename ImageType::IndexType> &indices);

  /**
  Image Intensity using the itk::RescaleIntensityImageFilter to the 0-255 range
  */
  ImageTypeFloat3D::Pointer RescaleImageIntensity(ImageTypeFloat3D::Pointer image);

  template<class ImageType>
  std::tuple<VectorDouble, VectorDouble, VectorDouble, VectorDouble, VectorDouble> GetAllFeaturesPerImagePerROI(typename ImageType::Pointer image, typename ImageType::Pointer mask, std::string modality);

  PerfusionMapType CombinePerfusionDataAndApplyExistingPerfusionModel(PerfusionMapType PerfusionDataMap, VariableSizeMatrixType TransformationMatrix, VariableLengthVectorType MeanVector);

  void ReadAllTheModelParameters(std::string modeldirectory,
    VariableSizeMatrixType &PCA_PERF,
    VariableSizeMatrixType &PCA_T1,
    VariableSizeMatrixType &PCA_T1CE,
    VariableSizeMatrixType &PCA_T2,
    VariableSizeMatrixType &PCA_FL,
    VariableSizeMatrixType &PCA_T1T1CE,
    VariableSizeMatrixType &PCA_T2FL,
    VariableSizeMatrixType &PCA_AX,
    VariableSizeMatrixType &PCA_FA,
    VariableSizeMatrixType &PCA_RAD,
    VariableSizeMatrixType &PCA_TR,
    VariableSizeMatrixType &PCA_PH,
    VariableSizeMatrixType &PCA_PSR,
    VariableSizeMatrixType &PCA_RCBV,
    VariableSizeMatrixType &PCA_PC1,
    VariableSizeMatrixType &PCA_PC2,
    VariableSizeMatrixType &PCA_PC3,
    VariableSizeMatrixType &PCA_PC4,
    VariableSizeMatrixType &PCA_PC5,
    VariableSizeMatrixType &PCA_PC6,
    VariableSizeMatrixType &PCA_PC7,
    VariableSizeMatrixType &PCA_PC8,
    VariableSizeMatrixType &PCA_PC9,
    VariableSizeMatrixType &PCA_PC10,
    VariableLengthVectorType &Mean_PERF,
    VariableLengthVectorType &Mean_T1,
    VariableLengthVectorType &Mean_T1CE,
    VariableLengthVectorType &Mean_T2,
    VariableLengthVectorType &Mean_FL,
    VariableLengthVectorType &Mean_T1T1CE,
    VariableLengthVectorType &Mean_T2FL,
    VariableLengthVectorType &Mean_AX,
    VariableLengthVectorType &Mean_FA,
    VariableLengthVectorType &Mean_RAD,
    VariableLengthVectorType &Mean_TR,
    VariableLengthVectorType &Mean_PH,
    VariableLengthVectorType &Mean_PSR,
    VariableLengthVectorType &Mean_RCBV,
    VariableLengthVectorType &Mean_PC1,
    VariableLengthVectorType &Mean_PC2,
    VariableLengthVectorType &Mean_PC3,
    VariableLengthVectorType &Mean_PC4,
    VariableLengthVectorType &Mean_PC5,
    VariableLengthVectorType &Mean_PC6,
    VariableLengthVectorType &Mean_PC7,
    VariableLengthVectorType &Mean_PC8,
    VariableLengthVectorType &Mean_PC9,
    VariableLengthVectorType &Mean_PC10);

  /**
  \brief Train a new model based on the inputs
  */
  bool TrainNewModelOnGivenData(const std::vector<std::map<CAPTK::ImageModalityType, std::string>> qualifiedsubjects, const std::string &outputdirectory, bool useConventioalData, bool useDTIData, bool usePerfData, bool useDistData);

  VariableSizeMatrixType LoadPseudoProgressionTrainingData(const std::vector<std::map<CAPTK::ImageModalityType, std::string>> &trainingsubjects, std::vector<double> &traininglabels, std::string outputdirectory);

  VariableSizeMatrixType LoadPseudoProgressionTestingData(const std::vector<std::map<CAPTK::ImageModalityType, std::string>> &trainingsubjects, std::vector<double> &traininglabels, std::string outputdirectory, std::string modeldirectory);
  /**
  \brief Recurrence Estimation using existing model
  */
  bool PseudoProgressionEstimateOnExistingModel(std::vector<std::map<CAPTK::ImageModalityType, std::string>> qualifiedsubjects,
    const std::string &modeldirectory,
    const std::string &inputdirectory,
    const std::string &outputdirectory,
    bool useConventionalData, bool useDTIData, bool usePerfData, bool useDistData);

  VariableSizeMatrixType SelectModelFeatures(const VariableSizeMatrixType &ModelFeatures, const VectorDouble &selectedFeatures);

  template<class ImageType>
  typename ImageType::Pointer RescaleImageIntensity(const typename ImageType::Pointer &image);

  VectorDouble GetHistogramFeatures(std::vector<float> intensities, int m_Bins);
  VectorDouble GetHistogramFeaturesWhole(std::vector<float> intensities);

  template< class TImageTypeShape = TImageType >
  VectorDouble GetShapeFeatures(typename TImageTypeShape::Pointer mask);

  VectorDouble GetIntensityFeatures(std::vector<float> m_intensities);

  template<class ImageType>
  VectorDouble GetGLCMFeatures(typename ImageType::Pointer image, typename ImageType::Pointer mask);

  template<class ImageType>
  typename ImageType::Pointer MakeAdditionalModality(typename ImageType::Pointer image1, typename ImageType::Pointer image2);


  template<class ImageType>
  VectorDouble GetRunLengthFeatures(typename ImageType::Pointer image, typename ImageType::Pointer mask);

  template<class PerfusionImageType, class ImageType>
  VectorVectorDouble GetPerfusionFeatures(typename PerfusionImageType::Pointer image, typename ImageType::Pointer mask);

  NiftiDataManager mNiftiLocalPtr;
  OutputWritingManager mOutputLocalPtr;
  FeatureReductionClass mFeatureReductionLocalPtr;
  FeatureScalingClass mFeatureScalingLocalPtr;
  FeatureExtractionClass mFeatureExtractionLocalPtr;
  //SVMClassificationClass mClassificationLocalPtr;

  void Run()
  {
  }

private:
  std::string mCurrentOutputDir;
};

template<class ImageType>
typename ImageType::Pointer PseudoProgressionEstimator::RescaleImageIntensity(const typename ImageType::Pointer &image)
{
  typedef itk::RescaleIntensityImageFilter< ImageType, ImageType > RescaleFilterType;
  typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput(image);
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(255);
  rescaleFilter->Update();
  typename ImageType::Pointer outputimage = rescaleFilter->GetOutput();
  return outputimage;
}

template<class PerfusionImageType, class ImageType>
VariableSizeMatrixType PseudoProgressionEstimator::LoadPerfusionData(typename ImageType::Pointer maskImagePointerNifti, typename PerfusionImageType::Pointer perfImagePointerNifti, std::vector<typename ImageType::IndexType> &qualifiedIndices)
{
  VectorVectorDouble pIntensities;
  //----------------------find the voxels of the mask------------------------------
  typedef itk::ImageRegionIteratorWithIndex <ImageTypeFloat3D> IteratorType;
  IteratorType maskIt(maskImagePointerNifti, maskImagePointerNifti->GetLargestPossibleRegion());

  maskIt.GoToBegin();
  int mask_counter = 0;
  while (!maskIt.IsAtEnd())
  {
    if (maskIt.Get() == 1)
    {
      mask_counter++;
      qualifiedIndices.push_back(maskIt.GetIndex());
    }
    ++maskIt;
  }
  //--------------------------populate the covariance matrix --------------------------------
  typename ImageTypeFloat4D::RegionType region = perfImagePointerNifti->GetLargestPossibleRegion();
  VariableSizeMatrixType revisedcovariancematrix;
  revisedcovariancematrix.SetSize(mask_counter, region.GetSize()[3]);
  for (size_t i = 0; i < qualifiedIndices.size(); i++)
  {
    typename ImageTypeFloat4D::IndexType regionIndex;
    regionIndex[0] = qualifiedIndices[i][0];
    regionIndex[1] = qualifiedIndices[i][1];
    regionIndex[2] = qualifiedIndices[i][2];
    for (size_t j = 0; j < region.GetSize()[3]; j++)
    {
      regionIndex[3] = j;
      revisedcovariancematrix(i, j) = perfImagePointerNifti->GetPixel(regionIndex);
    }

  }
  return revisedcovariancematrix;

}

template<class ImageType>
std::tuple<VectorDouble, VectorDouble, VectorDouble, VectorDouble, VectorDouble> PseudoProgressionEstimator::GetAllFeaturesPerImagePerROI(typename ImageType::Pointer image, typename ImageType::Pointer mask, std::string modality)
{
  std::vector<typename ImageType::IndexType> roiIndices;
  typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
  IteratorType imIt(mask, mask->GetLargestPossibleRegion());
  imIt.GoToBegin();
  while (!imIt.IsAtEnd())
  {
    if (imIt.Get() == 1)
      roiIndices.push_back(imIt.GetIndex());
    ++imIt;
  }


  typedef std::vector<float> VectorFloat;
  VectorFloat ROIIntensities;
  for (unsigned int i = 0; i < roiIndices.size(); i++)
    ROIIntensities.push_back(std::round(image.GetPointer()->GetPixel(roiIndices[i])));

  VectorDouble HistogramFeatures = GetHistogramFeatures(ROIIntensities, 10);
  VectorDouble HistogramFeatures1 = GetHistogramFeaturesWhole(ROIIntensities);
  VectorDouble IntensityFeatures = GetIntensityFeatures(ROIIntensities);
  VectorDouble GLCMFeatures = GetGLCMFeatures<ImageType>(image, mask);
  VectorDouble GLRLMFeatures = GetRunLengthFeatures<ImageType>(image, mask);

  //VectorDouble IntensityFeatures;
  //IntensityFeatures.push_back(100);
  //IntensityFeatures.push_back(100);
  //IntensityFeatures.push_back(100);
  //IntensityFeatures.push_back(100);
  //IntensityFeatures.push_back(100);
  //IntensityFeatures.push_back(100);
  //IntensityFeatures.push_back(100);

  /* VectorDouble GLCMFeatures;
  GLCMFeatures.push_back(100);
  GLCMFeatures.push_back(100);
  GLCMFeatures.push_back(100);
  GLCMFeatures.push_back(100);
  GLCMFeatures.push_back(100);
  GLCMFeatures.push_back(100);
  GLCMFeatures.push_back(100);
  GLCMFeatures.push_back(100);

  VectorDouble GLRLMFeatures;
  GLRLMFeatures.push_back(100);
  GLRLMFeatures.push_back(100);
  GLRLMFeatures.push_back(100);
  GLRLMFeatures.push_back(100);
  GLRLMFeatures.push_back(100);
  GLRLMFeatures.push_back(100);
  GLRLMFeatures.push_back(100);
  GLRLMFeatures.push_back(100);
  GLRLMFeatures.push_back(100);
  GLRLMFeatures.push_back(100);*/

  std::tuple<VectorDouble, VectorDouble, VectorDouble, VectorDouble, VectorDouble> new_tuple(HistogramFeatures, IntensityFeatures, GLCMFeatures, GLRLMFeatures, HistogramFeatures1);
  return new_tuple;
}

template< typename TImageTypeShape >
VectorDouble PseudoProgressionEstimator::GetShapeFeatures(typename TImageTypeShape::Pointer mask)
{
  std::vector<double> features;

  typedef  short LabelType;
  typedef itk::Image< LabelType, TImageTypeShape::ImageDimension > OutputImageType;
  typedef itk::ShapeLabelObject< LabelType, TImageTypeShape::ImageDimension > ShapeLabelObjectType;
  typedef itk::LabelMap< ShapeLabelObjectType > LabelMapType;

  typedef itk::ConnectedComponentImageFilter < TImageTypeShape, OutputImageType > ConnectedComponentImageFilterType;
  typedef itk::LabelImageToShapeLabelMapFilter < OutputImageType, LabelMapType > I2LType;
  typename ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
  connected->SetInput(mask);
  connected->Update();

  typename I2LType::Pointer i2l = I2LType::New();
  i2l->SetInput(connected->GetOutput());
  i2l->SetComputePerimeter(true);
  i2l->Update();
  typename LabelMapType::Pointer labelMap = i2l->GetOutput();
  auto spacing = mask->GetSpacing();
  double voxvol = spacing[0] * spacing[1] * spacing[2];

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

  features.push_back(mean_ecc1 / eccentricity1.size());
  features.push_back(mean_elong / elongation.size());
  features.push_back(mean_perim);
  features.push_back(mean_round / roundness.size());
  features.push_back(mean_flat / flatness.size());

  return features;
}

template<class ImageType>
VectorDouble PseudoProgressionEstimator::GetGLCMFeatures(typename ImageType::Pointer image, typename ImageType::Pointer mask)
{
  double m_Bins = 16;
  using FeatureextractionImageType = typename ImageType::Pointer;
  using Image2CoOccuranceType = itk::Statistics::ScalarImageToCooccurrenceMatrixFilter<ImageType>;
  using HistogramType = typename Image2CoOccuranceType::HistogramType;
  using Hist2FeaturesType = itk::Statistics::HistogramToTextureFeaturesFilter<HistogramType>;
  using OffsetType = typename ImageType::OffsetType;
  using Offsets = OffsetType;
  using OffsetVector = itk::VectorContainer< unsigned char, OffsetType >;


  double inputRadius = 1;
  double inputDirections = 13;
  itk::Neighborhood< typename ImageType::PixelType, ImageType::ImageDimension > neighborhood;
  neighborhood.SetRadius(inputRadius);
  auto size = neighborhood.GetSize();
  auto directionsToCompute = 1;
  for (size_t sizeCount = 0; sizeCount < ImageType::ImageDimension; sizeCount++)
  {
    directionsToCompute *= size[sizeCount];
  }
  if (inputDirections < directionsToCompute)
  {
    directionsToCompute = inputDirections;
  }

  typename OffsetVector::Pointer offsets = OffsetVector::New();

  for (int d = 0; d < directionsToCompute; d++)
  {
    offsets->push_back(neighborhood.GetOffset(d));
  }


  std::vector<double> features;
  double contrast = 0, correl = 0, ener = 0, entro = 0, homo = 0, clustershade = 0, clusterprominance = 0, autocorr = 0;

  for (size_t i = 0; i < offsets->size(); i++)
  {
    auto glcmGenerator = Image2CoOccuranceType::New();
    glcmGenerator->SetNumberOfBinsPerAxis(m_Bins); //reasonable number of bins
    glcmGenerator->SetPixelValueMinMax(0, 255);
    glcmGenerator->SetMaskImage(mask);
    glcmGenerator->SetInput(image);
    auto featureCalc = Hist2FeaturesType::New();

    glcmGenerator->SetOffset(offsets->at(i));
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

  contrast = contrast / offsets->size();
  correl = correl / offsets->size();
  ener = ener / offsets->size();
  homo = homo / offsets->size();
  entro = entro / offsets->size();
  clusterprominance = clusterprominance / offsets->size();
  clustershade = clustershade / offsets->size();
  autocorr = autocorr / offsets->size();

  features.push_back(correl);
  features.push_back(contrast);
  features.push_back(entro);
  features.push_back(homo);
  features.push_back(clustershade);
  features.push_back(clusterprominance);
  features.push_back(autocorr);
  features.push_back(ener);

  // TODO: Sung to add his GLCM extraction code here
  //featurevec[std::string(IndividualFeaturesString[Correlation]) + "_Sung"] = 0;
  return features;
}

template<class ImageType>
VectorDouble PseudoProgressionEstimator::GetRunLengthFeatures(typename ImageType::Pointer image, typename ImageType::Pointer mask)
{
  using FeatureextractionImageType = typename ImageType::Pointer;
  using Image2CoOccuranceType = itk::Statistics::ScalarImageToCooccurrenceMatrixFilter<ImageType>;
  using HistogramType = typename Image2CoOccuranceType::HistogramType;
  using Hist2FeaturesType = itk::Statistics::HistogramToTextureFeaturesFilter<HistogramType>;
  using OffsetType = typename ImageType::OffsetType;
  using Offsets = OffsetType;
  using OffsetVector = itk::VectorContainer< unsigned char, OffsetType >;

  double m_Bins = 16;
  double inputRadius = 1;
  double inputDirections = 13;
  itk::Neighborhood< typename ImageType::PixelType, ImageType::ImageDimension > neighborhood;
  neighborhood.SetRadius(inputRadius);
  auto size = neighborhood.GetSize();
  auto directionsToCompute = 1;
  for (size_t sizeCount = 0; sizeCount < ImageType::ImageDimension; sizeCount++)
  {
    directionsToCompute *= size[sizeCount];
  }
  if (inputDirections < directionsToCompute)
  {
    directionsToCompute = inputDirections;
  }

  typename OffsetVector::Pointer offsets = OffsetVector::New();

  for (int d = 0; d < directionsToCompute; d++)
  {
    offsets->push_back(neighborhood.GetOffset(d));
  }
  //---------------------------------------------
  using HistogramFrequencyContainerType = itk::Statistics::DenseFrequencyContainer2;
  using RunLengthFilterType = itk::Statistics::ScalarImageToRunLengthFeaturesFilter< ImageType, HistogramFrequencyContainerType >;
  using RunLengthMatrixGenerator = typename RunLengthFilterType::RunLengthMatrixFilterType;
  using RunLengthFeatures = typename RunLengthFilterType::RunLengthFeaturesFilterType;
  using InternalRunLengthFeatureName = typename RunLengthFeatures::RunLengthFeatureName;

  typename  RunLengthMatrixGenerator::Pointer matrix_generator = RunLengthMatrixGenerator::New();
  matrix_generator->SetMaskImage(mask);
  matrix_generator->SetInput(image);
  matrix_generator->SetInsidePixelValue(1);
  matrix_generator->SetPixelValueMinMax(0, 255);
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

  std::vector<double> features;
  typename  OffsetVector::ConstIterator offsetIt;
  size_t offsetNum = 0, featureNum = 0;

  std::vector<double> tempfeatures/*(requestedFeatures->size(), 0)*/;
  tempfeatures.resize(requestedFeatures->size());
  for (offsetIt = offsets->Begin(); offsetIt != offsets->End(); offsetIt++, offsetNum++)
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
      tempfeatures[featureNum] = tempfeatures[featureNum] + runLengthMatrixCalculator->GetFeature((InternalRunLengthFeatureName)fnameIt.Value());
    }
  }
  for (size_t i = 0; i < tempfeatures.size(); i++)
  {
    tempfeatures[i] = tempfeatures[i] / offsets->size();
    features.push_back(tempfeatures[i]);
  }
  return features;
}

template<class ImageType>
typename ImageType::Pointer PseudoProgressionEstimator::MakeAdditionalModality(typename ImageType::Pointer image1, typename ImageType::Pointer image2)
{
  typename ImageType::Pointer SubtractedImage = ImageType::New();
  SubtractedImage->CopyInformation(image1);
  SubtractedImage->SetRequestedRegion(image1->GetLargestPossibleRegion());
  SubtractedImage->SetBufferedRegion(image1->GetBufferedRegion());
  SubtractedImage->Allocate();
  SubtractedImage->FillBuffer(0);

  typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
  IteratorType subIt(SubtractedImage, SubtractedImage->GetLargestPossibleRegion());
  IteratorType im1It(image1, image1->GetLargestPossibleRegion());
  IteratorType im2It(image2, image2->GetLargestPossibleRegion());


  im1It.GoToBegin();
  im2It.GoToBegin();
  subIt.GoToBegin();
  while (!im1It.IsAtEnd())
  {
    subIt.Set(im1It.Get() - im2It.Get());
    ++im1It;
    ++im2It;
    ++subIt;
  }
  return SubtractedImage;
}
#endif

