/**
\file  EGFRvIIIIndexPredictor.h

\brief The header file containing the SurvivaPredictor class, used to find Survival Prediction Index
Library Dependencies: ITK 4.7+ <br>

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software/license.html

*/


#ifndef _EGFRvIIIIndexPredictor_h_
#define _EGFRvIIIIndexPredictor_h_

#include "FeatureScalingClass.h"
#include "FeatureExtractionClass.h"
#include "cbicaLogging.h"
#include "CaPTkEnums.h"
#include "itkCSVArray2DFileReader.h"
#include "itkImageFileReader.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"

#ifdef APP_BASE_CAPTK_H
#include "ApplicationBase.h"
#endif

typedef itk::Image< float, 3 > ImageType;
typedef itk::CSVArray2DFileReader<double> CSVFileReaderType;
typedef vnl_matrix<double> MatrixType;

#define EGFR_SIZE_COMP 100
#define EGFR_NO_NEAR_DIST 100
#define EGFR_NO_OF_FEATURES 448

/**
\class EGFRvIIIIndexPredictor

\brief Calculates EGFRvIII Status

Reference:

@article{akbari2018invivo,
title={In vivo evaluation of EGFRvIII mutation in primary glioblastoma patients via complex multiparametric MRI signature},
author={H.Akbari, S.Bakas, J.M.Pisapia, M.P.Nasrallah, M.Rozycki, M.Martinez-Lage, J.D.Morrissette, N.Dahmane, D.M.O’Rourke, C.Davatzikos},
journal={Neuro-oncology},
volume={20},
number={8},
pages={1068-1079},
year={2018},
publisher={Society for Neuro-Oncology}
}

*/
class EGFRvIIIIndexPredictor
#ifdef APP_BASE_CAPTK_H
  : public ApplicationBase
#endif
{
public:
  //! Default constructor
  EGFRvIIIIndexPredictor()
  {
    mTrainedFile = "EGFRvIII_SVM_Model.xml";
    logger.UseNewFile(loggerFile);
  }

  //! Default destructor
  ~EGFRvIIIIndexPredictor() {};

  std::string mTrainedFile;
  FeatureExtractionClass mFeatureExtractionLocalPtr;
  FeatureScalingClass mFeatureScalingLocalPtr;
  cbica::Logging logger;
  /**
  \brief Calculates the features for given images of one subject

  \param t1ceImagePointer			          Pointer to T1CE image
  \param t2flairImagePointer		        Pointer to T2-FLAIR image
  \param t1ImagePointer			            Pointer to T1-weighted image
  \param t2ImagePointer			            Pointer to T2 image
  \param rcbvImagePointer			          Pointer to RCBV image
  \param psrImagePointer			          Pointer to Percent Signal Recovery (PSR) image
  \param phImagePointer			            Pointer to Peak Height (PH) image
  \param axImagePointer			            Pointer to Axial Diffusivity image
  \param faImagePointer			            Pointer to Fractional Anisotropy image
  \param radImagePointer			          Pointer to Radial Diffusivity image
  \param trImagePointer			            Pointer to Trace image
  \param labelImagePointer		          Pointer to image having segmentation labels in patient space
  \param atlasImagePointer		          Pointer to image having segmentation labels in atla space
  \param templateImagePointer9regions		Pointer to atlas (having 9 regions) to calculate location features
  \param templateImagePointer21regions	Pointer to atlas (having 21 regions) to calculate location features
  \param templateImagePointerAllregions	Pointer to atlas (having all regions) to calculate location features
  \param POSAtlasImagePointer       		Pointer to EGFR positive atlas
  \param NEGAtlasImagePointer       		Pointer to EGFR negative atlas
  \param configuration			            Configuration to calculate histogram features
  */

  template<class ImageType>
  VectorDouble  LoadTestData(const typename ImageType::Pointer &t1ceImagePointer,
    const typename ImageType::Pointer &t2flairImagePointer,
    const typename ImageType::Pointer &t1ImagePointer,
    const typename ImageType::Pointer &t2ImagePointer,
    const typename ImageType::Pointer &rcbvImagePointer,
    const typename ImageType::Pointer &psrImagePointer,
    const typename ImageType::Pointer &phImagePointer,
    const typename ImageType::Pointer &axImagePointer,
    const typename ImageType::Pointer &faImagePointer,
    const typename ImageType::Pointer &radImagePointer,
    const typename ImageType::Pointer &trImagePointer,
    const typename ImageType::Pointer &labelImagePointer,
    const typename ImageType::Pointer &atlasImagePointer,
    const typename ImageType::Pointer &templateImagePointer9regions,
    const typename ImageType::Pointer &templateImagePointer21regions,
    const typename ImageType::Pointer &templateImagePointerAllregions,
    const typename ImageType::Pointer &POSAtlasImagePointer,
    const typename ImageType::Pointer &NEGAtlasImagePointer,
    const VariableSizeMatrixType &configuration);

  template<class ImageType>
  std::vector<double> GetSpatialLocationFeaturesForEGFR(typename ImageType::Pointer labelImagePointer,
    typename ImageType::Pointer EGFRvIII_Neg_Atlas,
    typename ImageType::Pointer EGFRvIII_Pos_Atlas);

  /**
  \brief Calculates the statistcal features (mean, satndard deviation) for given vector

  \param intensities			Input vector having intensities from a region of an image
  */
  VectorDouble GetStatisticalFeatures(const VectorDouble &intensities);

  /**
  \brief Normalizes the input histogram features

  \param intensities			Input vector having intensities from a region of an image
  \param size             Size of the region. Number of voxels in the histogram
  */
  VectorDouble NormalizeHistogramFeatures(VectorDouble inputHistogramFeatures, const int size);

  /**
  \brief Calculates the histogram binning based features for a given intensity vector

  \param intensities			Input vector having intensities from a region of an image
  \param start				    Start of the histogram bins
  \param interval				  Interval between the centers of two hsitogram bins
  \param end					    End of the histogram bins
  */
  VectorDouble GetHistogramFeatures(const VectorDouble &intensities, const double &start, const double &interval, const double &end);

  /**
  \brief Calculates volumetric measures and their ratios

  \param edemaSize			Volume of edema in terms of number of voxels
  \param tuSize				  Volume of enhancing tumor in terms of number of voxels
  \param neSize				  Volume of non-enahancing tumro core in terms of number of voxels
  \param totalSize			Volume of whole brain in terms of number of voxels
  */
  VectorDouble GetVolumetricFeatures(const double &edemaSize, const double &tuSize, const double &neSize, const double &totalSize);

  /**
  \brief Removes the connected components smaller than an area threshold from the tumor

  \param labelImagePointer		Pointer to image having segmentation labels in patient space
  \parm areathreshold         Area threshold to remove smaller regions
  */
  template<class ImageType>
  std::vector<typename ImageType::Pointer> RevisedTumorArea(const typename ImageType::Pointer &labelImagePointer, const int areathreshold);

  /**
  \brief Calculates distance features. 1. Distance of tumor from ventricles, 2. Distance of edema from ventricles

  \param edemaImage		Pointer to image having edema label
  \param tumorImage		Pointer to image having (enhancing tumor + non-enhancing tumor) label
  \param ventImage		Pointer to image having ventricles label
  */
  template<class ImageType>
  VectorDouble GetDistanceFeatures3(const typename ImageType::Pointer &edemaImage, 
    const typename ImageType::Pointer &tumorImage, 
    const typename ImageType::Pointer &ventImage);

  /**
  \brief Calculates the distance map of an image by using the given label as reference seed point

  \param labelImagePointer		Pointer to image having segmentation labels in patient space
  \param label					      Label to use as a reference seed point
  */
  template<class ImageType>
  typename ImageType::Pointer GetDistanceMapWithLabel(const typename ImageType::Pointer &labelImagePointer, 
    const int &label);

  /**
  \brief Calculates the distance map of an image by using non-zero voxels of the image as reference seed point

  \param image					Pointer to image having reference seed points
  */
  template<class ImageType>
  typename ImageType::Pointer GetDistanceMap(const typename ImageType::Pointer &image);

  /**
  \brief Calculates spatial location features for given tumor images

  \param labelImagePointer		Pointer to image having segmentation labels in atlas space
  \param atlasImagePointer		Pointer to jakob atlas image
  */
  template<class ImageType>
  VectorDouble GetSpatialLocationFeaturesForFixedNumberOfRegions(const typename ImageType::Pointer &labelImagePointer,const typename ImageType::Pointer &atlasImagePointer,const int number);
  
  /**
  \brief Apply an exsiting model on new patients to predit survival

  \param modeldirectory			Path to the directory where model files are stored
  \param inputdirectory			Path to the directory having test data
  \param qualifiedsubjects	List of qualifeid subjects having all the data available to apply a model
  \param outputdirectory		Path to the output directory
  */
  VectorDouble EGFRvIIIPredictionOnExistingModel(const std::string &modeldirectory,
    const std::string &inputdirectory,
    const std::vector< std::map< CAPTK::ImageModalityType, std::string > > &qualifiedsubjects,
    const std::string &outputdirectory);

  /**
  \brief Extracts the selected features from input data 

  \param inputFeatures			Input features matrix
  \param SelectedFeatures		Indices of features to be extracted
  */
  VariableSizeMatrixType SelectModelFeatures(const VariableSizeMatrixType &InputFeatures,
    const VariableLengthVectorType &SelectedFeatures);

  /**
  \brief Calculates distance from the hyperplane of linear SVM classifier

  \param testData			Final features of the test patients
  \param filename		  Path to the model file
  */
  VariableLengthVectorType DistanceFunctionLinear(const VariableSizeMatrixType &testData, const std::string &filename);

  /**
  \brief Calculates transpose of given input matrix

  \param  inputmatrix Input data
  */
  VariableSizeMatrixType MatrixTranspose(const VariableSizeMatrixType &inputmatrix);

  /**
  \brief Rescales intensity of an image in the range of 0-255

  \param  image     Input image
  */
  template<class ImageType>
  typename ImageType::Pointer RescaleImageIntensity(const typename ImageType::Pointer &image);

  /**
  \brief Caps the intensity profile of an image 

  \param  image     Input image
  */
  template<class ImageType>
  typename ImageType::Pointer CapImageIntensityWithPercentile(const typename ImageType::Pointer &image);


  //functions specific to the training process
  ImageType::Pointer MakeAtlases(const std::vector<ImageType::Pointer> &SegmentationaAllPatients, const std::vector<int> &LabelsAllPatientsconst, const int atlas_no);
  
  bool SelectModelFeaturesForTrainingSVMFFSel(const VariableSizeMatrixType &scaledFeatureSet, 
    const VectorDouble &AllLabels,const std::string outputdirectory);
  
  std::vector<int> UpdateUnselectedFeatures(std::vector<int> SelectedFeatures, int featuresize);
  
  VectorDouble InternalCrossValidationResubstitution(VariableSizeMatrixType inputFeatures, 
    std::vector<double> inputLabels, 
    double cValue, double gValue, int kerneltype);
  
  VectorDouble CalculatePerformanceMeasures(VariableLengthVectorType predictedLabels, VectorDouble GivenLabels);

  bool SelectModelFeaturesForTrainingFromStudy(const VariableSizeMatrixType &scaledFeatureSet, 
    const VectorDouble &AllLabels,
    const std::string outputdirectory);
  template<class ImageType>

  std::vector<double> GetSpatialLocationFeaturesForEGFRTrainingWithLOO(typename ImageType::Pointer labelImagePointer,
    typename ImageType::Pointer EGFRvIII_Neg_Atlas,
    typename ImageType::Pointer EGFRvIII_Pos_Atlas, int label);

  int PrepareNewEGFRvIIIPredictionModel(
    const std::string &inputdirectory,
    const std::vector< std::map< CAPTK::ImageModalityType, std::string > > &qualifiedsubjects,
    const std::string &outputdirectory);

  void Run()
  {}

private:

};

template<class ImageType>
typename ImageType::Pointer EGFRvIIIIndexPredictor::RescaleImageIntensity(const typename ImageType::Pointer &image)
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

template<class ImageType>
typename ImageType::Pointer EGFRvIIIIndexPredictor::CapImageIntensityWithPercentile(const typename ImageType::Pointer &image)
{
  VectorDouble imageIntensities; 
  typename ImageType::Pointer interImage = ImageType::New();
  interImage->CopyInformation(image);
  interImage->SetRequestedRegion(image->GetLargestPossibleRegion());
  interImage->SetBufferedRegion(image->GetBufferedRegion());
  interImage->Allocate();
  interImage->FillBuffer(0);


  typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
  IteratorType imIt(image, image->GetLargestPossibleRegion());
  IteratorType iminterIt(interImage, interImage->GetLargestPossibleRegion());

  imIt.GoToBegin();
  while (!imIt.IsAtEnd())
  {
    imageIntensities.push_back(imIt.Get());
    ++imIt;
  }
  std::sort(imageIntensities.begin(), imageIntensities.end(), std::less<double>());
  float capValue = imageIntensities[std::floor(99.99*(imageIntensities.size()) / 100)];

  imIt.GoToBegin();
  iminterIt.GoToBegin();
  while (!imIt.IsAtEnd())
  {
    if (imIt.Get() > capValue)
      iminterIt.Set(capValue);
    else
      iminterIt.Set(imIt.Get());

    ++imIt;
    ++iminterIt;
  }
  return interImage;
}


template<class ImageType>
VectorDouble EGFRvIIIIndexPredictor::GetDistanceFeatures3(const typename ImageType::Pointer &edemaImage, const typename ImageType::Pointer &tumorImage, const typename ImageType::Pointer &ventImage)
{
  //create distance image
  VectorDouble VentEdemaSum;
  VectorDouble VentTumorSum;
  
  //get sum of images
  typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
  IteratorType tumorIt(tumorImage, tumorImage->GetLargestPossibleRegion());
  IteratorType ventIt(ventImage, ventImage->GetLargestPossibleRegion());
  IteratorType edemaIt(edemaImage, edemaImage->GetLargestPossibleRegion());

  //convert negative values to zeros
  tumorIt.GoToBegin();
  ventIt.GoToBegin();
  edemaIt.GoToBegin();

  while (!tumorIt.IsAtEnd())
  {
    if (edemaIt.Get() < 0)
      edemaIt.Set(0);
    if (ventIt.Get() < 0)
      ventIt.Set(0);
    if (tumorIt.Get() < 0)
      tumorIt.Set(0);

    ++tumorIt;
    ++ventIt;
    ++edemaIt;
  }
  tumorIt.GoToBegin();
  ventIt.GoToBegin();
  edemaIt.GoToBegin();

  while (!tumorIt.IsAtEnd())
  {
    VentEdemaSum.push_back(edemaIt.Get() + ventIt.Get());
    VentTumorSum.push_back(tumorIt.Get() + ventIt.Get());
    ++tumorIt;
    ++ventIt;
    ++edemaIt;
  }
  std::sort(VentEdemaSum.begin(), VentEdemaSum.end());
  std::sort(VentTumorSum.begin(), VentTumorSum.end());

  double sumVentEdema = 0;
  double sumVentTumor = 0;
  for (int i = 0; i < EGFR_NO_NEAR_DIST; i++)
  {
    sumVentEdema = sumVentEdema + VentEdemaSum[i];
    sumVentTumor = sumVentTumor + VentTumorSum[i];
  }
  VectorDouble DistanceFeatures;
  DistanceFeatures.push_back(sumVentTumor / EGFR_NO_NEAR_DIST); 
  DistanceFeatures.push_back(sumVentEdema / EGFR_NO_NEAR_DIST);
  return DistanceFeatures;
}

template<class ImageType>
typename ImageType::Pointer EGFRvIIIIndexPredictor::GetDistanceMap(const typename ImageType::Pointer &labelImagePointer)
{
  typedef itk::Image<unsigned char, 3>  UnsignedCharImageType;
  typedef itk::CastImageFilter< ImageType, UnsignedCharImageType > CastFilterType;
  typename CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(labelImagePointer);
  typedef  itk::SignedMaurerDistanceMapImageFilter< UnsignedCharImageType, ImageType  > SignedMaurerDistanceMapImageFilterType;
  typename  SignedMaurerDistanceMapImageFilterType::Pointer distanceMapImageFilter = SignedMaurerDistanceMapImageFilterType::New();
  distanceMapImageFilter->SetInput(castFilter->GetOutput());
  distanceMapImageFilter->Update();
  typename ImageType::Pointer DistanceMap = distanceMapImageFilter->GetOutput();
  return DistanceMap;
}

template<class ImageType>
typename ImageType::Pointer EGFRvIIIIndexPredictor::GetDistanceMapWithLabel(const typename ImageType::Pointer &labelImagePointer, const int &label1)
{
  typename ImageType::Pointer interImage = ImageType::New();
  interImage->CopyInformation(labelImagePointer);
  interImage->SetRequestedRegion(labelImagePointer->GetLargestPossibleRegion());
  interImage->SetBufferedRegion(labelImagePointer->GetBufferedRegion());
  interImage->Allocate();
  interImage->FillBuffer(0);
  typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
  IteratorType imIt(labelImagePointer, labelImagePointer->GetLargestPossibleRegion());
  IteratorType interIt(interImage, interImage->GetLargestPossibleRegion());
  imIt.GoToBegin();
  interIt.GoToBegin();
  while (!interIt.IsAtEnd())
  {
    typename ImageType::IndexType index = imIt.GetIndex();
    if (imIt.Get() == label1)
      interIt.Set(CAPTK::VOXEL_STATUS::ON);
    else
      interIt.Set(CAPTK::VOXEL_STATUS::OFF);
    ++interIt;
    ++imIt;
  }
  typedef itk::Image<unsigned char, 3>  UnsignedCharImageType;
  typedef itk::CastImageFilter< ImageType, UnsignedCharImageType > CastFilterType;
  typename CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(interImage);
  typedef  itk::SignedMaurerDistanceMapImageFilter< UnsignedCharImageType, ImageType  > SignedMaurerDistanceMapImageFilterType;
  typename  SignedMaurerDistanceMapImageFilterType::Pointer distanceMapImageFilter = SignedMaurerDistanceMapImageFilterType::New();
  distanceMapImageFilter->SetInput(castFilter->GetOutput());
  distanceMapImageFilter->Update();
  typename ImageType::Pointer distanceimage = distanceMapImageFilter->GetOutput();
  return distanceimage;
}
//---------------------------------------------------------------------------------------

template<class ImageType>
VectorDouble  EGFRvIIIIndexPredictor::LoadTestData(const typename ImageType::Pointer &t1ceImagePointer,
  const typename ImageType::Pointer &t2flairImagePointer,
  const typename ImageType::Pointer &t1ImagePointer,
  const typename ImageType::Pointer &t2ImagePointer,
  const typename ImageType::Pointer &rcbvImagePointer,
  const typename ImageType::Pointer &psrImagePointer,
  const typename ImageType::Pointer &phImagePointer,
  const typename ImageType::Pointer &axImagePointer,
  const typename ImageType::Pointer &faImagePointer,
  const typename ImageType::Pointer &radImagePointer,
  const typename ImageType::Pointer &trImagePointer,
  const typename ImageType::Pointer &labelImagePointer, 
  const typename ImageType::Pointer &atlasImagePointer,
  const typename ImageType::Pointer &templateImagePointer9regions,
  const typename ImageType::Pointer &templateImagePointer21regions,
  const typename ImageType::Pointer &templateImagePointerAllregions,
  const typename ImageType::Pointer &POSAtlasImagePointer,
  const typename ImageType::Pointer &NEGAtlasImagePointer,
  const VariableSizeMatrixType &configuration)
{
  std::vector<typename ImageType::IndexType> edemaIndices;
  std::vector<typename ImageType::IndexType> etumorIndices;
  std::vector<typename ImageType::IndexType> tumorIndices;
  std::vector<typename ImageType::IndexType> necoreIndices;
  std::vector<typename ImageType::IndexType> brainIndices;
  std::vector<typename ImageType::IndexType> ventIndices;

  VectorDouble tumorIntensitiesT1;
  VectorDouble tumorIntensitiesT2;
  VectorDouble tumorIntensitiesT1CE;
  VectorDouble tumorIntensitiesT2Flair;
  VectorDouble tumorIntensitiesAX;
  VectorDouble tumorIntensitiesFA;
  VectorDouble tumorIntensitiesRAD;
  VectorDouble tumorIntensitiesTR;
  VectorDouble tumorIntensitiesRCBV;
  VectorDouble tumorIntensitiesPSR;
  VectorDouble tumorIntensitiesPH;

  VectorDouble edemaIntensitiesT1;
  VectorDouble edemaIntensitiesT2;
  VectorDouble edemaIntensitiesT1CE;
  VectorDouble edemaIntensitiesT2Flair;
  VectorDouble edemaIntensitiesAX;
  VectorDouble edemaIntensitiesFA;
  VectorDouble edemaIntensitiesRAD;
  VectorDouble edemaIntensitiesTR;
  VectorDouble edemaIntensitiesRCBV;
  VectorDouble edemaIntensitiesPSR;
  VectorDouble edemaIntensitiesPH;

  VectorDouble necoreIntensitiesT1;
  VectorDouble necoreIntensitiesT2;
  VectorDouble necoreIntensitiesT1CE;
  VectorDouble necoreIntensitiesT2Flair;
  VectorDouble necoreIntensitiesAX;
  VectorDouble necoreIntensitiesFA;
  VectorDouble necoreIntensitiesRAD;
  VectorDouble necoreIntensitiesTR;
  VectorDouble necoreIntensitiesRCBV;
  VectorDouble necoreIntensitiesPSR;
  VectorDouble necoreIntensitiesPH;


  typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
  IteratorType imIt(labelImagePointer, labelImagePointer->GetLargestPossibleRegion());
  imIt.GoToBegin();
  while (!imIt.IsAtEnd())
  {
    if (imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::TUMOR)
      etumorIndices.push_back(imIt.GetIndex());
    else if (imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::NONENHANCING)
      necoreIndices.push_back(imIt.GetIndex());
    else if (imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::EDEMA)
      edemaIndices.push_back(imIt.GetIndex());
    else
    {
    }
    if (imIt.Get() > CAPTK::GLISTR_OUTPUT_LABELS::ALL)
      brainIndices.push_back(imIt.GetIndex());
    if (imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::TUMOR || imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::NONENHANCING)
      tumorIndices.push_back(imIt.GetIndex());
    ++imIt;
  }
  for (unsigned int i = 0; i < etumorIndices.size(); i++)
  {
    tumorIntensitiesT1.push_back(std::round(t1ImagePointer.GetPointer()->GetPixel(etumorIndices[i])));
    tumorIntensitiesT2.push_back(std::round(t2ImagePointer.GetPointer()->GetPixel(etumorIndices[i])));
    tumorIntensitiesT1CE.push_back(std::round(t1ceImagePointer.GetPointer()->GetPixel(etumorIndices[i])));
    tumorIntensitiesT2Flair.push_back(std::round(t2flairImagePointer.GetPointer()->GetPixel(etumorIndices[i])));
    tumorIntensitiesAX.push_back(std::round(axImagePointer.GetPointer()->GetPixel(etumorIndices[i])));
    tumorIntensitiesFA.push_back(std::round(faImagePointer.GetPointer()->GetPixel(etumorIndices[i])));
    tumorIntensitiesRAD.push_back(std::round(radImagePointer.GetPointer()->GetPixel(etumorIndices[i])));
    tumorIntensitiesTR.push_back(std::round(trImagePointer.GetPointer()->GetPixel(etumorIndices[i])));
    tumorIntensitiesRCBV.push_back(std::round(rcbvImagePointer.GetPointer()->GetPixel(etumorIndices[i])));
    tumorIntensitiesPSR.push_back(std::round(psrImagePointer.GetPointer()->GetPixel(etumorIndices[i])));
    tumorIntensitiesPH.push_back(phImagePointer.GetPointer()->GetPixel(etumorIndices[i]));
  }

  for (unsigned int i = 0; i < edemaIndices.size(); i++)
  {
    edemaIntensitiesT1.push_back(std::round(t1ImagePointer.GetPointer()->GetPixel(edemaIndices[i])));
    edemaIntensitiesT2.push_back(std::round(t2ImagePointer.GetPointer()->GetPixel(edemaIndices[i])));
    edemaIntensitiesT1CE.push_back(std::round(t1ceImagePointer.GetPointer()->GetPixel(edemaIndices[i])));
    edemaIntensitiesT2Flair.push_back(std::round(t2flairImagePointer.GetPointer()->GetPixel(edemaIndices[i])));
    edemaIntensitiesAX.push_back(std::round(axImagePointer.GetPointer()->GetPixel(edemaIndices[i])));
    edemaIntensitiesFA.push_back(std::round(faImagePointer.GetPointer()->GetPixel(edemaIndices[i])));
    edemaIntensitiesRAD.push_back(std::round(radImagePointer.GetPointer()->GetPixel(edemaIndices[i])));
    edemaIntensitiesTR.push_back(std::round(trImagePointer.GetPointer()->GetPixel(edemaIndices[i])));
    edemaIntensitiesRCBV.push_back(std::round(rcbvImagePointer.GetPointer()->GetPixel(edemaIndices[i])));
    edemaIntensitiesPSR.push_back(std::round(psrImagePointer.GetPointer()->GetPixel(edemaIndices[i])));
    edemaIntensitiesPH.push_back(phImagePointer.GetPointer()->GetPixel(edemaIndices[i]));
  }

  for (unsigned int i = 0; i < necoreIndices.size(); i++)
  {
    necoreIntensitiesT1.push_back(std::round(t1ImagePointer.GetPointer()->GetPixel(necoreIndices[i])));
    necoreIntensitiesT2.push_back(std::round(t2ImagePointer.GetPointer()->GetPixel(necoreIndices[i])));
    necoreIntensitiesT1CE.push_back(std::round(t1ceImagePointer.GetPointer()->GetPixel(necoreIndices[i])));
    necoreIntensitiesT2Flair.push_back(std::round(t2flairImagePointer.GetPointer()->GetPixel(necoreIndices[i])));
    necoreIntensitiesAX.push_back(std::round(axImagePointer.GetPointer()->GetPixel(necoreIndices[i])));
    necoreIntensitiesFA.push_back(std::round(faImagePointer.GetPointer()->GetPixel(necoreIndices[i])));
    necoreIntensitiesRAD.push_back(std::round(radImagePointer.GetPointer()->GetPixel(necoreIndices[i])));
    necoreIntensitiesTR.push_back(std::round(trImagePointer.GetPointer()->GetPixel(necoreIndices[i])));
    necoreIntensitiesRCBV.push_back(std::round(rcbvImagePointer.GetPointer()->GetPixel(necoreIndices[i])));
    necoreIntensitiesPSR.push_back(std::round(psrImagePointer.GetPointer()->GetPixel(necoreIndices[i])));
    necoreIntensitiesPH.push_back(phImagePointer.GetPointer()->GetPixel(necoreIndices[i]));
  }

  //------------------------------------------------------get histogram and statistics features-------------------
  VectorDouble ETBinsT1 = NormalizeHistogramFeatures(GetHistogramFeatures(tumorIntensitiesT1, configuration(0, 0), configuration(0, 1), configuration(0, 2)), tumorIntensitiesT1.size());
  VectorDouble ETBinsT2 = NormalizeHistogramFeatures(GetHistogramFeatures(tumorIntensitiesT2, configuration(1, 0), configuration(1, 1), configuration(1, 2)), tumorIntensitiesT2.size());
  VectorDouble ETBinsT1CE = NormalizeHistogramFeatures(GetHistogramFeatures(tumorIntensitiesT1CE, configuration(2, 0), configuration(2, 1), configuration(2, 2)), tumorIntensitiesT1CE.size());
  VectorDouble ETBinsFlair = NormalizeHistogramFeatures(GetHistogramFeatures(tumorIntensitiesT2Flair, configuration(3, 0), configuration(3, 1), configuration(3, 2)), tumorIntensitiesT2Flair.size());
  VectorDouble ETBinsAX = NormalizeHistogramFeatures(GetHistogramFeatures(tumorIntensitiesAX, configuration(4, 0), configuration(4, 1), configuration(4, 2)), tumorIntensitiesAX.size());
  VectorDouble ETBinsFA = NormalizeHistogramFeatures(GetHistogramFeatures(tumorIntensitiesFA, configuration(5, 0), configuration(5, 1), configuration(5, 2)), tumorIntensitiesFA.size());
  VectorDouble ETBinsRAD = NormalizeHistogramFeatures(GetHistogramFeatures(tumorIntensitiesRAD, configuration(6, 0), configuration(6, 1), configuration(6, 2)), tumorIntensitiesRAD.size());
  VectorDouble ETBinsTR = NormalizeHistogramFeatures(GetHistogramFeatures(tumorIntensitiesTR, configuration(7, 0), configuration(7, 1), configuration(7, 2)), tumorIntensitiesTR.size());
  VectorDouble ETBinsRCBV = NormalizeHistogramFeatures(GetHistogramFeatures(tumorIntensitiesRCBV, configuration(8, 0), configuration(8, 1), configuration(8, 2)), tumorIntensitiesRCBV.size());
  VectorDouble ETBinsPSR = NormalizeHistogramFeatures(GetHistogramFeatures(tumorIntensitiesPSR, configuration(9, 0), configuration(9, 1), configuration(9, 2)), tumorIntensitiesPSR.size());
  VectorDouble ETBinsPH = NormalizeHistogramFeatures(GetHistogramFeatures(tumorIntensitiesPH, configuration(10, 0), configuration(10, 1), configuration(10, 2)), tumorIntensitiesPH.size());

  VectorDouble EDBinsT1 = NormalizeHistogramFeatures(GetHistogramFeatures(edemaIntensitiesT1, configuration(11, 0), configuration(11, 1), configuration(11, 2)), edemaIntensitiesT1.size());
  VectorDouble EDBinsT2 = NormalizeHistogramFeatures(GetHistogramFeatures(edemaIntensitiesT2, configuration(12, 0), configuration(12, 1), configuration(12, 2)), edemaIntensitiesT2.size());
  VectorDouble EDBinsT1CE = NormalizeHistogramFeatures(GetHistogramFeatures(edemaIntensitiesT1CE, configuration(13, 0), configuration(13, 1), configuration(13, 2)), edemaIntensitiesT1CE.size());
  VectorDouble EDBinsFlair = NormalizeHistogramFeatures(GetHistogramFeatures(edemaIntensitiesT2Flair, configuration(14, 0), configuration(14, 1), configuration(14, 2)), edemaIntensitiesT2Flair.size());
  VectorDouble EDBinsAX = NormalizeHistogramFeatures(GetHistogramFeatures(edemaIntensitiesAX, configuration(15, 0), configuration(15, 1), configuration(15, 2)), edemaIntensitiesAX.size());
  VectorDouble EDBinsFA = NormalizeHistogramFeatures(GetHistogramFeatures(edemaIntensitiesFA, configuration(16, 0), configuration(16, 1), configuration(16, 2)), edemaIntensitiesFA.size());
  VectorDouble EDBinsRAD = NormalizeHistogramFeatures(GetHistogramFeatures(edemaIntensitiesRAD, configuration(17, 0), configuration(17, 1), configuration(17, 2)), edemaIntensitiesRAD.size());
  VectorDouble EDBinsTR = NormalizeHistogramFeatures(GetHistogramFeatures(edemaIntensitiesTR, configuration(18, 0), configuration(18, 1), configuration(18, 2)), edemaIntensitiesTR.size());
  VectorDouble EDBinsRCBV = NormalizeHistogramFeatures(GetHistogramFeatures(edemaIntensitiesRCBV, configuration(19, 0), configuration(19, 1), configuration(19, 2)), edemaIntensitiesRCBV.size());
  VectorDouble EDBinsPSR = NormalizeHistogramFeatures(GetHistogramFeatures(edemaIntensitiesPSR, configuration(20, 0), configuration(20, 1), configuration(20, 2)), edemaIntensitiesPSR.size());
  VectorDouble EDBinsPH = NormalizeHistogramFeatures(GetHistogramFeatures(edemaIntensitiesPH, configuration(21, 0), configuration(21, 1), configuration(21, 2)), edemaIntensitiesPH.size());

  VectorDouble NCRBinsT1 = NormalizeHistogramFeatures(GetHistogramFeatures(necoreIntensitiesT1, configuration(22, 0), configuration(22, 1), configuration(22, 2)), necoreIntensitiesT1.size());
  VectorDouble NCRBinsT2 = NormalizeHistogramFeatures(GetHistogramFeatures(necoreIntensitiesT2, configuration(23, 0), configuration(23, 1), configuration(23, 2)), necoreIntensitiesT2.size());
  VectorDouble NCRBinsT1CE = NormalizeHistogramFeatures(GetHistogramFeatures(necoreIntensitiesT1CE, configuration(24, 0), configuration(24, 1), configuration(24, 2)), necoreIntensitiesT1CE.size());
  VectorDouble NCRBinsFlair = NormalizeHistogramFeatures(GetHistogramFeatures(necoreIntensitiesT2Flair, configuration(25, 0), configuration(25, 1), configuration(25, 2)), necoreIntensitiesT2Flair.size());
  VectorDouble NCRBinsAX = NormalizeHistogramFeatures(GetHistogramFeatures(necoreIntensitiesAX, configuration(26, 0), configuration(26, 1), configuration(26, 2)), necoreIntensitiesAX.size());
  VectorDouble NCRBinsFA = NormalizeHistogramFeatures(GetHistogramFeatures(necoreIntensitiesFA, configuration(27, 0), configuration(27, 1), configuration(27, 2)), necoreIntensitiesFA.size());
  VectorDouble NCRBinsRAD = NormalizeHistogramFeatures(GetHistogramFeatures(necoreIntensitiesRAD, configuration(28, 0), configuration(28, 1), configuration(28, 2)), necoreIntensitiesRAD.size());
  VectorDouble NCRBinsTR = NormalizeHistogramFeatures(GetHistogramFeatures(necoreIntensitiesTR, configuration(29, 0), configuration(29, 1), configuration(29, 2)), necoreIntensitiesTR.size());
  VectorDouble NCRBinsRCBV = NormalizeHistogramFeatures(GetHistogramFeatures(necoreIntensitiesRCBV, configuration(30, 0), configuration(30, 1), configuration(30, 2)), necoreIntensitiesRCBV.size());
  VectorDouble NCRBinsPSR = NormalizeHistogramFeatures(GetHistogramFeatures(necoreIntensitiesPSR, configuration(31, 0), configuration(31, 1), configuration(31, 2)), necoreIntensitiesPSR.size());
  VectorDouble NCRBinsPH = NormalizeHistogramFeatures(GetHistogramFeatures(necoreIntensitiesPH, configuration(32, 0), configuration(32, 1), configuration(32, 2)), necoreIntensitiesPH.size());

  VectorDouble ETStatisticsT1 = GetStatisticalFeatures(tumorIntensitiesT1);
  VectorDouble ETStatisticsT2 = GetStatisticalFeatures(tumorIntensitiesT2);
  VectorDouble ETStatisticsT1CE = GetStatisticalFeatures(tumorIntensitiesT1CE);
  VectorDouble ETStatisticsFlair = GetStatisticalFeatures(tumorIntensitiesT2Flair);
  VectorDouble ETStatisticsAX = GetStatisticalFeatures(tumorIntensitiesAX);
  VectorDouble ETStatisticsFA = GetStatisticalFeatures(tumorIntensitiesFA);
  VectorDouble ETStatisticsRAD = GetStatisticalFeatures(tumorIntensitiesRAD);
  VectorDouble ETStatisticsTR = GetStatisticalFeatures(tumorIntensitiesTR);
  VectorDouble ETStatisticsRCBV = GetStatisticalFeatures(tumorIntensitiesRCBV);
  VectorDouble ETStatisticsPSR = GetStatisticalFeatures(tumorIntensitiesPSR);
  VectorDouble ETStatisticsPH = GetStatisticalFeatures(tumorIntensitiesPH);

  VectorDouble EDStatisticsT1 = GetStatisticalFeatures(edemaIntensitiesT1);
  VectorDouble EDStatisticsT2 = GetStatisticalFeatures(edemaIntensitiesT2);
  VectorDouble EDStatisticsT1CE = GetStatisticalFeatures(edemaIntensitiesT1CE);
  VectorDouble EDStatisticsFlair = GetStatisticalFeatures(edemaIntensitiesT2Flair);
  VectorDouble EDStatisticsAX = GetStatisticalFeatures(edemaIntensitiesAX);
  VectorDouble EDStatisticsFA = GetStatisticalFeatures(edemaIntensitiesFA);
  VectorDouble EDStatisticsRAD = GetStatisticalFeatures(edemaIntensitiesRAD);
  VectorDouble EDStatisticsTR = GetStatisticalFeatures(edemaIntensitiesTR);
  VectorDouble EDStatisticsRCBV = GetStatisticalFeatures(edemaIntensitiesRCBV);
  VectorDouble EDStatisticsPSR = GetStatisticalFeatures(edemaIntensitiesPSR);
  VectorDouble EDStatisticsPH = GetStatisticalFeatures(edemaIntensitiesPH);

  VectorDouble NCRStatisticsT1 = GetStatisticalFeatures(necoreIntensitiesT1);
  VectorDouble NCRStatisticsT2 = GetStatisticalFeatures(necoreIntensitiesT2);
  VectorDouble NCRStatisticsT1CE = GetStatisticalFeatures(necoreIntensitiesT1CE);
  VectorDouble NCRStatisticsFlair = GetStatisticalFeatures(necoreIntensitiesT2Flair);
  VectorDouble NCRStatisticsAX = GetStatisticalFeatures(necoreIntensitiesAX);
  VectorDouble NCRStatisticsFA = GetStatisticalFeatures(necoreIntensitiesFA);
  VectorDouble NCRStatisticsRAD = GetStatisticalFeatures(necoreIntensitiesRAD);
  VectorDouble NCRStatisticsTR = GetStatisticalFeatures(necoreIntensitiesTR);
  VectorDouble NCRStatisticsRCBV = GetStatisticalFeatures(necoreIntensitiesRCBV);
  VectorDouble NCRStatisticsPSR = GetStatisticalFeatures(necoreIntensitiesPSR);
  VectorDouble NCRStatisticsPH = GetStatisticalFeatures(necoreIntensitiesPH);

  //find connected components and remove components smaller than a certain area threshold
  std::vector<typename ImageType::Pointer> RevisedImages = RevisedTumorArea<ImageType>(labelImagePointer,20);
  typename ImageType::Pointer tumorImage = RevisedImages[0];
  typename ImageType::Pointer etumorImage = RevisedImages[1];
  typename ImageType::Pointer ncrImage = RevisedImages[2];


  //find indices of tumor sub-regions and calculate volumetric features in the patient space
  necoreIndices.clear();
  tumorIndices.clear();
  
  IteratorType ncrIt(ncrImage, ncrImage->GetLargestPossibleRegion());
  ncrIt.GoToBegin();
  while (!ncrIt.IsAtEnd())
  {
    if (ncrIt.Get() == CAPTK::VOXEL_STATUS::ON)
      necoreIndices.push_back(ncrIt.GetIndex());
    ++ncrIt;
  }
  IteratorType etIt(etumorImage, etumorImage->GetLargestPossibleRegion());
  etIt.GoToBegin();
  while (!etIt.IsAtEnd())
  {
    if (etIt.Get() == CAPTK::VOXEL_STATUS::ON)
      etumorIndices.push_back(etIt.GetIndex());
    ++etIt;
  }
  IteratorType tumorIt(tumorImage, tumorImage->GetLargestPossibleRegion());
  tumorIt.GoToBegin();
  while (!tumorIt.IsAtEnd())
  {
    if (tumorIt.Get() == CAPTK::VOXEL_STATUS::ON)
      tumorIndices.push_back(tumorIt.GetIndex());
    ++tumorIt;
  }
  VectorDouble VolumetricFeatures_Segmentation = GetVolumetricFeatures(edemaIndices.size(), etumorIndices.size(), necoreIndices.size(), brainIndices.size());

  //find indices of tumor sub-regions and calculate volumetric features in the atlas space
  etumorIndices.clear();
  necoreIndices.clear();
  tumorIndices.clear();
  edemaIndices.clear();
  brainIndices.clear();
  ventIndices.clear();
  
  typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
  IteratorType atlimIt(atlasImagePointer, atlasImagePointer->GetLargestPossibleRegion());
  atlimIt.GoToBegin();
  while (!atlimIt.IsAtEnd())
  {
    //if (atlimIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::TUMOR)
    //  etumorIndices.push_back(atlimIt.GetIndex());
    //else if (atlimIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::NONENHANCING)
    //  necoreIndices.push_back(atlimIt.GetIndex());
    if (atlimIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::EDEMA)
      edemaIndices.push_back(atlimIt.GetIndex());
    else if (atlimIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::VENT)
      ventIndices.push_back(atlimIt.GetIndex());
    else
    {
    }
    if (atlimIt.Get() > CAPTK::GLISTR_OUTPUT_LABELS::ALL)
      brainIndices.push_back(atlimIt.GetIndex());
    //if (atlimIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::TUMOR || atlimIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::NONENHANCING)
    //  tumorIndices.push_back(atlimIt.GetIndex());
    ++atlimIt;
  }
  RevisedImages = RevisedTumorArea<ImageType>(atlasImagePointer,10);
  tumorImage = RevisedImages[0];
  etumorImage = RevisedImages[1];
  ncrImage = RevisedImages[2];

  IteratorType ncrAtlasIt(ncrImage, ncrImage->GetLargestPossibleRegion());
  ncrAtlasIt.GoToBegin();
  while (!ncrAtlasIt.IsAtEnd())
  {
    if (ncrAtlasIt.Get() == CAPTK::VOXEL_STATUS::ON)
      necoreIndices.push_back(ncrIt.GetIndex());
    ++ncrAtlasIt;
  }
  IteratorType etAtlasIt(etumorImage, etumorImage->GetLargestPossibleRegion());
  etAtlasIt.GoToBegin();
  while (!etAtlasIt.IsAtEnd())
  {
    if (etAtlasIt.Get() == CAPTK::VOXEL_STATUS::ON)
      etumorIndices.push_back(etAtlasIt.GetIndex());
    ++etAtlasIt;
  }
  IteratorType tumorAtlasIt(tumorImage, tumorImage->GetLargestPossibleRegion());
  tumorAtlasIt.GoToBegin();
  while (!tumorAtlasIt.IsAtEnd())
  {
    if (tumorAtlasIt.Get() == CAPTK::VOXEL_STATUS::ON)
      tumorIndices.push_back(tumorAtlasIt.GetIndex());
    ++tumorAtlasIt;
  }
  VectorDouble VolumetricFeatures_Atlas = GetVolumetricFeatures(edemaIndices.size(), etumorIndices.size(), necoreIndices.size(), brainIndices.size());

  //calculate distance maps and use these maps to calculate distance features
  typename ImageType::Pointer  edemaDistanceMap = GetDistanceMapWithLabel<ImageType>(atlasImagePointer, CAPTK::GLISTR_OUTPUT_LABELS::EDEMA);
  typename ImageType::Pointer ventImage = ImageType::New();
  ventImage->CopyInformation(atlasImagePointer);
  ventImage->SetRequestedRegion(atlasImagePointer->GetLargestPossibleRegion());
  ventImage->SetBufferedRegion(atlasImagePointer->GetBufferedRegion());
  ventImage->Allocate();
  ventImage->FillBuffer(0);

  typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;

  IteratorType labelIt(templateImagePointerAllregions, templateImagePointerAllregions->GetLargestPossibleRegion());
  IteratorType ventIt(ventImage, ventImage->GetLargestPossibleRegion());
  labelIt.GoToBegin();
  ventIt.GoToBegin();
  while (!labelIt.IsAtEnd())
  {
    if (labelIt.Get() == CAPTK::ORIGINAL_ATLAS_LABELS::VENT1 || 
      labelIt.Get() == CAPTK::ORIGINAL_ATLAS_LABELS::VENT2 || 
      labelIt.Get() == CAPTK::ORIGINAL_ATLAS_LABELS::VENT3 || 
      labelIt.Get() == CAPTK::ORIGINAL_ATLAS_LABELS::VENT4)
      ventIt.Set(CAPTK::VOXEL_STATUS::ON);
    ++labelIt;
    ++ventIt;
  }
  typename ImageType::Pointer  ventDistanceMap = GetDistanceMap<ImageType>(ventImage);
  typename ImageType::Pointer  tumorDistanceMap = GetDistanceMap<ImageType>(tumorImage);
  VectorDouble DistanceFeatures = GetDistanceFeatures3<ImageType>(edemaDistanceMap, tumorDistanceMap, ventDistanceMap);

  //spatial location features
  std::vector<double> spatialLocationFeatures_4 = GetSpatialLocationFeaturesForEGFR<ImageType>(atlasImagePointer, NEGAtlasImagePointer, POSAtlasImagePointer);
  VectorDouble spatialLocationFeatures_9        = GetSpatialLocationFeaturesForFixedNumberOfRegions<ImageType>(atlasImagePointer, templateImagePointer9regions,9);
  VectorDouble spatialLocationFeatures_21       = GetSpatialLocationFeaturesForFixedNumberOfRegions<ImageType>(atlasImagePointer, templateImagePointer21regions,21);

  //copy data from vectors to one final feature vector
  VectorDouble TestFeatures;
  //copy location features
  TestFeatures.insert(TestFeatures.end(), spatialLocationFeatures_4.begin(), spatialLocationFeatures_4.end());

  //copy volumetric features
  TestFeatures.insert(TestFeatures.end(), VolumetricFeatures_Segmentation.begin(), VolumetricFeatures_Segmentation.end());
  TestFeatures.insert(TestFeatures.end(), VolumetricFeatures_Atlas.begin(), VolumetricFeatures_Atlas.end());

  //copy distance features
  TestFeatures.insert(TestFeatures.end(), DistanceFeatures.begin(), DistanceFeatures.end());

  //copy statistics features
  TestFeatures.insert(TestFeatures.end(), ETStatisticsT1CE.begin(), ETStatisticsT1CE.end());
  TestFeatures.insert(TestFeatures.end(), ETStatisticsT1.begin(), ETStatisticsT1.end());
  TestFeatures.insert(TestFeatures.end(), ETStatisticsT2.begin(), ETStatisticsT2.end());
  TestFeatures.insert(TestFeatures.end(), ETStatisticsFlair.begin(), ETStatisticsFlair.end());
  TestFeatures.insert(TestFeatures.end(), ETStatisticsPH.begin(), ETStatisticsPH.end());
  TestFeatures.insert(TestFeatures.end(), ETStatisticsPSR.begin(), ETStatisticsPSR.end());
  TestFeatures.insert(TestFeatures.end(), ETStatisticsRCBV.begin(), ETStatisticsRCBV.end());
  TestFeatures.insert(TestFeatures.end(), ETStatisticsFA.begin(), ETStatisticsFA.end());
  TestFeatures.insert(TestFeatures.end(), ETStatisticsAX.begin(), ETStatisticsAX.end());
  TestFeatures.insert(TestFeatures.end(), ETStatisticsRAD.begin(), ETStatisticsRAD.end());
  TestFeatures.insert(TestFeatures.end(), ETStatisticsTR.begin(), ETStatisticsTR.end());

  TestFeatures.insert(TestFeatures.end(), NCRStatisticsT1CE.begin(), NCRStatisticsT1CE.end());
  TestFeatures.insert(TestFeatures.end(), NCRStatisticsT1.begin(), NCRStatisticsT1.end());
  TestFeatures.insert(TestFeatures.end(), NCRStatisticsT2.begin(), NCRStatisticsT2.end());
  TestFeatures.insert(TestFeatures.end(), NCRStatisticsFlair.begin(), NCRStatisticsFlair.end());
  TestFeatures.insert(TestFeatures.end(), NCRStatisticsPH.begin(), NCRStatisticsPH.end());
  TestFeatures.insert(TestFeatures.end(), NCRStatisticsPSR.begin(), NCRStatisticsPSR.end());
  TestFeatures.insert(TestFeatures.end(), NCRStatisticsRCBV.begin(), NCRStatisticsRCBV.end());
  TestFeatures.insert(TestFeatures.end(), NCRStatisticsFA.begin(), NCRStatisticsFA.end());
  TestFeatures.insert(TestFeatures.end(), NCRStatisticsAX.begin(), NCRStatisticsAX.end());
  TestFeatures.insert(TestFeatures.end(), NCRStatisticsRAD.begin(), NCRStatisticsRAD.end());
  TestFeatures.insert(TestFeatures.end(), NCRStatisticsTR.begin(), NCRStatisticsTR.end());

  TestFeatures.insert(TestFeatures.end(), EDStatisticsT1CE.begin(), EDStatisticsT1CE.end());
  TestFeatures.insert(TestFeatures.end(), EDStatisticsT1.begin(), EDStatisticsT1.end());
  TestFeatures.insert(TestFeatures.end(), EDStatisticsT2.begin(), EDStatisticsT2.end());
  TestFeatures.insert(TestFeatures.end(), EDStatisticsFlair.begin(), EDStatisticsFlair.end());
  TestFeatures.insert(TestFeatures.end(), EDStatisticsPH.begin(), EDStatisticsPH.end());
  TestFeatures.insert(TestFeatures.end(), EDStatisticsPSR.begin(), EDStatisticsPSR.end());
  TestFeatures.insert(TestFeatures.end(), EDStatisticsRCBV.begin(), EDStatisticsRCBV.end());
  TestFeatures.insert(TestFeatures.end(), EDStatisticsFA.begin(), EDStatisticsFA.end());
  TestFeatures.insert(TestFeatures.end(), EDStatisticsAX.begin(), EDStatisticsAX.end());
  TestFeatures.insert(TestFeatures.end(), EDStatisticsRAD.begin(), EDStatisticsRAD.end());
  TestFeatures.insert(TestFeatures.end(), EDStatisticsTR.begin(), EDStatisticsTR.end());

  //copy histogram features
  TestFeatures.insert(TestFeatures.end(), ETBinsT1CE.begin(), ETBinsT1CE.end());
  TestFeatures.insert(TestFeatures.end(), EDBinsT1CE.begin(), EDBinsT1CE.end());
  
  //doing it differently for this feature type as 8,9,10 number bins need to be excluded
  for (unsigned int histoIndex = 0; histoIndex<7; histoIndex++)
    TestFeatures.push_back(NCRBinsT1CE[histoIndex]);

  TestFeatures.insert(TestFeatures.end(), ETBinsT1.begin(), ETBinsT1.end());
  TestFeatures.insert(TestFeatures.end(), EDBinsT1.begin(), EDBinsT1.end());
  TestFeatures.insert(TestFeatures.end(), NCRBinsT1.begin(), NCRBinsT1.end());
  TestFeatures.insert(TestFeatures.end(), ETBinsT2.begin(), ETBinsT2.end());
  TestFeatures.insert(TestFeatures.end(), EDBinsT2.begin(), EDBinsT2.end());
  TestFeatures.insert(TestFeatures.end(), NCRBinsT2.begin(), NCRBinsT2.end());
  TestFeatures.insert(TestFeatures.end(), ETBinsFlair.begin(), ETBinsFlair.end());
  TestFeatures.insert(TestFeatures.end(), EDBinsFlair.begin(), EDBinsFlair.end());
  TestFeatures.insert(TestFeatures.end(), NCRBinsFlair.begin(), NCRBinsFlair.end());
  TestFeatures.insert(TestFeatures.end(), ETBinsPH.begin(), ETBinsPH.end());
  TestFeatures.insert(TestFeatures.end(), EDBinsPH.begin(), EDBinsPH.end());
  TestFeatures.insert(TestFeatures.end(), NCRBinsPH.begin(), NCRBinsPH.end());
  TestFeatures.insert(TestFeatures.end(), ETBinsPSR.begin(), ETBinsPSR.end());
  TestFeatures.insert(TestFeatures.end(), EDBinsPSR.begin(), EDBinsPSR.end());
  TestFeatures.insert(TestFeatures.end(), NCRBinsPSR.begin(), NCRBinsPSR.end());
  TestFeatures.insert(TestFeatures.end(), ETBinsRCBV.begin(), ETBinsRCBV.end());
  TestFeatures.insert(TestFeatures.end(), EDBinsRCBV.begin(), EDBinsRCBV.end());
  TestFeatures.insert(TestFeatures.end(), NCRBinsRCBV.begin(), NCRBinsRCBV.end());
  TestFeatures.insert(TestFeatures.end(), ETBinsFA.begin(), ETBinsFA.end());
  TestFeatures.insert(TestFeatures.end(), EDBinsFA.begin(), EDBinsFA.end());
  TestFeatures.insert(TestFeatures.end(), NCRBinsFA.begin(), NCRBinsFA.end());
  TestFeatures.insert(TestFeatures.end(), ETBinsAX.begin(), ETBinsAX.end());
  TestFeatures.insert(TestFeatures.end(), EDBinsAX.begin(), EDBinsAX.end());
  TestFeatures.insert(TestFeatures.end(), NCRBinsAX.begin(), NCRBinsAX.end());
  TestFeatures.insert(TestFeatures.end(), ETBinsRAD.begin(), ETBinsRAD.end());
  TestFeatures.insert(TestFeatures.end(), EDBinsRAD.begin(), EDBinsRAD.end());
  TestFeatures.insert(TestFeatures.end(), NCRBinsRAD.begin(), NCRBinsRAD.end());
  TestFeatures.insert(TestFeatures.end(), ETBinsTR.begin(), ETBinsTR.end());
  TestFeatures.insert(TestFeatures.end(), EDBinsTR.begin(), EDBinsTR.end());
  TestFeatures.insert(TestFeatures.end(), NCRBinsTR.begin(), NCRBinsTR.end());

  //copy spatial location features
  TestFeatures.insert(TestFeatures.end(), spatialLocationFeatures_9.begin(), spatialLocationFeatures_9.end());

  //doing it differently for this feature type as ROI number 17 and 21 need to be excluded
  for (unsigned int roiIndex = 0; roiIndex < 21; roiIndex++)
  {
    if (roiIndex == 16 || roiIndex == 20)
      continue;
    else
      TestFeatures.push_back(spatialLocationFeatures_21[roiIndex]);
  }
  return TestFeatures;
}

template<class ImageType>
std::vector<typename ImageType::Pointer> EGFRvIIIIndexPredictor::RevisedTumorArea(const typename ImageType::Pointer &labelImagePointer,const int areathreshold)
{
  std::vector<typename ImageType::Pointer> RevisedImages;

  typename ImageType::Pointer tumorImage = ImageType::New();
  tumorImage->CopyInformation(labelImagePointer);
  tumorImage->SetRequestedRegion(labelImagePointer->GetLargestPossibleRegion());
  tumorImage->SetBufferedRegion(labelImagePointer->GetBufferedRegion());
  tumorImage->Allocate();
  tumorImage->FillBuffer(0);

  typename ImageType::Pointer etumorImage = ImageType::New();
  etumorImage->CopyInformation(labelImagePointer);
  etumorImage->SetRequestedRegion(labelImagePointer->GetLargestPossibleRegion());
  etumorImage->SetBufferedRegion(labelImagePointer->GetBufferedRegion());
  etumorImage->Allocate();
  etumorImage->FillBuffer(0);

  typename ImageType::Pointer ncrImage = ImageType::New();
  ncrImage->CopyInformation(labelImagePointer);
  ncrImage->SetRequestedRegion(labelImagePointer->GetLargestPossibleRegion());
  ncrImage->SetBufferedRegion(labelImagePointer->GetBufferedRegion());
  ncrImage->Allocate();
  ncrImage->FillBuffer(0);


  typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
  IteratorType imIt(labelImagePointer, labelImagePointer->GetLargestPossibleRegion());
  IteratorType ncrIt(ncrImage, ncrImage->GetLargestPossibleRegion());
  IteratorType etIt(etumorImage, etumorImage->GetLargestPossibleRegion());
  IteratorType tumorIt(tumorImage, tumorImage->GetLargestPossibleRegion());
  imIt.GoToBegin();
  tumorIt.GoToBegin();
  etIt.GoToBegin();
  ncrIt.GoToBegin();

  while (!tumorIt.IsAtEnd())
  {
    if (imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::TUMOR || imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::NONENHANCING)
      tumorIt.Set(CAPTK::VOXEL_STATUS::ON);
    else
      tumorIt.Set(CAPTK::VOXEL_STATUS::OFF);

    if (imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::TUMOR)
      etIt.Set(CAPTK::VOXEL_STATUS::ON);
    else
      etIt.Set(CAPTK::VOXEL_STATUS::OFF);

    if (imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::NONENHANCING)
      ncrIt.Set(CAPTK::VOXEL_STATUS::ON);
    else
      ncrIt.Set(CAPTK::VOXEL_STATUS::OFF);

    ++tumorIt;
    ++etIt;
    ++imIt;
    ++ncrIt;
  }

  typedef itk::Image< unsigned short, 3 > OutputImageType;

  typedef itk::ConnectedComponentImageFilter <ImageType, OutputImageType> ConnectedComponentImageFilterType;
  typename ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
  connected->FullyConnectedOn();
  connected->SetInput(tumorImage);
  connected->Update();
  OutputImageType::Pointer labeledImage = connected->GetOutput();

  connected->GetObjectCount();
  std::vector<int> sizes;
  typedef itk::ImageRegionIteratorWithIndex <OutputImageType> OutputIteratorType;
  OutputIteratorType lbimIt(labeledImage, labeledImage->GetLargestPossibleRegion());

  for (unsigned int i = 0; i < connected->GetObjectCount(); i++)
  {
    int counter = 0;
    lbimIt.GoToBegin();
    while (!lbimIt.IsAtEnd())
    {
      if (lbimIt.Get() == i + 1)
        counter++;
      ++lbimIt;
    }
    sizes.push_back(counter);
  }
  for (unsigned int i = 0; i < connected->GetObjectCount(); i++)
  {
    if (sizes[i] <= areathreshold)
    {
      lbimIt.GoToBegin();
      tumorIt.GoToBegin();
      while (!lbimIt.IsAtEnd())
      {
        if (lbimIt.Get() == i + 1)
          tumorIt.Set(CAPTK::VOXEL_STATUS::OFF);
        ++lbimIt;
        ++tumorIt;
      }
    }
  }
  tumorIt.GoToBegin();
  etIt.GoToBegin();
  ncrIt.GoToBegin();

  while (!tumorIt.IsAtEnd())
  {
    etIt.Set(etIt.Get()*tumorIt.Get());
    ncrIt.Set(ncrIt.Get()*tumorIt.Get());

    ++tumorIt;
    ++etIt;
    ++ncrIt;
  }
  RevisedImages.push_back(tumorImage);
  RevisedImages.push_back(etumorImage);
  RevisedImages.push_back(ncrImage);

  return RevisedImages;
}

template<class ImageType>
VectorDouble EGFRvIIIIndexPredictor::GetSpatialLocationFeaturesForFixedNumberOfRegions(const typename ImageType::Pointer &labelImagePointer, const typename ImageType::Pointer &jacobtemplateImagePointer,const int numberofregions)
{
  std::vector<typename ImageType::Pointer> RevisedImages = RevisedTumorArea<ImageType>(labelImagePointer,10);
  typename ImageType::Pointer tumorImage = RevisedImages[0];
  typename ImageType::Pointer etumorImage = RevisedImages[1];
  typename ImageType::Pointer ncrImage = RevisedImages[2];


  typename ImageType::Pointer localizeImage = ImageType::New();
  localizeImage->CopyInformation(labelImagePointer);
  localizeImage->SetRequestedRegion(labelImagePointer->GetLargestPossibleRegion());
  localizeImage->SetBufferedRegion(labelImagePointer->GetBufferedRegion());
  localizeImage->Allocate();
  localizeImage->FillBuffer(0);

  //mulitply revised tumor image with the atlas ROI image
  typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
  IteratorType tumorIt(tumorImage, tumorImage->GetLargestPossibleRegion());
  IteratorType localizeIt(localizeImage, localizeImage->GetLargestPossibleRegion());
  IteratorType atlasIt(jacobtemplateImagePointer, jacobtemplateImagePointer->GetLargestPossibleRegion());

  atlasIt.GoToBegin();
  tumorIt.GoToBegin();
  localizeIt.GoToBegin();

  while (!tumorIt.IsAtEnd())
  {
    localizeIt.Set(tumorIt.Get()*atlasIt.Get());
    ++localizeIt;
    ++tumorIt;
    ++atlasIt;
  }
  //find number of voxels in pre-defined number of ROIs
  VectorDouble location;
  int tumorSize = 0;
  for (int i = 0; i < numberofregions; i++)
  {
    int counter = 0;
    localizeIt.GoToBegin();
    while (!localizeIt.IsAtEnd())
    {
      if (localizeIt.Get() == i + 1)
      {
        counter++;
        tumorSize++;
      }
      ++localizeIt;
    }
    location.push_back(counter);
  }
  //find percentages in predefined number of ROIs
  for (int i = 0; i < numberofregions; i++)
    location[i] = (location[i] * 100) / tumorSize;

  return location;
}

template<class ImageType>
std::vector<double> EGFRvIIIIndexPredictor::GetSpatialLocationFeaturesForEGFR(typename ImageType::Pointer labelImagePointer,
  typename ImageType::Pointer EGFRvIII_Neg_Atlas,
  typename ImageType::Pointer EGFRvIII_Pos_Atlas)
{
  typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;

  //find max and min 
  IteratorType negAtlasIt(EGFRvIII_Neg_Atlas, EGFRvIII_Neg_Atlas->GetLargestPossibleRegion());
  IteratorType posAtlasIt(EGFRvIII_Pos_Atlas, EGFRvIII_Pos_Atlas->GetLargestPossibleRegion());
  negAtlasIt.GoToBegin();
  posAtlasIt.GoToBegin();
  double maxposvalue = 0;
  double maxnegvalue = 0;
  while (!negAtlasIt.IsAtEnd())
  {
    if (negAtlasIt.Get() > maxnegvalue)
      maxnegvalue = negAtlasIt.Get();
    if (posAtlasIt.Get() > maxposvalue)
      maxposvalue = posAtlasIt.Get();
    ++negAtlasIt;
    ++posAtlasIt;
  }

  //multiply with 100
  negAtlasIt.GoToBegin();
  posAtlasIt.GoToBegin();
  while (!negAtlasIt.IsAtEnd())
  {
    negAtlasIt.Set(negAtlasIt.Get() * 100 / maxnegvalue);
    posAtlasIt.Set(posAtlasIt.Get() * 100 / maxposvalue);
    ++negAtlasIt;
    ++posAtlasIt;
  }
  //make new tumor with one ony;
  typename ImageType::Pointer tumorImage = ImageType::New();
  tumorImage->CopyInformation(labelImagePointer);
  tumorImage->SetRequestedRegion(labelImagePointer->GetLargestPossibleRegion());
  tumorImage->SetBufferedRegion(labelImagePointer->GetBufferedRegion());
  tumorImage->Allocate();
  tumorImage->FillBuffer(0);

  IteratorType tumorIt(tumorImage, tumorImage->GetLargestPossibleRegion());
  IteratorType imIt(labelImagePointer, labelImagePointer->GetLargestPossibleRegion());
  imIt.GoToBegin();
  tumorIt.GoToBegin();

  while (!tumorIt.IsAtEnd())
  {
    if (imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::TUMOR || imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::NONENHANCING)
      tumorIt.Set(CAPTK::VOXEL_STATUS::ON);
    else
      tumorIt.Set(CAPTK::VOXEL_STATUS::OFF);
    ++tumorIt;
    ++imIt;
  }

  //make new tumor with one ony;

  tumorIt.GoToBegin();
  negAtlasIt.GoToBegin();
  posAtlasIt.GoToBegin();
  std::vector<double> probNEGAtlas;
  std::vector<double> probPOSAtlas;

  while (!tumorIt.IsAtEnd())
  {
    if (tumorIt.Get() == CAPTK::VOXEL_STATUS::ON)
    {
//      if (negAtlasIt.Get() > 0)
        probNEGAtlas.push_back(negAtlasIt.Get());
//      if (posAtlasIt.Get() > 0)
        probPOSAtlas.push_back(posAtlasIt.Get());
    }
    ++tumorIt;
    ++negAtlasIt;
    ++posAtlasIt;
  }

  double sumPOS = 0;
  double sumNEG = 0;
  double average_POS_Atlas = 0;
  double average_NEG_Atlas = 0;
  double max_POS_Atlas=0;
  double max_NEG_Atlas=0;

  for (int i = 0; i < probPOSAtlas.size(); i++)
    sumPOS = sumPOS + probPOSAtlas[i];

  for (int i = 0; i < probNEGAtlas.size(); i++)
    sumNEG = sumNEG + probNEGAtlas[i];

  if (probPOSAtlas.size() > 0)
  {
    average_POS_Atlas = sumPOS / probPOSAtlas.size();
    max_POS_Atlas = *std::max_element(std::begin(probPOSAtlas), std::end(probPOSAtlas));
  }
  if (probNEGAtlas.size() > 0)
  {
    average_NEG_Atlas = sumNEG / probNEGAtlas.size();
    max_NEG_Atlas = *std::max_element(std::begin(probNEGAtlas), std::end(probNEGAtlas));
  }

  //previous implementation of location features
  std::vector<double> SpatialFeatures;
  SpatialFeatures.push_back(average_NEG_Atlas-average_POS_Atlas);
  SpatialFeatures.push_back(max_NEG_Atlas- max_POS_Atlas);
  SpatialFeatures.push_back(average_NEG_Atlas/average_POS_Atlas);
  SpatialFeatures.push_back(max_NEG_Atlas/max_POS_Atlas);

  //std::vector<double> SpatialFeatures;
  //SpatialFeatures.push_back(average_POS_Atlas - average_NEG_Atlas);
  //SpatialFeatures.push_back(max_POS_Atlas - max_NEG_Atlas);
  //SpatialFeatures.push_back(average_POS_Atlas / average_NEG_Atlas);
  //SpatialFeatures.push_back(max_POS_Atlas / max_NEG_Atlas);

  std::cout << SpatialFeatures[0] << " " << SpatialFeatures[1] << " " << SpatialFeatures[2] << " " << SpatialFeatures[3] << std::endl;

  return SpatialFeatures;
}

template<class ImageType>
std::vector<double> EGFRvIIIIndexPredictor::GetSpatialLocationFeaturesForEGFRTrainingWithLOO(typename ImageType::Pointer labelImagePointer,
  typename ImageType::Pointer EGFRvIII_Neg_Atlas,
  typename ImageType::Pointer EGFRvIII_Pos_Atlas, int label)
{
  typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;

  //make new tumor with one ony;
  typename ImageType::Pointer tumorImage = ImageType::New();
  tumorImage->CopyInformation(labelImagePointer);
  tumorImage->SetRequestedRegion(labelImagePointer->GetLargestPossibleRegion());
  tumorImage->SetBufferedRegion(labelImagePointer->GetBufferedRegion());
  tumorImage->Allocate();
  tumorImage->FillBuffer(0);

  IteratorType tumorIt(tumorImage, tumorImage->GetLargestPossibleRegion());
  IteratorType imIt(labelImagePointer, labelImagePointer->GetLargestPossibleRegion());
  imIt.GoToBegin();
  tumorIt.GoToBegin();

  while (!tumorIt.IsAtEnd())
  {
    if (imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::TUMOR || imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::NONENHANCING)
      tumorIt.Set(CAPTK::VOXEL_STATUS::ON);
    else
      tumorIt.Set(CAPTK::VOXEL_STATUS::OFF);
    ++tumorIt;
    ++imIt;
  }

  //--------------------------------------------------
  IteratorType negAtlasIt(EGFRvIII_Neg_Atlas, EGFRvIII_Neg_Atlas->GetLargestPossibleRegion());
  IteratorType posAtlasIt(EGFRvIII_Pos_Atlas, EGFRvIII_Pos_Atlas->GetLargestPossibleRegion());
  tumorIt.GoToBegin();
  negAtlasIt.GoToBegin();
  posAtlasIt.GoToBegin();
  while (!tumorIt.IsAtEnd())
  {
    if (tumorIt.Get() == CAPTK::VOXEL_STATUS::ON && label == 0)
      negAtlasIt.Set(negAtlasIt.Get() - 1);
    else if (tumorIt.Get() == CAPTK::VOXEL_STATUS::ON && label == 1)
      posAtlasIt.Set(posAtlasIt.Get() - 1);
    else
    {
    }
    ++tumorIt;
    ++negAtlasIt;
    ++posAtlasIt;
  }
  //------------------------------------------------------   

  //find max and min 
  negAtlasIt.GoToBegin();
  posAtlasIt.GoToBegin();
  double maxposvalue = 0;
  double maxnegvalue = 0;
  while (!negAtlasIt.IsAtEnd())
  {
    if (negAtlasIt.Get() > maxnegvalue)
      maxnegvalue = negAtlasIt.Get();
    if (posAtlasIt.Get() > maxposvalue)
      maxposvalue = posAtlasIt.Get();
    ++negAtlasIt;
    ++posAtlasIt;
  }

  //multiply with 100
  negAtlasIt.GoToBegin();
  posAtlasIt.GoToBegin();
  while (!negAtlasIt.IsAtEnd())
  {
    negAtlasIt.Set(negAtlasIt.Get() * 100 / maxnegvalue);
    posAtlasIt.Set(posAtlasIt.Get() * 100 / maxposvalue);
    ++negAtlasIt;
    ++posAtlasIt;
  }


  tumorIt.GoToBegin();
  negAtlasIt.GoToBegin();
  posAtlasIt.GoToBegin();
  std::vector<double> probNEGAtlas;
  std::vector<double> probPOSAtlas;

  while (!tumorIt.IsAtEnd())
  {
    if (tumorIt.Get() == CAPTK::VOXEL_STATUS::ON)
    {
      probNEGAtlas.push_back(negAtlasIt.Get());
      probPOSAtlas.push_back(posAtlasIt.Get());
    }
    ++tumorIt;
    ++negAtlasIt;
    ++posAtlasIt;
  }
  double max_POS_Atlas = *std::max_element(std::begin(probPOSAtlas), std::end(probPOSAtlas));
  double max_NEG_Atlas = *std::max_element(std::begin(probNEGAtlas), std::end(probNEGAtlas));

  double sumPOS = 0;
  double sumNEG = 0;

  for (int i = 0; i < probPOSAtlas.size(); i++)
  {
    sumPOS = sumPOS + probPOSAtlas[i];
    sumNEG = sumNEG + probNEGAtlas[i];
  }
  double average_POS_Atlas = sumPOS / probPOSAtlas.size();
  double average_NEG_Atlas = sumNEG / probNEGAtlas.size();

  std::vector<double> SpatialFeatures;
  SpatialFeatures.push_back(average_NEG_Atlas - average_POS_Atlas);
  SpatialFeatures.push_back(max_NEG_Atlas - max_POS_Atlas);
  SpatialFeatures.push_back(average_NEG_Atlas / average_POS_Atlas);
  SpatialFeatures.push_back(max_NEG_Atlas / max_POS_Atlas);

  return SpatialFeatures;
}
#endif







