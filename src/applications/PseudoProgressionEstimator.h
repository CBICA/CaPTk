/**
\file PseudoProgressionEstimator.h

This file holds the declaration of the class PseudoProgressionEstimator.

https://www.med.upenn.edu/sbia/software/ <br>
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
#include "CaPTkClassifierUtils.h"
#include "cbicaLogging.h"
#include "itkEnhancedScalarImageToRunLengthFeaturesFilter.h"
#include "itkRoundImageFilter.h"

#define RECURRENCE_MODEL_G 0.5
#define RECURRENCE_MODEL_RHO 0.0896
#define PSEUDOPROGRESSION_NO_OF_FEATURES 1041

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
    mPseudoTrainedFile = "Pseudo_SVM_Model.xml";
    mRecurrenceTrainedFile = "Recurrence_SVM_Model.xml";
    logger.UseNewFile(loggerFile);
  };

  //! Default Destructor
  ~PseudoProgressionEstimator();

  //! Gets last error
  std::string mLastEncounteredError;
  std::string mPseudoTrainedFile;
  std::string mRecurrenceTrainedFile;
  cbica::Logging logger;

  void WriteCSVFiles(VariableSizeMatrixType inputdata, std::string filepath);
  void WriteCSVFiles(vtkSmartPointer<vtkTable> inputdata, std::string filepath);
  void WriteCSVFiles(VectorVectorDouble inputdata, std::string filepath);
  void WriteCSVFiles(VariableLengthVectorType inputdata, std::string filepath);
  void WriteCSVFiles(std::vector<double> inputdata, std::string filepath);


  VariableSizeMatrixType GetModelSelectedFeatures(VariableSizeMatrixType & ScaledFeatureSetAfterAddingLabel, VariableLengthVectorType & psuSelectedFeatures);
  void WritePCAOutputs(std::string suffix, std::string outputdirectory, const VariableLengthVectorType mean, const VariableSizeMatrixType coefs);

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

  VectorDouble CombineEstimates(const VariableLengthVectorType &estimates1, const VariableLengthVectorType &estimates2);

  VariableLengthVectorType DistanceFunction(const VariableSizeMatrixType &testData, const std::string &filename, const double &rho, const double &bestg);
  /**
  \brief Get the size of the Feature Vector

  */
  int GetFeatureVectorSize(bool &useConvetionalData, bool &useDTIData, bool &usePerfData, bool &useDistData);

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

  template <typename T>
  std::vector<size_t> sort_indexes(const std::vector<T> &v);

  VariableSizeMatrixType SelectModelFeatures(const VariableSizeMatrixType &ModelFeatures);
  VectorDouble EffectSizeFeatureSelection(const VariableSizeMatrixType training_features, std::vector<double> target);

  template <class ImageType>
  typename ImageType::Pointer ReadNiftiImage(const std::string &filename);
  template<class ImageType>
  typename ImageType::Pointer RescaleImageIntensity(const typename ImageType::Pointer &image);

  VectorDouble GetHistogramFeatures(std::vector<float> intensities, int m_Bins);
  VectorDouble GetHistogramFeaturesWhole(std::vector<float> intensities);

  template< class TImageTypeShape = TImageType >
  VectorDouble GetShapeFeatures(typename TImageTypeShape::Pointer mask);

  VectorDouble GetIntensityFeatures(std::vector<float> m_intensities);

  template<class ImageType>
  VectorDouble GetGLCMFeatures(typename ImageType::Pointer image, typename ImageType::Pointer mask,double minvalue,double maxvalue);

  template<class ImageType>
  typename ImageType::Pointer MakeAdditionalModality(typename ImageType::Pointer image1, typename ImageType::Pointer image2);


  template<class ImageType>
  VectorDouble GetRunLengthFeatures(typename ImageType::Pointer image, typename ImageType::Pointer mask, double minvalue, double maxvalue);

  template<class PerfusionImageType, class ImageType>
  VectorVectorDouble GetPerfusionFeatures(typename PerfusionImageType::Pointer image, typename ImageType::Pointer mask);

  /**

  \brief Postprocessing of a recurrence map
  */
  template<class ImageType>
  VectorDouble RecurrenceMapPostprocessing(VectorDouble result, std::vector<typename ImageType::IndexType> indices, typename ImageType::Pointer RecurrenceProbabilityMap, typename ImageType::Pointer edemaMap);



  /**
  \brief Recurrence Estimation using existing model on given subject
  */
  template<class ImageType>
  void PseudoProgressionEstimateOnGivenSubject(typename ImageType::Pointer edemaMask,
    typename ImageType::Pointer tumorMask,
    typename ImageType::Pointer FinalT1CEImagePointer,
    typename ImageType::Pointer FinalT2FlairImagePointer,
    typename ImageType::Pointer FinalT1ImagePointer,
    typename ImageType::Pointer FinalT2ImagePointer,
    std::vector<typename ImageType::Pointer> FinalPerfusionImagePointer,
    std::vector<typename ImageType::Pointer> FinalDTIImagePointer,
    int imagetype, bool conDataPresent, bool dtiDataPresent, bool perfusionDataPresent, bool distanceDataPresent,
    bool useOtherModalities, std::string t1cebasefilename,
    const VectorVectorDouble &nearRegionIndices,
    const VectorVectorDouble &farRegionIndices);

  template<class ImageType>
  void PseudoProgressionEstimateOnGivenSubjectUsingExistingModel(typename ImageType::Pointer edemaMask,
    typename ImageType::Pointer tumorMask,
    typename ImageType::Pointer FinalT1CEImagePointer,
    typename ImageType::Pointer FinalT2FlairImagePointer,
    typename ImageType::Pointer FinalT1ImagePointer,
    typename ImageType::Pointer FinalT2ImagePointer,
    std::vector<typename ImageType::Pointer> FinalPerfusionImagePointer,
    std::vector<typename ImageType::Pointer> FinalDTIImagePointer,
    int imagetype, bool conDataPresent, bool dtiDataPresent, bool perfusionDataPresent, bool distanceDataPresent,
    bool useOtherModalities, std::string t1cebasefilename,
    const VectorVectorDouble &nearRegionIndices,
    const VectorVectorDouble &farRegionIndices, std::string modeldirectory);

  NiftiDataManager mNiftiLocalPtr;
  OutputWritingManager mOutputLocalPtr;
  FeatureReductionClass mFeatureReductionLocalPtr;
  FeatureScalingClass mFeatureScalingLocalPtr;
  FeatureExtractionClass mFeatureExtractionLocalPtr;
  //SVMClassificationClass mClassificationLocalPtr;

  void Run(const std::string &modeldirectory, const std::string &inputdirectory, const std::string &outputdirectory, bool useConvData, bool useDTIData, bool usePerfData, bool useDistData)
  {
    //		    RecurrenceEstimateOnExistingModel(modeldirectory, inputdirectory, outputdirectory, useT1Data, useT1CEData, useT2Data, useT2FlairData, useDTIData, usePerfData, useDistData);
  }

  template<class ImageType>
  void Run(typename ImageType::Pointer edemaMask,
    typename ImageType::Pointer tumorMask,
    typename ImageType::Pointer FinalT1CEImagePointer,
    typename ImageType::Pointer FinalT2FlairImagePointer,
    typename ImageType::Pointer FinalT1ImagePointer,
    typename ImageType::Pointer FinalT2ImagePointer,
    std::vector<typename ImageType::Pointer> FinalPerfusionImagePointer,
    std::vector<typename ImageType::Pointer> FinalDTIImagePointer,
    int imagetype, bool conDataPresent, bool dtiDataPresent, bool perfusionDataPresent, bool distanceDataPresent,
    bool useOtherModalities, const std::string &t1cebasefilename,
    const VectorVectorDouble &nearRegionIndices,
    const VectorVectorDouble &farRegionIndices, const std::string &outputDirName)
  {
    mCurrentOutputDir = outputDirName;
    PseudoProgressionEstimateOnGivenSubject<ImageType>(edemaMask, tumorMask, FinalT1CEImagePointer, FinalT2FlairImagePointer, FinalT1ImagePointer, FinalT2ImagePointer,
      FinalPerfusionImagePointer, FinalDTIImagePointer,
      imagetype, conDataPresent, dtiDataPresent, perfusionDataPresent, distanceDataPresent,
      useOtherModalities, t1cebasefilename, nearRegionIndices, farRegionIndices);
  }

  template <class TImageType>
  typename TImageType::Pointer GetImageWithLabels(std::vector<double> labels, typename TImageType::Pointer inputimage)
  {
    typename TImageType::Pointer output = TImageType::New();
    output->CopyInformation(inputimage);
    output->SetRequestedRegion(inputimage->GetLargestPossibleRegion());
    output->SetBufferedRegion(inputimage->GetBufferedRegion());
    output->Allocate();
    output->FillBuffer(0);

    typedef itk::ImageRegionIteratorWithIndex< TImageType > IteratorType;
    IteratorType maskIt(inputimage, inputimage->GetLargestPossibleRegion());
    IteratorType outputIt(output, output->GetLargestPossibleRegion());
    maskIt.GoToBegin();
    outputIt.GoToBegin();
    while (!maskIt.IsAtEnd())
    {
      for (int j = 0; j < labels.size(); j++)
      {
        if (maskIt.Get() == labels[j])
        {
          outputIt.Set(CAPTK::VOXEL_STATUS::ON);
        }
      }
      ++maskIt;
      ++outputIt;
    }
    return output;
  }


  template<class ImageType>
  void RunLoadedSubjectOnExistingModel(typename ImageType::Pointer edemaMask,
    typename ImageType::Pointer tumorMask,
    typename ImageType::Pointer FinalT1CEImagePointer,
    typename ImageType::Pointer FinalT2FlairImagePointer,
    typename ImageType::Pointer FinalT1ImagePointer,
    typename ImageType::Pointer FinalT2ImagePointer,
    std::vector<typename ImageType::Pointer> FinalPerfusionImagePointer,
    std::vector<typename ImageType::Pointer> FinalDTIImagePointer,
    int imagetype, bool conDataPresent, bool dtiDataPresent, bool perfusionDataPresent, bool distanceDataPresent,
    bool useOtherModalities, const std::string &t1cebasefilename,
    const VectorVectorDouble &nearRegionIndices,
    const VectorVectorDouble &farRegionIndices, const std::string &outputDirName, const std::string &modeldirectory)
  {
    mCurrentOutputDir = outputDirName;
    PseudoProgressionEstimateOnGivenSubjectUsingExistingModel<ImageType>(edemaMask, tumorMask, FinalT1CEImagePointer, FinalT2FlairImagePointer, FinalT1ImagePointer, FinalT2ImagePointer,
      FinalPerfusionImagePointer, FinalDTIImagePointer,
      imagetype, conDataPresent, dtiDataPresent, perfusionDataPresent, distanceDataPresent,
      useOtherModalities, t1cebasefilename, nearRegionIndices, farRegionIndices, modeldirectory);
  }


  std::string mRecurrenceMapFileName;
private:
  std::string mCurrentOutputDir;

  //feature header
  std::vector<std::string> FeatureLabels = { "Eccentricity","Elongation","Perimeter","Roundedness","Flatness",
"T1_Bins_1","T1_Bins_2","T1_Bins_3","T1_Bins_4","T1_Bins_5","T1_Bins_6","T1_Bins_7","T1_Bins_8","T1_Bins_9",
"T1_Bins_10","T1_Intensity_Min","T1_Intensity_Max","T1_Intensity_Mean","T1_Intensity_Variance","T1_Intensity_Std",
"T1_Intensity_Skew","T1_Intensity_Kurtosis","T1_GLCM_Correlation","T1_GLCM_Contrast","T1_GLCM_Entropy",
"T1_GLCM_Homogeneity","T1_GLCM_ClusterShade","T1_GLCM_ClusterProminence","T1_GLCM_AutoCorrelation",
"T1_GLCM_Energy","T1_GLRLM_ShortRunEmphasis","T1_GLRLM_LongRunEmphasis","T1_GLRLM_GLNonUniformity",
"T1_GLRLM_RLNonUniformity","T1_GLRLM_LowGreyLevelRunEmphasis","T1_GLRLM_HighGreyLevelRunEmphasis",
"T1_GLRLM_ShortRunLowGreyLevelEmphasis","T1_GLRLM_ShortRunHighGreyLevelEmphasis","T1_GLRLM_LongRunLowGreyLevelEmphasis",
"T1_GLRLM_LongRunHighGreyLevelEmphasis","T1CE_Bins_1","T1CE_Bins_2","T1CE_Bins_3","T1CE_Bins_4","T1CE_Bins_5",
"T1CE_Bins_6","T1CE_Bins_7","T1CE_Bins_8","T1CE_Bins_9","T1CE_Bins_10","T1CE_Intensity_Min","T1CE_Intensity_Max",
"T1CE_Intensity_Mean","T1CE_Intensity_Variance","T1CE_Intensity_Std","T1CE_Intensity_Skew","T1CE_Intensity_Kurtosis",
"T1CE_GLCM_Correlation","T1CE_GLCM_Contrast","T1CE_GLCM_Entropy","T1CE_GLCM_Homogeneity","T1CE_GLCM_ClusterShade",
"T1CE_GLCM_ClusterProminence","T1CE_GLCM_AutoCorrelation","T1CE_GLCM_Energy","T1CE_GLRLM_ShortRunEmphasis",
"T1CE_GLRLM_LongRunEmphasis","T1CE_GLRLM_GLNonUniformity","T1CE_GLRLM_RLNonUniformity","T1CE_GLRLM_LowGreyLevelRunEmphasis",
"T1CE_GLRLM_HighGreyLevelRunEmphasis","T1CE_GLRLM_ShortRunLowGreyLevelEmphasis","T1CE_GLRLM_ShortRunHighGreyLevelEmphasis",
"T1CE_GLRLM_LongRunLowGreyLevelEmphasis","T1CE_GLRLM_LongRunHighGreyLevelEmphasis","T2_Bins_1","T2_Bins_2","T2_Bins_3",
"T2_Bins_4","T2_Bins_5","T2_Bins_6","T2_Bins_7","T2_Bins_8","T2_Bins_9","T2_Bins_10","T2_Intensity_Min","T2_Intensity_Max",
"T2_Intensity_Mean","T2_Intensity_Variance","T2_Intensity_Std","T2_Intensity_Skew","T2_Intensity_Kurtosis","T2_GLCM_Correlation",
"T2_GLCM_Contrast","T2_GLCM_Entropy","T2_GLCM_Homogeneity","T2_GLCM_ClusterShade","T2_GLCM_ClusterProminence","T2_GLCM_AutoCorrelation",
"T2_GLCM_Energy","T2_GLRLM_ShortRunEmphasis","T2_GLRLM_LongRunEmphasis","T2_GLRLM_GLNonUniformity","T2_GLRLM_RLNonUniformity",
"T2_GLRLM_LowGreyLevelRunEmphasis","T2_GLRLM_HighGreyLevelRunEmphasis","T2_GLRLM_ShortRunLowGreyLevelEmphasis",
"T2_GLRLM_ShortRunHighGreyLevelEmphasis","T2_GLRLM_LongRunLowGreyLevelEmphasis","T2_GLRLM_LongRunHighGreyLevelEmphasis",
"FL_Bins_1","FL_Bins_2","FL_Bins_3","FL_Bins_4","FL_Bins_5","FL_Bins_6","FL_Bins_7","FL_Bins_8","FL_Bins_9","FL_Bins_10",
"FL_Intensity_Min","FL_Intensity_Max","FL_Intensity_Mean","FL_Intensity_Variance","FL_Intensity_Std","FL_Intensity_Skew",
"FL_Intensity_Kurtosis","FL_GLCM_Correlation","FL_GLCM_Contrast","FL_GLCM_Entropy","FL_GLCM_Homogeneity","FL_GLCM_ClusterShade",
"FL_GLCM_ClusterProminence","FL_GLCM_AutoCorrelation","FL_GLCM_Energy","FL_GLRLM_ShortRunEmphasis","FL_GLRLM_LongRunEmphasis",
"FL_GLRLM_GLNonUniformity","FL_GLRLM_RLNonUniformity","FL_GLRLM_LowGreyLevelRunEmphasis","FL_GLRLM_HighGreyLevelRunEmphasis",
"FL_GLRLM_ShortRunLowGreyLevelEmphasis","FL_GLRLM_ShortRunHighGreyLevelEmphasis","FL_GLRLM_LongRunLowGreyLevelEmphasis",
"FL_GLRLM_LongRunHighGreyLevelEmphasis","T1TC_Bins_1","T1TC_Bins_2","T1TC_Bins_3","T1TC_Bins_4","T1TC_Bins_5","T1TC_Bins_6",
"T1TC_Bins_7","T1TC_Bins_8","T1TC_Bins_9","T1TC_Bins_10","T1TC_Intensity_Min","T1TC_Intensity_Max","T1TC_Intensity_Mean",
"T1TC_Intensity_Variance","T1TC_Intensity_Std","T1TC_Intensity_Skew","T1TC_Intensity_Kurtosis","T1TC_GLCM_Correlation",
"T1TC_GLCM_Contrast","T1TC_GLCM_Entropy","T1TC_GLCM_Homogeneity","T1TC_GLCM_ClusterShade","T1TC_GLCM_ClusterProminence",
"T1TC_GLCM_AutoCorrelation","T1TC_GLCM_Energy","T1TC_GLRLM_ShortRunEmphasis","T1TC_GLRLM_LongRunEmphasis","T1TC_GLRLM_GLNonUniformity",
"T1TC_GLRLM_RLNonUniformity","T1TC_GLRLM_LowGreyLevelRunEmphasis","T1TC_GLRLM_HighGreyLevelRunEmphasis",
"T1TC_GLRLM_ShortRunLowGreyLevelEmphasis","T1TC_GLRLM_ShortRunHighGreyLevelEmphasis","T1TC_GLRLM_LongRunLowGreyLevelEmphasis",
"T1TC_GLRLM_LongRunHighGreyLevelEmphasis","T2FL_Bins_1","T2FL_Bins_2","T2FL_Bins_3","T2FL_Bins_4","T2FL_Bins_5","T2FL_Bins_6",
"T2FL_Bins_7","T2FL_Bins_8","T2FL_Bins_9","T2FL_Bins_10","T2FL_Intensity_Min","T2FL_Intensity_Max","T2FL_Intensity_Mean",
"T2FL_Intensity_Variance","T2FL_Intensity_Std","T2FL_Intensity_Skew","T2FL_Intensity_Kurtosis","T2FL_GLCM_Correlation",
"T2FL_GLCM_Contrast","T2FL_GLCM_Entropy","T2FL_GLCM_Homogeneity","T2FL_GLCM_ClusterShade","T2FL_GLCM_ClusterProminence",
"T2FL_GLCM_AutoCorrelation","T2FL_GLCM_Energy","T2FL_GLRLM_ShortRunEmphasis","T2FL_GLRLM_LongRunEmphasis","T2FL_GLRLM_GLNonUniformity",
"T2FL_GLRLM_RLNonUniformity","T2FL_GLRLM_LowGreyLevelRunEmphasis","T2FL_GLRLM_HighGreyLevelRunEmphasis","T2FL_GLRLM_ShortRunLowGreyLevelEmphasis",
"T2FL_GLRLM_ShortRunHighGreyLevelEmphasis","T2FL_GLRLM_LongRunLowGreyLevelEmphasis","T2FL_GLRLM_LongRunHighGreyLevelEmphasis",
"AX_Bins_1","AX_Bins_2","AX_Bins_3","AX_Bins_4","AX_Bins_5","AX_Bins_6","AX_Bins_7","AX_Bins_8","AX_Bins_9","AX_Bins_10",
"AX_Intensity_Min","AX_Intensity_Max","AX_Intensity_Mean","AX_Intensity_Variance","AX_Intensity_Std","AX_Intensity_Skew",
"AX_Intensity_Kurtosis","AX_GLCM_Correlation","AX_GLCM_Contrast","AX_GLCM_Entropy","AX_GLCM_Homogeneity","AX_GLCM_ClusterShade",
"AX_GLCM_ClusterProminence","AX_GLCM_AutoCorrelation","AX_GLCM_Energy","AX_GLRLM_ShortRunEmphasis","AX_GLRLM_LongRunEmphasis",
"AX_GLRLM_GLNonUniformity","AX_GLRLM_RLNonUniformity","AX_GLRLM_LowGreyLevelRunEmphasis","AX_GLRLM_HighGreyLevelRunEmphasis",
"AX_GLRLM_ShortRunLowGreyLevelEmphasis","AX_GLRLM_ShortRunHighGreyLevelEmphasis","AX_GLRLM_LongRunLowGreyLevelEmphasis",
"AX_GLRLM_LongRunHighGreyLevelEmphasis","FA_Bins_1","FA_Bins_2","FA_Bins_3","FA_Bins_4","FA_Bins_5","FA_Bins_6","FA_Bins_7",
"FA_Bins_8","FA_Bins_9","FA_Bins_10","FA_Intensity_Min","FA_Intensity_Max","FA_Intensity_Mean","FA_Intensity_Variance",
"FA_Intensity_Std","FA_Intensity_Skew","FA_Intensity_Kurtosis","FA_GLCM_Correlation","FA_GLCM_Contrast","FA_GLCM_Entropy",
"FA_GLCM_Homogeneity","FA_GLCM_ClusterShade","FA_GLCM_ClusterProminence","FA_GLCM_AutoCorrelation","FA_GLCM_Energy",
"FA_GLRLM_ShortRunEmphasis","FA_GLRLM_LongRunEmphasis","FA_GLRLM_GLNonUniformity","FA_GLRLM_RLNonUniformity",
"FA_GLRLM_LowGreyLevelRunEmphasis","FA_GLRLM_HighGreyLevelRunEmphasis","FA_GLRLM_ShortRunLowGreyLevelEmphasis",
"FA_GLRLM_ShortRunHighGreyLevelEmphasis","FA_GLRLM_LongRunLowGreyLevelEmphasis","FA_GLRLM_LongRunHighGreyLevelEmphasis",
"RAD_Bins_1","RAD_Bins_2","RAD_Bins_3","RAD_Bins_4","RAD_Bins_5","RAD_Bins_6","RAD_Bins_7","RAD_Bins_8","RAD_Bins_9",
"RAD_Bins_10","RAD_Intensity_Min","RAD_Intensity_Max","RAD_Intensity_Mean","RAD_Intensity_Variance","RAD_Intensity_Std",
"RAD_Intensity_Skew","RAD_Intensity_Kurtosis","RAD_GLCM_Correlation","RAD_GLCM_Contrast","RAD_GLCM_Entropy","RAD_GLCM_Homogeneity",
"RAD_GLCM_ClusterShade","RAD_GLCM_ClusterProminence","RAD_GLCM_AutoCorrelation","RAD_GLCM_Energy","RAD_GLRLM_ShortRunEmphasis",
"RAD_GLRLM_LongRunEmphasis","RAD_GLRLM_GLNonUniformity","RAD_GLRLM_RLNonUniformity","RAD_GLRLM_LowGreyLevelRunEmphasis",
"RAD_GLRLM_HighGreyLevelRunEmphasis","RAD_GLRLM_ShortRunLowGreyLevelEmphasis","RAD_GLRLM_ShortRunHighGreyLevelEmphasis",
"RAD_GLRLM_LongRunLowGreyLevelEmphasis","RAD_GLRLM_LongRunHighGreyLevelEmphasis","TR_Bins_1","TR_Bins_2","TR_Bins_3",
"TR_Bins_4","TR_Bins_5","TR_Bins_6","TR_Bins_7","TR_Bins_8","TR_Bins_9","TR_Bins_10","TR_Intensity_Min","TR_Intensity_Max",
"TR_Intensity_Mean","TR_Intensity_Variance","TR_Intensity_Std","TR_Intensity_Skew","TR_Intensity_Kurtosis","TR_GLCM_Correlation",
"TR_GLCM_Contrast","TR_GLCM_Entropy","TR_GLCM_Homogeneity","TR_GLCM_ClusterShade","TR_GLCM_ClusterProminence","TR_GLCM_AutoCorrelation",
"TR_GLCM_Energy","TR_GLRLM_ShortRunEmphasis","TR_GLRLM_LongRunEmphasis","TR_GLRLM_GLNonUniformity","TR_GLRLM_RLNonUniformity",
"TR_GLRLM_LowGreyLevelRunEmphasis","TR_GLRLM_HighGreyLevelRunEmphasis","TR_GLRLM_ShortRunLowGreyLevelEmphasis",
"TR_GLRLM_ShortRunHighGreyLevelEmphasis","TR_GLRLM_LongRunLowGreyLevelEmphasis","TR_GLRLM_LongRunHighGreyLevelEmphasis",
"PH_Bins_1","PH_Bins_2","PH_Bins_3","PH_Bins_4","PH_Bins_5","PH_Bins_6","PH_Bins_7","PH_Bins_8","PH_Bins_9","PH_Bins_10",
"PH_Intensity_Min","PH_Intensity_Max","PH_Intensity_Mean","PH_Intensity_Variance","PH_Intensity_Std","PH_Intensity_Skew",
"PH_Intensity_Kurtosis","PH_GLCM_Correlation","PH_GLCM_Contrast","PH_GLCM_Entropy","PH_GLCM_Homogeneity","PH_GLCM_ClusterShade",
"PH_GLCM_ClusterProminence","PH_GLCM_AutoCorrelation","PH_GLCM_Energy","PH_GLRLM_ShortRunEmphasis","PH_GLRLM_LongRunEmphasis",
"PH_GLRLM_GLNonUniformity","PH_GLRLM_RLNonUniformity","PH_GLRLM_LowGreyLevelRunEmphasis","PH_GLRLM_HighGreyLevelRunEmphasis",
"PH_GLRLM_ShortRunLowGreyLevelEmphasis","PH_GLRLM_ShortRunHighGreyLevelEmphasis","PH_GLRLM_LongRunLowGreyLevelEmphasis",
"PH_GLRLM_LongRunHighGreyLevelEmphasis","PS_Bins_1","PS_Bins_2","PS_Bins_3","PS_Bins_4","PS_Bins_5","PS_Bins_6","PS_Bins_7",
"PS_Bins_8","PS_Bins_9","PS_Bins_10","PS_Intensity_Min","PS_Intensity_Max","PS_Intensity_Mean","PS_Intensity_Variance",
"PS_Intensity_Std","PS_Intensity_Skew","PS_Intensity_Kurtosis","PS_GLCM_Correlation","PS_GLCM_Contrast","PS_GLCM_Entropy",
"PS_GLCM_Homogeneity","PS_GLCM_ClusterShade","PS_GLCM_ClusterProminence","PS_GLCM_AutoCorrelation","PS_GLCM_Energy",
"PS_GLRLM_ShortRunEmphasis","PS_GLRLM_LongRunEmphasis","PS_GLRLM_GLNonUniformity","PS_GLRLM_RLNonUniformity",
"PS_GLRLM_LowGreyLevelRunEmphasis","PS_GLRLM_HighGreyLevelRunEmphasis","PS_GLRLM_ShortRunLowGreyLevelEmphasis",
"PS_GLRLM_ShortRunHighGreyLevelEmphasis","PS_GLRLM_LongRunLowGreyLevelEmphasis","PS_GLRLM_LongRunHighGreyLevelEmphasis",
"RCBV_Bins_1","RCBV_Bins_2","RCBV_Bins_3","RCBV_Bins_4","RCBV_Bins_5","RCBV_Bins_6","RCBV_Bins_7","RCBV_Bins_8","RCBV_Bins_9",
"RCBV_Bins_10","RCBV_Intensity_Min","RCBV_Intensity_Max","RCBV_Intensity_Mean","RCBV_Intensity_Variance","RCBV_Intensity_Std",
"RCBV_Intensity_Skew","RCBV_Intensity_Kurtosis","RCBV_GLCM_Correlation","RCBV_GLCM_Contrast","RCBV_GLCM_Entropy",
"RCBV_GLCM_Homogeneity","RCBV_GLCM_ClusterShade","RCBV_GLCM_ClusterProminence","RCBV_GLCM_AutoCorrelation","RCBV_GLCM_Energy",
"RCBV_GLRLM_ShortRunEmphasis","RCBV_GLRLM_LongRunEmphasis","RCBV_GLRLM_GLNonUniformity","RCBV_GLRLM_RLNonUniformity",
"RCBV_GLRLM_LowGreyLevelRunEmphasis","RCBV_GLRLM_HighGreyLevelRunEmphasis","RCBV_GLRLM_ShortRunLowGreyLevelEmphasis",
"RCBV_GLRLM_ShortRunHighGreyLevelEmphasis","RCBV_GLRLM_LongRunLowGreyLevelEmphasis","RCBV_GLRLM_LongRunHighGreyLevelEmphasis",
"PCA1_Bins_1","PCA1_Bins_2","PCA1_Bins_3","PCA1_Bins_4","PCA1_Bins_5","PCA1_Bins_6","PCA1_Bins_7","PCA1_Bins_8","PCA1_Bins_9",
"PCA1_Bins_10","PCA1_Intensity_Min","PCA1_Intensity_Max","PCA1_Intensity_Mean","PCA1_Intensity_Variance","PCA1_Intensity_Std",
"PCA1_Intensity_Skew","PCA1_Intensity_Kurtosis","PCA1_GLCM_Correlation","PCA1_GLCM_Contrast","PCA1_GLCM_Entropy",
"PCA1_GLCM_Homogeneity","PCA1_GLCM_ClusterShade","PCA1_GLCM_ClusterProminence","PCA1_GLCM_AutoCorrelation","PCA1_GLCM_Energy",
"PCA1_GLRLM_ShortRunEmphasis","PCA1_GLRLM_LongRunEmphasis","PCA1_GLRLM_GLNonUniformity","PCA1_GLRLM_RLNonUniformity",
"PCA1_GLRLM_LowGreyLevelRunEmphasis","PCA1_GLRLM_HighGreyLevelRunEmphasis","PCA1_GLRLM_ShortRunLowGreyLevelEmphasis",
"PCA1_GLRLM_ShortRunHighGreyLevelEmphasis","PCA1_GLRLM_LongRunLowGreyLevelEmphasis","PCA1_GLRLM_LongRunHighGreyLevelEmphasis",
"PCA2_Bins_1","PCA2_Bins_2","PCA2_Bins_3","PCA2_Bins_4","PCA2_Bins_5","PCA2_Bins_6","PCA2_Bins_7","PCA2_Bins_8","PCA2_Bins_9",
"PCA2_Bins_10","PCA2_Intensity_Min","PCA2_Intensity_Max","PCA2_Intensity_Mean","PCA2_Intensity_Variance","PCA2_Intensity_Std",
"PCA2_Intensity_Skew","PCA2_Intensity_Kurtosis","PCA2_GLCM_Correlation","PCA2_GLCM_Contrast","PCA2_GLCM_Entropy",
"PCA2_GLCM_Homogeneity","PCA2_GLCM_ClusterShade","PCA2_GLCM_ClusterProminence","PCA2_GLCM_AutoCorrelation","PCA2_GLCM_Energy",
"PCA2_GLRLM_ShortRunEmphasis","PCA2_GLRLM_LongRunEmphasis","PCA2_GLRLM_GLNonUniformity","PCA2_GLRLM_RLNonUniformity",
"PCA2_GLRLM_LowGreyLevelRunEmphasis","PCA2_GLRLM_HighGreyLevelRunEmphasis","PCA2_GLRLM_ShortRunLowGreyLevelEmphasis",
"PCA2_GLRLM_ShortRunHighGreyLevelEmphasis","PCA2_GLRLM_LongRunLowGreyLevelEmphasis","PCA2_GLRLM_LongRunHighGreyLevelEmphasis",
"PCA3_Bins_1","PCA3_Bins_2","PCA3_Bins_3","PCA3_Bins_4","PCA3_Bins_5","PCA3_Bins_6","PCA3_Bins_7","PCA3_Bins_8","PCA3_Bins_9",
"PCA3_Bins_10","PCA3_Intensity_Min","PCA3_Intensity_Max","PCA3_Intensity_Mean","PCA3_Intensity_Variance","PCA3_Intensity_Std",
"PCA3_Intensity_Skew","PCA3_Intensity_Kurtosis","PCA3_GLCM_Correlation","PCA3_GLCM_Contrast","PCA3_GLCM_Entropy",
"PCA3_GLCM_Homogeneity","PCA3_GLCM_ClusterShade","PCA3_GLCM_ClusterProminence","PCA3_GLCM_AutoCorrelation","PCA3_GLCM_Energy",
"PCA3_GLRLM_ShortRunEmphasis","PCA3_GLRLM_LongRunEmphasis","PCA3_GLRLM_GLNonUniformity","PCA3_GLRLM_RLNonUniformity",
"PCA3_GLRLM_LowGreyLevelRunEmphasis","PCA3_GLRLM_HighGreyLevelRunEmphasis","PCA3_GLRLM_ShortRunLowGreyLevelEmphasis",
"PCA3_GLRLM_ShortRunHighGreyLevelEmphasis","PCA3_GLRLM_LongRunLowGreyLevelEmphasis","PCA3_GLRLM_LongRunHighGreyLevelEmphasis",
"PCA4_Bins_1","PCA4_Bins_2","PCA4_Bins_3","PCA4_Bins_4","PCA4_Bins_5","PCA4_Bins_6","PCA4_Bins_7","PCA4_Bins_8","PCA4_Bins_9",
"PCA4_Bins_10","PCA4_Intensity_Min","PCA4_Intensity_Max","PCA4_Intensity_Mean","PCA4_Intensity_Variance","PCA4_Intensity_Std",
"PCA4_Intensity_Skew","PCA4_Intensity_Kurtosis","PCA4_GLCM_Correlation","PCA4_GLCM_Contrast","PCA4_GLCM_Entropy","PCA4_GLCM_Homogeneity",
"PCA4_GLCM_ClusterShade","PCA4_GLCM_ClusterProminence","PCA4_GLCM_AutoCorrelation","PCA4_GLCM_Energy","PCA4_GLRLM_ShortRunEmphasis",
"PCA4_GLRLM_LongRunEmphasis","PCA4_GLRLM_GLNonUniformity","PCA4_GLRLM_RLNonUniformity","PCA4_GLRLM_LowGreyLevelRunEmphasis",
"PCA4_GLRLM_HighGreyLevelRunEmphasis","PCA4_GLRLM_ShortRunLowGreyLevelEmphasis","PCA4_GLRLM_ShortRunHighGreyLevelEmphasis",
"PCA4_GLRLM_LongRunLowGreyLevelEmphasis","PCA4_GLRLM_LongRunHighGreyLevelEmphasis","PCA5_Bins_1","PCA5_Bins_2","PCA5_Bins_3",
"PCA5_Bins_4","PCA5_Bins_5","PCA5_Bins_6","PCA5_Bins_7","PCA5_Bins_8","PCA5_Bins_9","PCA5_Bins_10","PCA5_Intensity_Min",
"PCA5_Intensity_Max","PCA5_Intensity_Mean","PCA5_Intensity_Variance","PCA5_Intensity_Std","PCA5_Intensity_Skew",
"PCA5_Intensity_Kurtosis","PCA5_GLCM_Correlation","PCA5_GLCM_Contrast","PCA5_GLCM_Entropy","PCA5_GLCM_Homogeneity",
"PCA5_GLCM_ClusterShade","PCA5_GLCM_ClusterProminence","PCA5_GLCM_AutoCorrelation","PCA5_GLCM_Energy","PCA5_GLRLM_ShortRunEmphasis",
"PCA5_GLRLM_LongRunEmphasis","PCA5_GLRLM_GLNonUniformity","PCA5_GLRLM_RLNonUniformity","PCA5_GLRLM_LowGreyLevelRunEmphasis",
"PCA5_GLRLM_HighGreyLevelRunEmphasis","PCA5_GLRLM_ShortRunLowGreyLevelEmphasis","PCA5_GLRLM_ShortRunHighGreyLevelEmphasis",
"PCA5_GLRLM_LongRunLowGreyLevelEmphasis","PCA5_GLRLM_LongRunHighGreyLevelEmphasis","PCA6_Bins_1","PCA6_Bins_2","PCA6_Bins_3",
"PCA6_Bins_4","PCA6_Bins_5","PCA6_Bins_6","PCA6_Bins_7","PCA6_Bins_8","PCA6_Bins_9","PCA6_Bins_10","PCA6_Intensity_Min",
"PCA6_Intensity_Max","PCA6_Intensity_Mean","PCA6_Intensity_Variance","PCA6_Intensity_Std","PCA6_Intensity_Skew",
"PCA6_Intensity_Kurtosis","PCA6_GLCM_Correlation","PCA6_GLCM_Contrast","PCA6_GLCM_Entropy","PCA6_GLCM_Homogeneity",
"PCA6_GLCM_ClusterShade","PCA6_GLCM_ClusterProminence","PCA6_GLCM_AutoCorrelation","PCA6_GLCM_Energy","PCA6_GLRLM_ShortRunEmphasis",
"PCA6_GLRLM_LongRunEmphasis","PCA6_GLRLM_GLNonUniformity","PCA6_GLRLM_RLNonUniformity","PCA6_GLRLM_LowGreyLevelRunEmphasis",
"PCA6_GLRLM_HighGreyLevelRunEmphasis","PCA6_GLRLM_ShortRunLowGreyLevelEmphasis","PCA6_GLRLM_ShortRunHighGreyLevelEmphasis",
"PCA6_GLRLM_LongRunLowGreyLevelEmphasis","PCA6_GLRLM_LongRunHighGreyLevelEmphasis","PCA7_Bins_1","PCA7_Bins_2","PCA7_Bins_3",
"PCA7_Bins_4","PCA7_Bins_5","PCA7_Bins_6","PCA7_Bins_7","PCA7_Bins_8","PCA7_Bins_9","PCA7_Bins_10","PCA7_Intensity_Min",
"PCA7_Intensity_Max","PCA7_Intensity_Mean","PCA7_Intensity_Variance","PCA7_Intensity_Std","PCA7_Intensity_Skew","PCA7_Intensity_Kurtosis",
"PCA7_GLCM_Correlation","PCA7_GLCM_Contrast","PCA7_GLCM_Entropy","PCA7_GLCM_Homogeneity","PCA7_GLCM_ClusterShade",
"PCA7_GLCM_ClusterProminence","PCA7_GLCM_AutoCorrelation","PCA7_GLCM_Energy","PCA7_GLRLM_ShortRunEmphasis","PCA7_GLRLM_LongRunEmphasis",
"PCA7_GLRLM_GLNonUniformity","PCA7_GLRLM_RLNonUniformity","PCA7_GLRLM_LowGreyLevelRunEmphasis","PCA7_GLRLM_HighGreyLevelRunEmphasis",
"PCA7_GLRLM_ShortRunLowGreyLevelEmphasis","PCA7_GLRLM_ShortRunHighGreyLevelEmphasis","PCA7_GLRLM_LongRunLowGreyLevelEmphasis",
"PCA7_GLRLM_LongRunHighGreyLevelEmphasis","PCA8_Bins_1","PCA8_Bins_2","PCA8_Bins_3","PCA8_Bins_4","PCA8_Bins_5","PCA8_Bins_6",
"PCA8_Bins_7","PCA8_Bins_8","PCA8_Bins_9","PCA8_Bins_10","PCA8_Intensity_Min","PCA8_Intensity_Max","PCA8_Intensity_Mean",
"PCA8_Intensity_Variance","PCA8_Intensity_Std","PCA8_Intensity_Skew","PCA8_Intensity_Kurtosis","PCA8_GLCM_Correlation",
"PCA8_GLCM_Contrast","PCA8_GLCM_Entropy","PCA8_GLCM_Homogeneity","PCA8_GLCM_ClusterShade","PCA8_GLCM_ClusterProminence",
"PCA8_GLCM_AutoCorrelation","PCA8_GLCM_Energy","PCA8_GLRLM_ShortRunEmphasis","PCA8_GLRLM_LongRunEmphasis","PCA8_GLRLM_GLNonUniformity",
"PCA8_GLRLM_RLNonUniformity","PCA8_GLRLM_LowGreyLevelRunEmphasis","PCA8_GLRLM_HighGreyLevelRunEmphasis",
"PCA8_GLRLM_ShortRunLowGreyLevelEmphasis","PCA8_GLRLM_ShortRunHighGreyLevelEmphasis","PCA8_GLRLM_LongRunLowGreyLevelEmphasis",
"PCA8_GLRLM_LongRunHighGreyLevelEmphasis","PCA9_Bins_1","PCA9_Bins_2","PCA9_Bins_3","PCA9_Bins_4","PCA9_Bins_5","PCA9_Bins_6",
"PCA9_Bins_7","PCA9_Bins_8","PCA9_Bins_9","PCA9_Bins_10","PCA9_Intensity_Min","PCA9_Intensity_Max","PCA9_Intensity_Mean",
"PCA9_Intensity_Variance","PCA9_Intensity_Std","PCA9_Intensity_Skew","PCA9_Intensity_Kurtosis","PCA9_GLCM_Correlation",
"PCA9_GLCM_Contrast","PCA9_GLCM_Entropy","PCA9_GLCM_Homogeneity","PCA9_GLCM_ClusterShade","PCA9_GLCM_ClusterProminence",
"PCA9_GLCM_AutoCorrelation","PCA9_GLCM_Energy","PCA9_GLRLM_ShortRunEmphasis","PCA9_GLRLM_LongRunEmphasis","PCA9_GLRLM_GLNonUniformity",
"PCA9_GLRLM_RLNonUniformity","PCA9_GLRLM_LowGreyLevelRunEmphasis","PCA9_GLRLM_HighGreyLevelRunEmphasis",
"PCA9_GLRLM_ShortRunLowGreyLevelEmphasis","PCA9_GLRLM_ShortRunHighGreyLevelEmphasis","PCA9_GLRLM_LongRunLowGreyLevelEmphasis",
"PCA9_GLRLM_LongRunHighGreyLevelEmphasis","PCA10_Bins_1","PCA10_Bins_2","PCA10_Bins_3","PCA10_Bins_4","PCA10_Bins_5",
"PCA10_Bins_6","PCA10_Bins_7","PCA10_Bins_8","PCA10_Bins_9","PCA10_Bins_10","PCA10_Intensity_Min","PCA10_Intensity_Max",
"PCA10_Intensity_Mean","PCA10_Intensity_Variance","PCA10_Intensity_Std","PCA10_Intensity_Skew","PCA10_Intensity_Kurtosis",
"PCA10_GLCM_Correlation","PCA10_GLCM_Contrast","PCA10_GLCM_Entropy","PCA10_GLCM_Homogeneity","PCA10_GLCM_ClusterShade",
"PCA10_GLCM_ClusterProminence","PCA10_GLCM_AutoCorrelation","PCA10_GLCM_Energy","PCA10_GLRLM_ShortRunEmphasis",
"PCA10_GLRLM_LongRunEmphasis","PCA10_GLRLM_GLNonUniformity","PCA10_GLRLM_RLNonUniformity","PCA10_GLRLM_LowGreyLevelRunEmphasis",
"PCA10_GLRLM_HighGreyLevelRunEmphasis","PCA10_GLRLM_ShortRunLowGreyLevelEmphasis","PCA10_GLRLM_ShortRunHighGreyLevelEmphasis",
"PCA10_GLRLM_LongRunLowGreyLevelEmphasis","PCA10_GLRLM_LongRunHighGreyLevelEmphasis","T1_PCA_1","T1_PCA_2","T1_PCA_3",
"T1_PCA_4","T1_PCA_5","T1_PCA_6","T1_PCA_7","T1_PCA_8","T1_PCA_9","T1_PCA_10","T1CE_PCA_1","T1CE_PCA_2","T1CE_PCA_3",
"T1CE_PCA_4","T1CE_PCA_5","T1CE_PCA_6","T1CE_PCA_7","T1CE_PCA_8","T1CE_PCA_9","T1CE_PCA_10","T1T1CE_PCA_1","T1T1CE_PCA_2",
"T1T1CE_PCA_3","T1T1CE_PCA_4","T1T1CE_PCA_5","T1T1CE_PCA_6","T1T1CE_PCA_7","T1T1CE_PCA_8","T1T1CE_PCA_9","T1T1CE_PCA_10",
"T2_PCA_1","T2_PCA_2","T2_PCA_3","T2_PCA_4","T2_PCA_5","T2_PCA_6","T2_PCA_7","T2_PCA_8","T2_PCA_9","T2_PCA_10","FL_PCA_1",
"FL_PCA_2","FL_PCA_3","FL_PCA_4","FL_PCA_5","FL_PCA_6","FL_PCA_7","FL_PCA_8","FL_PCA_9","FL_PCA_10","T2FL_PCA_1","T2FL_PCA_2",
"T2FL_PCA_3","T2FL_PCA_4","T2FL_PCA_5","T2FL_PCA_6","T2FL_PCA_7","T2FL_PCA_8","T2FL_PCA_9","T2FL_PCA_10","AX_PCA_1",
"AX_PCA_2","AX_PCA_3","AX_PCA_4","AX_PCA_5","AX_PCA_6","AX_PCA_7","AX_PCA_8","AX_PCA_9","AX_PCA_10","FA_PCA_1","FA_PCA_2",
"FA_PCA_3","FA_PCA_4","FA_PCA_5","FA_PCA_6","FA_PCA_7","FA_PCA_8","FA_PCA_9","FA_PCA_10","RAD_PCA_1","RAD_PCA_2","RAD_PCA_3",
"RAD_PCA_4","RAD_PCA_5","RAD_PCA_6","RAD_PCA_7","RAD_PCA_8","RAD_PCA_9","RAD_PCA_10","TR_PCA_1","TR_PCA_2","TR_PCA_3",
"TR_PCA_4","TR_PCA_5","TR_PCA_6","TR_PCA_7","TR_PCA_8","TR_PCA_9","TR_PCA_10","PH_PCA_1","PH_PCA_2","PH_PCA_3","PH_PCA_4",
"PH_PCA_5","PH_PCA_6","PH_PCA_7","PH_PCA_8","PH_PCA_9","PH_PCA_10","PSR_PCA_1","PSR_PCA_2","PSR_PCA_3","PSR_PCA_4",
"PSR_PCA_5","PSR_PCA_6","PSR_PCA_7","PSR_PCA_8","PSR_PCA_9","PSR_PCA_10","RCBV_PCA_1","RCBV_PCA_2","RCBV_PCA_3",
"RCBV_PCA_4","RCBV_PCA_5","RCBV_PCA_6","RCBV_PCA_7","RCBV_PCA_8","RCBV_PCA_9","RCBV_PCA_10","PCA1_PCA_1","PCA1_PCA_2",
"PCA1_PCA_3","PCA1_PCA_4","PCA1_PCA_5","PCA1_PCA_6","PCA1_PCA_7","PCA1_PCA_8","PCA1_PCA_9","PCA1_PCA_10","PCA2_PCA_1",
"PCA2_PCA_2","PCA2_PCA_3","PCA2_PCA_4","PCA2_PCA_5","PCA2_PCA_6","PCA2_PCA_7","PCA2_PCA_8","PCA2_PCA_9","PCA2_PCA_10",
"PCA3_PCA_1","PCA3_PCA_2","PCA3_PCA_3","PCA3_PCA_4","PCA3_PCA_5","PCA3_PCA_6","PCA3_PCA_7","PCA3_PCA_8","PCA3_PCA_9",
"PCA3_PCA_10","PCA4_PCA_1","PCA4_PCA_2","PCA4_PCA_3","PCA4_PCA_4","PCA4_PCA_5","PCA4_PCA_6","PCA4_PCA_7","PCA4_PCA_8",
"PCA4_PCA_9","PCA4_PCA_10","PCA5_PCA_1","PCA5_PCA_2","PCA5_PCA_3","PCA5_PCA_4","PCA5_PCA_5","PCA5_PCA_6","PCA5_PCA_7",
"PCA5_PCA_8","PCA5_PCA_9","PCA5_PCA_10","PCA6_PCA_1","PCA6_PCA_2","PCA6_PCA_3","PCA6_PCA_4","PCA6_PCA_5","PCA6_PCA_6",
"PCA6_PCA_7","PCA6_PCA_8","PCA6_PCA_9","PCA6_PCA_10","PCA7_PCA_1","PCA7_PCA_2","PCA7_PCA_3","PCA7_PCA_4","PCA7_PCA_5",
"PCA7_PCA_6","PCA7_PCA_7","PCA7_PCA_8","PCA7_PCA_9","PCA7_PCA_10","PCA8_PCA_1","PCA8_PCA_2","PCA8_PCA_3","PCA8_PCA_4",
"PCA8_PCA_5","PCA8_PCA_6","PCA8_PCA_7","PCA8_PCA_8","PCA8_PCA_9","PCA8_PCA_10","PCA9_PCA_1","PCA9_PCA_2","PCA9_PCA_3",
"PCA9_PCA_4","PCA9_PCA_5","PCA9_PCA_6","PCA9_PCA_7","PCA9_PCA_8","PCA9_PCA_9","PCA9_PCA_10","PCA10_PCA_1","PCA10_PCA_2",
"PCA10_PCA_3","PCA10_PCA_4","PCA10_PCA_5","PCA10_PCA_6","PCA10_PCA_7","PCA10_PCA_8","PCA10_PCA_9","PCA10_PCA_10" };

};

template<class ImageType>
void PseudoProgressionEstimator::PseudoProgressionEstimateOnGivenSubject(typename ImageType::Pointer edemaMask,
  typename ImageType::Pointer tumorMask,
  typename ImageType::Pointer FinalT1CEImagePointer,
  typename ImageType::Pointer FinalT2FlairImagePointer,
  typename ImageType::Pointer FinalT1ImagePointer,
  typename ImageType::Pointer FinalT2ImagePointer,
  std::vector<typename ImageType::Pointer> FinalPerfusionImagePointer,
  std::vector<typename ImageType::Pointer> FinalDTIImagePointer,
  int imagetype, bool conDataPresent, bool dtiDataPresent, bool perfusionDataPresent, bool distanceDataPresent,
  bool useOtherModalities, std::string t1cebasefilename,
  const VectorVectorDouble &nearRegionIndices,
  const VectorVectorDouble &farRegionIndices)
{
  //	//------------------------------------------training data formulation------------------------------------------
  //	mOutputLocalPtr.SetOutputDirectoryPath(mCurrentOutputDir);
  //	VectorVectorDouble mNearIntensities;
  //	VectorVectorDouble mFarIntensities;
  //	VectorVectorDouble pNearIntensities;
  //	VectorVectorDouble pFarIntensities;
  //	VectorDouble tNearIntensities;
  //	VectorDouble tFarIntensities;
  //
  //	VectorVectorDouble perfusionIntensities;
  //	VectorVectorDouble otherIntensities;
  //	VectorDouble distanceIntensities;
  //	VectorVectorDouble fNearIntensities;
  //	VectorVectorDouble fFarIntensities;
  //#ifdef APP_BASE_CAPTK_H
  //	messageUpdate("Recurrence Estimation");
  //	progressUpdate(0);
  //#endif
  //	mNiftiLocalPtr.LoadTrainingData(tumorMask, nearRegionIndices, farRegionIndices, FinalT1CEImagePointer, FinalT2FlairImagePointer, FinalT1ImagePointer, FinalT2ImagePointer, FinalPerfusionImagePointer, FinalDTIImagePointer, pNearIntensities, pFarIntensities, mNearIntensities, mFarIntensities, tNearIntensities, tFarIntensities, NIfTI, conDataPresent, dtiDataPresent, perfusionDataPresent, distanceDataPresent);
  //
  //	for (unsigned int i = 0; i < mNearIntensities.size(); i++)
  //	{
  //		if (perfusionDataPresent)
  //			perfusionIntensities.push_back(pNearIntensities[i]);
  //		if (useOtherModalities)
  //			otherIntensities.push_back(mNearIntensities[i]);
  //		if (distanceDataPresent)
  //			distanceIntensities.push_back(tNearIntensities[i]);
  //	}
  //	for (unsigned int i = 0; i < mFarIntensities.size(); i++)
  //	{
  //		if (perfusionDataPresent)
  //			perfusionIntensities.push_back(pFarIntensities[i]);
  //		if (useOtherModalities)
  //			otherIntensities.push_back(mFarIntensities[i]);
  //		if (distanceDataPresent)
  //			distanceIntensities.push_back(tFarIntensities[i]);
  //	}
  //#ifdef APP_BASE_CAPTK_H
  //	progressUpdate(10);
  //#endif
  //
  //	VariableLengthVectorType perfMeanVector;
  //	vtkSmartPointer<vtkTable> reducedPerfusionFeatures;
  //	//------------------------------------------reduce perfusion intensities------------------------------------------
  //	if (perfusionDataPresent)
  //	{
  //		reducedPerfusionFeatures = mFeatureReductionLocalPtr.GetDiscerningPerfusionTimePoints(perfusionIntensities);
  //	}
  //	//-----------------------------------develope final near and far vectors------------------------------------------
  //	for (unsigned int i = 0; i < mNearIntensities.size(); i++)
  //	{
  //		VectorDouble cIntensityVectorPerSub;
  //		if (perfusionDataPresent)
  //			for (int j = 0; j < NO_OF_PCS; j++)
  //				cIntensityVectorPerSub.push_back(reducedPerfusionFeatures->GetValue(i, j).ToDouble());
  //
  //		if (useOtherModalities)
  //			for (unsigned int j = 0; j < otherIntensities[i].size(); j++)
  //				cIntensityVectorPerSub.push_back(otherIntensities[i][j]);
  //
  //		if (distanceDataPresent)
  //			cIntensityVectorPerSub.push_back(distanceIntensities[i]);
  //
  //		fNearIntensities.push_back(cIntensityVectorPerSub);
  //	}
  //
  //	for (unsigned int i = mNearIntensities.size(); i < mNearIntensities.size() + mFarIntensities.size(); i++)
  //	{
  //		VectorDouble cIntensityVectorPerSub;
  //		if (perfusionDataPresent)
  //			for (int j = 0; j < NO_OF_PCS; j++)
  //				cIntensityVectorPerSub.push_back(reducedPerfusionFeatures->GetValue(i, j).ToDouble());
  //
  //		if (useOtherModalities)
  //			for (unsigned int j = 0; j < otherIntensities[i].size(); j++)
  //				cIntensityVectorPerSub.push_back(otherIntensities[i][j]);
  //
  //		if (distanceDataPresent)
  //			cIntensityVectorPerSub.push_back(distanceIntensities[i]);
  //		fFarIntensities.push_back(cIntensityVectorPerSub);
  //	}
  //	mFeatureExtractionLocalPtr.FormulateTrainingData(fNearIntensities, fFarIntensities);
  //	VariableSizeMatrixType TrainingData = mFeatureExtractionLocalPtr.GetTrainingData();
  //
  //#ifdef APP_BASE_CAPTK_H
  //	progressUpdate(20);
  //#endif
  //
  //	//typedef vnl_matrix<double> MatrixType;
  //	//MatrixType data;
  //	//data.set_size(120, 45);
  //	//for (int i = 0; i < perfusionIntensities.size(); i++)
  //	//	for (int j = 0; j < 45; j++)
  //	//		data(i, j) = perfusionIntensities[i][j];
  //	//typedef itk::CSVNumericObjectFileWriter<double, 120, 45> WriterType;
  //	//WriterType::Pointer writer = WriterType::New();
  //	//writer->SetFileName("perfusionData.csv");
  //	//writer->SetInput(&data);
  //	//writer->Write();
  //
  //
  //	//data.set_size(120,15);
  //	//for (unsigned int i = 0; i < TrainingData.Rows(); i++)
  //	//  for (unsigned int j = 0; j < TrainingData.Cols(); j++)
  //	//		data(i, j) = TrainingData[i][j];
  //	//writer->SetFileName("tData.csv");
  //	//writer->SetInput(&data);
  //	//writer->Write();
  //
  //
  //	VariableSizeMatrixType ScaledTrainingData = mFeatureScalingLocalPtr.ScaleGivenTrainingFeatures(TrainingData);
  //
  //	//data.set_size(120, 15);
  //	//for (unsigned int i = 0; i < ScaledTrainingData.Rows(); i++)
  //	//  for (unsigned int j = 0; j < ScaledTrainingData.Cols(); j++)
  //	//		data(i, j) = ScaledTrainingData[i][j];
  //	//writer->SetFileName("sData.csv");
  //	//writer->SetInput(&data);
  //	//writer->Write();
  //
  //	////------------------------------------------training process---------------------------------------------------
  //	//FILE *t;
  //	//t = fopen("TrainingData.txt", "w");
  //
  //	//for (int i = 0; i < ScaledTrainingData.Rows(); i++)
  //	//{
  //	//	fprintf(t, "%f ", ScaledTrainingData[i][ScaledTrainingData.Cols() - 1]);
  //	//	for (int j = 0; j < ScaledTrainingData.Cols() - 1; j++)
  //	//		fprintf(t, "%d:%lf ", j + 1, ScaledTrainingData[i][j]);
  //	//	fprintf(t, "\n");
  //	//}
  //	//fclose(t);
  //
  //	int size = GetFeatureVectorSize(conDataPresent, dtiDataPresent, perfusionDataPresent, distanceDataPresent);
  //	mOutputLocalPtr.SaveModelResults(ScaledTrainingData, mFeatureScalingLocalPtr.GetMeanVector(), mFeatureScalingLocalPtr.GetStdVector(), mFeatureReductionLocalPtr.GetPerfusionMeanVector(), mFeatureReductionLocalPtr.GetPCATransformationMatrix(),
  //		conDataPresent, dtiDataPresent, perfusionDataPresent, distanceDataPresent, size);
  //#ifdef APP_BASE_CAPTK_H
  //	progressUpdate(30);
  //#endif
  //
  //	VectorDouble trainDistances = trainOpenCVSVM(ScaledTrainingData, mOutputLocalPtr.mOutputDirectoryPath + "/" + mTrainedModelNameXML, true, Recurrence);
  //	//  mClassificationLocalPtr.Training(ScaledTrainingData, mOutputLocalPtr.mOutputDirectoryPath);
  //	//if (mClassificationLocalPtr->GetLastEncounteredError() != "")
  //	//  mLastEncounteredError = mClassificationLocalPtr->GetLastEncounteredError();
  //#ifdef APP_BASE_CAPTK_H
  //	progressUpdate(40);
  //#endif
  //	////------------------------------------------testing data formulation---------------------------------------------------
  //
  //	const typename ImageType::RegionType region = FinalT1CEImagePointer->GetLargestPossibleRegion();
  //	typename ImageType::Pointer RecurrenceProbabilityMap = ImageType::New();
  //	RecurrenceProbabilityMap->SetRegions(region);
  //	RecurrenceProbabilityMap->Allocate();
  //	RecurrenceProbabilityMap->SetSpacing(FinalT1CEImagePointer->GetSpacing());
  //	RecurrenceProbabilityMap->SetOrigin(FinalT1CEImagePointer->GetOrigin());
  //	RecurrenceProbabilityMap->SetDirection(FinalT1CEImagePointer->GetDirection());
  //
  //	//VectorVectorDouble perfusionIntensities;
  //	//VectorVectorDouble otherIntensities;
  //	//VectorDouble distanceIntensities;
  //	std::vector<typename ImageType::IndexType> testindices;
  //	testindices = mNiftiLocalPtr.LoadTestData(FinalT1CEImagePointer, FinalT2FlairImagePointer, FinalT1ImagePointer, FinalT2ImagePointer, FinalPerfusionImagePointer, FinalDTIImagePointer, tumorMask, edemaMask, perfusionIntensities, otherIntensities, distanceIntensities, DICOM, conDataPresent, dtiDataPresent, perfusionDataPresent, distanceDataPresent);
  //
  //
  //
  //	VectorVectorDouble reducedTestPerfusionFeatures;
  //	if (perfusionDataPresent)
  //		reducedTestPerfusionFeatures = mFeatureReductionLocalPtr.ApplyPCAOnTestData(perfusionIntensities);
  //
  //	int NumberOfPCs = 5;
  //	VectorVectorDouble globaltestintensities;
  //
  //	for (unsigned int k = 0; k < testindices.size(); k++)
  //	{
  //		VectorDouble inten;
  //
  //		if (perfusionDataPresent)
  //			for (int j = 0; j < NumberOfPCs; j++)
  //				inten.push_back(reducedTestPerfusionFeatures[k][j]);
  //
  //		if (useOtherModalities)
  //			for (unsigned int j = 0; j < otherIntensities[0].size(); j++)
  //				inten.push_back(otherIntensities[k][j]);
  //
  //		if (distanceDataPresent)
  //			inten.push_back(distanceIntensities[k]);
  //
  //		if (inten.size()>0)
  //			globaltestintensities.push_back(inten);
  //	}
  //	VariableSizeMatrixType TestingData = mFeatureExtractionLocalPtr.FormulateTestData(globaltestintensities);
  //	VariableSizeMatrixType ScaledTestingData = mFeatureScalingLocalPtr.ScaleGivenTestingFeatures(TestingData);
  //
  //#ifdef APP_BASE_CAPTK_H
  //	progressUpdate(60);
  //#endif
  //
  //	//qSleep(5000);
  //	//typedef vnl_matrix<double> MatrixTypeT;
  //	//MatrixTypeT dataT;
  //	//dataT.set_size(1517, 45);
  //	//for (int ii = 0; ii < perfusionIntensities.size(); ii++)
  //	//	for (int jj = 0; jj < 45; jj++)
  //	//		dataT(ii, jj) = reducedTestPerfusionFeatures[ii][jj];
  //	//typedef itk::CSVNumericObjectFileWriter<double, 1517, 45> WriterTypeT;
  //	//WriterTypeT::Pointer writerT = WriterTypeT::New();
  //	//writerT->SetFileName("perfusionDataTest.csv");
  //	//writerT->SetInput(&dataT);
  //	//writerT->Write();
  //
  //
  //	//dataT.set_size(1517, 15);
  //	//for (int ii = 0; ii < TestingData.Rows(); ii++)
  //	//	for (int jj = 0; jj < TestingData.Cols(); jj++)
  //	//		dataT(ii, jj) = TestingData[ii][jj];
  //	//writerT->SetFileName("tDataTest.csv");
  //	//writerT->SetInput(&dataT);
  //	//writerT->Write();
  //
  //	//dataT.set_size(1517, 15);
  //	//for (int ii = 0; ii < ScaledTestingData.Rows(); ii++)
  //	//	for (int jj = 0; jj < ScaledTestingData.Cols(); jj++)
  //	//		dataT(ii, jj) = ScaledTestingData[ii][jj];
  //	//writerT->SetFileName("sDataTest.csv");
  //	//writerT->SetInput(&dataT);
  //	//writerT->Write();
  //
  //	try
  //	{
  //		VectorDouble result = testOpenCVSVM(ScaledTestingData, mOutputLocalPtr.mOutputDirectoryPath + "/" + mTrainedModelNameXML);
  //		//--------------------------------------------------------
  //		VectorDouble positiveDistances;
  //		VectorDouble negativeDistances;
  //		for (size_t x = 0; x < trainDistances.size(); x++)
  //		{
  //			if (trainDistances[x] > 0)
  //				positiveDistances.push_back(trainDistances[x]);
  //			else if (trainDistances[x] < 0)
  //				negativeDistances.push_back(trainDistances[x]);
  //		}
  //		std::sort(positiveDistances.begin(), positiveDistances.end());
  //		std::sort(negativeDistances.rbegin(), negativeDistances.rend());
  //
  //		double positivePercentileCutOff = positiveDistances[positiveDistances.size() - 1];
  //		double negativePercentileCutOff = negativeDistances[negativeDistances.size() - 1];
  //
  //		for (size_t x = 0; x < positiveDistances.size(); x++)
  //		{
  //			double percentile = (100 * (x + 0.5)) / positiveDistances.size();
  //			if (percentile > 99)
  //			{
  //				positivePercentileCutOff = result[x];
  //				break;
  //			}
  //		}
  //		for (size_t x = 0; x < negativeDistances.size(); x++)
  //		{
  //			double percentile = (100 * (x + 0.5)) / negativeDistances.size();
  //			if (percentile > 99)
  //			{
  //				negativePercentileCutOff = result[x];
  //				break;
  //			}
  //		}
  //		for (size_t x = 0; x < result.size(); x++)
  //		{
  //			if (result[x] > positivePercentileCutOff)
  //				result[x] = positivePercentileCutOff;
  //			if (result[x] < negativePercentileCutOff)
  //				result[x] = negativePercentileCutOff;
  //		}
  //		double min = 0;
  //		double max = positivePercentileCutOff;
  //		for (size_t x = 0; x < result.size(); x++)
  //		{
  //			if (result[x] > 0)
  //				result[x] = (result[x] - min) / (max - min);
  //		}
  //		min = negativePercentileCutOff;
  //		max = 0;
  //		for (size_t x = 0; x < result.size(); x++)
  //		{
  //			if (result[x] < 0)
  //				result[x] = (result[x] - min) / (max - min);
  //		}
  //		for (int x = 0; x < result[x]; x++)
  //			result[x] = (result[x] + 1) / 2;
  //
  //		//--------------------------------------------------------
  //		//mClassificationLocalPtr.SetModelFileName(mOutputLocalPtr.mOutputDirectoryPath + "/FinalModelFile.model");
  //		//VectorVectorDouble result = mClassificationLocalPtr.Testing(ScaledTestingData, false, mOutputLocalPtr.mOutputDirectoryPath +  "/FinalModelFile.model");
  //#ifdef APP_BASE_CAPTK_H
  //		progressUpdate(70);
  //#endif
  //
  //		for (unsigned int x = 0; x < result.size(); x++)
  //			RecurrenceProbabilityMap->SetPixel(testindices[x], result[x]);
  //
  //#ifdef APP_BASE_CAPTK_H
  //		progressUpdate(80);
  //#endif
  //		typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
  //		IteratorType RecIt(RecurrenceProbabilityMap, RecurrenceProbabilityMap->GetLargestPossibleRegion());
  //		IteratorType EdeIt(edemaMask, edemaMask->GetLargestPossibleRegion());
  //		RecIt.GoToBegin();
  //		EdeIt.GoToBegin();
  //		while (!RecIt.IsAtEnd())
  //		{
  //			typename ImageType::IndexType index = RecIt.GetIndex();
  //			if (EdeIt.Get() != VOXEL_STATUS::ON)
  //				RecIt.Set(0);
  //			++RecIt;
  //			++EdeIt;
  //		}
  //		//if (x_index >= mBorderStartX && x_index <= mBorderEndX && y_index >= mBorderStartY && y_index <= mBorderEndY && z_index >= 80 && z_index <= 120)
  //	}
  //	catch (itk::ExceptionObject & excp)
  //	{
  //		std::string str(excp.GetDescription());
  //	}
  //#ifdef APP_BASE_CAPTK_H
  //	progressUpdate(90);
  //#endif
  //	//------------------------------------------Writing final output--------------------------------------------------
  //	mRecurrenceMapFileName = "RecurrenceMap.nii.gz";
  //	mOutputLocalPtr.WriteRecurrenceOutput<ImageType>(RecurrenceProbabilityMap, t1cebasefilename, mRecurrenceMapFileName);
  //
  //#ifdef APP_BASE_CAPTK_H
  //	progressUpdate(100);
  //#endif
}


template<class ImageType>
void PseudoProgressionEstimator::PseudoProgressionEstimateOnGivenSubjectUsingExistingModel(typename ImageType::Pointer edemaMask,
  typename ImageType::Pointer tumorMask,
  typename ImageType::Pointer FinalT1CEImagePointer,
  typename ImageType::Pointer FinalT2FlairImagePointer,
  typename ImageType::Pointer FinalT1ImagePointer,
  typename ImageType::Pointer FinalT2ImagePointer,
  std::vector<typename ImageType::Pointer> FinalPerfusionImagePointer,
  std::vector<typename ImageType::Pointer> FinalDTIImagePointer,
  int imagetype, bool conDataPresent, bool dtiDataPresent, bool perfusionDataPresent, bool distanceDataPresent,
  bool useOtherModalities, std::string t1cebasefilename,
  const VectorVectorDouble &nearRegionIndices,
  const VectorVectorDouble &farRegionIndices, const std::string modeldirectory)
{
  cbica::Logging(loggerFile, "mCurrentOutputDir::" + mCurrentOutputDir);
  mOutputLocalPtr.SetOutputDirectoryPath(mCurrentOutputDir);
  VariableSizeMatrixType pca_coefficients;
  VariableLengthVectorType pca_mean;
  VariableLengthVectorType mean;
  VariableLengthVectorType stds;

  //check for the presence of model file
  if (!cbica::fileExists(modeldirectory + "/" + mPseudoTrainedFile) || !cbica::fileExists(modeldirectory + "/" + mRecurrenceTrainedFile))
  {
    mLastEncounteredError = "SVM model file is not present in the directory: " + modeldirectory;
    return;
  }

  //check for the file having number of modalities
  //if (!cbica::fileExists(modeldirectory + "/Recurrence_SVM_Modalities.csv"))
  //{
  //	mLastEncounteredError = "File having record of modalities is not present in the directory: " + modeldirectory;
  //	return;
  //}
  //else
  //{
  //	vnl_matrix<double> modalities = mOutputLocalPtr.ReadNumberOfModalities(modeldirectory + "/Recurrence_SVM_Modalities.csv");
  //	std::string message = "Model supports following modalities:\n";
  //	if (modalities[0][0] == 1)
  //		message += "Conventional\n";
  //	if (modalities[0][1] == 1)
  //		message += "DTI\n";
  //	if (modalities[0][2] == 1)
  //		message += "Perfusion\n";
  //	if (modalities[0][3] == 1)
  //		message += "Distance trnaform";
  //	if ((useConventionalrData == true && modalities[0][0] == 0) || (useDTIData == true && modalities[0][1] == 0) || (usePerfData == true && modalities[0][2] == 0) ||
  //		(useDistData == true && modalities[0][3] == 0))
  //	{
  //		mLastEncounteredError = message;
  //		return;
  //	}
  //}

  //check for the presence of z-score record
  if (!cbica::fileExists(modeldirectory + "/Recurrence_ZScore_Mean.csv") || !cbica::fileExists(modeldirectory + "/Recurrence_ZScore_Std.csv"))
  {
    mLastEncounteredError = "Z-score record is not present in the directory: " + modeldirectory;
    return;
  }

  //check for the presence of PCA related record
  //if (usePerfData)
  //{
  if (!cbica::fileExists(modeldirectory + "/Recurrence_COEF.csv") || !cbica::fileExists(modeldirectory + "/Recurrence_MR.csv"))
  {
    mLastEncounteredError = "PCA parameters are not present in the directory: " + modeldirectory;
    return;
  }
  else
  {
    mOutputLocalPtr.ReadModelParameters(modeldirectory + "/Recurrence_ZScore_Mean.csv", modeldirectory + "/Recurrence_ZScore_Std.csv", modeldirectory + "/Recurrence_COEF.csv", modeldirectory + "/Recurrence_MR.csv", mean, stds, pca_coefficients, pca_mean);
    mFeatureReductionLocalPtr.SetParameters(pca_coefficients, pca_mean);
    mFeatureScalingLocalPtr.SetParameters(mean, stds);
  }
  //}
  //else
  //{
  //	mOutputLocalPtr.ReadModelParameters(modeldirectory + "/Recurrence_ZScore_Mean.csv", modeldirectory + "/Recurrence_ZScore_Std.csv", mean, stds);
  //	mFeatureScalingLocalPtr.SetParameters(mean, stds);
  //}
  cbica::Logging(loggerFile, "Model reading finished.");


  ////std::string modelFileName = "";
  ////std::string modelMeanName = "";
  ////std::string modelStdName = "";
  ////std::string modelCoefName = "";
  ////std::string modelMRName = "";

  ////if (QFile(QApplication::applicationFilePath() + "/../data/recurrence/Recurrence_SVM_Model.csv").exists() &&
  ////	QFile(QApplication::applicationFilePath() + "/../data/recurrence/Recurrence_SVM_Modalities.csv").exists() &&
  ////	QFile(QApplication::applicationFilePath() + "/../data/recurrence/Recurrence_ZScore_Mean.csv").exists() &&
  ////	QFile(QApplication::applicationFilePath() + "/../data/recurrence/Recurrence_ZScore_Std.csv").exists() &&
  ////	QFile(QApplication::applicationFilePath() + "/../data/recurrence/Recurrence_MR.csv").exists() &&
  ////	QFile(QApplication::applicationFilePath() + "/../data/recurrence/Recurrence_COEF.csv").exists())
  ////{
  ////	modelFileName = QApplication::applicationFilePath().toStdString() + "/../data/recurrence/" + mTrainedModelNameCSV;
  ////	modelMeanName = QApplication::applicationFilePath().toStdString() + "/../data/recurrence/Recurrence_ZScore_Mean.csv";
  ////	modelStdName = QApplication::applicationFilePath().toStdString() + "/../data/recurrence/Recurrence_ZScore_Std.csv";
  ////	modelCoefName = QApplication::applicationFilePath().toStdString() + "/../data/recurrence/Recurrence_COEF.csv";
  ////	modelMRName = QApplication::applicationFilePath().toStdString() + "/../data/recurrence/Recurrence_MR.csv";
  ////}
  ////else
  ////{
  ////	cbica::Logging(loggerFile, "Model files not present.");
  ////	return;
  ////}

  //////check for the presence of model file
  ////if (!cbica::fileExists("../data/recurrence/" + mTrainedModelNameCSV))
  ////{
  ////		mLastEncounteredError = "SVM model file is not present in the directory: " + modeldirectory;
  ////		return;
  ////}
  //////check for the file having number of modalities
  ////if (!cbica::fileExists("../data/recurrence/Recurrence_SVM_Modalities.csv"))
  ////{
  ////	mLastEncounteredError = "File having record of modalities is not present in the directory: " + modeldirectory;
  ////	return;
  ////}
  ////else
  ////{
  ////	vnl_matrix<double> modalities = mOutputLocalPtr.ReadNumberOfModalities("../data/recurrence/Recurrence_SVM_Modalities.csv");
  ////	std::string message = "Model supports following modalities:\n";
  ////	if (modalities[0][0] == 1)
  ////		message += "Conventional\n";
  ////	if (modalities[0][1] == 1)
  ////		message += "DTI\n";
  ////	if (modalities[0][2] == 1)
  ////		message += "Perfusion\n";
  ////	if (modalities[0][3] == 1)
  ////		message += "Distance trnaform";
  ////	if ((useConventionalrData == true && modalities[0][0] == 0) || (useDTIData == true && modalities[0][1] == 0) || (usePerfData == true && modalities[0][2] == 0) ||
  ////		(useDistData == true && modalities[0][3] == 0))
  ////	{
  ////		mLastEncounteredError = message;
  ////		return;
  ////	}
  ////}
  ////check for the presence of z-score record
  ////if (!cbica::fileExists("../data/recurrence/Recurrence_ZScore_Mean.csv") || !cbica::fileExists("../data/recurrence/Recurrence_ZScore_Std.csv"))
  ////{
  ////	mLastEncounteredError = "Z-score record is not present in the directory: " + modeldirectory;
  ////	return;
  ////}

  //////check for the presence of PCA related record
  ////if (usePerfData)
  ////{
  ////	if (!cbica::fileExists("../data/recurrence/Recurrence_COEF.csv") || !cbica::fileExists("../data/recurrence/Recurrence_MR.csv"))
  ////	{
  ////		mLastEncounteredError = "PCA parameters are not present in the directory: " + modeldirectory;
  ////		return;
  ////	}
  ////	else
  ////	{
  ////std::string data = cbica::normalizePath("../data/recurrence/Recurrence_ZScore_Mean.csv");

  ////cbica::Logging(loggerFile, "Path::" + data);


  ////mOutputLocalPtr.ReadModelParameters(modelMeanName, modelStdName, modelCoefName, modelMRName, mean, stds, pca_coefficients, pca_mean);
  //cbica::normalizePath("../data/recurrence/Recurrence_ZScore_Mean.csv"), 
  //cbica::normalizePath("../data/recurrence/Recurrence_ZScore_Std.csv"), 
  //cbica::normalizePath("../data/recurrence/Recurrence_COEF.csv"), 
  //cbica::normalizePath("../data/recurrence/Recurrence_MR.csv"), mean, stds, pca_coefficients, pca_mean);

  //mFeatureReductionLocalPtr.SetParameters(pca_coefficients, pca_mean);
  //mFeatureScalingLocalPtr.SetParameters(mean, stds);
  ////	}
  ////}
  ////else
  ////{
  ////	mOutputLocalPtr.ReadModelParameters("../data/recurrence/Recurrence_ZScore_Mean.csv", "../data/recurrence/Recurrence_ZScore_Std.csv", mean, stds);
  ////	mFeatureScalingLocalPtr.SetParameters(mean, stds);
  ////}


  ////------------------------------------------testing data formulation---------------------------------------------------





  VectorVectorDouble perfusionIntensities;
  VectorVectorDouble otherIntensities;
  VectorDouble distanceIntensities;
  std::vector<typename ImageType::IndexType> testindices;

  //typedef ImageTypeFloat3D OutputImageType;
  //OutputImageType::Pointer dilatedEdema;
  //typedef itk::BinaryBallStructuringElement<ImageType::PixelType, 3> StructuringElementType;
  //StructuringElementType structuringElement;
  //structuringElement.SetRadius(1);
  //structuringElement.CreateStructuringElement();
  //typedef itk::BinaryDilateImageFilter <ImageType, ImageType, StructuringElementType> BinaryDilateImageFilterType;
  //BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
  //dilateFilter->SetInput(edemaMask);
  //dilateFilter->SetKernel(structuringElement);
  //dilateFilter->SetDilateValue(VOXEL_STATUS::ON);
  //dilateFilter->Update();
  //dilatedEdema = dilateFilter->GetOutput();


  testindices = mNiftiLocalPtr.LoadTestData(FinalT1CEImagePointer, FinalT2FlairImagePointer, FinalT1ImagePointer, FinalT2ImagePointer, FinalPerfusionImagePointer, FinalDTIImagePointer, tumorMask, edemaMask, perfusionIntensities, otherIntensities, distanceIntensities, CAPTK::ImageExtension::DICOM, conDataPresent, dtiDataPresent, perfusionDataPresent, distanceDataPresent);

  if (testindices.empty())
  {
    cbica::Logging(loggerFile, "{RecurrDebug}testindices is empty");
  }

  VectorVectorDouble reducedTestPerfusionFeatures;
  reducedTestPerfusionFeatures = mFeatureReductionLocalPtr.ApplyPCAOnTestData(perfusionIntensities);

  if (reducedTestPerfusionFeatures.empty())
  {
    cbica::Logging(loggerFile, "{RecurrDebug}reducedTestPerfusionFeatures is empty");
  }

  int NumberOfPCs = 5;
  VectorVectorDouble globaltestintensities;

  for (unsigned int k = 0; k < testindices.size(); k++)
  {
    VectorDouble inten;

    for (int j = 0; j < NumberOfPCs; j++)
      inten.push_back(reducedTestPerfusionFeatures[k][j]);

    for (unsigned int j = 0; j < otherIntensities[0].size(); j++)
      inten.push_back(otherIntensities[k][j]);

    inten.push_back(distanceIntensities[k]);

    if (!inten.empty())
      globaltestintensities.push_back(inten);
  }

  if (globaltestintensities.empty())
  {
    cbica::Logging(loggerFile, "{RecurrDebug}globaltestintensities is empty");
  }

  VariableSizeMatrixType TestingData = mFeatureExtractionLocalPtr.FormulateTestData(globaltestintensities);

  //typedef vnl_matrix<double> MatrixType;
  //MatrixType data;
  //data.set_size(98721, 15);

  //for (unsigned int i = 0; i < TestingData.Rows(); i++)
  //	for (unsigned int j = 0; j < TestingData.Cols(); j++)
  //		data(i, j) = TestingData[i][j];
  //typedef itk::CSVNumericObjectFileWriter<double, 98721, 15> WriterType;
  //WriterType::Pointer writer = WriterType::New();
  //writer->SetFileName("testdata_loaded.csv");
  //writer->SetInput(&data);
  //try
  //{
  //	writer->Write();
  //}
  //catch (itk::ExceptionObject & excp)
  //{
  //	std::cerr << "Error: " << excp.what();
  //}



  VariableSizeMatrixType ScaledTestingData = mFeatureScalingLocalPtr.ScaleGivenTestingFeatures(TestingData);
  //for (int i = 0; i < TestingData.Rows(); i++)
  //	for (int j = 0; j < TestingData.Cols(); j++)
  //		data(i, j) = ScaledTestingData[i][j];
  //writer->SetFileName("scaledtestdata_loaded.csv");
  //writer->SetInput(&data);
  //try
  //{
  //	writer->Write();
  //}
  //catch (itk::ExceptionObject & excp)
  //{
  //}



  if (TestingData.Cols() == 0)
  {
    cbica::Logging(loggerFile, "{RecurrDebug}TestingData.Cols is empty");
  }

  if (ScaledTestingData.Cols() == 0)
  {
    cbica::Logging(loggerFile, "{RecurrDebug}ScaledTestingData.Cols is empty");
  }

#ifdef APP_BASE_CAPTK_H
  progressUpdate(60);
#endif

  //qSleep(5000);
  //typedef vnl_matrix<double> MatrixTypeT;
  //MatrixTypeT dataT;
  //dataT.set_size(48356, 45);
  //for (int ii = 0; ii < perfusionIntensities.size(); ii++)
  //	for (int jj = 0; jj < 45; jj++)
  //		dataT(ii, jj) = perfusionIntensities[ii][jj];
  //typedef itk::CSVNumericObjectFileWriter<double, 48356, 45> WriterTypeT;
  //WriterTypeT::Pointer writerT = WriterTypeT::New();
  //writerT->SetFileName("perfusionDataTest.csv");
  //writerT->SetInput(&dataT);
  //writerT->Write();


  //dataT.set_size(48356, 15);
  //for (int ii = 0; ii < TestingData.Rows(); ii++)
  //	for (int jj = 0; jj < TestingData.Cols(); jj++)
  //		dataT(ii, jj) = TestingData[ii][jj];
  //writerT->SetFileName("tDataTest.csv");
  //writerT->SetInput(&dataT);
  //writerT->Write();

  //dataT.set_size(48356, 15);
  //for (int ii = 0; ii < ScaledTestingData.Rows(); ii++)
  //	for (int jj = 0; jj < ScaledTestingData.Cols(); jj++)
  //		dataT(ii, jj) = ScaledTestingData[ii][jj];
  //writerT->SetFileName("sDataTest.csv");
  //writerT->SetInput(&dataT);
  //writerT->Write();

  typename ImageType::RegionType region = FinalT1CEImagePointer->GetLargestPossibleRegion();
  typename ImageType::Pointer RecProbabilityMap = ImageType::New();
  RecProbabilityMap->SetRegions(region);
  RecProbabilityMap->Allocate();
  RecProbabilityMap->SetSpacing(FinalT1CEImagePointer->GetSpacing());
  RecProbabilityMap->SetOrigin(FinalT1CEImagePointer->GetOrigin());
  RecProbabilityMap->SetDirection(FinalT1CEImagePointer->GetDirection());


  typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
  IteratorType imIt(edemaMask, edemaMask->GetLargestPossibleRegion());
  IteratorType RecIt(RecProbabilityMap, RecProbabilityMap->GetLargestPossibleRegion());
  imIt.GoToBegin(); RecIt.GoToBegin();

  while (!imIt.IsAtEnd())
  {
    if (imIt.Get() == 0)
      RecIt.Set(0);

    ++imIt;
    ++RecIt;
  }
  cbica::Logging(loggerFile, "Before testing.");
  VectorDouble result_modified;

  try
  {
    if (cbica::fileExists(modeldirectory + "/" + mPseudoTrainedFile))
    {
      //cbica::Logging(loggerFile, "Before testing 1.");
      VariableLengthVectorType result;
      result = DistanceFunction(ScaledTestingData, modeldirectory + "/" + mPseudoTrainedFile, RECURRENCE_MODEL_RHO, RECURRENCE_MODEL_G);
      for (unsigned int index = 0; index < result.Size(); index++)
      {
        RecProbabilityMap->SetPixel(testindices[index], result[index] * 1);
        result_modified.push_back(result[index]);
      }

    }
    else
    {
      //cbica::Logging(loggerFile, "Before testing 2.");
      VectorDouble result;
      result = testOpenCVSVM(ScaledTestingData, modeldirectory + "/" + mPseudoTrainedFile);
      for (unsigned int index = 0; index < result.size(); index++)
        RecProbabilityMap->SetPixel(testindices[index], result[index] * 1);
      result_modified = result;
    }
  }
  catch (itk::ExceptionObject & excp)
  {
    cbica::Logging(loggerFile, "Error caught during testing: " + std::string(excp.GetDescription()));
    exit(EXIT_FAILURE);
  }

  VectorDouble result_revised = RecurrenceMapPostprocessing<ImageType>(result_modified, testindices, RecProbabilityMap, edemaMask);
  for (unsigned int index = 0; index < result_modified.size(); index++)
    RecProbabilityMap->SetPixel(testindices[index], result_revised[index] * 1);

  //averaging filter
  typedef itk::MeanImageFilter<ImageType, ImageType > FilterType;
  typename FilterType::Pointer meanFilter = FilterType::New();
  typename FilterType::InputSizeType radius;
  radius.Fill(1);
  meanFilter->SetRadius(radius);
  meanFilter->SetInput(RecProbabilityMap);
  typename ImageType::Pointer RevisedRecurrenceMap = meanFilter->GetOutput();

  //------------------------------------------Writing final output--------------------------------------------------
  mRecurrenceMapFileName = "RecurrenceMap.nii.gz";
  mOutputLocalPtr.WriteRecurrenceOutput<ImageType>(RevisedRecurrenceMap, t1cebasefilename, mRecurrenceMapFileName);



  //	const typename ImageType::RegionType region = FinalT1CEImagePointer->GetLargestPossibleRegion();
  //	typename ImageType::Pointer RecurrenceProbabilityMap = ImageType::New();
  //	RecurrenceProbabilityMap->SetRegions(region);
  //	RecurrenceProbabilityMap->Allocate();
  //	RecurrenceProbabilityMap->SetSpacing(FinalT1CEImagePointer->GetSpacing());
  //	RecurrenceProbabilityMap->SetOrigin(FinalT1CEImagePointer->GetOrigin());
  //	RecurrenceProbabilityMap->SetDirection(FinalT1CEImagePointer->GetDirection());
  //	try
  //	{
  //		if (cbica::fileExists(modeldirectory + "/" + mTrainedModelNameCSV))
  //		{
  //			VariableLengthVectorType result;
  //			result = DistanceFunction(ScaledTestingData, modeldirectory + "/" + mTrainedModelNameCSV, RECURRENCE_MODEL_RHO, RECURRENCE_MODEL_G);
  //			for (unsigned int index = 0; index < result.Size(); index++)
  //				RecurrenceProbabilityMap->SetPixel(testindices[index], result[index] * 1);
  //			if (RecurrenceProbabilityMap.IsNull())
  //			{
  //				cbica::Logging(loggerFile, "{RecurrDebug}RecurrenceProbabilityMap is Null");
  //			}
  //		}
  //		else
  //		{
  //			VectorDouble result;
  //			result = testOpenCVSVM(ScaledTestingData, modeldirectory + "/" + mTrainedModelNameXML);
  //			for (unsigned int index = 0; index < result.size(); index++)
  //				RecurrenceProbabilityMap->SetPixel(testindices[index], result[index] * 1);	
  //		}
  //
  //		double minvalue = 0;
  //		//for (unsigned int index = 0; index < result.Size(); index++)
  //		//{
  //		//	if (result[index] < minvalue)
  //		//		minvalue = result[index];
  //		//}
  //
  //		//dataT.set_size(48356, 1);
  //		//for (int ii = 0; ii < result.Size(); ii++)
  //		//		dataT(ii, 0) = result[ii];
  //		//writerT->SetFileName("results.csv");
  //		//writerT->SetInput(&dataT);
  //		//writerT->Write();
  //
  //
  //		//VectorDouble result = testOpenCVSVM(ScaledTestingData, mOutputLocalPtr.mOutputDirectoryPath + "/" + mTrainedModelNameXML);
  //		//--------------------------------------------------------
  //		//VectorDouble positiveDistances;
  //		//VectorDouble negativeDistances;
  //		//for (size_t x = 0; x < trainDistances.size(); x++)
  //		//{
  //		//	if (trainDistances[x] > 0)
  //		//		positiveDistances.push_back(trainDistances[x]);
  //		//	else if (trainDistances[x] < 0)
  //		//		negativeDistances.push_back(trainDistances[x]);
  //		//}
  //		//std::sort(positiveDistances.begin(), positiveDistances.end());
  //		//std::sort(negativeDistances.rbegin(), negativeDistances.rend());
  //
  //		//double positivePercentileCutOff = positiveDistances[positiveDistances.size() - 1];
  //		//double negativePercentileCutOff = negativeDistances[negativeDistances.size() - 1];
  //
  //		//for (size_t x = 0; x < positiveDistances.size(); x++)
  //		//{
  //		//	double percentile = (100 * (x + 0.5)) / positiveDistances.size();
  //		//	if (percentile > 99)
  //		//	{
  //		//		positivePercentileCutOff = result[x];
  //		//		break;
  //		//	}
  //		//}
  //		//for (size_t x = 0; x < negativeDistances.size(); x++)
  //		//{
  //		//	double percentile = (100 * (x + 0.5)) / negativeDistances.size();
  //		//	if (percentile > 99)
  //		//	{
  //		//		negativePercentileCutOff = result[x];
  //		//		break;
  //		//	}
  //		//}
  //		//for (size_t x = 0; x < result.size(); x++)
  //		//{dd
  //		//	if (result[x] > positivePercentileCutOff)
  //		//		result[x] = positivePercentileCutOff;
  //		//	if (result[x] < negativePercentileCutOff)
  //		//		result[x] = negativePercentileCutOff;
  //		//}
  //		//double min = 0;
  //		//double max = positivePercentileCutOff;
  //		//for (size_t x = 0; x < result.size(); x++)
  //		//{
  //		//	if (result[x] > 0)
  //		//		result[x] = (result[x] - min) / (max - min);
  //		//}
  //		//min = negativePercentileCutOff;
  //		//max = 0;
  //		//for (size_t x = 0; x < result.size(); x++)
  //		//{
  //		//	if (result[x] < 0)
  //		//		result[x] = (result[x] - min) / (max - min);
  //		//}
  //		//for (int x = 0; x < result[x]; x++)
  //		//	result[x] = (result[x] + 1) / 2;
  //
  //		////--------------------------------------------------------
  //		////mClassificationLocalPtr.SetModelFileName(mOutputLocalPtr.mOutputDirectoryPath + "/FinalModelFile.model");
  //		////VectorVectorDouble result = mClassificationLocalPtr.Testing(ScaledTestingData, false, mOutputLocalPtr.mOutputDirectoryPath +  "/FinalModelFile.model");
  //#ifdef APP_BASE_CAPTK_H
  //		progressUpdate(70);
  //#endif
  //
  //		//for (unsigned int x = 0; x < result.size(); x++)
  //		//	RecurrenceProbabilityMap->SetPixel(testindices[x], result[x]);
  //
  //#ifdef APP_BASE_CAPTK_H
  //		progressUpdate(80);
  //#endif
  //
  //
  //		//typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
  //		//IteratorType RecIt(RecurrenceProbabilityMap, RecurrenceProbabilityMap->GetLargestPossibleRegion());
  //		//IteratorType EdeIt(edemaMask, edemaMask->GetLargestPossibleRegion());
  //
  //		//RecIt.GoToBegin();
  //		//while (!RecIt.IsAtEnd())
  //		//{
  //		//	if (RecIt.Get() < minvalue)
  //		//		minvalue = RecIt.Get();
  //		//	++RecIt;
  //		//}
  //
  //
  //		//RecIt.GoToBegin();
  //		//EdeIt.GoToBegin();
  //		//while (!RecIt.IsAtEnd())
  //		//{
  //		//	typename ImageType::IndexType index = RecIt.GetIndex();
  //		//	if (EdeIt.Get() != VOXEL_STATUS::ON)
  //		//		RecIt.Set(minvalue);
  //		//	++RecIt;
  //		//	++EdeIt;
  //		//}
  //		//if (x_index >= mBorderStartX && x_index <= mBorderEndX && y_index >= mBorderStartY && y_index <= mBorderEndY && z_index >= 80 && z_index <= 120)
  //	}
  //	catch (itk::ExceptionObject & excp)
  //	{
  //		std::string str(excp.GetDescription());
  //	}
#ifdef APP_BASE_CAPTK_H
  progressUpdate(90);
#endif

#ifdef APP_BASE_CAPTK_H
  progressUpdate(100);
#endif
}




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
template <class ImageType>
typename ImageType::Pointer PseudoProgressionEstimator::ReadNiftiImage(const std::string &filename)
{
  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  typename ImageReaderType::Pointer reader = ImageReaderType::New();
  reader->SetFileName(filename);

  try
  {
    reader->Update();
  }
  catch (itk::ExceptionObject& e)
  {
    std::cerr << "Error caught: " << e.what() << "\n";
    //cbica::Logging(loggerFile, "Error caught during testing: " + std::string(e.GetDescription()));
    exit(EXIT_FAILURE);
  }


  return reader->GetOutput();
}

template<class ImageType>
VectorDouble PseudoProgressionEstimator::RecurrenceMapPostprocessing(VectorDouble result, std::vector<typename ImageType::IndexType> testindices, typename ImageType::Pointer RecurrenceProbabilityMap, typename ImageType::Pointer edemaMap)
{
  VectorDouble revised_result;
  VectorDouble positiveDistances;
  VectorDouble negativeDistances;
  for (size_t x = 0; x < result.size(); x++)
  {
    if (result[x] > 0)
      positiveDistances.push_back(result[x]);
    else
      negativeDistances.push_back(result[x]);
  }
  std::sort(positiveDistances.begin(), positiveDistances.end());
  std::sort(negativeDistances.rbegin(), negativeDistances.rend());

  double positivePercentileCutOff = positiveDistances[positiveDistances.size() - 1];
  double negativePercentileCutOff = negativeDistances[negativeDistances.size() - 1];

  for (size_t x = 0; x < positiveDistances.size(); x++)
  {
    double percentile = (100 * (x + 0.5)) / positiveDistances.size();
    if (percentile > 99)
    {
      positivePercentileCutOff = positiveDistances[x];
      break;
    }
  }
  for (size_t x = 0; x < negativeDistances.size(); x++)
  {
    double percentile = (100 * (x + 0.5)) / negativeDistances.size();
    if (percentile > 99)
    {
      negativePercentileCutOff = negativeDistances[x];
      break;
    }
  }
  for (size_t x = 0; x < result.size(); x++)
  {
    if (result[x] > positivePercentileCutOff)
      result[x] = positivePercentileCutOff;
    if (result[x] < negativePercentileCutOff)
      result[x] = negativePercentileCutOff;
  }
  double min = 0;
  double max = positivePercentileCutOff;
  for (size_t x = 0; x < result.size(); x++)
  {
    if (result[x] > 0)
      result[x] = (result[x] - min) / (max - min);
  }
  min = negativePercentileCutOff;
  max = 0;
  for (size_t x = 0; x < result.size(); x++)
  {
    if (result[x] < 0)
      result[x] = (result[x] - min) / (max - min);
  }
  for (int x = 0; x < result[x]; x++)
    result[x] = (result[x] + 1) / 2;

  revised_result = result;
  return revised_result;



  ////eroding edema
  //typename ImageType::Pointer dilatedEdema;
  //typedef itk::BinaryBallStructuringElement<typename ImageType::PixelType, 3> StructuringElementType;
  //StructuringElementType structuringElement;
  //structuringElement.SetRadius(2);
  //structuringElement.CreateStructuringElement();
  //typedef itk::BinaryErodeImageFilter <ImageType, ImageType, StructuringElementType> BinaryErodeImageFilterType;
  //typename BinaryErodeImageFilterType::Pointer dilateFilter = BinaryErodeImageFilterType::New();
  //dilateFilter->SetInput(edemaMap);
  //dilateFilter->SetKernel(structuringElement);
  //dilateFilter->SetErodeValue(1);
  //dilateFilter->Update();
  //dilatedEdema = dilateFilter->GetOutput();

  ////find voxels which are not in eroded in edema
  //typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
  //IteratorType RecIt(RecurrenceProbabilityMap, RecurrenceProbabilityMap->GetLargestPossibleRegion());
  //IteratorType EdeIt(edemaMap, edemaMap->GetLargestPossibleRegion());
  //IteratorType DilIt(dilatedEdema, dilatedEdema->GetLargestPossibleRegion());
  //int counter1 = 0;
  //int counter2 = 0; int counter3 = 0;
  //int counter4 = 0;

  //RecIt.GoToBegin();
  //EdeIt.GoToBegin();
  //DilIt.GoToBegin();
  //while (!RecIt.IsAtEnd())
  //{
  //	if (EdeIt.Get() == 1)
  //		counter1++;

  //	if (DilIt.Get() == 0)
  //		counter2++;


  //	if (DilIt.Get() < 0)
  //		counter3++;


  //	if (DilIt.Get() == 1)
  //		counter4++;



  //	if (EdeIt.Get() > 0 && DilIt.Get() != 1)
  //	{
  //		if (RecIt.Get() > 0.7)
  //			RecIt.Set(0.1);
  //	}
  //	++RecIt;
  //	++EdeIt;
  //	++DilIt;
  //}
  //cbica::WriteImage<ImageType>(RecurrenceProbabilityMap, "BeforeMean.nii.gz");

  //typename ImageType::Pointer meanRecurrenceMap;
  //typedef itk::MeanImageFilter<ImageType, ImageType> filterType;
  //typename filterType::Pointer meanFilter = filterType::New();
  //meanFilter->SetInput(RecurrenceProbabilityMap);
  //meanFilter->SetRadius(2);
  //meanRecurrenceMap = meanFilter->GetOutput();
  //cbica::WriteImage<ImageType>(meanRecurrenceMap, "AfterMean.nii.gz");

  //RecurrenceIt.GoToBegin();
  //EdeIt.GoToBegin();
  //while (!RecurrenceIt.IsAtEnd())
  //{
  //	if (EdeIt.Get() == 0)
  //		RecurrenceIt.Set(minvalue);
  //	++RecurrenceIt;
  //	++EdeIt;
  //}


  //averaging filter
  //typedef itk::MeanImageFilter<ImageType, ImageType > FilterType;
  //FilterType::Pointer meanFilter = FilterType::New();
  //FilterType::InputSizeType radius;
  //radius.Fill(2);
  //meanFilter->SetRadius(radius);
  //meanFilter->SetInput(RecurrenceProbabilityMap);
  //typename ImageType::Pointer RevisedRecurrenceMap = meanFilter->GetOutput();





  //for (int i = 0; i < indices.size();i++)
  //	RevisedRecurrenceMap->SetPixel(indices[i],minvalue
  //IteratorType EdeIt(edemaMask, edemaMask->GetLargestPossibleRegion());
  //ErodedEdeIt.GoToBegin();
  //EdeIt.GoToBegin();
  //while (!ErodedEdeIt.IsAtEnd())
  //{
  //	++ErodedEdeIt;
  //	++EdeIt;
  //}
  //	return meanRecurrenceMap;
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

  double minvalue = *min_element(ROIIntensities.begin(), ROIIntensities.end());
  double maxvalue = *max_element(ROIIntensities.begin(), ROIIntensities.end());

  VectorDouble HistogramFeatures = GetHistogramFeatures(ROIIntensities, 10);
  VectorDouble HistogramFeatures1 = GetHistogramFeatures(ROIIntensities, 20);
  VectorDouble IntensityFeatures = GetIntensityFeatures(ROIIntensities);
  typename ImageType::Pointer isotropicImageGLCM = cbica::ResampleImage< ImageType >(image, 1.0, "Linear");
  typename ImageType::Pointer isotropicMaskGLCM  = cbica::ResampleImage< ImageType >(mask, 1.0, "Nearest");
  auto roundingFilter = itk::RoundImageFilter<ImageType,ImageType>::New();
  roundingFilter->SetInput(isotropicMaskGLCM);
  roundingFilter->Update();
  isotropicMaskGLCM = roundingFilter->GetOutput();
  VectorDouble GLCMFeatures = GetGLCMFeatures<ImageType>(isotropicImageGLCM, isotropicMaskGLCM,minvalue,maxvalue);

  typename ImageType::Pointer isotropicImageGLRLM = cbica::ResampleImage< ImageType >(image, 1.0, "Linear");
  typename ImageType::Pointer isotropicMaskGLRLM = cbica::ResampleImage< ImageType >(mask, 1.0, "Nearest");
  roundingFilter->SetInput(isotropicMaskGLRLM);
  roundingFilter->Update();
  isotropicMaskGLRLM = roundingFilter->GetOutput();
  VectorDouble GLRLMFeatures = GetRunLengthFeatures<ImageType>(isotropicImageGLRLM, isotropicMaskGLRLM, minvalue, maxvalue);

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
VectorDouble PseudoProgressionEstimator::GetGLCMFeatures(typename ImageType::Pointer image, typename ImageType::Pointer mask, double minvalue, double maxvalue)
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
  auto centerIndex = neighborhood.GetCenterNeighborhoodIndex();

  for (int d = directionsToCompute - 1; d >= 0; d--)
  {
    if (d != static_cast<int>(centerIndex))
    {
      offsets->push_back(neighborhood.GetOffset(d));
    }
  }


  std::vector<double> features;
  double contrast = 0, correl = 0, ener = 0, entro = 0, homo = 0, clustershade = 0, clusterprominance = 0, autocorr = 0;

  for (size_t i = 0; i < offsets->size(); i++)
  {
    auto glcmGenerator = Image2CoOccuranceType::New();
    glcmGenerator->SetNumberOfBinsPerAxis(m_Bins); //reasonable number of bins
    glcmGenerator->SetPixelValueMinMax(minvalue,maxvalue);
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
  return features;
}

template<class ImageType>
VectorDouble PseudoProgressionEstimator::GetRunLengthFeatures(typename ImageType::Pointer image, typename ImageType::Pointer mask, double minvalue, double maxvalue)
{
  double m_Bins = 16;
  using HistogramFrequencyContainerType = itk::Statistics::DenseFrequencyContainer2;
  using RunLengthFilterType = itk::Statistics::EnhancedScalarImageToRunLengthFeaturesFilter< ImageType, HistogramFrequencyContainerType >;
  using RunLengthMatrixGenerator = typename RunLengthFilterType::RunLengthMatrixFilterType;
  using RunLengthFeatures = typename RunLengthFilterType::RunLengthFeaturesFilterType;
  using OffsetType = typename ImageType::OffsetType;
  using OffsetVector = itk::VectorContainer< unsigned char, OffsetType >;

  //offset calculation
  double inputDirections = 13;
  itk::Neighborhood< typename ImageType::PixelType, ImageType::ImageDimension > neighborhood;
  neighborhood.SetRadius(1);
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
  auto centerIndex = neighborhood.GetCenterNeighborhoodIndex();
  for (int d = directionsToCompute - 1; d >= 0; d--)
  {
    if (d != static_cast<int>(centerIndex))
    {
      offsets->push_back(neighborhood.GetOffset(d));
    }
  }

  //matrix generation
  typename  RunLengthMatrixGenerator::Pointer matrix_generator = RunLengthMatrixGenerator::New();
  matrix_generator->SetInput(image);
  matrix_generator->SetMaskImage(mask);
  matrix_generator->SetInsidePixelValue(1);
  matrix_generator->SetPixelValueMinMax(minvalue,maxvalue);

  matrix_generator->SetDistanceValueMinMax(0,10);
  matrix_generator->SetNumberOfBinsPerAxis(m_Bins);
  typename  RunLengthFeatures::Pointer runLengthMatrixCalculator = RunLengthFeatures::New();
  typename  RunLengthFeatures::Pointer runLengthFeaturesCalculator = RunLengthFeatures::New();
  typename  OffsetVector::ConstIterator offsetIt;
  size_t offsetNum = 0;

  double sre = 0, lre = 0, gln = 0, glnn = 0, rln = 0, rlnn = 0, rp = 0, lglre = 0, hglre = 0, srlgle = 0, srhgle = 0, lrlgle = 0, lrhgle = 0;
  std::vector<double> features;

  for (offsetIt = offsets->Begin(); offsetIt != offsets->End(); offsetIt++, offsetNum++)
  {
    matrix_generator->SetOffset(offsetIt.Value());
    matrix_generator->Update();
    
    //auto outputglcmGenerator = matrix_generator->GetOutput();
    //auto iterator = outputglcmGenerator->Begin();
    //while (iterator != outputglcmGenerator->End())
    //{
    //  std::cout << iterator.GetIndex() << " " << iterator.GetFrequency() << std::endl;
    //  ++iterator;
    //}

    runLengthFeaturesCalculator->SetInput(matrix_generator->GetOutput());
    runLengthFeaturesCalculator->Update();


    sre += runLengthFeaturesCalculator->GetShortRunEmphasis();
    lre += runLengthFeaturesCalculator->GetLongRunEmphasis();
    gln += runLengthFeaturesCalculator->GetGreyLevelNonuniformity();
    rln += runLengthFeaturesCalculator->GetRunLengthNonuniformity();
    lglre += runLengthFeaturesCalculator->GetLowGreyLevelRunEmphasis();
    hglre += runLengthFeaturesCalculator->GetHighGreyLevelRunEmphasis();
    srlgle += runLengthFeaturesCalculator->GetShortRunLowGreyLevelEmphasis();
    srhgle += runLengthFeaturesCalculator->GetShortRunHighGreyLevelEmphasis();
    lrlgle += runLengthFeaturesCalculator->GetLongRunLowGreyLevelEmphasis();
    lrhgle += runLengthFeaturesCalculator->GetLongRunHighGreyLevelEmphasis();
  }
  sre /= offsets->size();
  lre /= offsets->size();
  gln /= offsets->size();
  rln /= offsets->size();
  lglre /= offsets->size();
  hglre /= offsets->size();
  srlgle /= offsets->size();
  srhgle /= offsets->size();
  lrlgle /= offsets->size();
  lrhgle /= offsets->size();

  features.push_back(sre);
  features.push_back(lre);
  features.push_back(gln);
  features.push_back(rln);
  features.push_back(lglre);
  features.push_back(hglre);
  features.push_back(srlgle);
  features.push_back(srhgle);
  features.push_back(lrlgle);
  features.push_back(lrhgle);

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

