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
#include "CaPTkClassifierUtils.h"
#include "cbicaLogging.h"

#define RECURRENCE_MODEL_G 0.5
#define RECURRENCE_MODEL_RHO 0.0896

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
  VectorDouble GetGLCMFeatures(typename ImageType::Pointer image, typename ImageType::Pointer mask);

  template<class ImageType>
  typename ImageType::Pointer MakeAdditionalModality(typename ImageType::Pointer image1, typename ImageType::Pointer image2);


  template<class ImageType>
  VectorDouble GetRunLengthFeatures(typename ImageType::Pointer image, typename ImageType::Pointer mask);

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

