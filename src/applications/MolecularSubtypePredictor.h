/**
\file  MolecularSubtypePredictor.h

\brief The header file containing the SurvivaPredictor class, used to find Survival Prediction Index
Library Dependencies: ITK 4.7+ <br>

https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html

*/


#ifndef _MolecularSubtypePredictor_h_
#define _MolecularSubtypePredictor_h_

//#include "CAPTk.h"
#include "NiftiDataManager.h"
#include "FeatureReductionClass.h"
#include "FeatureScalingClass.h"
#include "TrainingModule.h"
#include "FeatureExtractionClass.h"
#include "itkCSVArray2DFileReader.h"
#include "itkConnectedComponentImageFilter.h"
#include "CaPTkEnums.h"
#include "cbicaLogging.h"
#include "CaPTkGUIUtils.h"
#include "CaPTkUtils.h"

#ifdef APP_BASE_CAPTK_H
#include "ApplicationBase.h"
#endif

typedef itk::Image< float, 3 > ImageType;
typedef itk::CSVArray2DFileReader<double> CSVFileReaderType;
typedef vnl_matrix<double> MatrixType;

#define SURVIVAL_SIZE_COMP 100
#define SURVIVAL_NO_NEAR_DIST 100
#define MOLECULAR_NO_OF_FEATURES 161

/**
\class MolecularSubtypePredictor

\brief Calculates Survival Prediction Index

Reference:

@article{macyszyn2015imaging,
title={Imaging patterns predict patient survival and molecular subtype in glioblastoma via machine learning techniques},
author={Macyszyn, Luke and Akbari, Hamed and Pisapia, Jared M and Da, Xiao and Attiah, Mark and Pigrish, Vadim and Bi, Yingtao and Pal, Sharmistha and Davuluri, Ramana V and Roccograndi, Laura and others},
journal={Neuro-oncology},
volume={18},
number={3},
pages={417--425},
year={2015},
publisher={Society for Neuro-Oncology}
}

*/
class MolecularSubtypePredictor
#ifdef APP_BASE_CAPTK_H
	: public ApplicationBase
#endif
{
public:
	//! Default constructor
	MolecularSubtypePredictor()
	{
		mProneuralTrainedFile = "ProneuralModelFile.xml";
		mNeuralTrainedFile = "NeuralModelFile.xml";
		mClassicalTrainedFile = "ClassicalModelFile.xml";
		mMessenchymalTrainedFile = "MessenchymalModelFile.xml";
		mLastErrorMessage = "";
		logger.UseNewFile(loggerFile);
	}

	//! Default destructor
	~MolecularSubtypePredictor() {};

	std::string mProneuralTrainedFile, mNeuralTrainedFile, mClassicalTrainedFile, mMessenchymalTrainedFile;

	FeatureExtractionClass mFeatureExtractionLocalPtr;
	FeatureScalingClass mFeatureScalingLocalPtr;
	cbica::Logging logger;

	std::string mLastErrorMessage;

	//template<class ImageType>
	//typename ImageType::Pointer CalculateSERFeature(const typename ImageType::Pointer &enhancement1, const typename ImageType::Pointer &enhancement2);

	//template<class ImageType>
	//typename ImageType::Pointer CalculatePEFeature(const typename ImageType::Pointer &enhancement1, const typename ImageType::Pointer &enhancement2);

	//template<class ImageType>
	//typename ImageType::Pointer CalculateWOSFeature(const typename ImageType::Pointer &enhancement1, const typename ImageType::Pointer &enhancement2, const double &timepoint1, , const double &timepoint2 );

	//template<class ImageType>
	//typename ImageType::Pointer CalculateWISFeature(const typename ImageType::Pointer &enhancement1, const typename ImageType::Pointer &enhancement2, const double &timepoint1, , const double &timepoint2);


	/**
	\brief Calculates the features for given images of one subject

	\param t1ceImagePointer			Pointer to T1CE image
	\param t2flairImagePointer		Pointer to T2-FLAIR image
	\param t1ImagePointer			Pointer to T1-weighted image
	\param t2ImagePointer			Pointer to T2 image
	\param rcbvImagePointer			Pointer to RCBV image
	\param psrImagePointer			Pointer to Percent Signal Recovery (PSR) image
	\param phImagePointer			Pointer to Peak Height (PH) image
	\param axImagePointer			Pointer to Axial Diffusivity image
	\param faImagePointer			Pointer to Fractional Anisotropy image
	\param radImagePointer			Pointer to Radial Diffusivity image
	\param trImagePointer			Pointer to Trace image
	\param labelImagePointer		Pointer to image having segmentation labels in patient space
	\param atlasImagePointer		Pointer to image having segmentation labels in atla space
	\param templateImagePointer		Pointer to atlas to calculate location features
	\param configuration			Configurations for calculation of histogram features
	*/

	template<class ImageType>
	double SurvivalEstimateOnGivenSubject(typename ImageType::Pointer LabelMap,
		typename ImageType::Pointer AtlasMap,
		typename ImageType::Pointer TemplateMap,
		typename ImageType::Pointer FinalT1CEImagePointer,
		typename ImageType::Pointer FinalT2FlairImagePointer,
		typename ImageType::Pointer FinalT1ImagePointer,
		typename ImageType::Pointer FinalT2ImagePointer,
		typename ImageType::Pointer PHImagePointer,
		typename ImageType::Pointer PSRImagePointer,
		typename ImageType::Pointer RCBVImagePointer,
		typename ImageType::Pointer AXImagePointer,
		typename ImageType::Pointer FAImagePointer,
		typename ImageType::Pointer RADImagePointer,
		typename ImageType::Pointer TRImagePointer, int age, std::string modeldirectory);


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
		const typename ImageType::Pointer &templateImagePointer,
		const VariableSizeMatrixType &configuration);

	/**
	\brief Calculates the statistcal features (mean, satndard deviation) for given vector

	\param intensities			Input vector having intensities from a region of an image
	*/
	VectorDouble GetStatisticalFeatures(const VectorDouble &intensities);



	/**
	\brief Calculates the histogram binning based features for given vector

	\param intensities			Input vector having intensities from a region of an image
	\param start				Start of the histogram bins
	\param interval				Interval between the centers of two hsitogram bins
	\param end					End of the histogram bins
	*/
	VectorDouble GetHistogramFeatures(const VectorDouble &intensities, const double &start, const double &interval, const double &end);

	/**
	\brief Calculates volumetric measures adn their ratios

	\param edemaSize			Volume of edema in terms of number of voxels
	\param tuSize				Volume of enhancing tumor in terms of number of voxels
	\param neSize				Volume of non-enahancing tumro core in terms of number of voxels
	\param totalSize			Volume of whoe brain in terms of number of voxels
	*/
	VectorDouble GetVolumetricFeatures(const double &edemaSize, const double &tuSize, const double &neSize, const double &totalSize);

	/**
	\brief Removes the connected components smaller than an area threshold from the tumor

	\param labelImagePointer		Pointer to image having segmentation labels in patient space
	*/
	template<class ImageType>
	std::vector<typename ImageType::Pointer> RevisedTumorArea(const typename ImageType::Pointer &labelImagePointer);

	/**
	\brief Calculates the distance features. 1. Distanc eof tumro from ventricles, 2. Distance of edema from ventricles

	\param edemaImage		Pointer to image having edema label
	\param tumorImage		Pointer to image having (enhancing tumor + non-enahncign tumor) label
	\param ventImage		Pointer to image having ventricles label
	*/
	template<class ImageType>
	VectorDouble GetDistanceFeatures3(const typename ImageType::Pointer &edemaImage, const typename ImageType::Pointer &tumorImage, const typename ImageType::Pointer &ventImage);

	/**
	\brief Calculates the distance map of an image by keeping the given label as reference seed point

	\param labelImagePointer		Pointer to image having segmentation labels in patient space
	\param label					Label to use as a reference seed point
	*/

	template<class ImageType>
	typename ImageType::Pointer GetDistanceMapWithLabel(const typename ImageType::Pointer &labelImagePointer, const int &label1);

	/**
	\brief Calculates the distance map of an image by keeping the non-zero voxels of image as reference seed point

	\param image					Pointer to image having reference seed points
	*/
	template<class ImageType>
	typename ImageType::Pointer GetDistanceMap(const typename ImageType::Pointer &image);


	/**
	\brief Calculates the spatial location features for given tumor images

	\param labelImagePointer		Pointer to image having segmentation labels in patient space
	\param atlasImagePointer		Pointer to image having segmentation labels in atlas space
	*/
	template<class ImageType>
	VectorDouble GetSpatialLocationFeatures(const typename ImageType::Pointer &labelImagePointer,
		const typename ImageType::Pointer &atlasImagePointer);






	/**
	\brief Prepares a new survival prediction model

	\param inputdirectory			Path to the directory having training data
	\param qualifiedsubjects		List of qualifeid subjects having all the data avaialble to train a model
	\param outputdirectory			Path to the output directory
	*/
	int PrepareNewMolecularPredictionModel(
		const std::string &inputdirectory,
		const std::vector< std::map< CAPTK::ImageModalityType, std::string > > &qualifiedsubjects,
		const std::string &outputdirectory);

	/**
	\brief Apply an exsiting model on new patients to predit survival

	\param modeldirectory			Path to the directory where model fiels are stored
	\param inputdirectory			Path to the directory having test data
	\param qualifiedsubjects		List of qualifeid subjects having all the data avaialble to apply a model
	\param outputdirectory			Path to the output directory
	*/
	VectorDouble MolecularSubtypePredictionOnExistingModel(const std::string &modeldirectory,
		const std::string &inputdirectory,
		const std::vector< std::map< CAPTK::ImageModalityType, std::string > > &qualifiedsubjects,
		const std::string &outputdirectory);



  VariableSizeMatrixType SelectModelFeatures(const VariableSizeMatrixType &ModelFeatures, const VariableLengthVectorType &selectedFeatures);
  VariableSizeMatrixType SelectModelFeatures(const VariableSizeMatrixType &ModelFeatures, const VectorDouble &selectedFeatures);
  
  template<class ImageType>
	typename ImageType::Pointer RemoveSmallerComponentsFromTumor(const typename ImageType::Pointer &etumorImage, const typename ImageType::Pointer &ncrImage);

	VariableLengthVectorType DistanceFunction(const VariableSizeMatrixType &testData, const std::string &filename, const double &rho, const double &bestg);

	VectorDouble CombineEstimates(const VectorDouble &result_proneural,
		const VectorDouble &result_neural,
		const VectorDouble &result_messenchymal,
		const VectorDouble &result_classical);


	template <class ImageType>
	typename ImageType::Pointer ReadNiftiImage(const std::string &filename);
	template<class ImageType>
	typename ImageType::Pointer RescaleImageIntensity(const typename ImageType::Pointer &image);


	void Run()
	{

	}

private:

};

template<class ImageType>
typename ImageType::Pointer MolecularSubtypePredictor::RescaleImageIntensity(const typename ImageType::Pointer &image)
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
typename ImageType::Pointer MolecularSubtypePredictor::ReadNiftiImage(const std::string &filename)
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
VectorDouble MolecularSubtypePredictor::GetDistanceFeatures3(const typename ImageType::Pointer &edemaImage, const typename ImageType::Pointer &tumorImage, const typename ImageType::Pointer &ventImage)
{
	//-----------------------create distance image----------------------------
	VectorDouble VentEdemaSum;
	VectorDouble VentTumorSum;
	//-------------------------get sum of images----------------------------------
	typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
	IteratorType tumorIt(tumorImage, tumorImage->GetLargestPossibleRegion());
	IteratorType ventIt(ventImage, ventImage->GetLargestPossibleRegion());
	IteratorType edemaIt(edemaImage, edemaImage->GetLargestPossibleRegion());

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
	for (int i = 0; i < SURVIVAL_NO_NEAR_DIST; i++)
	{
		sumVentEdema = sumVentEdema + VentEdemaSum[i];
		sumVentTumor = sumVentTumor + VentTumorSum[i];
	}
	VectorDouble DistanceFeatures;
	DistanceFeatures.push_back(sumVentEdema / SURVIVAL_NO_NEAR_DIST);
	DistanceFeatures.push_back(sumVentTumor / SURVIVAL_NO_NEAR_DIST);

	return DistanceFeatures;
}
//--------------------------------------------------------------------------------------
template<class ImageType>
typename ImageType::Pointer MolecularSubtypePredictor::GetDistanceMap(const typename ImageType::Pointer &labelImagePointer)
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
typename ImageType::Pointer MolecularSubtypePredictor::GetDistanceMapWithLabel(const typename ImageType::Pointer &labelImagePointer, const int &label1)
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
VectorDouble  MolecularSubtypePredictor::LoadTestData(const typename ImageType::Pointer &t1ceImagePointer,
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
	const typename ImageType::Pointer &templateImagePointer,
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


	//VectorDouble GlistrFeatures = GetGlistrFeatures(parametersNames);


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
		else if (imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::VENT)
			ventIndices.push_back(imIt.GetIndex());
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
	VectorDouble TumorBinsT1 = GetHistogramFeatures(tumorIntensitiesT1, configuration(0, 0), configuration(0, 1), configuration(0, 2));
	VectorDouble TumorBinsT2 = GetHistogramFeatures(tumorIntensitiesT2, configuration(1, 0), configuration(1, 1), configuration(1, 2));
	VectorDouble TumorBinsT1CE = GetHistogramFeatures(tumorIntensitiesT1CE, configuration(2, 0), configuration(2, 1), configuration(2, 2));
	VectorDouble TumorBinsT2Flair = GetHistogramFeatures(tumorIntensitiesT2Flair, configuration(3, 0), configuration(3, 1), configuration(3, 2));
	VectorDouble TumorBinsAX = GetHistogramFeatures(tumorIntensitiesAX, configuration(4, 0), configuration(4, 1), configuration(4, 2));
	VectorDouble TumorBinsFA = GetHistogramFeatures(tumorIntensitiesFA, configuration(5, 0), configuration(5, 1), configuration(5, 2));
	VectorDouble TumorBinsRAD = GetHistogramFeatures(tumorIntensitiesRAD, configuration(6, 0), configuration(6, 1), configuration(6, 2));
	VectorDouble TumorBinsTR = GetHistogramFeatures(tumorIntensitiesTR, configuration(7, 0), configuration(7, 1), configuration(7, 2));
	VectorDouble TumorBinsRCBV = GetHistogramFeatures(tumorIntensitiesRCBV, configuration(8, 0), configuration(8, 1), configuration(8, 2));
	VectorDouble TumorBinsPSR = GetHistogramFeatures(tumorIntensitiesPSR, configuration(9, 0), configuration(9, 1), configuration(9, 2));
	VectorDouble TumorBinsPH = GetHistogramFeatures(tumorIntensitiesPH, configuration(10, 0), configuration(10, 1), configuration(10, 2));

	VectorDouble edemaBinsT1 = GetHistogramFeatures(edemaIntensitiesT1, configuration(11, 0), configuration(11, 1), configuration(11, 2));
	VectorDouble edemaBinsT2 = GetHistogramFeatures(edemaIntensitiesT2, configuration(12, 0), configuration(12, 1), configuration(12, 2));
	VectorDouble edemaBinsT1CE = GetHistogramFeatures(edemaIntensitiesT1CE, configuration(13, 0), configuration(13, 1), configuration(13, 2));
	VectorDouble edemaBinsT2Flair = GetHistogramFeatures(edemaIntensitiesT2Flair, configuration(14, 0), configuration(14, 1), configuration(14, 2));
	VectorDouble edemaBinsAX = GetHistogramFeatures(edemaIntensitiesAX, configuration(15, 0), configuration(15, 1), configuration(15, 2));
	VectorDouble edemaBinsFA = GetHistogramFeatures(edemaIntensitiesFA, configuration(16, 0), configuration(16, 1), configuration(16, 2));
	VectorDouble edemaBinsRAD = GetHistogramFeatures(edemaIntensitiesRAD, configuration(17, 0), configuration(17, 1), configuration(17, 2));
	VectorDouble edemaBinsTR = GetHistogramFeatures(edemaIntensitiesTR, configuration(18, 0), configuration(18, 1), configuration(18, 2));
	VectorDouble edemaBinsRCBV = GetHistogramFeatures(edemaIntensitiesRCBV, configuration(19, 0), configuration(19, 1), configuration(19, 2));
	VectorDouble edemaBinsPSR = GetHistogramFeatures(edemaIntensitiesPSR, configuration(20, 0), configuration(20, 1), configuration(20, 2));
	VectorDouble edemaBinsPH = GetHistogramFeatures(edemaIntensitiesPH, configuration(21, 0), configuration(21, 1), configuration(21, 2));

	VectorDouble necoreBinsT1 = GetHistogramFeatures(necoreIntensitiesT1, configuration(22, 0), configuration(22, 1), configuration(22, 2));
	VectorDouble necoreBinsT2 = GetHistogramFeatures(necoreIntensitiesT2, configuration(23, 0), configuration(23, 1), configuration(23, 2));
	VectorDouble necoreBinsT1CE = GetHistogramFeatures(necoreIntensitiesT1CE, configuration(24, 0), configuration(24, 1), configuration(24, 2));
	VectorDouble necoreBinsT2Flair = GetHistogramFeatures(necoreIntensitiesT2Flair, configuration(25, 0), configuration(25, 1), configuration(25, 2));
	VectorDouble necoreBinsAX = GetHistogramFeatures(necoreIntensitiesAX, configuration(26, 0), configuration(26, 1), configuration(26, 2));
	VectorDouble necoreBinsFA = GetHistogramFeatures(necoreIntensitiesFA, configuration(27, 0), configuration(27, 1), configuration(27, 2));
	VectorDouble necoreBinsRAD = GetHistogramFeatures(necoreIntensitiesRAD, configuration(28, 0), configuration(28, 1), configuration(28, 2));
	VectorDouble necoreBinsTR = GetHistogramFeatures(necoreIntensitiesTR, configuration(29, 0), configuration(29, 1), configuration(29, 2));
	VectorDouble necoreBinsRCBV = GetHistogramFeatures(necoreIntensitiesRCBV, configuration(30, 0), configuration(30, 1), configuration(30, 2));
	VectorDouble necoreBinsPSR = GetHistogramFeatures(necoreIntensitiesPSR, configuration(31, 0), configuration(31, 1), configuration(31, 2));
	VectorDouble necoreBinsPH = GetHistogramFeatures(necoreIntensitiesPH, configuration(32, 0), configuration(32, 1), configuration(32, 2));

	VectorDouble TumorStatisticsT1 = GetStatisticalFeatures(tumorIntensitiesT1);
	VectorDouble TumorStatisticsT2 = GetStatisticalFeatures(tumorIntensitiesT2);
	VectorDouble TumorStatisticsT1CE = GetStatisticalFeatures(tumorIntensitiesT1CE);
	VectorDouble TumorStatisticsT2Flair = GetStatisticalFeatures(tumorIntensitiesT2Flair);
	VectorDouble TumorStatisticsAX = GetStatisticalFeatures(tumorIntensitiesAX);
	VectorDouble TumorStatisticsFA = GetStatisticalFeatures(tumorIntensitiesFA);
	VectorDouble TumorStatisticsRAD = GetStatisticalFeatures(tumorIntensitiesRAD);
	VectorDouble TumorStatisticsTR = GetStatisticalFeatures(tumorIntensitiesTR);
	VectorDouble TumorStatisticsRCBV = GetStatisticalFeatures(tumorIntensitiesRCBV);
	VectorDouble TumorStatisticsPSR = GetStatisticalFeatures(tumorIntensitiesPSR);
	VectorDouble TumorStatisticsPH = GetStatisticalFeatures(tumorIntensitiesPH);

	VectorDouble edemaStatisticsT1 = GetStatisticalFeatures(edemaIntensitiesT1);
	VectorDouble edemaStatisticsT2 = GetStatisticalFeatures(edemaIntensitiesT2);
	VectorDouble edemaStatisticsT1CE = GetStatisticalFeatures(edemaIntensitiesT1CE);
	VectorDouble edemaStatisticsT2Flair = GetStatisticalFeatures(edemaIntensitiesT2Flair);
	VectorDouble edemaStatisticsAX = GetStatisticalFeatures(edemaIntensitiesAX);
	VectorDouble edemaStatisticsFA = GetStatisticalFeatures(edemaIntensitiesFA);
	VectorDouble edemaStatisticsRAD = GetStatisticalFeatures(edemaIntensitiesRAD);
	VectorDouble edemaStatisticsTR = GetStatisticalFeatures(edemaIntensitiesTR);
	VectorDouble edemaStatisticsRCBV = GetStatisticalFeatures(edemaIntensitiesRCBV);
	VectorDouble edemaStatisticsPSR = GetStatisticalFeatures(edemaIntensitiesPSR);
	VectorDouble edemaStatisticsPH = GetStatisticalFeatures(edemaIntensitiesPH);

	VectorDouble necoreStatisticsT1 = GetStatisticalFeatures(necoreIntensitiesT1);
	VectorDouble necoreStatisticsT2 = GetStatisticalFeatures(necoreIntensitiesT2);
	VectorDouble necoreStatisticsT1CE = GetStatisticalFeatures(necoreIntensitiesT1CE);
	VectorDouble necoreStatisticsT2Flair = GetStatisticalFeatures(necoreIntensitiesT2Flair);
	VectorDouble necoreStatisticsAX = GetStatisticalFeatures(necoreIntensitiesAX);
	VectorDouble necoreStatisticsFA = GetStatisticalFeatures(necoreIntensitiesFA);
	VectorDouble necoreStatisticsRAD = GetStatisticalFeatures(necoreIntensitiesRAD);
	VectorDouble necoreStatisticsTR = GetStatisticalFeatures(necoreIntensitiesTR);
	VectorDouble necoreStatisticsRCBV = GetStatisticalFeatures(necoreIntensitiesRCBV);
	VectorDouble necoreStatisticsPSR = GetStatisticalFeatures(necoreIntensitiesPSR);
	VectorDouble necoreStatisticsPH = GetStatisticalFeatures(necoreIntensitiesPH);

	//-------------------------------Find connected components and get volumetric features -----------------------------
	std::vector<typename ImageType::Pointer> RevisedImages = RevisedTumorArea<ImageType>(labelImagePointer);
	typename ImageType::Pointer tumorImage = RevisedImages[0];
	typename ImageType::Pointer etumorImage = RevisedImages[1];
	typename ImageType::Pointer ncrImage = RevisedImages[2];

	etumorIndices.clear();
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
	VectorDouble VolumetricFeatures = GetVolumetricFeatures(edemaIndices.size(), etumorIndices.size(), necoreIndices.size(), brainIndices.size());
	VectorDouble spatialLocationFeatures = GetSpatialLocationFeatures<ImageType>(atlasImagePointer, templateImagePointer);

	//--------------------------------Alternate function for distance features----------------------------------------------  
	typename ImageType::Pointer  edemaDistanceMap = GetDistanceMapWithLabel<ImageType>(labelImagePointer, CAPTK::GLISTR_OUTPUT_LABELS::EDEMA);

	typename ImageType::Pointer ventImage = ImageType::New();
	ventImage->CopyInformation(labelImagePointer);
	ventImage->SetRequestedRegion(labelImagePointer->GetLargestPossibleRegion());
	ventImage->SetBufferedRegion(labelImagePointer->GetBufferedRegion());
	ventImage->Allocate();
	ventImage->FillBuffer(0);

	typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;

	IteratorType labelIt(labelImagePointer, labelImagePointer->GetLargestPossibleRegion());
	IteratorType ventIt(ventImage, ventImage->GetLargestPossibleRegion());
	labelIt.GoToBegin();
	ventIt.GoToBegin();
	while (!labelIt.IsAtEnd())
	{
		if (labelIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::VENT)
			ventIt.Set(CAPTK::VOXEL_STATUS::ON);
		++labelIt;
		++ventIt;
	}
	typename ImageType::Pointer  ventDistanceMap = GetDistanceMap<ImageType>(ventImage);
	typename ImageType::Pointer  tumorDistanceMap = GetDistanceMap<ImageType>(tumorImage);

	VectorDouble DistanceFeatures = GetDistanceFeatures3<ImageType>(edemaDistanceMap, tumorDistanceMap, ventDistanceMap);




	//copy data from vectors to one final feature vector
	VectorDouble TestFeatures;
	TestFeatures.insert(TestFeatures.end(), VolumetricFeatures.begin(), VolumetricFeatures.end());
	TestFeatures.insert(TestFeatures.end(), DistanceFeatures.begin(), DistanceFeatures.end());
	TestFeatures.insert(TestFeatures.end(), TumorStatisticsT1CE.begin(), TumorStatisticsT1CE.end());
	TestFeatures.insert(TestFeatures.end(), TumorStatisticsT1.begin(), TumorStatisticsT1.end());
	TestFeatures.insert(TestFeatures.end(), TumorStatisticsT2.begin(), TumorStatisticsT2.end());
	TestFeatures.insert(TestFeatures.end(), TumorStatisticsT2Flair.begin(), TumorStatisticsT2Flair.end());
	TestFeatures.insert(TestFeatures.end(), TumorStatisticsPH.begin(), TumorStatisticsPH.end());
	TestFeatures.insert(TestFeatures.end(), TumorStatisticsPSR.begin(), TumorStatisticsPSR.end());
	TestFeatures.insert(TestFeatures.end(), TumorStatisticsRCBV.begin(), TumorStatisticsRCBV.end());
	TestFeatures.insert(TestFeatures.end(), TumorStatisticsFA.begin(), TumorStatisticsFA.end());
	TestFeatures.insert(TestFeatures.end(), TumorStatisticsAX.begin(), TumorStatisticsAX.end());
	TestFeatures.insert(TestFeatures.end(), TumorStatisticsRAD.begin(), TumorStatisticsRAD.end());
	TestFeatures.insert(TestFeatures.end(), TumorStatisticsTR.begin(), TumorStatisticsTR.end());
	TestFeatures.insert(TestFeatures.end(), TumorBinsT1CE.begin(), TumorBinsT1CE.end());
	TestFeatures.insert(TestFeatures.end(), edemaBinsT1CE.begin(), edemaBinsT1CE.end());
	TestFeatures.insert(TestFeatures.end(), necoreBinsT1CE.begin(), necoreBinsT1CE.end());
	TestFeatures.insert(TestFeatures.end(), TumorBinsT1.begin(), TumorBinsT1.end());
	TestFeatures.insert(TestFeatures.end(), edemaBinsT1.begin(), edemaBinsT1.end());
	TestFeatures.insert(TestFeatures.end(), necoreBinsT1.begin(), necoreBinsT1.end());
	TestFeatures.insert(TestFeatures.end(), TumorBinsT2.begin(), TumorBinsT2.end());
	TestFeatures.insert(TestFeatures.end(), edemaBinsT2.begin(), edemaBinsT2.end());
	TestFeatures.insert(TestFeatures.end(), necoreBinsT2.begin(), necoreBinsT2.end());
	TestFeatures.insert(TestFeatures.end(), TumorBinsT2Flair.begin(), TumorBinsT2Flair.end());
	TestFeatures.insert(TestFeatures.end(), edemaBinsT2Flair.begin(), edemaBinsT2Flair.end());
	TestFeatures.insert(TestFeatures.end(), necoreBinsT2Flair.begin(), necoreBinsT2Flair.end());
	TestFeatures.insert(TestFeatures.end(), TumorBinsPH.begin(), TumorBinsPH.end());
	TestFeatures.insert(TestFeatures.end(), edemaBinsPH.begin(), edemaBinsPH.end());
	TestFeatures.insert(TestFeatures.end(), necoreBinsPH.begin(), necoreBinsPH.end());
	TestFeatures.insert(TestFeatures.end(), TumorBinsPSR.begin(), TumorBinsPSR.end());
	TestFeatures.insert(TestFeatures.end(), edemaBinsPSR.begin(), edemaBinsPSR.end());
	TestFeatures.insert(TestFeatures.end(), necoreBinsPSR.begin(), necoreBinsPSR.end());
	TestFeatures.insert(TestFeatures.end(), TumorBinsRCBV.begin(), TumorBinsRCBV.end());
	TestFeatures.insert(TestFeatures.end(), edemaBinsRCBV.begin(), edemaBinsRCBV.end());
	TestFeatures.insert(TestFeatures.end(), necoreBinsRCBV.begin(), necoreBinsRCBV.end());
	TestFeatures.insert(TestFeatures.end(), TumorBinsFA.begin(), TumorBinsFA.end());
	TestFeatures.insert(TestFeatures.end(), edemaBinsFA.begin(), edemaBinsFA.end());
	TestFeatures.insert(TestFeatures.end(), necoreBinsFA.begin(), necoreBinsFA.end());
	TestFeatures.insert(TestFeatures.end(), TumorBinsAX.begin(), TumorBinsAX.end());
	TestFeatures.insert(TestFeatures.end(), edemaBinsAX.begin(), edemaBinsAX.end());
	TestFeatures.insert(TestFeatures.end(), necoreBinsAX.begin(), necoreBinsAX.end());
	TestFeatures.insert(TestFeatures.end(), TumorBinsRAD.begin(), TumorBinsRAD.end());
	TestFeatures.insert(TestFeatures.end(), edemaBinsRAD.begin(), edemaBinsRAD.end());
	TestFeatures.insert(TestFeatures.end(), necoreBinsRAD.begin(), necoreBinsRAD.end());
	TestFeatures.insert(TestFeatures.end(), TumorBinsTR.begin(), TumorBinsTR.end());
	TestFeatures.insert(TestFeatures.end(), edemaBinsTR.begin(), edemaBinsTR.end());
	TestFeatures.insert(TestFeatures.end(), necoreBinsTR.begin(), necoreBinsTR.end());
	TestFeatures.insert(TestFeatures.end(), spatialLocationFeatures.begin(), spatialLocationFeatures.end());
	return TestFeatures;
}
template<class ImageType>
std::vector<typename ImageType::Pointer> MolecularSubtypePredictor::RevisedTumorArea(const typename ImageType::Pointer &labelImagePointer)
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

		if (imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::TUMOR )
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
		if (sizes[i] < SURVIVAL_SIZE_COMP)
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
typename ImageType::Pointer MolecularSubtypePredictor::RemoveSmallerComponentsFromTumor(const typename ImageType::Pointer &etumorImage, const typename ImageType::Pointer &ncrImage)
{
	typename ImageType::Pointer tumorImage = ImageType::New();
	tumorImage->CopyInformation(etumorImage);
	tumorImage->SetRequestedRegion(etumorImage->GetLargestPossibleRegion());
	tumorImage->SetBufferedRegion(etumorImage->GetBufferedRegion());
	tumorImage->Allocate();
	tumorImage->FillBuffer(0);

	typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
	IteratorType tumorIt(tumorImage, tumorImage->GetLargestPossibleRegion());
	IteratorType ncrIt(ncrImage, ncrImage->GetLargestPossibleRegion());
	IteratorType etIt(etumorImage, etumorImage->GetLargestPossibleRegion());

	tumorIt.GoToBegin();
	etIt.GoToBegin();
	ncrIt.GoToBegin();

	while (!tumorIt.IsAtEnd())
	{
		if (etIt.Get() == CAPTK::VOXEL_STATUS::ON || ncrIt.Get() == CAPTK::VOXEL_STATUS::ON)
			tumorIt.Set(CAPTK::VOXEL_STATUS::ON);
		else
			tumorIt.Set(CAPTK::VOXEL_STATUS::OFF);

		++tumorIt;
		++etIt;
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
		if (sizes[i] < SURVIVAL_SIZE_COMP)
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
	return tumorImage;
}
template<class ImageType>
VectorDouble MolecularSubtypePredictor::GetSpatialLocationFeatures(const typename ImageType::Pointer &labelImagePointer, const typename ImageType::Pointer &jacobtemplateImagePointer)
{
	std::vector<typename ImageType::Pointer> RevisedImages = RevisedTumorArea<ImageType>(labelImagePointer);
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
	//find number of voxels in 9 ROIs
	VectorDouble location;
	int tumorSize = 0;
	for (unsigned int i = 0; i < 9; i++)
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
	//find percentages in 9 ROIs
	for (unsigned int i = 0; i < 9; i++)
		location[i] = (location[i] * 100) / tumorSize;

	return location;
}
//template<class ImageType>
//double MolecularSubtypePredictor::SurvivalEstimateOnGivenSubject(typename ImageType::Pointer LabelImagePointer,
//	typename ImageType::Pointer AtlasImagePointer,
//	typename ImageType::Pointer TemplateImagePointer,
//	typename ImageType::Pointer T1CEImagePointer,
//	typename ImageType::Pointer T2FlairImagePointer,
//	typename ImageType::Pointer T1ImagePointer,
//	typename ImageType::Pointer T2ImagePointer,
//	typename ImageType::Pointer PHImagePointer,
//	typename ImageType::Pointer PSRImagePointer,
//	typename ImageType::Pointer RCBVImagePointer,
//	typename ImageType::Pointer AXImagePointer,
//	typename ImageType::Pointer FAImagePointer,
//	typename ImageType::Pointer RADImagePointer,
//	typename ImageType::Pointer TRImagePointer, int age, std::string modeldirectory)
//{
//	VariableSizeMatrixType HistogramFeaturesConfigurations;
//	HistogramFeaturesConfigurations.SetSize(33, 3); //11 modalities*3 regions = 33 configurations*3 histogram features for each configuration
//
//	CSVFileReaderType::Pointer reader = CSVFileReaderType::New();
//	VectorDouble ages;
//	MatrixType dataMatrix;
//	try
//	{
//		reader->SetFileName(modeldirectory + "/Survival_HMFeatures_Configuration.csv");
//		reader->SetFieldDelimiterCharacter(',');
//		reader->HasColumnHeadersOff();
//		reader->HasRowHeadersOff();
//		reader->Parse();
//		dataMatrix = reader->GetArray2DDataObject()->GetMatrix();
//
//		for (unsigned int i = 0; i < dataMatrix.rows(); i++)
//			for (unsigned int j = 0; j < dataMatrix.cols(); j++)
//				HistogramFeaturesConfigurations(i, j) = dataMatrix(i, j);
//	}
//	catch (const std::exception& e1)
//	{
//		logger.WriteError("Cannot find the file 'Survival_HMFeatures_Configuration.csv' in the input directory. Error code : " + std::string(e1.what()));
//		exit(EXIT_FAILURE);
//	}
//	MatrixType meanMatrix;
//	VariableLengthVectorType mean;
//	VariableLengthVectorType stddevition;
//	try
//	{
//		reader->SetFileName(modeldirectory + "/Survival_ZScore_Mean.csv");
//		reader->SetFieldDelimiterCharacter(',');
//		reader->HasColumnHeadersOff();
//		reader->HasRowHeadersOff();
//		reader->Parse();
//		meanMatrix = reader->GetArray2DDataObject()->GetMatrix();
//
//		mean.SetSize(meanMatrix.size());
//		for (unsigned int i = 0; i < meanMatrix.size(); i++)
//			mean[i] = meanMatrix(i, 0);
//	}
//	catch (const std::exception& e1)
//	{
//		logger.WriteError("Cannot find the file 'mean.csv' in the model directory. Error code : " + std::string(e1.what()));
//		exit(EXIT_FAILURE);
//	}
//	MatrixType stdMatrix;
//	try
//	{
//		reader->SetFileName(modeldirectory + "/Survival_ZScore_Std.csv");
//		reader->SetFieldDelimiterCharacter(',');
//		reader->HasColumnHeadersOff();
//		reader->HasRowHeadersOff();
//		reader->Parse();
//		stdMatrix = reader->GetArray2DDataObject()->GetMatrix();
//
//		stddevition.SetSize(stdMatrix.size());
//		for (unsigned int i = 0; i < stdMatrix.size(); i++)
//			stddevition[i] = stdMatrix(i, 0);
//	}
//	catch (const std::exception& e1)
//	{
//		logger.WriteError("Cannot find the file 'std.csv' in the model directory. Error code : " + std::string(e1.what()));
//		exit(EXIT_FAILURE);
//	}
//	//----------------------------------------------------
//	VariableSizeMatrixType FeaturesOfAllSubjects;
//	FeaturesOfAllSubjects.SetSize(1, 161);
//
//	VectorDouble TestFeatures = LoadTestData<ImageType>(T1CEImagePointer, T2FlairImagePointer, T1ImagePointer, T2ImagePointer,
//		RCBVImagePointer, PSRImagePointer, PHImagePointer, AXImagePointer, FAImagePointer, RADImagePointer, TRImagePointer, LabelImagePointer, AtlasImagePointer, TemplateImagePointer, HistogramFeaturesConfigurations);
//	FeaturesOfAllSubjects(0, 0) = age;
//	for (unsigned int i = 1; i <= TestFeatures.size(); i++)
//		FeaturesOfAllSubjects(0, i) = TestFeatures[i - 1];
//
//	VariableSizeMatrixType ScaledTestingData = mFeatureScalingLocalPtr.ScaleGivenTestingFeatures(FeaturesOfAllSubjects, mean, stddevition);
//	VariableSizeMatrixType ScaledFeatureSetAfterAddingLabel;
//	ScaledFeatureSetAfterAddingLabel.SetSize(ScaledTestingData.Rows(), ScaledTestingData.Cols() + 1);
//	for (unsigned int i = 0; i < ScaledTestingData.Rows(); i++)
//	{
//		unsigned int j = 0;
//		for (j = 0; j < ScaledTestingData.Cols(); j++)
//			ScaledFeatureSetAfterAddingLabel(i, j) = ScaledTestingData(i, j);
//		ScaledFeatureSetAfterAddingLabel(i, j) = 0;
//	}
//	VariableSizeMatrixType SixModelSelectedFeatures = SelectSixMonthsModelFeatures(ScaledFeatureSetAfterAddingLabel);
//	VariableSizeMatrixType EighteenModelSelectedFeatures = SelectEighteenMonthsModelFeatures(ScaledFeatureSetAfterAddingLabel);
//	VectorDouble results;
//	try
//	{
//		VariableLengthVectorType result_6;
//		VariableLengthVectorType result_18;
//		if (cbica::fileExists(modeldirectory + "/Survival_SVM_Model6.csv") == true && cbica::fileExists(modeldirectory + "/Survival_SVM_Model18.csv") == true)
//		{
//			result_6 = DistanceFunction(SixModelSelectedFeatures, modeldirectory + "/Survival_SVM_Model6.csv", SURVIVAL_MODEL6_RHO, SURVIVAL_MODEL6_G);
//			result_18 = DistanceFunction(EighteenModelSelectedFeatures, modeldirectory + "/Survival_SVM_Model18.csv", SURVIVAL_MODEL18_RHO, SURVIVAL_MODEL18_G);
//		}
//		else if (cbica::fileExists(modeldirectory + "/Survival_SVM_Model6.xml") == true && cbica::fileExists(modeldirectory + "/Survival_SVM_Model18.xml") == true)
//		{
//		}
//		else
//		{
//			logger.WriteError("Error caught during testing: There is no exisitg model file in the model directory: " + modeldirectory);
//			exit(EXIT_FAILURE);
//		}
//
//		//Estimates=(abs(Estimates)<=1).*Estimates+(Estimates>1)-(Estimates<-1);
//		results = CombineEstimates(result_6, result_18);
//		//std::ofstream myfile;
//		//myfile.open(outputdirectory + "/results.csv");
//		//myfile << "SubjectName,Result\n";
//
//		//for (size_t i = 0; i < results.size(); i++)
//		//{
//		//	std::map< ImageModalityType, std::string > currentsubject = qualifiedsubjects[i];
//		//	myfile << static_cast<std::string>(currentsubject[IMAGE_TYPE_SUDOID]) + "," + std::to_string(result_6[i]) + "," + std::to_string(result_18[i]) + "," + std::to_string(results[i]) + "\n";
//		//}
//		//myfile.close();
//	}
//	catch (itk::ExceptionObject & excp)
//	{
//		logger.WriteError("Error caught during testing: " + std::string(excp.GetDescription()));
//		exit(EXIT_FAILURE);
//	}
//	return results[0];
//}


//template<class ImageType>
//typename ImageType::Pointer MolecularSubtypePredictor::CalculateSERFeature(const typename ImageType::Pointer &enhancement1, const typename ImageType::Pointer &enhancement2)
//{
//	typename ImageType::Pointer localizeImage = ImageType::New();
//	localizeImage->CopyInformation(enhancement1);
//	localizeImage->SetRequestedRegion(enhancement1->GetLargestPossibleRegion());
//	localizeImage->SetBufferedRegion(enhancement1->GetBufferedRegion());
//	localizeImage->Allocate();
//	localizeImage->FillBuffer(0);

//	typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
//	IteratorType et1It(enhancement1, enhancement1->GetLargestPossibleRegion());
//	IteratorType et2It(enhancement2, enhancement2->GetLargestPossibleRegion());
//	IteratorType localizeIt(localizeImage, localizeImage->GetLargestPossibleRegion());

//	et1It.GoToBegin();
//	et2It.GoToBegin();
//	localizeIt.GoToBegin();

//	while (!et1It.IsAtEnd())
//	{
//		localizeIt.Set(et1It.Get() / et2It.Get());
//		++localizeIt;
//		++et1It;
//		++et2It;
//	}
//	return localizeImage;
//}

//template<class ImageType>
//typename ImageType::Pointer MolecularSubtypePredictor::CalculatePEFeature(const typename ImageType::Pointer &enhancement1, const typename ImageType::Pointer &enhancement2)
//{
//	typename ImageType::Pointer localizeImage = ImageType::New();
//	localizeImage->CopyInformation(enhancement1);
//	localizeImage->SetRequestedRegion(enhancement1->GetLargestPossibleRegion());
//	localizeImage->SetBufferedRegion(enhancement1->GetBufferedRegion());
//	localizeImage->Allocate();
//	localizeImage->FillBuffer(0);

//	typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
//	IteratorType et1It(enhancement1, enhancement1->GetLargestPossibleRegion());
//	IteratorType et2It(enhancement2, enhancement2->GetLargestPossibleRegion());
//	IteratorType localizeIt(localizeImage, localizeImage->GetLargestPossibleRegion());

//	et1It.GoToBegin();
//	et2It.GoToBegin();
//	localizeIt.GoToBegin();

//	while (!et1It.IsAtEnd())
//	{
//		if (et1It.Get() > et2It.Get())
//			localizeIt.Set(et1It.Get());
//		else
//			localizeIt.Set(et2It.Get());
//		++localizeIt;
//		++et1It;
//		++et2It;
//	}
//	return localizeImage;

//}

//template<class ImageType>
//typename ImageType::Pointer MolecularSubtypePredictor::CalculateWOSFeature(const typename ImageType::Pointer &enhancement1, const typename ImageType::Pointer &enhancement1, const double &timepoint1, , const double &timepoint2)
//{
//	typename ImageType::Pointer localizeImage = ImageType::New();
//	localizeImage->CopyInformation(enhancement1);
//	localizeImage->SetRequestedRegion(enhancement1->GetLargestPossibleRegion());
//	localizeImage->SetBufferedRegion(enhancement1->GetBufferedRegion());
//	localizeImage->Allocate();
//	localizeImage->FillBuffer(0);

//	typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
//	IteratorType et1It(enhancement1, enhancement1->GetLargestPossibleRegion());
//	IteratorType et2It(enhancement2, enhancement2->GetLargestPossibleRegion());
//	IteratorType localizeIt(localizeImage, localizeImage->GetLargestPossibleRegion());

//	et1It.GoToBegin();
//	et2It.GoToBegin();
//	localizeIt.GoToBegin();

//	while (!et1It.IsAtEnd())
//	{
//		double val = (et1It.Get() - et2It.Get()) / (timepoint2-timepoint1);
//			localizeIt.Set(val);
//		++localizeIt;
//		++et1It;
//		++et2It;
//	}
//	return localizeImage;
//}



//template<class ImageType>
//typename ImageType::Pointer MolecularSubtypePredictor::CalculateKineticFeatures(const typename ImageType::Pointer &enhancement1, const typename ImageType::Pointer &enhancement2, const typename ImageType::Pointer &enhancement3, const double &timepoint1, , const double &timepoint2, const double &timepoint3)
//{
//	typename ImageType::Pointer enhancement1 = ImageType::New();
//	enhancement1->CopyInformation(input1);
//	enhancement1->SetRequestedRegion(input1->GetLargestPossibleRegion());
//	enhancement1->SetBufferedRegion(input1->GetBufferedRegion());
//	enhancement1->Allocate();
//	enhancement1->FillBuffer(0);
//	typename ImageType::Pointer enhancement2 = ImageType::New();
//	enhancement2->CopyInformation(input1);
//	enhancement2->SetRequestedRegion(input1->GetLargestPossibleRegion());
//	enhancement2->SetBufferedRegion(input1->GetBufferedRegion());
//	enhancement2->Allocate();
//	enhancement2->FillBuffer(0);

//	typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
//	IteratorType et1It(enhancement1, enhancement1->GetLargestPossibleRegion());
//	IteratorType et2It(enhancement2, enhancement2->GetLargestPossibleRegion());
//	IteratorType in1It(input1, input1->GetLargestPossibleRegion());
//	IteratorType in2It(input2, input2->GetLargestPossibleRegion());
//	IteratorType in3It(input3, input3->GetLargestPossibleRegion());

//	et1It.GoToBegin();
//	et2It.GoToBegin();
//	in1It.GoToBegin();
//	in2It.GoToBegin();
//	i3It.GoToBegin();

//	while (!et1It.IsAtEnd())
//	{
//		et1It.Set(in2It.Get() - in1It.Get()) ;
//		et3It.Set(in3It.Get() - in1It.Get());
//		++in1It;
//		++in2It;
//		++in3It;

//		++et1It;
//		++et2It;
//	}
//	return enhancement1;
//}
//template<class ImageType>
//typename ImageType::Pointer MolecularSubtypePredictor::CalculateWISFeature(const typename ImageType::Pointer &enhancement1, const typename ImageType::Pointer &enhancement1, const double &timepoint1, , const double &timepoint2)
//{
//}




#endif







