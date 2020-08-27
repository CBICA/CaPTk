/**
\file  ImagingSubtypePredictor.h

\brief The header file containing the SurvivaPredictor class, used to find Survival Prediction Index
Author: Saima Rathore
Library Dependecies: ITK 4.7+ <br>

https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2016 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/captk/license.html
*/

#ifndef _ImagingSubtypePredictor_h_
#define _ImagingSubtypePredictor_h_

//#include "CAPTk.h"
#include "NiftiDataManager.h"
#include "FeatureReductionClass.h"
#include "FeatureScalingClass.h"
#include "FeatureExtractionClass.h"
#include "cbicaLogging.h"
#include "cbicaITKUtilities.h"
#include "itkScalarImageToCooccurrenceMatrixFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkHistogramToTextureFeaturesFilter.h"

#ifdef APP_BASE_CAPTK_H
#include "ApplicationBase.h"
#endif

typedef itk::Image< float, 3 > ImageType;
typedef itk::CSVArray2DFileReader<double> CSVFileReaderType;
typedef vnl_matrix<double> MatrixType;
typedef ImageType::OffsetType OffsetType;
using OffsetVector = itk::VectorContainer < unsigned char, OffsetType >;

typedef itk::Statistics::ScalarImageToCooccurrenceMatrixFilter<ImageTypeFloat3D>Image2CoOccuranceType;
typedef Image2CoOccuranceType::HistogramType HistogramType;
typedef itk::Statistics::HistogramToTextureFeaturesFilter<HistogramType> Hist2FeaturesType1;

// pre-calculated values
#define SURVIVAL_MODEL6_RHO		-1.0927
#define SURVIVAL_MODEL6_G		0.0313
#define SURVIVAL_MODEL18_RHO	-0.2854
#define SURVIVAL_MODEL18_G		0.5

#define SURVIVAL_SIZE_COMP 100
#define SURVIVAL_NO_NEAR_DIST 100

/**
\class ImagingSubtypePredictor

\brief Calculates Survival Prediction Index

Reference:

@article{macyszyn2015imaging,
title={Imaging patterns predict patient survival and Imaging subtype in glioblastoma via machine learning techniques},
author={Macyszyn, Luke and Akbari, Hamed and Pisapia, Jared M and Da, Xiao and Attiah, Mark and Pigrish, Vadim and Bi, Yingtao and Pal, Sharmistha and Davuluri, Ramana V and Roccograndi, Laura and others},
journal={Neuro-oncology},
volume={18},
number={3},
pages={417--425},
year={2015},
publisher={Society for Neuro-Oncology}
}

*/
class ImagingSubtypePredictor
#ifdef APP_BASE_CAPTK_H
	: public ApplicationBase
#endif
{
public:
	//! Default constructor
	ImagingSubtypePredictor()
	{
		mSixTrainedFile = "SixModelFile.xml";
		mEighteenTrainedFile = "EighteenModelFile.xml";
	}

	//! Default destructor
	~ImagingSubtypePredictor() {};

	std::string mEighteenTrainedFile, mSixTrainedFile;
	FeatureExtractionClass mFeatureExtractionLocalPtr;
	FeatureScalingClass mFeatureScalingLocalPtr;
	cbica::Logging logger;

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
	VectorDouble GetTextureFeatures(const typename ImageType::Pointer &t1ceImagePointer, const typename ImageType::Pointer &labelImagePointer);

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
		const typename ImageType::Pointer &labelImagePointer);

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
	VectorDouble GetHistogramFeatures(const VectorDouble &intensities);

	/**
	\brief Calculates volumetric measures adn their ratios

	\param edemaSize			Volume of edema in terms of number of voxels
	\param tuSize				Volume of enhancing tumor in terms of number of voxels
	\param neSize				Volume of non-enahancing tumro core in terms of number of voxels
	\param totalSize			Volume of whoe brain in terms of number of voxels
	*/
	VectorDouble GetVolumetricFeatures(const double &edemaSize, const double &tuSize, const double &neSize, const double &totalSize);


	/**
	\brief Apply an exsiting model on new patients to predit survival

	\param modeldirectory			Path to the directory where model fiels are stored
	\param inputdirectory			Path to the directory having test data
	\param qualifiedsubjects		List of qualifeid subjects having all the data avaialble to apply a model
	\param outputdirectory			Path to the output directory
	*/
	template<class ImageType>
	double SubtypePredictionOnExistingModel(const typename ImageType::Pointer &t1ceImagePointer,
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
		const typename ImageType::Pointer &labelImagePointer);

	template <class ImageType>
	typename ImageType::Pointer ReadNiftiImage(const std::string &filename);
	template<class ImageType>
	typename ImageType::Pointer RescaleImageIntensity(const typename ImageType::Pointer &image);

	VectorDouble ScalingZeroToOne(const VectorDouble &data, const VariableLengthVectorType &minimum, const VariableLengthVectorType &maximum);

	VectorDouble CalculateEuclideanDistance(const VectorDouble &data, const VariableLengthVectorType &mean1, const VariableLengthVectorType &mean2, const VariableLengthVectorType &mean3);

	template <class TImageType = ImageTypeFloat3D, class Offsets >
	VectorDouble calculateGLCM(const typename TImageType::Pointer image, const typename TImageType::Pointer mask,
		typename Offsets::Pointer offset, const std::string &offset_select);


	void Run()
	{

	}

private:

};

template<class ImageType>
typename ImageType::Pointer ImagingSubtypePredictor::RescaleImageIntensity(const typename ImageType::Pointer &image)
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
typename ImageType::Pointer ImagingSubtypePredictor::ReadNiftiImage(const std::string &filename)
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
VectorDouble  ImagingSubtypePredictor::LoadTestData(const typename ImageType::Pointer &t1ceImagePointer,
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
	const typename ImageType::Pointer &labelImagePointer)
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
	VectorDouble TumorBinsT1 = GetHistogramFeatures(tumorIntensitiesT1);
	VectorDouble TumorBinsT2 = GetHistogramFeatures(tumorIntensitiesT2);
	VectorDouble TumorBinsT1CE = GetHistogramFeatures(tumorIntensitiesT1CE);
	VectorDouble TumorBinsT2Flair = GetHistogramFeatures(tumorIntensitiesT2Flair);
	VectorDouble TumorBinsAX = GetHistogramFeatures(tumorIntensitiesAX);
	VectorDouble TumorBinsFA = GetHistogramFeatures(tumorIntensitiesFA);
	VectorDouble TumorBinsRAD = GetHistogramFeatures(tumorIntensitiesRAD);
	VectorDouble TumorBinsTR = GetHistogramFeatures(tumorIntensitiesTR);
	VectorDouble TumorBinsRCBV = GetHistogramFeatures(tumorIntensitiesRCBV);
	VectorDouble TumorBinsPSR = GetHistogramFeatures(tumorIntensitiesPSR);
	VectorDouble TumorBinsPH = GetHistogramFeatures(tumorIntensitiesPH);

	VectorDouble edemaBinsT1 = GetHistogramFeatures(edemaIntensitiesT1);
	VectorDouble edemaBinsT2 = GetHistogramFeatures(edemaIntensitiesT2);
	VectorDouble edemaBinsT1CE = GetHistogramFeatures(edemaIntensitiesT1CE);
	VectorDouble edemaBinsT2Flair = GetHistogramFeatures(edemaIntensitiesT2Flair);
	VectorDouble edemaBinsAX = GetHistogramFeatures(edemaIntensitiesAX);
	VectorDouble edemaBinsFA = GetHistogramFeatures(edemaIntensitiesFA);
	VectorDouble edemaBinsRAD = GetHistogramFeatures(edemaIntensitiesRAD);
	VectorDouble edemaBinsTR = GetHistogramFeatures(edemaIntensitiesTR);
	VectorDouble edemaBinsRCBV = GetHistogramFeatures(edemaIntensitiesRCBV);
	VectorDouble edemaBinsPSR = GetHistogramFeatures(edemaIntensitiesPSR);
	VectorDouble edemaBinsPH = GetHistogramFeatures(edemaIntensitiesPH);

	VectorDouble necoreBinsT1 = GetHistogramFeatures(necoreIntensitiesT1);
	VectorDouble necoreBinsT2 = GetHistogramFeatures(necoreIntensitiesT2);
	VectorDouble necoreBinsT1CE = GetHistogramFeatures(necoreIntensitiesT1CE);
	VectorDouble necoreBinsT2Flair = GetHistogramFeatures(necoreIntensitiesT2Flair);
	VectorDouble necoreBinsAX = GetHistogramFeatures(necoreIntensitiesAX);
	VectorDouble necoreBinsFA = GetHistogramFeatures(necoreIntensitiesFA);
	VectorDouble necoreBinsRAD = GetHistogramFeatures(necoreIntensitiesRAD);
	VectorDouble necoreBinsTR = GetHistogramFeatures(necoreIntensitiesTR);
	VectorDouble necoreBinsRCBV = GetHistogramFeatures(necoreIntensitiesRCBV);
	VectorDouble necoreBinsPSR = GetHistogramFeatures(necoreIntensitiesPSR);
	VectorDouble necoreBinsPH = GetHistogramFeatures(necoreIntensitiesPH);


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

	//-------------------------------Find volumetric features -----------------------------
	VectorDouble VolumetricFeatures = GetVolumetricFeatures(edemaIndices.size(), etumorIndices.size(), necoreIndices.size(), brainIndices.size());

	//----------------------------------calculate texture features-----------------------------
	VectorDouble necoreTextureT1CE = GetTextureFeatures<ImageType>(t1ceImagePointer,labelImagePointer);

	//-------------------------------------------------------------------------------------
	//copy data from vectors to one final feature vector
	VectorDouble TestFeatures;
	
	TestFeatures.push_back(edemaStatisticsT1[0]);
	TestFeatures.push_back(edemaStatisticsT1CE[0]);
	TestFeatures.push_back(edemaStatisticsT2[0]);
	TestFeatures.push_back(edemaStatisticsT2Flair[0]);
	TestFeatures.push_back(edemaStatisticsAX[0]);
	TestFeatures.push_back(edemaStatisticsFA[0]);
	TestFeatures.push_back(edemaStatisticsRAD[0]);
	TestFeatures.push_back(edemaStatisticsTR[0]);
	TestFeatures.push_back(edemaStatisticsPH[0]);
	TestFeatures.push_back(edemaStatisticsRCBV[0]);
	TestFeatures.push_back(edemaStatisticsPSR[0]);

	TestFeatures.push_back(necoreStatisticsT1[0]);
	TestFeatures.push_back(necoreStatisticsT1CE[0]);
	TestFeatures.push_back(necoreStatisticsT2[0]);
	TestFeatures.push_back(necoreStatisticsT2Flair[0]);
	TestFeatures.push_back(necoreStatisticsAX[0]);
	TestFeatures.push_back(necoreStatisticsFA[0]);
	TestFeatures.push_back(necoreStatisticsRAD[0]);
	TestFeatures.push_back(necoreStatisticsTR[0]);
	TestFeatures.push_back(necoreStatisticsPH[0]);
	TestFeatures.push_back(necoreStatisticsRCBV[0]);
	TestFeatures.push_back(necoreStatisticsPSR[0]);

	TestFeatures.push_back(TumorStatisticsT1[0]);
	TestFeatures.push_back(TumorStatisticsT1CE[0]);
	TestFeatures.push_back(TumorStatisticsT2[0]);
	TestFeatures.push_back(TumorStatisticsT2Flair[0]);
	TestFeatures.push_back(TumorStatisticsAX[0]);
	TestFeatures.push_back(TumorStatisticsFA[0]);
	TestFeatures.push_back(TumorStatisticsRAD[0]);
	TestFeatures.push_back(TumorStatisticsTR[0]);
	TestFeatures.push_back(TumorStatisticsPH[0]);
	TestFeatures.push_back(TumorStatisticsRCBV[0]);
	TestFeatures.push_back(TumorStatisticsPSR[0]);

	TestFeatures.push_back(edemaStatisticsT1[1]);
	TestFeatures.push_back(edemaStatisticsT1CE[1]);
	TestFeatures.push_back(edemaStatisticsT2[1]);
	TestFeatures.push_back(edemaStatisticsT2Flair[1]);
	TestFeatures.push_back(edemaStatisticsAX[1]);
	TestFeatures.push_back(edemaStatisticsFA[1]);
	TestFeatures.push_back(edemaStatisticsRAD[1]);
	TestFeatures.push_back(edemaStatisticsTR[1]);
	TestFeatures.push_back(edemaStatisticsPH[1]);
	TestFeatures.push_back(edemaStatisticsRCBV[1]);
	TestFeatures.push_back(edemaStatisticsPSR[1]);

	TestFeatures.push_back(necoreStatisticsT1[1]);
	TestFeatures.push_back(necoreStatisticsT1CE[1]);
	TestFeatures.push_back(necoreStatisticsT2[1]);
	TestFeatures.push_back(necoreStatisticsT2Flair[1]);
	TestFeatures.push_back(necoreStatisticsAX[1]);
	TestFeatures.push_back(necoreStatisticsFA[1]);
	TestFeatures.push_back(necoreStatisticsRAD[1]);
	TestFeatures.push_back(necoreStatisticsTR[1]);
	TestFeatures.push_back(necoreStatisticsPH[1]);
	TestFeatures.push_back(necoreStatisticsRCBV[1]);
	TestFeatures.push_back(necoreStatisticsPSR[1]);

	TestFeatures.push_back(TumorStatisticsT1[1]);
	TestFeatures.push_back(TumorStatisticsT1CE[1]);
	TestFeatures.push_back(TumorStatisticsT2[1]);
	TestFeatures.push_back(TumorStatisticsT2Flair[1]);
	TestFeatures.push_back(TumorStatisticsAX[1]);
	TestFeatures.push_back(TumorStatisticsFA[1]);
	TestFeatures.push_back(TumorStatisticsRAD[1]);
	TestFeatures.push_back(TumorStatisticsTR[1]);
	TestFeatures.push_back(TumorStatisticsPH[1]);
	TestFeatures.push_back(TumorStatisticsRCBV[1]);
	TestFeatures.push_back(TumorStatisticsPSR[1]);




	TestFeatures.insert(TestFeatures.end(), edemaBinsT1.begin(), edemaBinsT1.end());
	TestFeatures.insert(TestFeatures.end(), necoreBinsT1.begin(), necoreBinsT1.end());
	TestFeatures.insert(TestFeatures.end(), TumorBinsT1.begin(), TumorBinsT1.end());

	TestFeatures.insert(TestFeatures.end(), edemaBinsT1CE.begin(), edemaBinsT1CE.end());
	TestFeatures.insert(TestFeatures.end(), necoreBinsT1CE.begin(), necoreBinsT1CE.end());
	TestFeatures.insert(TestFeatures.end(), TumorBinsT1CE.begin(), TumorBinsT1CE.end());

	TestFeatures.insert(TestFeatures.end(), edemaBinsT2.begin(), edemaBinsT2.end());
	TestFeatures.insert(TestFeatures.end(), necoreBinsT2.begin(), necoreBinsT2.end());
	TestFeatures.insert(TestFeatures.end(), TumorBinsT2.begin(), TumorBinsT2.end());

	TestFeatures.insert(TestFeatures.end(), edemaBinsT2Flair.begin(), edemaBinsT2Flair.end());
	TestFeatures.insert(TestFeatures.end(), necoreBinsT2Flair.begin(), necoreBinsT2Flair.end());
	TestFeatures.insert(TestFeatures.end(), TumorBinsT2Flair.begin(), TumorBinsT2Flair.end());

	TestFeatures.insert(TestFeatures.end(), edemaBinsAX.begin(), edemaBinsAX.end());
	TestFeatures.insert(TestFeatures.end(), necoreBinsAX.begin(), necoreBinsAX.end());
	TestFeatures.insert(TestFeatures.end(), TumorBinsAX.begin(), TumorBinsAX.end());

	TestFeatures.insert(TestFeatures.end(), edemaBinsFA.begin(), edemaBinsFA.end());
	TestFeatures.insert(TestFeatures.end(), necoreBinsFA.begin(), necoreBinsFA.end());
	TestFeatures.insert(TestFeatures.end(), TumorBinsFA.begin(), TumorBinsFA.end());

	TestFeatures.insert(TestFeatures.end(), edemaBinsRAD.begin(), edemaBinsRAD.end());
	TestFeatures.insert(TestFeatures.end(), necoreBinsRAD.begin(), necoreBinsRAD.end());
	TestFeatures.insert(TestFeatures.end(), TumorBinsRAD.begin(), TumorBinsRAD.end());

	TestFeatures.insert(TestFeatures.end(), edemaBinsTR.begin(), edemaBinsTR.end());
	TestFeatures.insert(TestFeatures.end(), necoreBinsTR.begin(), necoreBinsTR.end());
	TestFeatures.insert(TestFeatures.end(), TumorBinsTR.begin(), TumorBinsTR.end());


	TestFeatures.insert(TestFeatures.end(), edemaBinsPH.begin(), edemaBinsPH.end());
	TestFeatures.insert(TestFeatures.end(), necoreBinsPH.begin(), necoreBinsPH.end());
	TestFeatures.insert(TestFeatures.end(), TumorBinsPH.begin(), TumorBinsPH.end());


	TestFeatures.insert(TestFeatures.end(), edemaBinsRCBV.begin(), edemaBinsRCBV.end());
	TestFeatures.insert(TestFeatures.end(), necoreBinsRCBV.begin(), necoreBinsRCBV.end());
	TestFeatures.insert(TestFeatures.end(), TumorBinsRCBV.begin(), TumorBinsRCBV.end());

	TestFeatures.insert(TestFeatures.end(), edemaBinsPSR.begin(), edemaBinsPSR.end());
	TestFeatures.insert(TestFeatures.end(), necoreBinsPSR.begin(), necoreBinsPSR.end());
	TestFeatures.insert(TestFeatures.end(), TumorBinsPSR.begin(), TumorBinsPSR.end());

	return TestFeatures;
}


template <class ImageType>
double ImagingSubtypePredictor::SubtypePredictionOnExistingModel(const typename ImageType::Pointer &t1ceImagePointer,
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
	const typename ImageType::Pointer &labelImagePointer)
{
	std::string modelFileName = "";
	MatrixType meanMatrix;
	VariableLengthVectorType minimum;
	VariableLengthVectorType maximum;
	VariableLengthVectorType mean1;
	VariableLengthVectorType mean2;
	VariableLengthVectorType mean3;
	CSVFileReaderType::Pointer reader = CSVFileReaderType::New();
	try
	{
		reader->SetFileName(modelFileName);
		reader->SetFieldDelimiterCharacter(',');
		reader->HasColumnHeadersOff();
		reader->HasRowHeadersOff();
		reader->Parse();
		meanMatrix = reader->GetArray2DDataObject()->GetMatrix();

		minimum.SetSize(meanMatrix.rows());
		maximum.SetSize(meanMatrix.rows());
		mean1.SetSize(meanMatrix.rows());
		mean2.SetSize(meanMatrix.rows());
		mean3.SetSize(meanMatrix.rows());
		for (unsigned int i = 0; i < meanMatrix.size(); i++)
		{
			minimum[i] = meanMatrix(i, 0);
			maximum[i] = meanMatrix(i, 1);
			mean1[i] = meanMatrix(i, 2);
			mean2[i] = meanMatrix(i, 3);
			mean3[i] = meanMatrix(i, 4);
		}
	}
	catch (const std::exception& e1)
	{
		logger.WriteError("Cannot find the file 'mean.csv' in the model directory. Error code : " + std::string(e1.what()));
	}
	//----------------------------------------------------
	VectorDouble TestFeatures = LoadTestData<ImageType>(t1ceImagePointer, t2flairImagePointer, t1ImagePointer, t2ImagePointer,
		rcbvImagePointer, psrImagePointer, phImagePointer, axImagePointer, faImagePointer, radImagePointer, trImagePointer, labelImagePointer);

	TestFeatures = ScalingZeroToOne(TestFeatures,minimum,maximum);
	VectorDouble distances = CalculateEuclideanDistance(TestFeatures, mean1, mean2, mean3);
	double subtype=0.0;
	if (distances[0] < distances[1] && distances[0] < distances[2])
		subtype = 1;
	else if (distances[1] < distances[0] && distances[1] < distances[2])
		subtype = 2;
	else if (distances[2] < distances[0] && distances[2] < distances[1])
		subtype = 3;
	else
		subtype = 0;
	return subtype;

}
template<class ImageType>
VectorDouble ImagingSubtypePredictor::GetTextureFeatures(const typename ImageType::Pointer &t1ceImagePointer, const typename ImageType::Pointer &labelImagePointer)
{
	std::vector<typename ImageType::Pointer> RevisedImages;

	typename ImageType::Pointer edemaImage = ImageType::New();
	edemaImage->CopyInformation(labelImagePointer);
	edemaImage->SetRequestedRegion(labelImagePointer->GetLargestPossibleRegion());
	edemaImage->SetBufferedRegion(labelImagePointer->GetBufferedRegion());
	edemaImage->Allocate();
	edemaImage->FillBuffer(0);

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
	IteratorType edemaIt(edemaImage, edemaImage->GetLargestPossibleRegion());
	imIt.GoToBegin();
	edemaIt.GoToBegin();
	etIt.GoToBegin();
	ncrIt.GoToBegin();

	while (!etIt.IsAtEnd())
	{
		if (imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::EDEMA)
			edemaIt.Set(CAPTK::VOXEL_STATUS::ON);
		else
			edemaIt.Set(CAPTK::VOXEL_STATUS::OFF);

		if (imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::TUMOR)
			etIt.Set(CAPTK::VOXEL_STATUS::ON);
		else
			etIt.Set(CAPTK::VOXEL_STATUS::OFF);

		if (imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::NONENHANCING)
			ncrIt.Set(CAPTK::VOXEL_STATUS::ON);
		else
			ncrIt.Set(CAPTK::VOXEL_STATUS::OFF);

		++edemaIt;
		++etIt;
		++imIt;
		++ncrIt;
	}
	std::string offset_select = "Average";
	OffsetType offset;
	typedef itk::Neighborhood< float, 3 > NeighborhoodType;
	NeighborhoodType neighborhood;
	neighborhood.SetRadius(1);
	unsigned int centerIndex = neighborhood.GetCenterNeighborhoodIndex();
	auto offsets = OffsetVector::New();
	for (unsigned int d = 0; d < centerIndex; d++)
	{
		offset = neighborhood.GetOffset(d);
		offsets->push_back(neighborhood.GetOffset(d));
	}


	VectorDouble ncrFeatures = calculateGLCM<ImageType, OffsetVector>(t1ceImagePointer, ncrImage, offsets, offset_select);
	VectorDouble edemaFeatures = calculateGLCM<ImageType,OffsetVector>(t1ceImagePointer, edemaImage, offsets, offset_select);
	VectorDouble tuFeatures = calculateGLCM<ImageType, OffsetVector>(t1ceImagePointer, etumorImage, offsets, offset_select);

	VectorDouble features;
	return features;
}
template <class TImageType, typename Offsets >
VectorDouble ImagingSubtypePredictor::calculateGLCM(const typename TImageType::Pointer image, const typename TImageType::Pointer mask,
	typename Offsets::Pointer offset, const std::string &offset_select)
{
	VectorDouble featurevec;
	typename Offsets::Pointer offsets = Offsets::New();
	offsets = offset;

	typedef itk::ImageRegionIteratorWithIndex< TImageType > IteratorType;
	std::vector< typename TImageType::IndexType> index_vec;
	std::vector<double> nonzero_pix;
	IteratorType inputIt(mask, mask->GetLargestPossibleRegion());

	for (inputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt)
	{
		if (mask->GetPixel(inputIt.GetIndex()) != 0)
		{
			index_vec.push_back(inputIt.GetIndex());
			nonzero_pix.push_back(image->GetPixel(inputIt.GetIndex()));
		}
	}

	if (offset_select == "Average")
	{
		double contrast = 0, correl = 0, ener = 0, entro = 0, homo = 0, clustershade = 0, clusterprominance = 0, autocorr = 0;
		for (size_t i = 0; i < offsets->size(); i++)
		{
			Image2CoOccuranceType::Pointer glcmGenerator = Image2CoOccuranceType::New();
			glcmGenerator->SetOffset(offsets->at(i));
			glcmGenerator->SetNumberOfBinsPerAxis(16); //reasonable number of bins
			glcmGenerator->SetPixelValueMinMax(0,255);
			glcmGenerator->SetMaskImage(mask);
			glcmGenerator->SetInput(image);
			glcmGenerator->Update();
			auto featureCalc = Hist2FeaturesType1::New();
			featureCalc->SetInput(glcmGenerator->GetOutput());
			featureCalc->Update();

			contrast = contrast + featureCalc->GetFeature(Hist2FeaturesType1::Inertia);
			correl = correl + featureCalc->GetFeature(Hist2FeaturesType1::Correlation);
			ener = ener + featureCalc->GetFeature(Hist2FeaturesType1::Energy);
			homo = homo + featureCalc->GetFeature(Hist2FeaturesType1::InverseDifferenceMoment);
			entro = entro + featureCalc->GetFeature(Hist2FeaturesType1::Entropy);
			clustershade = clustershade + featureCalc->GetFeature(Hist2FeaturesType1::ClusterShade);
			clusterprominance = clusterprominance + featureCalc->GetFeature(Hist2FeaturesType1::ClusterProminence);
			autocorr = autocorr + featureCalc->GetFeature(Hist2FeaturesType1::HaralickCorrelation);
		}

		// calculate mean and variance
		double sum = std::accumulate(nonzero_pix.begin(), nonzero_pix.end(), 0.0);
		double pixelMean = sum / nonzero_pix.size();
		double var = 0;
		for (size_t x = 0; x < nonzero_pix.size(); x++)
		{
			var += pow((nonzero_pix[x] - pixelMean), 2);
		}
		double pixelVariance = (var / nonzero_pix.size());
		contrast = contrast / offsets->size();
		correl = correl / offsets->size();
		ener = ener / offsets->size();
		homo = homo / offsets->size();
		entro = entro / offsets->size();
		
		featurevec.push_back(ener);
		featurevec.push_back(contrast);
		featurevec.push_back(entro);
		featurevec.push_back(homo); 
		featurevec.push_back(correl);
		featurevec.push_back(pixelMean);
		featurevec.push_back(pixelVariance);
	}
	return featurevec;
}


#endif







