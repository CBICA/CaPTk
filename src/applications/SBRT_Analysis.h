#pragma once

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <float.h>
//#include <utility>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodIterator.h"

#include "itkScalarImageToCooccurrenceMatrixFilter.h"
#include "itkScalarImageToRunLengthFeaturesFilter.h"
#include "itkHistogramToTextureFeaturesFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVectorImage.h"
#include "itkVectorContainer.h"
#include "itkStatisticsImageFilter.h"
#include "itkImageToHistogramFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkShapeLabelObject.h"
#include "itkStatisticsLabelObject.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkLabelImageToStatisticsLabelMapFilter.h"
#include "itkStatisticsRelabelLabelMapFilter.h"
#include "itkLabelStatisticsImageFilter.h"

#include "itkScalarImageToTextureFeaturesFilter.h"
#include "itkScalarImageToRunLengthFeaturesFilter.h"

#include "vnl/vnl_math.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
//#include "vnl/algo/vnl_matrix_inverse.h"

//#include "opencv2/core.hpp"

#include "cbicaLogging.h"

/**
\class SBRT_Analysis
\brief SBRT_Analysis class - the main class structure enclosing all the functions 
		for predicting risks regarding survival and nodal failure 
		for early-stage non-small cell lung cancer patients treated with SBRT based on PET images.
*/

template < class TPixelType, unsigned int TDimension >
class SBRT_Analysis
{
public:
	typedef SBRT_Analysis			Self;

	typedef itk::Image< TPixelType, TDimension >		ScalarImageType;	
	typedef typename ScalarImageType::Pointer			ScalarImageTypePointer;
	typedef unsigned short LabelPixelType;
	typedef itk::Image< LabelPixelType, TDimension >	ScalarLabelImageType;
	typedef typename ScalarLabelImageType::Pointer 		ScalarLabelImageTypePointer;
	//typedef itk::Vector< TPixelType, TInputImageNum >   VectorType;
	//typedef itk::Image< VectorType,TDimension >			VectorImageType;
	//typedef typename VectorImageType::Pointer			VectorImageTypePointer;

	typedef itk::ImageRegionIterator< ScalarImageType >			ScalarIteratorType;
	typedef itk::NeighborhoodIterator< ScalarImageType >		ScalarNeighborhoodIteratorType;
	typedef itk::ImageRegionIterator< ScalarLabelImageType >	ScalarLabelImageIteratorType;
	//typedef itk::ImageRegionIterator< VectorImageType >			VectorIteratorType;
	//typedef itk::NeighborhoodIterator< VectorImageType >		VectorNeighborhoodIteratorType;
	typedef typename ScalarNeighborhoodIteratorType::RadiusType		ScalarNeighborhoodRadiusType;
	//typedef typename VectorNeighborhoodIteratorType::RadiusType		VectorNeighborhoodRadiusType;

	typedef typename ScalarImageType::IndexType			ScalarIndexType;
	typedef typename ScalarImageType::SizeType			ScalarSizeType;
	//typedef typename VectorImageType::IndexType			VectorIndexType;
	//typedef typename VectorImageType::SizeType			VectorSizeType;

	typedef itk::ImageFileReader< ScalarImageType >		ScalarImageReaderType;
	typedef itk::ImageFileWriter< ScalarImageType >		ScalarImageWriterType;
	typedef itk::ImageFileReader< ScalarLabelImageType >		ScalarLabelImageReaderType;
	typedef itk::ImageFileWriter< ScalarLabelImageType >		ScalarLabelImageWriterType;

	typedef itk::LabelStatisticsImageFilter< ScalarImageType, ScalarLabelImageType > LabelStatisticsImageFilterType;

	typedef itk::ConnectedComponentImageFilter< ScalarLabelImageType,ScalarLabelImageType > ConnectedComponentType;
	typedef itk::RelabelComponentImageFilter< ScalarLabelImageType,ScalarLabelImageType > RelabelType;

	typedef itk::ShapeLabelObject< LabelPixelType, TDimension >		ShapeLabelObjectType;
	typedef itk::LabelMap< ShapeLabelObjectType >         	LabelMapType;
	typedef itk::LabelImageToShapeLabelMapFilter< ScalarLabelImageType,LabelMapType > LabImg2LabMapType;

	typedef itk::StatisticsLabelObject< LabelPixelType, TDimension > StatisticsLabelObjectType;
	typedef itk::LabelMap< StatisticsLabelObjectType > 		StatisticsLabelMapType;
	typedef itk::LabelImageToStatisticsLabelMapFilter< ScalarLabelImageType, ScalarImageType, StatisticsLabelMapType > LabImg2StatLabMapType;
	typedef itk::StatisticsRelabelLabelMapFilter< StatisticsLabelMapType > StatisticsRelabelMapType;

	typedef itk::Statistics::ScalarImageToTextureFeaturesFilter<ScalarImageType> GlcmType;
	typedef itk::Statistics::ScalarImageToRunLengthFeaturesFilter<ScalarImageType> GlrlmType;

private:
	SBRT_Analysis(const Self &);
	void operator=(const Self &);

	ScalarImageTypePointer m_inputImage;
	ScalarImageTypePointer m_relabScalarImage;
	ScalarLabelImageTypePointer m_labelImage;
	ScalarLabelImageTypePointer m_relabImage;
	
	ScalarSizeType m_ImageSize;
	std::vector< float > m_ImageVoxelSize;
	
	unsigned int m_roiLabel;
	int m_numRoi;

	std::vector<float> m_shape; // future work: in case there are multiple ROIs
	std::vector<float> m_intensity;
	std::vector<float> m_glcm;
	std::vector<float> m_glrlm;
	std::vector<float> m_glszm;
	std::vector<std::string> m_shape_name;
	std::vector<std::string> m_intensity_name;
	std::vector<std::string> m_glcm_name;
	std::vector<std::string> m_glrlm_name;
	std::vector<std::string> m_glszm_name;

	std::vector< std::vector<float> > m_fea_all;
	std::vector<float> m_fea_vec;
	std::vector< std::vector<std::string> > m_feaName_all;

	//vnl_matrix<float> m_feaScaling;
	//vnl_matrix<float> m_metaFeaProj;
	//vnl_matrix<float> m_riskProj;
	std::vector<float> m_predRisk;

	//string m_outputBaseName;
	cbica::Logging *m_logger;
	bool m_logger_or_not;

public:
	SBRT_Analysis();
	~SBRT_Analysis();

    /**
    \brief Initialization function reads the input images.

    \param InputImageName refers to input PET image.
    \param InputLabelName refers to nodule labe images within which radiomic features will be extracted.
    */
	int Initialization( std::string InputImageName, std::string InputLabelName );

	/**
    \brief SetParameters function assigns the label value that corresponds to the nodule region in the mask image.

    \param roiLab refers to the nodule label.
    */
	void SetParameters( unsigned int roiLab );

	void SetLogger( std::string logName );
	
	/**
    \brief FeaExtraction function encloses all the functions that calculate different texture features.
    */
	void FeaExtraction();
	void CalcShapeFeature();
	void CalcIntensityFeature();
	void CalcGlcmFeature();
	void CalcGlrlmFeature();
	
	/**
    \brief OutputFeature function saves the texture features.
    \param oname is the file name to save the features.
    */
	void OutputFeature( std::string oname );

	/**
    \brief GetPredictedRisk function estimates the risks regarding survival and nodal failure based on the texture features.
    \param metaName refers to the model for meta-feature extraction.
    \param riskProjName refers to the model for risk prediction.
    */
	void GetPredictedRisk( std::string metaName, std::string riskProjName );

  /**
    \brief Get predicted Survival risk
    \return Survival risk value as a float
    */
  float GetSurivalRisk();

  /**
    \brief Get predicted Nodal Failure
    \return Nodal failure risk as a float
    */
  float GetNodalFailureRisk();
};

#include "SBRT_Analysis.hxx"


