#pragma once

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <float.h>
//#include <utility>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodIterator.h"

#include "itkVectorGradientMagnitudeImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkImageDuplicator.h"

#include "itkBinaryThresholdImageFunction.h"
#include "itkFloodFilledImageFunctionConditionalIterator.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
//#include "itkBinaryContourImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"

#include "vnl/vnl_math.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
//#include "vnl/algo/vnl_matrix_inverse.h"

#include "opencv2/core.hpp"

#include "cbicaLogging.h"

using namespace std;

/**
\class SBRT_LungField
\brief SBRT_LungField class - the main class structure enclosing all the functions for lung field segmentation based on PET/CT images.

*/

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
class SBRT_LungField
{
public:
	typedef SBRT_LungField			Self;

	typedef itk::Image< TPixelType, TDimension >		ScalarImageType;	
	typedef typename ScalarImageType::Pointer			ScalarImageTypePointer;
	typedef int 			LabelPixelType;
	typedef itk::Image< LabelPixelType, TDimension>				ScalarLabelImageType;
	typedef typename ScalarLabelImageType::Pointer 		ScalarLabelImageTypePointer;
	typedef itk::Vector< TPixelType, TInputImageNum >   VectorType;
	typedef itk::Image< VectorType,TDimension >			VectorImageType;
	typedef typename VectorImageType::Pointer			VectorImageTypePointer;

	typedef itk::ImageRegionIterator< ScalarImageType >			ScalarIteratorType;
	typedef itk::NeighborhoodIterator< ScalarImageType >		ScalarNeighborhoodIteratorType;
	typedef itk::ImageRegionIterator< ScalarLabelImageType >	ScalarLabelImageIteratorType;
	typedef itk::ImageRegionIterator< VectorImageType >			VectorIteratorType;
	typedef itk::NeighborhoodIterator< VectorImageType >		VectorNeighborhoodIteratorType;
	typedef typename ScalarNeighborhoodIteratorType::RadiusType		ScalarNeighborhoodRadiusType;
	typedef typename VectorNeighborhoodIteratorType::RadiusType		VectorNeighborhoodRadiusType;

	typedef typename ScalarImageType::IndexType			ScalarIndexType;
	typedef typename ScalarImageType::SizeType			ScalarSizeType;
	typedef typename VectorImageType::IndexType			VectorIndexType;
	typedef typename VectorImageType::SizeType			VectorSizeType;

	typedef itk::ImageFileReader< ScalarImageType >		ScalarImageReaderType;
	typedef itk::ImageFileWriter< ScalarImageType >		ScalarImageWriterType;
	typedef itk::ImageFileWriter< ScalarLabelImageType >		ScalarLabelImageWriterType;

	typedef itk::ExtractImageFilter< VectorImageType,VectorImageType > ExtractImageFilterType;
	typedef itk::VectorGradientMagnitudeImageFilter<VectorImageType,float> VectorGraMagFilterType;
	typedef itk::GradientMagnitudeImageFilter<ScalarImageType,ScalarImageType> ScalarGraMagFilterType;
	typedef itk::ImageDuplicator< ScalarImageType > DuplicatorType;

	typedef itk::LabelStatisticsImageFilter< ScalarImageType, ScalarLabelImageType > LabelStatisticsImageFilterType;

	typedef itk::BinaryThresholdImageFunction<ScalarImageType,float> BinaryFunctionType;
	typedef itk::FloodFilledImageFunctionConditionalIterator<ScalarImageType,BinaryFunctionType> FloodFilledIterator;

	typedef itk::BinaryThresholdImageFunction<ScalarLabelImageType,float> LabelBinaryFunctionType;
	typedef itk::FloodFilledImageFunctionConditionalIterator<ScalarLabelImageType,BinaryFunctionType> LabelFloodFilledIterator;

	typedef itk::ConnectedComponentImageFilter< ScalarLabelImageType,ScalarLabelImageType > ConnectedComponentType;
	typedef itk::RelabelComponentImageFilter< ScalarLabelImageType,ScalarLabelImageType > RelabelType;

private:
	SBRT_LungField(const Self &);
	void operator=(const Self &);

	VectorImageTypePointer m_InputImageVector;
	ScalarImageTypePointer m_Mask;
	ScalarImageTypePointer m_GraMagnitude;
	ScalarImageTypePointer m_OutputImage;
	ScalarLabelImageTypePointer m_labelingImage;
	ScalarLabelImageTypePointer m_KmeanLabelingImage;
	ScalarImageTypePointer m_voxelID;
	VectorNeighborhoodRadiusType m_VectorRadius;
	ScalarNeighborhoodRadiusType m_ScalarRadius;
	VectorSizeType m_ImageSize;
	vector< float > m_ImageVoxelSize;
	vnl_matrix< float > m_ImageFea;
	vnl_matrix< float > m_SvxFea;

  int m_maskAv;
	float m_compactness;
	float m_svxSize;
	float m_minSize;
	float m_step;
	int m_iterNum;
	int m_labNum;
	int m_validVxNum;
	int m_ptSz;
	int m_feaNum;
	vector<vector<float> > m_SeedCentroids;

	//string outputBaseName;
	cbica::Logging *m_logger;
	bool m_logger_or_not;

public:
	SBRT_LungField();
	~SBRT_LungField();
    
    /**
    \brief Initialization function reads and pre-processes the input image, and initialize the output images.

    \param InputImageName contains the names of input images.
    */
	int Initialization( vector<string> &InputImageName);
	
    /**
    \brief ReadMask function loads foreground mask image if provided.
    */
	void ReadMask( string maskName );

	/**
    \brief GetImageMask function generates foreground mask image automatically if it is not provided.
    */
	void GetImageMask();

	/**
    \brief SetParameters function sets parameters involved in the whole class.
    */
	void SetParameters(int maskAv, float svxSize,float compactness,float minSize, int iter, int ptSz );

	void SetLogger(string logName);
	
	/**
    \brief GenerateSeeds function generates centroids for the SLIC oversegmentation.
    */
	void GenerateSeeds();

	/**
    \brief FeaExtraction function extracts image features for each voxel, which are used for the SLIC oversegmentation.
    */
	void FeaExtraction();

	/**
    \brief PerformSLIC function performs the SLIC oversegmentation.
    \param compactness is the trade-off parameter for balancing the intensity distance and spatial distance for the SLIC oversegmentation.
    */
	void PerformSLIC(  float compactness );

	/**
    \brief EnforceSupervoxelLabelConnectivity function post-processes the initial SLIC oversegmentation so that isolated regions with the same label are eliminated.
    */
	void EnforceSupervoxelLabelConnectivity();

	/**
    \brief DoSupervoxelSegmentation function encloses all the fuctions than contribute to SLIC oversegmentation.
    */
	void DoSupervoxelSegmentation();

	/**
    \brief GetSvxlFea function extracts features for each super-voxel within the input image "imgName".
    */
	void GetSvxlFea( string imgName );

	/**
    \brief DoKmeans function performs K-means clustering at the super-voxel level.
    \param numK is the number of clusters
    \param refImg helps to identify the lung field region after K-means clustering.
    */
    void DoKmeans( int numK, string refImg );

    /**
    \brief WriteLabel function saves the lung field segmentation results.
    */
    void WriteLabel(string outBaseName);// , int maskAv );
};

#include "SBRT_LungField.hxx"
