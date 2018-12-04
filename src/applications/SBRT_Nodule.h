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
#include "itkMedianImageFilter.h"
#include "itkImageDuplicator.h"

#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryContourImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"

#include "vnl/vnl_math.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_sparse_matrix.h"
#include "vnl/vnl_vector.h"
#include "vnl/algo/vnl_matrix_inverse.h"
#include "vnl/vnl_random.h"

//#include "opencv2/core.hpp"

#include "cbicaLogging.h"

/**
\class SBRT_Nodule
\brief SBRT_Nodule class - the main class structure enclosing all the functions for lung nodule segmentation based on PET/CT images.

*/
template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
class SBRT_Nodule
{
public:
	typedef SBRT_Nodule			Self;

	typedef itk::Image< TPixelType, TDimension >		ScalarImageType;	
	typedef typename ScalarImageType::Pointer			ScalarImageTypePointer;
	typedef itk::Image< int, TDimension>				ScalarLabelImageType;
	typedef typename ScalarLabelImageType::Pointer 		ScalarLabelImageTypePointer;
	typedef itk::Vector< TPixelType, TInputImageNum >   VectorType;
	typedef itk::Image< VectorType,TDimension >			VectorImageType;
	typedef typename VectorImageType::Pointer			VectorImageTypePointer;
	typedef itk::Vector< TPixelType, 2 >   				LabValVectorType;
	typedef itk::Image< LabValVectorType,TDimension >	LabValImageType;
	typedef typename LabValImageType::Pointer			LabValImageTypePointer;

	typedef itk::ImageRegionIterator< ScalarImageType >			ScalarIteratorType;
	typedef itk::NeighborhoodIterator< ScalarImageType >		ScalarNeighborhoodIteratorType;
	typedef itk::ImageRegionIterator< ScalarLabelImageType >	ScalarLabelImageIteratorType;
	typedef itk::ImageRegionIterator< VectorImageType >			VectorIteratorType;
	typedef itk::NeighborhoodIterator< VectorImageType >		VectorNeighborhoodIteratorType;
	typedef itk::ImageRegionIterator< LabValImageType >			LabValIteratorType;
	typedef itk::NeighborhoodIterator< LabValImageType >		LabValNeighborhoodIteratorType;
	typedef typename ScalarNeighborhoodIteratorType::RadiusType		ScalarNeighborhoodRadiusType;
	typedef typename VectorNeighborhoodIteratorType::RadiusType		VectorNeighborhoodRadiusType;
	typedef typename LabValNeighborhoodIteratorType::RadiusType		LabValNeighborhoodRadiusType;

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
	typedef itk::MedianImageFilter<ScalarImageType,ScalarImageType> ScalarMedianFilterType;
	typedef itk::ImageDuplicator< ScalarImageType > DuplicatorType;

	//typedef itk::LabelStatisticsImageFilter< ScalarImageType, ScalarLabelImageType > LabelStatisticsImageFilterType;
	typedef itk::MinimumMaximumImageCalculator< ScalarImageType > ScalarMinMaxImageCalculatorType;
	typedef itk::ThresholdImageFilter< ScalarImageType >  ScalarThresholdImageFilterType;
	typedef itk::ConnectedComponentImageFilter< ScalarLabelImageType, ScalarLabelImageType > ScalarConnectedComponentImageFilterType;
	typedef itk::BinaryBallStructuringElement< TPixelType, TDimension > StructureElementType;
	typedef itk::BinaryDilateImageFilter< ScalarImageType, ScalarImageType, StructureElementType> BinaryDilateImageFilterType;
	typedef itk::BinaryContourImageFilter< ScalarImageType, ScalarImageType > BinaryContourImageFilterType;

private:
	SBRT_Nodule(const Self &);
	void operator=(const Self &);

	VectorImageTypePointer m_InputImageVector;
	std::vector<ScalarImageTypePointer> m_MfImage;
	//ScalarImageTypePointer m_MaskImage;
	ScalarImageTypePointer m_GraMagnitude;
	ScalarImageTypePointer m_OutputImage;
	ScalarImageTypePointer m_SeedImage;
	ScalarLabelImageTypePointer m_labelingImage;
	ScalarLabelImageTypePointer m_segImage;
	ScalarImageTypePointer m_voxelID;
	VectorNeighborhoodRadiusType m_VectorRadius;
	ScalarNeighborhoodRadiusType m_ScalarRadius;
	VectorSizeType m_ImageSize;
	std::vector< float > m_ImageVoxelSize;
	vnl_matrix< float > m_ImageFea;

	std::vector<float> m_modailityWt;
	float m_disThreshold;	//! the distance between the boundary of candidate nodule region and the background seeds
	float m_sigma; 			//! 
	float m_suvRatio;		//! threshold ratio (in SUV) used to identify nodule seeds

    int m_lf_lab;           //! label value for lung field
	int m_smoothIterNum;	//! number of iterations for median filtering of the input images
	int m_smoothR;			//! radius of neighborhood for median filtering
	int m_minNumFgSeed;		//! minimum number of nodule seeds (in voxel)
	int m_validVxNum;		//! number of foreground voxels that are involved in the segmentation
	int m_nei_r;			//! radius of neighborhood for label information propagation
	float m_alpha;			//! ratio parameter in label information propagation, 0<m_alpha<1
	int m_iterNum;			//! number of iterations for label information propagation

	std::string m_outputBaseName;
	cbica::Logging *m_logger;
	bool m_logger_or_not;

public:
	SBRT_Nodule();
	~SBRT_Nodule();

    /**
    \brief Initialization function reads and pre-processes the input image, and initialize the output images.

    \param InputImageName contains the names of input images.
    \param maskName is the file name of lung filed/foreground mask image within which the nodule segemtnation will be performed.
    */
	int Initialization( std::vector<std::string> &InputImageName, std::string maskName );
	
	/**
    \brief SetParameters function sets main parameters for generating seeds, creating graph, and label information propagation.
    */
	void SetParameters( int lfLab, std::string obasename, float suvR, int mnfs, std::vector<float> &mWt, float dThres, float wSigma, int iter, int siter, int sr );

	void SetLogger(std::string logName);

	/**
    \brief GenerateSeeds function generates nodule/background seeds that are used for label information propagation.
    */
	void GenerateSeeds();

	/**
    \brief ReadSeedImage function loads nodule/background seeds from optional input image where available.
    */
    void ReadSeedImage( std::string seedName );

    /**
    \brief PerformSegmentation function performs label information propagation for nodule segmentation.
    */
	void PerformSegmentation();

	/**
    \brief WriteLabel function saves the nodule segmentation results.
    */
	void WriteLabel(int seedAv);
};

#include "SBRT_Nodule.hxx"
