/*=========================================================================

  Program:   ITK General
  Language:  C++
  Date:      $Date: 2010/04/04 11:23:11 $
  Version:   $Revision: 1.22 $

  Author: Yuanjie Zheng (zheng.vision@gmail.com)
  Institution: PICSL

  This software is distributed WITHOUT ANY WARRANTY; without even 
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
  PURPOSE.

=========================================================================*/
#ifndef __BIASCORRECTION_ImageIO_h_
#define __BIASCORRECTION_ImageIO_h_

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImportImageFilter.h"
#include "itkImage.h"
#include "itkVector.h"

//typedef itk::Image<PixelType, ImageDimension> ImageType;

template <class PixelValType,unsigned int ImageDim,unsigned int PixelValDim>
bool readImage(const char* fn, PixelValType* &pData, unsigned int size[]);

//currently data can only save to a nifity file
template <class PixelValType,unsigned int ImageDim,unsigned int PixelValDim>
bool writeImage(const char* fn, PixelValType* pData, unsigned int size[]);

template <class PixelValType,unsigned int ImageDim,unsigned int PixelValDim>
bool readImage(const char* fn, PixelValType* &pData, unsigned int size[])
{
	//arbitrary variables
	unsigned int i,j;
	//all tpyes needed
	typedef itk::Vector< PixelValType,PixelValDim> PixelType;
	typedef itk::Image< PixelType, ImageDim >    ImageType;
	typedef itk::ImageFileReader< ImageType >  ReaderType;

	//read image
	typename ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( fn);
	try
	{
		reader->Update();
	}
	catch( itk::ExceptionObject & err ) 
    { 
		std::cerr << std::endl << "Input image can not be found !! Below is the information caught by ITK:" << std::endl << std::endl; 
		std::cerr << "ExceptionObject caught !" << std::endl; 
		std::cerr << err << std::endl; 
		return false;
    }

	//get image size
	typename ImageType::SizeType imageSize = 
	reader->GetOutput()->GetLargestPossibleRegion().GetSize();
	for(i=0;i<ImageDim;i++)
		size[i]=imageSize[i];

	
	//allocate space for the pData
	long nTotalValNum=1;
	for(i=0;i<ImageDim;i++)
		nTotalValNum*=size[i];
	nTotalValNum*=PixelValDim;

	pData=new PixelValType[nTotalValNum];

	//assign data value from reader to pData
	PixelType pixelVal;
	itk::ImageRegionConstIterator <ImageType> It( reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );

	long indVal=0;
	for(It.GoToBegin();!It.IsAtEnd();++It)
	{
			//pixelVal=It.GetPixel();
		pixelVal=It.Get();
		for(j=0;j<PixelValDim;j++)
			pData[indVal++]=pixelVal[j];

	}

	return true;
}

template <class PixelValType,unsigned int ImageDim,unsigned int PixelValDim, typename TImage>
bool writeImage(const char* fn, PixelValType* pData, unsigned int size[], float spacing[], typename TImage::DirectionType directionMatrix )
{
	//arbitrary variables
	unsigned int i,j;

	typedef float PixelType;//Terrible bug of ITK! from my experiments, PixelType can only be defined to unsigned char for saving a jpeg and bmp image,but can only be double for saving a nifty image. 
	typedef itk::Image< PixelType, ImageDim >    ImageType;
	typedef itk::ImageFileWriter< ImageType >  WriterType;
	typedef itk::ImportImageFilter< PixelType, ImageDim >   ImportFilterType;

	//Import filter
	typename ImportFilterType::Pointer importFilter = ImportFilterType::New();

	////set import filter
	//set region
	typename ImportFilterType::SizeType  imSize;
	for(i=0;i<ImageDim;i++)
		imSize[i]=size[i];

	typename ImportFilterType::IndexType start;
	start.Fill( 0 );

	typename ImportFilterType::RegionType region;
	region.SetIndex( start );
	region.SetSize(  imSize  );
	
	importFilter->SetRegion( region );
	
	//set origin
	double origin[ ImageDim ];
	for(i=0;i<ImageDim;i++)
		origin[i]=0.0;

	importFilter->SetOrigin( origin );

	//set spacing
	typename ImportFilterType::SpacingType  imSpacing;
	for(i=0;i<ImageDim;i++)
		imSpacing[i]=spacing[i];

	importFilter->SetSpacing( imSpacing );
	
	//set direction
	importFilter->SetDirection( directionMatrix );

	//copy data
	long nNumPixel=1;
	for(i=0;i<ImageDim;i++)
		nNumPixel*=size[i];
	PixelType * localBuffer = new PixelType[ nNumPixel ];
	
	PixelType pixelVal;
	long indVal=0;
	PixelValType val;
	for(i=0;i<nNumPixel;i++)
	{
		val=0;
		for(j=0;j<PixelValDim;j++)
			val+=pData[indVal++];
		val/=PixelValDim;
		pixelVal=static_cast<PixelType>(val);
		localBuffer[i]=pixelVal;
	}

	const bool importImageFilterWillOwnTheBuffer = true;
	importFilter->SetImportPointer( localBuffer, nNumPixel, 
                                  importImageFilterWillOwnTheBuffer );

	typedef itk::ImageFileWriter< ImageType > WriterType;
	typename WriterType::Pointer writer = WriterType::New();

	writer->SetFileName( fn );

	// Software Guide : BeginCodeSnippet
	writer->SetInput(  importFilter->GetOutput()  );

	try
    {
		writer->Update();
    }
	catch( itk::ExceptionObject & exp ) 
    {
		std::cerr << std::endl << "The image can not be written !! Below is the information caught by ITK:" << std::endl << std::endl; 
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << exp << std::endl;
		return false;
    }
	return true;
}

#endif // __BIASCORRECTION_ImageIO_h_
