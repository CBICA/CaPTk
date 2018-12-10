/*=========================================================================
  Program:   ITK General
  Language:  C++
  Date:      $Date: 2011/04/25 19:55:26 $
  ITK Version:   InsightToolkit-3.20.0

  Author: Yuanjie Zheng (zheng.vision@gmail.com)
  Institution: PICSL @ UPenn

  This software is distributed WITHOUT ANY WARRANTY; without even 
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
  PURPOSE.
=========================================================================*/
#include <stdio.h>

#include "ImageIOFunctions.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkBoundingBox.h"
#include "itkCastImageFilter.h"
//#include "itkGreyLevelRunLengthMatrixTextureCoefficientsCalculator.h"
#include "itkLabelStatisticsImageFilter.h"
//#include "itkMaskedScalarImageToGreyLevelRunLengthMatrixGenerator.h"

#include <string>
#include <vector>

#define isnan(x) ((x) != (x))

//predefine functions

void usage( int argc, char *argv[] );
void parse(int argc, char *argv[],char* &fn_input,char* &fn_output,unsigned int &winSZ,unsigned int &slidingDist);

template <unsigned int ImageDimension> int DownsampleImage( int argc, char *argv[] );

//maind code

int main( int argc, char *argv[] )
{
  	
	  usage(argc,argv);
	  switch( atoi( argv[1] ) )
	   {
	   case 2:
		 DownsampleImage<2>( argc, argv );
		 break;
	   case 3:
		 DownsampleImage<3>( argc, argv );
		 break;
	   default:
		  std::cerr << "Unsupported dimension" << std::endl;
		  exit( EXIT_FAILURE );
	   }
}


template <unsigned int ImageDimension>
int DownsampleImage( int argc, char *argv[] )
{

	//parse inputs
	char* inputImageFn=NULL; char* fn_output=NULL;
	char* fn_mask4Process=NULL;
	unsigned int winSZ; unsigned int numberOfBins;	unsigned int slidingDist;
	std::vector<int> offset_vectors[10]; int numVectors;

	parse(argc,argv,inputImageFn,fn_output,winSZ,slidingDist);
		  
	//Image type defination
	typedef float PixelType;
	typedef itk::Image<PixelType, ImageDimension> ImageType;
	typedef itk::ImageFileReader<ImageType> ReaderType;


    //get input image
    typename ReaderType::Pointer imageReader = ReaderType::New();
    imageReader->SetFileName( inputImageFn);
    imageReader->Update();
    typename ImageType::Pointer inputImage = imageReader->GetOutput();
    typename ImageType::SizeType inputImageSize = inputImage->GetLargestPossibleRegion().GetSize();
	typename itk::ImageRegionIterator<ImageType> ItIGlobal( inputImage,inputImage->GetLargestPossibleRegion() );
  	const typename ImageType::DirectionType& dirCosine = inputImage->GetDirection();
	
    //reset win size  
    unsigned int winHalfSZ=(unsigned int)floor(((double)winSZ)/2);
    winSZ=winHalfSZ*2+1;

	//get size of output image
    unsigned int outputImageSize[3];
    for(unsigned int i=0; i<ImageDimension; i++)
	    outputImageSize[i]=1+floor(((double)inputImageSize[i]-winSZ+1-1)/slidingDist);

    long nPixValNum=1;
	  for(unsigned int i=0;i<ImageDimension;i++)
		  nPixValNum*=outputImageSize[i];

    ///prepare for computing features
    float *outImage=new float[nPixValNum];

    //assign value to each pixel in downsampled image
	typename ImageType::SizeType imsz; imsz.Fill(winSZ);
	long pos=-1;
	for (ItIGlobal.GoToBegin(); !ItIGlobal.IsAtEnd(); ++ItIGlobal)
	{
		typename ImageType::IndexType thisidx=ItIGlobal.GetIndex();
		
		//if sliding point
		bool ifSlidingPoint=true;
		for(unsigned int i=0;i<ImageDimension;i++)
			if(thisidx[i]%slidingDist!=0)
				ifSlidingPoint=false;
		if(!ifSlidingPoint)
			continue;
		//if in domain
		bool ifInDomain=true;
		for(unsigned int i=0;i<ImageDimension;i++)
			if(thisidx[i]>inputImageSize[i]-winSZ)
				ifInDomain=false;
		if(!ifInDomain)
			continue;

		pos++;

		//show progress
		long howmany=(long)(nPixValNum*0.05);
		if(pos%howmany==0){
			int leftPercent=(int)((1-(double)pos/(double)nPixValNum)*100);
			std::cout << leftPercent << "% left ... ..." << std::endl;
		}

		//assign value
	    outImage[pos]=ItIGlobal.Get(); 

	}
	
	//save downsampled image
	//if(!writeImage<float,ImageDimension,1,ImageType>(fn_output,outImage,outputImageSize,dirCosine))
	//	delete[] outImage;
 
    delete[] outImage;
    return 0;
}


void usage( int argc, char *argv[] )
{
	if ( argc < 4 )
    {
		std::cerr << std::endl << "Infor: " <<std::endl << argv[0] << ": Downsample image by subsampling pixels based on sliding windows. Michael's version"<< std::endl << std::endl;
		std::cerr << "Usage: " <<std::endl << argv[0] << " imageDimension inputImageName outputImageName " << std::endl
		<< "with options: " << std::endl
		<< "-winsz<11> window size" << std::endl
		<< "-slidingdist<1> sliding distance" << std::endl << std::endl;

		std::cerr << "Example: " << std::endl << argv[0] << " 2 ./data/15582-CXJ.nii ./data/15582-CXJ-downsampled.nii  -slidingdist 3 -winsize 13"
		<< std::endl;
		exit( 1 );
    }

}

void parse(int argc, char *argv[],char* &fn_input,char* &fn_output,unsigned int &winSZ, unsigned int &slidingDist)
{
	//set default values
	unsigned int ImageDimension=atoi(argv[1]);
	fn_input=NULL;
	fn_output=NULL;
	winSZ=11;
    slidingDist=1;

	//parse
	fn_input=argv[2];
	fn_output=argv[3];

	int myargc=argc-4;
	char** myargv=argv+4;
	while((myargc > 1) && (myargv[0][0]=='-')){
		switch (myargv[0][1])
		{
		case 'w':
			winSZ=atoi(&myargv[1][0]);
			--myargc; --myargc;
			++myargv; ++myargv;
			break;
		case 's':
			slidingDist=atoi(&myargv[1][0]);
			--myargc; --myargc;
			++myargv; ++myargv;
			break;
		}
	}
}

