/*=========================================================================
  Program:   ITK General
  Language:  C++
  Date:      $Date: 2011/04/25 19:55:26 $
  ITK Version:   InsightToolkit-3.20.0

  Author: Yuanjie Zheng (zheng.vision@gmail.com)
  Institution: PICSL @ UPenn
  Note: Compute Fractal Dimension (FD) using differential box counting algorithm

  This software is distributed WITHOUT ANY WARRANTY; without even 
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
  PURPOSE.
=========================================================================*/
#include <stdio.h>

//#include "itkScalarToFractalImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkCastImageFilter.h"
//#include "itkGreyLevelCooccurrenceMatrixTextureCoefficientsCalculator.h"
#include "itkLabelStatisticsImageFilter.h"
//#include "itkMaskedScalarImageToGreyLevelCooccurrenceMatrixGenerator.h"
#include "ImageIOFunctions.h"

#include <string>
#include <vector>

#define isnan(x) ((x) != (x))

//predefine functions

void usage( int argc, char *argv[] );
void parse(int argc, char *argv[],char* &fn_input,char* &fn_output, char* &fn_mask4Process,unsigned int &radius);

template <unsigned int ImageDimension> int GenerateFDMeasuresWithDBC( int argc, char *argv[] );

float dot(const float * pVal1, const float * pVal2,unsigned int num);
float sum(const float * pVal, unsigned int num);

//maind code

int main( int argc, char *argv[] )
{
  	
	  usage(argc,argv);
	  switch( atoi( argv[1] ) )
	   {
	   case 2:
		 GenerateFDMeasuresWithDBC<2>( argc, argv );
		 break;
	   case 3:
		 GenerateFDMeasuresWithDBC<3>( argc, argv );
		 break;
	   default:
		  std::cerr << "Unsupported dimension" << std::endl;
		  exit( EXIT_FAILURE );
	   }
}


//const unsigned int ImageDimension = 2;
template <unsigned int ImageDimension>
int GenerateFDMeasuresWithDBC( int argc, char *argv[] )
{
	//parse inputs
	char* inputImageFn=NULL; char* fn_output=NULL;
	char* fn_mask4Process=NULL;
	unsigned int radius;
	parse(argc,argv,inputImageFn,fn_output,fn_mask4Process,radius);

	//Image type definition
	typedef float PixelType;

	//get input image
	PixelType* pInputImageData=NULL; unsigned int inputImageSize[3];
	if(!readImage<PixelType,ImageDimension,1>(inputImageFn, pInputImageData, inputImageSize))
		return 1;

	//pixel number
	long numPixel=1;
	for(unsigned int i=0;i<ImageDimension;i++)
		numPixel*=inputImageSize[i];

	//get mask4Process
	PixelType* pMaskImageData=NULL; unsigned int maskImageSize[3];
	FILE *fp=NULL;
	if(fn_mask4Process!=NULL)
        fp = fopen(fn_mask4Process,"r");
	if(!fp){//if mask not exist, create a mask using all voxels
		std::cout << "Mask not provided. All pixels will be processed!" << std::endl;
		pMaskImageData=new PixelType[numPixel];
		for(long l=0; l<numPixel; l++)
			pMaskImageData[l]=1;}
	else
		readImage<PixelType,ImageDimension,1>(fn_mask4Process, pMaskImageData, maskImageSize);
	
	//output image
	float *pFD=new float[numPixel];

	//preparation
	float *pB=new float[radius];
	float *logDiameters=new float[radius];
	for(unsigned int i=0; i<radius; i++)
		logDiameters[i]=log(double(i)*2+1);
	float Nxx=dot(logDiameters,logDiameters,radius)-sum(logDiameters,radius)*sum(logDiameters,radius)/radius;

	//repeat on each pixel
	long pos;
	for(unsigned int i1=radius; i1<inputImageSize[1]-radius; i1++)//number of rows
	{
		//show progress
		long howmany=(long)(inputImageSize[1]*0.05);
		if(i1%howmany==0){
			int leftPercent=(int)((1-(double)i1/(double)inputImageSize[1])*100);
			std::cout << leftPercent << "% left ... ..." << std::endl;
		}
		for(unsigned int i0=radius; i0<inputImageSize[0]-radius; i0++)//number of columns
		{//loop on each pixel
			pos=i0+i1*inputImageSize[0];

			//if in mask
			if(pMaskImageData[pos]<=0)
				continue;

			int maxDiameter=radius*2+1;
			for(int rd=1; rd<=radius; rd++){
				int diameter=rd*2+1;
				//max and min values
				PixelType maxValue = itk::NumericTraits<PixelType>::min();
				PixelType minValue = itk::NumericTraits<PixelType>::max();
				
				for(int mov1=-rd; mov1<=rd; mov1++){
					for(int mov0=-rd; mov0<=rd; mov0++){
						long tempos=(i0+mov0)+(i1+mov1)*inputImageSize[0];
						if ( pInputImageData[tempos] < minValue )
							minValue = pInputImageData[tempos];
						else if ( pInputImageData[tempos] > maxValue )
							maxValue = pInputImageData[tempos];
					}
				}
				
				PixelType F=floor((maxValue-minValue)/diameter+1);
				pB[radius-rd]=log(F*((double)maxDiameter*(double)maxDiameter/((double)diameter*(double)diameter)));
				
			}//loop on different radius
			float Nxy=dot(logDiameters,pB,radius)-sum(logDiameters,radius)*sum(pB,radius)/radius;
			pFD[pos]=Nxy/Nxx;
		}//loop on each pixel
	}

	//writeImage<float,ImageDimension,1>(fn_output,pFD,inputImageSize);

	delete[] pInputImageData; delete[] pMaskImageData; delete[] pFD; delete[] pB; delete[] logDiameters;
	return 0;
}

void usage( int argc, char *argv[] )
{
	if ( argc < 4 )
    {
		std::cerr << std::endl << "Infor: " <<std::endl << argv[0] << " Calculate image of fractal dimension with differential box counting (DBC) algorithm. "<< std::endl << std::endl;
		std::cerr << "Usage: " <<std::endl << argv[0] << " imageDimension inputImageName outputImageName " << std::endl
		<< "with options: " << std::endl
		<< "-radias<7>: radius" << std::endl
		<< "-template templateImageName" << std::endl
		<< std::endl;

		std::cerr << "Example: " << std::endl << argv[0] << " 2 ./data/15582-CXJ.nii ./result/15582-CXJ-fractal.nii -template ./data/15582-CXJ-template.nii -radius 8"
		<< std::endl << std::endl;
        std::cerr << "Note: " << " Output folder must be writable!"
		<< std::endl;
		exit( 1 );
    }

}

void parse(int argc, char *argv[],char* &fn_input,char* &fn_output, char* &fn_mask4Process,unsigned int &radius)
{
	//set default values
	unsigned int ImageDimension=atoi(argv[1]);
	fn_input=NULL;
	fn_output=NULL;
	fn_mask4Process=NULL;
	radius=7;


	//parse
	fn_input=argv[2];
	fn_output=argv[3];

	int myargc=argc-4;
	char** myargv=argv+4;
	while((myargc > 1) && (myargv[0][0]=='-')){
		switch (myargv[0][1])
		{
		case 't':
			fn_mask4Process=&myargv[1][0];
			--myargc; --myargc;
			++myargv; ++myargv;
			break;
		case 'r':
			radius=atoi(&myargv[1][0]);
			--myargc; --myargc;
			++myargv; ++myargv;
			break;
		}
	}
}

float dot(const float * pVal1, const float * pVal2,unsigned int num)
{
	float val=0;
	for(unsigned int i=0;i<num;i++)
		val+=pVal1[i]*pVal2[i];
	return val;
}
float sum(const float * pVal, unsigned int num)
{
	float val=0;
	for(unsigned int i=0;i<num;i++)
		val+=pVal[i];
	return val;
}