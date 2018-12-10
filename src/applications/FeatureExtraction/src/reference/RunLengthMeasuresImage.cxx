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
#include "itkScalarImageToRunLengthFeaturesFilter.h"

#include <string>
#include <vector>

#define isnan(x) ((x) != (x))

//predefine functions

void usage( int argc, char *argv[] );
void parse(int argc, char *argv[],char* &fn_input,char* &fn_output, char* &fn_mask4Process,unsigned int &winSZ,
		   unsigned int &numberOfBins,unsigned int &slidingDist,std::vector<int> *offsetvectors,int &numVectors);

template<class TValue> TValue Convert( std::string optionString );
template<class TValue> std::vector<TValue> ConvertVector( std::string optionString );
template <unsigned int ImageDimension> int GenerateRunLengthMeasureImage( int argc, char *argv[] );

//maind code

int main( int argc, char *argv[] )
{
  	
	  usage(argc,argv);
	  switch( atoi( argv[1] ) )
	   {
	   case 2:
		 GenerateRunLengthMeasureImage<2>( argc, argv );
		 break;
	   case 3:
		 GenerateRunLengthMeasureImage<3>( argc, argv );
		 break;
	   default:
		  std::cerr << "Unsupported dimension" << std::endl;
		  exit( EXIT_FAILURE );
	   }
}

//const unsigned int ImageDimension = 2;
template <unsigned int ImageDimension>
int GenerateRunLengthMeasureImage( int argc, char *argv[] )
{

	//parse inputs
	char* inputImageFn=NULL; char* fn_output=NULL;
	char* fn_mask4Process=NULL;
	unsigned int winSZ; unsigned int numberOfBins;	unsigned int slidingDist;
	std::vector<int> offset_vectors[10]; int numVectors;

	parse(argc,argv,inputImageFn,fn_output,fn_mask4Process,winSZ,numberOfBins,slidingDist,offset_vectors,numVectors);
		  
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
	const typename ImageType::SpacingType& inputImageSpacing = inputImage->GetSpacing();

    //get mask4Process
	FILE *fp=NULL;
	if(fn_mask4Process!=NULL)
        fp = fopen(fn_mask4Process,"r");
  
    typename ReaderType::Pointer mask4ProcessReader = ReaderType::New();
	typename ImageType::Pointer mask4Process = ImageType::New();
    if(!fp){//if mask not exist, create a mask using all voxels
		mask4Process->SetOrigin(inputImage->GetOrigin());
		mask4Process->SetSpacing(inputImage->GetSpacing());
		mask4Process->SetRegions(inputImage->GetLargestPossibleRegion());
		mask4Process->Allocate();
		mask4Process->FillBuffer( itk::NumericTraits<PixelType>::One);
    }
    else{//if mask exist, read the mask
	    mask4ProcessReader->SetFileName( fn_mask4Process );
	    mask4ProcessReader->Update();
	    mask4Process = mask4ProcessReader->GetOutput();
	    fclose(fp);
    }

	typename itk::ImageRegionIterator<ImageType> ItMGlobal( mask4Process,mask4Process->GetLargestPossibleRegion() );

   //reset win size  
    unsigned int winHalfSZ=(unsigned int)floor(((double)winSZ)/2);
    winSZ=winHalfSZ*2+1;

   //reset number of bins
    long maxnum=(long)(winSZ*winSZ)/3;
    numberOfBins=(numberOfBins>maxnum)?maxnum:numberOfBins;


	using HistogramFrequencyContainerType = itk::Statistics::DenseFrequencyContainer2;
	using RunLengthFilterType = itk::Statistics::ScalarImageToRunLengthFeaturesFilter< ImageType, HistogramFrequencyContainerType >;
	using RunLengthMatrixGenerator = typename RunLengthFilterType::RunLengthMatrixFilterType;
	using RunLengthFeatures = typename RunLengthFilterType::RunLengthFeaturesFilterType;
	using RunLengthFeatureName = typename RunLengthFilterType::RunLengthFeaturesFilterType;
	using InternalRunLengthFeatureName = typename RunLengthFeatureName::RunLengthFeatureName;
    //offfsets
    typename RunLengthMatrixGenerator::Pointer generator = RunLengthMatrixGenerator::New();
	typename RunLengthFilterType::FeatureNameVector::Pointer requestedFeatures = RunLengthFilterType::FeatureNameVector::New();
	typename RunLengthMatrixGenerator::OffsetVectorPointer offsets = RunLengthMatrixGenerator::OffsetVector::New();
	typename RunLengthMatrixGenerator::OffsetType offset;


	std::vector<int> vector;
	for(int i=0; i<numVectors; i++)
	{
		vector=offset_vectors[i];
		if( vector.size() > 0 && vector.size() != ImageDimension )
		  {
		  std::cerr << "Error:  offset size does not equal image dimension." << std::endl;
		  return EXIT_FAILURE;
		  }
		else
		  for ( unsigned int d = 0; d < ImageDimension; d++ )
			offset[d] = vector[d];

		offsets->push_back( offset );
	}
    //generator->SetOffsets( offsets );

    //local mask window
	typename ImageType::Pointer maskRoiIm = ImageType::New();
	typename ImageType::PointType origin;
	origin.Fill(0.0);
	maskRoiIm->SetOrigin(origin);
	maskRoiIm->SetSpacing(inputImage->GetSpacing());

	typename ImageType::SizeType size;
	size.Fill(winSZ);
	typename ImageType::RegionType region;
	region.SetSize(size);
	maskRoiIm->SetRegions(region);
	maskRoiIm->Allocate();
	maskRoiIm->FillBuffer( itk::NumericTraits<PixelType>::One);

    typename itk::ImageRegionIterator<ImageType> ItM( maskRoiIm,maskRoiIm->GetLargestPossibleRegion() );

	//get size of output image
    unsigned int outputImageSize[ImageDimension];
	float outputImageSpacing[ImageDimension];
    for(unsigned int i=0; i<ImageDimension; i++){
	    outputImageSize[i]=1+floor(((double)inputImageSize[i]-winSZ+1-1)/slidingDist);
		outputImageSpacing[i]=inputImageSize[i]/outputImageSize[i]*inputImageSpacing[i];
	}

    long nPixValNum=1;
	  for(unsigned int i=0;i<ImageDimension;i++)
		  nPixValNum*=outputImageSize[i];

    ///prepare for computing features
    float *psre=new float[nPixValNum]; float *plre=new float[nPixValNum]; float *pgln=new float[nPixValNum]; float *prln=new float[nPixValNum];
    float *prp=new float[nPixValNum]; float *plgre=new float[nPixValNum]; float *phgre=new float[nPixValNum]; 
	float *psrlge=new float[nPixValNum]; float *psrhge=new float[nPixValNum]; float *plrlge=new float[nPixValNum]; float *plrhge=new float[nPixValNum];

    //compute features for each pixel
	typename ImageType::SizeType imsz; imsz.Fill(winSZ);
	long pos=-1;
	for (ItIGlobal.GoToBegin(),ItMGlobal.GoToBegin(); !ItIGlobal.IsAtEnd(); ++ItIGlobal, ++ItMGlobal)
	{
		typename ImageType::IndexType thisidx=ItIGlobal.GetIndex();
		
		//cancel operations for pixels not in mask, and that will not be slided
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
		//if in mask
		typename ImageType::IndexType centidx=thisidx;
		for(unsigned int i=0;i<ImageDimension;i++)
			centidx[i]=centidx[i]+winHalfSZ;
		if(mask4Process->GetPixel( centidx )<=0)	
			continue;

		//show progress
		long howmany=(long)(nPixValNum*0.05);
		if(pos%howmany==0){
			int leftPercent=(int)((1-(double)pos/(double)nPixValNum)*100);
			std::cout << leftPercent << "% left ... ..." << std::endl;
		}

		//get roi of input image
	    typedef itk::RegionOfInterestImageFilter<ImageType,ImageType> ImageROIType;
	    typename ImageROIType::Pointer imageROI=ImageROIType::New();
	    typename ImageType::IndexType imstart=thisidx;

	    typename ImageType::RegionType desiredImageRegion;
	    desiredImageRegion.SetSize(imsz);
	    desiredImageRegion.SetIndex(imstart);
	    imageROI->SetRegionOfInterest(desiredImageRegion);
	    imageROI->SetInput(inputImage);
	    imageROI->Update();
	    typename ImageType::Pointer inputImageRoiIm=imageROI->GetOutput();

	    typename itk::ImageRegionIterator<ImageType> ItI( inputImageRoiIm,
		  inputImageRoiIm->GetLargestPossibleRegion() );

	  
	    ////////////////////////////////////////////////////////////////////////////////
	    //compute features here Note: ItI and ItM proved above
	    ////////////////////////////////////////////////////////////////////////////////
	    PixelType maxValue = itk::NumericTraits<PixelType>::min();
	    PixelType minValue = itk::NumericTraits<PixelType>::max();

	    for ( ItM.GoToBegin(), ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItM, ++ItI )
		{
			if ( ItM.Get() == itk::NumericTraits<PixelType>::One )
			  {
				if ( ItI.Get() < minValue )
					{	
						minValue = ItI.Get();
					}
				else if ( ItI.Get() > maxValue )
					{
						maxValue = ItI.Get();
					}
			  }

			if( isnan( ItI.Get() ) || std::isinf( ItI.Get() ) )
			  {
				ItM.Set( itk::NumericTraits<PixelType>::Zero );
			  }
		}

		generator->SetOffsets(offsets);
		generator->SetInput(inputImageRoiIm);
		generator->SetMaskImage( maskRoiIm );
		generator->SetNumberOfBinsPerAxis( numberOfBins );
		generator->SetPixelValueMinMax( minValue, maxValue );
		generator->SetDistanceValueMinMax( 0, 1.4142*(winSZ-1));
		generator->SetInsidePixelValue( 1 );
		generator->Update();

		//typedef itk::Statistics::HistogramToTextureFeaturesFilter<typename CooccurrenceMatrixGeneratorType::HistogramType> CalculatorType;
		//typename CalculatorType::Pointer calculator = CalculatorType::New();

		typename RunLengthFeatures::Pointer runLengthMatrixCalculator = RunLengthFeatures::New();
    auto temp = generator->GetOutput();
		runLengthMatrixCalculator->SetInput(generator->GetOutput());
		runLengthMatrixCalculator->Update();
		
		typedef float RealType;
		RealType sre = runLengthMatrixCalculator->GetShortRunEmphasis();
		RealType lre = runLengthMatrixCalculator->GetLongRunEmphasis();
		RealType gln = runLengthMatrixCalculator->GetGreyLevelNonuniformity();
		RealType rln = runLengthMatrixCalculator->GetRunLengthNonuniformity();
		RealType rp = static_cast<RealType>(runLengthMatrixCalculator->GetTotalNumberOfRuns()) / static_cast<RealType>(winSZ*winSZ);
		RealType lgre = runLengthMatrixCalculator->GetLowGreyLevelRunEmphasis();
		RealType hgre = runLengthMatrixCalculator->GetHighGreyLevelRunEmphasis();
		RealType srlge = runLengthMatrixCalculator->GetShortRunLowGreyLevelEmphasis();
		RealType srhge = runLengthMatrixCalculator->GetShortRunHighGreyLevelEmphasis();
		RealType lrlge = runLengthMatrixCalculator->GetLongRunLowGreyLevelEmphasis();
		RealType lrhge = runLengthMatrixCalculator->GetLongRunHighGreyLevelEmphasis();

	    psre[pos]=sre; plre[pos]=lre; pgln[pos]=gln; prln[pos]=rln; 
	    prp[pos]=rp; plgre[pos]=lgre; phgre[pos]=hgre;
		psrlge[pos]=srlge; psrhge[pos]=srhge; plrlge[pos]=lrlge; plrhge[pos]=lrhge; 

	}
	
	//save results
	char fullname[400];

	strcpy(fullname,fn_output); strcat(fullname,"_ShortRunEmphasis.nii");
	if(!writeImage<float,ImageDimension,1,ImageType>(fullname,psre,outputImageSize,outputImageSpacing,dirCosine))
		delete[] psre;
	strcpy(fullname,fn_output); strcat(fullname,"_LongRunEmphasis.nii");
	if(!writeImage<float,ImageDimension,1,ImageType>(fullname,plre,outputImageSize,outputImageSpacing,dirCosine))
		delete[] plre;
    strcpy(fullname,fn_output); strcat(fullname,"_GreyLevelNonuniformity.nii");
	if(!writeImage<float,ImageDimension,1,ImageType>(fullname,pgln,outputImageSize,outputImageSpacing,dirCosine))
		delete[] pgln;
	strcpy(fullname,fn_output); strcat(fullname,"_RunLengthNonuniformity.nii");
	if(!writeImage<float,ImageDimension,1,ImageType>(fullname,prln,outputImageSize,outputImageSpacing,dirCosine))
		delete[] prln;
	strcpy(fullname,fn_output); strcat(fullname,"_RunPercentage.nii");
	if(!writeImage<float,ImageDimension,1,ImageType>(fullname,prp,outputImageSize,outputImageSpacing,dirCosine))
		delete[] prp;
	strcpy(fullname,fn_output); strcat(fullname,"_LowGreyLevelRunEmphasis.nii");
	if(!writeImage<float,ImageDimension,1,ImageType>(fullname,plgre,outputImageSize,outputImageSpacing,dirCosine))
		delete[] plgre;
	strcpy(fullname,fn_output); strcat(fullname,"_HighGreyLevelRunEmphasis.nii");
	if(!writeImage<float,ImageDimension,1,ImageType>(fullname,phgre,outputImageSize,outputImageSpacing,dirCosine))
		delete[] phgre;
 
  delete[] psre; delete[] plre; delete[] pgln; delete[] prln;
  delete[] prp; delete[] plgre; delete[] phgre;
  return 0;
}

template<class TValue>
TValue Convert( std::string optionString )
			{
			TValue value;
			std::istringstream iss( optionString );
			iss >> value;
			return value;
			}

template<class TValue>
std::vector<TValue> ConvertVector( std::string optionString )
			{
			std::vector<TValue> values;
			std::string::size_type crosspos = optionString.find( 'x', 0 );

			if ( crosspos == std::string::npos )
					{
					values.push_back( Convert<TValue>( optionString ) );
					}
			else
					{
					std::string element = optionString.substr( 0, crosspos ) ;
					TValue value;
					std::istringstream iss( element );
					iss >> value;
					values.push_back( value );
					while ( crosspos != std::string::npos )
							{
							std::string::size_type crossposfrom = crosspos;
							crosspos = optionString.find( 'x', crossposfrom + 1 );
							if ( crosspos == std::string::npos )
									{
									element = optionString.substr( crossposfrom + 1, optionString.length() );
									}
							else
									{
									element = optionString.substr( crossposfrom + 1, crosspos ) ;
									}
							std::istringstream iss( element );
							iss >> value;
							values.push_back( value );
							}
					}
			return values;
}
void usage( int argc, char *argv[] )
{
	if ( argc < 4 )
    {
		std::cerr << std::endl << "Infor: " <<std::endl << argv[0] << " Calculate feature image of the run length measurements. "<< std::endl << std::endl;
		std::cerr << "Usage: " <<std::endl << argv[0] << " imageDimension inputImageName outputImageRootName " << std::endl
		<< "with options: " << std::endl
		<< "-template templateImageName" << std::endl
		<< "-winsz<11> window size" << std::endl
		<< "-numbins<50> number of bins" << std::endl
		<< "-slidingdist<1> sliding distance" << std::endl
		<< "-offset<1x1> nOffsets offset1 ... offsetn" << std::endl << std::endl;

		std::cerr << "Example: " << std::endl << argv[0] << " 2 ./data/15582-CXJ.nii ./result/15582-CXJ -template ./data/15582-CXJ-template.nii -winsze 11 -numbins 50 -slidingdist 3 -offset 2 1x0 1x1"
		<< std::endl << std::endl;
                std::cerr << "Note: " << " Output folder must be writable!"
		<< std::endl;
		exit( 1 );
    }

}

void parse(int argc, char *argv[],char* &fn_input,char* &fn_output, char* &fn_mask4Process,unsigned int &winSZ,unsigned int &numberOfBins,unsigned int &slidingDist,std::vector<int> *offsetvectors,int &numVectors)
{
	//set default values
	unsigned int ImageDimension=atoi(argv[1]);
	fn_input=NULL;
	fn_output=NULL;
	fn_mask4Process=NULL;
	winSZ=63; numberOfBins=50; slidingDist=63;

	numVectors=1;
	for(int d=0; d<ImageDimension; d++)
		offsetvectors[0].push_back(1);

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
		case 'w':
			winSZ=atoi(&myargv[1][0]);
			--myargc; --myargc;
			++myargv; ++myargv;
			break;
		case 'n':
			numberOfBins=atoi(&myargv[1][0]);
			--myargc; --myargc;
			++myargv; ++myargv;
			break;
		case 's':
			slidingDist=atoi(&myargv[1][0]);
			--myargc; --myargc;
			++myargv; ++myargv;
			break;
		case 'o':
			numVectors=atoi(&myargv[1][0]);
			for(int i=0;i<numVectors;i++)
				offsetvectors[i] = ConvertVector<int>( std::string( &myargv[2+i][0] ) );
			--myargc; --myargc;
			++myargv; ++myargv;
			for(int i=0;i<numVectors;i++){
				--myargc; ++myargv;
			}
			break;
		}
	}
}

