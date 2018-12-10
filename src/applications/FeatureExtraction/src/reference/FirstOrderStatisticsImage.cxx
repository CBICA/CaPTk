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

#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "ImageIOFunctions.h"


#define isnan(x) ((x) != (x))

using namespace std;

//predefine functions

void usage( int argc, char *argv[] );
void parse(int argc, char *argv[],char* &fn_input,char* &fn_output, char* &fn_mask4Process,unsigned int &winSZ,unsigned int &numberOfBins,unsigned int &slidingDist);


template <unsigned int ImageDimension>
int CalculateFirstOrderStatistics( int argc, char *argv[] );

//maind code

int main( int argc, char *argv[] )
{
  	
	  usage(argc,argv);
	  switch( atoi( argv[1] ) )
	   {
	   case 2:
		 CalculateFirstOrderStatistics<2>( argc, argv );
		 break;
	   case 3:
		 CalculateFirstOrderStatistics<3>( argc, argv );
		 break;
	   default:
		  std::cerr << "Unsupported dimension" << std::endl;
		  exit( EXIT_FAILURE );
	   }
}


template <unsigned int ImageDimension>
int CalculateFirstOrderStatistics( int argc, char *argv[] )
{

	//parse inputs
	char* inputImageFn=NULL; char* fn_output=NULL;
	char* fn_mask4Process=NULL;
	unsigned int winSZ; unsigned int numberOfBins;	unsigned int slidingDist;

	parse(argc,argv,inputImageFn,fn_output,fn_mask4Process,winSZ,numberOfBins,slidingDist);
		  
	//Image type defination
	typedef float PixelType;
	typedef itk::Image<PixelType, ImageDimension> ImageType;
	typedef itk::ImageFileReader<ImageType> ReaderType;

	//input image
	typename ReaderType::Pointer imageReader = ReaderType::New();
	imageReader->SetFileName( inputImageFn);
	imageReader->Update();
	typename ImageType::Pointer inputImage = imageReader->GetOutput();
	typename ImageType::SizeType inputImageSize = inputImage->GetLargestPossibleRegion().GetSize();
	typename itk::ImageRegionIterator<ImageType> ItIGlobal( inputImage,inputImage->GetLargestPossibleRegion() );
  	// Michael's modification
	const typename ImageType::DirectionType& dirCosine = inputImage->GetDirection();
	const typename ImageType::SpacingType& inputImageSpacing = inputImage->GetSpacing();

	//mask4Process
	FILE *fp=NULL;
	if(fn_mask4Process!=NULL)
        fp = fopen(fn_mask4Process,"r");

  
    typename ReaderType::Pointer mask4ProcessReader = ReaderType::New();
	typename ImageType::Pointer mask4Process = ImageType::New();
    if(!fp){//if mask not exist, create a mask using all voxels
		mask4Process->SetOrigin(inputImage->GetOrigin());
		mask4Process->SetSpacing(inputImage->GetSpacing());
		mask4Process->SetRegions(inputImage->GetLargestPossibleRegion());
		// mask4Process->SetDirection(inputImage->GetDirection());
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

   //win size  
  unsigned int winHalfSZ=(unsigned int)floor(((double)winSZ)/2);
  winSZ=winHalfSZ*2+1;

   //number of bins
  long maxnum=(long)(winSZ*winSZ)/3;
  numberOfBins=(numberOfBins>maxnum)?maxnum:numberOfBins;


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

  
  float *pmean=new float[nPixValNum];
  float *psigma=new float[nPixValNum];
  float *psum=new float[nPixValNum];
  float *pskewness=new float[nPixValNum];
  float *pkurtosis=new float[nPixValNum];
  float *pentropy=new float[nPixValNum];
  float *p5th=new float[nPixValNum];
  float *p95th=new float[nPixValNum];
  float *p5thmean=new float[nPixValNum];
  float *p95thmean=new float[nPixValNum];
  float *pmin=new float[nPixValNum];
  float *pmax=new float[nPixValNum];

   

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
		//cout << pos << endl;
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

		//local mask window
		typedef int MaskPixelType;
		typedef itk::Image<MaskPixelType, ImageDimension> MaskImageType;
		typename MaskImageType::Pointer maskRoiIm = MaskImageType::New();
		typename MaskImageType::PointType origin;
		origin.Fill(0.0);
		maskRoiIm->SetOrigin(inputImageRoiIm->GetOrigin());
		maskRoiIm->SetSpacing(inputImageRoiIm->GetSpacing());
		//maskRoiIm->SetRegions(inputImageRoiIm->GetLargestPossibleRegion());
		maskRoiIm->SetDirection(inputImageRoiIm->GetDirection());
		typename MaskImageType::SizeType size;
		size.Fill(winSZ);
		typename MaskImageType::RegionType region;
		region.SetSize(size);
		maskRoiIm->SetRegions(region);
		maskRoiIm->Allocate();
		maskRoiIm->FillBuffer(itk::NumericTraits<PixelType>::One);
		typename itk::ImageRegionIterator<MaskImageType> ItM(maskRoiIm, maskRoiIm->GetLargestPossibleRegion());
	    
	    ////////////////////////////////////////////////////////////////////////////////
	    //compute features here Note: ItI and ItM proved above
	    ////////////////////////////////////////////////////////////////////////////////
		 // maskRoiIm->FillBuffer( itk::NumericTraits<PixelType>::One);
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

				if( isnan( ItI.Get() ) || isinf( ItI.Get() ) )
				  {
					ItM.Set( itk::NumericTraits<PixelType>::Zero );
				  }
			}



		  PixelType mean; PixelType sigma; PixelType sum; PixelType variance; PixelType skewness; PixelType kurtosis;PixelType entropy = 0;

		  typedef itk::LabelStatisticsImageFilter<ImageType, MaskImageType> HistogramGeneratorType;
		  typename HistogramGeneratorType::Pointer stats = HistogramGeneratorType::New();
		  stats->SetInput( inputImageRoiIm );
		  stats->SetLabelInput( maskRoiIm );
		  stats->UseHistogramsOn();
		  stats->SetHistogramParameters( numberOfBins, minValue, maxValue );
		  stats->Update();

		  mean = stats->GetMean( 1 );
		  sum = stats->GetSum( 1 );
		  sigma = stats->GetSigma( 1 );
		  variance = sigma * sigma;

		  kurtosis = 0.0;
		  skewness = 0.0;

		  PixelType N = 0.0;
		  for ( ItI.GoToBegin(), ItM.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItM )
			{
			if ( ItM.Get() == itk::NumericTraits<PixelType>::One )
			  {
			  PixelType value = ItI.Get();

			  PixelType diff = value - mean;
			  skewness += ( diff * diff * diff );
			  kurtosis += ( diff * diff * diff * diff );

			  N += 1.0;
			  }
			}
		  skewness /= ( ( N - 1 ) * variance * sigma );
		  kurtosis /= ( ( N - 1 ) * variance * variance );

		double fifthPercentileValue = 0.0;
		double ninetyFifthPercentileValue = 0.0;

		double fifthPercentileMean = 0.0;
		double fifthN = 0.0;
		double ninetyFifthPercentileMean = 0.0;
		double ninetyFifthN = 0.0;

				typedef typename HistogramGeneratorType::HistogramType  HistogramType;
				const HistogramType *histogram = stats->GetHistogram( 1 );

				if( !histogram )
						{
						std::cerr << "ERROR:  No histogram created." << std::endl;
                        continue;
						}

				fifthPercentileValue = histogram->Quantile( 0, 0.05 );
				ninetyFifthPercentileValue = histogram->Quantile( 0, 0.95 );

				entropy = 0.0;
				for( unsigned int i = 0; i < histogram->Size(); i++ )
						{
						PixelType p = static_cast<PixelType>( histogram->GetFrequency( i, 0 )  )
								/ static_cast<PixelType>( histogram->GetTotalFrequency() );
						if ( p > 0 )
								{
								entropy += ( -p * vcl_log( p ) / vcl_log( 2.0 ) );
								}
						}

				for ( ItM.GoToBegin(), ItI.GoToBegin(); !ItM.IsAtEnd(); ++ItM, ++ItI )
						{
						PixelType value = ItI.Get();
						if ( value <= fifthPercentileValue )
								{
								fifthPercentileMean += value;
								fifthN++;
								}
						else if ( value >= ninetyFifthPercentileValue )
								{
								ninetyFifthPercentileMean += value;
								ninetyFifthN++;
								}
						}

				fifthPercentileMean /= fifthN;
				ninetyFifthPercentileMean /= ninetyFifthN;

			//assign each feature value
				//in total 12 features
		  pmean[pos]=mean; psigma[pos]=sigma; psum[pos]=sum; pskewness[pos]=skewness; pkurtosis[pos]=kurtosis; pentropy[pos]=entropy;
		  p5th[pos]=fifthPercentileValue; p95th[pos]=ninetyFifthPercentileValue; p5thmean[pos]=fifthPercentileMean; 
		  p95thmean[pos]=ninetyFifthPercentileMean; pmin[pos]=minValue; pmax[pos]=maxValue;
	}

  //save results
	char fullname[400];

	strcpy(fullname,fn_output); strcat(fullname,"_mean.nii");
	if(!writeImage<float,ImageDimension,1,ImageType>(fullname,pmean,outputImageSize,outputImageSpacing,dirCosine))
		delete[] pmean;

	strcpy(fullname,fn_output); strcat(fullname,"_sigma.nii");
	if(!writeImage<float,ImageDimension,1,ImageType>(fullname,psigma,outputImageSize,outputImageSpacing,dirCosine))
		delete[] psigma;

	strcpy(fullname,fn_output); strcat(fullname,"_sum.nii");
	if(!writeImage<float,ImageDimension,1,ImageType>(fullname,psum,outputImageSize,outputImageSpacing,dirCosine))
		delete[] psum;

	strcpy(fullname,fn_output); strcat(fullname,"_skewness.nii");
	if(!writeImage<float,ImageDimension,1,ImageType>(fullname,pskewness,outputImageSize,outputImageSpacing,dirCosine))
		delete[] pskewness;

	strcpy(fullname,fn_output); strcat(fullname,"_kurtosis.nii");
	if(!writeImage<float,ImageDimension,1,ImageType>(fullname,pkurtosis,outputImageSize,outputImageSpacing,dirCosine))
		delete[] pkurtosis;

	strcpy(fullname,fn_output); strcat(fullname,"_entropy.nii");
	if(!writeImage<float,ImageDimension,1,ImageType>(fullname,pentropy,outputImageSize,outputImageSpacing,dirCosine))
		delete[] pentropy;

	strcpy(fullname,fn_output); strcat(fullname,"_5th.nii");
	if(!writeImage<float,ImageDimension,1,ImageType>(fullname,p5th,outputImageSize,outputImageSpacing,dirCosine))
		delete[] p5th;

	strcpy(fullname,fn_output); strcat(fullname,"_95th.nii");
	if(!writeImage<float,ImageDimension,1,ImageType>(fullname,p95th,outputImageSize,outputImageSpacing,dirCosine))
		delete[] p95th;

	strcpy(fullname,fn_output); strcat(fullname,"_5thmean.nii");
	if(!writeImage<float,ImageDimension,1,ImageType>(fullname,p5thmean,outputImageSize,outputImageSpacing,dirCosine))
		delete[] p5thmean;

	strcpy(fullname,fn_output); strcat(fullname,"_95thmean.nii");
	if(!writeImage<float,ImageDimension,1,ImageType>(fullname,p95thmean,outputImageSize,outputImageSpacing,dirCosine))
		delete[] p95thmean;

	strcpy(fullname,fn_output); strcat(fullname,"_min.nii");
	if(!writeImage<float,ImageDimension,1,ImageType>(fullname,pmin,outputImageSize,outputImageSpacing,dirCosine))
		delete[] pmin;

	strcpy(fullname,fn_output); strcat(fullname,"_max.nii");
	if(!writeImage<float,ImageDimension,1,ImageType>(fullname,pmax,outputImageSize,outputImageSpacing,dirCosine))
		delete[] pmax;

  delete[] pmean; delete[] psigma; delete[] psum; delete[] pskewness; delete[] pkurtosis; delete[] pentropy;
  delete[] p5th; delete[] p95th; delete[] p5thmean; delete[] p95thmean; delete[] pmin; delete[] pmax;
  //system("pause");
  return 0;
}

void usage( int argc, char *argv[] )
{
	if ( argc < 4 )
    {
		std::cerr << std::endl << "Infor: " <<std::endl << argv[0] << " Calculate feature image for the 1st-order statistic features. "<< std::endl << std::endl;
		std::cerr << "Usage: " <<std::endl << argv[0] << " imageDimension inputImageName outputImageRootName " << std::endl
		<< "with options: " << std::endl
		<< "-template templateImageName" << std::endl
		<< "-winsz<11> window size" << std::endl
		<< "-numbins<50> number of bins" << std::endl;

		std::cerr << "Example: " << std::endl << argv[0] << " 2 ./data/15582-CXJ.nii ./result/15582-CXJ -template template.nii -winsze 11 -numbins 50"
		<< std::endl << std::endl;
                std::cerr << "Note: " << " Output folder must be writable!"
		<< std::endl;
		exit( 1 );
    }

}

void parse(int argc, char *argv[],char* &fn_input,char* &fn_output, char* &fn_mask4Process,unsigned int &winSZ,unsigned int &numberOfBins,unsigned int &slidingDist)
{
	//set default values
	unsigned int ImageDimension=atoi(argv[1]);
	fn_input=NULL;
	fn_output=NULL;
	fn_mask4Process=NULL;
	winSZ=63; numberOfBins=50; slidingDist=63;

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
			break;
		case 'w':
			winSZ=atoi(&myargv[1][0]);
			break;
		case 's':
			slidingDist=atoi(&myargv[1][0]);
			break;
		case 'n':
			numberOfBins=atoi(&myargv[1][0]);
			break;
		}
		--myargc; --myargc;
		++myargv; ++myargv;
	}
}

