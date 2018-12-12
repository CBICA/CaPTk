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
#include "iostream"
#include "vector"
#include "string.h"
#include "math.h"
#include "itkImage.h"
#include "itkConvolutionImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "ImageIOFunctions.h"
#include "itkAddImageFilter.h"
#include <itkPowImageFilter.h>
#include "itkWrapPadImageFilter.h"
#include "itkForwardFFTImageFilter.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkFFTShiftImageFilter.h"
#include <stdlib.h>
#define PI 3.14159265358979323846264


#define isnan(x) ((x) != (x))

using namespace std;

//predefine functions

void usage(int argc, char *argv[]);
void parse(int argc, char *argv[], char* &fn_input, char* &fn_output, char* &fn_mask4Process, unsigned int &winSZ, unsigned int &slidingDist, double &R, double &fmax, double &gamma_gabor, int &level, int &orientation);


template <unsigned int ImageDimension>
int CalculateGaborWavelet(int argc, char *argv[]);
vector<vector <double>> gabor_wavelet(int argc, char* argv[], double R0, double fmax0, double gamma_gabor0, int level0, int orientation0);

//maind code

int main(int argc, char *argv[])
{

	usage(argc, argv);
	switch (atoi(argv[1]))
	{
	case 2:
		CalculateGaborWavelet<2>(argc, argv);
		break;
	case 3:
		CalculateGaborWavelet<3>(argc, argv);
		break;
	default:
		std::cerr << "Unsupported dimension" << std::endl;
		exit(EXIT_FAILURE);
	}
		return 0;
}

template <unsigned int ImageDimension>
int CalculateGaborWavelet(int argc, char *argv[])
{

	//parse inputs
	char* inputImageFn = NULL; char* fn_output = NULL;
	char* fn_mask4Process = NULL;
	unsigned int winSZ;	unsigned int slidingDist;
	double R, fmax, gamma_gabor;
	int level, orientation;

	parse(argc, argv, inputImageFn, fn_output, fn_mask4Process, winSZ, slidingDist, R, fmax, gamma_gabor, level, orientation);

	typedef double PixelType;
	typedef itk::Image<PixelType, ImageDimension> ImageType;
	typedef itk::ImageFileReader<ImageType> ReaderType;

	//input image
	typename ReaderType::Pointer imageReader = ReaderType::New();
	imageReader->SetFileName(inputImageFn);
	imageReader->Update();
	typename ImageType::Pointer inputImage = imageReader->GetOutput();
	typename ImageType::SizeType inputImageSize = inputImage->GetLargestPossibleRegion().GetSize();
	typename itk::ImageRegionIterator<ImageType> ItIGlobal(inputImage, inputImage->GetLargestPossibleRegion());


	// Michael's modification
	const typename ImageType::DirectionType& dirCosine = inputImage->GetDirection();
	const typename ImageType::SpacingType& inputImageSpacing = inputImage->GetSpacing();

	//mask4Process
	FILE *fp = NULL;
	if (fn_mask4Process != NULL)
		fp = fopen(fn_mask4Process, "r");


	typename ReaderType::Pointer mask4ProcessReader = ReaderType::New();
	typename ImageType::Pointer mask4Process = ImageType::New();
	if (!fp){//if mask not exist, create a mask using all voxels
		mask4Process->SetOrigin(inputImage->GetOrigin());
		mask4Process->SetSpacing(inputImage->GetSpacing());
		mask4Process->SetRegions(inputImage->GetLargestPossibleRegion());
		// mask4Process->SetDirection(inputImage->GetDirection());
		mask4Process->Allocate();
		mask4Process->FillBuffer(itk::NumericTraits<PixelType>::One);
	}
	else{//if mask exist, read the mask
		mask4ProcessReader->SetFileName(fn_mask4Process);
		mask4ProcessReader->Update();
		mask4Process = mask4ProcessReader->GetOutput();
		fclose(fp);
	}
	typename itk::ImageRegionIterator<ImageType> ItMGlobal(mask4Process, mask4Process->GetLargestPossibleRegion());

	//win size  
	unsigned int winHalfSZ = (unsigned int)floor(((double)winSZ) / 2);
	winSZ = winHalfSZ * 2 + 1;

	//get size of output image
	unsigned int outputImageSize[ImageDimension];
	float outputImageSpacing[ImageDimension];
	for (unsigned int i = 0; i<ImageDimension; i++){
		outputImageSize[i] = 1 + floor(((double)inputImageSize[i] - winSZ + 1 - 1) / slidingDist);
		outputImageSpacing[i] = inputImageSize[i] / outputImageSize[i] * inputImageSpacing[i];
	}
	long nPixValNum = 1;
	for (unsigned int i = 0; i<ImageDimension; i++)
		nPixValNum *= outputImageSize[i];

	//local mask window
	typedef int MaskPixelType;
	typedef itk::Image<MaskPixelType, ImageDimension> MaskImageType;
	typename MaskImageType::Pointer maskRoiIm = MaskImageType::New();
	typename MaskImageType::PointType origin;
	origin.Fill(0.0);
	maskRoiIm->SetOrigin(origin);
	maskRoiIm->SetSpacing(inputImage->GetSpacing());
	typename MaskImageType::SizeType size;
	size.Fill(winSZ);
	typename MaskImageType::RegionType region;
	region.SetSize(size);
	maskRoiIm->SetRegions(region);
	maskRoiIm->Allocate();
	maskRoiIm->FillBuffer(itk::NumericTraits<PixelType>::One);
	typename itk::ImageRegionIterator<MaskImageType> ItM(maskRoiIm, maskRoiIm->GetLargestPossibleRegion());


	//compute features for each pixel
	typename ImageType::SizeType imsz; imsz.Fill(winSZ);
	long pos = -1;
	int count_lattice = 0;
	vector<vector<double>> gabor_features;
	vector<double> lattice_feature;
	for (ItIGlobal.GoToBegin(), ItMGlobal.GoToBegin(); !ItIGlobal.IsAtEnd(); ++ItIGlobal, ++ItMGlobal)
	{
		typename ImageType::IndexType thisidx = ItIGlobal.GetIndex();

		//if sliding point
		bool ifSlidingPoint = true;
		for (unsigned int i = 0; i<ImageDimension; i++)
			if (thisidx[i] % slidingDist != 0)
				ifSlidingPoint = false;
		if (!ifSlidingPoint)
			continue;
		//if in domain
		bool ifInDomain = true;
		for (unsigned int i = 0; i<ImageDimension; i++)
			if (thisidx[i]>inputImageSize[i] - winSZ)
				ifInDomain = false;
		if (!ifInDomain)
			continue;

		pos++;

		//if in mask
		typename ImageType::IndexType centidx = thisidx;
		for (unsigned int i = 0; i<ImageDimension; i++)
			centidx[i] = centidx[i] + winHalfSZ;
		if (mask4Process->GetPixel(centidx) <= 0)
			continue;

		count_lattice++;

		//show progress
		long howmany = (long)(nPixValNum*0.05);
		if (pos%howmany == 0)
		{
			int leftPercent = (int)((1 - (double)pos / (double)nPixValNum) * 100);
			std::cout << leftPercent << "% left ... ..." << std::endl;
		}

		//get roi of input image
		typedef itk::RegionOfInterestImageFilter<ImageType, ImageType> ImageROIType;
		typename ImageROIType::Pointer imageROI = ImageROIType::New();
		typename ImageType::IndexType imstart = thisidx;

		typename ImageType::RegionType desiredImageRegion;
		desiredImageRegion.SetSize(imsz);
		desiredImageRegion.SetIndex(imstart);
		imageROI->SetRegionOfInterest(desiredImageRegion);
		imageROI->SetInput(inputImage);
		imageROI->Update();
		typename ImageType::Pointer inputImageRoiIm = imageROI->GetOutput();

		typename itk::ImageRegionIteratorWithIndex<ImageType> ItI(inputImageRoiIm,
			inputImageRoiIm->GetLargestPossibleRegion());

		//Yifan's modification


		////////////////////////////////////////////////////////////////////////////////
		//compute features here Note: ItI and ItM proved above
		////////////////////////////////////////////////////////////////////////////////		
		lattice_feature.clear();
		//default u and v are from "https://www.mathworks.com/matlabcentral/fileexchange/44630-gabor-feature-extraction"
		for (int u = 1; u < level+1; u++)
		{
			for (int v = 1; v < orientation+1; v++)
			{
				//create a gabor filter with different parameters
				vector<vector <double>> result = gabor_wavelet(argc, argv, R,fmax, gamma_gabor, u, v);
				typename ImageType::RegionType region = inputImageRoiIm->GetLargestPossibleRegion();
				typename ImageType::SizeType size = region.GetSize();
				//result image
				typename ImageType::Pointer result_matrix = ImageType::New();
				typename ImageType::RegionType region_new;
				typename ImageType::IndexType start;
				start.Fill(0);
				//try to ignore the boundary effect, we can do the convolution without padding the image
				typename ImageType::SizeType size_new_matrix;
				size_new_matrix[0] = size[0] - R + 1;
				size_new_matrix[1] = size[1] - R + 1;
				region_new.SetSize(size_new_matrix);
				result_matrix->SetRegions(region_new);
				result_matrix->Allocate();
				result_matrix->FillBuffer(itk::NumericTraits< double >::Zero);

				vector<vector <double>> image_matrix;
				vector<double> temp2;
				for (int j = 0; j < size[1]; j++)
				{
					temp2.clear();
					for (int i = 0; i < size[0]; i++)
					{
						typename ImageType::IndexType index;
						index[0] = i;
						index[1] = j;
						temp2.push_back(inputImageRoiIm->GetPixel(index));
					}
					image_matrix.push_back(temp2);
				}

				typename ImageType::IndexType pixelIndex;

				for (int i = 0; i < size_new_matrix[1]; i++)
				{
					for (int j = 0; j < size_new_matrix[0]; j++)
					{
						double sum = 0;
						for (int k = 0; k < R; k++)
						{
							for (int l = 0; l < R; l++)
							{
								sum = sum + result[k][l] * image_matrix[i + k][j + l];
							}
						}
						pixelIndex[0] = j;
						pixelIndex[1] = i;
						result_matrix->SetPixel(pixelIndex, sum);
					}
				}
				
				typedef itk::StatisticsImageFilter<ImageType> StatisticsImageFilterType;
				typename StatisticsImageFilterType::Pointer statisticsImageFilter
					= StatisticsImageFilterType::New();
				statisticsImageFilter->SetInput(result_matrix);
				statisticsImageFilter->Update();
				lattice_feature.push_back(statisticsImageFilter->GetMean());
			}
		}

		gabor_features.push_back(lattice_feature);
	}

	
	for (int i = 0; i < gabor_features[0].size(); i++)
	{
		char fullname[300];
		//char *buffer1, buffer2[2];
		int integer_part = i / orientation + 1;
		int residual = i % int(orientation)+1;
		//string path = "C:/folder/branches/Libracxx/data/New_folder_power/";
		strcpy(fullname, fn_output);
		strcat(fullname, "_scale_");
		string s1 = std::to_string(integer_part);
		char const *pchar1 = s1.c_str();
		strcat(fullname, pchar1);
		strcat(fullname, "_orientation_");
		string s2 = std::to_string(residual);
		char const *pchar2 = s2.c_str();
		strcat(fullname, pchar2);
		strcat(fullname, ".nii");
		//float *temp = new float[nPixValNum];

		typename ImageType::Pointer result_image = ImageType::New();
		typename ImageType::RegionType result_region;
		typename ImageType::SizeType result_size;
		typename ImageType::IndexType result_index;
		result_size[0] = floor((inputImageSize[0] - winSZ) / slidingDist) + 1;
		result_size[1] = floor((inputImageSize[1] - winSZ) / slidingDist) + 1;
		result_region.SetSize(result_size);
		result_index.Fill(0);
		result_region.SetIndex(result_index);
		result_image->SetRegions(result_region);
    result_image->SetSpacing(outputImageSpacing);
    result_image->SetDirection(inputImage->GetDirection());
		result_image->Allocate();
		result_image->FillBuffer(itk::NumericTraits< float >::Zero);

		for (int k = 0; k < result_size[0]; k++)
		{
			for (int l = 0; l < result_size[1]; l++)
			{
				typename ImageType::IndexType temp_index;
				temp_index[0] = k;
				temp_index[1] = l;
				result_image->SetPixel(temp_index, gabor_features[l*result_size[0] + k][i]);
			}
		}
		
		typedef itk::ImageFileWriter< ImageType  > WriterType;
		typename WriterType::Pointer writer = WriterType::New();
		writer->SetFileName(fullname);
		writer->SetInput(result_image);
	
		try
		{
			writer->Update();
		}
		catch (itk::ExceptionObject & error)
		{
			std::cerr << "Error: " << error << std::endl;
			return EXIT_FAILURE;
		}

	}
	
	return 0;
}

void usage(int argc, char *argv[])
{
	if (argc < 3)
	{
		std::cerr << std::endl << "Infor: " << std::endl << argv[0] << " Calculate feature image for the 1st-order statistic features. " << std::endl << std::endl;
		std::cerr << "Usage: " << std::endl << argv[0] << " imageDimension inputImageName outputImageRootName " << std::endl
			<< "with options: " << std::endl
			<< "-template templateImageName" << std::endl
			<< "-winsz<11> window size" << std::endl
			<< "-numbins<50> number of bins" << std::endl;

		std::cerr << "Example: " << std::endl << argv[0] << " 2 ./data/15582-CXJ.nii ./result/15582-CXJ -template template.nii -winsze 11 -numbins 50"
			<< std::endl << std::endl;
		std::cerr << "Note: " << " Output folder must be writable!"
			<< std::endl;
		exit(1);
	}

}

void parse(int argc, char *argv[], char* &fn_input, char* &fn_output, char* &fn_mask4Process, unsigned int &winSZ, unsigned int &slidingDist, double &R, double &fmax, double &gamma_gabor, int &level, int &orientation)
{
	//set default values, all these default parameters are from "https://www.mathworks.com/matlabcentral/fileexchange/44630-gabor-feature-extraction"
	unsigned int ImageDimension = atoi(argv[1]);
	fn_input = NULL;
	fn_output = NULL;
	fn_mask4Process = NULL;
	winSZ = 20; slidingDist = 500;
	R = 16; fmax = 0.25; gamma_gabor = sqrt(2); 
	level = 4; orientation = 8;

	//parse
	fn_input = argv[2];
	fn_output = argv[3];


	int myargc = argc - 4;
	char** myargv = argv + 4;
	while ((myargc > 1) && (myargv[0][0] == '-')){
		switch (myargv[0][1])
		{
		case 't':
			fn_mask4Process = &myargv[1][0];
			break;
		case 'w':
			winSZ = atoi(&myargv[1][0]);
			break;
		case 's':
			slidingDist = atoi(&myargv[1][0]);
			break;
		case 'R':
			R = atoi(&myargv[1][0]);
			break;
		case 'f':
			fmax = atof(&myargv[1][0]);
			break;
		case 'g':
			gamma_gabor = atof(&myargv[1][0]);
			break;
		case 'l':
			level = atoi(&myargv[1][0]);
			break;
		case 'o':
			orientation = atoi(&myargv[1][0]);
			break;
		}
		--myargc; --myargc;
		++myargv; ++myargv;
	}
	//std::cout << "winsize" << winSZ << "," << slidingDist;
}


//gabor mask
vector<vector <double>> gabor_wavelet(int argc, char* argv[], double R0, double fmax0, double gamma_gabor0, int level0, int orientation0)
{
	double fu, alpha, tetav;
	fu = fmax0 / pow(gamma_gabor0, level0 - 1);
	alpha = fu / gamma_gabor0;
	tetav = (orientation0 - 1 + 0.0) / 8 * PI;
	double xprime, yprime;
	vector<vector <double>> realGW;
	vector <double> realtemp;
	for (double i = -R0 / 2 + 1; i <= R0 / 2; i++)
	{
		realtemp.clear();
		double realpart, imagpart;
		for (double j = -R0 / 2 + 1; j <= R0 / 2; j++)
		{
			xprime = (i - 0.5)*cos(tetav) + (j - 0.5)*sin(tetav);
			yprime = -(i - 0.5)*sin(tetav) + (j - 0.5)*cos(tetav);
			realpart = fu*fu / (PI*gamma_gabor0*gamma_gabor0)*exp(-alpha*alpha*(xprime*xprime + yprime*yprime))*cos(2 * PI*fu*xprime);
			imagpart = fu*fu / (PI*gamma_gabor0*gamma_gabor0)*exp(-alpha*alpha*(xprime*xprime + yprime*yprime))*sin(2 * PI*fu*xprime);
			realtemp.push_back(realpart);

		}
		realGW.push_back(realtemp);;
	}
	return realGW;
}