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
void parse(int argc, char *argv[], char* &fn_input, char* &fn_output, char* &fn_mask4Process, unsigned int &winSZ, unsigned int &slidingDist, double &Eta, double &epsi, int &radius_edge);


template <unsigned int ImageDimension>
int CalculateEdgeEnhance(int argc, char *argv[]);
vector<vector<int>> index_matrix(int n1, int n2);
vector<vector <double>> gD(vector<vector <double>> f, int scale, int ox, int oy);
vector <double> gDerivative(int a, vector <double> x, vector <double> Gs, int scale);

//maind code

int main(int argc, char *argv[])
{

	usage(argc, argv);
	switch (atoi(argv[1]))
	{
		//cout << 1 << endl;
	case 2:
		CalculateEdgeEnhance<2>(argc, argv);
		//std::cout << r;
		break;
	case 3:
		CalculateEdgeEnhance<3>(argc, argv);
		break;
	default:
		std::cerr << "Unsupported dimension" << std::endl;
		exit(EXIT_FAILURE);
	}
	return 0;
}

//const unsigned int ImageDimension = 2;
template <unsigned int ImageDimension>
int CalculateEdgeEnhance(int argc, char *argv[])
{

	//parse inputs
	char* inputImageFn = NULL; char* fn_output = NULL;
	char* fn_mask4Process = NULL;
	unsigned int winSZ;	unsigned int slidingDist;
	double Eta, epsi;
	int radius_edge;

	parse(argc, argv, inputImageFn, fn_output, fn_mask4Process, winSZ, slidingDist, Eta, epsi, radius_edge);

	//std::cout << "winsize" << winSZ << "," << slidingDist;
	//Image type defination
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
	vector<double> edge_enhance;
	//vector<double> lattice_feature;
	for (ItIGlobal.GoToBegin(), ItMGlobal.GoToBegin(); !ItIGlobal.IsAtEnd(); ++ItIGlobal, ++ItMGlobal)
	{
		typename ImageType::IndexType thisidx = ItIGlobal.GetIndex();

		//cancel operations for pixels not in mask, and that will not be slided
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
		//all the details are in this report: https://staff.fnwi.uva.nl/r.vandenboomgaard/nldiffusionweb/nldiffusioncode.pdf
		//all of this algorithm are translated from matlab version;
		typename ImageType::SizeType temp_size;
		temp_size = inputImageRoiIm->GetLargestPossibleRegion().GetSize();
		vector<vector <double>> f;
		vector<double> temp2;
		for (int j = 0; j < temp_size[1]; j++)
		{
			temp2.clear();
			for (int i = 0; i < temp_size[0]; i++)
			{
				typename ImageType::IndexType index;
				index[0] = i;
				index[1] = j;
				temp2.push_back(inputImageRoiIm->GetPixel(index));
			}
			f.push_back(temp2);
		}
		//double Eta, eps;
		vector<vector <double>> fsx = gD(f, radius_edge, 1, 0);
		vector<vector <double>> fsy = gD(f, radius_edge, 0, 1);
		vector<vector <double>> Cindex;
		vector<double> temp;
		for (int i = 0; i < fsx.size(); i++)
		{
			temp.clear();
			for (int j = 0; j < fsx[0].size(); j++)
			{
				double beta2 = exp(-(fsx[i][j] * fsx[i][j] + fsy[i][j] * fsy[i][j]) / (Eta*Eta));
				double beta1 = 0.2 * beta2;
				double factor = sqrt(fsx[i][j] * fsx[i][j] + fsy[i][j] * fsy[i][j]);
				double A = (beta1 * fsx[i][j] * fsx[i][j] + beta2 * fsy[i][j] * fsy[i][j])/factor;
				double B = ((beta2 - beta1) * fsx[i][j] * fsy[i][j])/factor;
				double C = B;
				double D = (beta2 * fsx[i][j] * fsx[i][j] + beta1 * fsy[i][j] * fsy[i][j])/factor;
				temp.push_back(((A - D)*(A - D) + 4 * B*C) / ((A + D + epsi)*(A + D + epsi)));
			}
			Cindex.push_back(temp);
		}

		////////////////////////////////////////////////////////////////////////////////
		//compute features here Note: ItI and ItM proved above
		////////////////////////////////////////////////////////////////////////////////		
		double sum_final = 0;
		for (int i = 0; i < Cindex.size(); i++)
		{
			for (int j = 0; j < Cindex[0].size(); j++)
			{
				sum_final += Cindex[i][j];
			}
		}
		//laws_feature_each_lattice.push_back(mask_feature);
		//stats->SetInput(inputImageRoiIm);
		double average = sum_final / (Cindex.size()*Cindex[0].size());
		edge_enhance.push_back(average);
		//cout << count_lattice << endl;

	}
	
	char fullname[300];
	strcpy(fullname, fn_output);
	strcat(fullname, "_radius_");
	string s1 = std::to_string(radius_edge);
	char const *pchar1 = s1.c_str();
	strcat(fullname, pchar1);
	strcat(fullname, "_edge_enhance.nii");
	
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
			result_image->SetPixel(temp_index, edge_enhance[l*result_size[0] + k]);
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

void parse(int argc, char *argv[], char* &fn_input, char* &fn_output, char* &fn_mask4Process, unsigned int &winSZ, unsigned int &slidingDist, double &Eta, double &epsi, int &radius_edge)
{
	//set default values
	unsigned int ImageDimension = atoi(argv[1]);
	fn_input = NULL;
	fn_output = NULL;
	fn_mask4Process = NULL;
	winSZ = 63; slidingDist = 63;
	Eta = 10; epsi = 10; radius_edge = 5;

	//parse
	fn_input = argv[2];
	fn_output = argv[3];

	//std::cout << "winsize" << winSZ << "," << slidingDist;

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
		case 'r':
			radius_edge = atoi(&myargv[1][0]);
			break;
		case 'e':
			epsi= atof(&myargv[1][0]);
			break;
		case 'E':
			Eta = atof(&myargv[1][0]);
			break;
		}
		--myargc; --myargc;
		++myargv; ++myargv;
	}
	//std::cout << "winsize" << winSZ << "," << slidingDist;
}


vector<vector <double>> gD(vector<vector <double>> f, int scale, int ox, int oy)
{
	int J = ceil(3*scale);
	vector <double> x, Gs;
	double sum = 0;
	for (int i = -J; i < J + 1; i++)
	{
		x.push_back(i);
		Gs.push_back(exp(-i*i/(2*scale*scale)));
		sum += exp(-i*i / (2 * scale*scale));
	}
	for (int i = 0; i < 2 * J + 1; i++)
	{
		Gs[i] = Gs[i] / sum;
	}
	vector <double> Gsx = gDerivative(ox, x, Gs, scale);
	vector <double> Gsy = gDerivative(oy, x, Gs, scale);
	int N = f.size();
	int M = f[0].size();
	int K = (Gsx.size() - 1) / 2;
	int L = (Gsy.size() - 1) / 2;
	vector<int> iind;
	for (int i = 0; i < N + 2 * K; i++)
	{
		iind.push_back(min(max(i+1-K,0),N-1));
	}
	vector<int> jind;
	for (int i = 0; i < M + 2 * L; i++)
	{
		jind.push_back(min(max(i + 1 - L, 0), M-1));
	}
	vector<vector <double>> fwb;
	vector<double> temp;
	for (int i = 0; i < iind.size(); i++)
	{
		temp.clear();
		for (int j = 0; j < jind.size(); j++)
		{
			temp.push_back(f[iind[i]][jind[j]]);
		}
		fwb.push_back(temp);
	}
	vector<vector <double>> filter;
	for (int i = 0; i < Gsx.size(); i++)
	{
		temp.clear();
		for (int j = 0; j < Gsy.size(); j++)
		{
			temp.push_back(Gsx[i]*Gsy[j]);
		}
		filter.push_back(temp);
	}
	vector<vector <double>> conv_final;
	for (int i = 0; i < iind.size() - Gsx.size() + 1; i++)
	{
		temp.clear();
		for (int j = 0; j < jind.size() - Gsy.size() + 1; j++)
		{
			double sum = 0;
			for (int k = 0; k < Gsx.size(); k++)
			{
				for (int l = 0; l < Gsy.size(); l++)
				{
					sum = sum + filter[k][l] * fwb[i + k][j + l];
				}
			}
			temp.push_back(sum);
		}
		conv_final.push_back(temp);
	}
	return conv_final;
}


vector <double> gDerivative(int a, vector <double> x, vector <double> Gs, int scale)
{
	vector <double> r;
	for (int i = 0; i < x.size(); i++)
	{
		if (a == 0)
		{
			r.push_back(Gs[i]);
		}
		else 
		{
			r.push_back(-x[i]/(scale*scale)*Gs[i]);
		}
	}
	return r;
}