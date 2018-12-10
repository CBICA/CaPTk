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



#define isnan(x) ((x) != (x))

using namespace std;

//predefine functions

void usage(int argc, char *argv[]);
void parse(int argc, char *argv[], char* &fn_input, char* &fn_output, char* &fn_mask4Process, unsigned int &winSZ, unsigned int &slidingDist);


//template <unsigned int ImageDimension>
int CalculatePowerSpectrum(int argc, char *argv[]);
vector<vector<int>> index_matrix(int n1, int n2);

//maind code

int main(int argc, char *argv[])
{

	usage(argc, argv);
	
	CalculatePowerSpectrum(argc, argv);
	
}

const unsigned int ImageDimension = 2;
//template <unsigned int ImageDimension>
int CalculatePowerSpectrum(int argc, char *argv[])
{
	//parse inputs
	char* inputImageFn = NULL; char* fn_output = NULL;
	char* fn_mask4Process = NULL;
	unsigned int winSZ; unsigned int slidingDist;

	parse(argc, argv, inputImageFn, fn_output, fn_mask4Process, winSZ, slidingDist);

	//Image type defination
	typedef float PixelType;
	typedef itk::Image<PixelType, ImageDimension> ImageType;
	typedef itk::ImageFileReader<ImageType> ReaderType;

	//input image
	ReaderType::Pointer imageReader = ReaderType::New();
	imageReader->SetFileName(inputImageFn);
	imageReader->Update();
	ImageType::Pointer inputImage = imageReader->GetOutput();
	ImageType::SizeType inputImageSize = inputImage->GetLargestPossibleRegion().GetSize();
	itk::ImageRegionIterator<ImageType> ItIGlobal(inputImage, inputImage->GetLargestPossibleRegion());




	// Michael's modification
	ImageType::DirectionType dirCosine = inputImage->GetDirection();
	ImageType::SpacingType inputImageSpacing = inputImage->GetSpacing();

	//mask4Process
	FILE *fp = NULL;
	if (fn_mask4Process != NULL)
		fp = fopen(fn_mask4Process, "r");


	ReaderType::Pointer mask4ProcessReader = ReaderType::New();
	ImageType::Pointer mask4Process = ImageType::New();
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
	itk::ImageRegionIterator<ImageType> ItMGlobal(mask4Process, mask4Process->GetLargestPossibleRegion());

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
	MaskImageType::Pointer maskRoiIm = MaskImageType::New();
	MaskImageType::PointType origin;
	origin.Fill(0.0);
	maskRoiIm->SetOrigin(origin);
	maskRoiIm->SetSpacing(inputImage->GetSpacing());
	MaskImageType::SizeType size;
	size.Fill(winSZ);
	MaskImageType::RegionType region;
	region.SetSize(size);
	maskRoiIm->SetRegions(region);
	maskRoiIm->Allocate();
	maskRoiIm->FillBuffer(itk::NumericTraits<PixelType>::One);
	itk::ImageRegionIterator<MaskImageType> ItM(maskRoiIm, maskRoiIm->GetLargestPossibleRegion());


	//compute features for each pixel
	ImageType::SizeType imsz; imsz.Fill(winSZ);
	long pos = -1;
	int count_lattice = 0;
	vector <double> all_power_spectrum;
	for (ItIGlobal.GoToBegin(), ItMGlobal.GoToBegin(); !ItIGlobal.IsAtEnd(); ++ItIGlobal, ++ItMGlobal)
	{
		ImageType::IndexType thisidx = ItIGlobal.GetIndex();

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
		//cout << pos << endl;
		//if in mask
		ImageType::IndexType centidx = thisidx;
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
		ImageROIType::Pointer imageROI = ImageROIType::New();
		ImageType::IndexType imstart = thisidx;

		ImageType::RegionType desiredImageRegion;
		desiredImageRegion.SetSize(imsz);
		desiredImageRegion.SetIndex(imstart);
		imageROI->SetRegionOfInterest(desiredImageRegion);
		imageROI->SetInput(inputImage);
		imageROI->Update();
		ImageType::Pointer inputImageRoiIm = imageROI->GetOutput();

		itk::ImageRegionIteratorWithIndex<ImageType> ItI(inputImageRoiIm,
			inputImageRoiIm->GetLargestPossibleRegion());

		//Yifan's modification


		////////////////////////////////////////////////////////////////////////////////
		//compute features here Note: ItI and ItM proved above
		////////////////////////////////////////////////////////////////////////////////
		typedef itk::Image< PixelType, 2 > RealImageType;
		typedef itk::WrapPadImageFilter< RealImageType, RealImageType > PadFilterType;
		PadFilterType::Pointer padFilter = PadFilterType::New();
		padFilter->SetInput(inputImageRoiIm);
		PadFilterType::SizeType padding;
		// Input size is [48, 62, 42].  Pad to [48, 64, 48].
		RealImageType::SizeType size00 = inputImageRoiIm->GetLargestPossibleRegion().GetSize();
		vector <int> reference;
		//create an array with numbers are a multiplication by 2,3 and 5 because the itk filter can only computed on a window size with multiplication by 2,3 and 5; it depends on our lattice size;
		for (int i = 0; i < 10; i++)
		{
			for (int j = 0; j < 7; j++)
			{
				for (int k = 0; k < 5; k++)
				{
					int s = pow(2, i)*pow(3, j)*pow(5, k);
					if (s < 1000)
					{
						reference.push_back(s);
					}

				}
			}
		}
		sort(reference.begin(), reference.end());
		//padding the window into the size in the array
		int kk;
		for (int i = 0; i < reference.size(); i++)
		{
			kk = reference[i] - size00[0];
			if (kk >= 0) break;
		}
		padding[0] = kk;
		padding[1] = kk;
		//padding[2] = 0;
		padFilter->SetPadUpperBound(padding);
		padFilter->Update();
		RealImageType::RegionType region0 = padFilter->GetOutput()->GetLargestPossibleRegion();
		RealImageType::SizeType size0 = region0.GetSize();

		//use itk filter to calculate power spectrum: https://itk.org/ITKExamples/src/Filtering/FFT/ComputeImageSpectralDensity/Documentation.html
		typedef itk::ForwardFFTImageFilter< RealImageType > ForwardFFTFilterType;
		typedef ForwardFFTFilterType::OutputImageType ComplexImageType;
		ForwardFFTFilterType::Pointer forwardFFTFilter = ForwardFFTFilterType::New();
		forwardFFTFilter->SetInput(padFilter->GetOutput());
		forwardFFTFilter->Update();

		typedef itk::ComplexToModulusImageFilter< ComplexImageType, RealImageType >
			ComplexToModulusFilterType;
		ComplexToModulusFilterType::Pointer complexToModulusFilter
			= ComplexToModulusFilterType::New();
		complexToModulusFilter->SetInput(forwardFFTFilter->GetOutput());
		complexToModulusFilter->Update();
		typedef unsigned short                           OutputPixelType;
		typedef itk::Image< OutputPixelType, 2 > OutputImageType;
		typedef itk::IntensityWindowingImageFilter< RealImageType, OutputImageType >
			WindowingFilterType;
		WindowingFilterType::Pointer windowingFilter
			= WindowingFilterType::New();
		windowingFilter->SetInput(complexToModulusFilter->GetOutput());
		windowingFilter->SetWindowMinimum(0);
		windowingFilter->SetWindowMaximum(100000);
		windowingFilter->Update();

		typedef itk::FFTShiftImageFilter< OutputImageType, OutputImageType > FFTShiftFilterType;
		FFTShiftFilterType::Pointer fftShiftFilter = FFTShiftFilterType::New();
		fftShiftFilter->SetInput(windowingFilter->GetOutput());
		fftShiftFilter->Update();

		OutputImageType::Pointer power_image = fftShiftFilter->GetOutput();
		OutputImageType::RegionType region = power_image->GetLargestPossibleRegion();
		OutputImageType::SizeType size = region.GetSize();
		//cout << size[0] << "," << size[1] << endl;
		//float power_number = min(size[0], size[1])/2;
		float power_number = ((winSZ / 2) / 10) * 10;
		vector <int> count;
		vector <unsigned short> sum_power_density;
		for (int i = 0; i < power_number + 1; i++)
		{
			sum_power_density.push_back(0);
			count.push_back(0);
		}

		vector<vector <int>> refer_matrix = index_matrix(size[0], size[1]);

		//compute average power correspond to different radius; paper: Modelling the Power Spectra of Natural Images: Statistics and Information
		for (int i = 0; i < size[0]; i++)
		{
			for (int j = 0; j < size[1]; j++)
			{
				RealImageType::IndexType index;
				index[0] = i;
				index[1] = j;
				//cout << fftShiftFilter->GetOutput()->GetPixel(index) << endl;
				if (refer_matrix[i][j]>power_number)
				{
					count[power_number]++;
					sum_power_density[power_number] += fftShiftFilter->GetOutput()->GetPixel(index);
				}
				else
				{
					count[refer_matrix[i][j]]++;
					sum_power_density[refer_matrix[i][j]] += fftShiftFilter->GetOutput()->GetPixel(index);
				}
			}
		}
		//compute the average on different radius
		double sum_y = 0;
		double sum_x = 0;
		double sum_y_x = 0;
		double sum_x_x = 0;
		double n = power_number;
		for (int i = 1; i < n + 1; i++)
		{
			double power_y = (sum_power_density[i] + 0.0) / (count[i] + 0.0);
			//aver_power.push_back((sum_power_density[i] + 0.0) / (count[i] + 0.0));
			//dist.push_back(i);
			sum_y_x += log(power_y)*log(i);
			sum_y += log(power_y);
			sum_x_x += log(i)*log(i);
			sum_x += log(i);
		}
		//fit the line and get the slope
		double beta;
		beta = (sum_y_x - sum_y*sum_x / n) / (sum_x_x - sum_x*sum_x / n);
		all_power_spectrum.push_back(-beta);

	}
	//total window*25*5 matrix     first dim: different masks; second dim: different masks; third dim: 5 laws measures;    In total: 125 features e.g. mask1_mean; mask15_std ...
	string filename = "_power_spectrum.nii";
	string path = fn_output + filename;
	char *cstrname = new char[path.length() + 1];
	strcpy(cstrname, path.c_str());
	//float *temp = new float[nPixValNum];

	ImageType::Pointer result_image = ImageType::New();
	ImageType::RegionType result_region;
	ImageType::SizeType result_size;
	ImageType::IndexType result_index;
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
			ImageType::IndexType temp_index;
			temp_index[0] = k;
			temp_index[1] = l;
			result_image->SetPixel(temp_index, all_power_spectrum[l*result_size[0] + k]);
		}
	}
	//for (int j = 0; j < nPixValNum; j++)
	//	temp[j] = all_laws_features[j][i];
	//writeImage<float, ImageDimension, 1, ImageType>(cstrname, temp, outputImageSize, outputImageSpacing, dirCosine);
	typedef itk::ImageFileWriter< ImageType  > WriterType;
	WriterType::Pointer writer = WriterType::New();
	cout << cstrname << endl;
	writer->SetFileName(cstrname);
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
	//if (!writeImage<float, ImageDimension, 1, ImageType>(cstrname, temp, outputImageSize, outputImageSpacing, dirCosine))
	//	delete[] temp;
	delete[] cstrname;



	//save results

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

void parse(int argc, char *argv[], char* &fn_input, char* &fn_output, char* &fn_mask4Process, unsigned int &winSZ, unsigned int &slidingDist)
{
	//set default values
	unsigned int ImageDimension = atoi(argv[1]);
	fn_input = NULL;
	fn_output = NULL;
	fn_mask4Process = NULL;
	winSZ = 63; slidingDist = 63;

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
		}
		--myargc; --myargc;
		++myargv; ++myargv;
	}
}

vector<vector<int>> index_matrix(int n1, int n2)
{
	vector<vector<int>> index;
	int center[2] = { floor(n1 / 2) + 1, floor(n2 / 2) + 1 };
	for (float i = 0; i < n1; i = i + 1)
	{
		vector<int> temp;
		temp.clear();
		for (float j = 0; j < n2; j = j + 1)
		{
			temp.push_back(max(abs(i - center[0]), abs(j - center[1])));
			//cout << ceil(sqrt(i*i + j*j)) << ",";
		}
		//cout << endl;
		index.push_back(temp);
	}
	return index;
}
