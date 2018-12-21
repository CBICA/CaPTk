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
#define PI 3.141592653589793



#define isnan(x) ((x) != (x))

using namespace std;

//predefine functions

void usage(int argc, char *argv[]);
void parse(int argc, char *argv[], char* &fn_input, char* &fn_output, char* &fn_mask4Process, unsigned int &winSZ, unsigned int &slidingDist, unsigned int &neighbours, unsigned int &radius, unsigned int &LBPstyle);


//template <unsigned int ImageDimension>
int CalculateLocalBinaryPattern(int argc, char *argv[]);
vector<vector<int>> index_matrix(int n1, int n2);

unsigned int rotateLeft(unsigned int i, unsigned int samples);
int NumberOfSetBits(int i);
vector<int> LUT(unsigned int samples, int type);

//maind code

int main(int argc, char *argv[])
{

	usage(argc, argv);

	CalculateLocalBinaryPattern(argc, argv);
	
}
//only defined in 2D space
const unsigned int ImageDimension = 2;

int CalculateLocalBinaryPattern(int argc, char *argv[])
{
	//parse inputs
	char* inputImageFn = NULL; char* fn_output = NULL;
	char* fn_mask4Process = NULL;
	unsigned int winSZ; unsigned int slidingDist;
	unsigned int neighbours; unsigned int radius; unsigned int LBPstyle;

	parse(argc, argv, inputImageFn, fn_output, fn_mask4Process, winSZ, slidingDist, neighbours, radius, LBPstyle);

	//Image type defination
	typedef double PixelType;
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
	double outputImageSpacing[ImageDimension];
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
	vector<vector <double>> all_lbp;
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
		//std::cout << pos << endl;
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
		vector<vector <double>> image_matrix;
		vector<double> temp2;
		for (int j = 0; j < winSZ; j++)
		{
			temp2.clear();
			for (int i = 0; i < winSZ; i++)
			{
				ImageType::IndexType index;
				index[0] = i;
				index[1] = j;
				//cout << image->GetPixel(index) << endl;
				temp2.push_back(inputImageRoiIm->GetPixel(index));
			}
			image_matrix.push_back(temp2);
		}

		////////////////////////////////////////////////////////////////////////////////
		//compute features here Note: ItI and ItM proved above
		////////////////////////////////////////////////////////////////////////////////
		int lbp_matrix_size = winSZ - 2 * radius;
		//int lbp_pattern[61][61] = {0};
		int** lbp_pattern = new int*[lbp_matrix_size];
		for (int i = 0; i < lbp_matrix_size; ++i)
			lbp_pattern[i] = new int[lbp_matrix_size]();


		for (auto n = 0; n<neighbours; n++) 
		{
			// sample points
			double x = static_cast<double>(radius)* cos(2.0*PI*n / static_cast<double>(neighbours));
			double y = static_cast<double>(radius)* -sin(2.0*PI*n / static_cast<double>(neighbours));
			// relative indices
			int fx = static_cast<int>(floor(x));
			int fy = static_cast<int>(floor(y));
			int cx = static_cast<int>(ceil(x));
			int cy = static_cast<int>(ceil(y));
			// fractional part
			double ty = y - fy;
			double tx = x - fx;
			// set interpolation weights
			double w1 = (1 - tx) * (1 - ty);
			double w2 = tx  * (1 - ty);
			double w3 = (1 - tx) *      ty;
			double w4 = tx  *      ty;
			// iterate through your data
			for (auto i = radius; i < winSZ - radius; i++) 
			{
				for (auto j = radius; j < winSZ - radius; j++)
				{
					double t = w1*image_matrix[i + fy][j + fx] + w2*image_matrix[i + fy][j + cx] + w3*image_matrix[i + cy][j + fx] + w4*image_matrix[i + cy][j + cx];
					lbp_pattern[i - radius][j - radius] += ((t > image_matrix[i][j]) && (abs(t - image_matrix[i][j]) > std::numeric_limits<double>::epsilon())) << n;
					//	std::cout << temp.at<float>(0, i - radius, j - radius);
				}
			}

		}

		vector<int> table = LUT(neighbours, LBPstyle);
		int hist_size = table.back();
		table.pop_back();
		//float *hist_lbp = new float[hist_size]{0};
		vector<double> hist_lbp(hist_size, 0.0);
		//cout << hist_size << endl;
		for (auto i = 0; i < lbp_matrix_size; i++){
			for (auto j = 0; j < lbp_matrix_size; j++){
				//	std::cout << temp.at<float>(0, i, j);
				int temp;
				int a = lbp_pattern[i][j];
				temp = table[lbp_pattern[i][j]];
				hist_lbp[temp] = hist_lbp[temp] + 1;
				//	std::cout << temp1.at<float>(0, i, j);
			}
		}
		all_lbp.push_back(hist_lbp);
		for (int i = 0; i < winSZ - 2 * radius; ++i) 
		{
			delete[] lbp_pattern[i];
		}
		delete[] lbp_pattern;

	}

	
	for (int i = 0; i < all_lbp[0].size(); i++)
	{
		char fullname[300];
		char *buffer1;
		//string path = "C:/folder/branches/Libracxx/data/New_folder_power/";
		strcpy(fullname, fn_output);
		strcat(fullname, "_");
		string s1 = std::to_string(i);
		char const *pchar1 = s1.c_str();
		strcat(fullname, pchar1);
		strcat(fullname, ".nii");
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
		result_image->Allocate();
		result_image->FillBuffer(itk::NumericTraits< double >::Zero);

		for (int k = 0; k < result_size[0]; k++)
		{
			for (int l = 0; l < result_size[1]; l++)
			{
				ImageType::IndexType temp_index;
				temp_index[0] = k;
				temp_index[1] = l;
				result_image->SetPixel(temp_index, all_lbp[l*result_size[0] + k][i]);
			}
		}

		
		typedef itk::ImageFileWriter< ImageType  > WriterType;
		WriterType::Pointer writer = WriterType::New();
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

void parse(int argc, char *argv[], char* &fn_input, char* &fn_output, char* &fn_mask4Process, unsigned int &winSZ, unsigned int &slidingDist, unsigned int &neighbours, unsigned int &radius, unsigned int &LBPstyle)
{
	//set default values
	unsigned int ImageDimension = atoi(argv[1]);
	fn_input = NULL;
	fn_output = NULL;
	fn_mask4Process = NULL;
	winSZ = 63; slidingDist = 63;
	neighbours = 8;
	radius = 1;
	LBPstyle = 2;

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
		case 'n':
			neighbours = atoi(&myargv[1][0]);
			break;
		case 'r':
			radius = atoi(&myargv[1][0]);
			break;
		case 'L':
			LBPstyle = atoi(&myargv[1][0]);
			break;
		}
		--myargc; --myargc;
		++myargv; ++myargv;
	}
}

unsigned int rotateLeft(unsigned int i, unsigned int samples) {
	unsigned int bg = ((i & (1 << (samples - 1))) >> (samples - 1)); // bitget(r,samples)
	unsigned int bs = (i << 1) & ((int)pow(2., (int)samples) - 1); // bitshift(r, 1, samples)
	unsigned int j = (bs + bg) & ((int)pow(2., (int)samples) - 1); // bitset( bs, 1, bg )
	return j;
}

int NumberOfSetBits(int i) {
	i = i - ((i >> 1) & 0x55555555);
	i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
	return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

//samples represent the number of neighbors, default = 8;
//case 0: original LBP; case 1: uniform LBP; case 2: rotation invariant LBP; case 3: uniform + rotation invariant LBP; default = 2
//last number of this return vector is the number of total patterns for each case;
//paper reference: Rotated Local Binary Pattern (RLBP): Rotation invariant texture descriptor
vector<int> LUT(unsigned int samples, int type)  //LUT: LBP lookup table, it create the pattern code of LBP pattern, different types have different patterns;
{
	vector<int> table;
	int newMax = 0; //number of patterns in the resulting LBP code
	int index = 0;
	if (type == 0) 
	{
		newMax = static_cast<int>(pow(2., (int)samples));
		for (int i = 0; i < newMax; i++) 
		{
			table.push_back(i);
		}
	}

	else if (type == 1) 
	{
		// Uniform 2
		newMax = samples * (samples - 1) + 3;
		for (unsigned int i = 0; i < pow(2., (int)(samples)); i++) 
		{
			unsigned int j = rotateLeft(i, samples);
			int numt = NumberOfSetBits(i ^ j);
			if (numt <= 2)
			{
				table.push_back(index);
				index = index + 1;
			}
			else 
			{
				table.push_back(newMax - 1);
			}
		}
	}
	else if (type == 2) 
	{
		long N = (int)pow(2., (int)samples);
		// Rotation Invariant
		int * tmpMap = new int[N];
		memset((void *)tmpMap, -1, N);

		for (long i = 0; i < N; i++) 
		{
			tmpMap[i] = -1;

			unsigned long rm = i;
			unsigned long r = i;
			for (unsigned int j = 1; j <= samples - 1; j++) 
			{
				r = rotateLeft(r, samples);
				if (r < rm)
					rm = r;
			}
			if (tmpMap[rm] < 0) 
			{
				tmpMap[rm] = newMax;
				newMax = newMax + 1;
			}
			table.push_back(tmpMap[rm]);
			//		std::cout << table->at(rm);
		}

	}
	else if (type == 3) 
	{
		// Rotation invariant uniform 2
		newMax = samples + 2;
		for (unsigned int i = 0; i <= pow(2., (int)samples) - 1; i++) 
		{
			unsigned int j = rotateLeft(i, samples); //bitset( bitshift( i, 1, samples ), 1, bitget( i, samples ) ); // rotate left
			unsigned int numt = NumberOfSetBits(i ^ j); //sum(bitget(bitxor(i,j),1:samples));
			if (numt <= 2)
			{
				table.push_back(NumberOfSetBits(i));
				//		std::cout << table->at(i);
			}
			else
			{
				table.push_back(samples + 1);
				//	std::cout << table->at(i);
			}
		}
	}
	table.push_back(newMax);
	return table;
}