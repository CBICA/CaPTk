/*=========================================================================
Program:   ITK General
Language:  C++
Date:      $Date: 2017/09/06 16:41:20 $
ITK Version:   InsightToolkit-4.12.0

Author: Yifan Hu (huyifan1989@gmail.com)
Institution: CBIG @ UPenn

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


#define isnan(x) ((x) != (x))

using namespace std;

//predefine functions

void usage(int argc, char *argv[]);
void parse(int argc, char *argv[], char* &fn_input, char* &fn_output, char* &fn_mask4Process, unsigned int &winSZ, unsigned int &slidingDist);


template <unsigned int ImageDimension>
int CalculateLawsStatistics(int argc, char *argv[]);
vector<vector<vector <int>>> laws_masks();
vector<vector<vector <int>>> laws_masks_3();

//maind code

int main(int argc, char *argv[])
{

	usage(argc, argv);
	switch (atoi(argv[1]))
	{
	case 2:
		CalculateLawsStatistics<2>(argc, argv);
		break;
	case 3:
		CalculateLawsStatistics<3>(argc, argv);
		break;
	default:
		std::cerr << "Unsupported dimension" << std::endl;
		exit(EXIT_FAILURE);
	}
}

template <unsigned int ImageDimension>
int CalculateLawsStatistics(int argc, char *argv[])
{

	//parse inputs
	char* inputImageFn = NULL; char* fn_output = NULL;
	char* fn_mask4Process = NULL;
	unsigned int winSZ; unsigned int slidingDist;

	parse(argc, argv, inputImageFn, fn_output, fn_mask4Process, winSZ, slidingDist);

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

	// Yifan's modification
	vector<vector <double>> all_laws_features(nPixValNum, vector<double>(125, 0));

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
		//cout << pos;
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
		typename ImageType::SizeType temp_size;		
		temp_size = inputImageRoiIm->GetLargestPossibleRegion().GetSize();
		vector<vector <double>> image_matrix;
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
			image_matrix.push_back(temp2);
		}
		

			////////////////////////////////////////////////////////////////////////////////
			//compute features here Note: ItI and ItM proved above
			////////////////////////////////////////////////////////////////////////////////
		typename ImageType::Pointer result_matrix = ImageType::New();
		typename ImageType::RegionType region_new;
		typename ImageType::IndexType start;
		start.Fill(0);
		typename ImageType::SizeType size_new_matrix;
		size_new_matrix[0] = winSZ - 5 + 1;
		size_new_matrix[1] = winSZ - 5 + 1;
		region_new.SetSize(size_new_matrix);
		result_matrix->SetRegions(region_new);
		result_matrix->Allocate();
		result_matrix->FillBuffer(itk::NumericTraits< double >::Zero);
		//here we can also choose mask size=3 by using laws_masks_3 function
		vector<vector<vector <int>>> mask_matrix = laws_masks();
		int mask_num = mask_matrix.size();
		typename ImageType::IndexType pixelIndex;
		//25 masks(9 masks when using mask size=3)
		//feature for each mask
		for (int m = 0; m < mask_num; m++)
		{
			int i, j;
			//convolution step
			for (i = 0; i < size_new_matrix[1]; i++)
			{
				for (j = 0; j < size_new_matrix[0]; j++)
				{

					double sum = 0;
					for (int k = 0; k < 5; k++)
					{
						for (int l = 0; l < 5; l++)
						{
							sum = sum + mask_matrix[m][k][l] * image_matrix[i + k][j + l];
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
			double origin_mean = statisticsImageFilter->GetMean();
			double origin_sigma = statisticsImageFilter->GetSigma();
			//mean and standard deviation
			all_laws_features[pos][5 * m] = origin_mean;
			all_laws_features[pos][5 * m + 1] = origin_sigma;
			typedef itk::AddImageFilter <ImageType, ImageType, ImageType> AddImageFilterType;
			typename AddImageFilterType::Pointer addImageFilter = AddImageFilterType::New();
			addImageFilter->SetInput(result_matrix);
			addImageFilter->SetConstant2(-statisticsImageFilter->GetMean());
			addImageFilter->Update();
			typename ImageType::Pointer result_matrix1 = addImageFilter->GetOutput();
			//power of value: 3 and 4
			typedef itk::PowImageFilter< ImageType, ImageType, ImageType> PowImageFilterType;
			typename PowImageFilterType::Pointer powImageFilter1 = PowImageFilterType::New();
			typename PowImageFilterType::Pointer powImageFilter2 = PowImageFilterType::New();
			//skewness
			powImageFilter1->SetInput(result_matrix1);
			powImageFilter1->SetConstant2(3);
			powImageFilter1->Update();
			typename ImageType::Pointer skewness_matrix = powImageFilter1->GetOutput();
			statisticsImageFilter->SetInput(skewness_matrix);
			statisticsImageFilter->Update();
			if (origin_sigma==0)
				all_laws_features[pos][5 * m + 2] = 0;
			else
				all_laws_features[pos][5 * m + 2] = (statisticsImageFilter->GetMean()) / origin_sigma;
			//kurtosis
			powImageFilter2->SetInput(result_matrix1);
			powImageFilter2->SetConstant2(4);
			powImageFilter2->Update();
			typename ImageType::Pointer kurtosis_matrix = powImageFilter2->GetOutput();
			statisticsImageFilter->SetInput(kurtosis_matrix);
			statisticsImageFilter->Update();
			if (origin_sigma == 0)
				all_laws_features[pos][5 * m + 3] = 0;
			else
				all_laws_features[pos][5 * m + 3] = (statisticsImageFilter->GetMean()) / pow(origin_sigma, 4) - 3;
			//entropy
			typename PowImageFilterType::Pointer powImageFilter3 = PowImageFilterType::New();
			powImageFilter3->SetInput(result_matrix);
			powImageFilter3->SetConstant2(2);
			powImageFilter3->Update();
			statisticsImageFilter->SetInput(powImageFilter3->GetOutput());
			statisticsImageFilter->Update();
			all_laws_features[pos][5 * m + 4] = statisticsImageFilter->GetMean();
		}
	}
	//total window*25*5 matrix     first dim: different masks; second dim: different masks; third dim: 5 laws measures;    In total: 125 features e.g. mask1_mean; mask15_std ...
	for (int i = 0; i < 125; i++)
	{
		char* measure[5] = { "_mean.nii", "_std.nii", "_skewness.nii", "_kurtosis.nii", "_entropy.nii"};
		char fullname[300];
		strcpy(fullname, fn_output);
		char buffer[2];
		int integer_part = i / 5 + 1;
		int residual = i % 5;
		strcat(fullname, "_mask_");
		string s1 = std::to_string(integer_part);
		char const *pchar1 = s1.c_str();
		strcat(fullname, pchar1);
		strcat(fullname, measure[residual]);

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
		result_image->FillBuffer(itk::NumericTraits< double >::Zero);

		for (int k = 0; k < result_size[0]; k++)
		{
			for (int l = 0; l < result_size[1]; l++)
			{
				typename ImageType::IndexType temp_index;
				temp_index[0] = k;
				temp_index[1] = l;
				result_image->SetPixel(temp_index, all_laws_features[l*result_size[0] + k][i]);
			}
		}

		typedef itk::ImageFileWriter< ImageType > WriterType;
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
	if (argc < 4)
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
//reference: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7724557
//this reference showed the best mask size, how to define masks and how to extract features
vector<vector<vector <int>>> laws_masks()
{
	int all_vectors[5][5] = { { 1, 4, 6, 4, 1 },        //Level detection
	{ -1, -2, 0, 2, 1 },                                //Edge detection
	{ -1, 0, 2, 0, -1 },                                //Spot detection
	{ -1, 2, 0, -2, 1 },                                //Wave detection
	{ 1, -4, 6, -4, 1 } };                              //Ripple detection
	
	//mask by level and edge could be considered as gradient
	vector<vector<vector<int>>> all_masks;

	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			vector<vector<int>> temp;
			int *vector1 = all_vectors[i];
			int *vector2 = all_vectors[j];
			for (int row = 0; row < 5; row++)
			{
				vector<int> temp_vector;
				for (int col = 0; col < 5; col++)
				{
					temp_vector.push_back(vector1[row] * vector2[col]);
				}
				temp.push_back(temp_vector);
				temp_vector.clear();
			}
			all_masks.push_back(temp);
		}
	}
	return all_masks;
}

vector<vector<vector <int>>> laws_masks_3()
{
	int all_vectors[3][3] = { { 1, 2, 1 },        //Level detection
	{ -1, 0, 1 },                                //Edge detection
	{ -1, 2, -1 }, };                            //Spot detection

	//mask by level and edge could be considered as gradient
	vector<vector<vector<int>>> all_masks;

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			vector<vector<int>> temp;
			int *vector1 = all_vectors[i];
			int *vector2 = all_vectors[j];
			for (int row = 0; row < 3; row++)
			{
				vector<int> temp_vector;
				for (int col = 0; col < 3; col++)
				{
					temp_vector.push_back(vector1[row] * vector2[col]);
				}
				temp.push_back(temp_vector);
				temp_vector.clear();
			}
			all_masks.push_back(temp);
		}
	}
	return all_masks;
}