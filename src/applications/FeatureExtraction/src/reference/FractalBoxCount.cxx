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


#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkConstantPadImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkImageDuplicator.h"
#include <itkMirrorPadImageFilter.h>
#include "itkRegionOfInterestImageFilter.h"
#include <itkExtractImageFilter.h>

//#include "commonheaders\cbicaCmdParser.h"
//#include "commonheaders\cbicaITKSafeImageIO.h"
//#include "commonheaders\cbicaUtilities.h"

#include <string>
#include <vector>
using namespace std;
#define isnan(x) ((x) != (x))


using DefaultImageType = itk::Image < float, 2 >;

void usage(int argc, char *argv[]);
void parse(int argc, char *argv[], char* &fn_input, char* &fn_output, char* &fn_mask4Process, unsigned int &winSZ, unsigned int &slidingDist);

//This feature can be only defined on 2D space.
int CalculateBoxCount(int argc, char *argv[]);
vector<double> dirichlet_energy_gradient(std::vector<int> u);
double log_with_zero(double x);
//maind code

//when fitting the regression line, we have log function on the expression, the log is defined as this in order to avoid inf value.
double log_with_zero(double x)
{
	if (x == 0)
		return -8;
	else
		return log(x);
}
int main(int argc, char *argv[])
{

	usage(argc, argv);
	CalculateBoxCount(argc, argv);
}

const unsigned int ImageDimension = 2;
int CalculateBoxCount(int argc, char *argv[])
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
	vector <double> box_counting_dimension;
	vector <double> minkovski_dimension;
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
        // this imageROI represents a lattice window
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
		ImageType::SizeType size = inputImageRoiIm->GetLargestPossibleRegion().GetSize();
		int width = 0;
        // supposed to be window size
		if (size[0] > size[1])
			width = size[0];
		else
			width = size[1];

		float p = std::log(width) / std::log(2);
		p = std::ceil(p);
		width = std::pow(2, p);

		typedef itk::StatisticsImageFilter<ImageType> StatisticsImageFilterType;

		typedef itk::MirrorPadImageFilter <ImageType, ImageType>
			MirrorPadImageFilterType;

		ImageType::SizeType lowerExtendRegion;
		lowerExtendRegion[0] = 0;
		lowerExtendRegion[1] = 0;

		ImageType::SizeType upperExtendRegion;
		upperExtendRegion[0] = width - size[0];
		upperExtendRegion[1] = width - size[1];

        // why mirror padding
		MirrorPadImageFilterType::Pointer padFilter
			= MirrorPadImageFilterType::New();
		
		padFilter->SetInput(inputImageRoiIm);
		padFilter->SetPadLowerBound(lowerExtendRegion);
		padFilter->SetPadUpperBound(upperExtendRegion);
		padFilter->Update();
		ImageType::Pointer paddedimage = padFilter->GetOutput();
		ImageType::SizeType newsize = paddedimage->GetLargestPossibleRegion().GetSize();
		StatisticsImageFilterType::Pointer statisticsImageFilter
			= StatisticsImageFilterType::New();
		statisticsImageFilter->SetInput(paddedimage);
		statisticsImageFilter->Update();


		vector<vector <double>> image_matrix;
		vector<double> temp2;
		for (int j = 0; j < newsize[1]; j++)
		{
			temp2.clear();
			for (int i = 0; i < newsize[0]; i++)
			{
				ImageType::IndexType index;
				index[0] = i;
				index[1] = j;
				temp2.push_back(paddedimage->GetPixel(index));
			}
			image_matrix.push_back(temp2);
		}
		//Definition of box-counting and Minkowski fractal dimension is in this paper "Fractal Analysis of Mammographic Parenchymal Patterns in Breast Cancer Risk Assessment", equation 2 for box-counting; equation 3 and 4 for minkowski;
		vector<double> area_sum;
		vector<double> minkovski_diff;
		vector<double> border_length;
		vector<double> border_length_minkov;
		for (int g = 0; g < p; g++)
		{
			//sum is for box-counting; sum1 is for minkowski
            // Michael: what is this size1?
			int size1 = std::pow(2, g);
			float sum = 0;
			float sum1 = 0;
			float square_sum1 = 0; float square_sum2 = 0;
			for (int i = 0; i < width - size1; i = i + size1)
			{
				for (int j = 0; j < width - size1; j = j + size1)
				{
					//box counting
					for (int k = 0; k < size1; k++)
					{
						for (int l = 0; l < size1; l++)
						{
							//square_sum1 and square_sum2 are part of box-counting equation 2 in the reference
							square_sum1 += abs(image_matrix[i + k][j + l] - image_matrix[i + k + size1][j + l]);
							square_sum2 += abs(image_matrix[i + k][j + l] - image_matrix[i + k][j + l + size1]);
						}	
					}
				}
			}
			//Minkovski
			if (size1 > 1)
			for (int i = 0; i < width; i++)
			{
				for (int j = 0; j < width; j++)
				{
					
						float temp_min = image_matrix[i][j]; float temp_max = image_matrix[i][j];
						for (int k = max(0,i-size1/2); k < min(i+size1/2,width); k++)
						{
							for (int l = max(0, j - size1 / 2); l < min(j + size1 / 2, width); l++)
							{
								//temp_min and temp_max is from equation 4 in the reference
								temp_min = temp_min< image_matrix[k][l] ? temp_min : image_matrix[k][l];
								temp_max = temp_max> image_matrix[k][l] ? temp_max : image_matrix[k][l];
							}
						}
						//equation 4
						sum1 += temp_max - temp_min;
				}
			}
			//equation 2
			sum = sum + size1*size1 + (abs(square_sum1) + abs(square_sum2)) * size1;
			//box counting
			area_sum.push_back(sum);
			border_length.push_back(size1+0.0);
			//Minkovski
			if (size1 > 1)
			{
				//equation 3
				minkovski_diff.push_back(sum1/(size1*size1*size1));
				border_length_minkov.push_back(1.0/size1);
			}
			//cout << sum << endl;
		}
		//box_counting
		double sum_y = 0;
		double sum_x = 0;
		double sum_y_x = 0;
		double sum_x_x = 0;
		//fit the line and get the slope which is the fractal dimension
		for (int i = 0; i < border_length.size(); i++)
		{
			sum_y_x += log_with_zero(area_sum[i]) * log_with_zero(border_length[i]);
			sum_y += log_with_zero(area_sum[i]);
			sum_x_x += log_with_zero(border_length[i])*log_with_zero(border_length[i]);
			sum_x += log_with_zero(border_length[i]);
		}

		double beta0;
		beta0 = (sum_y_x - sum_y*sum_x / border_length.size()) / (sum_x_x - sum_x*sum_x / border_length.size());
		box_counting_dimension.push_back(2 - beta0);
		//cout << 2-beta0 << endl;
		//second time for minkowski, variable with same meaning
		sum_y = 0;
		sum_x = 0;
		sum_y_x = 0;
		sum_x_x = 0;
		
		for (int i = 0; i < border_length_minkov.size(); i++)
		{
			sum_y_x += log_with_zero(minkovski_diff[i]) * log_with_zero(border_length_minkov[i]);
			sum_y += log_with_zero(minkovski_diff[i]);
			sum_x_x += log_with_zero(border_length_minkov[i])*log_with_zero(border_length_minkov[i]);
			sum_x += log_with_zero(border_length_minkov[i]);
		}

		double beta1;
		beta1 = (sum_y_x - sum_y*sum_x / border_length_minkov.size()) / (sum_x_x - sum_x*sum_x / border_length_minkov.size());
		minkovski_dimension.push_back(beta1);

	}

	string filename1 = "_box_counting.nii";
	string path1 = fn_output + filename1;
	char *cstrname1 = new char[path1.length() + 1];
	strcpy(cstrname1, path1.c_str());

	string filename2 = "_minkovski.nii";
	string path2 = fn_output + filename2;
	char *cstrname2 = new char[path2.length() + 1];
	strcpy(cstrname2, path2.c_str());

	ImageType::Pointer result_image1 = ImageType::New();
	ImageType::Pointer result_image2 = ImageType::New();
	ImageType::RegionType result_region;
	ImageType::SizeType result_size;
	ImageType::IndexType result_index;
	result_size[0] = floor((inputImageSize[0] - winSZ) / slidingDist) + 1;
	result_size[1] = floor((inputImageSize[1] - winSZ) / slidingDist) + 1;
	result_region.SetSize(result_size);
	result_index.Fill(0);
	result_region.SetIndex(result_index);
	result_image1->SetRegions(result_region);
  result_image1->SetSpacing(outputImageSpacing);
  result_image1->SetDirection(inputImage->GetDirection());
	result_image1->Allocate();
	result_image1->FillBuffer(itk::NumericTraits< float >::Zero);

	result_image2->SetRegions(result_region);
  result_image2->SetSpacing(outputImageSpacing);
  result_image2->SetDirection(inputImage->GetDirection());
	result_image2->Allocate();
	result_image2->FillBuffer(itk::NumericTraits< float >::Zero);

	for (int k = 0; k < result_size[0]; k++)
	{
		for (int l = 0; l < result_size[1]; l++)
		{
			ImageType::IndexType temp_index;
			temp_index[0] = k;
			temp_index[1] = l;
			result_image1->SetPixel(temp_index, box_counting_dimension[l*result_size[0] + k]);
			result_image2->SetPixel(temp_index, minkovski_dimension[l*result_size[0] + k]);
		}
	}

	typedef itk::ImageFileWriter< ImageType  > WriterType;
	WriterType::Pointer writer1 = WriterType::New();
	writer1->SetFileName(cstrname1);
	writer1->SetInput(result_image1);

	WriterType::Pointer writer2 = WriterType::New();
	writer2->SetFileName(cstrname2);
	writer2->SetInput(result_image2);

	try
	{
		writer1->Update();
		writer2->Update();
	}
	catch (itk::ExceptionObject & error)
	{
		std::cerr << "Error: " << error << std::endl;
		return EXIT_FAILURE;
	}
	delete[] cstrname1;
	delete[] cstrname2;

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
		}
		index.push_back(temp);
	}
	return index;
}











