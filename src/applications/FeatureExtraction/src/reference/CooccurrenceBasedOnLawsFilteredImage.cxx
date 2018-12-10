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
#include "itkCastImageFilter.h"
#include <itkDenseFrequencyContainer2.h>
#include "itkHistogramToTextureFeaturesFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkScalarImageToCooccurrenceMatrixFilter.h"
#include "ImageIOFunctions.h"
#include <itkSqrtImageFilter.h>



#define isnan(x) ((x) != (x))

using namespace std;

//predefine functions

void usage(int argc, char *argv[]);
void parse(int argc, char *argv[], char* &fn_input, char* &fn_output, char* &fn_mask4Process, unsigned int &winSZ, unsigned int &numberOfBins, unsigned int &slidingDist, std::vector<int> *offsetvectors, int &numVectors);


//template <unsigned int ImageDimension>
int CooccurrenceBasedOnLawsFilteredImage(int argc, char *argv[]);
itk::Image<double, 2>::Pointer mask_matrix_select(int a, int b);

//maind code

int main(int argc, char *argv[])
{

	usage(argc, argv);
	
	CooccurrenceBasedOnLawsFilteredImage(argc, argv);
	
}


//The idea of this feature is:
//1. Apply laws filter on the image. 2. Extract co-occurrence feature on the filtered images.
//Since laws feature is only defined on 2D image, now we only have 2D version for this feature.
const unsigned int ImageDimension = 2;
//template <unsigned int ImageDimension>
int CooccurrenceBasedOnLawsFilteredImage(int argc, char *argv[])
{

	//parse inputs
	char* inputImageFn = NULL; char* fn_output = NULL;
	char* fn_mask4Process = NULL;
	unsigned int winSZ; unsigned int slidingDist;
	unsigned int numberOfBins;
	std::vector<int> offset_vectors[10]; int numVectors;

	parse(argc, argv, inputImageFn, fn_output, fn_mask4Process, winSZ, numberOfBins, slidingDist, offset_vectors, numVectors);

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
	const ImageType::DirectionType& dirCosine = inputImage->GetDirection();
	const ImageType::SpacingType& inputImageSpacing = inputImage->GetSpacing();

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

	// Yifan's modification
	vector<vector <double>> all_features;
	vector<double> temp_feature;
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

		

		//Yifan's modification
		temp_feature.clear();
		for (int i = 0; i < 5; i++)
		{
			for (int j = 0; j <= i; j++)
			{
				//
				if (i == j)
				{
					//define masks 
					ImageType::Pointer masks = mask_matrix_select(i, j);
					////////////////////////////////////////////////////////////////////////////////
					//compute features here Note: ItI and ItM proved above
					////////////////////////////////////////////////////////////////////////////////
					typedef itk::ConvolutionImageFilter<ImageType> FilterType;
					FilterType::Pointer convolutionFilter = FilterType::New();
					convolutionFilter->SetInput(inputImageRoiIm);
					//convolve the image with mask
					convolutionFilter->SetKernelImage(masks);
					convolutionFilter->Update();


					itk::ImageRegionIteratorWithIndex<ImageType> ItI(convolutionFilter->GetOutput(),
						convolutionFilter->GetOutput()->GetLargestPossibleRegion());

					long maxnum = (long)(winSZ*winSZ) / 3;
					numberOfBins = (numberOfBins>maxnum) ? maxnum : numberOfBins;
					//offfsets
					typedef itk::Statistics::ScalarImageToCooccurrenceMatrixFilter<ImageType>  CooccurrenceMatrixGeneratorType;
					CooccurrenceMatrixGeneratorType::Pointer generator = CooccurrenceMatrixGeneratorType::New();
					CooccurrenceMatrixGeneratorType::OffsetVectorPointer offsets = CooccurrenceMatrixGeneratorType::OffsetVector::New();
					CooccurrenceMatrixGeneratorType::OffsetType offset;

					std::vector<int> vector;
					for (int i = 0; i<numVectors; i++)
					{
						vector = offset_vectors[i];
						if (vector.size() > 0 && vector.size() != ImageDimension)
						{
							std::cerr << "Error:  offset size does not equal image dimension." << std::endl;
							return EXIT_FAILURE;
						}
						else
							for (unsigned int d = 0; d < ImageDimension; d++)
								offset[d] = vector[d];

						offsets->push_back(offset);
					}
					generator->SetOffsets(offsets);

					//local mask window
					ImageType::Pointer maskRoiIm = ImageType::New();
					ImageType::PointType origin;
					origin.Fill(0.0);
					maskRoiIm->SetOrigin(origin);
					maskRoiIm->SetSpacing(inputImage->GetSpacing());

					ImageType::SizeType size;
					size.Fill(winSZ);
					ImageType::RegionType region;
					region.SetSize(size);
					maskRoiIm->SetRegions(region);
					maskRoiIm->Allocate();
					maskRoiIm->FillBuffer(itk::NumericTraits<PixelType>::One);

					itk::ImageRegionIterator<ImageType> ItM(maskRoiIm, maskRoiIm->GetLargestPossibleRegion());

					PixelType maxValue = itk::NumericTraits<PixelType>::min();
					PixelType minValue = itk::NumericTraits<PixelType>::max();

					for (ItM.GoToBegin(), ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItM, ++ItI)
					{
						if (ItM.Get() == itk::NumericTraits<PixelType>::One)
						{
							if (ItI.Get() < minValue)
							{
								minValue = ItI.Get();
							}
							else if (ItI.Get() > maxValue)
							{
								maxValue = ItI.Get();
							}
						}

						if (isnan(ItI.Get()) || std::isinf(ItI.Get()))
						{
							ItM.Set(itk::NumericTraits<PixelType>::Zero);
						}
					}
					//Input here is the convolved images
					generator->SetInput(convolutionFilter->GetOutput());
					generator->SetMaskImage(maskRoiIm);
					generator->SetNumberOfBinsPerAxis(numberOfBins);
					generator->SetPixelValueMinMax(minValue, maxValue);
					generator->Update();

					typedef itk::Statistics::HistogramToTextureFeaturesFilter<CooccurrenceMatrixGeneratorType::HistogramType> CalculatorType;
					CalculatorType::Pointer calculator = CalculatorType::New();
					calculator->SetInput(generator->GetOutput());
					calculator->Update();

					//typedef double RealType;
					//8 features are calculated
					temp_feature.push_back(calculator->GetEnergy());
					temp_feature.push_back(calculator->GetEntropy());
					temp_feature.push_back(calculator->GetCorrelation());
					temp_feature.push_back(calculator->GetInverseDifferenceMoment());
					temp_feature.push_back(calculator->GetInertia());
					temp_feature.push_back(calculator->GetClusterShade());
					temp_feature.push_back(calculator->GetClusterProminence());
					temp_feature.push_back(calculator->GetHaralickCorrelation());
				}
				//when i!=j, we compute the image by a combination of two transposed masks;
				//the pixel value of final image is sqrt(I1^2+I2^2)
				//this idea is from the definition of gradient on the image(Sobel operator: https://en.wikipedia.org/wiki/Sobel_operator); if I1 is gradient on x-axis (one of laws mask can be the filter to compute gradient), the I2 will be the gradient on y-axis. 
				//The sobel operator could be two of laws feature filters when mask size=3 
				else
				{
					ImageType::Pointer masks1 = mask_matrix_select(i, j);
					ImageType::Pointer masks2 = mask_matrix_select(j, i);
					////////////////////////////////////////////////////////////////////////////////
					//compute features here Note: ItI and ItM proved above
					////////////////////////////////////////////////////////////////////////////////
					typedef itk::ConvolutionImageFilter<ImageType> FilterType;
					FilterType::Pointer convolutionFilter1 = FilterType::New();
					convolutionFilter1->SetInput(inputImageRoiIm);
					convolutionFilter1->SetKernelImage(masks1);
					convolutionFilter1->Update();

					typedef itk::ConvolutionImageFilter<ImageType> FilterType;
					FilterType::Pointer convolutionFilter2 = FilterType::New();
					convolutionFilter2->SetInput(inputImageRoiIm);
					convolutionFilter2->SetKernelImage(masks2);
					convolutionFilter2->Update();

					typedef itk::PowImageFilter< ImageType, ImageType, ImageType> PowImageFilterType;
					PowImageFilterType::Pointer powImageFilter1 = PowImageFilterType::New();
					PowImageFilterType::Pointer powImageFilter2 = PowImageFilterType::New();
					//skewness
					powImageFilter1->SetInput(convolutionFilter1->GetOutput());
					powImageFilter1->SetConstant2(2);
					powImageFilter1->Update();

					powImageFilter2->SetInput(convolutionFilter2->GetOutput());
					powImageFilter2->SetConstant2(2);
					powImageFilter2->Update();

					typedef itk::AddImageFilter <ImageType, ImageType, ImageType> AddImageFilterType;
					AddImageFilterType::Pointer addImageFilter = AddImageFilterType::New();
					addImageFilter->SetInput1(powImageFilter1->GetOutput());
					addImageFilter->SetInput2(powImageFilter2->GetOutput());
					addImageFilter->Update();

					typedef itk::SqrtImageFilter<ImageType, ImageType> SqrtFilterType;
					SqrtFilterType::Pointer sqrtFilter = SqrtFilterType::New();
					sqrtFilter->SetInput(addImageFilter->GetOutput());
					sqrtFilter->Update();

					itk::ImageRegionIteratorWithIndex<ImageType> ItI(sqrtFilter->GetOutput(),
						sqrtFilter->GetOutput()->GetLargestPossibleRegion());

					long maxnum = (long)(winSZ*winSZ) / 3;
					numberOfBins = (numberOfBins>maxnum) ? maxnum : numberOfBins;
					//offsets
					typedef itk::Statistics::ScalarImageToCooccurrenceMatrixFilter<ImageType>  CooccurrenceMatrixGeneratorType;
					CooccurrenceMatrixGeneratorType::Pointer generator = CooccurrenceMatrixGeneratorType::New();
					CooccurrenceMatrixGeneratorType::OffsetVectorPointer offsets = CooccurrenceMatrixGeneratorType::OffsetVector::New();
					CooccurrenceMatrixGeneratorType::OffsetType offset;

					std::vector<int> vector;
					for (int i = 0; i<numVectors; i++)
					{
						vector = offset_vectors[i];
						if (vector.size() > 0 && vector.size() != ImageDimension)
						{
							std::cerr << "Error:  offset size does not equal image dimension." << std::endl;
							return EXIT_FAILURE;
						}
						else
							for (unsigned int d = 0; d < ImageDimension; d++)
								offset[d] = vector[d];

						offsets->push_back(offset);
					}
					generator->SetOffsets(offsets);

					//local mask window
					ImageType::Pointer maskRoiIm = ImageType::New();
					ImageType::PointType origin;
					origin.Fill(0.0);
					maskRoiIm->SetOrigin(origin);
					maskRoiIm->SetSpacing(sqrtFilter->GetOutput()->GetSpacing());

					ImageType::SizeType size;
					size.Fill(winSZ);
					ImageType::RegionType region;
					region.SetSize(size);
					maskRoiIm->SetRegions(region);
					maskRoiIm->Allocate();
					maskRoiIm->FillBuffer(itk::NumericTraits<PixelType>::One);

					itk::ImageRegionIterator<ImageType> ItM(maskRoiIm, maskRoiIm->GetLargestPossibleRegion());

					PixelType maxValue = itk::NumericTraits<PixelType>::min();
					PixelType minValue = itk::NumericTraits<PixelType>::max();

					for (ItM.GoToBegin(), ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItM, ++ItI)
					{
						if (ItM.Get() == itk::NumericTraits<PixelType>::One)
						{
							if (ItI.Get() < minValue)
							{
								minValue = ItI.Get();
							}
							else if (ItI.Get() > maxValue)
							{
								maxValue = ItI.Get();
							}
						}

						if (isnan(ItI.Get()) || std::isinf(ItI.Get()))
						{
							ItM.Set(itk::NumericTraits<PixelType>::Zero);
						}
					}
					
					

					generator->SetInput(sqrtFilter->GetOutput());
					generator->SetMaskImage(maskRoiIm);
					generator->SetNumberOfBinsPerAxis(numberOfBins);
					generator->SetPixelValueMinMax(minValue, maxValue);
					generator->Update();

					typedef itk::Statistics::HistogramToTextureFeaturesFilter<CooccurrenceMatrixGeneratorType::HistogramType> CalculatorType;
					CalculatorType::Pointer calculator = CalculatorType::New();
					calculator->SetInput(generator->GetOutput());
					calculator->Update();
					
					temp_feature.push_back(calculator->GetEnergy());
					temp_feature.push_back(calculator->GetEntropy());
					temp_feature.push_back(calculator->GetCorrelation());
					temp_feature.push_back(calculator->GetInverseDifferenceMoment());
					temp_feature.push_back(calculator->GetInertia());
					temp_feature.push_back(calculator->GetClusterShade());
					temp_feature.push_back(calculator->GetClusterProminence());
					temp_feature.push_back(calculator->GetHaralickCorrelation());
				}
			}

		}
		
		//convolutionFilter->GetOutput()
		//15 masks*8 features=120 features
		//feature for each mask
		all_features.push_back(temp_feature);
	}
	
	for (int i = 0; i < 120; i++)
	{
		char* measure[8] = { "_energy.nii", "_entropy.nii", "_correlation.nii", "_inversedifferencemoment.nii", "_inertia.nii", "_clustershade.nii", "_clusterprominence.nii", "_haralickcorrelation.nii"};
		char fullname[300];
		strcpy(fullname, fn_output);
		char buffer[2];
		int integer_part = i / 8 + 1;
		int residual = i % 8;
		strcat(fullname, "_mask_");
		string s1 = std::to_string(integer_part);
		char const *pchar1 = s1.c_str();
		strcat(fullname, pchar1);
		strcat(fullname, measure[residual]);

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
		result_image->FillBuffer(itk::NumericTraits< double >::Zero);

		for (int k = 0; k < result_size[0]; k++)
		{
			for (int l = 0; l < result_size[1]; l++)
			{
				ImageType::IndexType temp_index;
				temp_index[0] = k;
				temp_index[1] = l;
				result_image->SetPixel(temp_index, all_features[l*result_size[0] + k][i]);
			}
		}

		typedef itk::ImageFileWriter< ImageType > WriterType;
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

void parse(int argc, char *argv[], char* &fn_input, char* &fn_output, char* &fn_mask4Process, unsigned int &winSZ, unsigned int &numberOfBins, unsigned int &slidingDist, std::vector<int> *offsetvectors, int &numVectors)
{
	//set default values
	unsigned int ImageDimension = atoi(argv[1]);
	fn_input = NULL;
	fn_output = NULL;
	fn_mask4Process = NULL;
	winSZ = 63; slidingDist = 63;
	numberOfBins = 50;
	numVectors = 1;
	for (int d = 0; d<ImageDimension; d++)
		offsetvectors[0].push_back(11);
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
			numberOfBins = atoi(&myargv[1][0]);
			break;
		}
		--myargc; --myargc;
		++myargv; ++myargv;
	}
}

itk::Image<double, 2>::Pointer mask_matrix_select(int a, int b)
{
	int all_vectors[5][5] = { { 1, 4, 6, 4, 1 },
	{ -1, -2, 0, 2, 1 },
	{ -1, 0, 2, 0, -1 },
	{ -1, 2, 0, -2, 1 },
	{ 1, -4, 6, -4, 1 } };

	vector<vector<int>> masks;
	vector<int> temp;
	for (int i = 0; i < 5; i++)
	{
		temp.clear();
		for (int j = 0; j < 5; j++)
		{
			temp.push_back(all_vectors[a][i]*all_vectors[b][j]);
		}
		masks.push_back(temp);
	}
	
	typedef itk::Image<double, 2> ImageType;
	ImageType::Pointer image_temp_mask = ImageType::New();
	ImageType::RegionType region;
	ImageType::IndexType start;
	start.Fill(0);
	ImageType::SizeType size;
	size[0] = 5;
	size[1] = 5;
	region.SetSize(size);
	image_temp_mask->SetRegions(region);
	image_temp_mask->Allocate();
	image_temp_mask->FillBuffer(itk::NumericTraits< double >::Zero);
	for (int j = 0; j< 5; j++)
	{
		for (int k = 0; k<5; k++)
		{

			ImageType::IndexType pixelIndex;
			pixelIndex[0] = j;
			pixelIndex[1] = k;
			image_temp_mask->SetPixel(pixelIndex, masks[j][k]);
			ImageType::IndexType index;
			index[0] = j;
			index[1] = k;
		}
	}
	return image_temp_mask;
}