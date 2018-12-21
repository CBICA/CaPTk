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
#include "math.h"
#include "itkImage.h"
#include "itkConvolutionImageFilter.h"
#include "itkImageRegionIterator.h"



#define _USE_MATH_DEFINES

using namespace std;

#define isnan(x) ((x) != (x))

//predefine functions

itk::Image<float, 2>::Pointer convoltion_image_matrix(itk::Image<float, 2>::Pointer ImageI, vector<vector <float>> matrixA)
{
	typedef itk::Image<float, 2 > TImageType;
	TImageType::Pointer image = ImageI;
	TImageType::RegionType region = image->GetLargestPossibleRegion();
	TImageType::SizeType size = region.GetSize();

	TImageType::Pointer result_matrix = TImageType::New();
	TImageType::RegionType region_new;
	TImageType::IndexType start;
	start.Fill(0);
	TImageType::SizeType size_new_matrix;
	size_new_matrix[0] = size[0] - matrixA.size() + 1;
	size_new_matrix[1] = size[1] - matrixA[0].size() + 1;
	cout << size_new_matrix[0] << "," << size_new_matrix[1] << endl;
	region_new.SetSize(size_new_matrix);
	result_matrix->SetRegions(region_new);
	result_matrix->Allocate();
	result_matrix->FillBuffer(itk::NumericTraits< float >::Zero);

	for (int i = 0; i < size[0] - matrixA.size() + 1; i++)
	{
		for (int j = 0; j < size[1] - matrixA[0].size() + 1; j++)
		{
			cout << "(" << i << "," << j << ")";
			TImageType::SizeType cropSize1, cropSize2;
			cropSize1[0] = i;
			cropSize1[1] = j;
			cropSize2[0] = size[0] - (i + matrixA.size());
			cropSize2[1] = size[1] - (j + matrixA[0].size());
			typedef itk::CropImageFilter <TImageType, TImageType>
				CropImageFilterType;
			CropImageFilterType::Pointer cropFilter
				= CropImageFilterType::New();
			cropFilter->SetInput(image);
			//cropFilter->SetBoundaryCropSize(cropSize);
			cropFilter->SetUpperBoundaryCropSize(cropSize2);
			cropFilter->SetLowerBoundaryCropSize(cropSize1);
			cropFilter->Update();
			TImageType::Pointer croppedimage = cropFilter->GetOutput();
			TImageType::SizeType size2 = croppedimage->GetLargestPossibleRegion().GetSize();


			float sum = 0;
			typedef itk::ImageRegionIteratorWithIndex<TImageType> IteratorType;
			IteratorType it(croppedimage, croppedimage->GetLargestPossibleRegion());
			//cout << endl << it.Get() << "," << it.GetIndex();
			int index1 = it.GetIndex()[0];
			int index2 = it.GetIndex()[1];
			//cout << "," << matrixA[i][j];

			while (!it.IsAtEnd())
			{
				it.GetIndex();
				//std::cout << matrixA[it.GetIndex()[0] - index1][it.GetIndex()[1] - index2] << ',' << it.Get() << endl;
				sum = sum + it.Get()*matrixA[it.GetIndex()[1] - index2][it.GetIndex()[0] - index1];
				++it;

			}

			//cout << sum << endl;
			TImageType::IndexType pixelIndex;
			pixelIndex[0] = j;
			pixelIndex[1] = i;
			result_matrix->SetPixel(pixelIndex, sum);
		}
	}

	return result_matrix;
}
