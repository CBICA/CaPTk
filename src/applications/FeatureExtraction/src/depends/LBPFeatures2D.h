#pragma once
#include <iostream>
#include <iterator>
#include <stdlib.h>
#include <string>
#include <algorithm>
#include <exception>
#include <typeinfo>
#include <stdexcept>
#include <stdio.h>
#include <numeric>
#include <functional>
#include <vector>

#include <itkVectorImage.h>
#include "itkImageFileReader.h"
//#include "opencv/cv.h"
#include <itkImageRegionIteratorWithIndex.h>
#include "opencv2/core/core.hpp"
#define PI 3.14

//typedef  float           PixelType;
//static const unsigned int  Dimension = 3;
//typedef itk::Image< PixelType, Dimension > ImageType;
//
using ImageTypeFloat3D = itk::Image< float, 3 >;

//using FeatureBaseType1 = std::vector < std::tuple< std::string, std::string, double > >;


enum MappingType {
	none,
	u2 ,
	ri,
	riu2,
	hf
	
};

class LBPFeatures
{
public:
	LBPFeatures();
	~LBPFeatures();

	struct mapping_table{
		std::vector<int> table;
		int N;

	};
	mapping_table LUT(unsigned int samples, int type);
	std::vector<float> hist(std::vector<float>& src, int bins);

  template < class TImageType = ImageTypeFloat3D>
  std::vector<float> hist_mapping(typename TImageType::Pointer  roi, cv::Mat lbp_blk);

  template <class TImageType = ImageTypeFloat3D>
  void calculateLBP(typename TImageType::Pointer  inputImage, typename TImageType::Pointer  mask,
    typename TImageType::Pointer  roi, int radius, int neighbours, std::string modality, std::string label, /*FeatureBaseType1**/ std::map<std::string, double> &featurevec);
private:
  std::vector<double> lbp_feature_ri;
  std::vector<double> lbp_feature_riu2;
};

template < class TImageType >
std::vector<float> LBPFeatures::hist_mapping(typename TImageType::Pointer  roi, cv::Mat lbp_blk)
{
	std::vector<float> hist_map;
  typename TImageType::SizeType size = roi->GetLargestPossibleRegion().GetSize();
  itk::ImageRegionIteratorWithIndex< TImageType > inputIt(roi, roi->GetLargestPossibleRegion());
	for (unsigned int k = 0; k < size[2]; k++)
	{
    for (unsigned int j = 0; j < size[1]; j++)
		{
      for (unsigned int i = 0; i < size[0]; i++)
			{
			//	std::cout << lbp_blk.at<float>(i, j, k);
				if (roi->GetPixel(inputIt.GetIndex()) != 0)
				{
				//	std::cout << inputIt.GetIndex() << "&&" << roi->GetPixel(inputIt.GetIndex()).GetElement(0);

					hist_map.push_back(lbp_blk.at<float>(i, j, k));
				}
				++inputIt;

			}

		}
	}
	return  hist_map;
}

template < class TImageType >
void LBPFeatures::calculateLBP(typename TImageType::Pointer inputImage, typename TImageType::Pointer mask,
  typename TImageType::Pointer roi, int radius, int neighbours, std::string modality, std::string label, /*FeatureBaseType1**/ std::map<std::string, double> &featurevec)
{
  std::vector<std::string> lbp_featurename;
  std::vector<double> lbp_feature;
  typename TImageType::SizeType size = roi->GetLargestPossibleRegion().GetSize();
  int img_size[] = { (int)size[0], (int)size[1], (int)size[2] };

	cv::Mat imageblk(3, img_size, CV_32F, cv::Scalar(0.0));

	//int counter = 0;
  itk::ImageRegionIteratorWithIndex< TImageType > inputIt(roi, roi->GetLargestPossibleRegion());
	for (unsigned int k = 0; k< size[2]; k++)
	{
    for (unsigned int j = 0; j< size[1]; j++)
		{
      for (unsigned int i = 0; i< size[0]; i++)
			{
				imageblk.at<float>(i, j, k) = roi->GetPixel(inputIt.GetIndex());
				//std::cout << imageblk.at<float>(i, j, k) << ",,";
				++inputIt;
		  	}
		  }
  	}

	cv::Mat lbp_blk(3, img_size, CV_32F, cv::Scalar(0.0));
	cv::Mat lbp_blk_riu2(3, img_size, CV_32F, cv::Scalar(0.0));
	
	LBPFeatures::mapping_table map_table_ri = LUT(neighbours, 2);
	mapping_table map_table_riu2 = LUT(neighbours, 3);

  for (unsigned int x = 0; x < size[0]; x++){
   int slice_sx []= { 1, (int)size[1], (int)size[2] };
		cv::Mat slice = cv::Mat(3, slice_sx, CV_32F, cv::Scalar(0.0));
		imageblk.row(x).copyTo(slice);

		int lbp_sizex[] = { 1, (int)size[1] - 2,(int) size[2] - 2 };
		cv::Mat dst = cv::Mat(3, slice_sx, CV_32F, cv::Scalar(0.0));
		cv::Mat temp = cv::Mat(3, lbp_sizex, CV_32F, cv::Scalar(0.0));
		cv::Mat temp1 = cv::Mat(3, lbp_sizex, CV_32F, cv::Scalar(0.0));

    for (auto n = 0; n<neighbours; n++) {
			// sample points
			float x = static_cast<float>(radius)* cos(2.0*PI*n / static_cast<float>(neighbours));
			float y = static_cast<float>(radius)* -sin(2.0*PI*n / static_cast<float>(neighbours));
			// relative indices
			int fx = static_cast<int>(floor(x));
			int fy = static_cast<int>(floor(y));
			int cx = static_cast<int>(ceil(x));
			int cy = static_cast<int>(ceil(y));
			// fractional part
			float ty = y - fy;
			float tx = x - fx;
			// set interpolation weights
			float w1 = (1 - tx) * (1 - ty);
			float w2 = tx  * (1 - ty);
			float w3 = (1 - tx) *      ty;
			float w4 = tx  *      ty;
			// iterate through your data
      for (auto i = radius; i < slice.size[1] - radius; i++) {
        for (auto j = radius; j < slice.size[2] - radius; j++) {

					float t = w1*slice.at<float>(0, i + fy, j + fx) + w2*slice.at<float>(0, i + fy, j + cx) + w3*slice.at<float>(0, i + cy, j + fx) + w4*slice.at<float>(0, i + cy, j + cx);

					temp.at<float>(0, i - radius, j - radius) += ((t > slice.at<float>(0, i, j)) && (std::abs(t - slice.at<float>(0, i, j)) > std::numeric_limits<float>::epsilon())) << n;
				//	std::cout << temp.at<float>(0, i - radius, j - radius);
				}
			}

		}

		std::vector<int> table = map_table_ri.table;
		for (auto i = 0; i < temp.size[1]; i++){
      for (auto j = 0; j < temp.size[2]; j++){
			//	std::cout << temp.at<float>(0, i, j);
				temp1.at<float>(0, i, j) = table[temp.at<float>(0, i, j)];
			//	std::cout << temp1.at<float>(0, i, j);
			}
		}
		cv::Range ranges[3];
		ranges[0] = cv::Range::all();
		ranges[1] = cv::Range(1, size[1] - 1);
		ranges[2] = cv::Range(1, size[2] - 1);
		temp1.copyTo(dst(ranges));

		dst.copyTo(lbp_blk.row(x));
		table = map_table_riu2.table;
    for (auto i = 0; i < temp.size[1]; i++){
      for (auto j = 0; j < temp.size[2]; j++){
				temp.at<float>(0, i, j) = table[temp.at<float>(0, i, j)];
			}
		}
		temp.copyTo(dst(ranges));
		dst.copyTo(lbp_blk_riu2.row(x));

	}
	std::vector<float> hist_map;
	hist_map = hist_mapping< TImageType >(roi, lbp_blk);
	hist_map = hist(hist_map, map_table_ri.N);
  lbp_feature_ri.insert(lbp_feature_ri.end(),hist_map.begin(),hist_map.end());

	//write to file 
  hist_map = hist_mapping< TImageType >(roi, lbp_blk_riu2);
	hist_map = hist(hist_map, map_table_riu2.N);
  lbp_feature_riu2.insert(lbp_feature_riu2.end(), hist_map.begin(), hist_map.end());

  for (int y = 0; y <(int)size[1]; y++){
		int slice_sy[] = {(int) size[0], 1, (int)size[2] };
		cv::Mat slice = cv::Mat(3, slice_sy, CV_32F, cv::Scalar(0.0));
		imageblk.col(y).copyTo(slice);
		//std::cout << slice.at<float>(0, 0, 0) << ",," << slice.at<float>(0, 0, 1) << slice.at<float>(0, 0, 2) << slice.at<float>(3, 0, 2);
    int lbp_sizey[] = { (int)size[0] - 2, 1, (int)size[2] - 2 };
		cv::Mat dst = cv::Mat(3, slice_sy, CV_32F, cv::Scalar(0.0));
		cv::Mat temp = cv::Mat(3, lbp_sizey, CV_32F, cv::Scalar(0.0));
		cv::Mat temp1 = cv::Mat(3, lbp_sizey, CV_32F, cv::Scalar(0.0));

    for (int n = 0; n < (int)neighbours; n++) {
			// sample points
			float x = static_cast<float>(radius)* cos(2.0*PI*n / static_cast<float>(neighbours));
			float y = static_cast<float>(radius)* -sin(2.0*PI*n / static_cast<float>(neighbours));
			// relative indices
			int fx = static_cast<int>(floor(x));
			int fy = static_cast<int>(floor(y));
			int cx = static_cast<int>(ceil(x));
			int cy = static_cast<int>(ceil(y));
			// fractional part
			float ty = y - fy;
			float tx = x - fx;
			// set interpolation weights
			float w1 = (1 - tx) * (1 - ty);
			float w2 = tx  * (1 - ty);
			float w3 = (1 - tx) *      ty;
			float w4 = tx  *      ty;
			// iterate through your data
      for (auto i = radius; i < slice.size[0] - radius; i++) {
        for (auto j = radius; j < slice.size[2] - radius; j++) {

					float t = w1*slice.at<float>(i + fy, 0, j + fx) + w2*slice.at<float>(i + fy, 0, j + cx) + w3*slice.at<float>(i + cy, 0, j + fx) + w4*slice.at<float>(i + cy, 0, j + cx);

					temp.at<float>(i - radius, 0, j - radius) += ((t > slice.at<float>(i, 0, j)) && (std::abs(t - slice.at<float>(i, 0, j)) > std::numeric_limits<float>::epsilon())) << n;
				//	std::cout << temp.at<float>(i - radius, 0, j - radius);
				}
			}

		}

		std::vector<int> table = map_table_ri.table;
    for (auto i = 0; i < temp.size[0]; i++){
      for (auto j = 0; j < temp.size[2]; j++){
			//	std::cout << temp.at<float>(i, 0, j);
				temp1.at<float>(i, 0, j) = table[temp.at<float>(i, 0, j)];
			//	std::cout << temp1.at<float>(i, 0, j);
			}
		}
		cv::Range ranges[3];
		ranges[1] = cv::Range::all();
		ranges[0] = cv::Range(1, size[0] - 1);
		ranges[2] = cv::Range(1, size[2] - 1);
		temp1.copyTo(dst(ranges));
		dst.copyTo(lbp_blk.col(y));

		table = map_table_riu2.table;
    for (auto i = 0; i < temp.size[0]; i++){
      for (auto j = 0; j < temp.size[2]; j++){
			//	std::cout << temp.at<float>(i, 0, j);
				temp.at<float>(i, 0, j) = table[temp.at<float>(i, 0, j)];
			//	std::cout << temp.at<float>(i, 0, j);
			}
		}
		temp.copyTo(dst(ranges));
		dst.copyTo(lbp_blk_riu2.col(y));
	}

  hist_map = hist_mapping< TImageType >(roi, lbp_blk);
	hist_map = hist(hist_map, map_table_ri.N);
  lbp_feature_ri.insert(lbp_feature_ri.end(), hist_map.begin(), hist_map.end());
	//write to file 
  hist_map = hist_mapping< TImageType >(roi, lbp_blk_riu2);
	hist_map = hist(hist_map, map_table_riu2.N);
  lbp_feature_riu2.insert(lbp_feature_riu2.end(), hist_map.begin(), hist_map.end());


  for (int z = 0; z < (int)size[2]; z++){
		int slice_sz[] = { (int)size[0], (int)size[1], 1 };
		cv::Mat slice = cv::Mat(3, slice_sz, CV_32F, cv::Scalar(0.0));
		cv::Range ranges[3];
		ranges[0] = cv::Range::all();
		ranges[1] = cv::Range::all();
		ranges[2] = cv::Range(z, z + 1);

		slice = imageblk(ranges).clone();
	//	std::cout << slice.at<float>(0, 0, 0) << ",," << slice.at<float>(0, 1, 0) << slice.at<float>(0, 2, 0) << slice.at<float>(3, 2, 0);
    int lbp_sizez[] = { (int)size[0] - 2, (int)size[1] - 2, 1 };
		cv::Mat dst = cv::Mat(3, slice_sz, CV_32F, cv::Scalar(0.0));
		cv::Mat temp = cv::Mat(3, lbp_sizez, CV_32F, cv::Scalar(0.0));
		cv::Mat temp1 = cv::Mat(3, lbp_sizez, CV_32F, cv::Scalar(0.0));


    for (auto n = 0; n<neighbours; n++) {
			// sample points
			float x = static_cast<float>(radius)* cos(2.0*PI*n / static_cast<float>(neighbours));
			float y = static_cast<float>(radius)* -sin(2.0*PI*n / static_cast<float>(neighbours));
			// relative indices
			int fx = static_cast<int>(floor(x));
			int fy = static_cast<int>(floor(y));
			int cx = static_cast<int>(ceil(x));
			int cy = static_cast<int>(ceil(y));
			// fractional part
			float ty = y - fy;
			float tx = x - fx;
			// set interpolation weights
			float w1 = (1 - tx) * (1 - ty);
			float w2 = tx  * (1 - ty);
			float w3 = (1 - tx) *      ty;
			float w4 = tx  *      ty;
			// iterate through your data
      for (auto i = radius; i < slice.size[0] - radius; i++) {
        for (auto j = radius; j < slice.size[1] - radius; j++) {

					float t = w1*slice.at<float>(i + fy, j + fx, 0) + w2*slice.at<float>(i + fy, j + cx, 0) + w3*slice.at<float>(i + cy, j + fx, 0) + w4*slice.at<float>(i + cy, j + cx, 0);

					temp.at<float>(i - radius, j - radius, 0) += ((t > slice.at<float>(i, j, 0)) && (std::abs(t - slice.at<float>(i, j, 0)) > std::numeric_limits<float>::epsilon())) << n;
				//	std::cout << temp.at<float>(i - radius, j - radius, 0);
				}
			}

		}
		std::vector<int> table = map_table_ri.table;
    for (auto i = 0; i < temp.size[0]; i++){
      for (auto j = 0; j < temp.size[1]; j++){
				//std::cout << temp.at<float>(i, j, 0);
				temp1.at<float>(i, j, 0) = table[temp.at<float>(i, j, 0)];
			//	std::cout << temp1.at<float>(i, j, 0);
			}
		}
		cv::Range range[3];
		range[2] = cv::Range::all();
		range[0] = cv::Range(1, size[0] - 1);
		range[1] = cv::Range(1, size[1] - 1);
		temp1.copyTo(dst(range));
		dst.copyTo(lbp_blk(ranges));

		table = map_table_riu2.table;
    for (auto i = 0; i < temp.size[0]; i++){
      for (auto j = 0; j < temp.size[1]; j++){
			//	std::cout << temp.at<float>(i, j, 0);
				temp.at<float>(i, j, 0) = table[temp.at<float>(i, j, 0)];
			//	std::cout << temp.at<float>(i, j, 0);
			}
		}
		temp.copyTo(dst(range));
		dst.copyTo(lbp_blk_riu2(ranges));


	}

  hist_map = hist_mapping< TImageType >(roi, lbp_blk);
	hist_map = hist(hist_map, map_table_ri.N);
  lbp_feature_ri.insert(lbp_feature_ri.end(), hist_map.begin(), hist_map.end());
	//write to file 
  hist_map = hist_mapping< TImageType >(roi, lbp_blk_riu2);
	hist_map = hist(hist_map, map_table_riu2.N);
  lbp_feature_riu2.insert(lbp_feature_riu2.end(), hist_map.begin(), hist_map.end());
  std::stringstream ss;


  for (int i = 0; i <(int)( lbp_feature_ri.size()); i++)
  {
    lbp_featurename.push_back("ri" + std::to_string(i));
  }
  for (int i = 0; i <(int)lbp_feature_riu2.size(); i++)
  {
    lbp_featurename.push_back("riu2_" + std::to_string(i));
  }


  lbp_feature.insert(lbp_feature.end(),lbp_feature_ri.begin(),lbp_feature_ri.end());
  lbp_feature.insert(lbp_feature.end(), lbp_feature_riu2.begin(), lbp_feature_riu2.end());

  for (size_t i = 0; i < lbp_feature.size(); i++)
  {
    //std::string p = "_LBP2D_" + modality + "_" + label + "_" + lbp_featurename.at(i);
    //featurevec->push_back(std::tuple<std::string, std::string, double>("Texture", p, lbp_feature.at(i)));
    featurevec["LBP2D_" + std::to_string(i) + "_" + lbp_featurename[i]] = lbp_feature[i];
  }


}
