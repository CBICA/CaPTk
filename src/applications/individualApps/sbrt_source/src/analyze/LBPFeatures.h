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

#include "itkVectorImage.h"
#include "itkImageFileReader.h"
//#include "opencv/cv.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "opencv2/core/core.hpp"
#define PI 3.14
#include"generateTextureFeatures.h"

//typedef  float           PixelType;
//static const unsigned int  Dimension = 3;
//typedef itk::Image< PixelType, Dimension > ImageType;
//
typedef itk::Image< itk::Vector<PixelType, 2>, Dimension > MaskImageType;
typedef itk::ImageRegionIteratorWithIndex< MaskImageType > IteratorType_withmask;



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

	template<typename TPixel, unsigned int VDimension>
	std::vector<float> hist_mapping(itk::Image<TPixel, VDimension> *roi, cv::Mat lbp_blk);

	template<typename TPixel=float, unsigned int VDimension>
	void calculateLBP(itk::Image<itk::Vector<TPixel, 2>, VDimension> *inputImage, itk::Image<TPixel, VDimension> *mask,
    itk::Image<itk::Vector<TPixel, 2>, VDimension> *roi, int radius, int neighbours, std::vector<std::string>* lbp_featurename,
    std::vector<double>* lbp_feature);
private:
  std::vector<double> lbp_feature_ri;
  std::vector<double> lbp_feature_riu2;
};


template<typename TPixel, unsigned int VDimension>
std::vector<float> LBPFeatures::hist_mapping(itk::Image<TPixel, VDimension> *roi, cv::Mat lbp_blk)
{
	std::vector<float> hist_map;
	itk::Size<3> dim = roi->GetLargestPossibleRegion().GetSize();
  IteratorType_withmask inputIt(roi, roi->GetLargestPossibleRegion());
	for (unsigned int k = 0; k < dim[2]; k++)
	{
    for (unsigned int j = 0; j < dim[1]; j++)
		{
      for (unsigned int i = 0; i < dim[0]; i++)
			{
			//	std::cout << lbp_blk.at<float>(i, j, k);
				if (roi->GetPixel(inputIt.GetIndex())[1] != 0)
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

template<typename TPixel, unsigned int VDimension>
void LBPFeatures::calculateLBP(itk::Image<itk::Vector<TPixel, 2>, VDimension> *inputImage, itk::Image<TPixel, VDimension> *mask,
  itk::Image<itk::Vector<TPixel, 2>, VDimension> *roi, int radius, int neighbours, std::vector<std::string>* lbp_featurename,
  std::vector<double>* lbp_feature)
{
	itk::Size<Dimension> dim = roi->GetLargestPossibleRegion().GetSize();
  int img_size[] = { (int)dim[0], (int)dim[1], (int)dim[2] };

	cv::Mat imageblk(3, img_size, CV_32F, cv::Scalar(0.0));

	//int counter = 0;
	IteratorType_withmask inputIt(roi, roi->GetLargestPossibleRegion());
	for (unsigned int k = 0; k< dim[2]; k++)
	{
    for (unsigned int j = 0; j< dim[1]; j++)
		{
      for (unsigned int i = 0; i< dim[0]; i++)
			{
				imageblk.at<float>(i, j, k) = roi->GetPixel(inputIt.GetIndex())[0];
				//std::cout << imageblk.at<float>(i, j, k) << ",,";
				++inputIt;
		  	}
		  }
  	}

	cv::Mat lbp_blk(3, img_size, CV_32F, cv::Scalar(0.0));
	cv::Mat lbp_blk_riu2(3, img_size, CV_32F, cv::Scalar(0.0));
	// std::vector<cv::Mat>lbp_blk_riu2(3, cv::Mat(3, img_size, CV_32F, cv::Scalar(0.0)));
	
	LBPFeatures::mapping_table map_table_ri = LUT(neighbours, 2);
	mapping_table map_table_riu2 = LUT(neighbours, 3);

  for (unsigned int x = 0; x < dim[0]; x++){
   int slice_sx []= { 1, (int)dim[1], (int)dim[2] };
		cv::Mat slice = cv::Mat(3, slice_sx, CV_32F, cv::Scalar(0.0));
		imageblk.row(x).copyTo(slice);
	//	std::cout << slice.at<float>(0, 0, 0) << ",," << slice.at<float>(0, 1, 0) << slice.at<float>(0, 0, 2) << slice.at<float>(0, 3, 2);
		int lbp_sizex[] = { 1, (int)dim[1] - 2,(int) dim[2] - 2 };
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

					temp.at<float>(0, i - radius, j - radius) += ((t > slice.at<float>(0, i, j)) && (abs(t - slice.at<float>(0, i, j)) > std::numeric_limits<float>::epsilon())) << n;
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
		ranges[1] = cv::Range(1, dim[1] - 1);
		ranges[2] = cv::Range(1, dim[2] - 1);
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
	hist_map = hist_mapping(roi, lbp_blk);
	hist_map = hist(hist_map, map_table_ri.N);
  lbp_feature_ri.insert(lbp_feature_ri.end(),hist_map.begin(),hist_map.end());

	//write to file 
	hist_map = hist_mapping(roi, lbp_blk_riu2);
	hist_map = hist(hist_map, map_table_riu2.N);
  lbp_feature_riu2.insert(lbp_feature_riu2.end(), hist_map.begin(), hist_map.end());

  for (int y = 0; y <(int)dim[1]; y++){
		int slice_sy[] = {(int) dim[0], 1, (int)dim[2] };
		cv::Mat slice = cv::Mat(3, slice_sy, CV_32F, cv::Scalar(0.0));
		imageblk.col(y).copyTo(slice);
		//std::cout << slice.at<float>(0, 0, 0) << ",," << slice.at<float>(0, 0, 1) << slice.at<float>(0, 0, 2) << slice.at<float>(3, 0, 2);
    int lbp_sizey[] = { (int)dim[0] - 2, 1, (int)dim[2] - 2 };
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

					temp.at<float>(i - radius, 0, j - radius) += ((t > slice.at<float>(i, 0, j)) && (abs(t - slice.at<float>(i, 0, j)) > std::numeric_limits<float>::epsilon())) << n;
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
		ranges[0] = cv::Range(1, dim[0] - 1);
		ranges[2] = cv::Range(1, dim[2] - 1);
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

	hist_map = hist_mapping(roi, lbp_blk);
	hist_map = hist(hist_map, map_table_ri.N);
  lbp_feature_ri.insert(lbp_feature_ri.end(), hist_map.begin(), hist_map.end());
	//write to file 
	hist_map = hist_mapping(roi, lbp_blk_riu2);
	hist_map = hist(hist_map, map_table_riu2.N);
  lbp_feature_riu2.insert(lbp_feature_riu2.end(), hist_map.begin(), hist_map.end());


  for (int z = 0; z < (int)dim[2]; z++){
		int slice_sz[] = { (int)dim[0], (int)dim[1], 1 };
		cv::Mat slice = cv::Mat(3, slice_sz, CV_32F, cv::Scalar(0.0));
		cv::Range ranges[3];
		ranges[0] = cv::Range::all();
		ranges[1] = cv::Range::all();
		ranges[2] = cv::Range(z, z + 1);

		slice = imageblk(ranges).clone();
	//	std::cout << slice.at<float>(0, 0, 0) << ",," << slice.at<float>(0, 1, 0) << slice.at<float>(0, 2, 0) << slice.at<float>(3, 2, 0);
    int lbp_sizez[] = { (int)dim[0] - 2, (int)dim[1] - 2, 1 };
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

					temp.at<float>(i - radius, j - radius, 0) += ((t > slice.at<float>(i, j, 0)) && (abs(t - slice.at<float>(i, j, 0)) > std::numeric_limits<float>::epsilon())) << n;
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
		range[0] = cv::Range(1, dim[0] - 1);
		range[1] = cv::Range(1, dim[1] - 1);
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

	hist_map = hist_mapping(roi, lbp_blk);
	hist_map = hist(hist_map, map_table_ri.N);
  lbp_feature_ri.insert(lbp_feature_ri.end(), hist_map.begin(), hist_map.end());
	//write to file 
	hist_map = hist_mapping(roi, lbp_blk_riu2);
	hist_map = hist(hist_map, map_table_riu2.N);
  lbp_feature_riu2.insert(lbp_feature_riu2.end(), hist_map.begin(), hist_map.end());
  std::stringstream ss;


  for (int i = 0; i <(int)( lbp_feature_ri.size()); i++)
  {

    ss << i;
    std::string i_str = ss.str();
    lbp_featurename->push_back("lbp_ri" + i_str);
  }
  for (int i = 0; i <(int)lbp_feature_riu2.size(); i++)
  {
    ss << i;
    std::string i_str = ss.str();
    lbp_featurename->push_back("lbp_riu2_" + i_str);
  }


  lbp_feature->insert(lbp_feature->end(),lbp_feature_ri.begin(),lbp_feature_ri.end());
  lbp_feature->insert(lbp_feature->end(), lbp_feature_riu2.begin(), lbp_feature_riu2.end());

}