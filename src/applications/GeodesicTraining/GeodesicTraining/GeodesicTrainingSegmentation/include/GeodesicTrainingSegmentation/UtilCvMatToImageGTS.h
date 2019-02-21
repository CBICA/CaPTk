#ifndef H_CBICA_CV_MAT_TO_IMAGE_GTS
#define H_CBICA_CV_MAT_TO_IMAGE_GTS

#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkExceptionObject.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"

#include <opencv2/ml.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

namespace GeodesicTrainingSegmentation
{
	namespace CvMatToImageGTS
	{
		/**
		Fills an image with the contents of a cv::Mat
		@param image the itk image that will be filled
		@param mat the input cv::mat, only the first item of each row will be used
		*/
		template <typename TPixelType, unsigned int TDimensions>
		void FillImage(const typename itk::Image<TPixelType, TDimensions>::Pointer &image, const cv::Mat &mat)
		{
			itk::ImageRegionIteratorWithIndex<itk::Image<TPixelType, TDimensions>> iter_i(image, image->GetRequestedRegion());

			int i;
			for (iter_i.GoToBegin(), i = 0; !iter_i.IsAtEnd(); ++iter_i, i++)
			{
				iter_i.Set(mat.at<TPixelType>(i, 0));
			}
		}
	}
}

#endif // !H_CBICA_CV_MAT_TO_IMAGE_GTS
