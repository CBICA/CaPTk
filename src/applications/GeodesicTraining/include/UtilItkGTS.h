#ifndef H_CBICA_ITK_UTIL_GTS
#define H_CBICA_ITK_UTIL_GTS

#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkExceptionObject.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkNormalizeImageFilter.h"

namespace GeodesicTrainingSegmentation
{
	namespace ItkUtilGTS
	{
		/**Normalizes image values to be inside [0,255]*/
		template<class UImageType>
		typename UImageType::Pointer normalizeImage(typename UImageType::Pointer image, typename UImageType::PixelType max = 255)
		{
			auto filter = itk::RescaleIntensityImageFilter< UImageType, UImageType >::New();
			filter->SetInput(image);
			filter->SetOutputMinimum(0);
			filter->SetOutputMaximum(max);
			filter->Update();
			return filter->GetOutput();
		}

		/**Normalizes image values to be inside [0,255]*/
		template<class UImageType>
		void normalizeImageVector(std::vector<typename UImageType::Pointer>& images)
		{
			for (size_t i = 0; i < images.size(); i++) {
				images[i] = normalizeImage<UImageType>(images[i]);
			}
		}

		/**Normalizes image values to be inside [0,255]*/
		template<typename KeysType, class UImageType>
		void normalizeImageMap(std::map< KeysType, typename UImageType::Pointer>& map)
		{
			if (map.size() == 0) { return; }
			
			// Find a list of the different keys
			std::vector<KeysType> keys;
			keys.reserve(map.size());
			for (auto const& imap : map) {
				keys.push_back(imap.first);
			}

			// Normalize
			for (auto key : keys) {
				map[key] = normalizeImage<UImageType>(map[key]);
			}
		}

		/**Calculates the maximum pixel value of an image*/
		template<class UImageType>
		typename UImageType::PixelType getImageMaximum(typename UImageType::Pointer image)
		{
			auto imageCalculatorFilter = itk::MinimumMaximumImageCalculator<UImageType>::New();
			imageCalculatorFilter->SetImage(image);
			imageCalculatorFilter->ComputeMaximum();
			return imageCalculatorFilter->GetMaximum();
		}

		/**Changes the type of an image and rescales it to be inside [0,255]*/
		template<class TImageType, class UImageType>
		typename UImageType::Pointer castAndRescaleImage(typename TImageType::Pointer input)
		{
			typedef itk::CastImageFilter<TImageType, UImageType> CastFilterType;

			typename CastFilterType::Pointer castFilter = CastFilterType::New();
			castFilter->SetInput(input);

			auto filter = itk::RescaleIntensityImageFilter< UImageType, UImageType >::New();
			filter->SetInput(castFilter->GetOutput());
			filter->SetOutputMinimum(0);
			filter->SetOutputMaximum(255);
			filter->Update();
			return filter->GetOutput();
		}

		/**Create a zero image with exactly the same dimensions as another image*/
		template <class TImageType, class UImageType = TImageType>
		typename UImageType::Pointer initializeOutputImageBasedOn(typename TImageType::Pointer image)
		{
			typename UImageType::Pointer res = UImageType::New();

			res->SetRegions(image->GetLargestPossibleRegion());
			res->SetRequestedRegion(image->GetRequestedRegion());
			res->SetBufferedRegion(image->GetBufferedRegion());
			res->Allocate();
			res->FillBuffer(0);
			res->SetDirection(image->GetDirection());
			res->SetOrigin(image->GetOrigin());
			res->SetSpacing(image->GetSpacing());

			return res;
		}

		/**Normalize an image by setting its mean to 0 and variance to 1*/
		template<class UImageType, class VImageType = UImageType>
		typename VImageType::Pointer statisticalImageNormalization(typename UImageType::Pointer image, typename VImageType::PixelType normMax = 255)
		{
			auto normalizeFilter = itk::NormalizeImageFilter<UImageType, VImageType>::New();
			normalizeFilter->SetInput(image);
			normalizeFilter->Update();
			return normalizeImage<VImageType>(normalizeFilter->GetOutput(), normMax);
		}

		/**Normalize an image by setting its mean to 0 and variance to 1*/
		template<class UImageType, class VImageType = UImageType>
		void statisticalImageVectorNormalization(std::vector<typename UImageType::Pointer>& images, typename VImageType::PixelType normMax = 255)
		{
			for (size_t i = 0; i < images.size(); i++) {
				images[i] = statisticalImageNormalization<UImageType, VImageType>(images[i], normMax);
			}
		}

		/**Normalize an image by setting its mean to 0 and variance to 1*/
		template<typename KeysType, class UImageType, class VImageType = UImageType>
		void statisticalImageMapNormalization(std::map< KeysType, typename UImageType::Pointer>& map, typename VImageType::PixelType normMax = 255)
		{
			if (map.size() == 0) { return; }

			// Find a list of the different keys
			std::vector<KeysType> keys;
			keys.reserve(map.size());
			for (auto const& imap : map) {
				keys.push_back(imap.first);
			}

			// Statistical Normalize
			for (auto key : keys) {
				map[key] = statisticalImageNormalization<UImageType, VImageType>(map[key], normMax);
			}
		}

	}
}

#endif // !H_CBICA_ITK_UTIL_GTS