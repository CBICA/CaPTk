#ifndef H_CBICA_PROCESSING
#define H_CBICA_PROCESSING

#include <itkImageRegionIteratorWithIndex.h>
#include <itkNeighborhoodIterator.h>

#include <vector>
#include <set>
#include <cmath>

#include "cbicaITKSafeImageIO.h"
#include "cbicaITKImageInfo.h"

#include "SusanDenoising.h"
#include "UtilItkGTS.h"
#include "SvmSuite.h"

namespace GeodesicTrainingSegmentation
{
	/** This class handles preprocessing and postprocessing operations */
	template<typename PixelType, unsigned int Dimensions>
	class Processing
	{
	public:
		typedef int   AgdPixelType;
		typedef int   LabelsPixelType;
		typedef float PseudoProbPixelType;

		typedef itk::Image< AgdPixelType, Dimensions >        AgdImageType;
		typedef itk::Image< PixelType, Dimensions >           InputImageType;
		typedef itk::Image< LabelsPixelType, Dimensions >     LabelsImageType;
		typedef itk::Image< PseudoProbPixelType, Dimensions > PseudoProbImageType;

		typedef typename AgdImageType::Pointer                AgdImagePointer;
		typedef typename InputImageType::Pointer              InputImagePointer;
		typedef typename LabelsImageType::Pointer             LabelsImagePointer;
		typedef typename PseudoProbImageType::Pointer         PseudoProbImagePointer;

		void SetLimitPixels(bool limitPixels, int pixelLimit = 10000000)
		{
			m_limit_pixels = limitPixels;
			m_pixel_limit  = pixelLimit;
		}

		void SetDoStatisticalNormalization(bool doStatisticalNormalization, float imageToAgdMapsRatio = 6)
		{
			m_do_statistical_normalization = doStatisticalNormalization;
			m_image_to_agd_maps_ratio = imageToAgdMapsRatio;
		}

		void SetDoCurvatureAnisotropic(bool doCurvatureAnisotropic)
		{
			m_do_curvature_anisotropic = doCurvatureAnisotropic;
		}

		void SetDoSusanDenoising(bool doSusanDenoising)
		{
			m_do_susan_denoising = doSusanDenoising;
		}

		void SetVerbose(bool verbose)
		{
			m_verbose = verbose;
		}

		void SetSaveAll(bool saveAll)
		{
			m_save_all = saveAll;
		}

		void SetOutputFolder(std::string outputFolder)
		{
			m_output_folder = outputFolder;
		}

		void SetTimerEnabled(bool timerEnabled) 
		{
			m_timer_enabled = timerEnabled;
		}

		void SetNumberOfThreads(int numberOfThreads)
		{
			m_number_of_threads = numberOfThreads;
		}

		void PreProcess(std::vector<InputImagePointer>& inputImages, LabelsImagePointer& labels)
		{
			if (m_verbose) { std::cout << "Preprocessing:\n"; }

			m_original_input_image_size     = inputImages[0]->GetLargestPossibleRegion().GetSize();
			m_original_labels_image_size    = labels->GetLargestPossibleRegion().GetSize();
			m_original_input_image_spacing  = inputImages[0]->GetSpacing();
			m_original_labels_image_spacing = labels->GetSpacing();
			m_original_input_image_origin  = inputImages[0]->GetOrigin();
			m_original_labels_image_origin = labels->GetOrigin();

			/* ------ Find the pixel count of the original image ------ */
			
			unsigned long originalPixelCount = 1;
			for (unsigned int i=0; i < InputImageType::ImageDimension; i++) {
				originalPixelCount *= m_original_input_image_size[i];
			}

			/* ------ Print input information ------ */

			if (m_verbose)
			{
				std::cout << "\tImage size:     " << m_original_input_image_size << "\n";
				std::cout << "\tLabels size:    " << m_original_labels_image_size << "\n";
				std::cout << "\tImage spacing:  " << m_original_input_image_spacing << "\n";
				std::cout << "\tLabels spacing: " << m_original_labels_image_spacing << "\n";
				std::cout << "\tPixel count:    " << originalPixelCount << "\n";
				std::cout << "\tPixel limit:    " << ((m_limit_pixels)? 
					std::to_string(m_pixel_limit) : 
					"none") << "\n";
				std::cout << "\t---\n";
			}

			/* ------ Find out if spacing should be normalized ------ */

			// Find minimum, maximum and average spacing
			auto   minSpacing = m_original_input_image_spacing[0];
			auto   maxSpacing = m_original_input_image_spacing[0];
			double avgSpacing = m_original_input_image_spacing[0];
			for(int i=1; i<InputImageType::ImageDimension; i++)
			{
				if (m_original_input_image_spacing[i] > maxSpacing) {
					maxSpacing = m_original_input_image_spacing[i];
				} 
				if (m_original_input_image_spacing[i] < minSpacing) {
					minSpacing = m_original_input_image_spacing[i];
				} 
				avgSpacing += m_original_input_image_spacing[i];
			}
			avgSpacing /= InputImageType::ImageDimension;

			// This variable holds the normal image size if spacing will not be changed,
			// or holds the value of the new size that will result from (theoretically) changing the spacing.
			typename InputImageType::SizeType tSize = m_original_input_image_size;

			// Find if spacing should be normalized (i.e. there is a big difference across dimensions)
			bool shouldNormalizeSpacing = false;
			for(int i=0; i<InputImageType::ImageDimension; i++)
			{
				auto dimSpacingDiff = m_original_input_image_spacing[i] - minSpacing;
				if (dimSpacingDiff > 0.1 || dimSpacingDiff < -0.1)
				{
					shouldNormalizeSpacing = true;

					// Find how the size would look if we (theoretically) just normalized the spacing
					for (unsigned int i=0; i < InputImageType::ImageDimension; i++) {
						tSize[i] = m_original_input_image_size[i] * (
							static_cast<double>(m_original_input_image_spacing[i]) / minSpacing
						);
					}
					break;
				}
			}
			// Check for weird case where a 2D image is actually 3D (with one dimension with size 1)
			if (InputImageType::ImageDimension == 3)
			{
				bool sneaky2D = false;
				unsigned int occursInDim;
				for(int i=0; i<3; i++)
				{
					if (m_original_input_image_size[i] == 1)
					{
						sneaky2D = true;
						occursInDim = i;
						break;
					}
				}
				if (sneaky2D)
				{
					float spacing1 = -1, spacing2 = -1;
					for (int i=0; i<3; i++)
					{
						if (i == occursInDim) { continue; }
						if (spacing1 == -1) { spacing1 == m_original_input_image_spacing[i]; }
						else { spacing2 == m_original_input_image_spacing[i]; }
					}
					float dimSpacingDiff = spacing1 - spacing2;
					if (dimSpacingDiff > 0.1 || dimSpacingDiff < -0.1)
					{
						shouldNormalizeSpacing = true; // Keep the previous tSize for convinience
					}
					else {
						shouldNormalizeSpacing = false;
					}
				}
			}

			// Find the (theoretical) size that will result from potentially normalizing spacing
			unsigned long tPixelCount = 1;
			for(int i=0; i<InputImageType::ImageDimension; i++) {
				tPixelCount *= tSize[i];
			}

			// This variable holds the original spacing if spacing will not be changed
			// or the (theoretical) spacing if we just normalized the spacing
			typename InputImageType::SpacingType tSpacing;
			if (shouldNormalizeSpacing) {
				tSpacing.Fill(minSpacing);
			}
			else {
				tSpacing = m_original_input_image_spacing;
			}

			/* ------ Find out if size should be normalized (after potential spacing changes) ------ */
			
			bool shouldNormalizeSize = false;
			if (m_limit_pixels && tPixelCount > m_pixel_limit) 
			{
				shouldNormalizeSize = true; 
			}

			/* ------ Find target size and spacing ------ */

			typename InputImageType::SizeType targetSize;
			typename InputImageType::SpacingType targetSpacing = tSpacing;

			// Calculate target size if there is (theoretically) no size normalization
			for(int i=0; i<InputImageType::ImageDimension; i++) 
			{
				targetSize[i] = std::lround(
					1.0 * tSize[i] * m_original_labels_image_spacing[i] / tSpacing[i]
				);
				if (targetSize[i] < 1) { targetSize[i] = 1; }
			}

			if (shouldNormalizeSize) 
			{
				double ratio = 1.0 * m_pixel_limit / tPixelCount;

				typename InputImageType::SizeType targetSizeSmall;

				for(int i=0; i<InputImageType::ImageDimension; i++) 
				{
					targetSizeSmall[i] = std::lround(
						1.0 * targetSize[i] * std::pow(ratio, 1.0/InputImageType::ImageDimension)
					);

					targetSpacing[i] = std::lround(
						1.0 * tSpacing[i] * targetSize[i] / targetSizeSmall[i]
					);
				}

				targetSize = targetSizeSmall;
			}

			/* ------ Find target size and spacing for labels ------ */

			typename InputImageType::SizeType targetSizeLabels;
			typename InputImageType::SpacingType targetSpacingLabels;

			for(int i=0; i<InputImageType::ImageDimension; i++) 
			{
				targetSizeLabels[i] = targetSize[i];
				targetSpacingLabels[i] = targetSpacing[i];
			}

			/* ------ Resample the images and the labels image ------ */

			if (shouldNormalizeSize || shouldNormalizeSpacing)
			{			
				startTimer();
				for (auto& image : inputImages)
				{
					image = ItkUtilGTS::resampleNormalImage<InputImageType>(
						image, targetSize, targetSpacing
					);
				}
				stopTimerAndReport("Resampling images");

				startTimer();
				labels = ItkUtilGTS::resampleLabelsImage<LabelsImageType>(
					labels, targetSizeLabels, targetSpacingLabels, m_verbose
				);
				stopTimerAndReport("Resampling labels");
			}

			/* ------ Print image information if anything changed ------ */

			if ((shouldNormalizeSize || shouldNormalizeSpacing) && m_verbose)
			{
				std::cout << "\t After normalizing spacing and size:\n";
				std::cout << "\t\t Size:    " << labels->GetLargestPossibleRegion().GetSize() << "\n";
				std::cout << "\t\t Spacing: " << labels->GetSpacing() << "\n";
			}

			/* Save the intermediate images if saveall is set */
			
			if (m_save_all)
			{
				int y = 0;
				for (auto& image : inputImages) 
				{
					cbica::WriteImage<InputImageType>( image, m_output_folder + "/" +
						std::string("preprocessed_image_just_size_and_spacing") + std::to_string(++y) + ".nii.gz"
					);
				}		
				cbica::WriteImage<LabelsImageType>( labels, 
					m_output_folder + std::string("/") + "preprocessed_labels_just_size_and_spacing.nii.gz"
				);
			}

			/* ------ Normal filters after size and spacing are settled ------ */

			if (m_do_statistical_normalization)
			{
				startTimer();
				if (m_verbose) { std::cout << "\tNormalizing images\n"; }

				// Statistical Image Normalization is fast, no need for threads
				ItkUtilGTS::statisticalImageVectorNormalization<InputImageType>(
					inputImages, 
					std::lround(255 * m_image_to_agd_maps_ratio)
				);
				stopTimerAndReport("Preprocessing: Statistical image normalization");

				if (m_save_all)
				{		
					int w = 0;		
					for (auto& image : inputImages)
					{
						cbica::WriteImage<InputImageType>( image, m_output_folder + "/" +
							std::string("statistical_norm_image") + std::to_string(++w) + std::string(".nii.gz")
						);
					}
				}
			}

			// {
			// 	startTimer();
			// 	if (m_verbose) { std::cout << "\tBilateral filter\n"; }
			// 	for (auto& image : inputImages)
			// 	{
			// 		std::cout << "new\n";
			// 		image = ItkUtilGTS::BilateralImageFilter<InputImageType>(
			// 			image, true, 255 * m_image_to_agd_maps_ratio
			// 		);

			// 		cbica::WriteImage<InputImageType>( image, m_output_folder + "/" +
			// 			std::string("bilateral_norm_image_ignore") + std::string(".nii.gz")
			// 		);
			// 	}
			// 	stopTimerAndReport("Bilateral image filter");

			// 	if (m_save_all)
			// 	{		
			// 		int w = 0;		
			// 		for (auto& image : inputImages)
			// 		{
			// 			cbica::WriteImage<InputImageType>( image, m_output_folder + "/" +
			// 				std::string("bilateral_norm_image") + std::to_string(++w) + std::string(".nii.gz")
			// 			);
			// 		}
			// 	}				
			// }

			// {
			// 	startTimer();
			// 	if (m_verbose) { std::cout << "\tGADF filter\n"; }
			// 	for (auto& image : inputImages)
			// 	{
			// 		std::cout << "new\n";
			// 		image = ItkUtilGTS::gradientAnisotropicDiffusionImageFilter<InputImageType>(
			// 			image, true, 255 * m_image_to_agd_maps_ratio
			// 		);

			// 		cbica::WriteImage<InputImageType>( image, m_output_folder + "/" +
			// 			std::string("gadf_norm_image_ignore") + std::string(".nii.gz")
			// 		);
			// 	}
			// 	stopTimerAndReport("GADF image filter");

			// 	if (m_save_all)
			// 	{		
			// 		int w = 0;		
			// 		for (auto& image : inputImages)
			// 		{
			// 			cbica::WriteImage<InputImageType>( image, m_output_folder + "/" +
			// 				std::string("gadf_norm_image") + std::to_string(++w) + std::string(".nii.gz")
			// 			);
			// 		}
			// 	}				
			// }

			// Susan denoising
			if (m_do_susan_denoising) {
				startTimer();
				// TODO: Background threads
				if (m_verbose) { std::cout << "\tDenoising images\n"; }
				SusanDenoising susan;
				for (int i = 0; i < inputImages.size(); i++)
				{
					inputImages[i] = susan.Run<InputImageType>(inputImages[i]);
				}
				stopTimerAndReport("Preprocessing: Susan Denoising Filter");

				// Optionally save
				int z = 0;
				if (m_save_all)
				{				
					for (auto& image : inputImages) {
						cbica::WriteImage<InputImageType>( image, m_output_folder + "/" +
							std::string("susan_denoised_image") + std::to_string(++z) + std::string(".nii.gz")
						);
					}
				}
			}

			if (m_do_curvature_anisotropic)
			{
				startTimer();
				if (m_verbose) { std::cout << "\tApplying CADF... [Warning: This might take a long time]\n"; }
				
				int counterForThreadsVec = 0;
				std::vector<std::thread> threads(inputImages.size());
				int numberOfOpenThreads = 0;
				int oldestOpenThread = 0;

				for (int i = 0; i < inputImages.size(); i++)
				{
					if (numberOfOpenThreads == m_number_of_threads) {
						threads[oldestOpenThread].join();
						oldestOpenThread++;
						numberOfOpenThreads--;
					}

					numberOfOpenThreads++;
					threads[counterForThreadsVec++] = std::thread( 
						[&](int ii, std::vector<InputImagePointer>& inputImages) mutable 
						{
							inputImages[ii] = ItkUtilGTS::curvatureAnisotropicDiffusionImageFilter<InputImageType>(
								inputImages[ii]
							);
						},
						i,
						std::ref(inputImages)
					);
				}

				for (int i = oldestOpenThread; i < inputImages.size(); i++) {
					threads[i].join();
				}
				stopTimerAndReport("Preprocessing: Curvature Anisotropic Diffusion Image Filter");

				// Optionally save
				int z = 0;
				if (m_save_all)
				{				
					for (auto& image : inputImages) {
						// image = ItkUtilGTS::curvatureAnisotropicDiffusionImageFilter<InputImageType>(image);
						cbica::WriteImage<InputImageType>( image, m_output_folder + "/" +
							std::string("cadf_image") + std::to_string(++z) + std::string(".nii.gz")
						);
					}
				}
			}
		}

		void PostProcessLabelsImage(LabelsImagePointer& labels)
		{
			// Make the output the same as the original input
			if (m_verbose) { std::cout << "Postprocessing labels image...\n"; }

			if (m_save_all)
			{
				cbica::WriteImage<LabelsImageType>( labels, m_output_folder + "/" +
					std::string("labels_before_post_processing.nii.gz")
				);
			}

			startTimer();
			labels = ItkUtilGTS::changeImageSpacing<LabelsImageType>(
				labels, m_original_labels_image_spacing, true, 
				true, m_original_labels_image_size,
				true, m_original_labels_image_origin
			);
			stopTimerAndReport("Postprocessing: Changing size and spacing back");
			
			if (m_save_all)
			{
				cbica::WriteImage<LabelsImageType>( labels, m_output_folder + "/" +
					std::string("labels_res.nii.gz")
				);
			}
		}

		template <class TImageType>
		void PostProcessNormalImages(std::vector<typename TImageType::Pointer>& images)
		{
			for (auto& image : images) {
				image = ItkUtilGTS::changeImageSpacing<TImageType>(
					image, m_original_input_image_spacing, false, 
					true, m_original_input_image_size
				);	
			}
		}

		template <class TImageType>
		void PostProcessNormalImage(typename TImageType::Pointer& image)
		{
			if (m_verbose) { std::cout << "Postprocessing image\n"; }
			startTimer();
			std::vector<typename TImageType::Pointer> imageInVector;
			imageInVector.push_back(image);
			PostProcessNormalImages<TImageType>(imageInVector);
			stopTimerAndReport("Postprocessing: Changing a non-label image's size and spacing back");
		}

	private:
		typename InputImageType::SpacingType  m_original_input_image_spacing;
		typename LabelsImageType::SpacingType m_original_labels_image_spacing;

		typename InputImageType::SizeType  m_original_input_image_size;
		typename LabelsImageType::SizeType m_original_labels_image_size;

		typename InputImageType::PointType  m_original_input_image_origin;
		typename LabelsImageType::PointType m_original_labels_image_origin;

		bool  m_limit_pixels = true, m_do_curvature_anisotropic = false, 
		      m_do_susan_denoising = false, m_do_statistical_normalization = true,
		      m_verbose = false, m_save_all = false, m_timer_enabled = false;
		int   m_pixel_limit  = 10000000, m_number_of_threads = 16;
		float m_image_to_agd_maps_ratio = 6;
		std::string m_output_folder = cbica::getExecutablePath();
		SvmSuiteUtil::Timer m_timer;

		void startTimer() {
			if (m_timer_enabled) {
				m_timer.Reset();
			}
		}

		void stopTimerAndReport(std::string desc) {
			if (m_timer_enabled) {
				float diff = m_timer.Diff();

				std::ofstream timerFile;
				timerFile.open(m_output_folder + "/time_report.txt", std::ios_base::app); //append file
				timerFile << desc << ": " << diff << "s\n";
			}
		}

	};
}

#endif // ! H_CBICA_PROCESSING
