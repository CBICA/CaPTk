#ifndef  H_CBICA_GEODESIC_TRAINING_SEGMENTATION
#define  H_CBICA_GEODESIC_TRAINING_SEGMENTATION

#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkExceptionObject.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"
#include "itkScaleTransform.h"

#include "cbicaCmdParser.h"
#include "cbicaLogging.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaITKImageInfo.h"
#include "cbicaUtilities.h"

#include <iostream>
#include <fstream>
#include <typeinfo>
#include <exception>
#include <cmath>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <thread>
#include <mutex>

#include "UtilGTS.h"
#include "UtilItkGTS.h"
#include "UtilImageToCvMatGTS.h"
#include "UtilCvMatToImageGTS.h"
#include "OperationsSvmGTS.h"

#include "AdaptiveGeodesicDistance.h"
#include "SvmSuite.h"
#include "RandomForestSuite.h"
#include "Processing.h"

namespace GeodesicTrainingSegmentation
{
	enum MODE
	{
		SVM_PSEUDO,
		SVM_LABELS,
		AGD,
		LABELS_THRESHOLD,
		SEGMENT,
		REVERSE_GEOTRAIN,
		REVERSE_GEOTRAIN_SPORADIC,
		GEOTRAIN,
		GEOTRAIN_FULL,
		RF,
		AGD_RF,
		RF_AUTO,
		AGD_RF_AUTO,
		CHANGE_LABELS,
		GENERATE_CONFIG,
		CHECK_ACCURACY
	};

	enum MODALITY_MRI
	{
		FLAIR,
		T1,
		T1CE,
		T2
	};

	typedef int   AgdPixelType;
	typedef int   LabelsPixelType;
	typedef float PseudoProbPixelType;

	const LabelsPixelType DEFAULT_LABEL_TC = 1; // Tumor core
	const LabelsPixelType DEFAULT_LABEL_ET = 4; // Enhanced tumor
	const LabelsPixelType DEFAULT_LABEL_ED = 2; // Edema
	const LabelsPixelType DEFAULT_LABEL_HT = 3; // Healthy tissue

	const LabelsPixelType DEFAULT_LABEL_OF_INTEREST = 1;
	
	const std::string   DEFAULT_FILE_EXTENSION      = ".nii.gz";
	const float         DEFAULT_THRESHOLD           = 25;
	const int           DEFAULT_NUMBER_OF_THREADS   = 16;
	const int           MAX_THREADS_SVM_SUITE       = 32;
	const long          DEFAULT_MAX_IMAGE_SIZE_BEFORE_RESAMPLING = 10000000;

	const int           DEFAULT_MAX_SAMPLES_SVM_SUBSAMPLE      = 3000;
	const bool          DEFAULT_BALANCED_SUBSAMPLE             = true;
	const float         DEFAULT_INPUT_IMAGES_TO_AGD_MAPS_RATIO = 6;

	/** Class for all the Geodesic Training Segmentation operations */
	template<typename PixelType = float, unsigned int Dimensions = 3>
	class Coordinator
	{
	public:

		typedef itk::Image< AgdPixelType, Dimensions >        AgdImageType;
		typedef itk::Image< PixelType, Dimensions >           InputImageType;
		typedef itk::Image< LabelsPixelType, Dimensions >     LabelsImageType;
		typedef itk::Image< PseudoProbPixelType, Dimensions > PseudoProbImageType;

		typedef typename AgdImageType::Pointer                AgdImagePointer;
		typedef typename InputImageType::Pointer              InputImagePointer;
		typedef typename LabelsImageType::Pointer             LabelsImagePointer;
		typedef typename PseudoProbImageType::Pointer         PseudoProbImagePointer;

		/** Everything that the class does gets returned with this.
		 *  Not every field gets used every time
		 */
		typedef struct Result 
		{
			PseudoProbImagePointer posImage;
			PseudoProbImagePointer negImage;
			LabelsPixelType posLabel = 0;
			LabelsPixelType negLabel = 0;
			LabelsImagePointer labelsImage;
			InputImagePointer  segmentedFloatImage;
			AgdImagePointer    agdMapImage;

			std::map< LabelsPixelType, double > diceScore;
			std::map< LabelsPixelType, double > sensitivity;
			std::map< LabelsPixelType, int    > falsePositivesCount;

			double diceScoreAll;
			double sensitivityAll;
			int    falsePositivesCountAll;

			bool ok = true;
			std::string errorMessage = "";
		} Result;

		explicit Coordinator() {}

		virtual ~Coordinator() {}

		static_assert((Dimensions == 2 || Dimensions == 3), "2D or 3D Images supported");
		
		/** Main execution function */
		std::shared_ptr< Result > Execute()
		{
			std::shared_ptr< Result > gtsResult(new Result());

			// Check if everything that is needed is supplied
			if (!validState(gtsResult)) { return gtsResult; }

			message(
                std::string("-----------------------------------------\n") +
                std::string("Geodesic Training - Semantic Segmentation\n") +
                std::string("-----------------------------------------\n\n"), 
                "Executing"
			);	
			
			startTimer();

			//shrinkImagesIfNecessary();

			if (!m_changed_labels_map_manually_set && m_input_images_MRI.size() > 0) {
				// By default, remove healthy tissue after calculations
				m_change_labels_map[m_label_HT] = 0;
			}

			// // Normalize input (if applicable). Note: AGD maps are capped at 255
			// ItkUtilGTS::statisticalImageVectorNormalization<InputImageType>(m_input_images, std::lround(255 * m_image_to_agd_maps_ratio));
			// //ItkUtilGTS::statisticalImageMapNormalization<MODALITY_MRI, InputImageType>(m_input_images_MRI, std::lround(255 * m_image_to_agd_maps_ratio));
			m_processing.SetVerbose(m_verbose);
			m_processing.SetSaveAll(m_save_all);
			m_processing.SetOutputFolder(m_output_folder);
			m_processing.SetTimerEnabled(m_timer_enabled);
			m_processing.SetNumberOfThreads((m_max_threads)?32:m_number_of_threads);

			switch (m_mode) {
			case SVM_PSEUDO:
			{
				message("Number of modalities: " + std::to_string(m_input_images.size() + m_input_images_MRI.size()) + "\n");
				m_processing.PreProcess(m_input_images, m_labels_image);
				auto labelsCountMap = ParserGTS::CountsOfEachLabel<LabelsImageType>(m_labels_image);
				
				if (!validLabels(gtsResult, labelsCountMap, 2)) { return gtsResult; }
				
				addToInputImagesMRI(); // Nothing special happens for the different modalities here
				
				addCoordinatorMapsIfNecessary();

				// Construct pseudoprobability maps (~distance to hyperplane)
				svm<PixelType>(m_input_images, m_labels_image, gtsResult, labelsCountMap, true);

				m_processing.template PostProcessNormalImage<PseudoProbImageType>(gtsResult->posImage);
				m_processing.template PostProcessNormalImage<PseudoProbImageType>(gtsResult->negImage);

			} break;
			case SVM_LABELS:
			{
				message("Number of modalities: " + std::to_string(m_input_images.size() + m_input_images_MRI.size()) + "\n");
				m_processing.PreProcess(m_input_images, m_labels_image);
				auto labelsCountMap = ParserGTS::CountsOfEachLabel<LabelsImageType>(m_labels_image);

				if (!validLabels(gtsResult, labelsCountMap)) { return gtsResult; }

				addToInputImagesMRI(); // Nothing special happens for the different modalities here

				addCoordinatorMapsIfNecessary();

				// Construct a labeled image using SVM(s)
				svm<PixelType>(m_input_images, m_labels_image, gtsResult, labelsCountMap, false);
				if (!gtsResult->ok) { break; }

				m_processing.PostProcessLabelsImage(gtsResult->labelsImage);

				// Change labels if necessary
				changeLabels(gtsResult->labelsImage, m_change_labels_map);

				if (m_ground_truth_set) {
					checkAccuracyInRelationToGroundTruth(gtsResult->labelsImage, gtsResult);
				}

			} break;
			case AGD:
			{
				message("Number of modalities: " + std::to_string(m_input_images.size() + m_input_images_MRI.size()) + "\n");
				m_processing.PreProcess(m_input_images, m_labels_image);
				auto labelsCountMap = ParserGTS::CountsOfEachLabel<LabelsImageType>(m_labels_image);
				
				if (!validLabels(gtsResult, labelsCountMap, 1)) { return gtsResult; }

				addToInputImagesMRI(); // Nothing special happens for the different modalities here

				// Construct AGD maps
				gtsResult->agdMapImage = agd<PixelType>(m_input_images[0], m_label_of_interest, true);

				m_processing.template PostProcessNormalImage<AgdImageType>(gtsResult->agdMapImage);

				gtsResult->labelsImage = thresholdSingleImageToLabel<AgdImageType>(gtsResult->agdMapImage, m_threshold, m_label_of_interest);

				// Change labels if necessary
				changeLabels(gtsResult->labelsImage, m_change_labels_map);

				if (m_ground_truth_set) {
					checkAccuracyInRelationToGroundTruth(gtsResult->labelsImage, gtsResult);
				}

			} break;
			case LABELS_THRESHOLD:
			{
				// Construct labeled image. (For each pixel: if intensity <= threshold then label of interest, else label 0) 
				gtsResult->labelsImage = thresholdSingleImageToLabel<InputImageType>(m_input_images[0], m_threshold, m_label_of_interest);

			} break;
			case SEGMENT:
			{
				auto labelsCountMap = ParserGTS::CountsOfEachLabel<LabelsImageType>(m_labels_image);

				if (!validLabels(gtsResult, labelsCountMap, 1)) { return gtsResult; }

				addToInputImagesMRI(); // Nothing special happens for the different modalities here

				// Crop an image only at the labels of a particular class, according to labels image
				InputImagePointer segmentedImage = segmentSingleImage<InputImageType>(m_input_images[0], m_labels_image, m_label_of_interest);
				
				if (std::is_same<PixelType, float>::value) {
					// In case more pixel types are supported in the future
					gtsResult->segmentedFloatImage = segmentedImage;
				}

			} break;
			case REVERSE_GEOTRAIN_SPORADIC:
			{
				message("Number of modalities: " + std::to_string(m_input_images.size() + m_input_images_MRI.size()) + "\n");
				m_processing.PreProcess(m_input_images, m_labels_image);
				auto labelsCountMap = ParserGTS::CountsOfEachLabel<LabelsImageType>(m_labels_image);

				if (!validLabels(gtsResult, labelsCountMap)) { return gtsResult; }

				if(m_label_HT == 0 || labelsCountMap.count(m_label_HT) == 0)
				{
					gtsResult->ok = false;
					gtsResult->errorMessage = "Either the healthy tissue label was not set or there weren't any healthy tissue labels drawn";
					errorOccured(gtsResult->errorMessage);
					return gtsResult;
				}

				// Construct agd maps for each input image and input input images and agd maps to a SVM that produces labels
				std::vector< AgdImagePointer > agdImages;

				addToInputImagesMRI(); // No special operations for MRI images
				
				std::vector<LabelsPixelType> htKeyInList;
				htKeyInList.push_back(m_label_HT);

				if (m_input_images.size() != 0) {

					std::vector< AgdImagePointer > agdImagesAllKeys = agd<PixelType>(m_input_images, htKeyInList, true);
					agdImages.insert(agdImages.end(), agdImagesAllKeys.begin(), agdImagesAllKeys.end()); // Append
				}

				message("Converting AGD maps...", "Converting");
				for (auto agdImage : agdImages) {
					// Convert output of AGD to InputImageType
					auto agdImageCast = ItkUtilGTS::castAndRescaleImage< AgdImageType, InputImageType>(agdImage);

					// Apply negative filter so it is 0 outside of the brain (basically where there is no information - that is used to make it faster later)
					typedef itk::InvertIntensityImageFilter<InputImageType> InvertIntensityImageFilterType;
					typename InvertIntensityImageFilterType::Pointer invertIntensityFilter = InvertIntensityImageFilterType::New();
					invertIntensityFilter->SetInput(agdImageCast);
					invertIntensityFilter->SetMaximum(255);
					invertIntensityFilter->Update();

					//m_input_images.push_back(ItkUtilGTS::statisticalImageNormalization<InputImageType>(invertIntensityFilter->GetOutput(), VALUE));
					m_input_images.push_back(invertIntensityFilter->GetOutput());
					m_agd_maps_count++;
				}
				message("Converting AGD maps...", "Converting", 100);

				addCoordinatorMapsIfNecessary();

				// SVM_LABELS
				svm<PixelType>(m_input_images, m_labels_image, gtsResult, labelsCountMap, false);
				if (!gtsResult->ok) { break; }

				m_processing.PostProcessLabelsImage(gtsResult->labelsImage);

				// Change labels if necessary
				changeLabels(gtsResult->labelsImage, m_change_labels_map);

				if (m_ground_truth_set) {
					checkAccuracyInRelationToGroundTruth(gtsResult->labelsImage, gtsResult);
				}

			} break;
			case REVERSE_GEOTRAIN:
			{
				message("Number of modalities: " + std::to_string(m_input_images.size() + m_input_images_MRI.size()) + "\n");
				m_processing.PreProcess(m_input_images, m_labels_image);
				auto labelsCountMap = ParserGTS::CountsOfEachLabel<LabelsImageType>(m_labels_image);

				if (!validLabels(gtsResult, labelsCountMap)) { return gtsResult; }

				// Construct agd maps for each input image and input input images and agd maps to a SVM that produces labels
				std::vector< AgdImagePointer > agdImages;

				if (m_input_images.size() != 0) {
					auto keys = getMapKeyset(labelsCountMap);

					std::vector< AgdImagePointer > agdImagesAllKeys = agd<PixelType>(m_input_images, keys, true);
					agdImages.insert(agdImages.end(), agdImagesAllKeys.begin(), agdImagesAllKeys.end()); // Append
				}

				if (m_input_images_MRI.size() != 0) {
					// Special operations for MRI images
					std::vector< AgdImagePointer > agdImagesMRI = agdMRI<PixelType>(m_input_images_MRI, getMapKeyset(labelsCountMap), false, true);
					agdImages.insert(agdImages.end(), agdImagesMRI.begin(), agdImagesMRI.end()); // Append
					addToInputImagesMRI();
				}

				message("Converting AGD maps...", "Converting");
				for (auto agdImage : agdImages) {
					// Convert output of AGD to InputImageType
					auto agdImageCast = ItkUtilGTS::castAndRescaleImage< AgdImageType, InputImageType>(agdImage);

					// Apply negative filter so it is 0 outside of the brain (basically where there is no information - that is used to make it faster later)
					typedef itk::InvertIntensityImageFilter<InputImageType> InvertIntensityImageFilterType;
					typename InvertIntensityImageFilterType::Pointer invertIntensityFilter = InvertIntensityImageFilterType::New();
					invertIntensityFilter->SetInput(agdImageCast);
					invertIntensityFilter->SetMaximum(255);
					invertIntensityFilter->Update();

					//m_input_images.push_back(ItkUtilGTS::statisticalImageNormalization<InputImageType>(invertIntensityFilter->GetOutput(), VALUE));
					m_input_images.push_back(invertIntensityFilter->GetOutput());
					m_agd_maps_count++;
				}
				message("Converting AGD maps...", "Converting", 100);

				addCoordinatorMapsIfNecessary();

				// SVM_LABELS
				svm<PixelType>(m_input_images, m_labels_image, gtsResult, labelsCountMap, false);
				if (!gtsResult->ok) { break; }

				m_processing.PostProcessLabelsImage(gtsResult->labelsImage);

				// Change labels if necessary
				changeLabels(gtsResult->labelsImage, m_change_labels_map);

				if (m_ground_truth_set) {
					checkAccuracyInRelationToGroundTruth(gtsResult->labelsImage, gtsResult);
				}

			} break;
			case GEOTRAIN:
			{
				message("Number of modalities: " + std::to_string(m_input_images.size() + m_input_images_MRI.size()) + "\n");
				m_processing.PreProcess(m_input_images, m_labels_image);
				auto labelsCountMap = ParserGTS::CountsOfEachLabel<LabelsImageType>(m_labels_image);

				if (!validLabels(gtsResult, labelsCountMap, 2)) { return gtsResult; }

				addToInputImagesMRI(); // Nothing special happens for the different modalities here

				addCoordinatorMapsIfNecessary();

				// SVM_PSEUDO
				svm<PixelType>(m_input_images, m_labels_image, gtsResult, labelsCountMap, true);
				if (!gtsResult->ok) { break; }

				m_processing.template PostProcessNormalImage<PseudoProbImageType>(gtsResult->posImage);
				m_processing.template PostProcessNormalImage<PseudoProbImageType>(gtsResult->negImage);

				// AGD
				AgdImagePointer agdMapPos = agd<PseudoProbPixelType>(gtsResult->posImage, gtsResult->posLabel, true);
				AgdImagePointer agdMapNeg = agd<PseudoProbPixelType>(gtsResult->negImage, gtsResult->negLabel, true);

				// THRESHOLD
				gtsResult->labelsImage = thresholdSingleImageToLabel<AgdImageType>(agdMapPos, m_threshold, gtsResult->posLabel);
				thresholdSingleImageToLabel<AgdImageType>(agdMapNeg, m_threshold, gtsResult->negLabel);

				if (m_ground_truth_set) {
					checkAccuracyInRelationToGroundTruth(gtsResult->labelsImage, gtsResult);
				}

			} break;
			case GEOTRAIN_FULL:
			{
				message("Number of modalities: " + std::to_string(m_input_images.size() + m_input_images_MRI.size()) + "\n");
				m_processing.PreProcess(m_input_images, m_labels_image);
				auto labelsCountMap = ParserGTS::CountsOfEachLabel<LabelsImageType>(m_labels_image);

				if (!validLabels(gtsResult, labelsCountMap, 2)) { return gtsResult; }

				// Construct agd maps for each input image and input input images and agd maps to a SVM that produces labels
				std::vector< AgdImagePointer > agdImages;

				if (m_input_images.size() != 0) {
					auto keys = getMapKeyset(labelsCountMap);

					std::vector< AgdImagePointer > agdImagesAllKeys = agd<PixelType>(m_input_images, keys, true);
					agdImages.insert(agdImages.end(), agdImagesAllKeys.begin(), agdImagesAllKeys.end()); // Append
				}
				
				if (m_input_images_MRI.size() != 0) {
					// Special operations for MRI images
					std::vector< AgdImagePointer > agdImagesMRI = agdMRI<PixelType>(m_input_images_MRI, getMapKeyset(labelsCountMap), false, true);
					agdImages.insert(agdImages.end(), agdImagesMRI.begin(), agdImagesMRI.end()); // Append
					addToInputImagesMRI();
				}

				for (auto agdImage : agdImages) {
					// Convert output of AGD to InputImageType
					auto agdImageCast = ItkUtilGTS::castAndRescaleImage< AgdImageType, InputImageType>(agdImage);

					// Apply negative filter so it is 0 outside of the brain (basically where there is no information - that is used to make it faster later)
					typedef itk::InvertIntensityImageFilter<InputImageType> InvertIntensityImageFilterType;
					typename InvertIntensityImageFilterType::Pointer invertIntensityFilter = InvertIntensityImageFilterType::New();
					invertIntensityFilter->SetInput(agdImageCast);
					invertIntensityFilter->SetMaximum(255);
					invertIntensityFilter->Update();

					m_input_images.push_back(invertIntensityFilter->GetOutput());
					m_agd_maps_count++;
				}

				addCoordinatorMapsIfNecessary();

				// SVM_PSEUDO
				svm<PixelType>(m_input_images, m_labels_image, gtsResult, labelsCountMap, true);
				if (!gtsResult->ok) { break; }

				// AGD
				AgdImagePointer agdMapPos = agd<PseudoProbPixelType>(gtsResult->posImage, gtsResult->posLabel, true);
				AgdImagePointer agdMapNeg = agd<PseudoProbPixelType>(gtsResult->negImage, gtsResult->negLabel, true);

				gtsResult->labelsImage = thresholdSingleImageToLabel<AgdImageType>(agdMapPos, m_threshold, gtsResult->posLabel);
				
				m_processing.PostProcessLabelsImage(gtsResult->labelsImage);

				// THRESHOLD
				thresholdSingleImageToLabel<AgdImageType>(agdMapNeg, m_threshold, gtsResult->negLabel);

				if (m_ground_truth_set) {
					checkAccuracyInRelationToGroundTruth(gtsResult->labelsImage, gtsResult);
				}

			} break;
			case RF:
			{
				message("Number of modalities: " + std::to_string(m_input_images.size() + m_input_images_MRI.size()) + "\n");
				m_processing.PreProcess(m_input_images, m_labels_image);
				auto labelsCountMap = ParserGTS::CountsOfEachLabel<LabelsImageType>(m_labels_image);

				if (!validLabels(gtsResult, labelsCountMap)) { return gtsResult; }

				addToInputImagesMRI(); // Nothing special happens for the different modalities here

				addCoordinatorMapsIfNecessary();

				gtsResult->labelsImage = rf(false);

				m_processing.PostProcessLabelsImage(gtsResult->labelsImage);

				// Change labels if necessary
				changeLabels(gtsResult->labelsImage, m_change_labels_map);

				if (m_ground_truth_set) {
					checkAccuracyInRelationToGroundTruth(gtsResult->labelsImage, gtsResult);
				}
			} break;
			case AGD_RF:
			{
				message("Number of modalities: " + std::to_string(m_input_images.size() + m_input_images_MRI.size()) + "\n");
				m_processing.PreProcess(m_input_images, m_labels_image);
				auto labelsCountMap = ParserGTS::CountsOfEachLabel<LabelsImageType>(m_labels_image);

				if (!validLabels(gtsResult, labelsCountMap)) { return gtsResult; }

				// Construct agd maps for each input image and input input images and agd maps to a SVM that produces labels
				std::vector< AgdImagePointer > agdImages; //= agd<PixelType>(m_input_images, m_labels_list, false);

				if (m_input_images.size() != 0) {
					auto keys = getMapKeyset(labelsCountMap);

					std::vector< AgdImagePointer > agdImagesAllKeys = agd<PixelType>(m_input_images, keys, true);
					agdImages.insert(agdImages.end(), agdImagesAllKeys.begin(), agdImagesAllKeys.end()); // Append
				}

				if (m_input_images_MRI.size() != 0) {
					// Special operations for MRI images
					std::vector< AgdImagePointer > agdImagesMRI = agdMRI<PixelType>(m_input_images_MRI, getMapKeyset(labelsCountMap), false, true);
					agdImages.insert(agdImages.end(), agdImagesMRI.begin(), agdImagesMRI.end()); // Append
					addToInputImagesMRI();
				}

				for (auto agdImage : agdImages) {
					// Convert output of AGD to InputImageType
					auto agdImageCast = ItkUtilGTS::castAndRescaleImage< AgdImageType, InputImageType>(agdImage);

					// Apply negative filter so it is 0 outside of the brain (basically where there is no information - that is used to make it faster later)
					typedef itk::InvertIntensityImageFilter<InputImageType> InvertIntensityImageFilterType;
					typename InvertIntensityImageFilterType::Pointer invertIntensityFilter = InvertIntensityImageFilterType::New();
					invertIntensityFilter->SetInput(agdImageCast);
					invertIntensityFilter->SetMaximum(255);
					invertIntensityFilter->Update();

					m_input_images.push_back(invertIntensityFilter->GetOutput());
					m_agd_maps_count++;
				}

				addCoordinatorMapsIfNecessary();

				// RF
				gtsResult->labelsImage = rf(false);

				m_processing.PostProcessLabelsImage(gtsResult->labelsImage);

				// Change labels if necessary
				changeLabels(gtsResult->labelsImage, m_change_labels_map);

				if (m_ground_truth_set) {
					checkAccuracyInRelationToGroundTruth(gtsResult->labelsImage, gtsResult);
				}
			} break;
			case RF_AUTO:
			{
				message("Number of modalities: " + std::to_string(m_input_images.size() + m_input_images_MRI.size()) + "\n");
				m_processing.PreProcess(m_input_images, m_labels_image);
				auto labelsCountMap = ParserGTS::CountsOfEachLabel<LabelsImageType>(m_labels_image);

				if (!validLabels(gtsResult, labelsCountMap)) { return gtsResult; }

				addToInputImagesMRI(); // Nothing special happens for the different modalities here

				addCoordinatorMapsIfNecessary();

				gtsResult->labelsImage = rf(true);

				m_processing.PostProcessLabelsImage(gtsResult->labelsImage);

				// Change labels if necessary
				changeLabels(gtsResult->labelsImage, m_change_labels_map);

				if (m_ground_truth_set) {
					checkAccuracyInRelationToGroundTruth(gtsResult->labelsImage, gtsResult);
				}
			} break;
			case AGD_RF_AUTO:
			{
				message("Number of modalities: " + std::to_string(m_input_images.size() + m_input_images_MRI.size()) + "\n");
				m_processing.PreProcess(m_input_images, m_labels_image);
				auto labelsCountMap = ParserGTS::CountsOfEachLabel<LabelsImageType>(m_labels_image);

				if (!validLabels(gtsResult, labelsCountMap)) { return gtsResult; }

				// Construct agd maps for each input image and input input images and agd maps to a SVM that produces labels
				std::vector< AgdImagePointer > agdImages;

				if (m_input_images.size() != 0) {
					auto keys = getMapKeyset(labelsCountMap);

					std::vector< AgdImagePointer > agdImagesAllKeys = agd<PixelType>(m_input_images, keys, true);
					agdImages.insert(agdImages.end(), agdImagesAllKeys.begin(), agdImagesAllKeys.end()); // Append
				}

				if (m_input_images_MRI.size() != 0) {
					// Special operations for MRI images
					std::vector< AgdImagePointer > agdImagesMRI = agdMRI<PixelType>(m_input_images_MRI, getMapKeyset(labelsCountMap), false, true);
					agdImages.insert(agdImages.end(), agdImagesMRI.begin(), agdImagesMRI.end()); // Append
					addToInputImagesMRI();
				}

				for (auto agdImage : agdImages) {
					// Convert output of AGD to InputImageType
					auto agdImageCast = ItkUtilGTS::castAndRescaleImage< AgdImageType, InputImageType>(agdImage);

					// Apply negative filter so it is 0 outside of the brain (basically where there is no information - that is used to make it faster later)
					typedef itk::InvertIntensityImageFilter<InputImageType> InvertIntensityImageFilterType;
					typename InvertIntensityImageFilterType::Pointer invertIntensityFilter = InvertIntensityImageFilterType::New();
					invertIntensityFilter->SetInput(agdImageCast);
					invertIntensityFilter->SetMaximum(255);
					invertIntensityFilter->Update();

					m_input_images.push_back(invertIntensityFilter->GetOutput());
					m_agd_maps_count++;
				}

				addCoordinatorMapsIfNecessary();

				gtsResult->labelsImage = rf(true);

				m_processing.PostProcessLabelsImage(gtsResult->labelsImage);

				// Change labels if necessary
				changeLabels(gtsResult->labelsImage, m_change_labels_map);

				if (m_ground_truth_set) {
					checkAccuracyInRelationToGroundTruth(gtsResult->labelsImage, gtsResult);
				}
			} break;
			case CHANGE_LABELS:
			{
				gtsResult->labelsImage = changeLabels(m_labels_image, m_change_labels_map);
			} break;
			case GENERATE_CONFIG:
			{
				// Generate .yaml config from pretrained models (.xml) using their respective neighborhood radius and importance value
				std::vector<int> neighborhoodRadii(m_importance_values.size(), 1); // Neighborhood radius should be removed from SvmSuite in the future
				SvmSuite::generateConfig(m_pretrained_models_paths, neighborhoodRadii, m_importance_values, m_output_folder + "/config.yaml");
			} break;
			case CHECK_ACCURACY:
			{
				// changeLabels won't do anything if there is no need
				changeLabels(m_labels_image, m_change_labels_map);

				checkAccuracyInRelationToGroundTruth(m_labels_image, gtsResult);
			} break;
			}

			stopTimerAndReport("Total");
			message("", "Finished");
			
			return gtsResult;
		}

		//Setters

		void SetMode(MODE mode) {
			m_mode = mode;
		}
		void SetInputImages(std::vector< InputImagePointer > inputImages) {
			m_input_images = inputImages;
		}
		void SetInputImages(std::vector< std::string > inputImagesPaths) {
			if (inputImagesPaths.size() != 0) {
				std::vector< InputImagePointer > inputImages;

				for (auto inputImagePath : inputImagesPaths) {
					std::string extension = getFileExtension(inputImagePath);

					if (!m_file_extension_set_manually) {
						m_file_extension = extension;
					}

					inputImages.push_back(readImage<InputImageType>(inputImagePath));
				}

				m_input_images = inputImages; // Replace possible older values
			}
		}
		void SetInputImage(typename itk::Image< float, 3 >::Pointer inputImage) {
			m_input_images = std::vector< InputImagePointer >();
			m_input_images.push_back(inputImage);
		}
		void SetInputImage(std::string inputImagePath) {
			m_input_images = std::vector< InputImagePointer >();
			m_input_images.push_back(readImage<InputImageType>(inputImagePath));
		}
		void SetLabels(LabelsImagePointer labels) {
			m_labels_image = labels;
		}
		void SetLabels(std::string labelsPath) {
			if (labelsPath != "") {
				std::string extension = getFileExtension(labelsPath);

				if (!m_file_extension_set_manually) {
					m_file_extension = extension;
				}

				m_labels_image = readImage<LabelsImageType>(labelsPath);
			}
		}
		void SetOutputPath(std::string path) {
			m_output_folder = path;
			if (!cbica::isDir(m_output_folder)) {
				cbica::createDir(m_output_folder);
			}
		}
		void SetOutputPath(std::string path, std::string datasetName, std::string tag = "", bool includeDateTimeInTag = true) {
			if (datasetName != "") {
				if (!cbica::isDir(path)) {
					cbica::createDir(path);
				}
				m_output_folder = path + "\\" + datasetName;
				if (!cbica::isDir(m_output_folder)) {
					cbica::createDir(m_output_folder);
				}
				std::string subfolderName =
					((includeDateTimeInTag || tag == "") ?
						UtilGTS::currentDateTime()
						: "")
					+ ((tag != "") ?
					((includeDateTimeInTag) ?
						" "
						: "") + tag
						: "");

				//std::string subfolderName = current_date() + " " + current_time();
				//std::string subfolderName = "output";

				m_output_folder += "\\" + subfolderName;
				if (!cbica::isDir(m_output_folder)) {
					cbica::createDir(m_output_folder);
				}
			}
			else {
				this->SetOutputPath(path);
			}
		}
		void SetDoCoordinateMaps(bool doCoordinateMaps = true)
		{
			m_do_coordinate_maps = doCoordinateMaps;
		}
		void SetProcessing(bool onOrOff, float imageToAgdMapsRatio = 6, 
			bool doStatFilter = true, bool doCurvatureFilter = false, 
			bool limitPixels = true, int pixelLimit = 5000000)
		{
			m_process = onOrOff;

			if (m_process)
			{			
				m_processing.SetLimitPixels(limitPixels, pixelLimit);
				m_processing.SetDoStatisticalNormalization(doStatFilter, imageToAgdMapsRatio);
				m_processing.SetDoCurvatureAnisotropic(doCurvatureFilter);
			}
		}
		void SetSaveAll(bool saveAll) {
			m_save_all = saveAll;
		}
		void SaveOnlyNormalSegmentation(bool saveOnlySeg = true, std::string segName = "labels_res")
		{
			m_save_only_seg = saveOnlySeg;
			m_save_only_seg_name = segName;
		}
		void SetTimerEnabled(bool timerEnabled) {
			m_timer_enabled = timerEnabled;
		}
		void SetConfigFile(std::string configFilePath) {
			if (configFilePath != "") {
				m_config_file_path = configFilePath;
			}
		}
		void SetRfConfigFile(std::string rfConfigFilePath) {
			if (rfConfigFilePath != "") {
				m_rf_config_file_path = rfConfigFilePath;
			}
		}
		void SetThreshold(float threshold) {
			m_threshold = threshold;
		}
		void SetChangeLabelsMap(std::unordered_map< LabelsPixelType, LabelsPixelType > changeLabelsMap) {
			m_change_labels_map = changeLabelsMap;
			m_changed_labels_map_manually_set = true;
		}
		void SetImportanceValues(std::vector< double > importanceValues)
		{
			m_importance_values = importanceValues;
		}
		void SetPretrainedModelsPaths(std::vector< std::string > modelsPaths) {
			m_pretrained_models_paths = modelsPaths;
		}
		void SetNumberOfThreads(int numberOfThreads) {
			if (numberOfThreads > 0) {
				m_number_of_threads = numberOfThreads;
			}
		}
		void SetNumberOfThreadsMax(bool maxThreads) {
			m_max_threads = maxThreads;
		}
		void SetGroundTruth(std::string groundTruthPath, std::vector<LabelsPixelType> groundTruthSkip = std::vector<LabelsPixelType>())
		{
			if (groundTruthPath != "") {
				std::string extension = getFileExtension(groundTruthPath);

				if (!m_file_extension_set_manually) {
					m_file_extension = extension;
				}

				m_ground_truth_set = true;
				m_ground_truth = readImage<LabelsImageType>(groundTruthPath);
				m_ground_truth_skip = groundTruthSkip;
			}
		}
		void SetGroundTruth(LabelsImagePointer groundTruth, std::vector<LabelsPixelType> groundTruthSkip = std::vector<LabelsPixelType>())
		{
			m_ground_truth_set = true;
			m_ground_truth = groundTruth;
			m_ground_truth_skip = groundTruthSkip;
		}
		void SetOutputImageFileExtension(std::string fileExtension) {
			m_file_extension = fileExtension;
			m_file_extension_set_manually = true;
		}
		void SetInputImageMRI(MODALITY_MRI modality, InputImagePointer imageMRI) {
			m_input_images_MRI[modality] = imageMRI;
		}
		void SetInputImageMRI(MODALITY_MRI modality, std::string imagePathMRI) {
			if (imagePathMRI != "") {
				std::string extension = getFileExtension(imagePathMRI);

				if (!m_file_extension_set_manually) {
					m_file_extension = extension;
				}

				m_input_images_MRI[modality] = readImage<InputImageType>(imagePathMRI);
			}
		}
		void SetTumorCoreLabelMRI(LabelsPixelType labelTC) {
			m_label_TC = labelTC;
		}
		void SetEnhancedTumorLabelMRI(LabelsPixelType labelET) {
			m_label_ET = labelET;
		}
		void SetEdemaLabelMRI(LabelsPixelType labelED) {
			m_label_ED = labelED;
		}
		void SetHealthyTissueLabelMRI(LabelsPixelType labelHT) {
			m_label_HT = labelHT;
		}
		void SetLabelOfInterest(LabelsPixelType labelOfInterest) {
			m_label_of_interest = labelOfInterest;
		}
		void SetVerbose(bool verbose) {
			m_verbose = verbose;
		}
		void SetSubsampling(bool subsample, int maxSamples = 3000) {
			m_subsample = subsample;
			m_max_samples_svm_subsample = maxSamples;
		}
		void SetBalancedSubsampling(bool balancedSubsample) {
			m_balanced_subsample = balancedSubsample;
		}
		// void SetInputImageToAgdMapsRatio(float ratio) {
		// 	m_image_to_agd_maps_ratio = ratio;
		// }

	private:
		LabelsImagePointer                                     m_labels_image;
		std::vector< InputImagePointer >                       m_input_images;
		InputImagePointer                                      m_reference_image;
		GeodesicTrainingSegmentation::MODE                     m_mode = GeodesicTrainingSegmentation::MODE::REVERSE_GEOTRAIN;
		std::vector< double >                                  m_importance_values;       // For mode:  generateconfig
		std::unordered_map< LabelsPixelType, LabelsPixelType > m_change_labels_map;
		std::vector< std::string >                             m_pretrained_models_paths; // For mode:  generateconfig
		std::mutex                                             m_agd_mutex;
		LabelsImagePointer                                     m_ground_truth;
		std::vector<LabelsPixelType>                           m_ground_truth_skip = std::vector<LabelsPixelType>();
		std::map< MODALITY_MRI, InputImagePointer >            m_input_images_MRI;        // For MRI
		SvmSuiteUtil::Timer                                    m_timer;                   // For timer
		int                                                    m_number_of_threads = 16, m_agd_maps_count = 0,
                                                               m_max_samples_svm_subsample = DEFAULT_MAX_SAMPLES_SVM_SUBSAMPLE;
		float                                                  m_threshold = DEFAULT_THRESHOLD,
                                                               m_image_to_agd_maps_ratio = DEFAULT_INPUT_IMAGES_TO_AGD_MAPS_RATIO;
		bool            m_save_all = false, m_timer_enabled = false, m_subsample = true, m_balanced_subsample = DEFAULT_BALANCED_SUBSAMPLE,
                        m_file_extension_set_manually = false, m_verbose = false, m_were_images_shrunk = false, m_save_only_seg = false,
                        m_ground_truth_set = false, m_max_threads = false, m_changed_labels_map_manually_set = false, 
                        m_do_coordinate_maps = true, m_process = true;
		std::string     m_config_file_path = "", m_rf_config_file_path = "", m_save_only_seg_name = "",
                        m_file_extension = DEFAULT_FILE_EXTENSION, m_output_folder = cbica::getExecutablePath();
		LabelsPixelType m_label_TC = DEFAULT_LABEL_TC, m_label_ET = DEFAULT_LABEL_ET,
                        m_label_ED = DEFAULT_LABEL_ED, m_label_HT = DEFAULT_LABEL_HT,
                        m_label_of_interest = DEFAULT_LABEL_OF_INTEREST;
        GeodesicTrainingSegmentation::Processing<PixelType, Dimensions> m_processing;

		// For SVM

		template<typename TPixelType>
		void svm(std::vector < typename itk::Image<TPixelType, Dimensions>::Pointer > images, LabelsImagePointer labels, 
                 const std::shared_ptr<Result>& gtsResult, std::unordered_map<LabelsPixelType, int> labelsCountMap, bool predictFlags)
		{
			SvmManagerGTS<TPixelType, Dimensions> svmManagerGTS;

			if (m_config_file_path != "") {
				svmManagerGTS.AddSvmsFromConfig(m_config_file_path);
			}
			else {
				svmManagerGTS.AddSvmDescriptions(SvmSuite::SvmDescription::GetDefaultSvmDescriptions());
			}

			svmManagerGTS.SetOutputPath(m_output_folder);
			svmManagerGTS.SetTimerEnabled(m_timer_enabled);
			svmManagerGTS.SetSavingModelsEnabled(m_save_all);
			svmManagerGTS.SetVerbose(true);
			svmManagerGTS.SetSavingModelsEnabled(m_save_all);
			svmManagerGTS.SetNumberOfThreads((m_max_threads) ? MAX_THREADS_SVM_SUITE : m_number_of_threads);
			//svmManagerGTS.SetInputNormalization(true); // Old way

			std::cout << "Converting input images to matrices...";
			//auto data = ParserGTS::NormalizedParse<TPixelType, Dimensions>(images, labels, false); // New way
			//ParserGTS::ScaleSomeOfTheColumns(data->trainingMat,
			//	data->trainingMat.cols - m_agd_maps_count,
			//	data->trainingMat.cols - 1, 
			//	1 / m_image_to_agd_maps_ratio
			//);
			auto data = ParserGTS::Parse<TPixelType, Dimensions>(images, labels, false);
			std::cout << "finished\n";
			
			//cv::Mat sampleIdx;

			if (!m_balanced_subsample) {
				svmManagerGTS.SetSubsampling(m_subsample, m_max_samples_svm_subsample);
			}
			else {
				if (m_subsample && data->trainingMat.rows > m_max_samples_svm_subsample) {
					// /*sampleIdx = */CreateBalancedSubsample
					std::string errorMessageIfApplicable;
					if (!CreateBalancedSubsample(data, errorMessageIfApplicable, labelsCountMap, 
						m_max_samples_svm_subsample))
					{
						gtsResult->ok = false;
						gtsResult->errorMessage = errorMessageIfApplicable; 
						return;
					}
				}
			}

			message("", "SVM Operations");
			auto result = svmManagerGTS.TrainAndTestGTS(data, m_input_images, m_labels_image, predictFlags/*, sampleIdx*/);
			message("", "SVM Operations", 100);

			svmManagerGTS.GenerateConfigFromBestValues();

			if (predictFlags) {
				// Swap so that label of interest is always the pos
				if (result->negLabel == m_label_of_interest) {
					PseudoProbImagePointer tempImage = result->negImage;
					LabelsPixelType tempLabel = result->negLabel;
					result->negImage = result->posImage;
					result->negLabel = result->posLabel;
					result->posImage = tempImage;
					result->posLabel = tempLabel;
				}

				writeImage< PseudoProbImageType >(result->posImage, "pseudo_class" + std::to_string(result->posLabel));
				writeImage< PseudoProbImageType >(result->negImage, "pseudo_class" + std::to_string(result->negLabel));
			}
			else {
				writeImage< LabelsImageType >(result->labelsImage, "labels_res");
			}

			gtsResult->posImage = result->posImage;
			gtsResult->negImage = result->negImage;
			gtsResult->posLabel = result->posLabel;
			gtsResult->negLabel = result->negLabel;
			gtsResult->labelsImage = result->labelsImage;
		}

		// For AGD

		template<typename TPixelType>
		AgdImagePointer agd(typename itk::Image< TPixelType, Dimensions >::Pointer inputImage,
			LabelsPixelType agdLabel, bool saveResults = false)
		{
			std::vector < typename itk::Image< TPixelType, Dimensions >::Pointer > tempInputVector;
			tempInputVector.push_back(inputImage);

			return agd<TPixelType>(tempInputVector, agdLabel, saveResults)[0];
		}

		template<typename TPixelType>
		std::vector< AgdImagePointer > agd(std::vector< typename itk::Image< TPixelType, Dimensions >::Pointer > inputImages,
			std::vector< LabelsPixelType> agdLabels, bool saveResults = false)
		{
			message("AGD...", "AGD Operations");
			typedef itk::Image<TPixelType, Dimensions> AgdOriginalImageType;
			typedef itk::Image<AgdPixelType, 3> AgdImage3DType; // The algorithm implementation needs this type
			typedef AgdImage3DType::Pointer AgdImage3DPointer;

			std::vector< AgdImagePointer > agdResults;

			//int counterForFileName = 1;
			int counterForThreadsVec = 0;
			std::vector<std::thread> threads(inputImages.size() * agdLabels.size());
			int numberOfOpenThreads = 0;
			int oldestOpenThread = 0;

			for (int i = 0; i < inputImages.size(); i++)
			{
				for (int j = 0; j < agdLabels.size(); j++)
				{
					if (!m_max_threads && numberOfOpenThreads == m_number_of_threads) {
						threads[oldestOpenThread].join();
						oldestOpenThread++;
						numberOfOpenThreads--;
					}
					//std::cout << "Starting new thread for agd.\n";
					numberOfOpenThreads++;
					threads[counterForThreadsVec++] = std::thread(&Coordinator<PixelType, Dimensions>::agdThreadJob<TPixelType>, this,
						std::ref(agdResults), std::to_string(i + 1), inputImages[i], agdLabels[j], saveResults
					);
				}
			}

			for (int i = oldestOpenThread; i < inputImages.size() * agdLabels.size(); i++) {
				threads[i].join();
			}

			message("AGD...", "AGD Operations", 100);
			return agdResults;
		}

		template<typename TPixelType>
		std::vector< AgdImagePointer > agd(std::vector< typename itk::Image< TPixelType, Dimensions >::Pointer > inputImages,
			LabelsPixelType agdLabel, bool saveResults = false)
		{
			message("AGD...", "AGD Operations");
			typedef itk::Image<TPixelType, Dimensions> AgdOriginalImageType;
			typedef itk::Image<AgdPixelType, 3> AgdImage3DType; // The algorithm implementation needs this type
			typedef AgdImage3DType::Pointer AgdImage3DPointer;

			std::vector< AgdImagePointer > agdResults;

			//int counterForFileName = 1;
			int counterForThreadsVec = 0;
			std::vector<std::thread> threads(inputImages.size());
			int numberOfOpenThreads = 0;
			int oldestOpenThread = 0;

			for (int i = 0; i < inputImages.size(); i++)
			{
				if (!m_max_threads && numberOfOpenThreads == m_number_of_threads) {
					threads[oldestOpenThread].join();
					oldestOpenThread++;
					numberOfOpenThreads--;
				}
				//std::cout << "Starting new thread for agd.\n";
				numberOfOpenThreads++;
				threads[counterForThreadsVec++] = std::thread(&Coordinator<PixelType, Dimensions>::agdThreadJob<TPixelType>, this,
					std::ref(agdResults), std::to_string(i + 1), inputImages[i], agdLabel, saveResults
				);
			}

			for (int i = oldestOpenThread; i < inputImages.size(); i++) {
				threads[i].join();
			}

			message("AGD...", "AGD Operations", 100);
			return agdResults;
		}

		template<typename TPixelType>
		std::vector< AgdImagePointer > agdMRI(std::map< MODALITY_MRI, typename itk::Image< TPixelType, Dimensions >::Pointer > &inputImagesMRI,
			std::vector<LabelsPixelType> allLabels, bool saveResults = true, bool asSvmInput = false)
		{
			message("AGD (MRI)...", "AGD Operations");
			typedef itk::Image<TPixelType, Dimensions> AgdOriginalImageType;
			typedef itk::Image<AgdPixelType, Dimensions> AgdImageType; // The algorithm implementation needs this type
			typedef typename AgdImageType::Pointer AgdImagePointer;
			
			std::vector< AgdImagePointer > agdResults;

			std::map< MODALITY_MRI, std::vector<LabelsPixelType> > epicenterForEachModality;

			bool atLeastOneEpicenterSet = false;

			if (asSvmInput) {
				epicenterForEachModality[FLAIR] = {  };
				epicenterForEachModality[T1]    = {  };
				epicenterForEachModality[T1CE]  = {  };
				epicenterForEachModality[T2]    = {  };

				if ( std::find(allLabels.begin(), allLabels.end(), m_label_TC) != allLabels.end() ) {
					epicenterForEachModality[FLAIR].push_back(m_label_TC);
					epicenterForEachModality[T1].push_back(m_label_TC);
					epicenterForEachModality[T1CE].push_back(m_label_TC);
					epicenterForEachModality[T2].push_back(m_label_TC);
					atLeastOneEpicenterSet = true;
				}
				if ( std::find(allLabels.begin(), allLabels.end(), m_label_ET) != allLabels.end() ) {
					epicenterForEachModality[T1CE].push_back(m_label_ET);
					atLeastOneEpicenterSet = true;
				}
				if ( std::find(allLabels.begin(), allLabels.end(), m_label_ED) != allLabels.end() ) {
					epicenterForEachModality[FLAIR].push_back(m_label_ED);
					atLeastOneEpicenterSet = true;
				}
				if ( std::find(allLabels.begin(), allLabels.end(), m_label_HT) != allLabels.end() ) {
					epicenterForEachModality[FLAIR].push_back(m_label_HT);
					atLeastOneEpicenterSet = true;
				}

				for (LabelsPixelType label : allLabels) {
					if (label != m_label_TC && label != m_label_ET && label != m_label_ED && label != m_label_HT) {
						epicenterForEachModality[FLAIR].push_back(label);
						epicenterForEachModality[T1].push_back(label);
						epicenterForEachModality[T1CE].push_back(label);
						epicenterForEachModality[T2].push_back(label);
						atLeastOneEpicenterSet = true;
					}
				}
				
				/*epicenterForEachModality[FLAIR] = { m_label_TC, m_label_ED, m_label_HT };
				epicenterForEachModality[T1] = { m_label_TC };
				epicenterForEachModality[T1CE] = { m_label_TC, m_label_ET };
				epicenterForEachModality[T2] = { m_label_TC };*/

				//epicenterForEachModality[FLAIR] = { 1, 2, 3/*, 4*/ }; // class 2 is the edema 
				//epicenterForEachModality[T1]    = { 1/*, 2, 3, 4*/ }; // class 3 is the healthy tissue
				//epicenterForEachModality[T1CE]  = { 1, /*2, 3,*/ 4 }; // class 4 is the tumor border
				//epicenterForEachModality[T2]    = { 1/*, 2, 3, 4*/ }; // class 2 is the edema
			}
			else {
				// RF input
				if (m_label_TC != 0) {
					epicenterForEachModality[FLAIR] = { m_label_TC };
					epicenterForEachModality[T1]    = {  };
					epicenterForEachModality[T1CE]  = { m_label_TC };
					epicenterForEachModality[T2]    = {  };
					atLeastOneEpicenterSet = true;
				}

				for (LabelsPixelType label : allLabels) {
					if (label != m_label_TC && label != m_label_ET && label != m_label_ED && label != m_label_HT) {
						epicenterForEachModality[FLAIR].push_back(label);
						epicenterForEachModality[T1].push_back(label);
						epicenterForEachModality[T1CE].push_back(label);
						epicenterForEachModality[T2].push_back(label);
						atLeastOneEpicenterSet = true;
					}
				}

				//epicenterForEachModality[FLAIR] = { 1/*, 2, 3, 4*/ }; // class 2 is the edema 
				//epicenterForEachModality[T1]    = { /*1, 2, 3, 4*/ }; // class 3 is the healthy tissue
				//epicenterForEachModality[T1CE]  = { 1, /*2, 3, 4*/ }; // class 4 is the tumor border
				//epicenterForEachModality[T2]    = { /*1, 2, 3, 4*/ }; // class 2 is the edema
			}

			if (!atLeastOneEpicenterSet) {
				// Not really applicable anymore, but is kept in case anything changes in the future
				if (asSvmInput) {
					message("Can't perform AGD because there are no labels for TC, ET, ED or HT. Will continue without.\n", "");
				}
				else {
					message("Can't perform AGD because there are no labels for TC. Will continue without.\n", "");
				}
				return agdResults;
			}
			
			auto keys = getMapKeyset(inputImagesMRI);
			int threadsNumber = 0;

			for (auto key : keys) {
				threadsNumber += epicenterForEachModality[key].size();
			}

			//int counterForFileName = 1;
			int counterForThreadsVec = 0;
			std::vector<std::thread> threads(threadsNumber);
			int numberOfOpenThreads = 0;
			int oldestOpenThread = 0;

			for (auto key : keys)
			{
				for (LabelsPixelType epicenter : epicenterForEachModality[key])
				{
					if (!m_max_threads && numberOfOpenThreads == m_number_of_threads) {
						threads[oldestOpenThread].join();
						oldestOpenThread++;
						numberOfOpenThreads--;
					}
					//std::cout << "Starting new thread for agd MRI.\n";
					numberOfOpenThreads++;
					threads[counterForThreadsVec++] = std::thread(&Coordinator<PixelType, Dimensions>::agdThreadJob<TPixelType>, this,
						std::ref(agdResults), getModalityName(key), inputImagesMRI[key], epicenter, saveResults
					);
				}
			}

			for (int i = oldestOpenThread; i < threadsNumber; i++) {
				threads[i].join();
			}

			message("AGD (MRI)...", "AGD Operations", 100);
			return agdResults;
		}

		template<typename TPixelType>
		void agdThreadJob(std::vector< AgdImagePointer > &outputVec, std::string imageName,
			typename itk::Image< TPixelType, Dimensions >::Pointer inputImage, LabelsPixelType agdLabel, bool saveResults = true)
		{
			if (agdLabel == 0) {
				return;
			}

			typedef itk::Image<TPixelType, Dimensions> AgdOriginalImageType;
			typedef itk::Image<AgdPixelType, Dimensions> AgdImageType; // The algorithm implementation needs this type
			typedef typename AgdImageType::Pointer AgdImagePointer;

			// Normalize input to [0,255] (pseudoprob maps will be in [0,1], converting them to int would have been bad
			auto agdNormInputImage = ItkUtilGTS::normalizeImage<AgdOriginalImageType>(inputImage);

			/*std::cout << "Saving agd norm input image...";
			writeImage<AgdOriginalImageType>(agdNormInputImage, "agd" + std::to_string(counterForFileName)
				+ "_l" + std::to_string(agdLabels[i]) + "_normoriginput.nii.gz"
				);
			std::cout << "finished\n";*/

			// Convert input to AgdImageType
			//AgdImage3DPointer agdInput;
			//if (Dimensions < 3) {
			//	auto temp = castImage< AgdOriginalImageType, AgdImageType >(agdNormInputImage);
			//	
			//	typedef itk::JoinImageFilter<AgdImageType, AgdImageType3DType> JoinImageFilterType;
			//	
			//	JoinImageFilterType::Pointer joinFilter = JoinImageFilterType::New();
			//	joinFilter->SetInput(agdNormInputImage, agdInput);
			//	agdInput = joinFilter->GetOutput();
			//}
			//else {
			AgdImagePointer agdInput = ItkUtilGTS::castAndRescaleImage< AgdOriginalImageType, AgdImageType >(agdNormInputImage);
			//}

			/*std::cout << "Saving agd input image...";
			writeImage<AgdImageType>(agdInput, "agd" + std::to_string(counterForFileName)
				+ "_l" + std::to_string(agdLabels[i]) + "_input.nii.gz"
				);
			std::cout << "finished\n";*/

			// Get raw result of AGD algorithm (not in [0,255])
			AgdImagePointer agdOutRaw = AdaptiveGeodesicDistance::Run<AgdPixelType, Dimensions>(
				agdInput, agdInput, m_labels_image, agdLabel, false, true
			);

			// Normalize results to [0,255]
			//AgdImage3DPointer agdOutNorm = normalizeImage<AgdImage3DType>(agdOutRaw);

			// Convert output to correct dimensions
			//AgdImagePointer agdOut = castImage< AgdImage3DType, AgdImageType >(agdOutNorm);
			//AgdImagePointer agdOut = agdOutNorm;
			AgdImagePointer agdOut = ItkUtilGTS::normalizeImage<AgdImageType>(agdOutRaw);

			if (saveResults) {
				// Write results to file
				//std::cout << "Saving agd image...";
				writeImage<AgdImageType>(agdOut, "agd_" + imageName + "_class" + std::to_string(agdLabel));
				//std::cout << "finished\n";
			}

			std::lock_guard<std::mutex> lg(m_agd_mutex);
			outputVec.push_back(agdOut);
		}

		// For LABELSTHRES

		template<class UImageType>
		LabelsImagePointer thresholdSingleImageToLabel(typename UImageType::Pointer image, double threshold, LabelsPixelType labelOfInterest)
		{
			message("Thresholding...", "Thresholding");
			LabelsImagePointer res = ItkUtilGTS::initializeOutputImageBasedOn<UImageType, LabelsImageType>(image);

			itk::ImageRegionIterator<UImageType> iter_i(image, image->GetRequestedRegion());
			itk::ImageRegionIterator<LabelsImageType> iter_r(res, res->GetRequestedRegion());

			for (iter_i.GoToBegin(), iter_r.GoToBegin(); !iter_i.IsAtEnd(); ++iter_i, ++iter_r)
			{
				if (iter_i.Get() <= threshold) {
					iter_r.Set(labelOfInterest);
				}
				else {
					iter_r.Set(0);
				}
			}

			message("Thresholding...", "Thresholding", 100);
			writeImage<LabelsImageType>(res, "labelsthres_class" + std::to_string(labelOfInterest) + "_t" + std::to_string(m_threshold));

			return res;
		}

		// For SEGMENT

		template<class UImageType>
		typename UImageType::Pointer 
		segmentSingleImage(typename UImageType::Pointer image, LabelsImagePointer labels, LabelsPixelType labelOfInterest)
		{
			typename UImageType::Pointer res = ItkUtilGTS::initializeOutputImageBasedOn<UImageType, UImageType>(image);

			itk::ImageRegionIterator<UImageType> iter_i(image, image->GetRequestedRegion());
			itk::ImageRegionIterator<LabelsImageType> iter_l(labels, labels->GetRequestedRegion());
			itk::ImageRegionIterator<UImageType> iter_r(res, res->GetRequestedRegion());

			for (iter_i.GoToBegin(), iter_l.GoToBegin(), iter_r.GoToBegin(); !iter_i.IsAtEnd(); ++iter_i, ++iter_l, ++iter_r)
			{
				if (iter_l.Get() == labelOfInterest)
				{
					iter_r.Set(iter_i.Get());
				}
				else {
					iter_r.Set(0);
				}
			}
			writeImage<UImageType>(res, "segm");

			return res;
		}

		// For CHANGE_LABELS

		LabelsImagePointer changeLabels(LabelsImagePointer labels, std::unordered_map< int, int > &changeLabelsMap)
		{
			// isZeroIncludedInMap is for making it faster (usually 0 is not included) and map::find is time consuming
			bool isZeroIncludedInMap = false;

			if (changeLabelsMap.find(0) != changeLabelsMap.end()) {
				isZeroIncludedInMap = true;
			}

			if (changeLabelsMap.size() != 0) {
				message("Renaming labels...", "Changing labels");

				int howManyChanged = 0;

				itk::ImageRegionIteratorWithIndex< LabelsImageType > iter_l(labels, labels->GetRequestedRegion());

				for (iter_l.GoToBegin(); !iter_l.IsAtEnd(); ++iter_l)
				{
					LabelsPixelType prev_label = iter_l.Get();

					if ((prev_label == 0) && (!isZeroIncludedInMap)) {
						continue;
					}

					if (changeLabelsMap.find(prev_label) != changeLabelsMap.end()) {
						iter_l.Set(changeLabelsMap[prev_label]);
						howManyChanged++;
					}
				}

				writeImage<LabelsImageType>(labels, "labels_res_renamed");
				message("Renaming labels...", "Changing labels", 100);
				message("\t(Changed " + std::to_string(howManyChanged) + " labels)\n", "");
			}
			return labels;
		}

		// For RF

		LabelsImagePointer rf(bool trainAutoFlag)
		{
			message("Converting input images to matrices...", "Converting");
			auto data = ParserGTS::Parse<PixelType, Dimensions>(m_input_images, m_labels_image, false);
			message("Converting input images to matrices...", "Converting", 100);
			
			RFSuite::Manager rfManager;
			rfManager.SetTrainDataFromMats(data->trainingMat, data->labelsMat);
			rfManager.SetOutputPath(m_output_folder);
			rfManager.SetSaveAll(m_save_all);
			rfManager.SetVerbose(true);

			if (m_rf_config_file_path != "") {
				rfManager.SetParametersFromConfig(m_rf_config_file_path);
			}
			
			message("", "Training");
			if (!trainAutoFlag) {
				rfManager.Train();
			}
			else {
				rfManager.TrainAuto();
			}
			message("", "Training", 100);
			
			auto resMatPtr = rfManager.Test(data->testingMat, data->testingMat);

			message("Converting output matrix to image...", "Converting");
			LabelsImagePointer resImage = ItkUtilGTS::initializeOutputImageBasedOn<LabelsImageType>(m_labels_image);
			CvMatToImageGTS::FillImage<LabelsPixelType, Dimensions>(resImage, *resMatPtr);
			writeImage<LabelsImageType>(resImage, "res_rf");
			message("Converting output matrix to image...", "Converting", 100);
			
			return resImage;
		}

		// FOR CHECK_ACCURACY

		void checkAccuracyInRelationToGroundTruth(LabelsImagePointer labels, const std::shared_ptr<Result>& gtsResult)
		{
			// Convert vector to set for faster lookup
			
			std::set<LabelsPixelType> groundTruthSkip;
			groundTruthSkip.insert(0);

			for (LabelsPixelType skip : m_ground_truth_skip) {
				if (skip != 0) {
					message("Will skip label " + std::to_string(skip) + " at accuracy check.\n", "");
				}
				groundTruthSkip.insert(skip);
			}
			
			auto diffImage = ItkUtilGTS::initializeOutputImageBasedOn<LabelsImageType>(labels);

			message("Checking accuracy...", "Checking accuracy");

			itk::ImageRegionIteratorWithIndex< LabelsImageType > iter_l(labels, labels->GetRequestedRegion());
			itk::ImageRegionIteratorWithIndex< LabelsImageType > iter_d(diffImage, diffImage->GetRequestedRegion());
			itk::ImageRegionIteratorWithIndex< LabelsImageType > iter_g(m_ground_truth, m_ground_truth->GetRequestedRegion());

			LabelsPixelType lb;
			LabelsPixelType gt;

			std::map< LabelsPixelType, int > agreed, p1, t1, falsePositivesCount;

			bool foundValueLabels, foundValueGroundTruth;

			for (iter_l.GoToBegin(), iter_d.GoToBegin(), iter_g.GoToBegin(); !iter_l.IsAtEnd(); ++iter_l, ++iter_d, ++iter_g)
			{
				lb = iter_l.Get();
				gt = iter_g.Get();

				foundValueLabels      = (groundTruthSkip.find(lb) == groundTruthSkip.end());
				foundValueGroundTruth = (groundTruthSkip.find(gt) == groundTruthSkip.end());

				// Increment counters if necessary
				if (foundValueLabels) {
					incrementMapForCheckAccuracy(p1, lb);
				}
				if (foundValueGroundTruth) {
					incrementMapForCheckAccuracy(t1, gt);
				}
				if (foundValueLabels && foundValueGroundTruth && (lb == gt)) {
					incrementMapForCheckAccuracy(agreed, gt);
				}
				if (foundValueLabels && !foundValueGroundTruth) {
					incrementMapForCheckAccuracy(falsePositivesCount, lb);
				}
				/*if (foundValueLabels || foundValueGroundTruth) {
					correctSum += (lb == gt);
					allSum++;
				}*/

				// For saving the differences
				if (m_save_all) {
					/*if (foundValueLabels && !foundValueGroundTruth) {
						// False positive
						iter_l.Set(7); // Nothing particularly special about 7, it's just a high value
					}
					else*/ if (!foundValueLabels && !foundValueGroundTruth) {
						// Irrelevant voxel
						iter_d.Set(*groundTruthSkip.begin());
					}
					else if (foundValueLabels && foundValueGroundTruth && (lb == gt)) {
						// Voxel that the two segmentations agreed upon
						iter_d.Set(*groundTruthSkip.begin());
					}
					else {
						// Voxels where there was an actual difference
						iter_d.Set(5+gt); // Sets what it should have been + 5 (to differentiate false positives) 
					}
				}
			}

			message("Checking accuracy...", "Checking accuracy", 100);

			// Save the differences
			writeImage<LabelsImageType>(diffImage, "labels_diff_to_gt");

			// Make sure that each label exists everywhere

			for (auto p1k : getMapKeyset(p1))
			{
				if (t1.find(p1k) == t1.end()) {
					t1[p1k] = 0;
				}
				if (agreed.find(p1k) == agreed.end()) {
					agreed[p1k] = 0;
				}
				if (falsePositivesCount.find(p1k) == falsePositivesCount.end()) {
					falsePositivesCount[p1k] = 0;
				}
			}

			for (auto t1k : getMapKeyset(t1))
			{
				if (p1.find(t1k) == p1.end()) {
					p1[t1k] = 0;
				}
				if (agreed.find(t1k) == agreed.end()) {
					agreed[t1k] = 0;
				}
				if (falsePositivesCount.find(t1k) == falsePositivesCount.end()) {
					falsePositivesCount[t1k] = 0;
				}
			}

			// Get the different labels
			auto keys = getMapKeyset(p1);

			// Evaluation metrics for each label
			std::map< LabelsPixelType, double > diceScore, sensitivity;

			for (auto key : keys)
			{
				int valP1       = p1[key];
				int valT1       = t1[key];
				int valAgreed   = agreed[key];

				if (valP1 + valT1 != 0) {
					diceScore[key] = (static_cast<double>(valAgreed) * 2 * 100) / (valP1 + valT1);
				}
				else {
					std::cerr << "Error: checkAccuracyInRelationToGroundTruth: DSC DIV0\n";
				}

				if (valT1 != 0) {
					sensitivity[key] = static_cast<double>(valAgreed) * 100 / (valT1);
				}
				else {
					std::cerr << "Error: checkAccuracyInRelationToGroundTruth: Sens DIV0\n";
				}


			}

			// Raw results

			message("\n", "");
			
			for (auto key : keys) {
				if (m_verbose) {
					std::cout << "|P1| [LABEL " << key << "] (Output labels size):    " << p1[key] << "\n";
					std::cout << "|T1| [LABEL " << key << "] (Ground truth size):     " << t1[key] << "\n";
					std::cout << "|P1| /\\ |T1| [LABEL " << key << "] (Agreed upon):   " << agreed[key] << "\n\n";
					//std::cout << "Number of false positives [LABEL " << key << "]:    " << falsePositivesCount[key] << "\n\n";
				}

				if (m_save_all) {
					std::ofstream accuracyFile;
					accuracyFile.open(m_output_folder + "/accuracy_report.txt", std::ios_base::app); //append file

					accuracyFile << "|P1| [LABEL " << key << "] (Output labels size):    " << p1[key] << "\n";
					accuracyFile << "|T1| [LABEL " << key << "] (Ground truth size):     " << t1[key] << "\n";
					accuracyFile << "|P1| /\\ |T1| [LABEL " << key << "] (Agreed upon):    " << agreed[key] << "\n";
					accuracyFile << "Number of false positives [LABEL " << key << "]:    " << falsePositivesCount[key] << "\n\n";

					std::ofstream accuracyFileCompact;
					accuracyFileCompact.open(m_output_folder + "/accuracy_report_compact.txt", std::ios_base::app); //append file

					accuracyFileCompact << key << " " << p1[key] << "\n";
					accuracyFileCompact << key << " " << t1[key] << "\n";
					accuracyFileCompact << key << " " << agreed[key] << "\n";
					accuracyFileCompact << key << " " << falsePositivesCount[key] << "\n";
				}
			}

			// Evaluation metrics for all the labels combined
			
			double diceScoreAll = 0, sensitivityAll = 0;
			int p1All = 0, t1All = 0, agreedAll = 0, falsePositivesAll = 0;

			for (auto key : keys) {
				p1All += p1[key];
				t1All += t1[key];
				agreedAll += agreed[key];
				falsePositivesAll += falsePositivesCount[key];
			}

			diceScoreAll = (static_cast<double>(agreedAll) * 2 * 100) / (p1All + t1All);
			sensitivityAll = static_cast<double>(agreedAll) * 100 / (t1All);

			// Raw results for all

			if (m_verbose) {
				std::cout << "|P1| [ALL] (Output labels size):        " << p1All << "\n";
				std::cout << "|T1| [ALL] (Ground truth size):         " << t1All << "\n";
				std::cout << "|P1| /\\ |T1| [ALL] (Agreed upon):       " << agreedAll << "\n";
				std::cout << "Number of false positives [ALL]:        " << falsePositivesAll << "\n\n";
			}

			if (m_save_all) {
				std::ofstream accuracyFile;
				accuracyFile.open(m_output_folder + "/accuracy_report.txt", std::ios_base::app); //append file

				accuracyFile << "|P1| [ALL] (Output labels size):    " << p1All << "\n";
				accuracyFile << "|T1| [ALL] (Ground truth size):     " << t1All << "\n";
				accuracyFile << "|P1| /\\ |T1| [ALL] (Agreed upon):   " << agreedAll << "\n";
				accuracyFile << "Number of false positives [ALL]:    " << falsePositivesAll << "\n\n";

				std::ofstream accuracyFileCompact;
				accuracyFileCompact.open(m_output_folder + "/accuracy_report_compact.txt", std::ios_base::app); //append file

				accuracyFileCompact << 0 << " " << p1All << "\n";
				accuracyFileCompact << 0 << " " << t1All << "\n";
				accuracyFileCompact << 0 << " " << agreedAll << "\n";
				accuracyFileCompact << 0 << " " << falsePositivesAll << "\n";
			}

			// DSC For each label

			for (auto key : keys) {
				if (m_verbose) {
					std::cout << "DSC [" << key << "]   (%): " << diceScore[key] << "\n";
				}

				if (m_save_all) {
					std::ofstream accuracyFile;
					accuracyFile.open(m_output_folder + "/accuracy_report.txt", std::ios_base::app); //append file

					accuracyFile << "DSC [" << key << "]   (%): " << diceScore[key] << "\n";

					std::ofstream accuracyFileCompact;
					accuracyFileCompact.open(m_output_folder + "/accuracy_report_compact.txt", std::ios_base::app); //append file
				
					accuracyFileCompact << key << " " << diceScore[key] << "\n";
				}
			}

			// DSC For all

			if (m_verbose) {
				std::cout << "DSC [ALL] (%): " << diceScoreAll << "\n";
			}

			if (m_save_all) {
				std::ofstream accuracyFile;
				accuracyFile.open(m_output_folder + "/accuracy_report.txt", std::ios_base::app); //append file

				accuracyFile << "DSC [ALL] (%): " << diceScoreAll << "\n";

				std::ofstream accuracyFileCompact;
				accuracyFileCompact.open(m_output_folder + "/accuracy_report_compact.txt", std::ios_base::app); //append file

				accuracyFileCompact << 0 << " " << diceScoreAll << "\n";
			}


			// Sensitivity For each label
			
			if (m_verbose) {
				std::cout << "\n";
			}

			for (auto key : keys) {
				if (m_verbose) {
					std::cout << "Sens [" << key << "]   (%): " << sensitivity[key] << "\n";
				}

				if (m_save_all) {
					std::ofstream accuracyFile;
					accuracyFile.open(m_output_folder + "/accuracy_report.txt", std::ios_base::app); //append file

					accuracyFile << "Sens [" << key << "]   (%): " << sensitivity[key] << "\n";

					std::ofstream accuracyFileCompact;
					accuracyFileCompact.open(m_output_folder + "/accuracy_report_compact.txt", std::ios_base::app); //append file

					accuracyFileCompact << key << " " << sensitivity[key] << "\n";
				}
			}

			// Sensitivity For all

			if (m_verbose) {
				std::cout << "Sens [ALL] (%): " << sensitivityAll << "\n";
			}

			if (m_save_all) {
				std::ofstream accuracyFile;
				accuracyFile.open(m_output_folder + "/accuracy_report.txt", std::ios_base::app); //append file

				accuracyFile << "Sens [ALL] (%): " << sensitivityAll << "\n";

				std::ofstream accuracyFileCompact;
				accuracyFileCompact.open(m_output_folder + "/accuracy_report_compact.txt", std::ios_base::app); //append file

				accuracyFileCompact << 0 << " " << sensitivityAll << "\n";
			}
			
			// Update gtsResult

			gtsResult->diceScore = diceScore;
			gtsResult->sensitivity = sensitivity;
			gtsResult->falsePositivesCount = falsePositivesCount;
			gtsResult->diceScoreAll = diceScoreAll;
			gtsResult->sensitivityAll = sensitivityAll;
			gtsResult->falsePositivesCountAll = falsePositivesAll;

			message("", "Checking accuracy", 100);
		}

		void incrementMapForCheckAccuracy(std::map< LabelsPixelType, int > &map, LabelsPixelType key)
		{
			if (map.find(key) == map.end()) {
				map[key] = 1;
			}
			else {
				map[key] += 1;
			}
		}

		// Util

		void addCoordinatorMapsIfNecessary()
		{
			if (m_do_coordinate_maps)
			{				
				for (const InputImagePointer& image : GenerateCoordinateMaps<InputImageType>(m_input_images[0]))
				{
					m_input_images.push_back(image);
				}
			}
		}

		/** Coordinate maps means that for each dimension, a 3D map is generated
		 *  where the pixel values are distance to the leftmost pixels of that dimension. */
		template <class ImageType>
		static std::vector<typename ImageType::Pointer>
		GenerateCoordinateMaps(typename ImageType::Pointer reference)
		{
			typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;

			std::vector<typename ImageType::Pointer> outputVector;
			std::vector<Iterator> outputIterators;

			for (int i=0; i<ImageType::ImageDimension; i++)
			{
			  outputVector.push_back(ItkUtilGTS::initializeOutputImageBasedOn<ImageType>(reference));
			  outputIterators.push_back(
			    Iterator(outputVector[i], outputVector[i]->GetLargestPossibleRegion())
			  );
			  // outputIterators[i].GoToBegin();
			}

			Iterator iter_i(reference, reference->GetLargestPossibleRegion()); 

			// These needs to be initialized beforehand for optimization
			itk::Index<ImageType::ImageDimension> index;
			int i;

			for (iter_i.GoToBegin(); !iter_i.IsAtEnd(); 
			     ++iter_i)
			{
			  if (iter_i.Get() == 0) { continue; } // Save some time
			  index = iter_i.GetIndex();
			  
			  for (i=0; i<ImageType::ImageDimension; i++)
			  {
			    outputIterators[i].SetIndex(index);
			    // std::cout << index << "," << outputIterators[i].GetIndex() << "\n";
			    outputIterators[i].Set(index[i]);
			    // ++outputIterators[i];
			  }
			}
			return outputVector;
		}

		void shrinkImagesIfNecessary()
		{
			auto mriImagesKeys = getMapKeyset(m_input_images_MRI);

			if (m_input_images.size() != 0) {
				m_reference_image = m_input_images[0];
			}
			else {
				m_reference_image = m_input_images_MRI[mriImagesKeys[0]];
			}

			typename InputImageType::RegionType region = m_reference_image->GetLargestPossibleRegion();
			typename InputImageType::SizeType size = region.GetSize();
			typename InputImageType::SpacingType spacing = m_reference_image->GetSpacing();

			long totalSize = 1;

			itk::Index< Dimensions > centralPixel;
			std::string imageSizeString = "Image dimensions: ";
			for (int i = 0; i < Dimensions; i++) {
				totalSize *= size[i];
				imageSizeString += std::to_string(size[i]) + " ";
				centralPixel[i]  = size[i] / 2;
			}

			imageSizeString += "(Total: " + std::to_string(totalSize) + ")\n";
			message(imageSizeString);

			if (totalSize <= DEFAULT_MAX_IMAGE_SIZE_BEFORE_RESAMPLING) { return; }
			m_were_images_shrunk = true;

			itk::Point< double, Dimensions > centralPoint;
			for (int i = 0; i < Dimensions; i++) {
				centralPoint[i] = centralPixel[i];
			}

			using ScaleTransformType = itk::ScaleTransform< double, Dimensions >;
			typename ScaleTransformType::Pointer scaleTransform = ScaleTransformType::New();

			typename ScaleTransformType::ParametersType parameters = scaleTransform->GetParameters();
			float scale = 0.5;
			for (int i = 0; i < Dimensions; i++) {
				parameters[i] = scale; //???
			}

			scaleTransform->SetParameters(parameters);
			scaleTransform->SetCenter(centralPoint);

			// For images

			using LinearInterpolatorType = itk::LinearInterpolateImageFunction< InputImageType, double >;
			typename LinearInterpolatorType::Pointer interpolator = LinearInterpolatorType::New();

			using ResampleFilterType = itk::ResampleImageFilter< InputImageType, InputImageType >;
			typename ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();

			resampleFilter->SetTransform(scaleTransform);
			resampleFilter->SetInterpolator(interpolator);
			resampleFilter->SetSize(size);
			resampleFilter->SetOutputSpacing(spacing);
			
			for (size_t i = 0; i < m_input_images.size(); i++)
			{
				resampleFilter->SetInput(m_input_images[i]);
				m_input_images[i] = resampleFilter->GetOutput();
			}

			for (auto key : mriImagesKeys)
			{
				resampleFilter->SetInput(m_input_images_MRI[key]);
				m_input_images_MRI[key] = resampleFilter->GetOutput();
			}
			
			// For labels

			using LinearInterpolatorLabelsType = itk::LinearInterpolateImageFunction< LabelsImageType, double >;
			typename LinearInterpolatorLabelsType::Pointer interpolatorLabels = LinearInterpolatorLabelsType::New();

			using ResampleFilterLabelsType = itk::ResampleImageFilter< LabelsImageType, LabelsImageType >;
			typename ResampleFilterLabelsType::Pointer resampleFilterLabels = ResampleFilterLabelsType::New();

			resampleFilterLabels->SetTransform(scaleTransform);
			resampleFilterLabels->SetInterpolator(interpolatorLabels);
			resampleFilterLabels->SetSize(size);
			resampleFilterLabels->SetOutputSpacing(spacing);

			resampleFilterLabels->SetInput(m_labels_image);
			resampleFilterLabels->Update();
			m_labels_image = resampleFilterLabels->GetOutput();
			
			// Remove these
			if (m_input_images.size() != 0) { writeImage<InputImageType>(m_input_images[0], "downsampled_example"); }
			else { writeImage<InputImageType>(m_input_images_MRI[mriImagesKeys[0]], "downsampled_example"); }
		}

		template <class TImageType>
		typename TImageType::Pointer unshrinkImageIfNecessary(typename TImageType::Pointer image, typename TImageType::Pointer referenceImage)
		{

		}

		void addToInputImagesMRI()
		{
			auto keys = getMapKeyset(m_input_images_MRI);
			
			for (auto key : keys) {
				m_input_images.push_back(m_input_images_MRI[key]);
			}
		}

		bool validState(std::shared_ptr<Result> gtsResult) {
			switch (m_mode) {
			case SVM_PSEUDO:       // Same as GEOTRAIN
			case SVM_LABELS:       // Same as GEOTRAIN
			case REVERSE_GEOTRAIN: // Same as GEOTRAIN
			case GEOTRAIN:
			{
				if (m_input_images.size() + m_input_images_MRI.size() == 0 || m_labels_image == nullptr) {
					std::string errorMessage = "Need to supply input images and labels image.";
					errorOccured(errorMessage);
					gtsResult->errorMessage = errorMessage;
					gtsResult->ok = false;
					return false;
				}
			} break;
			case RF:              // Same as AGD_RF
			case AGD_RF:
			{
				// TODO: fill this
			} break;
			case RF_AUTO:
			case AGD_RF_AUTO:
			{
				// TODO: fill this
			} break;
			case AGD:
			{
				if (m_input_images.size() + m_input_images_MRI.size() == 0 || m_labels_image == nullptr) {
					std::string errorMessage = "Need to supply input images and labels image and a list containing the label of interest of each image for mode agd";
					errorOccured(errorMessage);
					gtsResult->errorMessage = errorMessage;
					gtsResult->ok = false;
					return false;
				}
				/*if (m_input_images.size() != m_labels_list.size()) {
					errorOccured("Number of input images must be equal to the number of agdlabels for mode agd");
					return false;
				}*/
			} break;
			case LABELS_THRESHOLD:
			{
				if (m_input_images.size() == 0) {
					std::string errorMessage = "Need to supply input images and list of labels for values<=threshold for mode labelsthres";
					errorOccured(errorMessage);
					gtsResult->errorMessage = errorMessage;
					gtsResult->ok = false;
					return false;
				}
			} break;
			case SEGMENT:
			{
				if ((m_input_images.size() != 1) || (m_labels_image == nullptr))
				{
					std::string errorMessage = "Need to supply one input image, one labels image and one label of interest for mode segment";
					errorOccured(errorMessage);
					gtsResult->errorMessage = errorMessage;
					gtsResult->ok = false;
					return false;
				}
			} break;
			case CHANGE_LABELS:
			{
				if (m_labels_image == nullptr) {
					std::string errorMessage = "Need to supply labels image for mode changelabels";
					errorOccured(errorMessage);
					gtsResult->errorMessage = errorMessage;
					gtsResult->ok = false;
					return false;
				}
			} break;
			case GENERATE_CONFIG:
			{
				if ((m_pretrained_models_paths.size() == 0) || (m_importance_values.size() == 0)) {
					std::string errorMessage = "Need to supply lists of input models and improtance values for each model for mode generateconfig.";
					errorOccured(errorMessage);
					gtsResult->errorMessage = errorMessage;
					gtsResult->ok = false;
					return false;
				}
				if (m_pretrained_models_paths.size() != m_importance_values.size()) {
					std::string errorMessage = "# of input models and # of importance values must be equal for mode generateconfig.";
					errorOccured(errorMessage);
					gtsResult->errorMessage = errorMessage;
					gtsResult->ok = false;
					return false;
				}
			} break;
			case CHECK_ACCURACY:
			{
				if ((m_labels_image == nullptr) || (m_ground_truth == nullptr)) {
					std::string errorMessage = "Need to supply labels image and ground truth image for mode accuracy";
					errorOccured(errorMessage);
					gtsResult->errorMessage = errorMessage;
					gtsResult->ok = false;
					return false;
				}
			} break;
			}

			return true;
		}

		/** Checks if the labels input image is valid. Also changes m_label_of_interest/m_label_TC/etc.. if necessary */
		bool validLabels(std::shared_ptr<Result> gtsResult, std::unordered_map<LabelsPixelType, int> labelsCountMap, int limitedNumberOfClasses = 0)
		{
			message("Labels detected:      ");
			auto keys = getMapKeyset(labelsCountMap);
			std::sort(keys.begin(), keys.end());

			for (LabelsPixelType key : keys) {
				message(std::to_string(key) + " ");
			}
			message("\n");

			if (limitedNumberOfClasses == 1) {
				if (labelsCountMap.size() == 0) {
					std::string errorMessage = "No labels were found in the image";
					errorOccured(errorMessage);
					gtsResult->errorMessage = errorMessage;
					gtsResult->ok = false;
					return false;
				}

				if (labelsCountMap.find(m_label_of_interest) == labelsCountMap.end()) {
					message("Label of interest (" + std::to_string(m_label_of_interest) + ") was not found.\n", "");
					auto keys = getMapKeyset(labelsCountMap);
					m_label_of_interest = labelsCountMap[keys[0]];
					message("\tUsing " + std::to_string(m_label_of_interest) + " instead.\n"
					      + "\tPlease provide a different one if you don't want this one to be used\n", "");
				}
			}
			else if (limitedNumberOfClasses == 2) {
				if (labelsCountMap.size() != 2) {
					std::string errorMessage = "There should only be two (non-zero) labels for 2-class modes";
					errorOccured(errorMessage);
					gtsResult->errorMessage = errorMessage;
					gtsResult->ok = false;
					return false;
				}

				if (labelsCountMap.find(m_label_of_interest) == labelsCountMap.end()) {
					message("Label of interest (" + std::to_string(m_label_of_interest) + ") was not found.\n", "");
					auto keys = getMapKeyset(labelsCountMap);
					m_label_of_interest = labelsCountMap[keys[0]];
					message("\tUsing " + std::to_string(m_label_of_interest) + " instead.\n"
					      + "\tPlease provide a different one if you don't want this one to be used\n", "");
				}

				if (m_input_images_MRI.size() != 0) {
					// For MRI. For each region set the label to zero if they were not included in the sample
					if (labelsCountMap.find(m_label_HT) == labelsCountMap.end()) {
						m_label_HT = 0;
					}
					if (labelsCountMap.find(m_label_TC) == labelsCountMap.end()) {
						m_label_TC = 0;
					}
					if (labelsCountMap.find(m_label_ET) == labelsCountMap.end()) {
						m_label_ET = 0;
					}
					if (labelsCountMap.find(m_label_ED) == labelsCountMap.end()) {
						m_label_ED = 0;
					}

					if (m_label_TC == 0 && m_label_ET == 0 && m_label_ED == 0 && m_label_HT == 0) {
						message("Note: No labels were found for TC/ET/ED/HT even though images were supplied through the special modalities input.\n");
						message("      Will proceed with performing AGD on all labels (if applicable)\n");
						addToInputImagesMRI();
						m_input_images_MRI.clear();
					}
				}
			}
			else {
				// n-classes
				if (labelsCountMap.size() < 2) {
					std::string errorMessage = 
						std::string("There should be at least two (non-zero) labels for multiclass modes. ") + 
						std::string("Please remember to draw a label for the background tissue too.");
					errorOccured(errorMessage);
					gtsResult->errorMessage = errorMessage;
					gtsResult->ok = false;
					return false;
				}
				
				if (m_input_images_MRI.size() != 0) {
					// For MRI. For each region set the label to zero if they were not included in the sample
					if (labelsCountMap.find(m_label_HT) == labelsCountMap.end()) {
						m_label_HT = 0;
					}
					if (labelsCountMap.find(m_label_TC) == labelsCountMap.end()) {
						m_label_TC = 0;
					}
					if (labelsCountMap.find(m_label_ET) == labelsCountMap.end()) {
						m_label_ET = 0;
					}
					if (labelsCountMap.find(m_label_ED) == labelsCountMap.end()) {
						m_label_ED = 0;
					}

					if (m_label_TC == 0 && m_label_ET == 0 && m_label_ED == 0 && m_label_HT == 0) {
						message("Note: No labels were found for TC/ET/ED/HT even though images were supplied through the special modalities input.\n");
						message("      Will proceed with performing AGD on all labels (if applicable)\n");
						addToInputImagesMRI();
						m_input_images_MRI.clear();
					}
				}
			}

			return true;
		}

		void errorOccured(std::string msg) {
			std::cerr << "GeodesicTrainingSegmentation error: " << msg << std::endl;
		}

		std::string getFileExtension(std::string fName) {
			std::string extension = cbica::getFilenameExtension(fName);
			std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);

			return extension;
		}

		template< typename U, typename V >
		std::vector< U > getMapKeyset(std::map< U, V > map)
		{
			// Find a list of the different keys
			std::vector< U > res;

			res.reserve(map.size());
			for (auto const& imap : map) {
				res.push_back(imap.first);
			}

			return res;
		}

		template< typename U, typename V >
		std::vector< U > getMapKeyset(std::unordered_map< U, V > map)
		{
			// Find a list of the different keys
			std::vector< U > res;

			res.reserve(map.size());
			for (auto const& imap : map) {
				res.push_back(imap.first);
			}

			return res;
		}

		std::string getModalityName(MODALITY_MRI modality) {
			switch (modality) {
				case FLAIR:
					return "flair";
				case T1:
					return "t1";
				case T1CE:
					return "t1ce";
				case T2:
					return "t2";
			}
			return "";
		}

		// I/O

		template<class TImageType>
		typename TImageType::Pointer readImage(std::string filename)
		{
			return cbica::ReadImage<TImageType>(filename);

			//if (((m_file_extension == ".nii.gz") || (m_file_extension == ".nii") 
			//	|| (m_file_extension == ".dcm") || (m_file_extension == ".dicom")) && (Dimensions > 2)) {
			//	return cbica::ReadImage<TImageType>(filename);
			//}
			//else {
			//	typedef itk::ImageFileReader<TImageType> ReaderType;
			//	typename ReaderType::Pointer reader = ReaderType::New();

			//	reader->SetFileName(filename);
			//	reader->Update();
			//	
			//	typename TImageType::Pointer image;
			//	image->Graft(reader->GetOutput());
			//	return image;
			//}
		}

		template<class TImageType>
		void writeImage(typename TImageType::Pointer image, std::string filename)
		{
			bool isOutputSegmentation = (filename == "labels_res" || filename == "labels_res_renamed");
			
			if (m_save_only_seg)
			{
				filename = m_save_only_seg_name;
			}

			if (m_save_all || (m_save_only_seg && isOutputSegmentation)) {
				std::string fileFullPath = m_output_folder + "/" + filename + m_file_extension;
				
				//message("Writing image to file...", "Writing to file");
				cbica::WriteImage<TImageType>(image, fileFullPath);
				//message("Writing image to file...", "Writing to file", 100);
				
				//std::string fileFullPath = m_output_folder + "/" + filename;

				//if (((m_file_extension == ".nii.gz") || (m_file_extension == ".nii") 
				//	|| (m_file_extension == ".dcm") || (m_file_extension == ".dicom")) && (Dimensions > 2)) {
				//	cbica::WriteImage<TImageType>(image, fileFullPath);
				//}
				//else {
				//	// Unknown image type
				//	typedef typename itk::ImageFileWriter<TImageType> WriterType;
				//	typename WriterType::Pointer writer = WriterType::New();
				//	writer->SetInput(image);
				//	writer->SetFileName(fileFullPath);
				//	writer->Update();
				//}
			}
		}

		/** This function is meant to be overriden if someone who uses the library wishes to */
		virtual void progressUpdate(std::string message, int progress) { }

		void message(std::string message, std::string shortMessage = "", int progress = -1)
		{
			if (m_verbose) {
				if (message != "") {
					std::string printMessage = (progress != -1) ? "\r" : "";
					printMessage += message;

					if (progress == 100) {
						printMessage += "finished\n";
					}
					else if (progress != -1) {
						printMessage += " [" + std::to_string(progress) + "%]";
					}
					
					std::cout << printMessage;
				}
			}

			// Update progress
			if (shortMessage != "") {
				if (progress == -1) {
					progressUpdate("GTS: " + shortMessage, 0);
				}
				else {
					progressUpdate("GTS: " + shortMessage, progress);
				}
			}
		}

		// For timer

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

#endif
