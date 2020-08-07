/**
\file  cbicaITKUtilities.h

\brief Some basic utility functions.

Dependecies: ITK (module_review, module_skullstrip enabled), OpenMP

https://www.cbica.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2015 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.cbica.upenn.edu/sbia/software/license.html

*/
#pragma once

#include <algorithm>
#include <functional>
#include <cmath>
#include <array>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "itkImage.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"

#include "itkHistogramMatchingImageFilter.h"
#include "itkAdaptiveHistogramEqualizationImageFilter.h"

#include "itkConnectedThresholdImageFilter.h"
#include "itkOrientImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkAddImageFilter.h"

//#include "itkMultiResolutionPDEDeformableRegistration.h"
//#include "itkDemonsRegistrationFilter.h"
//#include "itkDiffeomorphicDemonsRegistrationFilter.h"
//#include "itkFastSymmetricForcesDemonsRegistrationFilter.h"
//#include "itkSymmetricForcesDemonsRegistrationFilter.h"

#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLabelOverlapMeasuresImageFilter.h"
//#include "itkWindowedSincInterpolateImageFunction.h"

#include "itkTestingComparisonImageFilter.h"

#include "itkStripTsImageFilter.h"
#include "itkMaskImageFilter.h"

#include "cbicaUtilities.h"
#include "cbicaITKImageInfo.h"
//#include "cbicaITKSafeImageIO.h"

#include "HausdorffDistance.h"

#include "gdcmMD5.h"
#include "gdcmReader.h"

#include "DicomIOManager.h"

#include "itkNaryFunctorImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"

#include "itkJoinSeriesImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkStatisticsImageFilter.h"

using ImageTypeFloat3D = itk::Image< float, 3 >;
//unsigned int RmsCounter = 0;
//double MaxRmsE[4] = { 0.8, 0.75, 0.4, 0.2 };

enum DeformRegType
{
  Demons, DiffeomorphicDemons, SymmetricForcesDemons, FastSymmetricForcesDemons
};

enum InterpolatorType
{
  Linear, NearestNeighbor, BSpline
};

//! Helper class for computation
template <class TInputImage, class TOutputImage>
class NaryLabelVotingFunctor
{
public:
  typedef NaryLabelVotingFunctor<TInputImage, TOutputImage> Self;
  typedef typename TInputImage::PixelType InputPixelType;
  typedef typename TOutputImage::PixelType OutputPixelType;
  typedef std::vector<OutputPixelType> LabelArray;

  NaryLabelVotingFunctor(const LabelArray &labels)
    : m_LabelArray(labels), m_Size(labels.size()) {}

  NaryLabelVotingFunctor() : m_Size(0) {}


  OutputPixelType operator() (const std::vector<InputPixelType> &pix)
  {
    InputPixelType best_val = pix[0];
    int best_index = 0;
    for (int i = 1; i < m_Size; i++)
      if (pix[i] > best_val)
      {
        best_val = pix[i];
        best_index = i;
      }

    return m_LabelArray[best_index];
  }

  bool operator != (const Self &other)
  {
    return other.m_LabelArray != m_LabelArray;
  }

protected:
  LabelArray m_LabelArray;
  int m_Size;
};

/*
\namespace cbica
\brief Namespace for differentiating functions written for internal use
*/
namespace cbica
{
  //! The image type
  enum ImageModalityType
  {
    IMAGE_TYPE_UNDEFINED = 0, IMAGE_TYPE_T1, IMAGE_TYPE_T1CE, IMAGE_TYPE_T2,
    IMAGE_TYPE_T2FLAIR, IMAGE_TYPE_AX, IMAGE_TYPE_FA, IMAGE_TYPE_RAD, IMAGE_TYPE_TR,
    IMAGE_TYPE_PERFUSION, IMAGE_TYPE_DTI, IMAGE_TYPE_RECURRENCE_OUTPUT, IMAGE_TYPE_PP, IMAGE_TYPE_CT,
    IMAGE_TYPE_PET, IMAGE_TYPE_PSR, IMAGE_TYPE_PH, IMAGE_TYPE_RCBV, IMAGE_TYPE_SEG,
    IMAGE_TYPE_ATLAS, IMAGE_TYPE_PARAMS, IMAGE_TYPE_SUDOID, IMAGE_TYPE_NEAR, IMAGE_TYPE_FAR,
    IMAGE_MAMMOGRAM, IMAGE_TYPE_FEATURES
  };

  //! The modality strings that are used in the GUI 
  static const char ImageModalityString[ImageModalityType::IMAGE_TYPE_FEATURES + 1][15] =
  {
    "DEF", "T1", "T1Gd", "T2",
    "FLAIR", "DTI_AX", "DTI_FA", "DTI_RAD", "DTI_TR",
    "PERFUSION", "DTI", "REC", "PP", "CT",
    "PET", "pSR", "PH", "RCBV", "SEG",
    "ATLAS", "PARAMS", "SUDOID", "NEAR", "FAR",
    "FFDM", "FEAT"
  };

  /**
  \brief Guess Image Type

  \param str String to guess
  \return deduced type
  */
  inline int guessImageType(const std::string &fileName)
  {
    int ImageSubType = ImageModalityType::IMAGE_TYPE_UNDEFINED;
    std::string fileName_wrap = fileName;
    std::transform(fileName_wrap.begin(), fileName_wrap.end(), fileName_wrap.begin(), ::tolower);
    auto ext = cbica::getFilenameExtension(fileName_wrap, false);

    // lambda function to check the different string combinations in a file for modality check
    auto modalityCheckerFunction = [&](std::string baseModalityStringToCheck)
    {
      if (
        (fileName_wrap.find(baseModalityStringToCheck + ext) != std::string::npos) ||
        (fileName_wrap.find("_" + baseModalityStringToCheck) != std::string::npos) ||
        (fileName_wrap.find(baseModalityStringToCheck + "_") != std::string::npos) ||
        (fileName_wrap.find("." + baseModalityStringToCheck + ".") != std::string::npos) ||
        (fileName_wrap.find(baseModalityStringToCheck) != std::string::npos)
        )
      {
        return true;
      }
    };

    // using the lambda to figure out the modality
    if (
      modalityCheckerFunction("t1ce") ||
      modalityCheckerFunction("t1gd") ||
      modalityCheckerFunction("t1gad") ||
      modalityCheckerFunction("t1-ce") ||
      modalityCheckerFunction("t1-gd") ||
      modalityCheckerFunction("t1-gad")
      )
    {
      ImageSubType = ImageModalityType::IMAGE_TYPE_T1CE;
    }
    else if (modalityCheckerFunction("t1"))
    {
      ImageSubType = ImageModalityType::IMAGE_TYPE_T1;
    }
    else if (modalityCheckerFunction("t2"))
    {
      if ((fileName_wrap.find("flair") != std::string::npos)) // if there is "flair" present in this case, simply return flair
      {
        ImageSubType = ImageModalityType::IMAGE_TYPE_T2FLAIR;
      }
      else
      {
        ImageSubType = ImageModalityType::IMAGE_TYPE_T2;
      }
    }
    else if (modalityCheckerFunction("flair"))
    {
      ImageSubType = ImageModalityType::IMAGE_TYPE_T2FLAIR;
    }
    else if (modalityCheckerFunction("dti") || modalityCheckerFunction("b0"))
    {
      ImageSubType = ImageModalityType::IMAGE_TYPE_DTI;
    }
    else if (modalityCheckerFunction("radial") || modalityCheckerFunction("rad"))
    {
      ImageSubType = ImageModalityType::IMAGE_TYPE_RAD;
    }
    else if (modalityCheckerFunction("axial") || modalityCheckerFunction("ax"))
    {
      ImageSubType = ImageModalityType::IMAGE_TYPE_AX;
    }
    else if (modalityCheckerFunction("fractional") || modalityCheckerFunction("fa"))
    {
      ImageSubType = ImageModalityType::IMAGE_TYPE_FA;
    }
    else if (modalityCheckerFunction("trace") || modalityCheckerFunction("tr"))
    {
      ImageSubType = ImageModalityType::IMAGE_TYPE_TR;
    }
    else if (modalityCheckerFunction("rcbv"))
    {
      ImageSubType = ImageModalityType::IMAGE_TYPE_RCBV;
    }
    else if (modalityCheckerFunction("psr"))
    {
      ImageSubType = ImageModalityType::IMAGE_TYPE_PSR;
    }
    else if (modalityCheckerFunction("ph"))
    {
      ImageSubType = ImageModalityType::IMAGE_TYPE_PH;
    }
    else if (modalityCheckerFunction("perf") || modalityCheckerFunction("dsc"))
    {
      ImageSubType = ImageModalityType::IMAGE_TYPE_PERFUSION;
    }
    else if (modalityCheckerFunction("ct2pet") || modalityCheckerFunction("ct"))
    {
      ImageSubType = ImageModalityType::IMAGE_TYPE_CT;
    }
    else if (modalityCheckerFunction("pet"))
    {
      ImageSubType = ImageModalityType::IMAGE_TYPE_PET;
    }
    else if (
      modalityCheckerFunction("labelmap") ||
      modalityCheckerFunction("label-map") ||
      modalityCheckerFunction("segmentation") ||
      modalityCheckerFunction("annotation") ||
      modalityCheckerFunction("label") ||
      modalityCheckerFunction("roi"))
    {
      ImageSubType = ImageModalityType::IMAGE_TYPE_SEG;
    }
    return ImageSubType;
  }

  /**
  \brief Calculate and preserve the mask indeces

  \param inputModalitiesAndImages A collection of images which are stored in a per-modality basis (each entry corresponds to a subject, whose entries contain different modalities)
  \return A collection of indeces which constitute the non-zero locations per modality (each entry corresponds to a subject, which contains the locations of non-zero pixel values for all modalities)
  */
  template< class TImageType = ImageTypeFloat3D >
  std::vector< std::vector< typename TImageType::IndexType > > CreateMaskIndeces(const std::vector< std::vector< typename TImageType::Pointer > > &inputModalitiesAndImages)
  {
    std::vector< std::vector< typename TImageType::IndexType > > returnMaskIndeces;
    returnMaskIndeces.resize(inputModalitiesAndImages.size()); // pre-allocate data for speed

    // start data processing
    // made parallel for efficiency
    int threads = omp_get_max_threads(); // obtain maximum number of threads available on machine  
    //threads > inputModalitiesAndImages.size() ? threads = inputModalitiesAndImages.size() : threads = threads;
    //#pragma omp parallel for num_threads(threads)
    for (int i = 0; i < inputModalitiesAndImages.size(); i++)
    {
      std::vector< float > means;
      std::vector< typename TImageType::IndexType > tempIndeces;
      typename TImageType::SizeType size = inputModalitiesAndImages[i][0]->GetLargestPossibleRegion().GetSize();
      size_t totalImageSize = size[0] * size[1] * size[2];
      means.resize(totalImageSize);
      //std::fill(means.begin(), means.end(), 0);
      means.assign(totalImageSize, 0);

      for (size_t j = 0; j < inputModalitiesAndImages[i].size(); j++)
      {
        std::vector< float > tempVec;
        itk::ImageRegionIterator< TImageType > it(inputModalitiesAndImages[i][j], inputModalitiesAndImages[i][j]->GetLargestPossibleRegion());
        it.GoToBegin();

        while (!it.IsAtEnd())
        {
          tempVec.push_back(it.Get());
          tempIndeces.push_back(it.GetIndex());
          ++it;
        }
        if (tempVec.size() == means.size())
        {
          std::transform(means.begin(), means.end(), tempVec.begin(), means.begin(), std::plus< float >()); // add tempVec to means dector
        }
        else
        {
          std::cerr << "Mean vector calculation error.\n";
          exit(EXIT_FAILURE);
        }
      } // loop over all subjects in each modality

      //std::transform(means.begin(), means.end(), means.begin(), std::bind1st(std::divides< float >(), means.size())); // divide entire means vector by its size
      std::vector< typename TImageType::IndexType > tempMaskIndeces;
      for (size_t j = 0; j < means.size(); j++)
      {
        if (means[j] > 0)
        {
          tempMaskIndeces.push_back(tempIndeces[j]); // store indeces of non-zero mean values
        }
      }
      returnMaskIndeces[i] = tempMaskIndeces;
    } // loop over all modalities

    return returnMaskIndeces;
  }

  /**
  \brief Get Pixel Values of specified indeces of input Image

  \param inputImage The input image in itk::Image format
  \param indeced The indeces from which pixel values need to be extracted
  \return Vector of values whose data type is the same as image type
  */
  template < typename TImageType = ImageTypeFloat3D >
  std::vector< typename TImageType::PixelType > GetPixelValuesFromIndeces(const typename TImageType::Pointer inputImage, const std::vector< typename TImageType::IndexType > &indeces)
  {
    std::vector< typename TImageType::PixelType > returnVector;
    returnVector.resize(indeces.size()); // pre-allocation done for speed

    typedef itk::ImageRegionIterator< TImageType > IteratorType;
    IteratorType imageIterator(inputImage, inputImage->GetBufferedRegion());

    // made parallel for efficiency
    int threads = omp_get_max_threads(); // obtain maximum number of threads available on machine  
    //threads > returnVector.size() ? threads = returnVector.size() : threads = threads;
    //#pragma omp parallel for num_threads(threads)
    for (int i = 0; i < returnVector.size(); i++)
    {
      imageIterator.SetIndex(indeces[i]);
      returnVector[i] = imageIterator.Get();
    }

    return returnVector;
  }

  /**
  \brief Wrap of GetPixelValues
  */
  template < typename TImageType = ImageTypeFloat3D >
  std::vector< typename TImageType::PixelType > ExtractPixelValuesFromIndeces(const typename TImageType::Pointer inputImage, const std::vector< typename TImageType::IndexType > &indeces)
  {
    return GetPixelValuesFromIndeces< TImageType >(inputImage, indeces);
  }

  /**
  \brief Check properties of 2 images to see if they are defined in the same space.
  */
  template< typename TImageType >
  inline bool ImageSanityCheck(
    const typename TImageType::Pointer image1, 
    const typename TImageType::Pointer image2,
    const float nifti2dicomTolerance = 0.0,
    const float nifti2dicomOriginTolerance = 0.0)
  {
    auto size_1 = image1->GetLargestPossibleRegion().GetSize();
    auto size_2 = image2->GetLargestPossibleRegion().GetSize();

    auto origin_1 = image1->GetOrigin();
    auto origin_2 = image2->GetOrigin();

    auto spacing_1 = image1->GetSpacing();
    auto spacing_2 = image2->GetSpacing();

    auto directions_1 = image1->GetDirection();
    auto directions_2 = image2->GetDirection();

    if (directions_1.ColumnDimensions != directions_2.ColumnDimensions)
    {
      std::cerr << "Column dimension mismatch for directions.\n";
      return false;
    }
    if (directions_1.RowDimensions != directions_2.RowDimensions)
    {
      std::cerr << "Row dimension mismatch for directions.\n";
      return false;
    }
    // check with tolerance
    for (size_t i = 0; i < directions_1.RowDimensions; i++)
    {
      for (size_t j = 0; j < directions_1.ColumnDimensions; j++)
      {
        if (directions_1[i][j] != directions_2[i][j])
        {
          auto percentageDifference = std::abs(directions_1[i][j] - directions_2[i][j]) * 100 / std::abs(directions_1[i][j]);
          if (percentageDifference > nifti2dicomTolerance)
          {
            std::cerr << "Direction mismatch > " << nifti2dicomTolerance << 
              "% at location '[" << i << "," << j << "]' of direction matrix.\n";

            std::cout << "Direction matrix for input 1:\n" << directions_1 << "\n" <<
              "Direction matrix for input 2:\n" << directions_2 << "\n";
            return false;
          }
          else
          {
            std::cout << "Ignoring direction difference of '" <<
              percentageDifference << "%' in location '[" <<
              i << "," << j << "]' of direction matrix.\n";
          }
        }
      }
    }

    for (size_t d = 0; d < TImageType::ImageDimension; d++)
    {
      if (size_1[d] != size_2[d])
      {
        std::cerr << "Size mismatch at dimension '" << d << "'\n";
        return false;
      }
      if (origin_1[d] != origin_2[d])
      {
        auto percentageDifference = std::abs(origin_1[d] - origin_2[d]) * 100 / std::abs(origin_1[d]);
        if (percentageDifference > nifti2dicomOriginTolerance)
        {
          std::cerr << "Origin mismatch > " << percentageDifference <<
            percentageDifference << "%' in axis '" << d << "'.\n";

          std::cout << "Origin for input 1:\n" << origin_1 << "\n" <<
            "Origin for input 2:\n" << origin_2 << "\n";
          return false;
        }
        else
        {
          std::cout << "Ignoring origin difference of '" <<
            percentageDifference << "%' in axis '" << d << "'.\n";
        }
      }
      if (spacing_1[d] != spacing_2[d])
      {
        auto percentageDifference = std::abs(spacing_1[d] - spacing_2[d]) * 100;
        percentageDifference /= spacing_1[d];
        if (percentageDifference > 0.0001)
        {
          std::cerr << "Spacing mismatch at dimension '" << d << "'\n";
          return false;
        }
        else
        {
          std::cout << "Ignoring spacing difference of '" <<
            percentageDifference << "%' in dimension '" <<
            d << "'\n";
        }
      }
    }

    return true;
  }

  /**
  \brief Check properties of 2 images to see if they are defined in the same space.

  Checks are done based on cbica::ImageInfo class
  */
  inline bool ImageSanityCheck(const std::string &image1, const std::string &image2, bool FourDImageCheck = false)
  {
    auto imageInfo1 = cbica::ImageInfo(image1);
    auto imageInfo2 = cbica::ImageInfo(image2);

    auto dims = imageInfo1.GetImageDimensions();

    if (FourDImageCheck)
    {
      dims = 3; // this is for 4D images only
    }
    else // do the check when FourDImageCheck is disabled
    {
      if (imageInfo1.GetImageDimensions() != imageInfo2.GetImageDimensions())
      {
        std::cout << "The dimensions of the image_1 (" << image1 << ") and image_2 (" << image2 << ") doesn't match.\n";
        return false;
      }
    }

    // check size, spacing and origin information as well

    auto imageSize1 = imageInfo1.GetImageSize();
    auto imageSize2 = imageInfo2.GetImageSize();

    auto imageSpacing1 = imageInfo1.GetImageSpacings();
    auto imageSpacing2 = imageInfo2.GetImageSpacings();

    auto imageOrigin1 = imageInfo1.GetImageOrigins();
    auto imageOrigin2 = imageInfo2.GetImageOrigins();

    auto imageDirs1 = imageInfo1.GetImageDirections();
    auto imageDirs2 = imageInfo1.GetImageDirections();

    for (size_t d = 0; d < dims; d++)
    {
      if (imageSize1[d] != imageSize2[d])
      {
        std::cout << "The size in dimension[" << d << "] of the image_1 (" << image1 << ") and image_2 (" << image2 << ") doesn't match.\n";
        return false;
      }
      if (imageSpacing1[d] != imageSpacing2[d])
      {
        auto percentageDifference = std::abs(imageSpacing1[d] - imageSpacing2[d]) * 100;
        percentageDifference /= imageSpacing1[d];
        if (percentageDifference > 0.0001)
        {
          std::cerr << "Spacing mismatch at dimension '" << d << "'\n";
          return false;
        }
        else
        {
          std::cout << "Ignoring spacing difference of '" <<
            percentageDifference << "%' in dimension '" <<
            d << "'\n";
        }
      }
      if (imageOrigin1[d] != imageOrigin2[d])
      {
        std::cout << "The origin in dimension[" << d << "] of the image_1 (" << image1 << ") and image_2 (" << image2 << ") doesn't match.\n";
        return false;
      }

      if (imageDirs1[d].size() != imageDirs2[d].size())
      {
        std::cout << "The direction in dimension[" << d << "] of the image_1 (" << image1 << ") and image_2 (" << image2 << ") doesn't match.\n";
        return false;
      }
      else
      {
        for (size_t i = 0; i < imageDirs1[d].size(); i++)
        {
          if (imageDirs1[d][i] != imageDirs2[d][i])
          {
            auto percentageDifference = std::abs(imageDirs1[d][i] - imageDirs2[d][i]) * 100 / imageDirs1[d][i];
            if (percentageDifference > 0.0001)
            {
              std::cerr << "Direction mismatch at dimension '" << d << "'\n";
              return false;
            }
            else
            {
              std::cout << "Ignoring direction difference of '" <<
                percentageDifference << "%' in dimension '" <<
                d << "'\n";
            }
          }
        }
      }
    }

    return true;
  }

  /**
  \brief This function returns a joined N-D image with an input of a vector of (N-1)-D images

  Uses the itk::JoinSeriesImageFilter to accomplish this

  \param inputImage The vector of images from which the larger image is to be extracted
  \param newSpacing The spacing in the new dimension
  */
  template< class TInputImageType, class TOutputImageType >
  typename TOutputImageType::Pointer GetJoinedImage(std::vector< typename TInputImageType::Pointer > &inputImages, double newSpacing = 1.0)
  {
    if (TOutputImageType::ImageDimension - 1 != TInputImageType::ImageDimension)
    {
      std::cerr << "Only works when input and output image dimensions are N and (N+1), respectively.\n";
      //return typename TOutputImageType::New();
      exit(EXIT_FAILURE);
    }
    auto joinFilter = /*typename*/ itk::JoinSeriesImageFilter< TInputImageType, TOutputImageType >::New();
    joinFilter->SetSpacing(newSpacing);

    for (size_t N = 0; N < inputImages.size(); N++)
    {
      if (!ImageSanityCheck< TInputImageType >(inputImages[0], inputImages[N]))
      {
        std::cerr << "Image Sanity check failed in index '" << N << "'\n";
        //return typename TOutputImageType::New();
        exit(EXIT_FAILURE);
      }
      joinFilter->SetInput(N, inputImages[N]);
    }
    try
    {
      joinFilter->Update();
    }
    catch (const std::exception& e)
    {
      std::cerr << "Joining failed: " << e.what() << "\n";
    }
    return joinFilter->GetOutput();
  }

  /**
  \brief This function returns a vector of (N-1)-D images with an input of an N-D image

  Uses the itk::ExtractImageFilter to accomplish this

  \param inputImage The larger image series from which the sub-images in the specified axis are extracted
  \param axisToExtract The axis along with the images are to be extracted from; defaults to TInputImageType::ImageDimension - for extraction along Z, use '3'
  \param directionsCollapseIdentity Whether direction cosines are to be normalized to identity or not; defaults to not
  */
  template< class TInputImageType, class TOutputImageType >
  std::vector< typename TOutputImageType::Pointer > GetExtractedImages(typename TInputImageType::Pointer inputImage,
    int axisToExtract = TInputImageType::ImageDimension - 1, bool directionsCollapseIdentity = false)
  {
    std::vector<typename TOutputImageType::Pointer> returnImages;

    if (TOutputImageType::ImageDimension != TInputImageType::ImageDimension - 1)
    {
      std::cerr << "Only works when input and output image dimensions are N and (N-1), respectively.\n";
      return returnImages;
    }
    // set the sub-image properties
    auto imageSize = inputImage->GetLargestPossibleRegion().GetSize();
    auto regionSize = imageSize;
    regionSize[axisToExtract] = 0;
    returnImages.resize(imageSize[axisToExtract]);

    typename TInputImageType::IndexType regionIndex;
    regionIndex.Fill(0);

    // loop through time points
    for (size_t i = 0; i < imageSize[axisToExtract]; i++)
    {
      regionIndex[axisToExtract] = i;
      typename TInputImageType::RegionType desiredRegion(regionIndex, regionSize);
      auto extractor = /*typename*/ itk::ExtractImageFilter< TInputImageType, TOutputImageType >::New();
      extractor->SetExtractionRegion(desiredRegion);
      extractor->SetInput(inputImage);
      if (directionsCollapseIdentity)
      {
        extractor->SetDirectionCollapseToIdentity();
      }
      else
      {
        extractor->SetDirectionCollapseToSubmatrix();
      }
      try
      {
        extractor->Update();
      }
      catch (const std::exception& e)
      {
        std::cerr << "Extracting failed: " << e.what() << "\n";
      }
      auto temp = extractor->GetOutput();
      //temp->DisconnectPipeline(); // ensure a hard copy is done 
      returnImages[i] = temp;
    }
    return returnImages;
  }

  ///**
  //\brief Get MD5 sum of a supplied file

  //\param fileName The input file
  //\return The MD5 checksum
  //*/
  //inline std::string GetMD5Sum(const std::string &fileName)
  //{
  //  gdcm::MD5 md5Computer;
  //  char digStr[1024/*MAX_PATH*/];
  //  md5Computer.ComputeFile(fileName.c_str(), digStr);
  //  return std::string(digStr);
  //}

  ///**
  //\brief Wrap of GetMD5Sum()
  //*/
  //inline std::string ComputeMD5Sum(const std::string &fileName)
  //{
  //  return GetMD5Sum(fileName);
  //}

  /**
  \brief Get the indeces of the image which are not zero
  */
  template < typename TImageType = ImageTypeFloat3D >
  std::vector< typename TImageType::IndexType > GetIndexFromNonZeroPixels(const typename TImageType::Pointer inputImage, const std::string valuesToExclude = "0")
  {
    std::vector< typename TImageType::IndexType > returnVector;

    itk::ImageRegionConstIterator< TImageType > iterator(inputImage, inputImage->GetLargestPossibleRegion());
    for (iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator)
    {
      if (iterator.Get() != 0)
      {
        returnVector.push_back(iterator.GetIndex());
      }
    }

    return returnVector;
  }

  /**
  \brief Get the indeces of the image which are not zero

  \param inputImage The input image on which the matching needs to be done
  \param referenceImage The reference image based on which the
  \param numberOfMatchPoints Governs the number of quantile values to be matched
  \param numberOfHistogramLevels Sets the number of bins used when creating histograms of the source and reference images
  */
  template < typename TImageType = ImageTypeFloat3D >
  typename TImageType::Pointer GetHistogramMatchedImage(const typename TImageType::Pointer inputImage, const typename TImageType::Pointer referenceImage,
    const int numberOfMatchPoints = 40, const int numberOfHistogramLevels = 100)
  {
    auto filter = itk::HistogramMatchingImageFilter< TImageType, TImageType >::New();
    filter->SetInput(inputImage);
    filter->SetReferenceImage(referenceImage);
    if (numberOfHistogramLevels != 100)
    {
      filter->SetNumberOfHistogramLevels(numberOfHistogramLevels);
    }
    filter->ThresholdAtMeanIntensityOn();
    filter->SetNumberOfMatchPoints(numberOfMatchPoints);
    filter->Update();

    return filter->GetOutput();
  }

  /**
  \brief Get the indeces of the image which are not zero

  \param inputImage The input image on which the matching needs to be done
  \param referenceImage The reference image based on which the
  \param alpha Ranges between 0-1; with 1 giving result same as input image and lower values behaving as unsharp filters; default = 0.3
  \param beta Ranges between 0-1; with 1 giving result same as input image and lower values behaving as unsharp filters; default = 0.3
  \param radius Ranges between 1-10 with default = 1
  */
  template < typename TImageType = ImageTypeFloat3D >
  typename TImageType::Pointer GetAdaptiveHistogramEqualizedImage(const typename TImageType::Pointer inputImage, const typename TImageType::Pointer referenceImage,
    const float alpha = 0.3, const float beta = 0.3, const float radius = 1, const int numberOfHistogramLevels = 100)
  {
    auto filter = itk::AdaptiveHistogramEqualizationImageFilter< TImageType, TImageType >::New();
    filter->SetInput(inputImage);
    filter->SetAlpha(alpha);
    filter->SetBeta(beta);
    filter->SetRadius(radius);
    filter->Update();

    return filter->GetOutput();
  }

  /**
  \brief Get result of Image comparison between 2 images

  This runs itk::Testing::ComparisonImageFilter inside so the inputs are identical. Always updates the largest possible region.

  \param referenceImage The reference image for comparison
  \param checkImage The image to check
  \param differenceThreshold The minimum number of different pixels among both images; default is 0
  \param toleranceRadius The maximum distance to look for a matching pixel; default is 0
  \param numberOfPixelsTolerance The maximum maximum number of pixels that can be different; default is 10
  \param averageIntensityDifference The maximum allowed average intensity difference between both images; default is 0
  \return True if images are similar in accordance with passed parameters
  */
  template< class TImageType = ImageTypeFloat3D >
  bool GetResultOfImageComparasion(const typename TImageType::Pointer referenceImage, const typename TImageType::Pointer checkImage,
    const typename TImageType::PixelType differenceThreshold = 0, const unsigned int toleranceRadius = 0,
    const unsigned long long numberOfPixelsTolerance = 10, const typename TImageType::PixelType averageIntensityDifference = 0)
  {
    // check if sizes are different - that is a clear indicator that the images are NOT similar
    if (referenceImage->GetLargestPossibleRegion().GetSize() != checkImage->GetLargestPossibleRegion().GetSize())
    {
      return false;
    }

    // initialize the comparator
    auto diff = itk::Testing::ComparisonImageFilter< TImageType, TImageType >::New();
    diff->SetValidInput(referenceImage);
    diff->SetTestInput(checkImage);
    diff->SetDifferenceThreshold(differenceThreshold);
    diff->SetToleranceRadius(toleranceRadius);
    diff->UpdateLargestPossibleRegion();

    // check for different conditions 
    if (static_cast<typename TImageType::PixelType>(diff->GetTotalDifference()) > averageIntensityDifference)
    {
      // if there is an appreciable intensity difference between the images, check the number of difference pixels
      if (diff->GetNumberOfPixelsWithDifferences() > numberOfPixelsTolerance)
      {
        return false;
      }
    }

    return true;
  }

  /**
  \brief Perform the deformable registration

  \param movingImage The moving image for registration
  \param referenceImage The reference image for registration
  \param multiResLevels Number of multi-resolution levels for registration, defaults to 5
  \param iterationStart Start size of iteration for first multiResLevel, defaults to 10
  \param iterationStep Step size of the iterations increasing over each MultiResLevel, defaults to 10
  \param iterationEnd End size of iteration for first multiResLevel, defaults to 50
  \param regType The type of registration to perform, defaults to 'Demons'
  \param interpolatorType The type of interpolator to use, defaults to 'Linear'
  */
  /*template< class TImageType = ImageTypeFloat3D >
  typename TImageType::Pointer GetDeformRegisteredImage(const typename TImageType::Pointer movingImage, const typename TImageType::Pointer referenceImage,
  const unsigned int multiResLevels = 5,
  const unsigned int iterationStart = 10, const unsigned int iterationStep = 10, const unsigned int iterationEnd = 50,
  const int regType = Demons, const int interpolatorType = Linear)
  {
  // do histogram matchin of the 2 images
  auto movingImage_histoMatched = GetHistogramMatchedImage< TImageType >(movingImage, referenceImage, 40, 1024);

  // setup the displacement field
  using VectorPixelType = itk::Vector< float, TImageType::ImageDimension >;
  using DisplacementFieldType = itk::Image< VectorPixelType, TImageType::ImageDimension >;

  auto multiRes = typename itk::MultiResolutionPDEDeformableRegistration< TImageType, TImageType, DisplacementFieldType >::New();

  // set the registration type
  switch (regType)
  {
  case DiffeomorphicDemons:
  {
  auto filter = typename itk::DiffeomorphicDemonsRegistrationFilter< TImageType, TImageType, DisplacementFieldType >::New();
  filter->SetStandardDeviations(1.0);
  multiRes->SetRegistrationFilter(filter);
  break;
  }
  case SymmetricForcesDemons:
  {
  auto filter = typename itk::SymmetricForcesDemonsRegistrationFilter< TImageType, TImageType, DisplacementFieldType >::New();
  filter->SetStandardDeviations(1.0);
  multiRes->SetRegistrationFilter(filter);
  break;
  }
  case FastSymmetricForcesDemons:
  {
  auto filter = typename itk::FastSymmetricForcesDemonsRegistrationFilter< TImageType, TImageType, DisplacementFieldType >::New();
  filter->SetStandardDeviations(1.0);
  multiRes->SetRegistrationFilter(filter);
  break;
  }
  default: // does Demons
  {
  auto filter = typename itk::DemonsRegistrationFilter< TImageType, TImageType, DisplacementFieldType >::New();
  filter->SetStandardDeviations(1.0);
  multiRes->SetRegistrationFilter(filter);
  break;
  }
  }

  // set the different parameters of the multi resolution registration
  multiRes->SetNumberOfLevels(multiResLevels);
  multiRes->SetFixedImage(referenceImage);
  multiRes->SetMovingImage(movingImage_histoMatched);

  // set up the iterations
  std::vector< unsigned int > iterations_vector;
  for (size_t i = iterationStart; i <= iterationEnd; i += iterationStep)
  {
  iterations_vector.push_back(i);
  }
  unsigned int *iterations_array = &iterations_vector[0];
  multiRes->SetNumberOfIterations(iterations_array);

  // start the regisration
  try
  {
  multiRes->Update();
  }
  catch (itk::ExceptionObject & excp)
  {
  std::cerr << excp << std::endl;
  return referenceImage;
  }

  // warp the moving image
  auto warper = typename itk::WarpImageFilter< TImageType, TImageType, DisplacementFieldType >::New();
  warper->SetInput(movingImage);
  warper->SetOutputSpacing(referenceImage->GetSpacing());
  warper->SetOutputOrigin(referenceImage->GetOrigin());
  warper->SetOutputDirection(referenceImage->GetDirection());
  warper->SetDisplacementField(multiRes->GetOutput());

  // set up the interpolator type
  switch (interpolatorType)
  {
  case NearestNeighbor:
  {
  auto interpolator = typename itk::NearestNeighborInterpolateImageFunction< TImageType, double >::New();
  warper->SetInterpolator(interpolator);
  break;
  }
  case BSpline:
  {
  auto interpolator = typename itk::BSplineInterpolateImageFunction< TImageType, double >::New();
  warper->SetInterpolator(interpolator);
  break;
  }
  default: // linear by default
  {
  auto interpolator = typename itk::LinearInterpolateImageFunction< TImageType, double >::New();
  warper->SetInterpolator(interpolator);
  break;
  }
  }

  // perform the warping
  try
  {
  warper->Update();
  }
  catch (itk::ExceptionObject & excp)
  {
  std::cerr << excp << std::endl;
  return referenceImage;
  }

  return warper->GetOutput();
  }*/

  /**
  \brief Get the image orientation

  \param inputImage The input image
  \param desiredOrientation The desired orientation to conver the image to
  \return A pair of string (which represents the orientation) and an itk::Image which represents the inputImage in RAI form
  */
  template< class TImageType = ImageTypeFloat3D >
  std::pair< std::string, typename TImageType::Pointer > GetImageOrientation(const typename TImageType::Pointer inputImage, const std::string &desiredOrientation = "RAI")
  {
    if (TImageType::ImageDimension != 3)
    {
      std::cerr << "This function is only defined for 3D and 4D images.\n";
      exit(EXIT_FAILURE);
    }
    auto orientFilter = itk::OrientImageFilter< TImageType, TImageType >::New();
    orientFilter->SetInput(inputImage);
    orientFilter->UseImageDirectionOn();
    orientFilter->SetDirectionTolerance(0);
    orientFilter->SetCoordinateTolerance(0);

    auto desiredOrientation_wrap = desiredOrientation;
    std::transform(desiredOrientation_wrap.begin(), desiredOrientation_wrap.end(), desiredOrientation_wrap.begin(), ::toupper);

    std::map< std::string, itk::SpatialOrientation::ValidCoordinateOrientationFlags > orientationMap;
    orientationMap["Axial"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI;
    orientationMap["Coronal"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RSA;
    orientationMap["Sagittal"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ASL;
    orientationMap["RIP"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIP;
    orientationMap["LIP"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LIP;
    orientationMap["RSP"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RSP;
    orientationMap["LSP"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LSP;
    orientationMap["RIA"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIA;
    orientationMap["LIA"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LIA;
    orientationMap["RSA"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RSA;
    orientationMap["LSA"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LSA;
    orientationMap["IRP"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IRP;
    orientationMap["ILP"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ILP;
    orientationMap["SRP"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SRP;
    orientationMap["SLP"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SLP;
    orientationMap["IRA"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IRA;
    orientationMap["ILA"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ILA;
    orientationMap["SRA"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SRA;
    orientationMap["SLA"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SLA;
    orientationMap["RPI"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RPI;
    orientationMap["LPI"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPI;
    orientationMap["RAI"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI;
    orientationMap["LAI"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAI;
    orientationMap["RPS"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RPS;
    orientationMap["LPS"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPS;
    orientationMap["RAS"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAS;
    orientationMap["LAS"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAS;
    orientationMap["PRI"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PRI;
    orientationMap["PLI"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PLI;
    orientationMap["ARI"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ARI;
    orientationMap["ALI"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ALI;
    orientationMap["PRS"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PRS;
    orientationMap["PLS"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PLS;
    orientationMap["ARS"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ARS;
    orientationMap["ALS"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ALS;
    orientationMap["IPR"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IPR;
    orientationMap["SPR"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SPR;
    orientationMap["IAR"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IAR;
    orientationMap["SAR"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SAR;
    orientationMap["IPL"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IPL;
    orientationMap["SPL"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SPL;
    orientationMap["IAL"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IAL;
    orientationMap["SAL"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SAL;
    orientationMap["PIR"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PIR;
    orientationMap["PSR"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PSR;
    orientationMap["AIR"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_AIR;
    orientationMap["ASR"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ASR;
    orientationMap["PIL"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PIL;
    orientationMap["PSL"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PSL;
    orientationMap["AIL"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_AIL;
    orientationMap["ASL"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ASL;

    // set the desired orientation and update
    orientFilter->SetDesiredCoordinateOrientation(orientationMap[desiredOrientation_wrap]);
    orientFilter->Update();
    auto outputImage = orientFilter->GetOutput();

    std::string returnString;

    for (auto it = orientationMap.begin(); it != orientationMap.end(); ++it)
    {
      if (it->second == orientFilter->GetGivenCoordinateOrientation())
      {
        returnString = it->first;
      }
    }
    if (returnString.empty())
    {
      returnString = "Unknown";
    }

    return std::make_pair(returnString, outputImage);
  }

  /**
  \brief Get skull stripped image

  Templated over InputImageType, AtlasImageType and AtlasLabelType

  \param inputImage The input image on which to run the skull stripping
  \param atlasImage The atlas image
  \param atlasLabelImage The atlas label
  */
  template< class TImageType = ImageTypeFloat3D, class TAtlasImageType = TImageType, class TAtlasLabelType = TImageType >
  typename TImageType::Pointer GetSkullStrippedImage(const typename TImageType::Pointer inputImage,
    const typename TAtlasImageType::Pointer atlasImage, const typename TAtlasLabelType::Pointer atlasLabelImage)
  {
    // skull stripping initialization
    auto skullStripper = itk::StripTsImageFilter< TImageType, TAtlasImageType, TAtlasLabelType>::New();
    skullStripper->SetInput(inputImage);
    skullStripper->SetAtlasImage(atlasImage);
    skullStripper->SetAtlasBrainMask(atlasLabelImage);

    // actually do the skull stripping
    try
    {
      skullStripper->Update();
    }
    catch (itk::ExceptionObject &exception)
    {
      std::cerr << "Exception caught: " << exception << "\n";
      return inputImage;
    }

    // apply the generated mask
    auto masker = itk::MaskImageFilter< TImageType, TAtlasLabelType, TImageType >::New();
    masker->SetInput(inputImage);
    masker->SetMaskImage(skullStripper->GetOutput());

    try
    {
      masker->Update();
    }
    catch (itk::ExceptionObject &exception)
    {
      std::cerr << "Exception caught: " << exception << "\n";
      return inputImage;
    }

    return masker->GetOutput();
  }

  /**
  \brief Get the distance between 2 indeces of an itk::Image
  */
  template< class TImageType = ImageTypeFloat3D >
  float GetDistanceBetweenIndeces(const typename TImageType::IndexType point1, const typename TImageType::IndexType point2)
  {
    float currentDist = 0.0;
    for (size_t i = 0; i < TImageType::ImageDimension; i++)
    {
      currentDist += powf(point1[i] - point2[i], 2);
    }
    currentDist = sqrtf(currentDist);

    return currentDist;
  }

  /**
  \brief Get the distance between 2 itk::P of an itk::Image
  */
  inline float GetDistanceBetweenIndeces(const float* point1, const float* point2)
  {
    float currentDist = 0.0;
    for (size_t i = 0; i < 3; i++)
    {
      currentDist += powf(point1[i] - point2[i], 2);
    }
    currentDist = sqrtf(currentDist);

    return currentDist;
  }

  /**
  \brief Get the maximum distance and corresponding coordinate from a seed point in a label map

  \param inputLabelMap The label map on which to do the calculation
  \param indexForComputation The index of the seed point from where to do the distance measurements
  \param realCoordinatesPassed Bool which denotes if indexForComputation is real (e.g. the tumorPoints used by GLISTR) or an image index

  \return The maximum distance and the corrensponding index
  */
  template< class TImageType = ImageTypeFloat3D >
  std::pair< float, typename TImageType::IndexType > GetMaxDistanceInLabelMap(const typename TImageType::Pointer inputLabelMap,
    const typename TImageType::IndexType indexForComputation,
    bool realCoordinateInput = false, bool realCoordinateOutput = false)
  {
    auto indexToUse = indexForComputation;
    if (realCoordinateInput)
    {
      for (size_t i = 0; i < TImageType::ImageDimension; i++)
      {
        // gets the index of the point in question
        indexToUse[i] = std::abs((indexToUse[i] * inputLabelMap->GetSpacing()[i]) + inputLabelMap->GetOrigin()[i]);
      }
    }
    // setup the connected component segmentation
    auto connectedComponentFilter = itk::ConnectedThresholdImageFilter< TImageType, TImageType >::New();
    connectedComponentFilter->SetInput(inputLabelMap);
    connectedComponentFilter->SetSeed(indexToUse);
    connectedComponentFilter->SetReplaceValue(1);
    // we only want the selected voxel value to be segmented and nothing else
    auto currentPixelVal = inputLabelMap->GetPixel(indexToUse);
    connectedComponentFilter->SetLower(currentPixelVal);
    connectedComponentFilter->SetUpper(currentPixelVal);
    connectedComponentFilter->Update();

    itk::ImageRegionConstIterator <TImageType> iterator(connectedComponentFilter->GetOutput(), connectedComponentFilter->GetOutput()->GetLargestPossibleRegion());

    float maxDist = 0;
    typename TImageType::IndexType index_maxDist;

    // iterate through the whole image and find maximum distance
    for (iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator)
    {
      if (iterator.Get() > 0)
      {
        auto currentIndex = iterator.GetIndex();
        float currentDist = 0.0;// TBD  gcc is unable to deduce sutable type. Please fix this -> GetDistanceBetweenIndeces(currentIndex, indexToUse);
        currentDist = GetDistanceBetweenIndeces<TImageType>(currentIndex, indexToUse);

        if (currentDist > maxDist)
        {
          maxDist = currentDist;
          index_maxDist = currentIndex;
        }
      }
    }

    if (realCoordinateOutput)
    {
      for (size_t i = 0; i < TImageType::ImageDimension; i++)
      {
        // gets the index of the point in question
        index_maxDist[i] = (index_maxDist[i] * inputLabelMap->GetSpacing()[i]) + inputLabelMap->GetOrigin()[i];
      }
    }

    return std::make_pair(maxDist, index_maxDist);
  }

  /**
  \brief Create an empty (optionally pass a value) ITK image based on an input image with same properties

  \param inputImage The image to base the output on
  \param value The value to populate the new image with; defaults to '0'
  */
  template< class TImageType = ImageTypeFloat3D >
  typename TImageType::Pointer CreateImage(const typename TImageType::Pointer inputImage, const typename TImageType::PixelType value = 0)
  {
    typename TImageType::Pointer new_image = TImageType::New();
    new_image->SetLargestPossibleRegion(inputImage->GetLargestPossibleRegion());
    new_image->SetRequestedRegion(inputImage->GetRequestedRegion());
    new_image->SetBufferedRegion(inputImage->GetBufferedRegion());
    new_image->SetDirection(inputImage->GetDirection());
    new_image->SetOrigin(inputImage->GetOrigin());
    new_image->SetSpacing(inputImage->GetSpacing());
    new_image->Allocate();
    new_image->FillBuffer(value);
    return new_image;
  }

  /**
  \brief Create an empty (optionally pass a value) ITK image based on an input image with same properties

  \param inputImage The image to base the output on
  \param oldValues Values separated by 'x'
  \param newValues Values separated by 'x'
  */
  template< class TImageType = ImageTypeFloat3D >
  typename TImageType::Pointer ChangeImageValues(const typename TImageType::Pointer inputImage, const std::string &oldValues, const std::string &newValues)
  {
    auto oldValues_split = cbica::stringSplit(oldValues, "x");
    auto newValues_split = cbica::stringSplit(newValues, "x");
    if (oldValues_split.size() != newValues_split.size())
    {
      std::cerr << "Change values needs the old and new values to be of same size, for example '-cv 1x2,2x3.\n";
      return nullptr;
    }

    itk::ImageRegionConstIterator< TImageType > iterator(inputImage, inputImage->GetBufferedRegion());
    auto outputImage = inputImage;
    outputImage->DisconnectPipeline();
    itk::ImageRegionIterator< TImageType > outputIterator(outputImage, outputImage->GetBufferedRegion());
    outputIterator.GoToBegin();
    for (iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator, ++outputIterator)
    {
      for (size_t i = 0; i < oldValues_split.size(); i++)
      {
        if (iterator.Get() == std::atof(oldValues_split[i].c_str()))
        {
          outputIterator.Set(std::atof(newValues_split[i].c_str()));
        }
      }
    }

    return outputImage;
  }

  /**
  \brief Get distances in world coordinates across axes for an image

  \param inputImage
  */
  template< typename TImageType = ImageTypeFloat3D >
  itk::Vector< float, TImageType::ImageDimension > GetDistances(const typename TImageType::Pointer inputImage)
  {
    itk::Vector< float, TImageType::ImageDimension > distances;
    itk::Point< float, TImageType::ImageDimension > start_worldCoordinates, end_worldCoordinates;

    typename TImageType::IndexType start_image, end_image;

    auto size = inputImage->GetBufferedRegion().GetSize();

    for (size_t i = 0; i < TImageType::ImageDimension; i++)
    {
      start_image[i] = 0;
      end_image[i] = size[i] - 1;
    }

    inputImage->TransformIndexToPhysicalPoint(start_image, start_worldCoordinates);
    inputImage->TransformIndexToPhysicalPoint(end_image, end_worldCoordinates);

    for (size_t i = 0; i < TImageType::ImageDimension; i++)
    {
      distances[i] = std::abs(end_worldCoordinates[i] - start_worldCoordinates[i]); // real world image span along each axis
    }

    return distances;
  }

  /**
  \brief Resample an image to an isotropic resolution using the specified output spacing vector

  This filter uses the example https://itk.org/Wiki/ITK/Examples/ImageProcessing/ResampleImageFilter as a base while processing time-stamped images as well
  \param inputImage The input image to process
  \param outputSpacing The output spacing, always isotropic
  \param interpolator The type of interpolator to use; can be Linear, BSpline or NearestNeighbor
  \return The resized image
  */
  template< class TImageType = ImageTypeFloat3D >
  typename TImageType::Pointer ResampleImage(const typename TImageType::Pointer inputImage, const typename TImageType::SpacingType outputSpacing,
    typename TImageType::SizeType outputSize, const std::string interpolator = "Linear")
  {
    std::string interpolator_wrap = interpolator;
    std::transform(interpolator_wrap.begin(), interpolator_wrap.end(), interpolator_wrap.begin(), ::tolower);

    auto resampler = itk::ResampleImageFilter< TImageType, TImageType >::New();
    resampler->SetInput(inputImage);
    resampler->SetSize(outputSize);
    resampler->SetOutputSpacing(outputSpacing);
    resampler->SetOutputOrigin(inputImage->GetOrigin());
    resampler->SetOutputDirection(inputImage->GetDirection());
    resampler->SetOutputStartIndex(inputImage->GetLargestPossibleRegion().GetIndex());
    resampler->SetTransform(itk::IdentityTransform< double, TImageType::ImageDimension >::New());
    if (interpolator_wrap == "bspline")
    {
      auto interpolatorFunc = itk::BSplineInterpolateImageFunction< TImageType, double >::New();
      resampler->SetInterpolator(interpolatorFunc);
    }
    else if (interpolator_wrap.find("bicubic") != std::string::npos)
    {
      auto interpolatorFunc = itk::BSplineInterpolateImageFunction< TImageType >::New();
      interpolatorFunc->SetSplineOrder(3);
      resampler->SetInterpolator(interpolatorFunc);
    }
    else if (interpolator_wrap.find("nearest") != std::string::npos)
    {
      // if nearest label is found, do something special
      if (interpolator_wrap.find("nearestlabel") != std::string::npos)
      {
        ///
        // this situation can only happen in case of a mask image
        // The label image assumed to be an image of shortsC
        typedef itk::Image<short, TImageType::ImageDimension> LabelImageType;
        /// histogram calculation from ITK -- for texture feature pipeline
        auto caster = itk::CastImageFilter< TImageType, LabelImageType >::New();
        caster->SetInput(inputImage); //original input binary mask, not the masked image
        caster->Update();

        auto moving = caster->GetOutput();

        // Scan the unique labels in the image
        std::set<short> label_set;
        short *labels = moving->GetBufferPointer();
        int n_pixels = moving->GetPixelContainer()->Size();

        // Get the list of unique pixels
        short last_pixel = 0;
        for (int j = 0; j < n_pixels; j++)
        {
          short pixel = labels[j];
          if (last_pixel != pixel)
          {
            label_set.insert(pixel);
            last_pixel = pixel;
            //if (label_set.size() > 1000)
            //  throw GreedyException("Label wise interpolation not supported for image %s "
            //    "which has over 1000 distinct labels", filename);
          }
        }

        // Turn this set into an array
        std::vector<short> label_array(label_set.begin(), label_set.end());

        // Create a N-way voting filter
        typedef NaryLabelVotingFunctor<TImageType, LabelImageType> VotingFunctor;
        VotingFunctor vf(label_array);

        typedef itk::NaryFunctorImageFilter<TImageType, LabelImageType, VotingFunctor> VotingFilter;
        typename VotingFilter::Pointer fltVoting = VotingFilter::New();
        fltVoting->SetFunctor(vf);

        // Create a mini-pipeline of streaming filters
        for (size_t j = 0; j < label_array.size(); j++)
        {
          // Set up a threshold filter for this label
          typedef itk::BinaryThresholdImageFilter<LabelImageType, TImageType> ThresholdFilterType;
          typename ThresholdFilterType::Pointer fltThreshold = ThresholdFilterType::New();
          fltThreshold->SetInput(moving);
          fltThreshold->SetLowerThreshold(label_array[j]);
          fltThreshold->SetUpperThreshold(label_array[j]);
          fltThreshold->SetInsideValue(1.0);
          fltThreshold->SetOutsideValue(0.0);

          // Set up a smoothing filter for this label
          typedef itk::SmoothingRecursiveGaussianImageFilter<TImageType, TImageType> SmootherType;
          typename SmootherType::Pointer fltSmooth = SmootherType::New();
          fltSmooth->SetInput(fltThreshold->GetOutput());

          // Work out the sigmas for the filter
          //if (r_param.images[i].interp.sigma.physical_units)
          //{
          //  fltSmooth->SetSigma(r_param.images[i].interp.sigma.sigma);
          //}
          //else
          {
            typename SmootherType::SigmaArrayType sigma_array;
            auto default_sigma = std::sqrt(3);
            for (unsigned int d = 0; d < TImageType::ImageDimension; d++)
              sigma_array[d] = /*r_param.images[i].interp.sigma.sigma*/default_sigma * moving->GetSpacing()[d];
            fltSmooth->SetSigmaArray(sigma_array);
          }

          //// TODO: we should really be coercing the output into a vector image to speed up interpolation!
          //typedef FastWarpCompositeImageFilter<TImageType, TImageType, VectorImageType> InterpFilter;
          //typename InterpFilter::Pointer fltInterp = InterpFilter::New();
          //fltInterp->SetMovingImage(fltSmooth->GetOutput());
          //fltInterp->SetDeformationField(warp);
          //fltInterp->SetUsePhysicalSpace(true);

          //fltInterp->Update();

          // Add to the voting filter
          fltVoting->SetInput(j, fltSmooth->GetOutput());
        }

        // TODO: test out streaming!
        // Run this big pipeline
        fltVoting->Update();

        auto caster_back = itk::CastImageFilter< LabelImageType, TImageType >::New();
        caster_back->SetInput(fltVoting->GetOutput()); //original input binary mask, not the masked image
        caster_back->Update();
        return caster_back->GetOutput();
        ///
      }
      else
      {
        auto interpolatorFunc = itk::NearestNeighborInterpolateImageFunction< TImageType, double >::New();
        resampler->SetInterpolator(interpolatorFunc);
      }
    }
    //else if (interpolator_wrap.find("window") != std::string::npos)
    //{
    //  constexpr unsigned int WindowRadius = 5; // pass as parameter
    //  auto interpolatorFunc = itk::WindowedSincInterpolateImageFunction< TImageType, 
    //    WindowRadius,  // pass as parameter
    //    itk::Function::HammingWindowFunction< WindowRadius >, // pass as parameter
    //    itk::ConstantBoundaryCondition< TImageType >, // pass as parameter
    //    double  >::New();
    //  resampler->SetInterpolator(interpolatorFunc);
    //}
    else
    {
      auto interpolatorFunc = itk::LinearInterpolateImageFunction< TImageType, double >::New();
      resampler->SetInterpolator(interpolatorFunc);
    }
    resampler->UpdateLargestPossibleRegion();

    return resampler->GetOutput();
  }

  /**
  \brief Resample an image to an isotropic resolution using the specified output spacing vector

  This filter uses the example https://itk.org/Wiki/ITK/Examples/ImageProcessing/ResampleImageFilter as a base while processing time-stamped images as well
  \param inputImage The input image to process
  \param outputSpacing The desired output spacing
  \param interpolator The type of interpolator to use; can be Linear, BSpline or NearestNeighbor
  \return The resized image
  */
  template< class TImageType = ImageTypeFloat3D >
  typename TImageType::Pointer ResampleImage(const typename TImageType::Pointer inputImage,
    const itk::Vector< double, TImageType::ImageDimension > outputSpacing,
    const std::string interpolator = "Linear")
  {
    auto outputSize = inputImage->GetLargestPossibleRegion().GetSize();
    auto inputSpacing = inputImage->GetSpacing();
    if (TImageType::ImageDimension != 4)
    {
      for (size_t i = 0; i < TImageType::ImageDimension; i++)
      {
        outputSize[i] = std::round(outputSize[i] * inputSpacing[i] / outputSpacing[i]);
      }
    }
    else // preserve all time points of a time series image
    {
      for (size_t i = 0; i < 3; i++)
      {
        outputSize[i] = std::round(outputSize[i] * inputSpacing[i] / outputSpacing[i]);
      }
    }

    return ResampleImage< TImageType >(inputImage, outputSpacing, outputSize, interpolator);

  }

  /**
  \brief Resample an image to an isotropic resolution using the specified output spacing

  This filter uses the example https://itk.org/Wiki/ITK/Examples/ImageProcessing/ResampleImageFilter as a base while processing time-stamped images as well
  \param inputImage The input image to process
  \param outputSpacing The output spacing, always isotropic
  \param interpolator The type of interpolator to use; can be Linear, BSpline or NearestNeighbor
  \return The resized image
  */
  template< class TImageType = ImageTypeFloat3D >
  typename TImageType::Pointer ResampleImage(const typename TImageType::Pointer inputImage,
    const float outputSpacing = 1.0, const std::string interpolator = "Linear")
  {
    auto outputSize = inputImage->GetLargestPossibleRegion().GetSize();
    auto inputSpacing = inputImage->GetSpacing();
    auto outputSpacingVector = inputSpacing;
    if (TImageType::ImageDimension != 4)
    {
      outputSpacingVector.Fill(outputSpacing);
      for (size_t i = 0; i < TImageType::ImageDimension; i++)
      {
        outputSize[i] = std::round(outputSize[i] * inputSpacing[i] / outputSpacing);
      }
    }
    else // preserve all time points of a time series image
    {
      for (size_t i = 0; i < 3; i++)
      {
        outputSpacingVector[i] = outputSpacing;
        outputSize[i] = std::round(outputSize[i] * inputSpacing[i] / outputSpacing);
      }
      outputSpacingVector[3] = inputSpacing[3];
    }

    return ResampleImage< TImageType >(inputImage, outputSpacingVector, outputSize, interpolator);

  }

  /**
  \brief Resize an input image by a factor (expressed as a percentage)

  This filter uses the example https://itk.org/Wiki/ITK/Examples/ImageProcessing/ResampleImageFilter as a base while processing time-stamped images as well
  \param inputImage The input image to process
  \param resizeFactor The resize factor; can be greater than 100 (which causes an expanded image to be written) but can never be less than zero
  \param interpolator The type of interpolator to use; can be Linear, BSpline, BiCubic or NearestNeighbor
  \return The resized image
  */
  template< class TImageType = ImageTypeFloat3D >
  typename TImageType::Pointer ResizeImage(const typename TImageType::Pointer inputImage, const size_t resizeFactor, const std::string &interpolator = "Linear")
  {
    const float factor = static_cast<float>(resizeFactor) / 100;
    auto outputSize = inputImage->GetLargestPossibleRegion().GetSize();
    auto outputSpacing = inputImage->GetSpacing();
    if (TImageType::ImageDimension != 4)
    {
      for (size_t i = 0; i < TImageType::ImageDimension; i++)
      {
        outputSize[i] = outputSize[i] * factor;
        outputSpacing[i] = outputSpacing[i] / factor;
      }
    }
    else // preserve all time points of a time series image
    {
      for (size_t i = 0; i < 3; i++)
      {
        outputSize[i] = outputSize[i] * factor;
        outputSpacing[i] = outputSpacing[i] / factor;
      }
    }

    return ResampleImage < TImageType >(inputImage, outputSpacing, outputSize, interpolator);

  }

  /**
  \brief Get the unique values in an image

  \param inputImage The input image
  \param sortResult Whether the output should be sorted in ascending order or not, defaults to true
  */
  template< class TImageType = ImageTypeFloat3D >
  std::vector< typename TImageType::PixelType > GetUniqueValuesInImage(typename TImageType::Pointer inputImage, 
    bool sortResult = true)
  {
    itk::ImageRegionConstIterator< TImageType > iterator(inputImage, inputImage->GetBufferedRegion());

    std::vector< typename TImageType::PixelType > uniqueValues;

    for (iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator)
    {
      auto currentValue = iterator.Get();
      if (std::find(uniqueValues.begin(), uniqueValues.end(), currentValue) == uniqueValues.end())
      {
        uniqueValues.push_back(currentValue);
      }
    }

    if (sortResult)
    {
      std::sort(uniqueValues.begin(), uniqueValues.end(), std::less< typename TImageType::PixelType >());
    }

    return uniqueValues;
  }

  /**
  \brief Get the unique value labels in an image

  \param inputImage The input image
  */
  template< class TImageType = ImageTypeFloat3D >
  std::map< int, typename TImageType::Pointer > GetUniqueLabelImagessFromImage(typename TImageType::Pointer inputImage)
  {
    /// fix this up
    itk::ImageRegionConstIterator< TImageType > iterator(inputImage, inputImage->GetBufferedRegion());

    std::vector< typename TImageType::PixelType > uniqueValues;
    std::map< int, typename TImageType::Pointer > returnImages;

    for (iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator)
    {
      auto currentValue = iterator.Get();

      // if a new value has been found, initialize a new image
      if (std::find(uniqueValues.begin(), uniqueValues.end(), currentValue) == uniqueValues.end())
      {
        uniqueValues.push_back(currentValue);
        returnImages[currentValue] = CreateImage< TImageType >(inputImage);
      }

      returnImages[currentValue]->SetPixel(iterator.GetIndex(), 1);
    }

    return returnImages;
  }

  /**
  \brief Get Non-zero indeces of image
  */
  template< class TImageType = ImageTypeFloat3D >
  std::vector< typename TImageType::IndexType > GetNonZeroIndeces(typename TImageType::Pointer inputImage)
  {
    std::vector< typename TImageType::IndexType > outputIndeces;
    itk::ImageRegionConstIterator< TImageType > iterator(inputImage, inputImage->GetBufferedRegion());

    for (iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator)
    {
      if (iterator.Get() != 0)
      {
        outputIndeces.push_back(iterator.GetIndex());
      }
    }
    return outputIndeces;
  }

  /**
  \brief Get hausdorff distance between 2 labels

  \param inputLabel_1 The first label file
  \param inputLabel_2 The second label file
  \param percentile The percentile value for hausdorff; defaults to 0.95
  \return Hausdorff distance
  */
  template < typename TImageType = ImageTypeFloat3D >
  float GetHausdorffDistance(const typename TImageType::Pointer input_1, const typename TImageType::Pointer input_2, float percentile = 0.95)
  {
    // sanity check for percentile
    if (percentile > 1)
    {
      percentile = percentile / 100.0;
    }
    auto filter = HausdorffDistanceImageToImageMetric< TImageType, TImageType >::New();
    filter->SetFixedImage(input_1);
    filter->SetMovingImage(input_2);
    filter->SetPercentile(percentile);

    return filter->GetValue();
  }

  /**
  \brief Get sensitivity and specificity between 2 labels

  \param inputLabel_1 The first label file
  \param inputLabel_2 The second label file
  \return Map of metrics
  */
  template < typename TImageType = ImageTypeFloat3D >
  std::map< std::string, float > GetSensitivityAndSpecificity(const typename TImageType::Pointer input_1, const typename TImageType::Pointer input_2)
  {
    std::map< std::string, float > returnStruct;

    itk::ImageRegionConstIterator< TImageType > inputIterator(input_1, input_1->GetBufferedRegion()),
      outputIterator(input_2, input_2->GetBufferedRegion());

    std::vector< float > inputVector_1, inputVector_2;
    // iterate through the entire input image and if the label value matches the input value,
    // put '1' in the corresponding location of the output
    for (inputIterator.GoToBegin(); !inputIterator.IsAtEnd(); ++inputIterator)
    {
      outputIterator.SetIndex(inputIterator.GetIndex());
      inputVector_1.push_back(inputIterator.Get());
      inputVector_2.push_back(outputIterator.Get());
    }

    auto temp_roc = cbica::ROC_Values(inputVector_1, inputVector_2);

    returnStruct["Sensitivity"] = temp_roc["Sensitivity"];
    returnStruct["Specificity"] = temp_roc["Specificity"];
    returnStruct["Accuracy"] = temp_roc["Accuracy"];
    returnStruct["Precision"] = temp_roc["Precision"];

    return returnStruct;
  }

  /**
  \brief Get the statistics between 2 labels

  \param inputLabel_1 The first label file
  \param inputLabel_2 The second label file
  \return Map of various statistics and corresponding values
  */
  template < typename TImageType = ImageTypeFloat3D >
  std::map< std::string, double > GetLabelStatistics(const typename TImageType::Pointer inputLabel_1,
    const typename TImageType::Pointer inputLabel_2)
  {
    std::map< std::string, double > returnMap;

    auto uniqueLabels = GetUniqueValuesInImage< TImageType >(inputLabel_1);
    auto uniqueLabelsRef = GetUniqueValuesInImage< TImageType >(inputLabel_2);

    // sanity check
    if (uniqueLabels.size() != uniqueLabelsRef.size())
    {
      std::cerr << "The number of unique labels in input and reference image are not consistent.\n";
      return returnMap;
    }
    else
    {
      for (size_t i = 0; i < uniqueLabels.size(); i++)
      {
        if (uniqueLabels[i] != uniqueLabelsRef[i])
        {
          std::cerr << "The label values in input and reference image are not consistent.\n";
          return returnMap;
        }
      }
    }

    auto similarityFilter = itk::LabelOverlapMeasuresImageFilter< TImageType >::New();

    similarityFilter->SetSourceImage(inputLabel_1);
    similarityFilter->SetTargetImage(inputLabel_2);
    similarityFilter->Update();

    returnMap["Overlap_Overall"] = similarityFilter->GetTotalOverlap();
    //std::cout << "=== Entire Masked Area ===\n";
    returnMap["Jaccard_Overall"] = similarityFilter->GetUnionOverlap();
    returnMap["Dice_Overall"] = similarityFilter->GetMeanOverlap();
    returnMap["VolumeSimilarity_Overall"] = similarityFilter->GetVolumeSimilarity();
    returnMap["FalseNegativeError_Overall"] = similarityFilter->GetFalseNegativeError();
    returnMap["FalsePositiveError_Overall"] = similarityFilter->GetFalsePositiveError();

    if (uniqueLabels.size() > 2) // basically if there is something more than 0 and 1
    {
      //std::cout << "=== Individual Labels ===\n";
      //std::cout << "Property,Value\n";
      for (size_t i = 0; i < uniqueLabels.size(); i++)
      {
        auto uniqueLabels_string = std::to_string(uniqueLabels[i]);
        returnMap["Overlap_Label" + uniqueLabels_string] = similarityFilter->GetTargetOverlap(uniqueLabels[i]);
        returnMap["Jaccard_Label" + uniqueLabels_string] = similarityFilter->GetUnionOverlap(uniqueLabels[i]);
        returnMap["Dice_Label" + uniqueLabels_string] = similarityFilter->GetMeanOverlap(uniqueLabels[i]);
        returnMap["VolumeSimilarity_Label" + uniqueLabels_string] = similarityFilter->GetVolumeSimilarity(uniqueLabels[i]);
        returnMap["FalseNegativeError_Label" + uniqueLabels_string] = similarityFilter->GetFalseNegativeError(uniqueLabels[i]);
        returnMap["FalsePositiveError_Label" + uniqueLabels_string] = similarityFilter->GetFalsePositiveError(uniqueLabels[i]);
      }
    }

    // overall stats
    {
      auto temp_roc = GetSensitivityAndSpecificity< TImageType >(inputLabel_1, inputLabel_2);

      for (const auto &metric : temp_roc)
      {
        returnMap[metric.first + "_Overall"] = metric.second;
      }

      returnMap["Hausdorff95_Overall"] = GetHausdorffDistance< TImageType >(inputLabel_1, inputLabel_2, 0.95);
    }

    auto inputLabelsImages_1 = GetUniqueLabelImagessFromImage< TImageType >(inputLabel_1);
    auto inputLabelsImages_2 = GetUniqueLabelImagessFromImage< TImageType >(inputLabel_2);
    
    for (const auto &label : inputLabelsImages_1)
    {
      if (label.first != 0) // we don't care about the background value
      {
        auto valueString = std::to_string(label.first);

        auto temp_roc = GetSensitivityAndSpecificity< TImageType >(label.second, inputLabelsImages_2[label.first]);

        for (const auto &metric : temp_roc)
        {
          returnMap[metric.first + "_" + valueString] = metric.second;
        }

        returnMap["Hausdorff95_" + valueString] = GetHausdorffDistance< TImageType >(label.second, inputLabelsImages_2[label.first], 0.95);
      }

    }    

    return returnMap;
  }

  /**
  \brief Get the BraTS statistics between 2 brain labels

  Requires the following labesl to be initialized in both masks, otherwise the estimate is given as '0': 1,2,4

  \param inputLabel_1 The first brain label file
  \param inputLabel_2 The second brain label file
  \return Map of various statistics and corresponding values
  */
  template < typename TImageType = ImageTypeFloat3D >
  std::map< std::string, std::map< std::string, double > > GetBraTSLabelStatistics(const typename TImageType::Pointer inputLabel_1,
    const typename TImageType::Pointer inputLabel_2)
  {
    // return variable
    std::map< std::string, // the label string
      std::map< std::string, // the metric
      double > > // the value
      returnMap;

    auto uniqueLabels = GetUniqueValuesInImage< TImageType >(inputLabel_1);
    auto uniqueLabelsRef = GetUniqueValuesInImage< TImageType >(inputLabel_2);

    ///// this is not needed as we need to generate stats for all the BraTS values regardless
    //// sanity check
    //if (uniqueLabels.size() != uniqueLabelsRef.size())
    //{
    //  std::cerr << "The number of unique labels in input and reference image are not consistent.\n";
    //  return returnMap;
    //}
    //else
    //{
    //  for (size_t i = 0; i < uniqueLabels.size(); i++)
    //  {
    //    if (uniqueLabels[i] != uniqueLabelsRef[i])
    //    {
    //      std::cerr << "The label values in input and reference image are not consistent.\n";
    //      return returnMap;
    //    }
    //  }
    //}

    // remove background
    uniqueLabels.erase(uniqueLabels.begin());
    uniqueLabelsRef.erase(uniqueLabelsRef.begin());

    const std::vector< int > bratsValues = { 1,2,4 };
    std::map< int, std::string > bratsLabels;
    bratsLabels[1] = "NET";
    bratsLabels[2] = "ED";
    bratsLabels[4] = "ET";
    std::vector< int > missingLabels, missingLabelsRef;

    // check brats labels and populate missing labels to 
    for (size_t i = 0; i < bratsValues.size(); i++)
    {
      if (std::find(uniqueLabels.begin(), uniqueLabels.end(), bratsValues[i]) == uniqueLabels.end())
      {
        std::cerr << "Missing BraTS label in Label_1: " << bratsValues[i] << "\n";
        missingLabels.push_back(bratsValues[i]);
      }
      if (std::find(uniqueLabelsRef.begin(), uniqueLabelsRef.end(), bratsValues[i]) == uniqueLabelsRef.end())
      {
        std::cerr << "Missing BraTS label in Label_2: " << bratsValues[i] << "\n";
        missingLabelsRef.push_back(bratsValues[i]);
      }
    }

    auto inputLabelsImages_1 = GetUniqueLabelImagessFromImage< TImageType >(inputLabel_1);
    auto inputLabelsImages_2 = GetUniqueLabelImagessFromImage< TImageType >(inputLabel_2);

    // erase background label
    inputLabelsImages_1.erase(inputLabelsImages_1.find(0));
    inputLabelsImages_2.erase(inputLabelsImages_2.find(0));

    // populate empty masks for the missing labels
    for (size_t i = 0; i < missingLabels.size(); i++)
    {
      inputLabelsImages_1[missingLabels[i]] = CreateImage< TImageType >(inputLabel_1);
    }
    for (size_t i = 0; i < missingLabelsRef.size(); i++)
    {
      inputLabelsImages_2[missingLabelsRef[i]] = CreateImage< TImageType >(inputLabel_2);
    }

    // put individual labels in a single structure for easy processing
    std::map< std::string, typename TImageType::Pointer > regionsToCompare_1, regionsToCompare_2;

    for (const auto &label : inputLabelsImages_1)
    {
      regionsToCompare_1[bratsLabels[label.first]] = label.second;
    }
    for (const auto &label : inputLabelsImages_2)
    {
      regionsToCompare_2[bratsLabels[label.first]] = label.second;
    }

    // combine NET+ET to get TC and TC+ED to get WT
    {
      auto adder_1 = itk::AddImageFilter< TImageType >::New();
      adder_1->SetInput1(regionsToCompare_1["NET"]);
      adder_1->SetInput2(regionsToCompare_1["ET"]);
      adder_1->Update();

      auto adder_2 = itk::AddImageFilter< TImageType >::New();
      adder_2->SetInput1(regionsToCompare_2["NET"]);
      adder_2->SetInput2(regionsToCompare_2["ET"]);
      adder_2->Update();

      regionsToCompare_1["TC"] = adder_1->GetOutput();
      regionsToCompare_1["TC"]->DisconnectPipeline();

      regionsToCompare_2["TC"] = adder_2->GetOutput();
      regionsToCompare_2["TC"]->DisconnectPipeline();


      auto adder_1_WT = itk::AddImageFilter< TImageType >::New();
      adder_1_WT->SetInput1(regionsToCompare_1["ED"]);
      adder_1_WT->SetInput2(regionsToCompare_1["TC"]);
      adder_1_WT->Update();

      auto adder_2_WT = itk::AddImageFilter< TImageType >::New();
      adder_2_WT->SetInput1(regionsToCompare_2["ED"]);
      adder_2_WT->SetInput2(regionsToCompare_2["TC"]);
      adder_2_WT->Update();

      regionsToCompare_1["WT"] = adder_1_WT->GetOutput();
      regionsToCompare_1["WT"]->DisconnectPipeline();

      regionsToCompare_2["WT"] = adder_2_WT->GetOutput();
      regionsToCompare_2["WT"]->DisconnectPipeline();
    }

    // iterate over all regions (including missing brats regions in either image) and populate statistics
    for (const auto &region : regionsToCompare_1)
    {
      auto labelString = region.first; // the label for stats
      if (!labelString.empty())
      {
        // get the images to compare up front
        auto imageToCompare_1 = region.second;
        auto imageToCompare_2 = regionsToCompare_2[labelString];

        // in case one of the labels is missing, just put something
        auto stats_1 = itk::StatisticsImageFilter< TImageType >::New();
        stats_1->SetInput(imageToCompare_1);
        stats_1->Update();
        auto max_1 = stats_1->GetMaximum();
        auto stats_2 = itk::StatisticsImageFilter< TImageType >::New();
        stats_2->SetInput(imageToCompare_2);
        stats_2->Update();
        auto max_2 = stats_2->GetMaximum();

        auto similarityFilter = itk::LabelOverlapMeasuresImageFilter< TImageType >::New();

        similarityFilter->SetSourceImage(imageToCompare_1);
        similarityFilter->SetTargetImage(imageToCompare_2);
        similarityFilter->Update();

        returnMap[labelString]["Overlap"] = similarityFilter->GetTotalOverlap();
        returnMap[labelString]["Jaccard"] = similarityFilter->GetUnionOverlap();
        returnMap[labelString]["Dice"] = similarityFilter->GetMeanOverlap();
        if (std::isinf(returnMap[labelString]["Dice"]))
        {
          // this happens in the case where there is a label missing in both the reference and input annotations
          returnMap[labelString]["Dice"] = 1;
        }
        returnMap[labelString]["VolumeSimilarity"] = similarityFilter->GetVolumeSimilarity();
        returnMap[labelString]["FalseNegativeError"] = similarityFilter->GetFalseNegativeError();
        returnMap[labelString]["FalsePositiveError"] = similarityFilter->GetFalsePositiveError();

        auto temp_roc = GetSensitivityAndSpecificity< TImageType >(imageToCompare_1, imageToCompare_2);

        for (const auto &metric : temp_roc)
        {
          returnMap[labelString][metric.first] = metric.second;
        }

        if ((max_1 == 0) && (max_2 == 0))
        {
          returnMap[labelString]["Sensitivity"] = 1;
          returnMap[labelString]["Specificity"] = 1;
        }
        if (std::isnan(returnMap[labelString]["Sensitivity"]))
        {
          if (max_1 != max_2)
          {
            returnMap[labelString]["Sensitivity"] = 0;
          }
          else
          {
            returnMap[labelString]["Sensitivity"] = 1;
          }
        }
        if (std::isinf(returnMap[labelString]["Sensitivity"]))
        {
          returnMap[labelString]["Sensitivity"] = 1;
        }

        /// not used till implementation gets standardized
        //returnMap[labelString]["Hausdorff95"] = GetHausdorffDistance< TImageType >(imageToCompare_1, imageToCompare_2, 0.95);
        //returnMap[labelString]["Hausdorff99"] = GetHausdorffDistance< TImageType >(imageToCompare_1, imageToCompare_2, 0.99);
        bool hausdorffFound = true;
        std::string hausdorffExe = cbica::getExecutablePath() + "/Hausdorff95"
#if WIN32
          + ".exe"
#endif
          ;
        if (!cbica::isFile(hausdorffExe))
        {
          hausdorffExe = cbica::getExecutablePath() + "../hausdorff95/Hausdorff95"
#if WIN32
            + ".exe"
#endif
            ;
          if (!cbica::isFile(hausdorffExe))
          {
            std::cerr << "Could not find Hausdorff95 executable, so not computing this metric.\n";
            hausdorffFound = false;
          }
        }

        if (hausdorffFound)
        {
          if ((max_1 == 0) || (max_2 == 0))
          {
            // this is the case where one of the labels is missing
            returnMap[labelString]["Hausdorff95"] = NAN;
          }
          else
          {
            auto tempDir = cbica::createTmpDir();
            auto file_1 = tempDir + "/mask_1.nii.gz";
            auto file_2 = tempDir + "/mask_2.nii.gz";
            auto writer = itk::ImageFileWriter< TImageType >::New();
            writer->SetInput(imageToCompare_1);
            writer->SetFileName(file_1);
            try
            {
              writer->Write();
            }
            catch (itk::ExceptionObject &e)
            {
              std::cerr << "Error occurred while trying to write the image '" << file_1 << "': " << e.what() << "\n";
            }
            writer->SetInput(imageToCompare_2);
            writer->SetFileName(file_2);
            try
            {
              writer->Write();
            }
            catch (itk::ExceptionObject &e)
            {
              std::cerr << "Error occurred while trying to write the image '" << file_2 << "': " << e.what() << "\n";
            }
            std::array< char, 128 > buffer;
            std::string result;
            FILE *pPipe;
#if WIN32
#define POPEN _popen
#define PCLOSE _pclose
#else
#define POPEN popen
#define PCLOSE pclose
#endif
            pPipe = POPEN((hausdorffExe + " -gt " + file_1 + " -m " + file_2).c_str(), "r");
            if (!pPipe)
            {
              std::cerr << "Couldn't start command.\n";
            }
            while (fgets(buffer.data(), 128, pPipe) != NULL)
            {
              result += buffer.data();
            }
            auto returnCode = PCLOSE(pPipe);
            // remove "\n"
            result.pop_back();
            result.pop_back();
            returnMap[labelString]["Hausdorff95"] = std::atof(result.c_str());
            cbica::removeDirectoryRecursively(tempDir, true);
          }
          // in case a label is not defined, use the longest diagonal
          if (std::isnan(returnMap[labelString]["Hausdorff95"]) || std::isinf(returnMap[labelString]["Hausdorff95"]))
          {
            // correct prediction for missing label
            if ((max_1 == 0) && (max_2 == 0))
            {
              returnMap[labelString]["Hausdorff95"] = 0;
            }
            else
            {
              auto size = imageToCompare_1->GetLargestPossibleRegion().GetSize();
              auto diag_plane_squared = std::pow(size[0], 2) + std::pow(size[1], 2);
              auto diag_cube = std::sqrt(std::pow(size[2], 2) + diag_plane_squared);
              returnMap[labelString]["Hausdorff95"] = diag_cube;
            }
          }
        } // end hausdorff found
      } // end labelString check
    }

    return returnMap;
  }

}
