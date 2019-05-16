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

//#include "itkMultiResolutionPDEDeformableRegistration.h"
//#include "itkDemonsRegistrationFilter.h"
//#include "itkDiffeomorphicDemonsRegistrationFilter.h"
//#include "itkFastSymmetricForcesDemonsRegistrationFilter.h"
//#include "itkSymmetricForcesDemonsRegistrationFilter.h"

#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
//#include "itkWindowedSincInterpolateImageFunction.h"

#include "itkTestingComparisonImageFilter.h"

#include "itkStripTsImageFilter.h"
#include "itkMaskImageFilter.h"

#include "cbicaUtilities.h"
#include "cbicaITKImageInfo.h"

#include "gdcmMD5.h"
#include "gdcmReader.h"

#include "DicomIOManager.h"

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

/*
\namespace cbica
\brief Namespace for differentiating functions written for internal use
*/
namespace cbica
{
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
    auto diff = typename itk::Testing::ComparisonImageFilter< TImageType, TImageType >::New();
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
  \brief Check properties of 2 images to see if they are defined in the same space.
  */
  template< typename TImageType >
  inline bool ImageSanityCheck(const typename TImageType::Pointer image1, const typename TImageType::Pointer image2)
  {
    auto size_1 = image1->GetLargestPossibleRegion().GetSize();
    auto size_2 = image2->GetLargestPossibleRegion().GetSize();

    auto origin_1 = image1->GetOrigin();
    auto origin_2 = image2->GetOrigin();

    auto spacing_1 = image1->GetSpacing();
    auto spacing_2 = image2->GetSpacing();

    for (size_t i = 0; i < TImageType::ImageDimension; i++)
    {
      if (size_1[i] != size_2[i])
      {
        std::cerr << "Size mismatch at dimension '" << i << "'\n";
        return false;
      }
      if (origin_1[i] != origin_2[i])
      {
        std::cerr << "Origin mismatch at dimension '" << i << "'\n";
        return false;
      }
      if (spacing_1[i] != spacing_2[i])
      {
        std::cerr << "Spacing mismatch at dimension '" << i << "'\n";
        return false;
      }
    }

    return true;
  }

  /**
  \brief Check properties of 2 images to see if they are defined in the same space.

  Checks are done based on cbica::ImageInfo class
  */
  inline bool ImageSanityCheck(const std::string &image1, const std::string &image2)
  {
    auto imageInfo1 = cbica::ImageInfo(image1);
    auto imageInfo2 = cbica::ImageInfo(image2);

    if (imageInfo1.GetImageDimensions() != imageInfo2.GetImageDimensions())
    {
      std::cout << "The dimensions of the image_1 (" << image1 << ") and image_2 (" << image2 << ") doesn't match.\n";
      return false;
    }

    // check size, spacing and origin information as well
    auto dims = imageInfo1.GetImageDimensions();

    auto imageSize1 = imageInfo1.GetImageSize();
    auto imageSize2 = imageInfo2.GetImageSize();

    auto imageSpacing1 = imageInfo1.GetImageSpacings();
    auto imageSpacing2 = imageInfo2.GetImageSpacings();

    auto imageOrigin1 = imageInfo1.GetImageOrigins();
    auto imageOrigin2 = imageInfo2.GetImageOrigins();

    for (size_t d = 0; d < dims; d++)
    {
      if (imageSize1[d] != imageSize2[d])
      {
        std::cout << "The size in dimension[" << d << "] of the image_1 (" << image1 << ") and image_2 (" << image2 << ") doesn't match.\n";
        return false;
      }
      if (imageSpacing1[d] != imageSpacing2[d])
      {
        std::cout << "The spacing in dimension[" << d << "] of the image_1 (" << image1 << ") and image_2 (" << image2 << ") doesn't match.\n";
        return false;
      }
      if (imageOrigin1[d] != imageOrigin2[d])
      {
        std::cout << "The origin in dimension[" << d << "] of the image_1 (" << image1 << ") and image_2 (" << image2 << ") doesn't match.\n";
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
  \return A pair of string (which represents the orientation) and an itk::Image which represents the inputImage in RAI form
  */
  template< class TImageType = ImageTypeFloat3D >
  std::pair< std::string, typename TImageType::Pointer > GetImageOrientation(const typename TImageType::Pointer inputImage)
  {
    using namespace itk::SpatialOrientation;
    auto orientFilter = itk::OrientImageFilter< TImageType, TImageType >::New();
    orientFilter->SetInput(inputImage);
    orientFilter->SetDesiredCoordinateOrientation(ITK_COORDINATE_ORIENTATION_RAI);
    orientFilter->Update();

    std::string returnString;

    switch (orientFilter->GetGivenCoordinateOrientation())
    {
    case ITK_COORDINATE_ORIENTATION_RIP:
    {
      returnString = "RIP";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_LIP:
    {
      returnString = "LIP";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_RSP:
    {
      returnString = "RSP";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_LSP:
    {
      returnString = "LSP";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_RIA:
    {
      returnString = "RIA";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_LIA:
    {
      returnString = "LIA";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_LSA:
    {
      returnString = "LSA";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_IRP:
    {
      returnString = "IRP";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_ILP:
    {
      returnString = "ILP";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_SRP:
    {
      returnString = "SRP";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_SLP:
    {
      returnString = "SLP";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_IRA:
    {
      returnString = "IRA";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_ILA:
    {
      returnString = "ILA";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_SRA:
    {
      returnString = "SRA";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_SLA:
    {
      returnString = "SLA";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_RPI:
    {
      returnString = "RPI";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_LPI:
    {
      returnString = "LPI";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_RAI:
    {
      returnString = "RAI";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_LAI:
    {
      returnString = "LAI";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_RPS:
    {
      returnString = "RPS";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_LPS:
    {
      returnString = "LPS";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_RAS:
    {
      returnString = "RAS";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_LAS:
    {
      returnString = "LAS";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_PRI:
    {
      returnString = "PRI";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_PLI:
    {
      returnString = "PLI";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_ARI:
    {
      returnString = "ARI";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_ALI:
    {
      returnString = "ALI";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_PRS:
    {
      returnString = "PRS";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_PLS:
    {
      returnString = "PLS";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_ARS:
    {
      returnString = "ARS";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_ALS:
    {
      returnString = "ALS";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_IPR:
    {
      returnString = "IPR";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_SPR:
    {
      returnString = "SPR";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_IAR:
    {
      returnString = "IAR";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_SAR:
    {
      returnString = "SAR";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_IPL:
    {
      returnString = "IPL";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_SPL:
    {
      returnString = "SPL";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_IAL:
    {
      returnString = "IAL";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_SAL:
    {
      returnString = "SAL";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_PIR:
    {
      returnString = "PIR";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_PSR:
    {
      returnString = "PSR";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_AIR:
    {
      returnString = "AIR";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_ASR:
    {
      returnString = "ASR";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_PIL:
    {
      returnString = "PIL";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_PSL:
    {
      returnString = "PSL";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_AIL:
    {
      returnString = "AIL";
      break;
    }
    case ITK_COORDINATE_ORIENTATION_ASL:
    {
      returnString = "ASL";
      break;
    }
    default:
    {
      returnString = "UNKNOWN";
      break;
    }
    }

    return std::make_pair(returnString, orientFilter->GetOutput());
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
  typename TImageType::Pointer ResampleImage(const typename TImageType::Pointer inputImage, const itk::Vector< double, TImageType::ImageDimension > &outputSpacing, const std::string interpolator = "Linear")
  {
    auto outputSize = inputImage->GetLargestPossibleRegion().GetSize();
    auto outputSpacingVector = outputSpacing;
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

    return ResampleImage< TImageType >(inputImage, outputSpacingVector, outputSize, interpolator);

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
      auto interpolatorFunc = itk::NearestNeighborInterpolateImageFunction< TImageType, double >::New();
      resampler->SetInterpolator(interpolatorFunc);
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
  \brief Resample an image to an isotropic resolution using the specified output spacing

  This filter uses the example https://itk.org/Wiki/ITK/Examples/ImageProcessing/ResampleImageFilter as a base while processing time-stamped images as well
  \param inputImage The input image to process
  \param outputSpacing The output spacing, always isotropic
  \param interpolator The type of interpolator to use; can be Linear, BSpline or NearestNeighbor
  \return The resized image
  */
  template< class TImageType = ImageTypeFloat3D >
  typename TImageType::Pointer ResampleImage(const typename TImageType::Pointer inputImage, const float outputSpacing = 1.0, const std::string interpolator = "Linear")
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
  std::vector< typename TImageType::PixelType > GetUniqueValuesInImage(typename TImageType::Pointer inputImage, bool sortResult = true)
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

}
