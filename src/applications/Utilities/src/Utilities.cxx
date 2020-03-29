#include "cbicaCmdParser.h"
#include "cbicaLogging.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"
#include "DicomMetadataReader.h"

#include "itkBoundingBox.h"
#include "itkPointSet.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkJoinSeriesImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkLabelOverlapMeasuresImageFilter.h"
#include "itkHausdorffDistanceImageFilter.h"

#include "cbicaDCMQIWrapper.h"

//! Detail the available algorithms to make it easier to initialize
enum AvailableAlgorithms
{
  None,
  Resize,
  Resample,
  SanityCheck,
  Information,
  Casting,
  UniqueValues,
  TestComparison,
  BoundingBox,
  CreateMask,
  ChangeValue,
  DicomLoadTesting,
  Dicom2Nifti,
  Nifti2Dicom,
  Nifti2DicomSeg,
  OrientImage,
  ThresholdAbove,
  ThresholdBelow,
  ThresholdAboveAndBelow,
  ThresholdOtsu,
  ThresholdBinary,
  ConvertFormat,
  Image2World,
  World2Image,
  ImageStack2Join,
  JoinedImage2Stack,
  LabelSimilarity,
  Hausdorff
};

int requestedAlgorithm = 0;

std::string inputImageFile, inputMaskFile, outputImageFile, targetImageFile, resamplingInterpolator = "LINEAR", dicomSegJSON, orientationDesired, coordinateToTransform, referenceMaskForSimilarity;
std::string dicomFolderPath;
size_t resize = 100;
int testRadius = 0, testNumber = 0;
float testThresh = 0.0, testAvgDiff = 0.0, lowerThreshold = 1, upperThreshold = std::numeric_limits<float>::max();
std::string changeOldValues, changeNewValues, resamplingResolution_full = "1.0,1.0,1.0", resamplingReference;
float resamplingResolution = 1.0, thresholdAbove = 0.0, thresholdBelow = 0.0, thresholdOutsideValue = 0.0, thresholdInsideValue = 1.0;
float imageStack2JoinSpacing = 1.0;
int joinedImage2stackedAxis;

bool uniqueValsSort = true, boundingBoxIsotropic = true;

template< class TImageType >
int algorithmsRunner()
{
  if (requestedAlgorithm == Resize)
  {
    auto outputImage = cbica::ResizeImage< TImageType >(cbica::ReadImage< TImageType >(inputImageFile), resize, resamplingInterpolator);
    cbica::WriteImage< TImageType >(outputImage, outputImageFile);

    std::cout << "Resizing by a factor of " << resize << "% completed.\n";
  }
  else if (requestedAlgorithm == Resample)
  {
    auto inputImage = cbica::ReadImage< TImageType >(inputImageFile);
    itk::Vector< double, TImageType::ImageDimension > outputSpacing = inputImage->GetSpacing();
    if (!resamplingReference.empty())
    {
      auto imageInfo_1 = cbica::ImageInfo(inputImageFile);
      auto imageInfo_2 = cbica::ImageInfo(resamplingReference);
      if (imageInfo_1.GetImageDimensions() != imageInfo_2.GetImageDimensions())
      {
        std::cerr << "The resampling reference image needs to be the same dimension as the input image.\n";
        return EXIT_FAILURE;
      }
      outputSpacing = cbica::ReadImage< TImageType >(resamplingReference)->GetSpacing();
    }
    else
    {
      auto resolution_split = cbica::stringSplit(resamplingResolution_full, ",");
      if (resolution_split.size() == 1)
      {
        std::cout << "Isotropic resultion of '" << resolution_split[0] << "' has been selected.\n";
        for (size_t d = 1; d < TImageType::ImageDimension; d++)
        {
          resolution_split.push_back(resolution_split[0]);
        }
      }
      if (resolution_split.size() != TImageType::ImageDimension)
      {
        std::cerr << "The resampling resolution needs to be the same dimension as the input image.\n";
        return EXIT_FAILURE;
      }
      for (size_t d = 0; d < TImageType::ImageDimension; d++)
      {
        outputSpacing[d] = std::atof(resolution_split[d].c_str());
      }
    }
    auto outputImage = cbica::ResampleImage< TImageType >(inputImage, outputSpacing, resamplingInterpolator);
    cbica::WriteImage< TImageType >(outputImage, outputImageFile);

    std::cout << "Resampled image to a resolution of '" << outputSpacing << "' using interpolator '" << resamplingInterpolator << "'.\n";
  }
  else if (requestedAlgorithm == UniqueValues)
  {
    bool sort = true;
    if (uniqueValsSort == 0)
    {
      sort = false;
    }
    auto uniqueValues = cbica::GetUniqueValuesInImage< TImageType >(cbica::ReadImage < TImageType >(inputImageFile), sort);

    if (!uniqueValues.empty())
    {
      std::cout << "Unique values:\n";
      for (size_t i = 0; i < uniqueValues.size(); i++)
      {
        std::cout << cbica::to_string_precision(uniqueValues[i]) << "\n";
      }
    }
  }
  else if (requestedAlgorithm == CreateMask)
  {
    auto thresholder = itk::BinaryThresholdImageFilter< TImageType, TImageType >::New();
    thresholder->SetInput(cbica::ReadImage< TImageType >(inputImageFile));
    thresholder->SetLowerThreshold(lowerThreshold);
    thresholder->SetUpperThreshold(upperThreshold);
    thresholder->SetOutsideValue(0);
    thresholder->SetInsideValue(1);
    thresholder->Update();

    cbica::WriteImage< TImageType >(thresholder->GetOutput(), outputImageFile);
    std::cout << "Create Mask completed.\n";
  }
  else if (requestedAlgorithm == DicomLoadTesting)
  {
    auto readDicomImage = cbica::ReadImage< TImageType >(dicomFolderPath);
    if (!readDicomImage)
    {
      std::cout << "Dicom Load Failed" << std::endl;
      return EXIT_FAILURE;
    }
    auto inputNiftiImage = cbica::ReadImage< TImageType >(inputImageFile);

    //cbica::WriteImage< TImageType>(inputNiftiImage, "readNifti.nii.gz");

    bool differenceFailed = false;

    auto diffFilter = itk::Testing::ComparisonImageFilter< TImageType, TImageType >::New();
    diffFilter->SetValidInput(inputNiftiImage);
    diffFilter->SetTestInput(readDicomImage);
    auto size = inputNiftiImage->GetBufferedRegion().GetSize();
    diffFilter->VerifyInputInformationOn();
    diffFilter->SetDifferenceThreshold(testThresh);
    diffFilter->SetToleranceRadius(testRadius);
    diffFilter->UpdateLargestPossibleRegion();

    size_t totalSize = 1;
    for (size_t d = 0; d < TImageType::ImageDimension; d++)
    {
      totalSize *= size[d];
    }

    const double averageIntensityDifference = diffFilter->GetTotalDifference();
    const unsigned long numberOfPixelsWithDifferences = diffFilter->GetNumberOfPixelsWithDifferences();

    std::cout << "Total Voxels in Image       : " << totalSize << "\n";
    std::cout << "Number of Difference Voxels : " << diffFilter->GetNumberOfPixelsWithDifferences() << "\n";
    std::cout << "Percentage of Diff Voxels   : " << diffFilter->GetNumberOfPixelsWithDifferences() * 100 / totalSize << "\n";
    std::cout << "Minimum Intensity Difference: " << diffFilter->GetMinimumDifference() << "\n";
    std::cout << "Maximum Intensity Difference: " << diffFilter->GetMaximumDifference() << "\n";
    std::cout << "Average Intensity Difference: " << diffFilter->GetMeanDifference() << "\n";
    std::cout << "Overall Intensity Difference: " << diffFilter->GetTotalDifference() << "\n";

    int numberOfPixelsTolerance = 70;
    if (averageIntensityDifference > 0.0)
    {
      if (static_cast<int>(numberOfPixelsWithDifferences) >
        numberOfPixelsTolerance)
      {
        differenceFailed = true;
      }
      else
      {
        differenceFailed = false;
      }
    }
    else
    {
      differenceFailed = false;
    }
    return differenceFailed;
    //return EXIT_FAILURE;
  }
  else if (requestedAlgorithm == Dicom2Nifti)
  {
    auto readDicomImage = cbica::ReadImage< TImageType >(inputImageFile);
    if (!readDicomImage)
    {
      std::cout << "Dicom Load Failed" << std::endl;
      return EXIT_FAILURE;
    }

    cbica::WriteImage< TImageType>(readDicomImage, outputImageFile);

    if (!targetImageFile.empty())
    {
      if (!cbica::ImageSanityCheck< TImageType >(readDicomImage, cbica::ReadImage< TImageType >(targetImageFile)))
      {
        std::cerr << "Input image and target image physical space mismatch.\n";
        return EXIT_FAILURE;
      }
      else
      {
        auto diffFilter = itk::Testing::ComparisonImageFilter< TImageType, TImageType >::New();
        auto inputImage = cbica::ReadImage< TImageType >(inputImageFile);
        auto size = inputImage->GetBufferedRegion().GetSize();
        diffFilter->SetValidInput(cbica::ReadImage< TImageType >(targetImageFile));
        diffFilter->SetTestInput(inputImage);
        diffFilter->VerifyInputInformationOn();
        diffFilter->SetDifferenceThreshold(testThresh);
        diffFilter->SetToleranceRadius(testRadius);
        diffFilter->UpdateLargestPossibleRegion();

        size_t totalSize = 1;
        for (size_t d = 0; d < TImageType::ImageDimension; d++)
        {
          totalSize *= size[d];
        }

        std::cout << "Total Voxels/Pixels in Image: " << totalSize << "\n";
        std::cout << "Number of Difference Voxels : " << diffFilter->GetNumberOfPixelsWithDifferences() << "\n";
        std::cout << "Percentage of Diff Voxels   : " << diffFilter->GetNumberOfPixelsWithDifferences() * 100 / totalSize << "\n";
        std::cout << "Minimum Intensity Difference: " << diffFilter->GetMinimumDifference() << "\n";
        std::cout << "Maximum Intensity Difference: " << diffFilter->GetMaximumDifference() << "\n";
        std::cout << "Average Intensity Difference: " << diffFilter->GetMeanDifference() << "\n";
        std::cout << "Overall Intensity Difference: " << diffFilter->GetTotalDifference() << "\n";
      }
    }
  }
  else if (requestedAlgorithm == Nifti2Dicom)
  {
    auto referenceDicom = targetImageFile;
    cbica::WriteDicomImageFromReference< TImageType >(referenceDicom, cbica::ReadImage< TImageType >(inputImageFile), outputImageFile);
    if (cbica::exists(outputImageFile))
    {
      std::cout << "Finished writing the DICOM file.\n";
    }
  }
  else if (requestedAlgorithm == Nifti2DicomSeg)
  {
    if (TImageType::ImageDimension != 3)
    {
      std::cerr << "NIfTI to DICOM-Seg conversion can only be done for 3D images.\n";
      return EXIT_SUCCESS;
    }
    auto actualDicomReferenceDir = targetImageFile;

    if (!cbica::isDir(targetImageFile)) // the reference dicom image series
    {
      if (!cbica::IsDicom(targetImageFile))
      {
        std::cerr << "The supplied reference DICOM file, '" << targetImageFile << "' was not detected as a DICOM file. Please try again.\n";
      }
      else
      {
        actualDicomReferenceDir = cbica::getFilenamePath(targetImageFile);
      }
      std::cerr << "The supplied reference DICOM directory, '" << targetImageFile << "' was not detected as a directory. Please try again.\n";
    }
    cbica::ConvertNiftiToDicomSeg(inputImageFile, actualDicomReferenceDir, dicomSegJSON, outputImageFile);
  }
  else if (requestedAlgorithm == ChangeValue)
  {
    auto outputImage = cbica::ChangeImageValues< TImageType >(cbica::ReadImage< TImageType >(inputImageFile), changeOldValues, changeNewValues);

    if (!outputImage.IsNull())
    {
      cbica::WriteImage< TImageType >(outputImage, outputImageFile);
      std::cout << "Create Mask completed.\n";
      return EXIT_SUCCESS;
    }
    else
    {
      return EXIT_FAILURE;
    }

  }
  else if (requestedAlgorithm == Casting)
  {
    if (targetImageFile == "uchar")
    {
      using DefaultPixelType = unsigned char;
      using CurrentImageType = itk::Image< DefaultPixelType, TImageType::ImageDimension >;
      cbica::WriteImage< CurrentImageType >(cbica::ReadImage< CurrentImageType >(inputImageFile), outputImageFile);
    }
    else if (targetImageFile == "char")
    {
      using DefaultPixelType = char;
      using CurrentImageType = itk::Image< DefaultPixelType, TImageType::ImageDimension >;
      cbica::WriteImage< CurrentImageType >(cbica::ReadImage< CurrentImageType >(inputImageFile), outputImageFile);
    }
    else if (targetImageFile == "ushort")
    {
      using DefaultPixelType = unsigned short;
      using CurrentImageType = itk::Image< DefaultPixelType, TImageType::ImageDimension >;
      cbica::WriteImage< CurrentImageType >(cbica::ReadImage< CurrentImageType >(inputImageFile), outputImageFile);
    }
    else if (targetImageFile == "short")
    {
      using DefaultPixelType = short;
      using CurrentImageType = itk::Image< DefaultPixelType, TImageType::ImageDimension >;
      cbica::WriteImage< CurrentImageType >(cbica::ReadImage< CurrentImageType >(inputImageFile), outputImageFile);
    }
    else if (targetImageFile == "uint")
    {
      using DefaultPixelType = unsigned int;
      using CurrentImageType = itk::Image< DefaultPixelType, TImageType::ImageDimension >;
      cbica::WriteImage< CurrentImageType >(cbica::ReadImage< CurrentImageType >(inputImageFile), outputImageFile);
    }
    else if (targetImageFile == "int")
    {
      using DefaultPixelType = int;
      using CurrentImageType = itk::Image< DefaultPixelType, TImageType::ImageDimension >;
      cbica::WriteImage< CurrentImageType >(cbica::ReadImage< CurrentImageType >(inputImageFile), outputImageFile);
    }
    else if (targetImageFile == "ulong")
    {
      using DefaultPixelType = unsigned long;
      using CurrentImageType = itk::Image< DefaultPixelType, TImageType::ImageDimension >;
      cbica::WriteImage< CurrentImageType >(cbica::ReadImage< CurrentImageType >(inputImageFile), outputImageFile);
    }
    else if (targetImageFile == "long")
    {
      using DefaultPixelType = long;
      using CurrentImageType = itk::Image< DefaultPixelType, TImageType::ImageDimension >;
      cbica::WriteImage< CurrentImageType >(cbica::ReadImage< CurrentImageType >(inputImageFile), outputImageFile);
    }
    else if (targetImageFile == "ulonglong")
    {
      using DefaultPixelType = unsigned long long;
      using CurrentImageType = itk::Image< DefaultPixelType, TImageType::ImageDimension >;
      cbica::WriteImage< CurrentImageType >(cbica::ReadImage< CurrentImageType >(inputImageFile), outputImageFile);
    }
    else if (targetImageFile == "longlong")
    {
      using DefaultPixelType = long long;
      using CurrentImageType = itk::Image< DefaultPixelType, TImageType::ImageDimension >;
      cbica::WriteImage< CurrentImageType >(cbica::ReadImage< CurrentImageType >(inputImageFile), outputImageFile);
    }
    else if (targetImageFile == "float")
    {
      using DefaultPixelType = float;
      using CurrentImageType = itk::Image< DefaultPixelType, TImageType::ImageDimension >;
      cbica::WriteImage< CurrentImageType >(cbica::ReadImage< CurrentImageType >(inputImageFile), outputImageFile);
    }
    else if (targetImageFile == "double")
    {
      using DefaultPixelType = double;
      using CurrentImageType = itk::Image< DefaultPixelType, TImageType::ImageDimension >;
      cbica::WriteImage< CurrentImageType >(cbica::ReadImage< CurrentImageType >(inputImageFile), outputImageFile);
    }
    else
    {
      std::cerr << "Undefined pixel type cast requested, cannot process.\n";
      return EXIT_FAILURE;
    }

    std::cout << "Casting completed.\n";
  }
  else if (requestedAlgorithm == TestComparison)
  {
    auto diffFilter = itk::Testing::ComparisonImageFilter< TImageType, TImageType >::New();
    auto inputImage = cbica::ReadImage< TImageType >(inputImageFile);
    auto size = inputImage->GetBufferedRegion().GetSize();
    diffFilter->SetValidInput(cbica::ReadImage< TImageType >(targetImageFile));
    diffFilter->SetTestInput(inputImage);
    if (cbica::ImageSanityCheck(inputImageFile, targetImageFile))
    {
      diffFilter->VerifyInputInformationOn();
      diffFilter->SetDifferenceThreshold(testThresh);
      diffFilter->SetToleranceRadius(testRadius);
      diffFilter->UpdateLargestPossibleRegion();

      size_t totalSize = 1;
      for (size_t d = 0; d < TImageType::ImageDimension; d++)
      {
        totalSize *= size[d];
      }

      std::cout << "Total Voxels/Pixels in Image: " << totalSize << "\n";
      std::cout << "Number of Difference Voxels : " << diffFilter->GetNumberOfPixelsWithDifferences() << "\n";
      std::cout << "Percentage of Diff Voxels   : " << diffFilter->GetNumberOfPixelsWithDifferences() * 100 / totalSize << "\n";
      std::cout << "Minimum Intensity Difference: " << diffFilter->GetMinimumDifference() << "\n";
      std::cout << "Maximum Intensity Difference: " << diffFilter->GetMaximumDifference() << "\n";
      std::cout << "Average Intensity Difference: " << diffFilter->GetMeanDifference() << "\n";
      std::cout << "Overall Intensity Difference: " << diffFilter->GetTotalDifference() << "\n";
    }
    else
    {
      std::cerr << "Images are in different spaces (size/origin/spacing mismatch).\n";
      return EXIT_FAILURE;
    }
  }
  else if (requestedAlgorithm == BoundingBox)
  {
    if (!cbica::ImageSanityCheck(inputImageFile, targetImageFile))
    {
      std::cerr << "Input image and mask are in different spaces, cannot compute bounding box.\n";
      return EXIT_FAILURE;
    }

    auto inputImage = cbica::ReadImage< TImageType >(inputImageFile);
    auto maskImage = cbica::ReadImage< TImageType >(targetImageFile);
    auto outputImage = cbica::CreateImage< TImageType >(maskImage);
    auto size = maskImage->GetBufferedRegion().GetSize();

    auto nonZeroIndeces = cbica::GetNonZeroIndeces< TImageType >(maskImage);

    using PointSetType = itk::PointSet< typename TImageType::PixelType, TImageType::ImageDimension >;
    using PointIdentifier = typename PointSetType::PointIdentifier;
    using PointType = typename PointSetType::PointType;

    auto pointSet = PointSetType::New();
    auto points = pointSet->GetPoints();

    for (size_t i = 0; i < nonZeroIndeces.size(); i++)
    {
      PointType currentPoint;
      for (size_t d = 0; d < TImageType::ImageDimension; d++)
      {
        currentPoint[d] = nonZeroIndeces[i][d];
      }
      points->InsertElement(i, currentPoint);
    }

    auto boundingBoxCalculator = itk::BoundingBox< PointIdentifier, TImageType::ImageDimension, typename TImageType::PixelType >::New();
    boundingBoxCalculator->SetPoints(points);
    boundingBoxCalculator->ComputeBoundingBox();

    auto max = boundingBoxCalculator->GetMaximum();
    auto min = boundingBoxCalculator->GetMinimum();
    auto center = boundingBoxCalculator->GetCenter();

    std::vector< float > distances;

    for (size_t d = 0; d < TImageType::ImageDimension; d++)
    {
      distances.push_back(std::abs(max[d] - min[d]));
    }

    size_t maxDist = 0;
    for (size_t i = 1; i < TImageType::ImageDimension; i++)
    {
      if (distances[maxDist] < distances[i])
      {
        maxDist = i;
      }
    }

    auto minFromCenter = center, maxFromCenter = center;
    for (size_t d = 0; d < TImageType::ImageDimension; d++)
    {
      size_t isoTropicCheck = maxDist;
      if (!boundingBoxIsotropic)
      {
        isoTropicCheck = d;
      }
      minFromCenter[d] = std::round(minFromCenter[d] - distances[isoTropicCheck] / 2);
      if (minFromCenter[d] < 0)
      {
        minFromCenter[d] = 0;
      }

      maxFromCenter[d] = std::round(maxFromCenter[d] + distances[isoTropicCheck] / 2);
      if (maxFromCenter[d] > size[d])
      {
        maxFromCenter[d] = size[d];
      }
    }

    itk::ImageRegionConstIterator< TImageType > iterator(inputImage, inputImage->GetBufferedRegion());
    itk::ImageRegionIterator< TImageType > outputIterator(outputImage, outputImage->GetBufferedRegion());

    switch (TImageType::ImageDimension)
    {
    case 2:
    {
      for (size_t x = minFromCenter[0]; x < maxFromCenter[0]; x++)
      {
        for (size_t y = minFromCenter[1]; x < maxFromCenter[1]; x++)
        {
          typename TImageType::IndexType currentIndex;
          currentIndex[0] = x;
          currentIndex[1] = y;

          iterator.SetIndex(currentIndex);
          outputIterator.SetIndex(currentIndex);
          outputIterator.Set(iterator.Get());
        }
      }
      break;
    }
    case 3:
    {
      for (size_t x = minFromCenter[0]; x < maxFromCenter[0]; x++)
      {
        for (size_t y = minFromCenter[1]; y < maxFromCenter[1]; y++)
        {
          for (size_t z = minFromCenter[2]; z < maxFromCenter[2]; z++)
          {
            typename TImageType::IndexType currentIndex;
            currentIndex[0] = x;
            currentIndex[1] = y;
            currentIndex[2] = z;

            iterator.SetIndex(currentIndex);
            outputIterator.SetIndex(currentIndex);
            outputIterator.Set(iterator.Get());
          }
        }
      }
      break;
    }
    default:
      std::cerr << "Images other than 2D or 3D are not supported, right now.\n";
      return EXIT_FAILURE;
    }

    cbica::WriteImage< TImageType >(outputImage, outputImageFile);
  }
  else if (requestedAlgorithm == ThresholdAbove)
  {
    auto thresholder = itk::ThresholdImageFilter< TImageType >::New();
    thresholder->SetInput(cbica::ReadImage< TImageType >(inputImageFile));
    thresholder->SetOutsideValue(thresholdOutsideValue);
    thresholder->ThresholdAbove(thresholdAbove);
    thresholder->Update();
    cbica::WriteImage< TImageType >(thresholder->GetOutput(), outputImageFile);
  }
  else if (requestedAlgorithm == ThresholdBelow)
  {
    auto thresholder = itk::ThresholdImageFilter< TImageType >::New();
    thresholder->SetInput(cbica::ReadImage< TImageType >(inputImageFile));
    thresholder->SetOutsideValue(thresholdOutsideValue);
    thresholder->ThresholdBelow(thresholdBelow);
    thresholder->Update();
    cbica::WriteImage< TImageType >(thresholder->GetOutput(), outputImageFile);
  }
  else if (requestedAlgorithm == ThresholdAboveAndBelow)
  {
    auto thresholder = itk::ThresholdImageFilter< TImageType >::New();
    thresholder->SetInput(cbica::ReadImage< TImageType >(inputImageFile));
    thresholder->SetOutsideValue(thresholdOutsideValue);
    thresholder->ThresholdOutside(thresholdBelow, thresholdAbove);
    thresholder->Update();
    cbica::WriteImage< TImageType >(thresholder->GetOutput(), outputImageFile);
  }
  else if (requestedAlgorithm == ThresholdBinary)
  {
    auto thresholder = itk::BinaryThresholdImageFilter< TImageType, TImageType >::New();
    thresholder->SetInput(cbica::ReadImage< TImageType >(inputImageFile));
    thresholder->SetOutsideValue(thresholdOutsideValue);
    thresholder->SetInsideValue(thresholdInsideValue);
    thresholder->SetLowerThreshold(thresholdBelow);
    thresholder->SetUpperThreshold(thresholdAbove);
    thresholder->Update();
    cbica::WriteImage< TImageType >(thresholder->GetOutput(), outputImageFile);
  }
  else if (requestedAlgorithm == ThresholdOtsu)
  {
    auto thresholder = itk::OtsuThresholdImageFilter< TImageType, TImageType >::New();
    thresholder->SetInput(cbica::ReadImage< TImageType >(inputImageFile));
    if (!inputMaskFile.empty())
    {
      thresholder->SetMaskImage(cbica::ReadImage< TImageType >(inputMaskFile));
    }
    thresholder->SetOutsideValue(thresholdOutsideValue);
    thresholder->SetInsideValue(thresholdInsideValue);
    thresholder->Update();
    std::cout << "Otsu Threshold Value: " << thresholder->GetThreshold() << "\n";
    cbica::WriteImage< TImageType >(thresholder->GetOutput(), outputImageFile);
  }
  else if (requestedAlgorithm == ConvertFormat)
  {
    cbica::WriteImage< TImageType >(cbica::ReadImage< TImageType >(inputImageFile), outputImageFile);
  }
  else if (requestedAlgorithm == Image2World)
  {
    auto inputImage = cbica::ReadImage< TImageType >(inputImageFile);
    auto coordinate_split = cbica::stringSplit(coordinateToTransform, ",");
    if (coordinate_split.size() != TImageType::ImageDimension)
    {
      std::cerr << "Please provide the coordinate to transform in the same format as the input image (i.e., N coordinates for a ND image).\n";
      return EXIT_FAILURE;
    }
    typename TImageType::IndexType indexToConvert;
    for (size_t d = 0; d < TImageType::ImageDimension; d++)
    {
      indexToConvert[d] = std::atoi(coordinate_split[d].c_str());
    }
    typename TImageType::PointType output;
    inputImage->TransformIndexToPhysicalPoint(indexToConvert, output);
    std::cout << indexToConvert << " ==> " << output << "\n";
  }
  else if (requestedAlgorithm == World2Image)
  {
    auto inputImage = cbica::ReadImage< TImageType >(inputImageFile);
    auto coordinate_split = cbica::stringSplit(coordinateToTransform, ",");
    if (coordinate_split.size() != TImageType::ImageDimension)
    {
      std::cerr << "Please provide the coordinate to transform in the same format as the input image (i.e., N coordinates for a ND image).\n";
      return EXIT_FAILURE;
    }
    typename TImageType::IndexType output;
    typename TImageType::PointType indexToConvert;
    for (size_t d = 0; d < TImageType::ImageDimension; d++)
    {
      indexToConvert[d] = std::atof(coordinate_split[d].c_str());
    }
    inputImage->TransformPhysicalPointToIndex(indexToConvert, output);
    std::cout << indexToConvert << " ==> " << output << "\n";
  }
  else if (requestedAlgorithm == LabelSimilarity)
  {
    // this filter only works on unsigned int type
    using DefaultImageType = itk::Image< unsigned int, TImageType::ImageDimension >;
    auto inputImage = cbica::ReadImage< DefaultImageType >(inputImageFile);
    auto referenceImage = cbica::ReadImage< DefaultImageType >(referenceMaskForSimilarity);
    auto uniqueLabels = cbica::GetUniqueValuesInImage< DefaultImageType >(inputImage);
    auto uniqueLabelsRef = cbica::GetUniqueValuesInImage< DefaultImageType >(referenceImage);

    // sanity check
    if (uniqueLabels.size() != uniqueLabelsRef.size())
    {
      std::cerr << "The number of unique labels in input and reference image are not consistent.\n";
      return EXIT_FAILURE;
    }
    else
    {
      for (size_t i = 0; i < uniqueLabels.size(); i++)
      {
        if (uniqueLabels[i] != uniqueLabelsRef[i])
        {
          std::cerr << "The label values in input and reference image are not consistent.\n";
          return EXIT_FAILURE;
        }
      }
    }

    auto similarityFilter = itk::LabelOverlapMeasuresImageFilter< DefaultImageType >::New();

    similarityFilter->SetSourceImage(inputImage);
    similarityFilter->SetTargetImage(referenceImage);
    similarityFilter->Update();

    //std::cout << "=== Entire Masked Area ===\n";
    std::cout << "Property,Value\n";
    std::cout << "TotalOverlap," << similarityFilter->GetTotalOverlap() << "\n";
    std::cout << "Union(Jaccard)_Overall," << similarityFilter->GetUnionOverlap() << "\n";
    std::cout << "Mean(DICE)_Overall," << similarityFilter->GetMeanOverlap() << "\n";
    std::cout << "VolumeSimilarity_Overall," << similarityFilter->GetVolumeSimilarity() << "\n";
    std::cout << "FalseNegativeError_Overall," << similarityFilter->GetFalseNegativeError() << "\n";
    std::cout << "FalsePositiveError_Overall," << similarityFilter->GetFalsePositiveError() << "\n";

    if (uniqueLabels.size() > 2) // basically if there is something more than 0 and 1
    {
      //std::cout << "=== Individual Labels ===\n";
      //std::cout << "Property,Value\n";
      for (size_t i = 0; i < uniqueLabels.size(); i++)
      {
        auto uniqueLabels_string = std::to_string(uniqueLabels[i]);
        std::cout << "LabelValue," << uniqueLabels_string << "\n";
        std::cout << "TargetOverlap_Label" + uniqueLabels_string + "," << similarityFilter->GetTargetOverlap(uniqueLabels[i]) << "\n";
        std::cout << "Union(Jaccard)_Label" + uniqueLabels_string + "," << similarityFilter->GetUnionOverlap(uniqueLabels[i]) << "\n";
        std::cout << "Mean(DICE)_Label" + uniqueLabels_string + "," << similarityFilter->GetMeanOverlap(uniqueLabels[i]) << "\n";
        std::cout << "VolumeSimilarity_Label" + uniqueLabels_string + "," << similarityFilter->GetVolumeSimilarity(uniqueLabels[i]) << "\n";
        std::cout << "FalseNegativeError_Label" + uniqueLabels_string + "," << similarityFilter->GetFalseNegativeError(uniqueLabels[i]) << "\n";
        std::cout << "FalsePositiveError_Label" + uniqueLabels_string + "," << similarityFilter->GetFalsePositiveError(uniqueLabels[i]) << "\n";
      }
    }
      return EXIT_SUCCESS;
  }
  else if (requestedAlgorithm == Hausdorff)
  {
    auto filter = itk::HausdorffDistanceImageFilter< TImageType, TImageType >::New();
    filter->SetInput1(cbica::ReadImage< TImageType >(inputImageFile));
    filter->SetInput2(cbica::ReadImage< TImageType >(referenceMaskForSimilarity));
    try
    {
      filter->Update();
    }
    catch (const std::exception&e)
    {
      std::cerr << "Hausdorff calculation error: " << e.what() << "\n";
      return EXIT_FAILURE;
    }   

    std::cout << "Hausdorff Distance: " << filter->GetHausdorffDistance() << "\n";
  }

  // if no other algorithm has been selected and mask & output files are present and in same space as input, apply it
  else if (cbica::isFile(inputMaskFile) && !outputImageFile.empty())
  {
    auto errorMessage = "The input mask and mask are not in the same space.\n";
    if (cbica::ImageSanityCheck(inputMaskFile, inputImageFile))
    {
      auto masker = itk::MaskImageFilter< TImageType, TImageType >::New();
      masker->SetInput(cbica::ReadImage< TImageType >(inputImageFile));
      masker->SetMaskImage(cbica::ReadImage< TImageType >(inputMaskFile));
      masker->Update();
      cbica::WriteImage< TImageType >(masker->GetOutput(), outputImageFile);
    }
    else if (TImageType::ImageDimension == 4)
    {
      // input image is 4D and mask is 3D
      if (!cbica::ImageSanityCheck(inputMaskFile, inputImageFile, true))
      {
        std::cerr << errorMessage;
        return EXIT_FAILURE;
      }

      auto inputImage = cbica::ReadImage< TImageType >(inputImageFile);
      using TBaseImageType = itk::Image< typename TImageType::PixelType, 3 >;
      auto maskImage = cbica::ReadImage< TBaseImageType >(inputMaskFile);
      auto inputImages = cbica::GetExtractedImages<
        TImageType, TBaseImageType >(
          inputImage);

      std::vector< typename TBaseImageType::Pointer > outputImages;
      outputImages.resize(inputImages.size());

      for (size_t i = 0; i < inputImages.size(); i++)
      {
        auto masker = itk::MaskImageFilter< TBaseImageType, TBaseImageType >::New();
        masker->SetInput(inputImages[i]);
        masker->SetMaskImage(maskImage);
        masker->Update();
        outputImages[i] = masker->GetOutput();
      }

      auto output = cbica::GetJoinedImage< TBaseImageType, TImageType >(outputImages, inputImage->GetSpacing()[3]);
      cbica::WriteImage< TImageType >(output, outputImageFile);
    }
    else
    {
      std::cerr << errorMessage;
      return EXIT_FAILURE;
    }
  }
  else
  {
    // no algorithm has been selected
  }
  
  return EXIT_SUCCESS;
}

//! helper function to image stack joining
template< class TImageType, class TOutputImageType >
int algorithmsRunner_imageStack2join(std::vector< std::string >& inputImageFiles)
{
  std::vector< typename TImageType::Pointer > inputImages;
  inputImages.resize(inputImageFiles.size());
  for (size_t i = 0; i < inputImageFiles.size(); i++)
  {
    inputImages[i] = cbica::ReadImage< TImageType >(inputImageFiles[i]);
  }
  auto output = cbica::GetJoinedImage< TImageType, TOutputImageType >(inputImages, imageStack2JoinSpacing);
  cbica::WriteImage< TOutputImageType >(output, outputImageFile);
  return EXIT_SUCCESS;
}

//! helper function to image stack extraction
template< class TImageType, class TOutputImageType >
int algorithmsRunner_join2imageStack()
{
  auto outputImages = cbica::GetExtractedImages< TImageType, TOutputImageType >(cbica::ReadImage< TImageType >(inputImageFile));
  std::string path, base, ext;
  cbica::splitFileName(outputImageFile, path, base, ext);
  cbica::createDir(path);

  for (size_t i = 0; i < outputImages.size(); i++)
  {
    cbica::WriteImage< TOutputImageType >(
      outputImages[i],
      cbica::normalizePath(path + "/" + base + "_" + std::to_string(i) + ".nii.gz")
      );
  }
  return EXIT_SUCCESS;
}


int main(int argc, char** argv)
{
  cbica::CmdParser parser(argc, argv, "Utilities");

  parser.addOptionalParameter("i", "inputImage", cbica::Parameter::STRING, "File or Dir", "Input Image (all CaPTk supported images) for processing", "Directory to a single series DICOM only");
  parser.addOptionalParameter("m", "maskImage", cbica::Parameter::STRING, "File or Dir", "Input Mask (all CaPTk supported images) for processing", "Directory to a single series DICOM only");
  parser.addOptionalParameter("o", "outputImage", cbica::Parameter::FILE, "NIfTI", "Output Image for processing");
  parser.addOptionalParameter("df", "dicomDirectory", cbica::Parameter::DIRECTORY, "none", "Absolute path of directory containing single dicom series");
  parser.addOptionalParameter("r", "resize", cbica::Parameter::INTEGER, "10-500", "Resize an image based on the resizing factor given", "Example: -r 150 resizes inputImage by 150%", "Defaults to 100, i.e., no resizing", "Resampling can be done on image with 100");
  parser.addOptionalParameter("rr", "resizeResolution", cbica::Parameter::STRING, "0-10", "[Resample] Resolution of the voxels/pixels to change to", "Resize value needs to be 100", "Defaults to " + resamplingResolution_full, "Use '-rf' for a reference file");
  parser.addOptionalParameter("rf", "resizeReference", cbica::Parameter::FILE, "NIfTI image", "[Resample] Reference image on which resampling is to be done", "Resize value needs to be 100", "Use '-ri' for resize resolution");
  parser.addOptionalParameter("ri", "resizeInterp", cbica::Parameter::STRING, "NEAREST:LINEAR:BSPLINE:BICUBIC", "[Resample] The interpolation type to use for resampling or resizing", "Defaults to LINEAR");
  parser.addOptionalParameter("s", "sanityCheck", cbica::Parameter::FILE, "NIfTI Reference", "Do sanity check of inputImage with the file provided in with this parameter", "Performs checks on size, origin & spacing", "Pass the target image after '-s'");
  parser.addOptionalParameter("inf", "information", cbica::Parameter::BOOLEAN, "true or false", "Output the information in inputImage", "If DICOM file is detected, the tags are written out");
  parser.addOptionalParameter("c", "cast", cbica::Parameter::STRING, "(u)char, (u)int, (u)long, (u)longlong, float, double", "Change the input image type", "Examples: '-c uchar', '-c float', '-c longlong'");
  parser.addOptionalParameter("un", "uniqueVals", cbica::Parameter::BOOLEAN, "true or false", "Output the unique values in the inputImage", "Pass value '1' for ascending sort or '0' for no sort", "Defaults to '1'");
  parser.addOptionalParameter("b", "boundingBox", cbica::Parameter::FILE, "NIfTI Mask", "Extracts the smallest bounding box around the mask file", "With respect to inputImage", "Writes to outputImage");
  parser.addOptionalParameter("bi", "boundingIso", cbica::Parameter::BOOLEAN, "Isotropic Box or not", "Whether the bounding box is Isotropic or not", "Defaults to true");
  parser.addOptionalParameter("tb", "testBase", cbica::Parameter::FILE, "NIfTI Reference", "Baseline image to compare inputImage with");
  parser.addOptionalParameter("tr", "testRadius", cbica::Parameter::INTEGER, "0-10", "Maximum distance away to look for a matching pixel", "Defaults to 0");
  parser.addOptionalParameter("tt", "testThresh", cbica::Parameter::FLOAT, "0-5", "Minimum threshold for pixels to be different", "Defaults to 0.0");
  parser.addOptionalParameter("cm", "createMask", cbica::Parameter::STRING, "N.A.", "Create a binary mask out of a provided (float) thresholds","Format: -cm lower,upper", "Output is 1 if value >= lower or <= upper", "Defaults to 1,Max");
  parser.addOptionalParameter("cv", "changeValue", cbica::Parameter::STRING, "N.A.", "Change the specified pixel/voxel value", "Format: -cv oldValue1xoldValue2,newValue1xnewValue2", "Can be used for multiple number of value changes", "Defaults to 3,4");
  parser.addOptionalParameter("d2n", "dicom2Nifti", cbica::Parameter::FILE, "NIfTI Reference", "If path to reference is present, then image comparison is done", "Use '-i' to pass input DICOM image", "Use '-o' to pass output image file");
  parser.addOptionalParameter("n2d", "nifi2dicom", cbica::Parameter::DIRECTORY, "DICOM Reference", "A reference DICOM is passed after this parameter", "The header information from the DICOM reference is taken to write output", "Use '-i' to pass input NIfTI image", "Use '-o' to pass output DICOM directory");
  parser.addOptionalParameter("ds", "dcmSeg", cbica::Parameter::DIRECTORY, "DICOM Reference", "A reference DICOM is passed after this parameter", "The header information from the DICOM reference is taken to write output", "Use '-i' to pass input NIfTI image", "Use '-o' to pass output DICOM file");
  parser.addOptionalParameter("dsJ", "dcmSegJSON", cbica::Parameter::FILE, "JSON file for Metadata", "The extra metadata needed to generate the DICOM-Seg object", "Use http://qiicr.org/dcmqi/#/seg to create it", "Use '-i' to pass input NIfTI segmentation image", "Use '-o' to pass output DICOM file");
  parser.addOptionalParameter("or", "orient", cbica::Parameter::STRING, "Desired 3 letter orientation", "The desired orientation of the image", "See the following for supported orientations (use last 3 letters only):", "https://itk.org/Doxygen/html/namespaceitk_1_1SpatialOrientation.html#a8240a59ae2e7cae9e3bad5a52ea3496e");
  parser.addOptionalParameter("thA", "threshAbove", cbica::Parameter::FLOAT, "Desired_Threshold", "The intensity ABOVE which pixels of the input image will be", "made to OUTSIDE_VALUE (use '-tOI')", "Generates a grayscale image");
  parser.addOptionalParameter("thB", "threshBelow", cbica::Parameter::FLOAT, "Desired_Threshold", "The intensity BELOW which pixels of the input image will be", "made to OUTSIDE_VALUE (use '-tOI')", "Generates a grayscale image");
  parser.addOptionalParameter("tAB", "threshAnB", cbica::Parameter::STRING, "Lower_Threshold,Upper_Threshold", "The intensities outside Lower and Upper thresholds will be", "made to OUTSIDE_VALUE (use '-tOI')", "Generates a grayscale image");
  parser.addOptionalParameter("thO", "threshOtsu", cbica::Parameter::BOOLEAN, "0-1", "Whether to do Otsu threshold", "Generates a binary image which has been thresholded using Otsu", "Use '-tOI' to set Outside and Inside Values", "Optional mask to localize Otsu search area");
  parser.addOptionalParameter("tBn", "thrshBinary", cbica::Parameter::STRING, "Lower_Threshold,Upper_Threshold", "The intensity BELOW and ABOVE which pixels of the input image will be", "made to OUTSIDE_VALUE (use '-tOI')", "Default for OUTSIDE_VALUE=0");
  parser.addOptionalParameter("tOI", "threshOutIn", cbica::Parameter::STRING, "Outside_Value,Inside_Value", "The values that will go inside and outside the thresholded region", "Defaults to '0,1', i.e., a binary output");
  parser.addOptionalParameter("i2w", "image2world", cbica::Parameter::STRING, "x,y,z", "The world coordinates that will be converted to image coordinates for the input image", "Example: '-i2w 10,20,30'");
  parser.addOptionalParameter("w2i", "world2image", cbica::Parameter::STRING, "i,j,k", "The image coordinates that will be converted to world coordinates for the input image", "Example: '-w2i 10.5,20.6,30.2'");
  parser.addOptionalParameter("j2e", "joined2extracted", cbica::Parameter::BOOLEAN, "0-1", "Axis to extract is always the final axis (axis '3' for a 4D image)", "The '-o' parameter can be used for output: '-o /path/to/extracted_'");
  parser.addOptionalParameter("e2j", "extracted2joined", cbica::Parameter::FLOAT, "0-10", "The spacing in the new direction", "Pass the folder containing all images in '-i'");
  parser.addOptionalParameter("ls", "labelSimilarity", cbica::Parameter::FILE, "NIfTI Reference", "Calculate similarity measures for 2 label maps", "Pass the reference map after '-ls' and the comparison will be done with '-i'", "For images with more than 2 labels, individual label stats are also presented");
  parser.addOptionalParameter("hd", "hausdorffDist", cbica::Parameter::FILE, "NIfTI Reference", "Calculate the Hausdorff Distance for the input image and", "the one passed after '-hd'");

  parser.addExampleUsage("-i C:/test.nii.gz -o C:/test_int.nii.gz -c int", "Cast an image pixel-by-pixel to a signed integer");
  parser.addExampleUsage("-i C:/test.nii.gz -o C:/test_75.nii.gz -r 75 -ri linear", "Resize an image by 75% using linear interpolation");
  parser.addExampleUsage("-i C:/test.nii.gz -inf", "Prints out image information to console (for DICOMs, this does a full dump of the tags)");
  parser.addExampleUsage("-i C:/test/1.dcm -o C:/test.nii.gz -d2n C:/test_reference.nii.gz", "DICOM to NIfTI conversion and do sanity check of the converted image with the reference image");
  parser.addExampleUsage("-i C:/test/1.nii.gz -o C:/test_dicom/ -n2d C:/referenceDICOM/", "NIfTI to DICOM conversion and do sanity check of the converted image with the reference image");
  parser.addExampleUsage("-i C:/test/1.nii.gz -o C:/test_dicom/ -ds C:/referenceDICOM/ -dsJ C:/dicomSeg.json", "NIfTI Segmentation to DICOM-Seg conversion using the supplied reference DICOM and JSON");
  parser.addExampleUsage("-i C:/test/1.nii.gz -o C:/output.nii.gz -or RAI", "Re-orient input image to RAI");
  parser.addExampleUsage("-i C:/test/1.nii.gz -o C:/output.nii.gz -thO 1", "Otsu Threshold");
  parser.addExampleUsage("-i C:/test/1.nii.gz -o C:/output.nii.gz -tBn 50,100", "Binary Threshold between 50 and 100 with default outside & inside values");
  parser.addExampleUsage("-i C:/test/1.nii.gz -o C:/output.nii.gz -tAB 50,100 -tOI -100,10000", "Above & Below Threshold between 50 and 100 with outside value -100 and inside value 10000");
  parser.addExampleUsage("-i C:/test/1.nii.gz -o C:/output -j2e 1", "Extract the joined image into its series");
  parser.addExampleUsage("-i C:/test/ -o C:/output.nii.gz -e2j 1.5", "Join the extracted images into a single image with spacing in the new dimension as 1.5");
  parser.addExampleUsage("-i C:/test/input.nii.gz -o C:/output.nii.gz -rr 1.0 -ri LINEAR", "Calculates an isotropic image from the input with spacing '1.0' in all dimensions using linear interpolation");
  parser.addExampleUsage("-i C:/test/input.nii.gz -o C:/output.nii.gz -rr 1.0,2.0,3.0 -ri LINEAR", "Calculates an anisotropic image from the input with spacing '1.0' in x, '2.0' in y and '3.0' in z using linear interpolation");
  parser.addExampleUsage("-i C:/test/input.nii.gz -o C:/output.nii.gz -rf C:/reference.nii.gz -ri LINEAR", "Calculates an isotropic image from the input with spacing from the reference image using linear interpolation");
  parser.addExampleUsage("-i C:/test/outputMask.nii.gz -l2s C:/referenceMask.nii.gz", "Calculates Total/Union/Mean Overlap (different DICE coefficients), Volume Similarity, False Positive/Negative Error for all labels");

  parser.addApplicationDescription("This application has various utilities that can be used for constructing pipelines around CaPTk's functionalities. Please add feature requests on the CaPTk GitHub page at https://github.com/CBICA/CaPTk.");
  
  if (parser.isPresent("i"))
  {
    parser.getParameterValue("i", inputImageFile);
  }
  if (parser.isPresent("m"))
  {
    parser.getParameterValue("m", inputMaskFile);
  }
  if (parser.isPresent("o"))
  {
    parser.getParameterValue("o", outputImageFile);
  }
  if(parser.isPresent("df"))
  {
    requestedAlgorithm = DicomLoadTesting;
    parser.getParameterValue("df", dicomFolderPath);
    parser.getParameterValue("i", inputImageFile);
    parser.getParameterValue("tt", testThresh);
  }
  else if (parser.isPresent("d2n"))
  {
    requestedAlgorithm = Dicom2Nifti;

	// check if the Nifti output filename has been specified?
	if (parser.isPresent("o"))
	{
		parser.getParameterValue("o", outputImageFile);
	}
	else
	{
		std::cerr << "DICOM2Nifti conversion requested but the output Nifti filename was not found. Please use the '-o' parameter.\n";
		return EXIT_FAILURE;
	}

	// check if the Dicom input has been specified?
	if (parser.isPresent("i"))
	{
		parser.getParameterValue("i", inputImageFile);
	}
	else
	{
		std::cerr << "DICOM2Nifti conversion requested but the input Dicom data was not found. Please use the '-i' parameter.\n";
		return EXIT_FAILURE;
	}
    parser.getParameterValue("d2n", targetImageFile);

  }
  else if (parser.isPresent("n2d"))
  {
    requestedAlgorithm = Nifti2Dicom;
    parser.getParameterValue("n2d", targetImageFile); // in this case, it is the DICOM reference file
    parser.getParameterValue("o", outputImageFile);
  }
  else if (parser.isPresent("ds"))
  {
    requestedAlgorithm = Nifti2DicomSeg;
    parser.getParameterValue("ds", targetImageFile); // in this case, it is the DICOM reference file
    parser.getParameterValue("o", outputImageFile);
    if (parser.isPresent("dsJ"))
    {
      parser.getParameterValue("dsJ", dicomSegJSON);
    }
    else
    {
      std::cerr << "DICOM-Seg conversion requested but the JSON file containing the metadata was not found. Please use the '-dsJ' parameter.\n";
      return EXIT_FAILURE;
    }
  }
  else if (parser.isPresent("r"))
  {
    requestedAlgorithm = Resize;
    parser.getParameterValue("r", resize);
    if (parser.isPresent("ri"))
    {
      parser.getParameterValue("ri", resamplingInterpolator);
    }
  }
  else if (parser.isPresent("rr"))
  {
    requestedAlgorithm = Resample;
    parser.getParameterValue("rr", resamplingResolution_full);
    if (parser.isPresent("ri"))
    {
      parser.getParameterValue("ri", resamplingInterpolator);
    }
  }
  else if (parser.isPresent("rf"))
  {
    requestedAlgorithm = Resample;
    parser.getParameterValue("rf", resamplingReference);
    if (parser.isPresent("ri"))
    {
      parser.getParameterValue("ri", resamplingInterpolator);
    }
  }
  else if (parser.isPresent("s"))
  {
    parser.getParameterValue("s", targetImageFile);
    requestedAlgorithm = SanityCheck;
  }
  else if (parser.isPresent("inf"))
  {
    requestedAlgorithm = Information;
  }
  else if (parser.isPresent("c"))
  {
    parser.getParameterValue("c", targetImageFile);
    requestedAlgorithm = Casting;
  }
  else if (parser.isPresent("un"))
  {
    parser.getParameterValue("un", uniqueValsSort);
    requestedAlgorithm = UniqueValues;
  }
  else if (parser.isPresent("or"))
  {
    parser.getParameterValue("or", orientationDesired);
    requestedAlgorithm = OrientImage;
  }
  else if (parser.isPresent("tb"))
  {
    parser.getParameterValue("tb", targetImageFile);
    requestedAlgorithm = TestComparison;
    if (parser.isPresent("tr"))
    {
      parser.getParameterValue("tr", testRadius);
    }
    if (parser.isPresent("tt"))
    {
      parser.getParameterValue("tt", testThresh);
    }
  }
  else if (parser.isPresent("b"))
  {
    parser.getParameterValue("b", targetImageFile);
    requestedAlgorithm = BoundingBox;
    if (parser.isPresent("bi"))
    {
      parser.getParameterValue("bi", boundingBoxIsotropic);
    }
  }
  else if (parser.isPresent("cm"))
  {
    std::string temp;
    parser.getParameterValue("cm", temp);
    if (!temp.empty())
    {
      auto temp2 = cbica::stringSplit(temp, ",");
      lowerThreshold = std::atof(temp2[0].c_str());
      if (temp2.size() == 2)
      {
        upperThreshold = std::atof(temp2[1].c_str());
      }
    }

    requestedAlgorithm = CreateMask;
  }
  else if (parser.isPresent("cv"))
  {
    std::string temp;
    parser.getParameterValue("cv", temp);
    if (!temp.empty())
    {
      auto temp2 = cbica::stringSplit(temp, ",");
      if (temp2.size() != 2)
      {
        std::cerr << "Change value needs 2 values in the format '-cv oldValue,newValue' to work.\n";
        return EXIT_FAILURE;
      }
      changeOldValues = temp2[0];
      changeNewValues = temp2[1];
    }

    requestedAlgorithm = ChangeValue;
  }

  /// common for all thresholding
  if (parser.isPresent("tOI"))
  {
    std::string lowerAndUpperVals;
    parser.getParameterValue("tOI", lowerAndUpperVals);
    auto temp = cbica::stringSplit(lowerAndUpperVals, ",");
    thresholdOutsideValue = std::atof(temp[0].c_str());
    if (temp.size() > 1)
    {
      thresholdInsideValue = std::atof(temp[1].c_str());
    }
  }
  ///

  else if (parser.isPresent("thA"))
  {
    requestedAlgorithm = ThresholdAbove;
    std::string thresholds;
    parser.getParameterValue("thA", thresholdAbove);
  }
  
  else if (parser.isPresent("thB"))
  {
    requestedAlgorithm = ThresholdBelow;
    std::string thresholds;
    parser.getParameterValue("thB", thresholdBelow);
  }
  
  else if (parser.isPresent("tAB"))
  {
    requestedAlgorithm = ThresholdAboveAndBelow;
    std::string thresholds;
    parser.getParameterValue("tAB", thresholds);
    auto temp = cbica::stringSplit(thresholds, ",");
    if (temp.size() < 2)
    {
      std::cerr << "Please provide 2 thresholds, for instance '-tAB 15,100'.\n";
    }
    thresholdBelow = std::atof(temp[0].c_str());
    thresholdAbove = std::atof(temp[1].c_str());
  }

  else if (parser.isPresent("thO"))
  {
    bool temp;
    parser.getParameterValue("thO", temp);
    if (temp)
    {
      requestedAlgorithm = ThresholdOtsu;
    }
  }

  else if (parser.isPresent("tBn"))
  {
    requestedAlgorithm = ThresholdBinary;
    std::string thresholds;
    parser.getParameterValue("tBn", thresholds);
    auto temp = cbica::stringSplit(thresholds, ",");
    if (temp.size() < 2)
    {
      std::cerr << "Please provide 2 thresholds, for instance '-tBn 15,100'.\n";
    }
    thresholdBelow = std::atof(temp[0].c_str());
    thresholdAbove = std::atof(temp[1].c_str());
  }

  else if (parser.isPresent("cov"))
  {
    // will not check if the flag is enabled or not; will simply go from 
    // inputImageFile -> outputImageFile
    requestedAlgorithm = ConvertFormat;
  }

  else if (parser.isPresent("i2w"))
  {
    requestedAlgorithm = Image2World;
    parser.getParameterValue("i2w", coordinateToTransform);
    if (coordinateToTransform.empty())
    {
      std::cerr << "Please provide the coordinate to transform in a comma-separated format: '-i2w 10,20,30'\n";
      return EXIT_FAILURE;
    }
  }

  else if (parser.isPresent("w2i"))
  {
    requestedAlgorithm = World2Image;
    parser.getParameterValue("w2i", coordinateToTransform);
    if (coordinateToTransform.empty())
    {
      std::cerr << "Please provide the coordinate to transform in a comma-separated format: '-w2i 10.5,20.6,30.8'\n";
      return EXIT_FAILURE;
    }
  }

  else if (parser.isPresent("j2e"))
  {
    requestedAlgorithm = JoinedImage2Stack;

    auto imageInfo_first = cbica::ImageInfo(inputImageFile);

    switch (imageInfo_first.GetImageDimensions())
    {
    case 3:
    {
      using TImageType = itk::Image<float, 3>;
      using TOImageType = itk::Image<float, 2>;
      algorithmsRunner_join2imageStack< TImageType, TOImageType >();
    }
    case 4:
    {
      using TImageType = itk::Image<float, 4>;
      using TOImageType = itk::Image<float, 3>;
      algorithmsRunner_join2imageStack< TImageType, TOImageType >();
    }
    default:
      break;
    }
    return EXIT_SUCCESS;
  }
  else if (parser.isPresent("e2j"))
  {
    requestedAlgorithm = ImageStack2Join;
    int pos;
    parser.compareParameter("e2j", pos);
    if (argc <= pos + 1)
    {
      std::cerr << "Please provide the spacing along the joined axis: '-e2j 1.0'.\n";
      return EXIT_FAILURE;
    }
    parser.getParameterValue("e2j", imageStack2JoinSpacing);

    auto imagesToJoin = cbica::filesInDirectory(inputImageFile);
    if (imagesToJoin.empty())
    {
      std::cerr << "Please pass a folder containing images the join using '-i'.\n";
      return EXIT_FAILURE;
    }
    auto imageInfo_first = cbica::ImageInfo(imagesToJoin[0]);
    for (size_t i = 1; i < imagesToJoin.size(); i++)
    {
      auto currentImageInfo = cbica::ImageInfo(imagesToJoin[i]);
      if (imageInfo_first.GetImageDimensions() != currentImageInfo.GetImageDimensions())
      {
        std::cerr << "The image '" << imagesToJoin[i] << "' has a different dimension with the first image in series; cannot proceed.\n";
        return EXIT_FAILURE;
      }      
    }

    switch (imageInfo_first.GetImageDimensions())
    {
    case 2:
    {
      using TImageType = itk::Image<float, 2>;
      using TOImageType = itk::Image<float, 3>;
      algorithmsRunner_imageStack2join< TImageType, TOImageType >(imagesToJoin);
    }
    case 3:
    {
      using TImageType = itk::Image<float, 3>;
      using TOImageType = itk::Image<float, 4>;
      algorithmsRunner_imageStack2join< TImageType, TOImageType >(imagesToJoin);
    }
    default:
      break;
    }
    return EXIT_SUCCESS;
  }
  else if (parser.isPresent("ls"))
  {
    requestedAlgorithm = LabelSimilarity;
    parser.getParameterValue("ls", referenceMaskForSimilarity);
  }
  else if (parser.isPresent("hd"))
  {
    requestedAlgorithm = Hausdorff;
    parser.getParameterValue("hd", referenceMaskForSimilarity);
  }

  // this doesn't need any template initialization
  if (requestedAlgorithm == SanityCheck)
  {
    if (cbica::ImageSanityCheck(inputImageFile, targetImageFile))
    {
      std::cout << "Images are in the same space.\n";
      return EXIT_SUCCESS;
    }
    else
    {
      std::cerr << "Images are in different spaces.\n";
      return EXIT_FAILURE;
    }
  }
  auto inputImageInfo = cbica::ImageInfo(inputImageFile);

  if (requestedAlgorithm == Information)
  {
    std::cout << "ITK Image information:.\n";
    auto dims = inputImageInfo.GetImageDimensions();
    auto size = inputImageInfo.GetImageSize();
    auto origin = inputImageInfo.GetImageOrigins();
    auto spacing = inputImageInfo.GetImageSpacings();
    auto directions = inputImageInfo.GetImageDirections();
    auto size_string = std::to_string(size[0]);
    auto origin_string = std::to_string(origin[0]);
    auto spacing_string = std::to_string(spacing[0]);
    auto directions_string = "[" + std::to_string(directions[0][0]) + "x" + std::to_string(directions[0][1]) + "x" + std::to_string(directions[0][2]);
    size_t totalSize = size[0];
    for (size_t i = 1; i < dims; i++)
    {
      size_string += "x" + std::to_string(size[i]);
      origin_string += "x" + std::to_string(origin[i]);
      spacing_string += "x" + std::to_string(spacing[i]);
      directions_string += ";" + std::to_string(directions[i][0]) + "x" + std::to_string(directions[i][1]) + "x" + std::to_string(directions[i][2]);
      totalSize *= size[i];
    }
    directions_string += "]";
    std::cout << "Property,Value\n";
    std::cout << "Dimensions," << dims << "\n";
    std::cout << "Size," << size_string << "\n";
    std::cout << "Total," << totalSize << "\n";
    std::cout << "Origin," << origin_string << "\n";
    std::cout << "Spacing," << spacing_string << "\n";
    std::cout << "Component," << inputImageInfo.GetComponentTypeAsString() << "\n";
    std::cout << "Pixel Type," << inputImageInfo.GetPixelTypeAsString() << "\n";
    std::cout << "Directions," << directions_string << "\n";

    if (cbica::IsDicom(inputImageFile)) // if dicom file
    {
      std::cout << "DICOM file detected, will print out all tags from here.\n";
      DicomMetadataReader reader;
      reader.SetFilePath(inputImageFile);
      bool readStatus = reader.ReadMetaData();
      if (!readStatus)
      {
        std::cout << "Could not read dicom image.\n";
        return EXIT_FAILURE;
      }

      auto readMap = reader.GetMetaDataMap();
      std::pair<std::string, std::string> labelValuePair;
      std::cout << "Tag" << "," << "Description" << "," << "Value" << std::endl;
      for (auto itr = readMap.begin(); itr != readMap.end(); ++itr)
      {
        labelValuePair = itr->second;
        std::cout << itr->first.c_str() << ","
          << labelValuePair.first.c_str() << ","
          << labelValuePair.second.c_str() << "\n";
      }
    }
    return EXIT_SUCCESS;
  }

  switch (inputImageInfo.GetImageDimensions())
  {
  case 2:
  {
    using ImageType = itk::Image< float, 2 >;
    return algorithmsRunner< ImageType >();

    break;
  }
  case 3:
  {
    using ImageType = itk::Image< float, 3 >;

    if (requestedAlgorithm == OrientImage) // this does not work for 2 or 4-D images
    {
      auto output = cbica::GetImageOrientation< ImageType >(cbica::ReadImage< ImageType >(inputImageFile), orientationDesired);
      std::cout << "Original Image Orientation: " << output.first << "\n";
      std::string path, base, ext;
      cbica::splitFileName(outputImageFile, path, base, ext);
      if (ext.find(".nii") != std::string::npos)
      {
        std::cerr << "WARNING: NIfTI files do NOT support orientation properly [https://github.com/InsightSoftwareConsortium/ITK/issues/1042].\n";
      }
      if (ext != ".mha")
      {
        auto tempOutputFile = path + "/" + base + ".mha"; // this is done to ensure NIfTI IO issues are taken care of
        cbica::WriteImage< ImageType >(output.second, tempOutputFile);
        auto reorientedInput = cbica::ReadImage< ImageType >(tempOutputFile);
        cbica::WriteImage< ImageType >(reorientedInput, outputImageFile);
        std::remove(tempOutputFile.c_str());
      }
      else
      {
        cbica::WriteImage< ImageType >(output.second, outputImageFile);
      }
      return EXIT_SUCCESS;
    }

    return algorithmsRunner< ImageType >();

    break;
  }
  case 4:
  {
    using ImageType = itk::Image< float, 4 >;

    if (requestedAlgorithm == OrientImage) // this does not work for 2 or 4-D images
    {
      using OrientationImageType = itk::Image< float, 3 >;
      auto imageSize = inputImageInfo.GetImageSize();

      auto inputImage = cbica::ReadImage< ImageType >(inputImageFile); 

      // set the sub-image properties
      typename ImageType::IndexType regionIndex;
      typename ImageType::SizeType regionSize;
      regionSize[0] = imageSize[0];
      regionSize[1] = imageSize[1];
      regionSize[2] = imageSize[2];
      regionSize[3] = 0;
      regionIndex[0] = 0;
      regionIndex[1] = 0;
      regionIndex[2] = 0;
      regionIndex[3] = 0;

      std::vector<ImageType::Pointer> OnePatientperfusionImages;

      auto joinFilter = itk::JoinSeriesImageFilter< OrientationImageType, ImageType >::New();

      bool originalOrientationOutput = false;
      // loop through time points
      for (size_t i = 0; i < imageSize[3]; i++)
      {
        regionIndex[3] = i;
        typename ImageType::RegionType desiredRegion(regionIndex, regionSize);
        auto filter = itk::ExtractImageFilter< ImageType, OrientationImageType >::New();
        filter->SetExtractionRegion(desiredRegion);
        filter->SetInput(inputImage);
        filter->SetDirectionCollapseToSubmatrix();
        filter->Update();
        auto CurrentTimePoint = filter->GetOutput();

        auto output = cbica::GetImageOrientation< OrientationImageType >(CurrentTimePoint, orientationDesired);
        if (!originalOrientationOutput)
        {
          std::cout << "Original Image Orientation: " << output.first << "\n";
          originalOrientationOutput = true;
        }

        joinFilter->SetInput(i, output.second);
      }

      joinFilter->Update();

      std::string path, base, ext;
      cbica::splitFileName(outputImageFile, path, base, ext);
      if (ext != ".mha")
      {
        auto tempOutputFile = path + "/" + base + ".mha"; // this is done to ensure NIfTI IO issues are taken care of
        cbica::WriteImage< ImageType >(joinFilter->GetOutput(), tempOutputFile);
        cbica::WriteImage< ImageType >(cbica::ReadImage< ImageType >(tempOutputFile), outputImageFile);
        std::remove(tempOutputFile.c_str());
      }
      else
      {
        cbica::WriteImage< ImageType >(joinFilter->GetOutput(), outputImageFile);
      }
      return EXIT_SUCCESS;
    }

    return algorithmsRunner< ImageType >();

    break;
  }
  default:
    std::cerr << "Supplied image has an unsupported dimension of '" << inputImageInfo.GetImageDimensions() << "'; only 2, 3 and 4 D images are supported.\n";
    return EXIT_FAILURE; // exiting here because no further processing should be done on the image
  }

  return EXIT_SUCCESS;
}