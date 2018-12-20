#include "cbicaCmdParser.h"
#include "cbicaLogging.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"

#include "ZScoreNormalizer.h"

#include "itkBoundingBox.h"
#include "itkPointSet.h"

//! Detail the available algorithms to make it easier to initialize
enum AvailableAlgorithms
{
  None,
  HistogramMatching,
  Resize,
  SanityCheck,
  Information,
  Casting,
  UniqueValues,
  TestComparison,
  BoundingBox,
  ZScoreNormalize
};

int requestedAlgorithm = 0;

std::string inputImageFile, outputImageFile, targetImageFile;
size_t resize = 100;
int histoMatchQuantiles = 40, histoMatchBins = 100, testRadius = 0, testNumber = 0;
float testThresh = 0.0, testAvgDiff = 0.0;

bool uniqueValsSort = true, boundingBoxIsotropic = true;

template< class TImageType >
int algorithmsRunner()
{
  if ((requestedAlgorithm == Resize) && (resize != 100))
  {
    auto outputImage = cbica::ResizeImage< TImageType >(cbica::ReadImage< TImageType >(inputImageFile), resize);
    cbica::WriteImage< TImageType >(outputImage, outputImageFile);

    std::cout << "Resizing by a factor of " << resize << "% completed.\n";
    return EXIT_SUCCESS;
  }

  if (requestedAlgorithm == UniqueValues)
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
    return EXIT_SUCCESS;
  }

  if (requestedAlgorithm == HistogramMatching)
  {
    cbica::WriteImage< TImageType >(
      cbica::GetHistogramMatchedImage< TImageType >(
        cbica::ReadImage< TImageType >(inputImageFile), cbica::ReadImage< TImageType >(targetImageFile), histoMatchQuantiles, histoMatchBins), outputImageFile);
    std::cout << "Histogram matching completed.\n";
    return EXIT_SUCCESS;
  }

  if (requestedAlgorithm == Casting)
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
    return EXIT_SUCCESS;
  }

  if (requestedAlgorithm == TestComparison)
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

      std::cout << "Minimum Intensity Difference: " << diffFilter->GetMinimumDifference() << "\n";
      std::cout << "Maximum Intensity Difference: " << diffFilter->GetMaximumDifference() << "\n";
      std::cout << "Average Intensity Difference: " << diffFilter->GetMeanDifference() << "\n";
      std::cout << "Overall Intensity Difference: " << diffFilter->GetTotalDifference() << "\n";
      
      size_t totalSize = 1;
      for (size_t d = 0; d < TImageType::ImageDimension; d++)
      {
        totalSize *= size[d];
      }
      
      std::cout << "Number of Difference Voxels : " << diffFilter->GetNumberOfPixelsWithDifferences() << "\n";
      std::cout << "Total Voxels in Image       : " << totalSize << "\n";
      std::cout << "Percentage of Diff Voxels   : " << diffFilter->GetNumberOfPixelsWithDifferences() * 100 / totalSize << "\n";
    }
    else
    {
      std::cerr << "Images are in different spaces (size/origin/spacing mismatch).\n";
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }

  if (requestedAlgorithm == BoundingBox)
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
        for (size_t y = minFromCenter[1]; x < maxFromCenter[1]; x++)
        {
          for (size_t z = minFromCenter[2]; z < maxFromCenter[2]; x++)
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
  
  if (requestedAlgorithm == ZScoreNormalize)
  {
    auto currentDataDir = cbica::normalizePath(inputImageFile); // this is a data directory for now and will be changed
    auto outputDir = cbica::normalizePath(outputImageFile);// this is a data directory for now and will be changed
    cbica::createDir(outputDir);
    std::string hggFL_file = outputDir + "/hgg_fl.cfg", hggT1C_file = outputDir + "/hgg_t1c.cfg", hggT1_file = outputDir + "/hgg_t1.cfg", hggT2_file = outputDir + "/hgg_t2.cfg",
      lggFL_file = outputDir + "/lgg_fl.cfg", lggT1C_file = outputDir + "/lgg_t1c.cfg", lggT1_file = outputDir + "/lgg_t1.cfg", lggT2_file = outputDir + "/lgg_t2.cfg",
      hggMask_file = outputDir + "/hgg_brainMask.cfg", lggMask_file = outputDir + "/lgg_brainMask.cfg";

    auto output_HGG = outputDir + "/HGG";
    auto output_LGG = outputDir + "/LGG";
    cbica::createDir(output_HGG);
    cbica::createDir(output_LGG);
    auto allFolders = cbica::subdirectoriesInDirectory(currentDataDir, true);

    for (size_t i = 0; i < allFolders.size(); i++)
    {
      auto filesInFolder = cbica::filesInDirectory(allFolders[i]);
      for (size_t j = 0; j < filesInFolder.size(); j++)
      {
        std::ofstream file_hggFL(hggFL_file, std::ofstream::out | std::ofstream::app), file_hggT1C(hggT1C_file, std::ofstream::out | std::ofstream::app), file_hggT1(hggT1_file, std::ofstream::out | std::ofstream::app), file_hggT2(hggT2_file, std::ofstream::out | std::ofstream::app),
          file_lggFL(lggFL_file, std::ofstream::out | std::ofstream::app), file_lggT1C(lggT1C_file, std::ofstream::out | std::ofstream::app), file_lggT1(lggT1_file, std::ofstream::out | std::ofstream::app), file_lggT2(lggT2_file, std::ofstream::out | std::ofstream::app),
          file_hggBM(hggMask_file, std::ofstream::out | std::ofstream::app), file_lggBM(lggMask_file, std::ofstream::out | std::ofstream::app);
        std::string inputBase, inputExt, inputPath;
        cbica::splitFileName(filesInFolder[j], inputPath, inputBase, inputExt);
        if (filesInFolder[j].find("_seg.nii.gz") == std::string::npos)
        {
          if (inputExt == ".nii.gz")
          {
            auto inputImage = cbica::ReadImage< TImageType >(filesInFolder[j]);
            ZScoreNormalizer< TImageType > normalizer;
            normalizer.SetInputImage(inputImage);
            normalizer.Update();
            auto outputImage = normalizer.GetOutput();

            if (allFolders[i].find("HGG") != std::string::npos)
            {
              cbica::WriteImage< TImageType >(outputImage, output_HGG + "/" + inputBase + inputExt);
              if (filesInFolder[j].find("_t1ce.nii.gz") != std::string::npos)
              {
                file_hggT1C << output_HGG + "/" + inputBase + inputExt << "\n";
              }
              else if (filesInFolder[j].find("_t1.nii.gz") != std::string::npos)
              {
                file_hggT1 << output_HGG + "/" + inputBase + inputExt << "\n";
              }
              else if (filesInFolder[j].find("_flair.nii.gz") != std::string::npos)
              {
                file_hggFL << output_HGG + "/" + inputBase + inputExt << "\n";
              }
              else if (filesInFolder[j].find("_t2.nii.gz") != std::string::npos)
              {
                file_hggT2 << output_HGG + "/" + inputBase + inputExt << "\n";
              }
            }
            else
            {
              cbica::WriteImage< TImageType >(outputImage, output_LGG + "/" + inputBase + inputExt);
              if (filesInFolder[j].find("_t1ce.nii.gz") != std::string::npos)
              {
                file_lggT1C << output_LGG + "/" + inputBase + inputExt << "\n";
              }
              else if (filesInFolder[j].find("_t1.nii.gz") != std::string::npos)
              {
                file_lggT1 << output_LGG + "/" + inputBase + inputExt << "\n";
              }
              else if (filesInFolder[j].find("_flair.nii.gz") != std::string::npos)
              {
                file_lggFL << output_LGG + "/" + inputBase + inputExt << "\n";
              }
              else if (filesInFolder[j].find("_t2.nii.gz") != std::string::npos)
              {
                file_lggT2 << output_LGG + "/" + inputBase + inputExt << "\n";
              }
            }
          }
        }
        else
        {
          auto currentFile = cbica::replaceString(filesInFolder[j], "_seg.nii.gz", "t1ce.nii.gz");

          auto thresholder = itk::BinaryThresholdImageFilter< TImageType, TImageType >::New();
          thresholder->SetInput(cbica::ReadImage< TImageType >(currentFile));
          thresholder->SetLowerThreshold(1);
          thresholder->SetUpperThreshold(std::numeric_limits<float>::max());
          thresholder->SetOutsideValue(0);
          thresholder->SetInsideValue(1);
          thresholder->Update();

          if (allFolders[i].find("HGG") != std::string::npos)
          {
            cbica::WriteImage< TImageType >(thresholder->GetOutput(), output_HGG + "/brainMask_" + inputBase + inputExt);
            file_hggBM << output_HGG + "/brainMask_" + inputBase + inputExt << "\n";
          }
          else
          {
            cbica::WriteImage< TImageType >(thresholder->GetOutput(), output_LGG + "/brainMask_" + inputBase + inputExt);
            file_lggBM << output_LGG + "/brainMask_" + inputBase + inputExt << "\n";
          }
        }
        file_hggT1C.close();
        file_hggT1.close();
        file_hggT2.close();
        file_hggFL.close();
        file_lggT1C.close();
        file_lggT1.close();
        file_lggT2.close();
        file_lggFL.close();
        file_hggBM.close();
        file_lggBM.close();
      }
    }
  }

  return EXIT_SUCCESS;
}


int main(int argc, char** argv)
{
  cbica::CmdParser parser(argc, argv);

  parser.addOptionalParameter("i", "inputImage", cbica::Parameter::FILE, "NIfTI", "Input Image for processing");
  parser.addOptionalParameter("o", "outputImage", cbica::Parameter::FILE, "NIfTI", "Output Image for processing");
  parser.addOptionalParameter("r", "resize", cbica::Parameter::INTEGER, "10-500", "Resize an image based on the resizing factor given", "Example: -r 150 resizes inputImage by 150%", "Defaults to 100, i.e., no resizing");
  parser.addOptionalParameter("s", "sanityCheck", cbica::Parameter::FILE, "NIfTI Reference", "Do sanity check of inputImage with the file provided in with this parameter", "Performs checks on size, origin & spacing",
    "Pass the target image after '-s'");
  parser.addOptionalParameter("inf", "information", cbica::Parameter::BOOLEAN, "true or false", "Output the information in inputImage");
  parser.addOptionalParameter("c", "cast", cbica::Parameter::STRING, "(u)char, (u)int, (u)long, (u)longlong, float, double", "Change the input image type", "Examples: '-c uchar', '-c float', '-c longlong'");
  parser.addOptionalParameter("un", "uniqueVals", cbica::Parameter::BOOLEAN, "true or false", "Output the unique values in the inputImage", "Pass value '1' for ascending sort or '0' for no sort", "Defaults to '1'");
  parser.addOptionalParameter("b", "boundingBox", cbica::Parameter::FILE, "NIfTI Mask", "Extracts the smallest bounding box around the mask file", "With respect to inputImage", "Writes to outputImage");
  parser.addOptionalParameter("bi", "boundingIso", cbica::Parameter::BOOLEAN, "Isotropic Box or not", "Whether the bounding box is Isotropic or not", "Defaults to true");
  parser.addOptionalParameter("tb", "testBase", cbica::Parameter::FILE, "NIfTI Reference", "Baseline image to compare inputImage with");
  parser.addOptionalParameter("tr", "testRadius", cbica::Parameter::INTEGER, "0-10", "Maximum distance away to look for a matching pixel", "Defaults to 0");
  parser.addOptionalParameter("tt", "testThresh", cbica::Parameter::FLOAT, "0-5", "Minimum threshold for pixels to be different", "Defaults to 0.0");
  parser.addOptionalParameter("hi", "histoMatch", cbica::Parameter::FILE, "NIfTI Target", "Match inputImage with the file provided in with this parameter", "Pass the target image after '-h'");
  parser.addOptionalParameter("hb", "hMatchBins", cbica::Parameter::INTEGER, "1-1000", "Number of histogram bins for histogram matching", "Only used for histoMatching", "Defaults to 100");
  parser.addOptionalParameter("hq", "hMatchQnts", cbica::Parameter::INTEGER, "1-1000", "Number of quantile values to match for histogram matching", "Only used for histoMatching", "Defaults to 40");
  parser.addOptionalParameter("utB", "unitTestBuffer", cbica::Parameter::STRING, "N.A.", "Buffer test of application");
  parser.addOptionalParameter("zn", "zScoreNorm", cbica::Parameter::BOOLEAN, "N.A.", "Z-Score normalization");

  /// unit testing
  if (parser.isPresent("utB"))
  {
    char buff[100];
    sprintf(buff, "This is a testing scenario.\n");
    buff[0] = '\0';
    cbica::sleep();
    return EXIT_SUCCESS;
  }
  /// unit testing

  if (parser.isPresent("i"))
  {
    parser.getParameterValue("i", inputImageFile);
  }
  if (parser.isPresent("o"))
  {
    parser.getParameterValue("o", outputImageFile);
  }
  if (parser.isPresent("hi"))
  {
    parser.getParameterValue("hi", targetImageFile);
    requestedAlgorithm = HistogramMatching;
    if (parser.isPresent("hb"))
    {
      parser.getParameterValue("hb", histoMatchBins);
    }
    if (parser.isPresent("hq"))
    {
      parser.getParameterValue("hq", histoMatchQuantiles);
    }
  }
  if (parser.isPresent("r"))
  {
    parser.getParameterValue("r", resize);
    requestedAlgorithm = Resize;
  }
  if (parser.isPresent("s"))
  {
    parser.getParameterValue("s", targetImageFile);
    requestedAlgorithm = SanityCheck;
  }
  if (parser.isPresent("inf"))
  {
    requestedAlgorithm = Information;
  }
  if (parser.isPresent("c"))
  {
    parser.getParameterValue("c", targetImageFile);
    requestedAlgorithm = Casting;
  }
  if (parser.isPresent("un"))
  {
    parser.getParameterValue("un", uniqueValsSort);
    requestedAlgorithm = UniqueValues;
  }
  if (parser.isPresent("zn"))
  {
    requestedAlgorithm = ZScoreNormalize;
  }
  if (parser.isPresent("tb"))
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
  if (parser.isPresent("b"))
  {
    parser.getParameterValue("b", targetImageFile);
    requestedAlgorithm = BoundingBox;
    if (parser.isPresent("bi"))
    {
      parser.getParameterValue("bi", boundingBoxIsotropic);
    }
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
  if (requestedAlgorithm == ZScoreNormalize)
  {
    using ImageType = itk::Image< float, 3 >;
    return algorithmsRunner< ImageType >();
    return EXIT_SUCCESS;
  }

  auto inputImageInfo = cbica::ImageInfo(inputImageFile);

  if (requestedAlgorithm == Information)
  {
    auto dims = inputImageInfo.GetImageDimensions();
    auto size = inputImageInfo.GetImageSize();
    auto origin = inputImageInfo.GetImageOrigins();
    auto spacing = inputImageInfo.GetImageSpacings();
    auto size_string = std::to_string(size[0]);
    auto origin_string = std::to_string(origin[0]);
    auto spacing_string = std::to_string(spacing[0]);
    size_t totalSize = size[0];
    for (size_t i = 1; i < dims; i++)
    {
      size_string += "x" + std::to_string(size[i]);
      origin_string += "x" + std::to_string(origin[i]);
      spacing_string += "x" + std::to_string(spacing[i]);
      totalSize *= size[i];
    }
    std::cout << "Dimensions: " << dims << "\n";
    std::cout << "Size      : " << size_string << "\n";
    std::cout << "Total     : " << totalSize << "\n";
    std::cout << "Origin    : " << origin_string << "\n";
    std::cout << "Spacing   : " << spacing_string << "\n";
    std::cout << "Component : " << inputImageInfo.GetComponentTypeAsString() << "\n";
    std::cout << "Pixel Type: " << inputImageInfo.GetPixelTypeAsString() << "\n";
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
    return algorithmsRunner< ImageType >();

    break;
  }
  case 4:
  {
    using ImageType = itk::Image< float, 4 >;
    return algorithmsRunner< ImageType >();

    break;
  }
  default:
    std::cerr << "Supplied image has an unsupported dimension of '" << inputImageInfo.GetImageDimensions() << "'; only 2, 3 and 4 D images are supported.\n";
    return EXIT_FAILURE; // exiting here because no further processing should be done on the image
  }

  return EXIT_SUCCESS;
}