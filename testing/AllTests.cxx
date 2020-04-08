#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkTestingComparisonImageFilter.h"
#include "itkExceptionObject.h"

#include "cbicaCmdParser.h"
#include "cbicaUtilities.h"
#include "cbicaITKSafeImageIO.h"

#include "GeodesicSegmentation.h"
#include "EGFRvIIISurrogateIndex.h"
#include "CaPTkEnums.h"
#include "CaPTkUtils.h"

int main(int argc, char** argv)
{
  cbica::CmdParser parser = cbica::CmdParser(argc, argv);

  parser.setExeName("CaPTk_Tests");
  parser.addOptionalParameter("b", "bufferTest", cbica::Parameter::NONE, "none", "Buffer Test");
  parser.addOptionalParameter("i", "inputFile", cbica::Parameter::FILE, ".nii.gz file", "Generic input file");
  parser.addOptionalParameter("geo", "geodesic", cbica::Parameter::DIRECTORY, "Path to 'TestData/Geodesic'", "Geodesic test");
  parser.addOptionalParameter("egfr", "egfrviii", cbica::Parameter::DIRECTORY, "Path to 'TestData/EGFRvIII'", "EGFRvIII test");
  parser.addOptionalParameter("recur", "recurrence", cbica::Parameter::DIRECTORY, "Path to 'TestData/Recurrence'", "Recurrence test");

  std::string dataDir;

  if (parser.isPresent("h"))
  {
    parser.echoHelp();
    return EXIT_SUCCESS;
  }

  if (parser.isPresent("u"))
  {
    parser.echoUsage();
    return EXIT_SUCCESS;
  }

  if (parser.isPresent("v"))
  {
    parser.echoVersion();
    return EXIT_SUCCESS;
  }

  if (parser.isPresent("b"))
  {
    char buff[100];
    sprintf(buff, "This is a testing scenario.\n");
    buff[0] = '\0';
    cbica::sleep();
    return EXIT_SUCCESS;
  }

  const int numberOfPixelsTolerance = 10; // number of pixels that are acceptable to have intensity differences
  std::string inputFile, drawingFile;

  if (parser.isPresent("geodesic"))
  {
    parser.getParameterValue("geodesic", dataDir);
    if (cbica::directoryExists(dataDir))
    {
      using ImageType = itk::Image< short, 3 >;
      auto inputImage = cbica::ReadImage< ImageType >(dataDir + "input_flair.nii.gz");
      auto mask = cbica::ReadImage< ImageType >(dataDir + "_drawing.nii.gz");
      auto outputToCheck = cbica::ReadImage< ImageType >(dataDir + "input_flair_geosMap.nii.gz");

      VectorVectorDouble Indices;
      typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
      IteratorType maskIt(mask, mask->GetLargestPossibleRegion());
      maskIt.GoToBegin();
      for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++maskIt)
      {
        if (maskIt.Get() != 0) // anything which is not the background is taken as input 
        {
          VectorDouble localIndex;
          localIndex.push_back(maskIt.GetIndex()[0]);
          localIndex.push_back(maskIt.GetIndex()[1]);
          localIndex.push_back(maskIt.GetIndex()[2]);

          Indices.push_back(localIndex);
        }
      }

      GeodesicSegmentation< ImageType > segmentationClass;
      auto geosOutput = segmentationClass.Run(inputImage, Indices);

      using DiffType = itk::Testing::ComparisonImageFilter<ImageType, ImageType>;
      auto checkerType = DiffType::New();
      checkerType->SetValidInput(outputToCheck);
      checkerType->SetTestInput(geosOutput);
      checkerType->SetDifferenceThreshold(10);
      checkerType->SetToleranceRadius(1);
      checkerType->UpdateLargestPossibleRegion();

      if ((checkerType->GetTotalDifference() > 0.0) &&
        (static_cast<int>(checkerType->GetNumberOfPixelsWithDifferences()) > numberOfPixelsTolerance))
      {
        cbica::Logging(loggerFile, "Geodesic test failed with parameters: averageIntensityDifference = '"
          + std::to_string(checkerType->GetTotalDifference()) + "' and numberOfPixelsWithDifferences = '"
          + std::to_string(checkerType->GetNumberOfPixelsWithDifferences()));
        return EXIT_FAILURE;
      }

    }
    else
    {
      cbica::Logging(loggerFile, "Geodesic testing data directory not found.");
      return EXIT_FAILURE;
    }
  }

  if (parser.isPresent("egfrviii"))
  {
    parser.getParameterValue("egfrviii", dataDir);
    if (cbica::directoryExists(dataDir))
    {
      using ImageType = itk::Image< float, 3 >;
      using ImageTypePerfusion = itk::Image< float, 4 >;
      std::vector<ImageType::Pointer> Perfusion_Registered;
      std::vector<ImageType::IndexType> nearIndices, farIndices;
      auto perfusionImage = cbica::ReadImage< ImageTypePerfusion >(dataDir + "DSC-MRI_data_changed_rescaled.nii.gz");
      auto maskImage = cbica::ReadImage< ImageType >(dataDir + "Near_Far_masks_changed_rescaled.nii.gz");
      typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
      IteratorType maskIt(maskImage, maskImage->GetLargestPossibleRegion());
      for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++maskIt)
      {
        if (maskIt.Get() == 1)
          nearIndices.push_back(maskIt.GetIndex());
        else if (maskIt.Get() == 2)
          farIndices.push_back(maskIt.GetIndex());
      }

      EGFRStatusPredictor egfrEstimator;
      auto EGFRStatusParams = egfrEstimator.PredictEGFRStatus< ImageType, ImageTypePerfusion >(perfusionImage, nearIndices, farIndices);

      if ((EGFRStatusParams[0] < 0.02) || (EGFRStatusParams[0] > 0.03)) // check for PHI
      {
        cbica::Logging(loggerFile, "EGFRvIII Test failure for PHI value.");
        return EXIT_FAILURE;
      }

      float peakHeightRatio = EGFRStatusParams[1] / EGFRStatusParams[2];
      if ((peakHeightRatio < 2) || (peakHeightRatio > 3)) // peak height ratio
      {
        cbica::Logging(loggerFile, "EGFRvIII Test failure for peak height rati.o");
        return EXIT_FAILURE;
      }

      if (EGFRStatusParams[3] != 54) // near voxels used
      {
        cbica::Logging(loggerFile, "EGFRvIII Test failure for number of near voxels used.");
        return EXIT_FAILURE;
      }

      if (EGFRStatusParams[4] != 46) // far voxels used
      {
        cbica::Logging(loggerFile, "EGFRvIII Test failure for number of far voxels used.");
        return EXIT_FAILURE;
      }
    }
    else
    {
      cbica::Logging(loggerFile, "EGFRvIII testing data directory not found.");
      return EXIT_FAILURE;
    }
  }

  if (parser.isPresent("recurrence"))
  {
    std::string dataDir = DATA_DIR;
    if (cbica::directoryExists(dataDir + "CaPTk_SampleData/Recurrence"))
    {
      dataDir = dataDir + "Recurrence/sub_1";
      auto filesForSub = cbica::filesInDirectory(dataDir);
      using ImageType = itk::Image< float, 3 >;
      using ImageTypePerfusion = itk::Image< float, 4 >;

      enum DTIModalities
      {
        AD, AX, B0, FA, TR
      };
      std::map< int, ImageType::Pointer > imagesWithModalities; // keeping this int instead of string because of efficiency
      imagesWithModalities[CAPTK::ImageModalityType::IMAGE_TYPE_T1] = ImageType::New();
      imagesWithModalities[CAPTK::ImageModalityType::IMAGE_TYPE_T1CE] = ImageType::New();
      imagesWithModalities[CAPTK::ImageModalityType::IMAGE_TYPE_T2] = ImageType::New();
      imagesWithModalities[CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR] = ImageType::New();
      imagesWithModalities[CAPTK::ImageModalityType::IMAGE_TYPE_DTI + AD] = ImageType::New();
      imagesWithModalities[CAPTK::ImageModalityType::IMAGE_TYPE_DTI + AX] = ImageType::New();
      imagesWithModalities[CAPTK::ImageModalityType::IMAGE_TYPE_DTI + B0] = ImageType::New();
      imagesWithModalities[CAPTK::ImageModalityType::IMAGE_TYPE_DTI + FA] = ImageType::New();
      imagesWithModalities[CAPTK::ImageModalityType::IMAGE_TYPE_DTI + TR] = ImageType::New();
      imagesWithModalities[CAPTK::ImageModalityType::IMAGE_TYPE_DTI + TR + 1] = ImageType::New();
      // read the different modalities
      for (size_t i = 0; i < filesForSub.size(); i++)
      {
        switch (guessImageType(filesForSub[i]))
        {
        case CAPTK::ImageModalityType::IMAGE_TYPE_T1:
        {
          imagesWithModalities[CAPTK::ImageModalityType::IMAGE_TYPE_T1] = cbica::ReadImage<ImageType>(dataDir + filesForSub[i]);
          break;
        }
        case CAPTK::ImageModalityType::IMAGE_TYPE_T1CE:
        {
          imagesWithModalities[CAPTK::ImageModalityType::IMAGE_TYPE_T1CE] = cbica::ReadImage<ImageType>(dataDir + filesForSub[i]);
          break;
        }
        case CAPTK::ImageModalityType::IMAGE_TYPE_T2:
        {
          imagesWithModalities[CAPTK::ImageModalityType::IMAGE_TYPE_T2] = cbica::ReadImage<ImageType>(dataDir + filesForSub[i]);
          break;
        }
        case CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR:
        {
          imagesWithModalities[CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR] = cbica::ReadImage<ImageType>(dataDir + filesForSub[i]);
          break;
        }
        case CAPTK::ImageModalityType::IMAGE_TYPE_DTI:
        {
          // do more processing for different sub modalities - AD, AX, B0, FA, TR
          if (filesForSub[i].find("AD") != std::string::npos)
          {
            imagesWithModalities[CAPTK::ImageModalityType::IMAGE_TYPE_DTI + AD] = cbica::ReadImage<ImageType>(dataDir + filesForSub[i]);
          }
          else if (filesForSub[i].find("AX") != std::string::npos)
          {
            imagesWithModalities[CAPTK::ImageModalityType::IMAGE_TYPE_DTI + AX] = cbica::ReadImage<ImageType>(dataDir + filesForSub[i]);
          }
          else if (filesForSub[i].find("B0") != std::string::npos)
          {
            imagesWithModalities[CAPTK::ImageModalityType::IMAGE_TYPE_DTI + B0] = cbica::ReadImage<ImageType>(dataDir + filesForSub[i]);
          }
          else if (filesForSub[i].find("FA") != std::string::npos)
          {
            imagesWithModalities[CAPTK::ImageModalityType::IMAGE_TYPE_DTI + FA] = cbica::ReadImage<ImageType>(dataDir + filesForSub[i]);
          }
          else if (filesForSub[i].find("TR") != std::string::npos)
          {
            imagesWithModalities[CAPTK::ImageModalityType::IMAGE_TYPE_DTI + TR] = cbica::ReadImage<ImageType>(dataDir + filesForSub[i]);
          }
          break;
        }
        default: // this is for the drawing
        {
          if (filesForSub[i].find("drawing") != std::string::npos)
          {
            imagesWithModalities[CAPTK::ImageModalityType::IMAGE_TYPE_DTI + TR + 1] = cbica::ReadImage<ImageType>(dataDir + filesForSub[i]);
          }
          break;
        }
        }
      }

    }
    else
    {
      cbica::Logging(loggerFile, "Recurrence testing data directory not found.");
      return EXIT_FAILURE;
    }

  }

  return EXIT_SUCCESS;
}
