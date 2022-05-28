#include "GeodesicSegmentation.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

#include "cbicaCmdParser.h"
#include "cbicaUtilities.h"
#include "cbicaITKSafeImageIO.h"

using VectorDouble = std::vector < double >;
int main(int argc, char **argv)
{
  auto parser = cbica::CmdParser(argc, argv, "GeodesicSegmentation");

  parser.addOptionalParameter("c", "inputCSV", cbica::Parameter::FILE, "csv File", "CSV containing list of input subjects and GLISTR/GLISTRboost labels", "If this is provided, only 'o' is needed as parameter");
  parser.addOptionalParameter("ci", "imageCSV", cbica::Parameter::STRING, "String", "CSV header of input image to process", "Required if 'c' is passed");
  parser.addOptionalParameter("cm", "imageCSV", cbica::Parameter::STRING, "String", "CSV header of mask correspoding to input image", "Required if 'c' is passed");
  //parser.addOptionalParameter("d", "dataDir", cbica::Parameter::FILE, "csv File", "CSV containing list of input subjects and GLISTR/GLISTRboost labels", "If this is provided, only 'o' is needed as parameter");
  //parser.addOptionalParameter("dr", "radiusCSV", cbica::Parameter::INTEGER, "0-10", "The eroding radius for the fused mask", "Default = 5");
  parser.addOptionalParameter("i", "image", cbica::Parameter::FILE, "NIfTI or DICOM", "Input image which needs to be segmented");
  parser.addOptionalParameter("m", "mask", cbica::Parameter::FILE, "NIfTI or DICOM", "Seed points in the image from where to start computation");
  parser.addRequiredParameter("o", "output", cbica::Parameter::STRING, "NIfTI or Directory Path", "Output File to save the results for single subject", "If 'c' is passed, then this expects a directory where results will be stored");
  parser.addOptionalParameter("n", "normalize", cbica::Parameter::BOOLEAN, "flag", "Normalize output map between 0-255 or not", "Defaults to true/1");
  parser.addOptionalParameter("t", "threshold", cbica::Parameter::INTEGER, "0-255", "Threshold distance for geodesic mask", "By default, full geodesic mask is written");
  //parser.exampleUsage("");
  parser.addExampleUsage("-i C:/input.nii.gz -m C:/mask.nii.gz -o C:/output.nii.gz",
    "The geodesic segmentation of input is calculated based on the seeded mask");
  parser.addApplicationDescription("Geodesic Distance based segmentation");

  // argc check, "u", "h" and "v" checks are done inside CmdParser

  //if (parser.isPresent("d"))
  //{
  //  std::string dataDir, dir_output;
  //  parser.getParameterValue("d", dataDir);
  //  parser.getParameterValue("o", dir_output);

  //  auto allFolders = cbica::subdirectoriesInDirectory(dataDir);

  //  for (size_t i = 0; i < allFolders.size(); i++)
  //  {
  //    using ImageType = itk::Image< float, 3 >;

  //    auto currentPath = dataDir + allFolders[i];
  //    auto allFiles = cbica::filesInDirectory(currentPath);
  //    std::string file_image, file_mask;


  //    //// check which file is which
  //    if (allFiles[0].find("_Flair.nii.gz") != std::string::npos)
  //    {
  //      file_image = currentPath + "/" + allFiles[0];
  //      file_mask = currentPath + "/" + allFiles[1];
  //    }
  //    else
  //    {
  //      file_image = currentPath + "/" + allFiles[1];
  //      file_mask = currentPath + "/" + allFiles[0];
  //    }

  //    //for (size_t j = 0; j < allFiles.size(); j++)
  //    //{
  //    //  if ((currentPath + allFiles[i]).find("_Flair.nii.gz") != std::string::npos)
  //    //  {
  //    //    file_image = currentPath + allFiles[i];
  //    //  }
  //    //  else if ((currentPath + allFiles[i]).find("_OT.nii.gz") != std::string::npos)
  //    //  {
  //    //    file_mask = currentPath + allFiles[i];
  //    //  }
  //    //}

  //    std::cout << "Reading and altering mask...\n";
  //    auto mask = cbica::ReadImage< ImageType >(file_mask);
  //    itk::ImageRegionIteratorWithIndex< ImageType > mask_orgIt(mask, mask->GetLargestPossibleRegion());
  //    for (mask_orgIt.GoToBegin(); !mask_orgIt.IsAtEnd(); ++mask_orgIt)
  //    {
  //      auto currentValue = mask_orgIt.Get();
  //      if ((currentValue == 1) || (currentValue == 2) || (currentValue == 4))
  //      {
  //        mask_orgIt.Set(1);
  //      }
  //      else
  //      {
  //        mask_orgIt.Set(0);
  //      }
  //    }

  //    // scale the input for consistency with captk
  //    auto resalceInput = itk::RescaleIntensityImageFilter< ImageType >::New();
  //    resalceInput->SetInput(cbica::ReadImage< ImageType >(file_image));
  //    resalceInput->SetOutputMinimum(0);
  //    resalceInput->SetOutputMaximum(255);
  //    resalceInput->Update();
  //    
  //    for (size_t rad = 5; rad <= 30; rad+=5)
  //    {
  //      auto currentOutput = dir_output + "/r" + std::to_string(rad);
  //      if (!cbica::directoryExists(currentOutput))
  //      {
  //        cbica::createDir(currentOutput);
  //      }

  //      // do erosion
  //      std::cout << "Doing mask erosion...\n";
  //      using StructuringElementType = itk::BinaryBallStructuringElement< ImageType::PixelType, 3 >;
  //      StructuringElementType structuringElement;
  //      structuringElement.SetRadius(rad);
  //      structuringElement.CreateStructuringElement();

  //      auto eroder = itk::BinaryErodeImageFilter< ImageType, ImageType, StructuringElementType >::New();
  //      eroder->SetInput(mask);
  //      eroder->SetRadius(rad);
  //      eroder->SetKernel(structuringElement);
  //      eroder->Update();

  //      cbica::WriteImage< ImageType >(eroder->GetOutput(), currentOutput + "/erodedMask.nii.gz");
  //      auto mask_eroded = eroder->GetOutput();

  //      std::cout << "Populating calculation indeces (1,2,4)...\n";
  //      // populate the indeces to process
  //      VectorVectorDouble IndicesToProcess;
  //      itk::ImageRegionIteratorWithIndex< ImageType > maskIt(mask_eroded, mask_eroded->GetLargestPossibleRegion());
  //      for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++maskIt)
  //      {
  //        auto currentValue = maskIt.Get();
  //        if (maskIt.Get() > 0)
  //        {
  //          auto currentIndex = maskIt.GetIndex();
  //          VectorDouble localIndex;
  //          localIndex.push_back(currentIndex[0]);
  //          localIndex.push_back(currentIndex[1]);
  //          localIndex.push_back(currentIndex[2]);
  //          IndicesToProcess.push_back(localIndex);
  //        }
  //      }

  //      if (!IndicesToProcess.empty())
  //      {
  //        std::cout << "Running Geodesic segmentation...\n";
  //        GeodesicSegmentation segmentation;
  //        auto geosOutput = segmentation.Run< ImageType >(resalceInput->GetOutput(), IndicesToProcess);
  //        cbica::WriteImage< ImageType >(geosOutput, currentOutput + "/rawgeos.nii.gz");

  //        std::cout << "Doing scaling of the Geodesic output...\n";
  //        auto rescaleFilter = itk::RescaleIntensityImageFilter< ImageType >::New();
  //        rescaleFilter->SetInput(geosOutput);
  //        rescaleFilter->SetOutputMinimum(0);
  //        rescaleFilter->SetOutputMaximum(255);
  //        rescaleFilter->Update();
  //        geosOutput = rescaleFilter->GetOutput();

  //        std::string currentOutputFile = currentOutput + "/" + cbica::getFilenameBase(file_image) + "_geos255.nii.gz";
  //        std::cout << "Writing Output File: '" + currentOutputFile + "'\n";
  //        cbica::WriteImage< ImageType >(geosOutput, currentOutputFile);
  //      }
  //      else
  //      {
  //        std::cerr << "Subject '" << allFiles[0] << "' didn't have any non-zero segmentation label.\n";
  //      }

  //    } // radius for-loop
  //  } // sub-folders for-loop
  //} // if-loop

  // check if CSV file is provided
  if (parser.isPresent("c"))
  {
    std::string file_csv, dir_output, header_image, header_mask;
    parser.getParameterValue("c", file_csv);

    parser.getParameterValue("o", dir_output);

    if (parser.isPresent("ci"))
    {
      parser.getParameterValue("ci", header_image);
    }
    else
    {
      std::cerr << "Image header is missing from CSV. Please check.\n";
      return EXIT_FAILURE;
    }

    if (parser.isPresent("cm"))
    {
      parser.getParameterValue("cm", header_mask);
    }
    else
    {
      std::cerr << "Mask header is missing from CSV. Please check.\n";
      return EXIT_FAILURE;
    }

    auto csvInfo = cbica::parseCSVFile(file_csv, header_image + "," + header_mask, "");

    for (size_t i = 0; i < csvInfo.size(); i++)
    {
      using ImageType = itk::Image< float, 3 >;

      std::string file_image = csvInfo[i].inputImages[0],
        file_mask = csvInfo[i].inputImages[1];

      std::string currentOutputFile = dir_output + "/" + cbica::getFilenameBase(file_image) + "_geos255.nii.gz";

      if (!cbica::fileExists(currentOutputFile))
      {
        auto mask = cbica::ReadImage< ImageType >(file_mask);

        // read and fuse the labels 1,2,4 into 1 and put everything else to 0
        VectorVectorDouble IndicesToProcess;
        itk::ImageRegionIteratorWithIndex< ImageType > maskIt(mask, mask->GetLargestPossibleRegion());
        for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++maskIt)
        {
          auto currentValue = maskIt.Get();
          if ((currentValue == 1) || (currentValue == 2) || (currentValue == 4))
          {
            auto currentIndex = maskIt.GetIndex();
            VectorDouble localIndex;
            localIndex.push_back(currentIndex[0]);
            localIndex.push_back(currentIndex[1]);
            localIndex.push_back(currentIndex[2]);
            IndicesToProcess.push_back(localIndex);
          }
        }

        std::cout << "Running Geodesic segmentation...\n";
        GeodesicSegmentation< ImageType > segmentation;
        auto geosOutput = segmentation.Run/*< ImageType >*/(cbica::ReadImage< ImageType >(file_image), IndicesToProcess);

        std::cout << "Doing scaling of the Geodesic output...\n";
        auto rescaleFilter = itk::RescaleIntensityImageFilter< ImageType >::New();
        rescaleFilter->SetInput(geosOutput);
        rescaleFilter->SetOutputMinimum(0);
        rescaleFilter->SetOutputMaximum(255);
        rescaleFilter->Update();
        geosOutput = rescaleFilter->GetOutput();

        std::cout << "Writing Output File: '" + currentOutputFile + "'\n";
        cbica::WriteImage< ImageType >(geosOutput, currentOutputFile);
      }
    }
  }
  // otherwise, do single subject processing
  else
  {
    std::string inputFile, drawingFile, outputFile;
    bool normalize = true;

    if (parser.isPresent("i"))
    {
      parser.getParameterValue("i", inputFile);
    }
    else
    {
      std::cerr << "Input Image is missing. Please check.\n";
      return EXIT_FAILURE;
    }

    if (parser.isPresent("m"))
    {
      parser.getParameterValue("m", drawingFile);
    }
    else
    {
      std::cerr << "Input mask is missing. Please check.\n";
      return EXIT_FAILURE;
    }

    if (parser.isPresent("o"))
    {
      parser.getParameterValue("o", outputFile);
    }
    else
    {
      std::cerr << "Output Image is missing. Please check.\n";
      return EXIT_FAILURE;
    }

    if (!cbica::fileExists(inputFile))
    {
      std::cerr << "Input image file not found :'" << inputFile << "'\n";
      return EXIT_FAILURE;
    }

    if (!cbica::fileExists(drawingFile))
    {
      std::cerr << "Input mask file not found :'" << drawingFile << "'\n";
      return EXIT_FAILURE;
    }

    auto tempExt = cbica::getFilenameExtension(outputFile, false);
    if ((tempExt != ".nii.gz") && (tempExt != ".nii"))
    {
      std::cerr << "Output file extension needs to be compatible with NIfTI format. Only '.nii.gz' or '.nii' allowed\n";
      return EXIT_FAILURE;
    }

    if (parser.isPresent("n"))
    {
      bool temp;
      parser.getParameterValue("n", temp);
      if (!temp)
      {
        normalize = false;
      }
    }

    unsigned short threshold = 0;
    if (parser.isPresent("t"))
    {
      int temp;
      parser.getParameterValue("t", temp);
      threshold = temp;
    }

    std::cout << "Reading inputs.\n";
    using ImageType = itk::Image< float, 3 >;
    using ImageTypeGeodesic = itk::Image < short, 3 >;
    auto inputImage = cbica::ReadImage< ImageType >(inputFile);
    auto mask = cbica::ReadImage< ImageType >(drawingFile);
    using CastFilterType = itk::CastImageFilter<ImageType, ImageTypeGeodesic>;
    auto castFilter = CastFilterType::New();
    castFilter->SetInput(inputImage);
    castFilter->Update();
    auto convertedInputImage = castFilter->GetOutput();
    using RescaleFilterType = itk::RescaleIntensityImageFilter<ImageTypeGeodesic, ImageTypeGeodesic>;
    auto rescaleFilter = RescaleFilterType::New();
    rescaleFilter->SetInput(convertedInputImage);
    rescaleFilter->SetOutputMinimum(0);
    rescaleFilter->SetOutputMaximum(255);
    rescaleFilter->Update();
    convertedInputImage = rescaleFilter->GetOutput();

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

    std::cout << "Performing Geodesic map estimation.\n";
    GeodesicSegmentation< ImageTypeGeodesic > segmentationClass;
    auto geosOutput = segmentationClass.Run/*< ImageType >*/(convertedInputImage, Indices);
    ImageType::PixelType maxVal;

    if (normalize)
    {
      std::cout << "Performing normalization.\n";
      maxVal = 255;
      auto rescaleFilter = itk::RescaleIntensityImageFilter< ImageTypeGeodesic >::New();
      rescaleFilter->SetInput(geosOutput);
      rescaleFilter->SetOutputMinimum(0);
      rescaleFilter->SetOutputMaximum(maxVal);
      rescaleFilter->Update();
      geosOutput = rescaleFilter->GetOutput();
    }
    else
    {
      auto calculator = itk::MinimumMaximumImageCalculator< ImageTypeGeodesic >::New();
      calculator->SetImage(geosOutput);
      calculator->ComputeMaximum();
      maxVal = calculator->GetMaximum();
    }

    if (threshold != 0)
    {
      std::cout << "Performing threshold with threshold = " << threshold << ".\n";
      auto thresholder = itk::BinaryThresholdImageFilter< ImageTypeGeodesic, ImageTypeGeodesic >::New();
      thresholder->SetInput(geosOutput);
      thresholder->SetOutsideValue(0);
      thresholder->SetInsideValue(1);
      thresholder->SetLowerThreshold(0);
      thresholder->SetUpperThreshold(threshold);
      thresholder->Update();
      geosOutput = thresholder->GetOutput();
    }

    std::cout << "Writing Output File:\n" + outputFile + "\n";
    cbica::WriteImage< ImageTypeGeodesic >(geosOutput, outputFile);
  }

  std::cout << "Finished Successfully.\n";
  return EXIT_SUCCESS;
}
