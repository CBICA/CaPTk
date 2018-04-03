#include <iostream>
#include <iterator>
#include <stdlib.h>
#include <string>
#include <algorithm>
#include <exception>
#include <typeinfo>
#include <stdexcept>
#include <stdio.h>
#include <numeric>
#include <functional>
#include <vector>

#include "cbicaCmdParser.h"
#include "cbicaLogging.h"

#include "Joint_Segm.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaUtilities.h"
//#include "cbicaPreProcessImage.h"
#include "itkImageFileReader.h"
#include "itkScalarImageToRunLengthFeaturesFilter.h"
#include "itkHistogramToTextureFeaturesFilter.h"
#include "itkVectorImage.h"
#include "itkDenseFrequencyContainer2.h"
#include "itkScalarImageToCooccurrenceMatrixFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVectorContainer.h"
#include "itkStatisticsImageFilter.h"
//#include "opencv/cv.h"
#include "opencv2/core/core.hpp"
#include "itkShiftScaleImageFilter.hxx"
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkConstantPadImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryMask3DMeshSource.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkBinaryImageToLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelSelectionLabelMapFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkOtsuMultipleThresholdsImageFilter.h"
#include "itkPasteImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkNormalizeImageFilter.h"

//! standard definitions 
typedef itk::Image< float, 3 > ImageType;
typedef ImageType::Pointer InputImagePointerType;
typedef itk::ImageFileReader< ImageType> ImageReaderType;
typedef ImageReaderType::Pointer InputImageFileReaderPointerType;
using namespace std;

int main(int argc, char** argv)
{
  cbica::CmdParser parser(argc, argv);

  parser.addOptionalParameter("d", "dataDir", cbica::Parameter::DIRECTORY, "none", "Parent directory where all data is present", "Not required if individual files are passed");
  parser.addRequiredParameter("o", "outputDir", cbica::Parameter::DIRECTORY, "none", "Output directory where all data is written", "Not required if individual files are passed");
  parser.addOptionalParameter("c", "ctFile", cbica::Parameter::FILE, "NIfTI", "CT image of the lung", "Not needed if '-d' is being used");
  parser.addOptionalParameter("p", "petFile", cbica::Parameter::FILE, "NIfTI", "PET image of the lung", "Not needed if '-d' is being used");
  parser.addOptionalParameter("m", "maskOut", cbica::Parameter::FILE, "NIfTI", "Mask output of the lung", "Not needed if '-d' is being used");
  //parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");

  cbica::Logging logger;
  std::string loggerFile;
  std::string datadirectory, ctImageFileName, petImageFileName, outputDir, maskImageFileName;
  bool loggerRequested = false;

  Joint_Segm segmenation;

  if (parser.isPresent("L"))
  {
    parser.getParameterValue("L", loggerFile);
    loggerRequested = true;
    logger.UseNewFile(loggerFile);
  }

  if (!cbica::isDir(outputDir))
  {
    cbica::createDir(outputDir);
  }

  if (parser.isPresent("d"))
  {
    // dataDir has been provided
    parser.getParameterValue("d", datadirectory);
    // ensure dataDir ends with a '/
    if (datadirectory[datadirectory.length() - 1] != '/')
    {
      datadirectory.append("/");
    }

    if (parser.isPresent("o"))
    {
      parser.getParameterValue("o", outputDir);
    }

    if (outputDir.empty())
    {
      logger.WriteError("OutputDir cannot be empty.");
      return EXIT_FAILURE;
    }

    if (!cbica::isDir(datadirectory))
    {
      logger.WriteError("The data directory is not valid");
      return EXIT_FAILURE;
    }

    std::vector< std::string > detectedSubjects = cbica::subdirectoriesInDirectory(datadirectory);
    if (detectedSubjects.empty())
    {
      std::cerr << "No images were found in the data directory '" + datadirectory + "'\n";
      return EXIT_FAILURE;
    }
    for (size_t subject = 0; subject < detectedSubjects.size(); subject++)
    {
      const std::string subjectUnderConsideration = datadirectory + detectedSubjects[subject] + "/";
      std::vector< std::string > filesPerSubject = cbica::filesInDirectory(subjectUnderConsideration, false);
      InputImagePointerType imgin = ImageType::New();
      InputImagePointerType PET_imgin = ImageType::New();


      for (size_t j = 0; j < filesPerSubject.size(); j++)
      {
        FileNameParts  fileUnderConsideration = FileNameParts(subjectUnderConsideration + filesPerSubject[j]);
        if (!cbica::fileExists(fileUnderConsideration.fullFileName))
        {
          //logger.WriteError("Aborting:   Could not read " + fileUnderConsideration.fullFileName);
          std::cerr << "Could not read the file: " + fileUnderConsideration.fullFileName << "\n";
          return EXIT_FAILURE;
        }
        if (fileUnderConsideration.base == "crop_Ct2Pet")
        {
          if (!cbica::fileExists(fileUnderConsideration.fullFileName))
          {
            logger.WriteError("Aborting:   Could not read " + fileUnderConsideration.fullFileName);
            return EXIT_FAILURE;
          }
          else
            imgin = cbica::ReadImage(fileUnderConsideration.fullFileName);
        }
        if (fileUnderConsideration.base == "crop_Pet")
        {
          if (!cbica::fileExists(fileUnderConsideration.fullFileName))
          {
            logger.WriteError("Aborting:   Could not read " + fileUnderConsideration.fullFileName);
            return EXIT_FAILURE;
          }
          PET_imgin = cbica::ReadImage(fileUnderConsideration.fullFileName);
        }

      }
      logger.Write("Images Read");
      maskImageFileName = outputDir + "/" + subjectUnderConsideration + "_mask.nii.gz";
      auto maskImage = segmenation.generateseeds< ImageType >(imgin, PET_imgin);
      cbica::WriteImage< ImageType >(maskImage, maskImageFileName);
    } // subjects for-loop
  } // dataDir if-loop
  else
  {
    // check for individual file types
    if (parser.isPresent("c"))
    {
      parser.getParameterValue("c", ctImageFileName);
    }
    if (parser.isPresent("p"))
    {
      parser.getParameterValue("p", petImageFileName);
    }
    if (parser.isPresent("m"))
    {
      parser.getParameterValue("m", maskImageFileName);
    }

    if (ctImageFileName.empty() || petImageFileName.empty() || maskImageFileName.empty())
    {
      logger.WriteError("One of the required files is missing. Please check");
      return EXIT_FAILURE;
    }

    if (!cbica::fileExists(ctImageFileName))
    {
      logger.WriteError("Could not read the file '" + ctImageFileName + "'");
      return EXIT_FAILURE;
    }
    if (!cbica::fileExists(petImageFileName))
    {
      logger.WriteError("Could not read the file '" + petImageFileName + "'");
      return EXIT_FAILURE;
    }
    //if (!cbica::fileExists(maskImageFileName))
    //{
    //  logger.WriteError("Could not read the file '" + maskImageFileName + "'");
    //  return EXIT_FAILURE;
    //}

    cout << " Reading input files..." << endl;
    InputImagePointerType imgin = cbica::ReadImage(ctImageFileName);
    InputImagePointerType PET_imgin = cbica::ReadImage(petImageFileName);
    // InputImagePointerType mask = cbica::ReadImage(maskImageFileName);

    cout << " Generating segmentation....." << endl;
    auto maskImage = segmenation.generateseeds< ImageType >(imgin, PET_imgin);
    cbica::WriteImage< ImageType >(maskImage, maskImageFileName);

  }
  cout << "All Done!" << endl;
  return EXIT_SUCCESS;
}