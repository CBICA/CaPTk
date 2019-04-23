#include "cbicaCmdParser.h"
#include "cbicaLogging.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"

#include "LibraPreprocess.h"

std::string inputImageFile, outputDir;

bool debugMode;

size_t resizingFactor = 100;

std::string findRelativeApplicationPath(const std::string appName)
{
  auto appName_wrap = appName;
  std::string winExt = 
#if WIN32
    winExt = ".exe";
#else
    winExt = "";
#endif
  
  if (appName.find("libra") != std::string::npos)
  {
#if WIN32
    winExt = ".bat";
#endif
  }
  auto currentApplicationPath = cbica::getExecutablePath() + "/";
  
  auto appName_path =
#ifndef __APPLE__
  cbica::normPath(currentApplicationPath + appName + winExt);
#else
  cbica::normPath(currentApplicationPath + "../Resources/bin/" + appName_wrap);
#endif  

  if (!cbica::isFile(appName_path))
  {
    std::cerr << "Please install CaPTk properly (LIBRA executable needs to be in the same location as current executable).\n";
    exit(EXIT_FAILURE);
  }
  return appName_path;
}

inline std::string getCaPTkDataDir()
{
  auto captk_currentApplicationPath = cbica::getExecutablePath();
  auto captk_dataDir = captk_currentApplicationPath + "../data/";
  if (!cbica::exists(captk_dataDir))
  {
    captk_dataDir = captk_currentApplicationPath + "../../data/";
    if (!cbica::exists(captk_dataDir))
    {
      captk_dataDir = captk_currentApplicationPath + "../Resources/data/";
      if (!cbica::exists(captk_dataDir))
      {
        std::cerr << "Data Directory not found. Please re-install CaPTk.\n";
        return "";
      }
    }
  }

  return cbica::normPath(captk_dataDir);
}

//template< class TImageType >
int algorithmsRunner()
{
  LibraPreprocess< LibraImageType > preprocessingObj;
  preprocessingObj.SetInputFileName(inputImageFile);
  preprocessingObj.SetResizingFactor(resizingFactor);
  if (debugMode)
  {
    preprocessingObj.EnableDebugMode();
  }
  preprocessingObj.Update();
  if (debugMode)
  {
    std::cout << "Done.\n";
  }

  auto outputFileName = outputDir + "/temp/" + cbica::getFilenameBase(inputImageFile) + "_preprocessed.nii.gz";
  
  cbica::WriteImage< LibraImageType >(preprocessingObj.GetOutputImage(), outputFileName);

  auto libraPath = findRelativeApplicationPath("libra");

  std::string command = libraPath + " " + inputImageFile + " " + outputDir + "/temp/" + cbica::getFilenameBase(inputImageFile) + " true true";
  std::cout << "Running LIBRA Single Image with command '" + command + "'\n";
  std::system(command.c_str());
  std::cout << "Done.\n";

  auto outputTotalMask = outputDir + "/temp/" + cbica::getFilenameBase(inputImageFile) + "/Result_Images/totalmask/totalmask.dcm";
  //auto outputTotalMaskImage = cbica::ReadImage< LibraImageType >(outputTotalMask);
  auto dicomReader = itk::ImageSeriesReader< LibraImageType >::New();
  dicomReader->SetImageIO(itk::GDCMImageIO::New());
  dicomReader->SetFileName(outputTotalMask);
  try
  {
    dicomReader->Update();
  }
  catch (itk::ExceptionObject & err)
  {
    std::cerr << "Error while loading DICOM image(s): " << err.what() << "\n";
  }
  auto outputTotalMaskImage = dicomReader->GetOutput();
  
  auto outputRelevantMaskImage = cbica::ChangeImageValues< LibraImageType >(outputTotalMaskImage, "2", "1");
  auto outputRelevantMaskImage_flipped = preprocessingObj.ApplyFlipToMaskImage(outputRelevantMaskImage);

  auto outputRelevantMaskFile = outputDir + "/temp/" + cbica::getFilenameBase(inputImageFile) + "_mask.nii.gz";
  cbica::WriteImage< LibraImageType >(outputRelevantMaskImage_flipped, outputRelevantMaskFile);

  auto featureExtractionPath = findRelativeApplicationPath("FeatureExtraction");

  auto currentDataDir = getCaPTkDataDir();
  auto latticeFeatureParamFilePath = getCaPTkDataDir() + "/featureExtraction/2_params_default_lattice.csv";
  if (!cbica::isFile(latticeFeatureParamFilePath))
  {
    std::cerr << "The default lattice parameter file, '2_params_default_lattice.csv' was not found in the data directory, '" << currentDataDir << "'; please check.\n";
    exit(EXIT_FAILURE);
  }

  command = featureExtractionPath + "-n Lattice -p " + latticeFeatureParamFilePath +
    " -o " + outputDir +
    " -i " + outputFileName + " -t MAM " +
    " -m " + outputRelevantMaskFile + " -l TT -r 1 -f 1";

  std::cout << "Running FeatureExtraction with command '" + command + "'\n";
  std::system(command.c_str());
  std::cout << "Done.\n";

  return EXIT_SUCCESS;
}


int main(int argc, char** argv)
{
  cbica::CmdParser parser(argc, argv);

  parser.addRequiredParameter("i", "inputImage", cbica::Parameter::FILE, "DICOM", "Input Image for processing");
  parser.addRequiredParameter("o", "outputDir", cbica::Parameter::DIRECTORY, "NIfTI", "Dir with write access", "All output files are written here");
  parser.addOptionalParameter("d", "debugMode", cbica::Parameter::BOOLEAN, "0 or 1", "Enabled debug mode", "Default: 0");
  parser.addOptionalParameter("r", "resize", cbica::Parameter::INTEGER, "0 - 100", "What resizing factor is to be applied", "Default: " + std::to_string(resizingFactor));

  parser.getParameterValue("i", inputImageFile);
  parser.getParameterValue("o", outputDir);

  if (parser.isPresent("d"))
  {
    parser.getParameterValue("d", debugMode);
  }
  if (parser.isPresent("r"))
  {
    parser.getParameterValue("r", resizingFactor);
  }
  auto inputImageInfo = cbica::ImageInfo(inputImageFile);

  switch (inputImageInfo.GetImageDimensions())
  {
  case 2:
  {
    //using ImageType = itk::Image< float, 2 >;
    return algorithmsRunner();

    break;
  }
  default:
    std::cerr << "Supplied image has an unsupported dimension of '" << inputImageInfo.GetImageDimensions() << "'; only 2 D images are supported.\n";
    return EXIT_FAILURE; // exiting here because no further processing should be done on the image
  }

  return EXIT_SUCCESS;
}