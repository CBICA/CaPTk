#include "cbicaCmdParser.h"
#include "cbicaLogging.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"

#include "LibraPreprocess.h"

#include "ZScoreNormalizer.h"
#include "FeatureExtraction.h"

#include <QApplication>
#include "CaPTkGUIUtils.h"
#include "ApplicationDownloadManager.h"
#include <QDebug>

std::string inputImageFile, outputDir;

bool debugMode;

size_t resizingFactor = 100;

void updateProgress(int progress)
{
    // qDebug() << "progress = " << progress << endl;
    std::cout << progress << " %" << "\t\r" << std::flush;
}

void appExit(bool isFinished)
{
    if (isFinished) {
        std::cout << "Installation successfully finished\n";
    }
    else {
        std::cout << "Installation failed\n";
    }
}

std::string findRelativeApplicationPath(const std::string appName)
{
  std::string winExt =
#if WIN32
    ".exe";
#else
    "";
#endif

  if (appName.find("libra") != std::string::npos)
  {
#if WIN32
    winExt = ".bat";
#endif
  }

  auto currentApplicationPath = cbica::normPath(cbica::getExecutablePath()) + "/";

  auto appName_path = cbica::normPath(currentApplicationPath +
#ifndef __APPLE__
    appName + winExt
#else
    "../Resources/bin/" + appName_wrap
#endif  
  );

  if (!cbica::isFile(appName_path))
  {
    std::cerr << "Please install CaPTk properly (LIBRA executable needs to be in the same location as current executable).\n";
    exit(EXIT_FAILURE);
  }
  return appName_path;
}

//template< class TImageType >
int algorithmsRunner()
{
  if (debugMode)
  {
    std::cout << "Starting pre-processing.\n";
  }
  if (!cbica::IsDicom(inputImageFile))
  {
    std::cerr << "The input image is not a DICOM image; please provide a DICOM image to continue.\n";
    return EXIT_FAILURE;
  }
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

  // auto libraPath = findRelativeApplicationPath("libra");
  //auto libraPath = cbica::normPath("C:/Projects/CaPTk_myFork/src/applications/individualApps/libra/libra.bat");
  std::string libraPath = getApplicationDownloadPath("libra");

  cbica::createDir(outputDir + "/temp");

  std::string command = libraPath + " " + inputImageFile + " " + outputDir + "/temp/" + cbica::getFilenameBase(inputImageFile)
#if WIN32
    + " true true"
#endif
    ;
  std::cout << "Running LIBRA Single Image with command '" + command + "'\n";
  std::system(command.c_str());
  auto outputTotalMask = outputDir + "/temp/" + cbica::getFilenameBase(inputImageFile) + "/Result_Images/totalmask/totalmask.dcm";
  if (cbica::isFile(outputTotalMask))
  {
    std::cout << "Done.\n";

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

    auto preprocessedImage = preprocessingObj.GetOutputImage();
    ZScoreNormalizer< LibraImageType > normalizer;
    normalizer.SetInputImage(preprocessingObj.GetOutputImage());
    normalizer.SetInputMask(outputRelevantMaskImage_flipped);
    normalizer.SetCutoffs(0, 0);
    normalizer.SetQuantiles(0, 0);
    normalizer.Update();

    auto outputFileName = outputDir + "/temp/" + cbica::getFilenameBase(inputImageFile) + "_preprocessed_normalized.nii.gz";

    cbica::WriteImage< LibraImageType >(normalizer.GetOutput(), outputFileName);

    auto outputRelevantMaskFile = outputDir + "/temp/" + cbica::getFilenameBase(inputImageFile) + "_mask.nii.gz";
    cbica::WriteImage< LibraImageType >(outputRelevantMaskImage_flipped, outputRelevantMaskFile);

    //auto featureExtractionPath = findRelativeApplicationPath("FeatureExtraction");

    auto currentDataDir = getCaPTkDataDir();
    auto latticeFeatureParamFilePath = currentDataDir + "/features/2_params_default_lattice.csv";
    if (!cbica::isFile(latticeFeatureParamFilePath))
    {
      std::cerr << "The default lattice parameter file, '2_params_default_lattice.csv' was not found in the data directory, '" << currentDataDir << "'; please check.\n";
      exit(EXIT_FAILURE);
    }

    std::cout << "Running CaPTk's FeatureExtraction.\n";

    std::vector< LibraImageType::Pointer > inputImages;
    inputImages.push_back(normalizer.GetOutput());
    FeatureExtraction< LibraImageType > features;
    features.SetPatientID("Lattice");
    features.SetInputImages(inputImages, "MAM");
    features.SetSelectedROIsAndLabels("1", "Breast");
    features.SetMaskImage(outputRelevantMaskImage_flipped);
    features.SetWriteFeatureMaps(true);
    features.SetValidMask();
    features.SetRequestedFeatures(latticeFeatureParamFilePath);
    features.SetOutputFilename(cbica::normPath(outputDir + "/features/output.csv"));
    features.SetVerticallyConcatenatedOutput(true);
    features.Update();

    std::cout << "Done.\n";

    return EXIT_SUCCESS;
  }
  else
  {
    std::cerr << "Libra did not succeed. Please recheck.\n";
    return EXIT_FAILURE;
  }
}


int main(int argc, char** argv)
{
  QApplication app(argc, argv);
  cbica::CmdParser parser(argc, argv);

  parser.addRequiredParameter("i", "inputImage", cbica::Parameter::FILE, "DICOM", "Input Image for processing");
  parser.addRequiredParameter("o", "outputDir", cbica::Parameter::DIRECTORY, "NIfTI", "Dir with write access", "All output files are written here");
  parser.addOptionalParameter("d", "debugMode", cbica::Parameter::BOOLEAN, "0 or 1", "Enabled debug mode", "Default: 0");
  parser.addOptionalParameter("r", "resize", cbica::Parameter::INTEGER, "0 - 100", "What resizing factor is to be applied", "Default: " + std::to_string(resizingFactor));

  parser.getParameterValue("i", inputImageFile);
  parser.getParameterValue("o", outputDir);

  parser.addApplicationDescription("This application uses LIBRA to extract the breast mask and then perform feature extraction");
  parser.addExampleUsage("-i C:/test/Case1.dcm -o C:/outputDir -d 1", "This command takes the input mammogram and estimates breast density");

  cbica::createDir(outputDir);
  inputImageFile = cbica::normPath(inputImageFile);
  outputDir = cbica::normPath(outputDir);

  if (parser.isPresent("d"))
  {
    parser.getParameterValue("d", debugMode);
  }
  if (parser.isPresent("r"))
  {
    parser.getParameterValue("r", resizingFactor);
  }
  //auto inputImageInfo = cbica::ImageInfo(inputImageFile);

  //switch (inputImageInfo.GetImageDimensions())
  //{
  //case 2:
  //{
    //using ImageType = itk::Image< float, 2 >;

  cbica::createDir(loggerFolder);
  cbica::createDir(downloadFolder);

  ApplicationDownloadManager appDownloadMngr;
  QObject::connect(&appDownloadMngr, &ApplicationDownloadManager::updateProgressSignal, updateProgress);
  QObject::connect(&appDownloadMngr, SIGNAL(extractResult(bool)), QApplication::instance(), SLOT(quit()));
  std::string libraPath = appDownloadMngr.getApplication("libra", true);

  if (!libraPath.empty()) { // libra is not present
    bool ret = algorithmsRunner(); 
    app.quit();
    return ret;
  }
  
  // else
  
  //  break;
  //}
  //default:
  //  std::cerr << "Supplied image has an unsupported dimension of '" << inputImageInfo.GetImageDimensions() << "'; only 2 D images are supported.\n";
  //  return EXIT_FAILURE; // exiting here because no further processing should be done on the image
  //}
  return app.exec();
}