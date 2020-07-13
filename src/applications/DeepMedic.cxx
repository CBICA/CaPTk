#include "ZScoreNormalizer.h"
#include "P1P2Normalizer.h"
#include "cbicaUtilities.h"
#include "cbicaCmdParser.h"
#include "cbicaLogging.h"
#include "cbicaITKUtilities.h"
#include "cbicaITKImageInfo.h"
#include "cbicaITKSafeImageIO.h"
#include "CaPTkGUIUtils.h"

#include "itkStatisticsImageFilter.h"
#include "itkBinaryFillholeImageFilter.h"
#include "itkJoinSeriesImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkPermuteAxesImageFilter.h"

#ifdef _WIN32
#include <direct.h>
#endif

//#include "CAPTk.h"

std::string inputImageFiles, inputMaskName, modelDirName, outputDirectory, outputFile, loggerFileIn;
std::vector< std::string > inputImageFilesVector;
float quantLower = 5, quantUpper = 95, cutOffLower = 3, cutOffUpper = 3;
bool maskProvided = false, debugMode = false, enableNormalization = true;
int inferenceType = 0;
float resamplingRes = 1;

enum InferenceTypes
{
  TumorSegmentation,
  SkullStripping,
  MaxType
};

template< class TImageType >
typename TImageType::Pointer HoleFillForSingleAxis(typename TImageType::Pointer input, size_t axisToIterate)
{
  auto size = input->GetLargestPossibleRegion().GetSize();
  auto origin = input->GetOrigin();
  typename TImageType::IndexType regionIndex;
  typename TImageType::SizeType regionSize;

  regionSize = size;
  regionSize[axisToIterate] = 0;
  regionIndex.Fill(0);

  itk::FixedArray< unsigned int, TImageType::ImageDimension > order;
  order[0] = 0;
  order[1] = 0;
  order[axisToIterate] = 1;
  auto permuter = itk::PermuteAxesImageFilter< TImageType >::New();
  permuter->SetInput(input);
  permuter->SetOrder(order);
  permuter->Update();
  cbica::WriteImage< TImageType >(permuter->GetOutput(), 
    outputDirectory + "/segm_permuter" + std::to_string(axisToIterate) + ".nii.gz");

  using TImageType2D = itk::Image< typename TImageType::PixelType, 2 >;
  auto extractor = itk::ExtractImageFilter< TImageType, TImageType2D >::New();
  extractor->SetInput(input);
  extractor->SetDirectionCollapseToIdentity(); // This is required.

  //auto joiner = itk::JoinSeriesImageFilter< itk::Image<float,2>, itk::Image<float, 3> >::New();
  auto joiner = itk::JoinSeriesImageFilter< TImageType2D, TImageType >::New();
  joiner->SetOrigin(input->GetOrigin()[axisToIterate]);
  joiner->SetSpacing(input->GetSpacing()[axisToIterate]);

  for (size_t i = 0; i < size[axisToIterate]; i++)
  {
    regionIndex[axisToIterate] = i;
    typename TImageType::RegionType desiredRegion(regionIndex, regionSize);
    extractor->SetExtractionRegion(desiredRegion);
    extractor->Update();

    auto debugImage = extractor->GetOutput();

    auto holeFiller = itk::BinaryFillholeImageFilter< TImageType2D >::New();
    holeFiller->SetInput(extractor->GetOutput());
    holeFiller->SetForegroundValue(1);
    holeFiller->Update();

    joiner->SetInput(i, holeFiller->GetOutput());
  }

  joiner->Update();
  return joiner->GetOutput();
}

template< class TImageType >
typename TImageType::Pointer HoleFillOnThreeAxes(typename TImageType::Pointer input)
{
  auto holeFiller = itk::BinaryFillholeImageFilter< TImageType >::New();
  holeFiller->SetInput(input);
  holeFiller->SetForegroundValue(1);
  holeFiller->SetFullyConnected(true);
  return holeFiller->GetOutput();

  auto output = HoleFillForSingleAxis< TImageType >(input, 2);
  if (cbica::ImageSanityCheck<TImageType>(input, output))
  {
    output = HoleFillForSingleAxis< TImageType >(output, 1);
    if (cbica::ImageSanityCheck<TImageType>(input, output))
    {
      output = HoleFillForSingleAxis< TImageType >(output, 0);
      if (cbica::ImageSanityCheck<TImageType>(input, output))
      {
        return output;
      }
    }
  }
  
  std::cerr << "Something went wrong with hole filling, please check.\n";
  exit(EXIT_FAILURE);
}

template< class TImageType >
void algorithmRunner()
{
  if (!cbica::isFile(modelDirName + "/modelConfig.txt"))
  {
    std::cerr << "'modelConfig.txt' was not found in the directory, please check.\n";
    return;
  }
  if (!cbica::isFile(modelDirName + "/model.ckpt.index"))
  {
    std::cerr << "'model.ckpt' was not found in the directory, please check.\n";
    return;
  }
  if (cbica::isFile(modelDirName + "/VERSION.yaml"))
  {
    if (!cbica::IsCompatible(modelDirName + "/VERSION.yaml"))
    {
      std::cerr << "The version of model is incompatible with this version of CaPTk.\n";
      return;
    }
  }
  else
  {
    std::cout << "!!! WARNING: Using non-CaPTk approved model at user's own risk...\n";
  }
  
  auto filesInDir = cbica::filesInDirectory(modelDirName);
  for (size_t i = 0; i < filesInDir.size(); i++)
  {
    if (filesInDir[i].find("model.ckpt.data") != std::string::npos) // find an appopriate checkpoint
    {
      break;
    }
  }

  std::vector< typename TImageType::Pointer > inputImages;
  inputImages.resize(inputImageFilesVector.size());

  for (size_t i = 0; i < inputImageFilesVector.size(); i++)
  {
    inputImages[i] = cbica::ReadImage< TImageType >(inputImageFilesVector[i]);
  }

  auto maskImage = cbica::CreateImage< TImageType >(inputImages[0]);

  if (maskProvided)
  {
    maskImage = cbica::ReadImage< TImageType >(inputMaskName);
  }

  // TBD: this requires cleanup
  if (modelDirName.find("tumor") != std::string::npos)
  {
    inferenceType = 0;
  }
  else if (modelDirName.find("skull") != std::string::npos)
  {
    inferenceType = 1;
  }

  if (enableNormalization)
  {
    if (inferenceType == TumorSegmentation)
    {
      std::cout << "=== Checking and rectifying (z-score) normalization status.\n";

      for (size_t i = 0; i < inputImages.size(); i++)
      {
        auto statsCalculator = itk::StatisticsImageFilter< TImageType >::New();
        statsCalculator->SetInput(inputImages[i]);
        statsCalculator->Update();
        if (statsCalculator->GetMean() != 0)
        {
          if (debugMode)
          {
            std::cout << "== Starting Normalization of image '" << i << "'.\n";
          }
          ZScoreNormalizer< TImageType > normalizer;
          normalizer.SetInputImage(inputImages[i]);
          if (maskProvided)
          {
            normalizer.SetInputMask(maskImage);
          }
          normalizer.SetCutoffs(cutOffLower, cutOffUpper);
          normalizer.SetQuantiles(quantLower, quantUpper);
          normalizer.Update();
          inputImages[i] = normalizer.GetOutput();

          if (debugMode)
          {
            std::cout << "== Done.\n";
          }
        }
      }
      std::cout << "=== Done.\n";
    }

    if (inferenceType == SkullStripping)
    {
      std::cout << "=== Starting P1P2Normalize.\n";

      for (size_t i = 0; i < inputImages.size(); i++)
      {
        P1P2Normalizer< TImageType > normalizer;
        normalizer.SetInputImage(inputImages[i]);
        inputImages[i] = normalizer.GetOutput();
      }
      std::cout << "=== Done.\n";
    }
  }

  if (inferenceType <= SkullStripping) // this is present as a sanity check in case we design a model that does not use isotropic resampling, which is unlikely at this point
  {
    if (resamplingRes > 0)
    {
      std::cout << "=== Starting resampling of images to isotropic resolution.\n";

      for (size_t i = 0; i < inputImages.size(); i++)
      {
        inputImages[i] = cbica::ResampleImage< TImageType >(inputImages[i], resamplingRes);
      }
      maskImage = cbica::ResampleImage< TImageType >(maskImage, resamplingRes, "Nearest"); // default is linear resampling to isotropic resolution of 1.0
      std::cout << "=== Done.\n";
    }
  }

  cbica::createDir(outputDirectory);
  std::cout << "Starting DeepMedic Segmentation.\n";

  std::vector< std::string > file_inputNormalized;
  file_inputNormalized.resize(inputImages.size());
  std::string files_forCommand;

  for (size_t i = 0; i < inputImages.size(); i++)
  {
    file_inputNormalized[i] = cbica::normPath(outputDirectory + "/" + std::to_string(i) + "_normalized.nii.gz");
    cbica::WriteImage< TImageType >(inputImages[i], file_inputNormalized[i]);
    files_forCommand += file_inputNormalized[i] + ",";
  }

  files_forCommand.pop_back(); // remove last ","

  auto dmExe = getApplicationPath("deepMedicInference");
  //std::string dmExe = "C:/Projects/CaPTk_myFork/src/applications/individualApps/deepmedic/deepMedicRun.exe"; // FOR DEBUGGING ONLY //

#ifdef _WIN32
  SetCurrentDirectory(cbica::getFilenamePath(dmExe).c_str());
#endif

  // order for 4-modality applications: t1,t1c,t2,fl

  auto fullCommand = dmExe +
    " -i " + files_forCommand +
    " -model " + cbica::normPath(modelDirName + "/modelConfig.txt") +
    " -load " + cbica::normPath(modelDirName + "/model.ckpt") +
    " -test " + cbica::normPath(getCaPTkDataDir() + "/deepMedic/configFiles/testApiConfig.txt") +
    " -o " + outputDirectory;

  if (debugMode)
  {
    std::cout << "Running the following command:\n" << fullCommand << "\n";
  }

  if (std::system(fullCommand.c_str()) != 0)
  {
    std::cerr << "DeepMedic exited with code !=0.\n";
    exit(EXIT_FAILURE);
  }

  auto outputImageFile_temp = outputDirectory + "/predictions/testApiSession/predictions/Segm.nii.gz";
  if (cbica::exists(outputImageFile_temp))
  {
    auto outputImage_temp = cbica::ReadImage< TImageType >(outputImageFile_temp);
    if (inferenceType == SkullStripping)
    {
      if (debugMode)
      {
        std::cout << "=== Performing hole-filling operation for skull stripping.\n";
      }
      auto outputImageWithHoles = outputImage_temp;

      auto holeFiller = itk::BinaryFillholeImageFilter< TImageType >::New();
      holeFiller->SetInput(outputImageWithHoles);
      holeFiller->SetForegroundValue(1);
      holeFiller->SetFullyConnected(false);
      holeFiller->Update();

      cbica::WriteImage< TImageType >(holeFiller->GetOutput(), outputImageFile_temp);
      if (debugMode)
      {
        std::cout << "=== Done.\n";
      }
    }
    else if (inferenceType == TumorSegmentation)
    {
      if (debugMode)
      {
        std::cout << "=== Changing the output label value from '3' to '4' for BraTS consistency.\n";
      }
      auto outputImageWithNewValues = cbica::ChangeImageValues< TImageType >(outputImage_temp, "3", "4");

      cbica::WriteImage< TImageType >(outputImageWithNewValues, outputImageFile_temp);
    }
    
    if (!outputFile.empty())
    {
      cbica::WriteImage< TImageType >(
        cbica::ReadImage< TImageType >(outputImageFile_temp),
        outputFile
        );
    }
  }

  return;
}

int main(int argc, char **argv)
{
  cbica::CmdParser parser(argc, argv, "DeepMedic");

  parser.addRequiredParameter("i", "inputImages", cbica::Parameter::STRING, "NIfTI files", "Input images provided in a comma-separated manner", "Should be in the same order as trained model", "Example: '-i /path/t1.nii.gz,/path/t1c.nii.gz'");
  parser.addOptionalParameter("m", "mask", cbica::Parameter::FILE, "", "The Optional input mask file.", "This is needed for normalization only");
  parser.addRequiredParameter("md", "modelDir", cbica::Parameter::DIRECTORY, "", "The trained model to use", "See examples in '${CaPTk_installDir}/data/deepMedic/saved_models'");
  parser.addRequiredParameter("o", "output", cbica::Parameter::DIRECTORY, "", "The output directory");
  parser.addOptionalParameter("zn", "zScoreNorm", cbica::Parameter::BOOLEAN, "N.A.", "Z-Score normalization", "Set to '0' if you are passing normalized images");
  parser.addOptionalParameter("zq", "zNormQuant", cbica::Parameter::FLOAT, "0-100", "The Lower-Upper Quantile range to remove", "Default: " + std::to_string(quantLower) + "," + std::to_string(quantUpper));
  parser.addOptionalParameter("zc", "zNormCut", cbica::Parameter::FLOAT, "0-10", "The Lower-Upper Cut-off (multiple of stdDev) to remove", "Default: " + std::to_string(cutOffLower) + "," + std::to_string(cutOffUpper));
  parser.addOptionalParameter("rr", "resizeResolution", cbica::Parameter::FLOAT, "0-10", "[Resample] Isotropic resampling resolution to change to", "Defaults to " + std::to_string(resamplingRes), "If '0' is passed, no resampling is done.");
  parser.addOptionalParameter("d", "debugModel", cbica::Parameter::BOOLEAN, "0-1", "Enable/disable debug mode for extra information on console");

  parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");
  parser.addApplicationDescription("This is a Deep Learning based inference engine based on DeepMedic (see documentation for details)");
  parser.addExampleUsage("-i c:/t1.nii.gz,c:/t1gc.nii.gz,c:/t2.nii.gz,c:/fl.nii.gz -o c:/output -md c:/CaPTk_install/data/deepMedic/saved_models/brainTumorSegmentation", "This does a tumor segmentation of the input structural data");
  parser.addExampleUsage("-i c:/t1.nii.gz,c:/t1gc.nii.gz,c:/t2.nii.gz,c:/fl.nii.gz -o c:/output -md c:/CaPTk_install/data/deepMedic/saved_models/skullStripping", "This does a skull stripping of the input structural data");
  parser.addExampleUsage("-i c:/t1_withSkull.nii.gz -o c:/output -md c:/CaPTk_install/data/deepMedic/saved_models/skullStripping_modalityAgnostic", "This does a skull-stripping using a single modality (can use either T1,T2-Gd,T2 or Flair)");

  // parameters to get from the command line
  cbica::Logging logger;

  if (parser.isPresent("L"))
  {
    parser.getParameterValue("L", loggerFileIn);
    logger.UseNewFile(loggerFileIn);
  }

  if (parser.isPresent("d"))
  {
    parser.getParameterValue("d", debugMode);
  }
  parser.getParameterValue("i", inputImageFiles);
  inputImageFilesVector = cbica::stringSplit(inputImageFiles, ",");
  if (inputImageFilesVector.size() == 1)
  {
    // just in case the user passes input files with a different delimiter
    inputImageFilesVector = cbica::stringSplit(inputImageFiles, "|");
    if (inputImageFilesVector.size() == 1)
    {
      // at this point, we assume the user wants to pass a single image and proceed
    }
  }
  
  parser.getParameterValue("o", outputDirectory);  
  // sanity check in case the user has passed a file
  if (!cbica::getFilenameExtension(outputDirectory, false).empty())
  {
    outputFile = outputDirectory;
    outputDirectory = cbica::getFilenamePath(outputFile, false);
  }

  if (parser.isPresent("m"))
  {
    parser.getParameterValue("m", inputMaskName);
    maskProvided = true;
  }

  if (parser.isPresent("md"))
  {
    parser.getParameterValue("md", modelDirName);
  }

  if (parser.isPresent("zn"))
  {
    parser.getParameterValue("zn", enableNormalization);
  }
  if (parser.isPresent("zc"))
  {
    std::string tempCutOff;
    parser.getParameterValue("zc", tempCutOff);
    auto temp = cbica::stringSplit(tempCutOff, ",");
    if (temp.size() == 2)
    {
      cutOffLower = std::atof(temp[0].c_str());
      cutOffUpper = std::atof(temp[1].c_str());

      if (cutOffUpper < cutOffLower)
      {
        std::swap(cutOffUpper, cutOffLower);
      }
    }
  }
  if (parser.isPresent("zq"))
  {
    std::string tempCutOff;
    parser.getParameterValue("zq", tempCutOff);
    auto temp = cbica::stringSplit(tempCutOff, ",");
    if (temp.size() == 2)
    {
      quantLower = std::atof(temp[0].c_str());
      quantUpper = std::atof(temp[1].c_str());

      if (quantUpper < quantLower)
      {
        std::swap(quantUpper, quantLower);
      }
    }
  }

  if (parser.isPresent("rr"))
  {
    parser.getParameterValue("rr", resamplingRes);
  }
  
  if (parser.isPresent("t"))
  {
    parser.getParameterValue("t", inferenceType);
  }

  if (!cbica::isFile(inputImageFilesVector[0]))
  {
    std::cerr << "Could not find the image file '" << inputImageFilesVector[0] << "' in the filesystem.\n";
    return EXIT_FAILURE;
  }
  // perform sanity checks of input images only if multiple images are passed
  if (inputImageFilesVector.size() > 1)
  {
    for (size_t i = 1; i < inputImageFilesVector.size(); i++)
    {
      if (!cbica::ImageSanityCheck(inputImageFilesVector[0], inputImageFilesVector[i]))
      {
        std::cerr << "The first image file '" << inputImageFilesVector[0] << "' has different physical dimensions with '" << inputImageFilesVector[i] << "', please register them before trying again.\n";
        return EXIT_FAILURE;
      }
    }
  }
  
  if (maskProvided)
  {
    if(!cbica::ImageSanityCheck(inputImageFilesVector[0], inputMaskName))
    {
      std::cerr << "The input image(s) and mask are in inconsistent spaces, please register them before trying.\n";
      return EXIT_FAILURE;
    }
  }

  auto imageInfo = cbica::ImageInfo(inputImageFilesVector[0]);

  switch (imageInfo.GetImageDimensions())
  {
  case 3:
  {
    using ImageType = itk::Image< float, 3 >;
    algorithmRunner< ImageType >();

    break;
  }
  default:
    std::cerr << "Only 2D images are supported right now.\n";
    return EXIT_FAILURE;
  }

  std::cout << "Finished successfully.\n";
  return EXIT_SUCCESS;
}