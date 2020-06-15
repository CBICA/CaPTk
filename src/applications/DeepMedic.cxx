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

std::string inputT1ce, inputT1, inputT2, inputFlair, inputMaskName, modelDirName, inputBVecName, outputDirectory, outputFile, loggerFileIn;
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
    std::cout << "Using non-CaPTk approved model at user's own risk.\n";
  }
  
  auto filesInDir = cbica::filesInDirectory(modelDirName);
  for (size_t i = 0; i < filesInDir.size(); i++)
  {
    if (filesInDir[i].find("model.ckpt.data") != std::string::npos) // find an appopriate checkpoint
    {
      break;
    }
  }

  auto t1cImg = cbica::ReadImage< TImageType >(inputT1ce);
  auto t1Img = cbica::ReadImage< TImageType >(inputT1);
  auto t2Img = cbica::ReadImage< TImageType >(inputT2);
  auto flImg = cbica::ReadImage< TImageType >(inputFlair);
  auto maskImage = cbica::CreateImage< TImageType >(t1cImg);

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

  // per-patient registration
  auto greedyExe = getApplicationPath("GreedyRegistration");
  if (!cbica::ImageSanityCheck< TImageType >(t1cImg, maskImage))
  {
    auto tempFile_input = outputDirectory + "/maskToT1gd_input.nii.gz";
    auto tempFile = outputDirectory + "/maskToT1gd.nii.gz";
    cbica::WriteImage< TImageType >(maskImage, tempFile_input);
    auto greedyCommand = greedyExe +
      " -i " + tempFile_input +
      " -f " + inputT1ce +
      " -t " + outputDirectory + "/tempMatrix.mat" +
      " -o " + tempFile + " -reg -trf -a -m MI -n 100x50x5"
      ;

    std::cout << "== Starting per-subject registration of Mask to T1-Ce using Greedy.\n";
    std::system(greedyCommand.c_str());
    maskImage = cbica::ReadImage< TImageType >(tempFile);
    std::cout << "== Done.\n";
  }
  if (!cbica::ImageSanityCheck< TImageType >(t1cImg, t1Img))
  {
    auto tempFile = outputDirectory + "/T1ToT1gd.nii.gz";
    auto greedyCommand = greedyExe +
      " -i " + inputT1 +
      " -f " + inputT1ce +
      " -t " + outputDirectory + "/tempMatrix.mat" +
      " -o " + tempFile + " -reg -trf -a -m MI -n 100x50x5"
      ;

    std::cout << "== Starting per-subject registration of T1 to T1-Ce using Greedy.\n";
    std::system(greedyCommand.c_str());
    t1Img = cbica::ReadImage< TImageType >(tempFile);
    std::cout << "== Done.\n";
  }
  if (!cbica::ImageSanityCheck< TImageType >(t1cImg, t2Img))
  {
    auto tempFile = outputDirectory + "/T2ToT1gd.nii.gz";
    auto greedyCommand = greedyExe +
      " -i " + inputT2 +
      " -f " + inputT1ce +
      " -t " + outputDirectory + "/tempMatrix.mat" +
      " -o " + tempFile + " -reg -trf -a -m MI -n 100x50x5"
      ;

    std::cout << "== Starting per-subject registration of T2 to T1-Ce using Greedy.\n";
    std::system(greedyCommand.c_str());
    t2Img = cbica::ReadImage< TImageType >(tempFile);
    std::cout << "== Done.\n";
  }
  if (!cbica::ImageSanityCheck< TImageType >(t1cImg, flImg))
  {
    auto tempFile = outputDirectory + "/FLToT1gd.nii.gz";
    auto greedyCommand = greedyExe +
      " -i " + inputFlair +
      " -f " + inputT1ce +
      " -t " + outputDirectory + "/tempMatrix.mat" +
      " -o " + tempFile + " -reg -trf -a -m MI -n 100x50x5"
      ;

    std::cout << "== Starting per-subject registration of T2-Flair to T1-Ce using Greedy.\n";
    std::system(greedyCommand.c_str());
    flImg = cbica::ReadImage< TImageType >(tempFile);
    std::cout << "== Done.\n";
  }

  if (enableNormalization)
  {
    if (inferenceType == TumorSegmentation)
    {
      std::cout << "=== Checking and rectifying (z-score) normalization status.\n";
      auto statsCalculator = itk::StatisticsImageFilter< TImageType >::New();
      statsCalculator->SetInput(t1cImg);
      statsCalculator->Update();
      if (statsCalculator->GetMean() != 0)
      {
        std::cout << "== Starting Normalization of T1CE image.\n";
        ZScoreNormalizer< TImageType > normalizer;
        normalizer.SetInputImage(t1cImg);
        if (maskProvided)
        {
          normalizer.SetInputMask(maskImage);
        }
        normalizer.SetCutoffs(cutOffLower, cutOffUpper);
        normalizer.SetQuantiles(quantLower, quantUpper);
        normalizer.Update();
        t1cImg = normalizer.GetOutput();
        std::cout << "== Done.\n";
      }

      statsCalculator->SetInput(t1Img);
      statsCalculator->Update();
      if (statsCalculator->GetMean() != 0)
      {
        std::cout << "== Starting Normalization of T1 image.\n";
        ZScoreNormalizer< TImageType > normalizer;
        normalizer.SetInputImage(t1Img);
        if (maskProvided)
        {
          normalizer.SetInputMask(maskImage);
        }
        normalizer.SetCutoffs(cutOffLower, cutOffUpper);
        normalizer.SetQuantiles(quantLower, quantUpper);
        normalizer.Update();
        t1Img = normalizer.GetOutput();
        std::cout << "== Done.\n";
      }

      statsCalculator->SetInput(t2Img);
      statsCalculator->Update();
      if (statsCalculator->GetMean() != 0)
      {
        std::cout << "== Starting Normalization of T2 image.\n";
        ZScoreNormalizer< TImageType > normalizer;
        normalizer.SetInputImage(t2Img);
        if (maskProvided)
        {
          normalizer.SetInputMask(maskImage);
        }
        normalizer.SetCutoffs(cutOffLower, cutOffUpper);
        normalizer.SetQuantiles(quantLower, quantUpper);
        normalizer.Update();
        t2Img = normalizer.GetOutput();
        std::cout << "== Done.\n";
      }

      statsCalculator->SetInput(flImg);
      statsCalculator->Update();
      if (statsCalculator->GetMean() != 0)
      {
        std::cout << "== Starting Normalization of T2-Flair image.\n";
        ZScoreNormalizer< TImageType > normalizer;
        normalizer.SetInputImage(flImg);
        if (maskProvided)
        {
          normalizer.SetInputMask(maskImage);
        }
        normalizer.SetCutoffs(cutOffLower, cutOffUpper);
        normalizer.SetQuantiles(quantLower, quantUpper);
        normalizer.Update();
        flImg = normalizer.GetOutput();
        std::cout << "== Done.\n";
      }
      std::cout << "=== Done.\n";
    }

    if (inferenceType == SkullStripping)
    {
      std::cout << "=== Starting P1P2Normalize.\n";

      P1P2Normalizer< TImageType > normalizer;
      normalizer.SetInputImage(t1cImg);
      t1cImg = normalizer.GetOutput();
      normalizer.SetInputImage(t1Img);
      t1Img = normalizer.GetOutput();
      normalizer.SetInputImage(t2Img);
      t2Img = normalizer.GetOutput();
      normalizer.SetInputImage(flImg);
      flImg = normalizer.GetOutput();
      std::cout << "=== Done.\n";
    }
  }

  if (inferenceType <= SkullStripping)
  {
    if (resamplingRes > 0)
    {
      std::cout << "=== Starting resampling of images to isotropic resolution.\n";
      t1cImg = cbica::ResampleImage< TImageType >(t1cImg, resamplingRes); // default is linear resampling to isotropic resolution of 1.0
      t1Img = cbica::ResampleImage< TImageType >(t1Img, resamplingRes); // default is linear resampling to isotropic resolution of 1.0
      flImg = cbica::ResampleImage< TImageType >(flImg, resamplingRes); // default is linear resampling to isotropic resolution of 1.0
      t2Img = cbica::ResampleImage< TImageType >(t2Img, resamplingRes); // default is linear resampling to isotropic resolution of 1.0
      maskImage = cbica::ResampleImage< TImageType >(maskImage, resamplingRes, "Nearest"); // default is linear resampling to isotropic resolution of 1.0
      std::cout << "=== Done.\n";
    }
  }

  cbica::createDir(outputDirectory);
  std::cout << "Starting DeepMedic Segmentation.\n";

  std::string file_t1ceNorm = cbica::normPath(outputDirectory + "/t1ce_normalized.nii.gz"),
    file_t1Norm = cbica::normPath(outputDirectory + "/t1_normalized.nii.gz"),
    file_t2Norm = cbica::normPath(outputDirectory + "/t2_normalized.nii.gz"),
    file_flNorm = cbica::normPath(outputDirectory + "/fl_normalized.nii.gz");
  cbica::WriteImage< TImageType >(t1cImg, file_t1ceNorm);
  cbica::WriteImage< TImageType >(t1Img, file_t1Norm);
  cbica::WriteImage< TImageType >(t2Img, file_t2Norm);
  cbica::WriteImage< TImageType >(flImg, file_flNorm);

  auto dmExe = getApplicationPath("deepMedicRun");
  //std::string dmExe = "C:/Projects/CaPTk_myFork/src/applications/individualApps/deepmedic/deepMedicRun.exe";

#ifdef _WIN32
  SetCurrentDirectory(cbica::getFilenamePath(dmExe).c_str());
#endif

  auto fullCommand = dmExe +
    " -t1 " + file_t1Norm +
    " -t1c " + file_t1ceNorm +
    " -t2 " + file_t2Norm +
    " -fl " + file_flNorm +
    " -model " + cbica::normPath(modelDirName + "/modelConfig.txt") +
    " -load " + cbica::normPath(modelDirName + "/model.ckpt") +
    " -test " + cbica::normPath(getCaPTkDataDir() + "/deepMedic/configFiles/testApiConfig.txt") +
    " -o " + outputDirectory;

  std::cout << "Running the following command:\n" << fullCommand << "\n";

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
      std::cout << "=== Performing hole-filling operation for skull stripping.\n";
      auto outputImageWithHoles = outputImage_temp;

      auto holeFiller = itk::BinaryFillholeImageFilter< TImageType >::New();
      holeFiller->SetInput(outputImageWithHoles);
      holeFiller->SetForegroundValue(1);
      holeFiller->SetFullyConnected(false);
      holeFiller->Update();

      cbica::WriteImage< TImageType >(holeFiller->GetOutput(), outputImageFile_temp);
      std::cout << "=== Done.\n";
    }
    else if (inferenceType == TumorSegmentation)
    {
      std::cout << "=== Changing the output label value from '3' to '4' for BraTS consistency.\n";
      auto outputImageWithNewValues = cbica::ChangeImageValues< TImageType >(outputImage_temp, "3", "4");

      cbica::WriteImage< TImageType >(outputImageWithNewValues, outputImageFile_temp);
    }
    
    outputImage_temp = cbica::ReadImage< TImageType >(outputImageFile_temp);

    {
      std::cout << "== Starting resampling of output segmentation back to patient space.\n";
      auto t1cImg_original = cbica::ReadImage< TImageType >(inputT1ce);
      auto resampledMask = cbica::ResampleImage< TImageType >(outputImage_temp,
        t1cImg_original->GetSpacing(),
        t1cImg_original->GetLargestPossibleRegion().GetSize(), "nearest");
      cbica::WriteImage< TImageType >(
        resampledMask,
        outputImageFile_temp
        );
      std::cout << "== Done.\n";
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

  parser.addRequiredParameter("t1c", "T1CE", cbica::Parameter::FILE, "", "The input T1CE or T1Gd image file.");
  parser.addRequiredParameter("t1", "T1", cbica::Parameter::FILE, "", "The input T1 image file.");
  parser.addRequiredParameter("fl", "FLAIR", cbica::Parameter::FILE, "", "The input T2-FLAIR image file.");
  parser.addRequiredParameter("t2", "T2", cbica::Parameter::FILE, "", "The input T2 image file.");
  parser.addOptionalParameter("m", "mask", cbica::Parameter::FILE, "", "The Optional input mask file.", "This is needed for normalization only");
  parser.addRequiredParameter("md", "modelDir", cbica::Parameter::DIRECTORY, "", "The trained model to use", "See examples in '${CaPTk_installDir}/data/deepMedic/saved_models'");
  parser.addRequiredParameter("o", "output", cbica::Parameter::DIRECTORY, "", "The output directory");
  parser.addOptionalParameter("zn", "zScoreNorm", cbica::Parameter::BOOLEAN, "N.A.", "Z-Score normalization", "Set to '0' if you are passing normalized images");
  parser.addOptionalParameter("zq", "zNormQuant", cbica::Parameter::FLOAT, "0-100", "The Lower-Upper Quantile range to remove", "Default: " + std::to_string(quantLower) + "," + std::to_string(quantUpper));
  parser.addOptionalParameter("zc", "zNormCut", cbica::Parameter::FLOAT, "0-10", "The Lower-Upper Cut-off (multiple of stdDev) to remove", "Default: " + std::to_string(cutOffLower) + "," + std::to_string(cutOffUpper));
  parser.addOptionalParameter("rr", "resizeResolution", cbica::Parameter::FLOAT, "0-10", "[Resample] Isotropic resampling resolution to change to", "Defaults to " + std::to_string(resamplingRes), "If '0' is passed, no resampling is done.");

  parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");
  parser.addApplicationDescription("This is a Deep Learning based inference engine based on DeepMedic (see documentation for details)");
  parser.addExampleUsage("-t1 c:/t1.nii.gz -t1c c:/t1gc.nii.gz -t2 c:/t2.nii.gz -fl c:/fl.nii.gz -o c:/output -md c:/CaPTk_install/data/deepMedic/saved_models/skullStripping", "This does a skull stripping of the input structural data");
  parser.addExampleUsage("-t1 c:/t1.nii.gz -t1c c:/t1gc.nii.gz -t2 c:/t2.nii.gz -fl c:/fl.nii.gz -o c:/output -md c:/CaPTk_install/data/deepMedic/saved_models/brainTumorSegmentation", "This does a tumor segmentation of the input structural data");

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
  parser.getParameterValue("t1c", inputT1ce);
  parser.getParameterValue("t1", inputT1);
  parser.getParameterValue("t2", inputT2);
  parser.getParameterValue("fl", inputFlair);

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

  auto imageInfo = cbica::ImageInfo(inputT1ce);
  
  if (!cbica::ImageSanityCheck(inputT1, inputT1ce))
  {
    std::cerr << "T1 and T1CE images are in inconsistent spaces, please register them before trying.\n";
    return EXIT_FAILURE;
  }
  if (!cbica::ImageSanityCheck(inputT2, inputT1ce))
  {
    std::cerr << "T2 and T1CE images are in inconsistent spaces, please register them before trying.\n";
    return EXIT_FAILURE;
  }
  if (!cbica::ImageSanityCheck(inputFlair, inputT1ce))
  {
    std::cerr << "T2-Flair and T1CE images are in inconsistent spaces, please register them before trying.\n";
    return EXIT_FAILURE;
  }

  if (maskProvided)
  {
    if(!cbica::ImageSanityCheck(inputT1ce, inputMaskName))
    {
      std::cerr << "The input image(s) and mask are in inconsistent spaces, please register them before trying.\n";
      return EXIT_FAILURE;
    }
  }

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