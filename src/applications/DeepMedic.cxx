#include "ZScoreNormalizer.h"
#include "cbicaUtilities.h"
#include "cbicaCmdParser.h"
#include "cbicaLogging.h"
#include "cbicaITKUtilities.h"
#include "cbicaITKImageInfo.h"
#include "cbicaITKSafeImageIO.h"
#include "CaPTkGUIUtils.h"

#ifdef _WIN32
#include <direct.h>
#endif

//#include "CAPTk.h"

std::string inputT1ce, inputT1, inputT2, inputFlair, inputMaskName, modelDirName, inputBVecName, outputDirectory, loggerFileIn;
float quantLower = 5, quantUpper = 95, cutOffLower = 3, cutOffUpper = 3;
bool maskProvided = false;

template< class TImageType >
void algorithmRunner()
{
  auto t1cImg = cbica::ReadImage< TImageType >(inputT1ce);
  auto t1Img = cbica::ReadImage< TImageType >(inputT1);
  auto t2Img = cbica::ReadImage< TImageType >(inputT2);
  auto flImg = cbica::ReadImage< TImageType >(inputFlair);
  auto maskImage = cbica::CreateImage< TImageType >(t1cImg);
  if (maskProvided)
  {
    maskImage = cbica::ReadImage< TImageType >(inputMaskName);
  }

  // per-patient registration
  auto greedyExe = getApplicationPath("GreedyRegistration");
  if (!cbica::ImageSanityCheck< TImageType >(t1cImg, maskProvided))
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

    std::cout << "Starting per-subject registration of Mask to T1-Ce using Greedy.\n";
    std::system(greedyCommand.c_str());
    maskImage = cbica::ReadImage< TImageType >(tempFile);
    std::cout << "Done.\n";
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

    std::cout << "Starting per-subject registration of T1 to T1-Ce using Greedy.\n";
    std::system(greedyCommand.c_str());
    t1Img = cbica::ReadImage< TImageType >(tempFile);
    std::cout << "Done.\n";
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

    std::cout << "Starting per-subject registration of T2 to T1-Ce using Greedy.\n";
    std::system(greedyCommand.c_str());
    t2Img = cbica::ReadImage< TImageType >(tempFile);
    std::cout << "Done.\n";
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

    std::cout << "Starting per-subject registration of FL to T1-Ce using Greedy.\n";
    std::system(greedyCommand.c_str());
    flImg = cbica::ReadImage< TImageType >(tempFile);
    std::cout << "Done.\n";
  }

  auto statsCalculator = itk::StatisticsImageFilter< TImageType >::New();
  statsCalculator->SetInput(t1cImg);
  statsCalculator->Update();
  if (statsCalculator->GetMean() != 0)
  {
    std::cout << "Starting Normalization of T1CE image.\n";
    ZScoreNormalizer< TImageType > normalizer;
    normalizer.SetInputImage(t1cImg);
    normalizer.SetInputMask(maskImage);
    normalizer.SetCutoffs(cutOffLower, cutOffUpper);
    normalizer.SetQuantiles(quantLower, quantUpper);
    normalizer.Update();
    t1cImg = normalizer.GetOutput();
    std::cout << "Finished Normalization of T1CE image.\n";
  }

  statsCalculator->SetInput(t1Img);
  statsCalculator->Update();
  if (statsCalculator->GetMean() != 0)
  {
    std::cout << "Starting Normalization of T1 image.\n";
    ZScoreNormalizer< TImageType > normalizer;
    normalizer.SetInputImage(t1Img);
    normalizer.SetInputMask(maskImage);
    normalizer.SetCutoffs(cutOffLower, cutOffUpper);
    normalizer.SetQuantiles(quantLower, quantUpper);
    normalizer.Update();
    t1Img = normalizer.GetOutput();
    std::cout << "Finished Normalization of T1 image.\n";
  }

  statsCalculator->SetInput(t2Img);
  statsCalculator->Update();
  if (statsCalculator->GetMean() != 0)
  {
    std::cout << "Starting Normalization of T2 image.\n";
    ZScoreNormalizer< TImageType > normalizer;
    normalizer.SetInputImage(t2Img);
    normalizer.SetInputMask(maskImage);
    normalizer.SetCutoffs(cutOffLower, cutOffUpper);
    normalizer.SetQuantiles(quantLower, quantUpper);
    normalizer.Update();
    t2Img = normalizer.GetOutput();
    std::cout << "Finished Normalization of T2 image.\n";
  }

  statsCalculator->SetInput(flImg);
  statsCalculator->Update();
  if (statsCalculator->GetMean() != 0)
  {
    std::cout << "Starting Normalization of T2-Flair image.\n";
    ZScoreNormalizer< TImageType > normalizer;
    normalizer.SetInputImage(flImg);
    normalizer.SetInputMask(maskImage);
    normalizer.SetCutoffs(cutOffLower, cutOffUpper);
    normalizer.SetQuantiles(quantLower, quantUpper);
    normalizer.Update();
    flImg = normalizer.GetOutput();
    std::cout << "Finished Normalization of T2-Flair image.\n";
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

  if (!cbica::isFile(modelDirName + "/modelConfig.txt"))
  {
    std::cerr << "'modelConfig.txt' was not found in the directory, please check.\n";
    return;
  }
  if (!cbica::isFile(modelDirName + "/model.ckpt"))
  {
    std::cerr << "'model.ckpt' was not found in the directory, please check.\n";
    return;
  }

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
    " -test " + cbica::normPath(getCaPTkDataDir() + "/deepmedic/configFiles/testApiConfig.txt") +
    " -o " + outputDirectory;

  if (std::system(fullCommand.c_str()) != 0)
  {
    std::cerr << "DeepMedic exited with code !=0.\n";
    exit(EXIT_FAILURE);
  }

  return;
}

int main(int argc, char **argv)
{
  cbica::CmdParser parser(argc, argv, "DeepMedic");

  parser.addRequiredParameter("t1c", "T1CE", cbica::Parameter::FILE, "", "The input T1CE or T1Gd image file.");
  parser.addRequiredParameter("t1", "T1", cbica::Parameter::FILE, "", "The input T1 image file.");
  parser.addRequiredParameter("fl", "FLAIR", cbica::Parameter::FILE, "", "The input T2-FLAIR image file.");
  parser.addRequiredParameter("t2", "FLAIR", cbica::Parameter::FILE, "", "The input T2 image file.");
  parser.addOptionalParameter("m", "mask", cbica::Parameter::FILE, "", "The Optional input mask file.", "This is needed for normalization only");
  parser.addOptionalParameter("md", "modelDir", cbica::Parameter::DIRECTORY, "", "The trained model to use", "Defaults to 'CaPTk_installDir/data/deepMedic/brainSegmentation'");
  parser.addRequiredParameter("o", "output", cbica::Parameter::DIRECTORY, "", "The output File.");

  parser.addOptionalParameter("ql", "quantLower", cbica::Parameter::FLOAT, "0-100", "The Lower Quantile range to remove", "This is needed for normalization only", "Default: 5");
  parser.addOptionalParameter("qu", "quantUpper", cbica::Parameter::FLOAT, "0-100", "The Upper Quantile range to remove", "This is needed for normalization only", "Default: 95");
  parser.addOptionalParameter("cl", "cutOffLower", cbica::Parameter::FLOAT, "0-10", "The Lower Cut-off (multiple of stdDev) to remove", "This is needed for normalization only", "Default: 3");
  parser.addOptionalParameter("cu", "cutOffUpper", cbica::Parameter::FLOAT, "0-10", "The Upper Cut-off (multiple of stdDev) to remove", "This is needed for normalization only", "Default: 3");
  parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");
  parser.exampleUsage("DeepMedic -t1 <input dir> -t1c <input dir> -t2 <input dir> -fl <input dir> -o <output dir>");

  // parameters to get from the command line
  cbica::Logging logger;

  if (parser.isPresent("L"))
  {
    parser.getParameterValue("L", loggerFileIn);
    logger.UseNewFile(loggerFileIn);
  }
  parser.getParameterValue("t1c", inputT1ce);
  parser.getParameterValue("t1", inputT1);
  parser.getParameterValue("t2", inputT2);
  parser.getParameterValue("fl", inputFlair);
  parser.getParameterValue("o", outputDirectory);

  if (parser.isPresent("m"))
  {
    parser.getParameterValue("m", inputMaskName);
    maskProvided = true;
  }

  if (parser.isPresent("md"))
  {
    parser.getParameterValue("md", modelDirName);
  }

  if (parser.isPresent("ql"))
  {
    parser.getParameterValue("ql", quantLower);
  }
  if (parser.isPresent("qu"))
  {
    parser.getParameterValue("qu", quantUpper);
  }

  if (parser.isPresent("cl"))
  {
    parser.getParameterValue("cl", cutOffLower);
  }
  if (parser.isPresent("cu"))
  {
    parser.getParameterValue("cu", cutOffUpper);
  }

  //std::cout << "Input File:" << inputFileName << std::endl;
  //if (!inputMaskName.empty())
  //{
  //  std::cout << "Input Mask:" << inputMaskName << std::endl;
  //}
  //std::cout << "Output File:" << outputFileName << std::endl;
  //std::cout << "Quant Lower:" << quantLower << std::endl;
  //std::cout << "Quant Upper:" << quantUpper << std::endl;
  //std::cout << "CutOff Lower:" << cutOffLower << std::endl;
  //std::cout << "CutOff Upper:" << cutOffUpper << std::endl;

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