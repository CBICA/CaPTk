#include "DeepMedicNormalizer.h"
#include "cbicaUtilities.h"
#include "cbicaCmdParser.h"
#include "cbicaLogging.h"
#include "cbicaITKUtilities.h"
#include "cbicaITKImageInfo.h"
#include "cbicaITKSafeImageIO.h"

std::string inputFileName, inputMaskName, inputBValName, inputBVecName, outputFileName, loggerFile;
float quantLower = 5, quantUpper = 95, cutOffLower = 3, cutOffUpper = 3;
bool maskProvided = false;

template< class TImageType >
void algorithmRunner()
{
  auto inputImage = cbica::ReadImage< TImageType >(inputFileName);
  auto maskImage = cbica::CreateImage< TImageType >(inputImage);
  if (maskProvided)
  {
    maskImage = cbica::ReadImage< TImageType >(inputMaskName);
  }

  std::cout << "Starting Normalization.\n";
  DeepMedicNormalizer< TImageType > normalizer;
  normalizer.SetInputImage(inputImage);
  normalizer.SetInputMask(maskImage);
  normalizer.SetCutoffs(cutOffLower, cutOffUpper);
  normalizer.SetQuantiles(quantLower, quantUpper);
  normalizer.Update();
  std::cout << "Finished Normalization.\n";

  cbica::WriteImage< TImageType >(normalizer.GetOutput(), outputFileName);

  return;
}

int main(int argc, char **argv)
{
  cbica::CmdParser parser(argc, argv);

  parser.addRequiredParameter("i", "input", cbica::Parameter::FILE, "", "The input image file.");
  parser.addOptionalParameter("m", "mask", cbica::Parameter::FILE, "", "The Optional input mask file.");
  parser.addRequiredParameter("o", "output", cbica::Parameter::FILE, "", "The output File.");

  parser.addOptionalParameter("ql", "quantLower", cbica::Parameter::FLOAT, "0-100", "The Lower Quantile range to remove", "Default: 5");
  parser.addOptionalParameter("qu", "quantUpper", cbica::Parameter::FLOAT, "0-100", "The Upper Quantile range to remove", "Default: 95");
  parser.addOptionalParameter("cl", "cutOffLower", cbica::Parameter::FLOAT, "0-10", "The Lower Cut-off (multiple of stdDev) to remove", "Default: 3");
  parser.addOptionalParameter("cu", "cutOffUpper", cbica::Parameter::FLOAT, "0-10", "The Upper Cut-off (multiple of stdDev) to remove", "Default: 3");
  parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");

  // parameters to get from the command line
  cbica::Logging logger;

  if (parser.isPresent("L"))
  {
    parser.getParameterValue("L", loggerFile);
    logger.UseNewFile(loggerFile);
  }
  parser.getParameterValue("i", inputFileName);
  parser.getParameterValue("o", outputFileName);

  if (parser.isPresent("m"))
  {
    parser.getParameterValue("m", inputMaskName);
    maskProvided = true;
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

  auto imageInfo = cbica::ImageInfo(inputFileName);
  if (maskProvided)
  {
    if(!cbica::ImageSanityCheck(inputFileName, inputMaskName))
    {
      return EXIT_FAILURE;
    }
  }

  switch (imageInfo.GetImageDimensions())
  {
  case 2:
  {
    using ImageType = itk::Image< float, 2 >;
    algorithmRunner< ImageType >();

    break;
  }
  case 3:
  {
    using ImageType = itk::Image< float, 3 >;
    algorithmRunner< ImageType >();

    break;
  }
  default:
    std::cerr << "Only 2D or 3D images are supported right now.\n";
    return EXIT_FAILURE;
  }


  std::cout << "Finished successfully.\n";
  return EXIT_SUCCESS;
}