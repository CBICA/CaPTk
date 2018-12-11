#include "PerfusionPCA.h"
#include "cbicaUtilities.h"
#include "cbicaCmdParser.h"
#include "cbicaITKSafeImageIO.h"

int main(int argc, char **argv)
{
  cbica::CmdParser parser = cbica::CmdParser(argc, argv, "PerfusionPCA");
  parser.addRequiredParameter("i", "input", cbica::Parameter::STRING, "", "The input DSC-MRI image.");
  parser.addRequiredParameter("m", "mask", cbica::Parameter::STRING, "", "The input mask.");
  parser.addRequiredParameter("n", "number of PCAs", cbica::Parameter::STRING, "", "The number of principal components.");
  parser.addRequiredParameter("o", "output", cbica::Parameter::STRING, "", "The output directory.");
  parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");

  parser.writeCWLFile(".", parser.getExeName(), false);
  parser.exampleUsage("");

  // parameters to get from the command line
  cbica::Logging logger;
  std::string loggerFile;
  bool loggerRequested = false;

  int tempPosition;

  double inputPCs = 0;
  std::string inputFileName, inputMaskName, outputDirectoryName;

  if ((argc < 1) || (parser.compareParameter("u", tempPosition)))
  {
    parser.echoUsage();
    return EXIT_SUCCESS;
  }
  if (parser.compareParameter("L", tempPosition))
  {
    loggerFile = argv[tempPosition + 1];
    loggerRequested = true;
    logger.UseNewFile(loggerFile);
  }
  if (parser.compareParameter("i", tempPosition))
  {
    inputFileName = argv[tempPosition + 1];
  }
  if (parser.compareParameter("m", tempPosition))
  {
    inputMaskName = argv[tempPosition + 1];
  }
  if (parser.compareParameter("n", tempPosition))
  {
	  inputPCs = atof(argv[tempPosition + 1]);
  }
  if (parser.compareParameter("o", tempPosition))
  {
    outputDirectoryName = argv[tempPosition + 1];
  }

  std::cout << "Input File:" << inputFileName << std::endl;
  std::cout << "Input Mask:" << inputMaskName << std::endl;
  std::cout << "Output Directory:" << outputDirectoryName << std::endl;

  if (!cbica::isFile(inputFileName))
  {
    std::cout << "The input file does not exist:" << inputFileName << std::endl;
    return EXIT_FAILURE;
  }
  if (!cbica::directoryExists(outputDirectoryName))
  {
    if (!cbica::createDirectory(outputDirectoryName))
      std::cout << "The output directory can not be created:" << outputDirectoryName << std::endl;
    return EXIT_FAILURE;
  }
  ImageTypeFloat4D::Pointer perfusionImage = cbica::ReadImage<ImageTypeFloat4D>(inputFileName);
  ImageTypeFloat3D::Pointer maskImage = cbica::ReadImage<ImageTypeFloat3D>(inputMaskName);

  PerfusionPCA object_pca;
  std::vector<ImageTypeFloat3D::Pointer> individual_pcs = object_pca.Run<ImageTypeFloat4D, ImageTypeFloat3D>(maskImage, perfusionImage);
    for (int index = 0; index < inputPCs; index++)
      cbica::WriteImage< ImageTypeFloat3D >(individual_pcs[index], outputDirectoryName + "/pca_" + std::to_string(index) + ".nii.gz");

  std::cout<<"principal components have been saved at the specified locations."; 
  std::cout << "Finished successfully.\n";

  return EXIT_SUCCESS;
}