#include "PerfusionDerivatives.h"
#include "cbicaUtilities.h"
#include "cbicaCmdParser.h"

int main(int argc, char **argv)
{
  cbica::CmdParser parser = cbica::CmdParser(argc, argv, "PerfusionDerivatives");
  parser.addRequiredParameter("i", "input", cbica::Parameter::STRING, "", "The input DSC-MRI image.");
  parser.addRequiredParameter("e", "echo time", cbica::Parameter::STRING, "", "The echo time in seconds.");
  parser.addRequiredParameter("p", "PSR", cbica::Parameter::STRING, "", "The Percent Signal Recovery image (1=YES, 0=NO, 1 (Default))");
  parser.addRequiredParameter("pH", "peekHeight", cbica::Parameter::STRING, "", "The Peak Height image (1=YES, 0=NO, 1 (Default))");
  parser.addRequiredParameter("r", "ap-RCBV", cbica::Parameter::STRING, "", "Automatially-extracted proxy to reletive cerebral volume image (1=YES, 0=NO, 1 (Default))");
  parser.addRequiredParameter("o", "output", cbica::Parameter::STRING, "", "The output directory.");
  parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");
  parser.exampleUsage("PerfusionDerivatives -i AAAC_PreOp_perf_pp.nii.gz -e 1 -o <output dir> -p 1 -r 1");

  // parameters to get from the command line
  cbica::Logging logger;
  std::string loggerFile;
  bool loggerRequested = false;

  int tempPosition;
  int psrPresent = 1;
  int phPresent = 1;
  int rcbvPresent = 1;

  double inputEchoName = 0;
  std::string inputFileName, outputDirectoryName;

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

  if (parser.compareParameter("e", tempPosition))
  {
	  inputEchoName = atof(argv[tempPosition + 1]);
  }
  if (parser.compareParameter("p", tempPosition))
  {
	  psrPresent = atoi(argv[tempPosition + 1]);
  }
  if (parser.compareParameter("h", tempPosition))
  {
	  phPresent = atoi(argv[tempPosition + 1]);
  }
  if (parser.compareParameter("r", tempPosition))
  {
	  rcbvPresent = atoi(argv[tempPosition + 1]);
  }

  if (parser.compareParameter("o", tempPosition))
  {
    outputDirectoryName = argv[tempPosition + 1];
  }

  // std::cout << "Input File:" << inputFileName << std::endl;
  // std::cout << "Output Directory:" << outputDirectoryName << std::endl;
  // cbica::Logging(loggerFile, "Input directory name: " + inputFileName + "\n");
	// cbica::Logging(loggerFile, "Output directory name: " + outputDirectoryName + "\n");

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
  if (psrPresent == 0 && phPresent ==0 && rcbvPresent==0)
  {
    std::cout << "Please select atleast one of the given three measures (PSR, PH, RCBV)." << std::endl;
    return EXIT_FAILURE;
  }
  PerfusionDerivatives objPerfusion;
  std::vector<typename ImageTypeFloat3D::Pointer> measures = objPerfusion.Run<ImageTypeFloat3D, ImageTypeFloat4D>(inputFileName, rcbvPresent, psrPresent, phPresent, inputEchoName);
  // std::cout << "Writing measures to the specified output directory.\n";
  for (int i = 0; i < measures.size(); i++)
  {
    typedef itk::ImageFileWriter< ImageTypeFloat3D > WriterType;
    WriterType::Pointer writer1 = WriterType::New();
    writer1->SetFileName(outputDirectoryName + "/Atlas_" + std::to_string(i + 1) + ".nii.gz");
    writer1->SetInput(measures[i]);
    writer1->Update();
  }
  // std::cout << "Finished successfully.\n";
  // std::cout << "\nPress any key to continue............\n";

  return EXIT_SUCCESS;
}