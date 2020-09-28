#include "PerfusionDerivatives.h"
#include "cbicaUtilities.h"
#include "cbicaCmdParser.h"

bool IsWithinRange(float v)
{
	return 0.0 <= v && v <= 100.0;
}

int main(int argc, char **argv)
{
  cbica::CmdParser parser = cbica::CmdParser(argc, argv, "PerfusionDerivatives");
  parser.addRequiredParameter("i", "input", cbica::Parameter::STRING, "", "The input DSC-MRI image.");
  parser.addOptionalParameter("bs", "baseline start", cbica::Parameter::FLOAT, "", "The baseline start threshold percentage (range[0-100], Default = 0)");
  parser.addOptionalParameter("be", "baseline end", cbica::Parameter::FLOAT, "", "The baseline end threshold percentage (range[0-100], Default = 20)");
  parser.addOptionalParameter("rs", "recovery start", cbica::Parameter::FLOAT, "", "The recovery start threshold percentage (range[0-100], Default = 66)");
  parser.addOptionalParameter("re", "recovery end", cbica::Parameter::FLOAT, "", "The recovery end threshold percentage (range[0-100], Default = 88)");
  parser.addOptionalParameter("p", "PSR", cbica::Parameter::STRING, "", "Output the percentage signal recovery image (1=YES, 0=NO, 1 (Default))");
  parser.addOptionalParameter("pH", "peakHeight", cbica::Parameter::STRING, "", "Output the peak height image (1=YES, 0=NO, 1 (Default))");
  parser.addOptionalParameter("r", "ap-RCBV", cbica::Parameter::STRING, "", "Output the automatially-extracted proxy to reletive cerebral volume image (1=YES, 0=NO, 1 (Default))");
  parser.addRequiredParameter("o", "output", cbica::Parameter::STRING, "", "The output directory.");
  parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");
  parser.addExampleUsage("-i input_image.nii.gz -o <output dir>", 
    "Calculates the perfusion derivates of the input image.");
  parser.addApplicationDescription("Perfusion Derivatives calculation based on specific parameters");

  // parameters to get from the command line
  cbica::Logging logger;
  std::string loggerFile;
  bool loggerRequested = false;

  int tempPosition;
  int psrPresent = 1;
  int phPresent = 1;
  int rcbvPresent = 1;

  float bst,bet,rst,ret;

  bool bstValid = false;
  bool betValid = false;
  bool rstValid = false;
  bool retValid = false;

  //auto IsValid = [](float v){ return 0.0 <= v && v <= 100.0; };

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
  if (parser.compareParameter("p", tempPosition))
  {
	  psrPresent = atoi(argv[tempPosition + 1]);
  }
  if (parser.compareParameter("pH", tempPosition))
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

  if (parser.isPresent("bs"))
  {
	  parser.getParameterValue("bs", bst);
	  bstValid = IsWithinRange(bst);
	  if (!bstValid)
	  {
		  std::cout << "Baseline start threshold value must be between 0 and 100";
		  return EXIT_FAILURE;
	  }
  }
  if (parser.isPresent("be"))
  {
	  parser.getParameterValue("be", bet);
	  betValid = IsWithinRange(bet);
	  if (!betValid)
	  {
		  std::cout << "Baseline end threshold value must be between 0 and 100";
		  return EXIT_FAILURE;
	  }
  }
  if (parser.isPresent("rs"))
  {
	  parser.getParameterValue("rs", rst);
	  rstValid = IsWithinRange(rst);
	  if (!rstValid)
	  {
		  std::cout << "Recovery start threshold value must be between 0 and 100";
		  return EXIT_FAILURE;
	  }
  }
  if (parser.isPresent("re"))
  {
	  parser.getParameterValue("re", ret);
	  retValid = IsWithinRange(ret);
	  if (!retValid)
	  {
		  std::cout << "Recovery end threshold value must be between 0 and 100";
		  return EXIT_FAILURE;
	  }
  }

  std::cout << "Input File:" << inputFileName << std::endl;
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

  if (psrPresent == 0 && phPresent ==0 && rcbvPresent==0)
  {
    std::cout << "Please select atleast one of the given three measures (PSR, PH, RCBV)." << std::endl;
    return EXIT_FAILURE;
  }
  PerfusionDerivatives objPerfusion;
  if (bstValid)
	  objPerfusion.SetBaselineStartPercentage(bst);
  if (betValid)
	  objPerfusion.SetBaselineEndPercentage(bet);
  if (rstValid)
	  objPerfusion.SetRecoveryStartPercentage(rst);
  if (retValid)
	  objPerfusion.SetRecoveryEndPercentage(ret);
  std::vector<typename ImageTypeFloat3D::Pointer> perfusionDerivatives = objPerfusion.Run<ImageTypeFloat3D, ImageTypeFloat4D>(inputFileName, rcbvPresent, psrPresent, phPresent, outputDirectoryName);

  //write perfusion derivatives
  if (perfusionDerivatives.size() > 0)
  {
    if (psrPresent==1)
      cbica::WriteImage<ImageTypeFloat3D>(perfusionDerivatives[0], outputDirectoryName + "/PSR.nii.gz");
    if (phPresent==1)
      cbica::WriteImage<ImageTypeFloat3D>(perfusionDerivatives[1], outputDirectoryName + "/PH.nii.gz");
    if (rcbvPresent==1)
      cbica::WriteImage<ImageTypeFloat3D>(perfusionDerivatives[2], outputDirectoryName + "/ap-RCBV.nii.gz");
    std::cout << "Perfusion derivatives written to the specified location." << std::endl;
    std::cout << "Finished successfully.\n";
  }
  else
  {
    std::cerr << "No derivates were calculated. Please make sure that the input image does not contain negative values.\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}