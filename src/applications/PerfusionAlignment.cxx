#include "PerfusionAlignment.h"
#include "cbicaUtilities.h"
#include "cbicaCmdParser.h"

int main(int argc, char **argv)
{
  cbica::CmdParser parser = cbica::CmdParser(argc, argv, "PerfusionAlignment");
  parser.addRequiredParameter("i", "input", cbica::Parameter::STRING, "", "The input DSC-MRI image.");
  parser.addRequiredParameter("d", "dicom file", cbica::Parameter::STRING, "", "The input dicom image.");
  parser.addRequiredParameter("o", "output", cbica::Parameter::STRING, "", "The output directory.");
  parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");
  parser.exampleUsage("PerfusionAlignment -i AAAC_PreOp_perf_pp.nii.gz -d AAAC_PreOp_perf_pp.dcm -o <output dir>");
  PerfusionAlignment objPerfusion;

  // parameters to get from the command line
  cbica::Logging logger;
  std::string loggerFile;
  bool loggerRequested = false;

  int tempPosition;
  double inputEchoName = 0;
  std::string inputFileName, inputDicomName, outputDirectoryName;

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
  if (parser.compareParameter("d", tempPosition))
  {
    inputDicomName = argv[tempPosition + 1];
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
  
  std::vector< /*typename*/ ImageTypeFloat3D::Pointer > PerfusionAlignment = objPerfusion.Run< ImageTypeFloat3D, ImageTypeFloat4D >(inputFileName,inputDicomName);
  std::cout << "Writing measures to the specified output directory.\n";

  if (PerfusionAlignment.size() == 0)
  {
  }
  else
  {
    //if (psrPresent == true)
    //  cbica::WriteImage< ImageTypeFloat3D >(PerfusionAlignment[0], outputDirectoryName + "/PSR.nii.gz");
    //if (phPresent == true)
    //  cbica::WriteImage< ImageTypeFloat3D >(PerfusionAlignment[1], outputDirectoryName + "/PH.nii.gz");
    //if (rcbvPresent == true)
    //  cbica::WriteImage< ImageTypeFloat3D >(PerfusionAlignment[2], outputDirectoryName + "/ap-RCBV.nii.gz");

    std::cout << "Perfusion derivatives have been saved at the specified locations.\n";
  }
  std::cout << "Finished successfully.\n";
  std::cout << "\nPress any key to continue............\n";

  return EXIT_SUCCESS;
}