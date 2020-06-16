#include "PerfusionAlignment.h"
#include "cbicaUtilities.h"
#include "cbicaCmdParser.h"

int main(int argc, char **argv)
{
  cbica::CmdParser parser = cbica::CmdParser(argc, argv, "PerfusionAlignment");
  parser.addRequiredParameter("i", "input", cbica::Parameter::STRING, "", "The input DSC-MRI image.");
  parser.addRequiredParameter("c", "t1ce file", cbica::Parameter::STRING, "", "The input T1 post-weighted image.");
  parser.addRequiredParameter("b", "time-points before drop", cbica::Parameter::STRING, "", "The number of time-points before the drop.");
  parser.addRequiredParameter("a", "time-points after drop", cbica::Parameter::STRING, "", "The number of time-points after the drop.");
  parser.addRequiredParameter("e", "echo time", cbica::Parameter::FLOAT, "", "Echo time.");

  parser.addRequiredParameter("o", "output", cbica::Parameter::STRING, "", "The output directory.");
  parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");
  //parser.exampleUsage("PerfusionAlignment -i AAAC_PreOp_perf_pp.nii.gz -d AAAC_PreOp_perf_pp.dcm -o <output dir>");
  parser.addExampleUsage("-i AAAC_PreOp_perf_pp.nii.gz -d AAAC_PreOp_perf_pp.dcm -o <output dir>", "Aligns the perfusion signal of the input image based on the time points");
  parser.addApplicationDescription("Perfusion Alignment of the input based based on specified time points");
  

  // parameters to get from the command line
  cbica::Logging logger;
  std::string loggerFile;
  bool loggerRequested = false;
  int tempPosition;
  int pointsbeforedrop, pointsafterdrop;
  double echotime;
  std::string inputFileName, inputDicomName, outputDirectoryName,inputt1ceName;

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
  if (parser.compareParameter("o", tempPosition))
  {
    outputDirectoryName = argv[tempPosition + 1];
  }
  if (parser.compareParameter("c", tempPosition))
  {
    inputt1ceName = argv[tempPosition + 1];
  }

  if (parser.compareParameter("b", tempPosition))
    pointsbeforedrop = atoi(argv[tempPosition + 1]);

  if (parser.compareParameter("e", tempPosition))
    echotime = atof(argv[tempPosition + 1]);

  if (parser.compareParameter("a", tempPosition))
    pointsafterdrop = atoi(argv[tempPosition + 1]);

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
  
  PerfusionAlignment objPerfusion;
  std::vector<double> OriginalCurve, RevisedCurve;
  std::vector<typename ImageTypeFloat3D::Pointer> PerfusionAlignment = objPerfusion.Run<ImageTypeFloat3D, ImageTypeFloat4D>(inputFileName,inputt1ceName, pointsbeforedrop,pointsafterdrop,OriginalCurve,RevisedCurve,echotime);
  //std::vector<typename ImageTypeFloat3D::Pointer> PerfusionAlignment = objPerfusion.Run<ImageTypeFloat3D, ImageTypeFloat4D>("//cbica-cifs/hasun/comp_space/180815_Henry_Ford/Protocols/5_SSFinal/2/2/2_perf_LPS_r_SSFinal.nii.gz", "W:/perf/MSh_PERF_AX-1001_echo1_I000001.dcm", "//cbica-cifs/hasun/comp_space/180815_Henry_Ford/Protocols/5_SSFinal/2/2/2_t1ce_LPS_r_SSFinal.nii.gz", 17, 40);
  for (int index = 0; index < PerfusionAlignment.size(); index++)
  {
    std::cout << "Writing time-point: " << index+1 << "/" << PerfusionAlignment.size() << std::endl;
    cbica::WriteImage<ImageTypeFloat3D>(PerfusionAlignment[index], outputDirectoryName + std::to_string(index+1+pointsbeforedrop) + ".nii.gz");
  }
  
  std::ofstream myfile;
  myfile.open(outputDirectoryName + "/original_curve.csv");
  for (unsigned int index1 = 0; index1 < OriginalCurve.size(); index1++)
    myfile << std::to_string(OriginalCurve[index1]) << "\n";
  myfile.close();

  myfile.open(outputDirectoryName + "/revised_curve.csv");
  for (unsigned int index1 = 0; index1 < RevisedCurve.size(); index1++)
    myfile << std::to_string(RevisedCurve[index1]) << "\n";
  myfile.close();


  std::cout << "Finished successfully.\n";
  std::cout << "\nPress any key to continue............\n";

  return EXIT_SUCCESS;
}

//-i //cbica-cifs/hasun/comp_space/180815_Henry_Ford/Protocols/5_SSFinal/2/2/2_perf_LPS_r_SSFinal.nii.gz -d W:/perf/MSh_PERF_AX-1001_echo1_I000001.dcm -c //cbica-cifs/hasun/comp_space/180815_Henry_Ford/Protocols/5_SSFinal/2/2/2_t1ce_LPS_r_SSFinal.nii.gz -b 17 -a 40 -o W:/perf