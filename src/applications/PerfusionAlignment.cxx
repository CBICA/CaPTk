#include "PerfusionAlignment.h"
#include "cbicaUtilities.h"
#include "cbicaCmdParser.h"

int main(int argc, char **argv)
{
  cbica::CmdParser parser = cbica::CmdParser(argc, argv, "PerfusionAlignment");
  parser.addRequiredParameter("i", "input", cbica::Parameter::STRING, "", "The input DSC-MRI image.");
  parser.addRequiredParameter("b", "time-points before drop", cbica::Parameter::STRING, "", "The number of time-points before the drop.");
  parser.addRequiredParameter("a", "time-points after drop", cbica::Parameter::STRING, "", "The number of time-points after the drop.");
  parser.addRequiredParameter("t", "time-domain resolution", cbica::Parameter::FLOAT, "", "The time-interval (spacing) between two consecutive volumes in time-domain (in seconds).");

  parser.addRequiredParameter("o", "output", cbica::Parameter::STRING, "", "The output directory.");
  parser.addOptionalParameter("s", "drop-scaling", cbica::Parameter::BOOLEAN, "", "Whether to scale the value of the drop for the curve? 1=yes, 0=no, default=0.");
  parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");
  //parser.exampleUsage("PerfusionAlignment -i AAAC_PreOp_perf_pp.nii.gz -d AAAC_PreOp_perf_pp.dcm -o <output dir>");
  parser.addExampleUsage("-i AAAC_PreOp_perf_pp.nii.gz -c AAAC_PreOp_t1ce_pp.nii.gz -b 15 -a 17 -t 2 -o <output dir>", "Aligns the perfusion signal of the input image based on the time points");
  parser.addApplicationDescription("Perfusion Alignment of the input based based on specified time points");

  // parameters to get from the command line
  cbica::Logging logger;
  std::string loggerFile;
  bool loggerRequested = false;
  int tempPosition;
  int pointsbeforedrop, pointsafterdrop;
  bool dropscaling = 0;
  double timeresolution;
  std::string inputFileName, inputDicomName, outputDirectoryName, inputt1ceName;

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

  if (parser.compareParameter("b", tempPosition))
    pointsbeforedrop = atoi(argv[tempPosition + 1]);

  if (parser.compareParameter("t", tempPosition))
    timeresolution = atof(argv[tempPosition + 1]);

  if (parser.compareParameter("a", tempPosition))
    pointsafterdrop = atoi(argv[tempPosition + 1]);

  if (parser.isPresent("s"))
  {
    parser.getParameterValue("s", dropscaling);
  }

  if (!cbica::isFile(inputFileName))
  {
    std::cout << "The input file does not exist:" << inputFileName << std::endl;
    return EXIT_FAILURE;
  }
  cbica::createDirectory(outputDirectoryName);

  PerfusionAlignment objPerfusion;
  std::vector<double> OriginalCurve, InterpolatedCurve, RevisedCurve, TruncatedCurve;
  std::vector<typename ImageTypeFloat3D::Pointer> PerfusionAlignment = objPerfusion.Run<ImageTypeFloat3D, ImageTypeFloat4D>(inputFileName, pointsbeforedrop, pointsafterdrop, OriginalCurve, InterpolatedCurve, RevisedCurve, TruncatedCurve, timeresolution, dropscaling);

  if (!PerfusionAlignment.empty())
  {
    auto joinedImage = cbica::GetJoinedImage< ImageTypeFloat3D, ImageTypeFloat4D >(PerfusionAlignment);
    cbica::WriteImage< ImageTypeFloat4D >(joinedImage, outputDirectoryName + "/perfusionAlignedImage.nii.gz");

    WriteCSVFiles(OriginalCurve, outputDirectoryName + "/original_curve.csv");
    WriteCSVFiles(InterpolatedCurve, outputDirectoryName + "/interpolated_curve.csv");
    WriteCSVFiles(RevisedCurve, outputDirectoryName + "/revised_curve.csv");
    WriteCSVFiles(TruncatedCurve, outputDirectoryName + "/truncated_curve.csv");

    std::cout << "Finished successfully.\n";
    std::cout << "\nPress any key to continue............\n";
  }
  else
  {
    std::cerr << "Something went wrong and CaPTk could not align the perfusion signal correctly. Please see prior messages for details.\n";
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
