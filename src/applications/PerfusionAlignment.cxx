#include "PerfusionAlignment.h"
#include "cbicaUtilities.h"
#include "cbicaCmdParser.h"

int main(int argc, char **argv)
{
  std::cout << "This functionality has been removed from this CaPTk release, \
and we are actively testing an optimized robust implementation that would enable \
generalization in multi-institutional data. We expect this to be released in our \
next patch release, expected in Q4 2020.\n";
  return EXIT_FAILURE;
  int tempPosition;
  int pointsbeforedrop, pointsafterdrop;
  bool dropscaling = 0;
  double timeresolution;
  size_t time_beforeDrop, time_afterDrop;
  float baseline = 300, stdDev = 10, time_inputPerfTime, time_outputPerfTime = 1.0, scale_maxIntensityBeforeDrop = 300, scale_intensityDropInMeanCurve = 100;
  cbica::CmdParser parser = cbica::CmdParser(argc, argv, "PerfusionAlignment");
  parser.addRequiredParameter("i", "input", cbica::Parameter::FILE, "File with read access", "The input DSC-MRI image.");
  parser.addOptionalParameter("m", "mask", cbica::Parameter::STRING, "File with read access", "The mask or the type of default masking", "1: otsu+stdDev-drop; 2: stdDev-volume+otsu");
  parser.addRequiredParameter("t1", "tInputPerfRep", cbica::Parameter::FLOAT, "0-100", "The input perfusion repetition time");
  parser.addOptionalParameter("t2", "tOutputPerfRep", cbica::Parameter::FLOAT, "0-100", "The output perfusion repetition time", "Defaults to: " + std::to_string(time_outputPerfTime));
  parser.addRequiredParameter("t3", "tBeforeDrop", cbica::Parameter::INTEGER, "0-100", "Number of timepoints BEFORE drop (starting index '1')");
  parser.addRequiredParameter("t4", "tAfterDrop", cbica::Parameter::INTEGER, "0-100", "Number of timepoints AFTER drop (starting index '1')");

  parser.addOptionalParameter("s1", "sMaxBeforeDrop", cbica::Parameter::FLOAT, "0-100", "The output perfusion repetition time", "Defaults to: " + std::to_string(scale_maxIntensityBeforeDrop));
  parser.addOptionalParameter("s2", "sDropInMeanCurve", cbica::Parameter::FLOAT, "0-100", "The intensity drop in mean curve inside ROI", "Defaults to: " + std::to_string(scale_intensityDropInMeanCurve));

  parser.addRequiredParameter("o", "output", cbica::Parameter::STRING, "Directory with write access", "The output directory.");


  //parser.addRequiredParameter("b", "timePointsBeforeDrop", cbica::Parameter::STRING, "0-100", "The number of time-points before the drop.");
  //parser.addRequiredParameter("a", "timePointsAfterDrop", cbica::Parameter::STRING, "0-100", "The number of time-points after the drop.");
  //parser.addRequiredParameter("t", "timeDomainResolution", cbica::Parameter::FLOAT, "0.01-100", "The time-interval (spacing) between two consecutive volumes in time-domain (in seconds).");

  parser.addOptionalParameter("d", "dropScaling", cbica::Parameter::BOOLEAN, "0-1", "Whether to scale the value of the drop for the curve? 1=yes, 0=no; defaults to '0'");
  parser.addOptionalParameter("s", "stdDev", cbica::Parameter::FLOAT, "0-100", "The standard deviation threshold of time series signal above which the", "location is considered to part of brain", "Defaults to " + std::to_string(stdDev));
  parser.addOptionalParameter("bl", "baseLine", cbica::Parameter::FLOAT, "0-1000", "The value of the baseline to which the output gets scaled", "Only used if '-d' is '1'", "Defaults to " + std::to_string(baseline));
  parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");
  //parser.exampleUsage("PerfusionAlignment -i AAAC_PreOp_perf_pp.nii.gz -d AAAC_PreOp_perf_pp.dcm -o <output dir>");
  parser.addExampleUsage("-i AAAC_PreOp_perf_pp.nii.gz -c AAAC_PreOp_t1ce_pp.nii.gz -b 15 -a 17 -t 2 -o <output dir>", "Aligns the perfusion signal of the input image based on the time points");
  parser.addApplicationDescription("Perfusion Alignment of the input based based on specified time points");

  // parameters to get from the command line
  cbica::Logging logger;
  std::string loggerFile;
  bool loggerRequested = false;
  std::string inputFileName, inputDicomName, outputDirectoryName, inputt1ceName;

  parser.getParameterValue("i", inputFileName);
  parser.getParameterValue("t1", time_inputPerfTime);
  parser.getParameterValue("t3", time_beforeDrop);
  parser.getParameterValue("t4", time_afterDrop);

  if (parser.isPresent("t2"))
  {
    parser.getParameterValue("t2", time_outputPerfTime);
    if (time_outputPerfTime < 0)
    {
      std::cerr << "The output time resolution cannot be less than 0, using default, instead.\n";
      time_outputPerfTime = 1;
    }
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

  if (parser.isPresent("d"))
  {
    parser.getParameterValue("d", dropscaling);
  }
  if (parser.isPresent("s"))
  {
    parser.getParameterValue("s", stdDev);
  }
  if (parser.isPresent("bl"))
  {
    parser.getParameterValue("bl", baseline);
  }

  if (!cbica::isFile(inputFileName))
  {
    std::cout << "The input file does not exist:" << inputFileName << std::endl;
    return EXIT_FAILURE;
  }
  cbica::createDirectory(outputDirectoryName);

  PerfusionAlignment objPerfusion;
  std::vector<double> OriginalCurve, InterpolatedCurve, RevisedCurve, TruncatedCurve;
  //std::vector<typename ImageTypeFloat3D::Pointer>  = 
  auto output =  objPerfusion.Run<ImageTypeFloat3D, ImageTypeFloat4D>(inputFileName, pointsbeforedrop, pointsafterdrop, OriginalCurve, InterpolatedCurve, RevisedCurve, TruncatedCurve, timeresolution, dropscaling, stdDev, baseline);

  auto PerfusionAlignment = output.first;
  auto calculatedMask = output.second;

  if (!PerfusionAlignment.empty())
  {
    auto joinedImage = cbica::GetJoinedImage< ImageTypeFloat3D, ImageTypeFloat4D >(PerfusionAlignment);
    cbica::WriteImage< ImageTypeFloat4D >(joinedImage, outputDirectoryName + "/perfusionAlignedImage.nii.gz");

    cbica::WriteImage< ImageTypeFloat3D >(calculatedMask, outputDirectoryName + "/calculatedMask.nii.gz");

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
