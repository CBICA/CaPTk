#include "PerfusionAlignment.h"
#include "cbicaUtilities.h"
#include "cbicaCmdParser.h"

int main(int argc, char **argv)
{
  size_t time_beforeDrop, time_afterDrop;
  float time_inputPerfTime, time_outputPerfTime = 1.0, scale_maxIntensityBeforeDrop = 300, scale_intensityDropInMeanCurve = 100;
  cbica::CmdParser parser = cbica::CmdParser(argc, argv, "PerfusionAlignment");
  parser.addRequiredParameter("i", "input", cbica::Parameter::FILE, "File with read access", "The input DSC-MRI image.");
  parser.addOptionalParameter("m", "mask", cbica::Parameter::STRING, "File with read access", "The brain mask corresponding to the input", "If not defined, it performs Otsu thresholding on input");
  parser.addRequiredParameter("t1", "tInputPerfRep", cbica::Parameter::FLOAT, "0-100", "The input perfusion repetition time");
  parser.addOptionalParameter("t2", "tOutputPerfRep", cbica::Parameter::FLOAT, "0-100", "The output perfusion repetition time", "Defaults to: " + std::to_string(time_outputPerfTime));
  parser.addRequiredParameter("t3", "tBeforeDrop", cbica::Parameter::INTEGER, "0-100", "Number of timepoints BEFORE drop (starting index '1')");
  parser.addRequiredParameter("t4", "tAfterDrop", cbica::Parameter::INTEGER, "0-100", "Number of timepoints AFTER drop (starting index '1')");

  parser.addOptionalParameter("s1", "sMaxBeforeDrop", cbica::Parameter::FLOAT, "0-100", "The maximum intensity before drop in mean curve inside ROI", "Defaults to: " + std::to_string(scale_maxIntensityBeforeDrop));
  parser.addOptionalParameter("s2", "sDropInMeanCurve", cbica::Parameter::FLOAT, "0-100", "The intensity drop in mean curve inside ROI", "Defaults to: " + std::to_string(scale_intensityDropInMeanCurve));

  parser.addRequiredParameter("o", "output", cbica::Parameter::STRING, "Directory with write access", "The output directory.");

  parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");
  parser.addExampleUsage("-i AAAC_PreOp_perf_pp.nii.gz -m brainMask.nii.gz -t1 2 -t3 15 -t4 30 -o C:/testPerf", "Aligns the perfusion signal of the input image based on the time points");
  parser.addApplicationDescription("Perfusion Alignment of the input based based on specified time points");

  // parameters to get from the command line
  cbica::Logging logger;
  std::string loggerFile;
  bool loggerRequested = false;
  std::string inputFileName, inputMaskFileName, inputDicomName, outputDirectoryName;

  parser.getParameterValue("i", inputFileName);
  parser.getParameterValue("t1", time_inputPerfTime);
  parser.getParameterValue("t3", time_beforeDrop);
  parser.getParameterValue("t4", time_afterDrop);
  parser.getParameterValue("o", outputDirectoryName);

  if (parser.isPresent("m"))
  {
    parser.getParameterValue("m", inputMaskFileName);
  }

  if (parser.isPresent("t2"))
  {
    parser.getParameterValue("t2", time_outputPerfTime);
    if (time_outputPerfTime < 0)
    {
      std::cerr << "The output time resolution cannot be less than 0, using default, instead.\n";
      time_outputPerfTime = 1;
    }
  }
  if (parser.isPresent("s1"))
  {
    parser.getParameterValue("s1", scale_maxIntensityBeforeDrop);
  }
  if (parser.isPresent("s2"))
  {
    parser.getParameterValue("s2", scale_intensityDropInMeanCurve);
  }
  if (parser.isPresent("L"))
  {
    parser.getParameterValue("L", loggerFile);
    loggerRequested = true;
    logger.UseNewFile(loggerFile);
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
  auto output = objPerfusion.Run<ImageTypeFloat3D, ImageTypeFloat4D>(inputFileName, inputMaskFileName, time_beforeDrop, time_afterDrop, OriginalCurve, InterpolatedCurve, RevisedCurve, TruncatedCurve,
    time_inputPerfTime, time_outputPerfTime, scale_maxIntensityBeforeDrop, scale_intensityDropInMeanCurve); // , dropscaling, stdDev, baseline);

  auto PerfusionAlignment = output.first;
  auto calculatedMask = output.second;

  if (!PerfusionAlignment.empty())
  {
    auto joinedImage = cbica::GetJoinedImage< ImageTypeFloat3D, ImageTypeFloat4D >(PerfusionAlignment, time_outputPerfTime);
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
