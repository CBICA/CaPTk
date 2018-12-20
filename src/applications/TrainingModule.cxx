#include "TrainingModule.h"
#include <ctime>
#include "cbicaCmdParser.h"
//#include "CAPTk.h"

void showVersionInfo()
{
  //cout << "Survival Predictor: version: " << VERSION_STRING_W << endl;
}
void showhelpinfo(const std::string exeName)
{
  std::cout << "Usage:   " << exeName << " [-option] [argument]" << std::endl;
  std::cout << "Option:  " << "-h show help information" << std::endl;
  std::cout << "         " << "-v show version information" << std::endl;
  std::cout << "         " << "-f feature file path" << std::endl;
  std::cout << "         " << "-l label file path" << std::endl;
  std::cout << "         " << "-o output directory path" << std::endl;
  std::cout << "Example: " << exeName << " -f ../features.csv -l ../targets.csv -o ../outputfolder" << std::endl;
}

int main(int argc, char *argv[])
{
  cbica::CmdParser parser = cbica::CmdParser(argc, argv, "TrainingModule");
  parser.addRequiredParameter("f", "features", cbica::Parameter::STRING, "", "The input file having features (*.csv).");
  parser.addRequiredParameter("l", "label", cbica::Parameter::STRING, "", "The input file having target labels (*.csv).");
  parser.addRequiredParameter("c", "classifier", cbica::Parameter::INTEGER, "", "The SVM kernel to be used in developing model (1=Linear, 2=RBF).");
  parser.addRequiredParameter("n", "configuration", cbica::Parameter::INTEGER, "", "The Configuration type, either Cross-validation or Split Train-Test (1=CrossValidtion, 2=TrainTest).");

  parser.addRequiredParameter("k", "configuration parameters", cbica::Parameter::INTEGER, "", "The number of folds for Crossvalidation (5/10) and the size of training set for TrainTest (k<n).");
  parser.addRequiredParameter("o", "output", cbica::Parameter::STRING, "", "The output direcory to write output");
  parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");

  // parameters to get from the command line
  cbica::Logging logger;
  std::string loggerFile;
  bool loggerRequested = false;

  int tempPosition;
  std::string inputFeaturesFile, inputLabelsFile, outputDirectoryName, toWrite;
  int classifierType;
  int foldType;
  int confType;
  if ((argc < 1) || (parser.compareParameter("u", tempPosition)))
  {
    parser.echoUsage();
    return EXIT_SUCCESS;
  }
  if (parser.compareParameter("h", tempPosition))
  {
    parser.echoHelp();
    return EXIT_SUCCESS;
  }
  if (parser.compareParameter("v", tempPosition))
  {
    parser.echoVersion();
    return EXIT_SUCCESS;
  }
  if (parser.compareParameter("L", tempPosition))
  {
    loggerFile = argv[tempPosition + 1];
    loggerRequested = true;
    logger.UseNewFile(loggerFile);
  }
  // should checks for nifti files be placed?
  if (parser.compareParameter("f", tempPosition))
  {
    inputFeaturesFile = argv[tempPosition + 1];
    //inputFeaturesFile = "W:/Projects/PSU/PSU_Selected_Features_WithoutScaling.csv";
    //inputFeaturesFile = "W:/Projects/PSU/REC_Selected_Features_NewSheet_WithoutScaling.csv";
  }
  if (parser.compareParameter("l", tempPosition))
  {
    inputLabelsFile = argv[tempPosition + 1];
    //inputLabelsFile = "W:/Projects/PSU/labels_pseudoprogression.csv";
    ////inputLabelsFile = "W:/Projects/PSU/labels_recurrence.csv";
  }
  if (parser.compareParameter("o", tempPosition))
  {
    outputDirectoryName = argv[tempPosition + 1];
  }
  if (parser.compareParameter("c", tempPosition))
  {
    classifierType = atoi(argv[tempPosition + 1]);
  }
  if (parser.compareParameter("k", tempPosition))
  {
    foldType = atoi(argv[tempPosition + 1]);
  }
  if (parser.compareParameter("n", tempPosition))
  {
    confType = atoi(argv[tempPosition + 1]);
  }
  TrainingModule mTrainingSimulator;
  if(mTrainingSimulator.Run(inputFeaturesFile, inputLabelsFile, outputDirectoryName, classifierType, foldType, confType)==true)
    std::cout << "Finished successfully!!!\n";
  else
    std::cout << "Encountered an error!!!\n";

  int a;
  std::cin >> a;
  return EXIT_SUCCESS;
}
