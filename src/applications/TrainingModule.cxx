#include "TrainingModule.h"
#include <ctime>
#include "cbicaCmdParser.h"
#include "CAPTk.h"

void showVersionInfo()
{
  //cout << "Survival Predictor: version: " << VERSION_STRING_W << endl;
}
void showhelpinfo(const std::string exeName)
{
  cout << "Usage:   " << exeName << " [-option] [argument]" << endl;
  cout << "Option:  " << "-h show help information" << endl;
  cout << "         " << "-v show version information" << endl;
  cout << "         " << "-f feature file path" << endl;
  cout << "         " << "-l label file path" << endl;
  cout << "         " << "-o output directory path" << endl;
  cout << "Example: " << exeName << " -f ../features.csv -l ../targets.csv -o ../outputfolder" << endl;
}

int main(int argc, char *argv[])
{
  cbica::CmdParser parser = cbica::CmdParser(argc, argv, "TrainingModule");
  parser.addRequiredParameter("f", "features", cbica::Parameter::STRING, "", "The input file having features (*.csv).");
  parser.addRequiredParameter("l", "label", cbica::Parameter::STRING, "", "The input file having target labels (*.csv).");
  parser.addRequiredParameter("c", "classifier", cbica::Parameter::INTEGER, "", "The SVM kernel to be used in developing model (1=Linear, 2=RBF).");
  parser.addRequiredParameter("k", "No. of folds", cbica::Parameter::INTEGER, "", "The number of folds to develop the model (5/10).");
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
  }

  if (parser.compareParameter("l", tempPosition))
  {
    inputLabelsFile = argv[tempPosition + 1];
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
  // std::cout << inputFeaturesFile << std::endl;
  // std::cout << inputLabelsFile << std::endl;
  // std::cout << "Classifier: " << std::to_string(classifierType) << " , Folds: " << std::to_string(foldType)<< std::endl;
  
  // cbica::Logging(loggerFile, inputFeaturesFile + "\n");
  // cbica::Logging(loggerFile, inputLabelsFile + "\n");
  // cbica::Logging(loggerFile, "Classifier: " + std::to_string(classifierType) + " , Folds: " + std::to_string(foldType) + "\n");

  TrainingModule mTrainingSimulator;
  mTrainingSimulator.Run(inputFeaturesFile, inputLabelsFile, outputDirectoryName,classifierType,foldType);


  std::cout << "Finished successfully\n";

  return EXIT_SUCCESS;
}
