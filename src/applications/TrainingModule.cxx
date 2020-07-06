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
  parser.addRequiredParameter("o", "output", cbica::Parameter::STRING, "", "The output direcory to write output");

  parser.addOptionalParameter("l", "label", cbica::Parameter::STRING, "", "The input file having target labels (*.csv).");
  parser.addOptionalParameter("c", "classifier", cbica::Parameter::INTEGER, "", "The SVM kernel to be used in developing model (1=Linear, 2=RBF).");
  parser.addOptionalParameter("s", "feature selection", cbica::Parameter::INTEGER, "", "The feature selection method to be used in developing model (1=EffectSize, 2=Correlation, 3=SVM FFS, 4=SVM RFE).");
  parser.addOptionalParameter("n", "configuration", cbica::Parameter::INTEGER, "", "The Configuration type, Cross-validation (n=1), Split Train-Test (n=2), Train only (n=3), and Test only (n=4).");
  parser.addOptionalParameter("k", "configuration parameters", cbica::Parameter::INTEGER, "", "The number of folds for Cross-validation (5/10) and the size of training set for TrainTest (k<n).");
  parser.addOptionalParameter("p", "hyperparameteroptimization", cbica::Parameter::INTEGER, "", "Whether parameters of the classifier need to be optimized or not during feature selection (1=yes, 0 =No)");
  parser.addOptionalParameter("r", "internalcrossvalidation", cbica::Parameter::INTEGER, "", "Internal cross-validation during feature selection (1=resubstitution, 2=5-fold)");

  parser.addOptionalParameter("m", "output", cbica::Parameter::STRING, "", "The model direcory (needed only when n=4)");
  parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");
  //parser.exampleUsage("TrainingModule -f features2.csv -l labels2.csv -c 1 -o <output dir> -k 5");
  parser.addExampleUsage(" -f features2.csv -l labels2.csv -c 1 -o <output dir> -k 5", 
    "Trains a new Linear SVM model based on the input features in 'feature2.csv' and corresponding labels in 'labels2.csv' with cross-validation of 5");
  parser.addApplicationDescription("Molecular Subtype Training and Prediction application");

  // parameters to get from the command line
  cbica::Logging logger;
  std::string loggerFile;
  bool loggerRequested = false;

  int tempPosition;
  std::string inputFeaturesFile, inputLabelsFile, outputDirectoryName, modelDirectoryName, toWrite;

  //assigning default values
  int classifierType=1;
  int featureselectionType=1;
  int optimizationType=0;
  int crossvalidationType=1;
  int foldType=10;
  int confType=1;

  TrainingModule mTrainingSimulator;
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
    std::cout << "Input Features File:"<<inputFeaturesFile << std::endl;
  }
  if (parser.compareParameter("l", tempPosition))
  {
    inputLabelsFile = argv[tempPosition + 1];
    std::cout << "Input Labels File:" << inputLabelsFile << std::endl;
  }
  if (parser.compareParameter("o", tempPosition))
  {
    outputDirectoryName = argv[tempPosition + 1];
    cbica::createDir(outputDirectoryName);
  }
  if (parser.compareParameter("m", tempPosition))
  {
    modelDirectoryName = argv[tempPosition + 1];
  }
  if (parser.compareParameter("c", tempPosition))
  {
    classifierType = atoi(argv[tempPosition + 1]);
  }
  if (parser.compareParameter("s", tempPosition))
  {
    featureselectionType = atoi(argv[tempPosition + 1]);
  }
  if (parser.compareParameter("k", tempPosition))
  {
    foldType = atoi(argv[tempPosition + 1]);
  }
  if (parser.compareParameter("n", tempPosition))
  {
    confType = atoi(argv[tempPosition + 1]);
  }
  if (parser.compareParameter("p", tempPosition))
  {
    optimizationType = atoi(argv[tempPosition + 1]);
  }
  if (parser.compareParameter("r", tempPosition))
  {
    crossvalidationType = atoi(argv[tempPosition + 1]);
  }
  //TrainingModule mTrainingSimulator;
  std::cout << "Calling function" << std::endl;

  if (confType == CAPTK::ClassificationConfigurationType::CONF_TYPE_SPLIT_TRAIN |
    confType == CAPTK::ClassificationConfigurationType::CONF_TYPE_DOUBLE |
    confType == CAPTK::ClassificationConfigurationType::CONF_TYPE_KFOLD_CV)
  {
    if (inputLabelsFile == "")
    {
      std::cout << "Please provide the class label file." << std::endl;
      return EXIT_FAILURE;
    }
  }
  else if (confType == CAPTK::ClassificationConfigurationType::CONF_TYPE_SPLIT_TEST)
  {
    if (modelDirectoryName == "")
    {
      std::cout << "Please provide the model directory name." << std::endl;
      return EXIT_FAILURE;
    }
  }
  if (mTrainingSimulator.Run(inputFeaturesFile, outputDirectoryName, inputLabelsFile, modelDirectoryName, classifierType, foldType, confType,featureselectionType,optimizationType, crossvalidationType) == true)
    std::cout << "Finished successfully!!!\n";
  else
    std::cout << "Encountered an error!!!\n";

  //int a;
  //std::cin >> a;
  return EXIT_SUCCESS;
}
