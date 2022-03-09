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
  parser.addRequiredParameter("o", "output", cbica::Parameter::STRING, "", "The directory where to create output");

  parser.addOptionalParameter("l", "label", cbica::Parameter::STRING, "", "The input file having target labels (*.csv).");
  parser.addOptionalParameter("c", "classifier", cbica::Parameter::INTEGER, "", 
      "The classifier to be used in developing the model. \
      (1=Linear SVM, 2=RBF SVM, 3=Polynomial SVM, 4=Sigmoid SVM, 5=Chi-squared SVM, \
      6=Intersection SVM, 7=Random Forest, 8=SGD-based SVM, 9=Boosted Trees). Default: Linear SVM");
  parser.addOptionalParameter("s", "featureselection", cbica::Parameter::INTEGER, "",
      "The feature selection method to be used in developing the model \
       (1=Effect-size FS, 2=Forward FS, 3=Recursive Feature Elimination, 4=Random-Forest-based FS, \
       5=RELIEF-F FS). Default: Effect-size FS");
  parser.addOptionalParameter("t", "task", cbica::Parameter::STRING, "", "Execution mode. Valid values are 'crossvalidate', 'train', 'test'.");
  parser.addOptionalParameter("k", "kfolds", cbica::Parameter::INTEGER, "", "The number of folds for Cross-validation (default 5).");
  parser.addOptionalParameter("p", "hyperparameteroptimization", cbica::Parameter::INTEGER, "", "Whether hyperparameters of the classifier need to be optimized or not during feature selection (1=yes, 2 =No) Default: No");
  parser.addOptionalParameter("r", "internalcrossvalidation", cbica::Parameter::INTEGER, "", "Internal cross-validation during feature selection (1=resubstitution, 2=5-fold)");
  parser.addOptionalParameter("x", "maxfeatures", cbica::Parameter::INTEGER, "", 
      "Maximum number of features to select. Up to this many features will be included.  \
      A value of 0 behaves differently between methods. For Forward and Effect-size FS and Recursive Feature Elimination, \
      produces the best performing feature set overall. For Random Forest and Relief-F, selects all features, but in the order of importance.");

  parser.addOptionalParameter("Cmax", "csearchmaximum", cbica::Parameter::FLOAT, "", "Log2 of the higher bound of the C hyperparameter search space (used for optimization of SVMs) Default: 5");
  parser.addOptionalParameter("Cmin", "csearchmaximum", cbica::Parameter::FLOAT, "", "Log2 of the lower bound of the C hyperparameter search space (used for optimization of SVMs) Default: -5");
  parser.addOptionalParameter("Gmax", "gsearchminimum", cbica::Parameter::FLOAT, "", "Log2 of the higher bound of the Gamma (G) hyperparameter search space (used for optimization of RBF, Polynomial, Sigmod and Chi-squared SVMs) Default 5");
  parser.addOptionalParameter("Gmin", "gsearchmaximum", cbica::Parameter::FLOAT, "", "Log2 of the lower bound of the Gamma (G) hyperparameter search space (used for optimization of RBF, Polynomial, Sigmod and Chi-squared SVMs) Default -5");
  // TODO: add grid-search parametrization for the other classification-problem SVM params (Coef0, Degree)
  // P, Nu may be necessary for the SVMs if we extend beyond the classification problem (C_SVC)

  parser.addOptionalParameter("RFeps", "rfepsilon", cbica::Parameter::FLOAT, "", "Random Forest: Determine a convergence threshold for termination. Default: 0.01");
  parser.addOptionalParameter("RFmaxiters", "rfmaxiterations", cbica::Parameter::INTEGER, "", "Random Forest: Determine maximum number of forest iterations before termination. Default: 50");

  parser.addOptionalParameter("m", "model", cbica::Parameter::STRING, "", "The model directory (needed only when n=3)");
  parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");
  //parser.exampleUsage("TrainingModule -f features2.csv -l labels2.csv -c 1 -o <output dir> -k 5");
  parser.addExampleUsage(" -f features2.csv -l labels2.csv -c 1 -s 1 -o <output dir> -k 5", 
    "Trains a new Linear SVM model based on the input features in 'feature2.csv' and corresponding labels in 'labels2.csv' with 5 cross-validation folds");
  parser.addApplicationDescription("Training Module");

  // parameters to get from the command line
  cbica::Logging logger;
  std::string loggerFile;
  bool loggerRequested = false;

  int tempPosition;
  std::string inputFeaturesFile, inputLabelsFile, outputDirectoryName, modelDirectoryName, toWrite;

  //assigning default values -- TODO: grab these from params defaults so they're in a central place
  int classifierType=1;
  int featureselectionType=1;
  int optimizationType=2;
  int crossvalidationType=1;
  int foldType=5;
  int confType=0;
  std::string confTypeString;

  double cMin = -5;
  double cMax = 5;
  double gMin = -5;
  double gMax = 5;

  TrainingModule mTrainingSimulator;
  TrainingModuleParameters params; // Track parameters to be passed to TrainingModule with this

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
    params.inputFeaturesFile = inputFeaturesFile;
    std::cout << "Input Features File:"<<inputFeaturesFile << std::endl;
  }
  if (parser.compareParameter("l", tempPosition))
  {
    inputLabelsFile = argv[tempPosition + 1];
    params.inputLabelsFile = inputLabelsFile;
    std::cout << "Input Labels File:" << inputLabelsFile << std::endl;
  }
  if (parser.compareParameter("o", tempPosition))
  {
    outputDirectoryName = argv[tempPosition + 1];
    params.outputDirectory = outputDirectoryName;
    cbica::createDir(outputDirectoryName);
  }
  if (parser.compareParameter("m", tempPosition))
  {
    modelDirectoryName = argv[tempPosition + 1];
    params.modelDirectory = modelDirectoryName;
    std::cout << "Model directory:" << inputLabelsFile << std::endl;
  }
  if (parser.compareParameter("c", tempPosition))
  {
    classifierType = atoi(argv[tempPosition + 1]);
    params.classifierType = (CAPTK::ClassifierType)classifierType;
  }
  if (parser.compareParameter("s", tempPosition))
  {
    featureselectionType = atoi(argv[tempPosition + 1]);
    params.featureSelectionType = (CAPTK::FeatureSelectionType)featureselectionType;
  }
  if (parser.compareParameter("k", tempPosition))
  {
    foldType = atoi(argv[tempPosition + 1]);
    params.folds = foldType;
  }
  if (parser.compareParameter("t", tempPosition))
  {
    confTypeString = std::string(argv[tempPosition + 1]);
    std::transform(confTypeString.begin(), confTypeString.end(), confTypeString.begin(), ::tolower);

    if (confTypeString == "cv" || confTypeString == "crossvalidate" || confTypeString == "crossvalidation")
    {
        confType = CAPTK::ClassificationConfigurationType::CONF_TYPE_KFOLD_CV;
    }
    else if (confTypeString == "train" || confTypeString == "training")
    {
        confType = CAPTK::ClassificationConfigurationType::CONF_TYPE_SPLIT_TRAIN;
    }
    else if (confTypeString == "test" || confTypeString == "testing" || confTypeString == "inference")
    {
        confType = CAPTK::ClassificationConfigurationType::CONF_TYPE_SPLIT_TEST;
    }
    else
    {
        confType = CAPTK::ClassificationConfigurationType::CONF_TYPE_UNDEFINED;
    }
    
    params.configurationType = (CAPTK::ClassificationConfigurationType)confType;
  }
  if (parser.compareParameter("p", tempPosition))
  {
    optimizationType = atoi(argv[tempPosition + 1]);
    params.optimizationType = (CAPTK::OptimizationType)optimizationType;
  }
  if (parser.compareParameter("r", tempPosition))
  {
    crossvalidationType = atoi(argv[tempPosition + 1]);
    params.crossValidationType = (CAPTK::CrossValidationType)crossvalidationType;
  }
  if (parser.compareParameter("cMin", tempPosition))
  {
      cMin = atof(argv[tempPosition + 1]);
      params.cMin = cMin;
  }
  if (parser.compareParameter("cMax", tempPosition))
  {
      cMax = atof(argv[tempPosition + 1]);
      params.cMax = cMax;
  }
  if (parser.compareParameter("gMin", tempPosition))
  {
      gMin = atof(argv[tempPosition + 1]);
      params.gMin = gMin;
  }
  if (parser.compareParameter("gMax", tempPosition))
  {
      gMax = atof(argv[tempPosition + 1]);
      params.gMax = gMax;
  }
  //TrainingModule mTrainingSimulator;
  std::cout << "Calling function" << std::endl;

  if (confType == CAPTK::ClassificationConfigurationType::CONF_TYPE_SPLIT_TRAIN |
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

    if (inputLabelsFile != "")
    {
        std::cout << "Using the provided labels file " + inputLabelsFile + " as ground truth to generate accuracy metrics." << std::endl;
        params.testPredictionsAgainstProvidedLabels = true;
    }

  }

  TrainingModuleResult result = mTrainingSimulator.Run(params);
  if (result.success == true)
  {
      std::cout << "Finished successfully!" << std::endl;
      return EXIT_SUCCESS;
  }
  else
  {
      std::cout << "Encountered an error. Please check logs for additional information." << std::endl;
      std::cout << result.message << std::endl;
      return EXIT_FAILURE;
  }
}
