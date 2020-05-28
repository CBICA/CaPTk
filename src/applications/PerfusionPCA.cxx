#include "PerfusionPCA.h"
#include "cbicaUtilities.h"
#include "cbicaCmdParser.h"
#include "cbicaITKSafeImageIO.h"

std::vector<std::map<CAPTK::ImageModalityType, std::string>> LoadQualifiedSubjectsFromGivenDirectoryForPCA(const std::string directoryname)
{
  std::map<CAPTK::ImageModalityType, std::string> OneQualifiedSubject;
  std::vector<std::map<CAPTK::ImageModalityType, std::string>> QualifiedSubjects;
  std::vector<std::string> subjectNames = cbica::subdirectoriesInDirectory(directoryname);
  std::sort(subjectNames.begin(), subjectNames.end());

  for (unsigned int sid = 0; sid < subjectNames.size(); sid++)
  {
    std::string subjectPath = directoryname + "/" + subjectNames[sid];
    
    std::string perfFilePath = "";
    std::string labelPath = "";

    std::vector<std::string> files;
    files = cbica::filesInDirectory(subjectPath + "", false);

    for (unsigned int i = 0; i < files.size(); i++)
    {
      std::string filePath = subjectPath + "/" + files[i], filePath_lower;
      std::string extension = cbica::getFilenameExtension(filePath, false);
      filePath_lower = filePath;
      std::transform(filePath_lower.begin(), filePath_lower.end(), filePath_lower.begin(), ::tolower);
      if ((filePath_lower.find("label") != std::string::npos || filePath_lower.find("segmentation") != std::string::npos)
        && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
        labelPath = subjectPath + "/" + files[i];
      else if ((filePath_lower.find("perf") != std::string::npos || filePath_lower.find("PERF") != std::string::npos || filePath_lower.find("DSC") != std::string::npos)
        && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
        perfFilePath = subjectPath +  "/" + files[i];
    }

    if (labelPath.empty() || perfFilePath.empty())
      continue;

    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION] = perfFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG] = labelPath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID] = subjectNames[sid];

    QualifiedSubjects.push_back(OneQualifiedSubject);
  }
  return QualifiedSubjects;
}


int main(int argc, char **argv)
{
  cbica::CmdParser parser = cbica::CmdParser(argc, argv, "PerfusionPCA");
  parser.addRequiredParameter("i", "input", cbica::Parameter::STRING, "", "The input directory.");
  parser.addRequiredParameter("t", "type", cbica::Parameter::STRING, "", "The option of preparing a new model (=0), and for testing on an existing model (=1)");
  parser.addRequiredParameter("n", "number of PCAs", cbica::Parameter::STRING, "", "The number of principal components.");
  parser.addOptionalParameter("m", "model", cbica::Parameter::STRING, "", "The directory having PCA models");
  parser.addRequiredParameter("o", "output", cbica::Parameter::STRING, "", "The output directory.");
  parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");
  //parser.exampleUsage("");
  parser.addExampleUsage("-t 0 -i C:/properly/formatted/inputDir -o C:/outputDir -n 5", "Trains a new model based on the samples in inputDir");
  parser.addExampleUsage("-t 1 -i C:/input -m C:/model -o C:/output -m C:/modelDir -n 5", "Tests an existing model for inputs in 'C:/input' based on 'C:/modelDir' ");
  parser.addApplicationDescription("Perfusion PCA Training and Prediction application");

  // parameters to get from the command line
  cbica::Logging logger;
  std::string loggerFile;
  bool loggerRequested = false;

  int tempPosition;
  int applicationType;
  applicationType = 0;

  double inputPCs = 0;
  std::string inputFileName, inputMaskName, outputDirectoryName, modelDirectoryName;

  if ((argc < 1) || (parser.compareParameter("u", tempPosition)))
  {
    parser.echoUsage();
    return EXIT_SUCCESS;
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
  if (parser.compareParameter("n", tempPosition))
  {
	  inputPCs = atof(argv[tempPosition + 1]);
  }
  if (parser.compareParameter("o", tempPosition))
  {
    outputDirectoryName = argv[tempPosition + 1];
  }
  if (parser.compareParameter("m", tempPosition))
  {
    modelDirectoryName = argv[tempPosition + 1];
  }
  if (parser.compareParameter("t", tempPosition))
  {
    applicationType = atoi(argv[tempPosition + 1]);
  }

  std::cout << "Input File:" << inputFileName << std::endl;
  std::cout << "Output Directory:" << outputDirectoryName << std::endl;

  if (!cbica::isDir(inputFileName))
  {
    std::cout << "The input directory does not exist:" << inputFileName << std::endl;
    return EXIT_FAILURE;
  }
  if (!cbica::directoryExists(outputDirectoryName))
  {
    if (!cbica::createDirectory(outputDirectoryName))
      std::cout << "The output directory can not be created:" << outputDirectoryName << std::endl;
    return EXIT_FAILURE;
  }
  std::vector<std::map<CAPTK::ImageModalityType, std::string>> QualifiedSubjects = LoadQualifiedSubjectsFromGivenDirectoryForPCA(inputFileName);

  if (QualifiedSubjects.size() == 0)
  {
    std::cout << "There is no subject with the required input in the given directory." << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Number of subjects with the required input: " << QualifiedSubjects.size() << std::endl;
  PerfusionPCA object_pca;
  if (applicationType == CAPTK::MachineLearningApplicationSubtype::TESTING)
  {
    std::cout << "Model directory name:" << modelDirectoryName << std::endl;
    if (!cbica::directoryExists(modelDirectoryName))
    {
      std::cout << "The model directory does not exist:" << modelDirectoryName << std::endl;
      return EXIT_FAILURE;
    }
    object_pca.ApplyExistingPCAModel(inputPCs, inputFileName, outputDirectoryName, QualifiedSubjects,modelDirectoryName);
  }
  else if (applicationType == CAPTK::MachineLearningApplicationSubtype::TRAINING)
    object_pca.TrainNewPerfusionModel(inputPCs, inputFileName, outputDirectoryName, QualifiedSubjects);
  else
  {
    parser.echoVersion();
    return EXIT_SUCCESS;
  }

  std::cout<<"principal components have been saved at the specified locations.\n"; 
  std::cout << "Finished successfully.\n";

  return EXIT_SUCCESS;
}