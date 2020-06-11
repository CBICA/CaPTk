#include "RecurrenceEstimator.h"
#include "cbicaUtilities.h"
#include "cbicaCmdParser.h"

#include "CaPTkGUIUtils.h"

//------------------Survival Prediction on existing model-----------------------
std::vector<std::map<CAPTK::ImageModalityType, std::string>> LoadQualifiedSubjectsFromGivenDirectoryForRecurrence(const CAPTK::MachineLearningApplicationSubtype type, const std::string &directoryname, const bool &useConventionalData, const bool &useDTIData, const bool &usePerfData, const bool &useDistData)
{
  std::map<CAPTK::ImageModalityType, std::string> OneQualifiedSubject;
  std::vector<std::map<CAPTK::ImageModalityType, std::string>> QualifiedSubjects;
  std::vector<std::string> subjectNames = cbica::subdirectoriesInDirectory(directoryname);
  std::sort(subjectNames.begin(), subjectNames.end());

  for (unsigned int sid = 0; sid < subjectNames.size(); sid++)
  {
    std::string subjectPath = directoryname + "/" + subjectNames[sid];
    std::string t1ceFilePath = "";
    std::string t1FilePath = "";
    std::string t2FilePath = "";
    std::string t2FlairFilePath = "";
    std::string axFilePath = "";
    std::string faFilePath = "";
    std::string radFilePath = "";
    std::string trFilePath = "";
    std::string perfFilePath = "";
    std::string labelPath = "";
    std::string nearFilePath = "";
    std::string farFilePath = "";

    std::vector<std::string> files;

    if (cbica::directoryExists(subjectPath + "/DRAWING"))
    {
      files = cbica::filesInDirectory(subjectPath + "/DRAWING", false);
      for (unsigned int i = 0; i < files.size(); i++)
      {
        std::string filePath = subjectPath + "/DRAWING/" + files[i];
        std::string extension = cbica::getFilenameExtension(filePath, false);
        if ((files[i].find("near") != std::string::npos || files[i].find("Infiltrated") != std::string::npos) && (extension == IMG_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          nearFilePath = subjectPath + "/DRAWING/" + files[i];

        if ((files[i].find("far") != std::string::npos || files[i].find("Pure") != std::string::npos) && (extension == IMG_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          farFilePath = subjectPath + "/DRAWING/" + files[i];
      }
    }

    if (cbica::directoryExists(subjectPath + "/SEGMENTATION"))
    {
      files = cbica::filesInDirectory(subjectPath + "/SEGMENTATION", false);
      for (unsigned int i = 0; i < files.size(); i++)
      {
        std::string filePath = subjectPath + "/SEGMENTATION/" + files[i];
        std::string extension = cbica::getFilenameExtension(filePath, false);
        if ((files[i].find("label-map") != std::string::npos || files[i].find("label") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          labelPath = subjectPath + "/SEGMENTATION/" + files[i];
      }
    }
    if (cbica::directoryExists(subjectPath + "/PERFUSION"))
    {
      files = cbica::filesInDirectory(subjectPath + "/PERFUSION", false);
      for (unsigned int i = 0; i < files.size(); i++)
      {
        std::string filePath = subjectPath + "/PERFUSION/" + files[i];
        std::string extension = cbica::getFilenameExtension(filePath, false);
        if ((files[i].find("perf") != std::string::npos || files[i].find("Perf") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          perfFilePath = subjectPath + "/PERFUSION/" + files[i];
      }
    }
    if (useConventionalData && cbica::directoryExists(subjectPath + "/CONVENTIONAL"))
    {
      files = cbica::filesInDirectory(subjectPath + "/CONVENTIONAL", false);
      for (unsigned int i = 0; i < files.size(); i++)
      {
        std::string filePath = subjectPath + "/CONVENTIONAL/" + files[i];
        std::string extension = cbica::getFilenameExtension(filePath, false);

        if ((files[i].find("t1ce") != std::string::npos || files[i].find("T1CE") != std::string::npos || files[i].find("T1ce") != std::string::npos || files[i].find("T1-gd") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          t1ceFilePath = subjectPath + "/CONVENTIONAL/" + files[i];
        else if ((files[i].find("t1") != std::string::npos || files[i].find("T1") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          t1FilePath = subjectPath + "/CONVENTIONAL/" + files[i];
        else if ((files[i].find("t2") != std::string::npos || files[i].find("T2") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          t2FilePath = subjectPath + "/CONVENTIONAL/" + files[i];
        else if ((files[i].find("flair") != std::string::npos || files[i].find("FLAIR") != std::string::npos || files[i].find("Flair") != std::string::npos || files[i].find("T2-Flair") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          t2FlairFilePath = subjectPath + "/CONVENTIONAL/" + files[i];
      }
    }


    if (useDTIData && cbica::directoryExists(subjectPath + "/DTI"))
    {
      files = cbica::filesInDirectory(subjectPath + "/DTI", false);
      for (unsigned int i = 0; i < files.size(); i++)
      {
        std::string filePath = subjectPath + "/DTI/" + files[i];
        std::string extension = cbica::getFilenameExtension(filePath, false);

        if ((files[i].find("Axial") != std::string::npos || files[i].find("axial") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          axFilePath = subjectPath + "/DTI/" + files[i];
        else if ((files[i].find("Fractional") != std::string::npos || files[i].find("fractional") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          faFilePath = subjectPath + "/DTI/" + files[i];
        else if ((files[i].find("Radial") != std::string::npos || files[i].find("radial") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          radFilePath = subjectPath + "/DTI/" + files[i];
        else if ((files[i].find("Trace") != std::string::npos || files[i].find("trace") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          trFilePath = subjectPath + "/DTI/" + files[i];
      }
    }
    if ((useConventionalData && t1FilePath.empty()) || (useConventionalData && t2FilePath.empty()) || (useConventionalData && t1ceFilePath.empty()) || (useConventionalData && t2FlairFilePath.empty()) ||
      (usePerfData && perfFilePath.empty()) || (useDTIData && axFilePath.empty()) || (useDTIData && faFilePath.empty()) || (useDTIData && radFilePath.empty()) || (useDTIData && trFilePath.empty()))
      continue;

    if ((nearFilePath.empty() || farFilePath.empty()) && type == CAPTK::MachineLearningApplicationSubtype::TRAINING)
      continue;

    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1] = t1FilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2] = t2FilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1CE] = t1ceFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR] = t2FlairFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_AX] = axFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_FA] = faFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_RAD] = radFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_TR] = trFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION] = perfFilePath;

    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG] = labelPath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID] = subjectNames[sid];

    if (type == CAPTK::MachineLearningApplicationSubtype::TRAINING)
    {
      OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_NEAR] = nearFilePath;
      OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_FAR] = farFilePath;
    }
    QualifiedSubjects.push_back(OneQualifiedSubject);

  }
  return QualifiedSubjects;
}


int RecurrencePredictionOnExistingModel(const std::string modeldirectory,
  const std::string inputdirectory,
  const std::string outputdirectory)
{
  std::cout << "Module loaded: Recurrence Prediction using Existing Model." << std::endl;
  std::vector<double> finalresult;
  std::vector<std::map<CAPTK::ImageModalityType, std::string>> QualifiedSubjects = LoadQualifiedSubjectsFromGivenDirectoryForRecurrence(CAPTK::MachineLearningApplicationSubtype::TESTING, inputdirectory, true, true, true, true);
  std::cout << "Number of subjects with required input: " << QualifiedSubjects.size() << std::endl;
  RecurrenceEstimator objRecurrenceEstimator;
  for (unsigned int subjectID = 0; subjectID < QualifiedSubjects.size(); subjectID++)
  {
    std::map<CAPTK::ImageModalityType, std::string> onesubject = QualifiedSubjects[subjectID];
    std::vector<std::map<CAPTK::ImageModalityType, std::string>> OneQualifiedSubject;
    OneQualifiedSubject.push_back(QualifiedSubjects[subjectID]);
    std::cout << "Patient's data loaded: " << subjectID + 1 << "/" << QualifiedSubjects.size() << ", Sudo ID: " << static_cast<std::string>(OneQualifiedSubject[0][CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]) << std::endl;
    objRecurrenceEstimator.RecurrenceEstimateOnExistingModel(OneQualifiedSubject, modeldirectory, inputdirectory, outputdirectory, true, true, true, true);
  }
  return EXIT_SUCCESS;
}
int PrepareNewRecurrencePredictionModel(const std::string inputdirectory, const std::string outputdirectory)
{
  std::cout << "Module loaded: Prepare Recurrence Prediction Model." << std::endl;
  std::vector<double> finalresult;
  std::vector<std::map<CAPTK::ImageModalityType, std::string>> QualifiedSubjects = LoadQualifiedSubjectsFromGivenDirectoryForRecurrence(CAPTK::MachineLearningApplicationSubtype::TRAINING, inputdirectory, true, true, true, true);
  RecurrenceEstimator objRecurrencePredictor;
  std::cout << "Number of subjects with required input: " << QualifiedSubjects.size() << std::endl;
  objRecurrencePredictor.TrainNewModelOnGivenData(QualifiedSubjects, outputdirectory, true, true, true, true);
  return EXIT_SUCCESS;
}
int main(int argc, char **argv)
{
  cbica::CmdParser parser = cbica::CmdParser(argc, argv, "RecurrenceEstimator");
  parser.addRequiredParameter("t", "type", cbica::Parameter::STRING, "", "The option of preparing a new model (=0), and for testing on an existing model (=1)");
  parser.addRequiredParameter("i", "input", cbica::Parameter::STRING, "", "The input directory having test subjects");
  parser.addOptionalParameter("m", "model", cbica::Parameter::STRING, "", "The directory having SVM models", "Penn Model: " + getAppropriateDownloadLink("RecurrenceEstimator", "Model"));
  parser.addRequiredParameter("o", "output", cbica::Parameter::STRING, "", "The output direcory to write output");
  parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");
  //parser.exampleUsage("RecurrenceEstimator -t 0 -i <input dir> -o <output dir>");
  parser.addExampleUsage("-t 0 -i C:/properly/formatted/inputDir -o C:/outputDir", "Trains a new model based on the samples in inputDir");
  parser.addExampleUsage("-t 1 -i C:/input -m C:/model -o C:/output", "Tests an existing model for inputs in 'C:/input' based on 'C:/model' ");
  parser.addApplicationDescription("Recurrence Estimator Training and Prediction application");

  //parser.addOptionalParameter("c", "T1ce usage", cbica::Parameter::STRING, "", "Whether to use conventional imaging for model building or not? 1 for YES, 2 for NO");
  //parser.addOptionalParameter("p", "Perfusion usage", cbica::Parameter::STRING, "", "Whether to use perfusion imaging for model building or not? 1 for YES, 2 for NO");
  //parser.addOptionalParameter("d", "DTI usage", cbica::Parameter::STRING, "", "Whether to use DTI imaging for model building or not? 1 for YES, 2 for NO");
  //parser.addOptionalParameter("s", "Distance usage", cbica::Parameter::STRING, "", "Whether to use distance measure for model building or not? 1 for YES, 2 for NO");

  // parameters to get from the command line
  cbica::Logging logger;
  std::string loggerFile;
  bool loggerRequested = false;

  int tempPosition;
  std::string inputDirectoryName, modelDirectoryName, outputDirectoryName, toWrite;
  int applicationType;
  //int useConventional;
  //int usePerfusion;
  //int useDTI;
  //int useDistance;
  
  if (parser.compareParameter("L", tempPosition))
  {
    loggerFile = argv[tempPosition + 1];
    loggerRequested = true;
    logger.UseNewFile(loggerFile);
  }
  if (parser.compareParameter("i", tempPosition))
  {
    inputDirectoryName = argv[tempPosition + 1];
  }

  if (parser.compareParameter("m", tempPosition))
  {
    modelDirectoryName = argv[tempPosition + 1];
  }

  if (parser.compareParameter("o", tempPosition))
  {
    outputDirectoryName = argv[tempPosition + 1];
  }
  if (parser.compareParameter("t", tempPosition))
  {
    applicationType = atoi(argv[tempPosition + 1]);
  }

  //if (parser.compareParameter("c", tempPosition))
  //{
  //	if (atoi(argv[tempPosition + 1]) == 1)
  //		useConventional = true;
  //	else if (atoi(argv[tempPosition + 1]) == 2)
  //		useConventional = false;
  //	else
  //	{
  //		std::cout << "Please specify whether to use conventional imaging." << std::endl;
  //		return EXIT_FAILURE;
  //	}
  //}
  //if (parser.compareParameter("d", tempPosition))
  //{
  //	if (atoi(argv[tempPosition + 1]) == 1)
  //		useDTI = true;
  //	else if (atoi(argv[tempPosition + 1]) == 2)
  //		useDTI = false;
  //	else
  //	{
  //		std::cout << "Please specify whether to use DTI imaging." << std::endl;
  //		return EXIT_FAILURE;
  //	}
  //}
  //if (parser.compareParameter("s", tempPosition))
  //{
  //	if (atoi(argv[tempPosition + 1]) == 1)
  //		useDistance = true;
  //	else if (atoi(argv[tempPosition + 1]) == 2)
  //		useDistance = false;
  //	else
  //	{
  //		std::cout << "Please specify whether to use Distance measure." << std::endl;
  //		return EXIT_FAILURE;
  //	}
  //}
  //if (parser.compareParameter("p", tempPosition))
  //{
  //	if (atoi(argv[tempPosition + 1]) == 1)
  //		usePerfusion = true;
  //	else if (atoi(argv[tempPosition + 1]) == 2)
  //		usePerfusion = false;
  //	else
  //	{
  //		std::cout << "Please specify whether to use Perfusion imaging." << std::endl;
  //		return EXIT_FAILURE;
  //	}
  //}
  if (applicationType == CAPTK::MachineLearningApplicationSubtype::TESTING && modelDirectoryName.empty())
  {
    std::cout << "Please specify a directory having model file." << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Input directory name:" << inputDirectoryName << std::endl;
  std::cout << "Output directory name:" << outputDirectoryName << std::endl;

  if (!cbica::directoryExists(outputDirectoryName))
  {
    if (!cbica::createDirectory(outputDirectoryName))
    {
      std::cout << "The output directory can not be created:" << outputDirectoryName << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (applicationType == CAPTK::MachineLearningApplicationSubtype::TESTING)
  {
    std::cout << "Model directory name:" << modelDirectoryName << std::endl;
    if (!cbica::directoryExists(modelDirectoryName))
    {
      std::cout << "The model directory does not exist:" << modelDirectoryName << std::endl;
      return EXIT_FAILURE;
    }
    if (cbica::isFile(modelDirectoryName + "/VERSION.yaml"))
    {
        if (!cbica::IsCompatible(modelDirectoryName + "/VERSION.yaml"))
        {
            std::cerr << "The version of model is incompatible with this version of CaPTk.\n";
            return EXIT_FAILURE;
        }
    }
    RecurrencePredictionOnExistingModel(modelDirectoryName, inputDirectoryName, outputDirectoryName);
  }
  else if (applicationType == CAPTK::MachineLearningApplicationSubtype::TRAINING)
    PrepareNewRecurrencePredictionModel(inputDirectoryName, outputDirectoryName);
  else
  {
    parser.echoVersion();
    return EXIT_SUCCESS;
  }

  std::cout << "Finished successfully\n";

  return EXIT_SUCCESS;
}