#include "PseudoProgressionEstimator.h"
#include "cbicaUtilities.h"
#include "cbicaCmdParser.h"
#include "CaPTkEnums.h"

//------------------Survival Prediction on existing model-----------------------

std::vector<std::map<CAPTK::ImageModalityType, std::string>>  LoadQualifiedSubjectsFromGivenDirectoryForPseudoProgression(const CAPTK::MachineLearningApplicationSubtype type, const std::string &directoryname, const bool &useConventionalData, const bool &useDTIData, const bool &usePerfData, const bool &useDistData)
{
  std::map<CAPTK::ImageModalityType, std::string> OneQualifiedSubject;
  std::vector<std::map<CAPTK::ImageModalityType, std::string>> QualifiedSubjects;
  std::vector<std::string> subjectNames = cbica::subdirectoriesInDirectory(directoryname);
  std::cout << subjectNames.size();
  std::sort(subjectNames.begin(), subjectNames.end());

  for (unsigned int sid = 0; sid < subjectNames.size(); sid++)
  {
    std::string subjectPath = directoryname + "/" + subjectNames[sid];

    std::string t1ceFilePath = "";
    std::string t1FilePath = "";
    std::string t2FilePath = "";
    std::string t2FlairFilePath = "";

    //std::string t1t1ceFilePath    = "";
    //std::string t2t2FlairFilePath = "";

    std::string axFilePath = "";
    std::string faFilePath = "";
    std::string radFilePath = "";
    std::string trFilePath = "";
    std::string phFilePath = "";
    std::string psrFilePath = "";
    std::string rcbvFilePath = "";
    std::string perfFilePath = "";
    std::string featuresFilePath = "";

    std::string labelPath = "";
    std::string atlasPath = "";

    std::vector<std::string> files;


    if (cbica::directoryExists(subjectPath + "/SEGMENTATION"))
    {
      files = cbica::filesInDirectory(subjectPath + "/SEGMENTATION", false);
      for (unsigned int i = 0; i < files.size(); i++)
      {
        std::string filePath = subjectPath + "/SEGMENTATION/" + files[i];
        std::string extension = cbica::getFilenameExtension(filePath, false);
        if ((files[i].find("label-map") != std::string::npos || files[i].find("label") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          labelPath = subjectPath + "/SEGMENTATION/" + files[i];
        else if ((files[i].find("atlas") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          atlasPath = subjectPath + "/SEGMENTATION/" + files[i];
      }
    }
    if (cbica::directoryExists(subjectPath + "/PERFUSION"))
    {
      files = cbica::filesInDirectory(subjectPath + "/PERFUSION", false);
      for (unsigned int i = 0; i < files.size(); i++)
      {
        std::string filePath = subjectPath + "/PERFUSION" + "/" + files[i], filePath_lower;
        std::string extension = cbica::getFilenameExtension(filePath, false);
        filePath_lower = filePath;
        std::transform(filePath_lower.begin(), filePath_lower.end(), filePath_lower.begin(), ::tolower);
        if ((filePath_lower.find("rcbv") != std::string::npos)
          && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          rcbvFilePath = subjectPath + "/PERFUSION" + "/" + files[i];
        else if ((filePath_lower.find("psr") != std::string::npos)
          && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          psrFilePath = subjectPath + "/PERFUSION" + "/" + files[i];
        else if ((filePath_lower.find("ph") != std::string::npos)
          && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          phFilePath = subjectPath + "/PERFUSION" + "/" + files[i];
        else if ((filePath_lower.find("perf") != std::string::npos)
          && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          perfFilePath = subjectPath + "/PERFUSION" + "/" + files[i];
      }
    }
    if (useConventionalData && cbica::directoryExists(subjectPath + "/CONVENTIONAL"))
    {
      files = cbica::filesInDirectory(subjectPath + "/CONVENTIONAL", false);
      for (unsigned int i = 0; i < files.size(); i++)
      {
        std::string filePath = subjectPath + "/CONVENTIONAL/" + files[i];
        std::string extension = cbica::getFilenameExtension(filePath, false);

        //if ((files[i].find("t1t1ce") != std::string::npos || files[i].find("T1T1CE") != std::string::npos || files[i].find("T1T1ce") != std::string::npos || files[i].find("T1-T1-gd") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
        //  t1t1ceFilePath = subjectPath + "/CONVENTIONAL/" + files[i];
        //else if ((files[i].find("t2t2flair") != std::string::npos || files[i].find("T2T2FLAIR") != std::string::npos || files[i].find("T2T2Flair") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
        //  t2t2FlairFilePath = subjectPath + "/CONVENTIONAL/" + files[i];
        //else
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
    if (cbica::fileExists(subjectPath + "/features.csv"))
      featuresFilePath = subjectPath + "/features.csv";

    if (labelPath.empty() || t1FilePath.empty() || t2FilePath.empty() || t1ceFilePath.empty() || t2FlairFilePath.empty() || rcbvFilePath.empty() || axFilePath.empty() || faFilePath.empty() || radFilePath.empty() || trFilePath.empty() || psrFilePath.empty() || phFilePath.empty())
      continue;

    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1] = t1FilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2] = t2FilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1CE] = t1ceFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR] = t2FlairFilePath;
    //OneQualifiedSubject[IMAGE_TYPE_T1T1CE]    = t1t1ceFilePath; 
    //OneQualifiedSubject[IMAGE_TYPE_T2FL]      = t2t2FlairFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_AX] = axFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_FA] = faFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_RAD] = radFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_TR] = trFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_PH] = phFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_PSR] = psrFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_RCBV] = rcbvFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION] = perfFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_FEATURES] = featuresFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG] = labelPath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID] = subjectNames[sid];
    QualifiedSubjects.push_back(OneQualifiedSubject);

    std::cout << subjectNames[sid] << std::endl;
  }
  return QualifiedSubjects;
}
int SurvivalPredictionOnExistingModel(const std::string modeldirectory,
  const std::string inputdirectory,
  const std::string outputdirectory)
{
  std::cout << "Module loaded: Pseudoprogression Estimation on Existing Model:" << std::endl;
  std::vector<double> finalresult;
  std::vector<std::map<CAPTK::ImageModalityType, std::string>> QualifiedSubjects = LoadQualifiedSubjectsFromGivenDirectoryForPseudoProgression(CAPTK::MachineLearningApplicationSubtype::TESTING, inputdirectory, true, true, true, true);
  std::cout << "Number of subjects with required input: " << QualifiedSubjects.size() << std::endl;
  PseudoProgressionEstimator objPseudoProgressionEstimator;
  bool data = objPseudoProgressionEstimator.PseudoProgressionEstimateOnExistingModel(QualifiedSubjects, modeldirectory, inputdirectory, outputdirectory, true, true, true, true);
  //for (unsigned int subjectID = 0; subjectID < QualifiedSubjects.size(); subjectID++)
  //{
  //	std::map<ImageModalityType, std::string> onesubject = QualifiedSubjects[subjectID];
  //	//std::cout << static_cast<std::string>(onesubject[IMAGE_TYPE_SUDOID]) << ": " << result[subjectID] << std::endl;
  //}
  return EXIT_SUCCESS;
}
int PrepareNewSurvivalPredictionModel(const std::string inputdirectory, const std::string outputdirectory)
{
  std::cout << "Module loaded: Prepare Pseudoprogression Prediction Model." << std::endl;
  std::vector<double> finalresult;
  std::vector<std::map<CAPTK::ImageModalityType, std::string>> QualifiedSubjects = LoadQualifiedSubjectsFromGivenDirectoryForPseudoProgression(CAPTK::MachineLearningApplicationSubtype::TRAINING, inputdirectory, true, true, true, true);
  PseudoProgressionEstimator objPseudoProgressionEstimator;
  std::cout << "Number of subjects with required input: " << QualifiedSubjects.size() << std::endl;
  if (QualifiedSubjects.size() == 0)
    std::cout << "No subject found with required input. Exiting...." << std::endl;
  else if (QualifiedSubjects.size() >0 && QualifiedSubjects.size() <= 20)
    std::cout << "There should be atleast 20 patients to build reliable pseudo-progression model. Exiting...." << std::endl;
  else
    objPseudoProgressionEstimator.TrainNewModelOnGivenData(QualifiedSubjects, outputdirectory, true, true, true, true);
  return EXIT_SUCCESS;
}
int main(int argc, char **argv)
{
  cbica::CmdParser parser = cbica::CmdParser(argc, argv, "PseudoProgressionEstimator");
  parser.addRequiredParameter("t", "type", cbica::Parameter::STRING, "", "The option of preparing a new model (=0), and for testing on an existing model (=1)");
  parser.addRequiredParameter("i", "input", cbica::Parameter::STRING, "", "The input directory having test subjects");
  parser.addOptionalParameter("m", "model", cbica::Parameter::STRING, "", "The directory having SVM models");
  parser.addRequiredParameter("o", "output", cbica::Parameter::STRING, "", "The output direcory to write output");
  parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");

  // parameters to get from the command line
  cbica::Logging logger;
  std::string loggerFile;
  bool loggerRequested = false;



  int tempPosition;
  std::string inputDirectoryName, modelDirectoryName, outputDirectoryName, toWrite;
  int applicationType;

  applicationType = 0;

  if (parser.compareParameter("L", tempPosition))
  {
    loggerFile = argv[tempPosition + 1];
    loggerRequested = true;
    logger.UseNewFile(loggerFile);
  }
  if (parser.compareParameter("i", tempPosition))
  {
    inputDirectoryName = argv[tempPosition + 1];
    //inputDirectoryName = "E:/SoftwareDevelopmentProjects/PseudoprogressionRelatedMaterial/TrainingData";
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
  if (applicationType == CAPTK::MachineLearningApplicationSubtype::TESTING && modelDirectoryName == "")
  {
    std::cout << "Please specify a directory having model file." << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Input directory name:" << inputDirectoryName << std::endl;
  std::cout << "Output directory name:" << outputDirectoryName << std::endl;


  if (!cbica::directoryExists(inputDirectoryName))
  {
    std::cout << "The input directory does not exist:" << inputDirectoryName << std::endl;
    return EXIT_FAILURE;
  }
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
    SurvivalPredictionOnExistingModel(modelDirectoryName, inputDirectoryName, outputDirectoryName);
  }
  else if (applicationType == CAPTK::MachineLearningApplicationSubtype::TRAINING)
    PrepareNewSurvivalPredictionModel(inputDirectoryName, outputDirectoryName);
  else
  {
    parser.echoVersion();
    return EXIT_SUCCESS;
  }

  std::cout << "Finished successfully\n";

  return EXIT_SUCCESS;
}