#include "MolecularSubtypePredictor.h"
#include "cbicaUtilities.h"
#include "cbicaCmdParser.h"
#include "CaPTkGUIUtils.h"

//------------------Molecular subtype prediction on existing model-----------------------
std::vector<std::map<CAPTK::ImageModalityType, std::string>> LoadQualifiedSubjectsFromGivenDirectory(const std::string directoryname)
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
    std::string rcbvFilePath = "";
    std::string psrFilePath = "";
    std::string phFilePath = "";
    std::string labelPath = "";
    std::string atlasPath = "";
    std::string parametersPath = "";
    std::string featureFilePath = "";

    std::vector<std::string> files;

    if (cbica::directoryExists(subjectPath + "/SEGMENTATION"))
    {
      files = cbica::filesInDirectory(subjectPath + "/SEGMENTATION", false);
      if (files.size() == 1)
      {
        labelPath = subjectPath + "/SEGMENTATION" + "/" + files[0];
      }
      else
      {
        for (unsigned int i = 0; i < files.size(); i++)
        {
          std::string filePath = subjectPath + "/SEGMENTATION" + "/" + files[i], filePath_lower;
          std::string extension = cbica::getFilenameExtension(filePath, false);
          filePath_lower = filePath;
          std::transform(filePath_lower.begin(), filePath_lower.end(), filePath_lower.begin(), ::tolower);
          if ((filePath_lower.find("atlas") != std::string::npos || filePath_lower.find("jakob_label") != std::string::npos)
            && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
            atlasPath = subjectPath + "/SEGMENTATION" + "/" + files[i];
          else if ((filePath_lower.find("segmentation") != std::string::npos || filePath_lower.find("label") != std::string::npos)
            && (extension == HDR_EXT || extension
            == NII_EXT || extension == NII_GZ_EXT))
            labelPath = subjectPath + "/SEGMENTATION" + "/" + files[i];
          else if ((filePath_lower.find("parameter") != std::string::npos)
            && (extension == PARAM_EXT))
            parametersPath = subjectPath + "/SEGMENTATION" + "/" + files[i];
        }
      }
    }

    if (cbica::directoryExists(subjectPath + "/CONVENTIONAL"))
    {
      files = cbica::filesInDirectory(subjectPath + "/CONVENTIONAL", false);
      for (unsigned int i = 0; i < files.size(); i++)
      {
        std::string filePath = subjectPath + "/CONVENTIONAL" + "/" + files[i], filePath_lower;
        std::string extension = cbica::getFilenameExtension(filePath, false);
        filePath_lower = filePath;
        std::transform(filePath_lower.begin(), filePath_lower.end(), filePath_lower.begin(), ::tolower);

        if ((filePath_lower.find("t1ce") != std::string::npos || filePath_lower.find("t1-gd") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          t1ceFilePath = subjectPath + "/CONVENTIONAL" + "/" + files[i];
        else if ((filePath_lower.find("t1") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          t1FilePath = subjectPath + "/CONVENTIONAL" + "/" + files[i];
        else if ((filePath_lower.find("t2") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          t2FilePath = subjectPath + "/CONVENTIONAL" + "/" + files[i];
        else if ((filePath_lower.find("flair") != std::string::npos || filePath_lower.find("t2-flair") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          t2FlairFilePath = subjectPath + "/CONVENTIONAL" + "/" + files[i];
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
      }
    }
    if (cbica::directoryExists(subjectPath + "/DTI"))
    {
      files = cbica::filesInDirectory(subjectPath + "/DTI", false);
      for (unsigned int i = 0; i < files.size(); i++)
      {
        std::string filePath = subjectPath + "/DTI/" + files[i], filePath_lower;
        std::string extension = cbica::getFilenameExtension(filePath, false);
        filePath_lower = filePath;
        std::transform(filePath_lower.begin(), filePath_lower.end(), filePath_lower.begin(), ::tolower);

        if ((filePath_lower.find("axial") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          axFilePath = subjectPath + "/DTI/" + files[i];
        else if ((filePath_lower.find("fractional") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          faFilePath = subjectPath + "/DTI/" + files[i];
        else if ((filePath_lower.find("radial") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          radFilePath = subjectPath + "/DTI/" + files[i];
        else if ((filePath_lower.find("trace") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          trFilePath = subjectPath + "/DTI/" + files[i];
      }
    }
    if (cbica::fileExists(subjectPath + "/features.csv"))
      featureFilePath = subjectPath + "/features.csv";

    if (labelPath.empty() || t1FilePath.empty() || t2FilePath.empty() || t1ceFilePath.empty() || t2FlairFilePath.empty() || rcbvFilePath.empty() || axFilePath.empty() || faFilePath
      == "" || radFilePath.empty() || trFilePath.empty() || psrFilePath.empty() || phFilePath.empty() || featureFilePath.empty())
      continue;

    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1] = t1FilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2] = t2FilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1CE] = t1ceFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR] = t2FlairFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_AX] = axFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_FA] = faFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_RAD] = radFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_TR] = trFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_PSR] = psrFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_PH] = phFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_RCBV] = rcbvFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG] = labelPath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_ATLAS] = atlasPath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_PARAMS] = parametersPath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_FEATURES] = featureFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID] = subjectNames[sid];

    QualifiedSubjects.push_back(OneQualifiedSubject);
  }
  return QualifiedSubjects;
}

int MolecularSubtypePredictionOnExistingModel(const std::string modeldirectory,
  const std::string inputdirectory,
  const std::string outputdirectory)
{
  std::cout << "Module loaded: Moelcular Subtype Prediction on Existing Model:" << std::endl;
  std::vector<double> finalresult;
  std::vector<std::map<CAPTK::ImageModalityType, std::string>> QualifiedSubjects = LoadQualifiedSubjectsFromGivenDirectory(inputdirectory);
  std::cout << "Number of subjects with required input: " << QualifiedSubjects.size() << std::endl;

  MolecularSubtypePredictor objMolecularSubtypePredictor;

  VectorDouble result = objMolecularSubtypePredictor.MolecularSubtypePredictionOnExistingModel(modeldirectory, inputdirectory, QualifiedSubjects, outputdirectory);

  for (unsigned int subjectID = 0; subjectID < QualifiedSubjects.size(); subjectID++)
  {
    std::map<CAPTK::ImageModalityType, std::string> onesubject = QualifiedSubjects[subjectID];
    if(result[subjectID]==1)
      std::cout << static_cast<std::string>(onesubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]) << ": Proneural" << std::endl;
    else if (result[subjectID] == 2)
      std::cout << static_cast<std::string>(onesubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]) << ": Neural" << std::endl;
    else if (result[subjectID] == 3)
      std::cout << static_cast<std::string>(onesubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]) << ": Mesenchymal" << std::endl;
    else
      std::cout << static_cast<std::string>(onesubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]) << ": Classical" << std::endl;

  }
  return EXIT_SUCCESS;
}


#ifdef BUILD_MSUBTYPE
int PrepareNewMolecularPredictionModel(const std::string inputdirectory,
  const std::string outputdirectory)
{
  std::cout << "Module loaded: Prepare Molecular Subtype Prediction Model." << std::endl;
  std::vector<double> finalresult;
  std::vector<std::map<CAPTK::ImageModalityType, std::string>> QualifiedSubjects = LoadQualifiedSubjectsFromGivenDirectory(inputdirectory);
  MolecularSubtypePredictor objMolecularSubtypePredictor;
  std::cout << "Number of subjects with required input: " << QualifiedSubjects.size() << std::endl;
  objMolecularSubtypePredictor.PrepareNewMolecularPredictionModel(inputdirectory, QualifiedSubjects, outputdirectory);
  return EXIT_SUCCESS;
}
#endif

int main(int argc, char **argv)
{
  std::cout << "This functionality has been removed from this CaPTk release, \
and we are actively testing an optimized robust implementation.\n";
  return EXIT_FAILURE;
#ifdef BUILD_MSUBTYPE
  cbica::CmdParser parser = cbica::CmdParser(argc, argv, "MolecularSubtypePredictor");
  parser.addRequiredParameter("t", "type", cbica::Parameter::STRING, "", "The option of preparing a new model (=0), and for testing on an existing model (=1)");
  parser.addRequiredParameter("i", "input", cbica::Parameter::STRING, "", "The input directory having test subjects");
  parser.addOptionalParameter("m", "model", cbica::Parameter::STRING, "", "The directory having SVM models", "Penn Model: " + getAppropriateDownloadLink("MolecularSubtypePredictor", "Model"));
  parser.addRequiredParameter("o", "output", cbica::Parameter::STRING, "", "The output direcory to write output");
  parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");
  //parser.exampleUsage("MolecularSubtypePredictor -t 0 -i <input dir> -o <output dir>");
  parser.addExampleUsage("-t 0 -i C:/properly/formatted/inputDir -o C:/outputDir", "Trains a new model based on the samples in inputDir");
  parser.addExampleUsage("-t 1 -i C:/input -m C:/model -o C:/output", "Tests an existing model for inputs in 'C:/input' based on 'C:/model' ");
  parser.addApplicationDescription("Molecular Subtype Training and Prediction application");

  // parameters to get from the command line
  cbica::Logging logger;
  std::string loggerFile;
  bool loggerRequested = false;

  int tempPosition;
  std::string inputDirectoryName, modelDirectoryName, outputDirectoryName, toWrite;
  int applicationType;

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

  std::cout << "Input directory name:" << inputDirectoryName << std::endl;
  std::cout << "Output directory name:" << outputDirectoryName << std::endl;

  if (applicationType == CAPTK::MachineLearningApplicationSubtype::TESTING)
  {
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
    else
      std::cout << "Model directory name:" << modelDirectoryName << std::endl;
    MolecularSubtypePredictionOnExistingModel(modelDirectoryName, inputDirectoryName, outputDirectoryName);
  }
  else if (applicationType == CAPTK::MachineLearningApplicationSubtype::TRAINING)
    PrepareNewMolecularPredictionModel(inputDirectoryName, outputDirectoryName);
  else
  {
    parser.echoVersion();
    return EXIT_SUCCESS;
  }

  std::cout << "Finished successfully\n";

  return EXIT_SUCCESS;
#else
  std::cout << "This application has been disabled.\n";
  return EXIT_SUCCESS;
#endif
}