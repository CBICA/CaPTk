#include "PopulationAtlases.h"
#include "cbicaUtilities.h"
#include "cbicaCmdParser.h"
#include "CaPTkEnums.h"

//------------------Survival Prediction on existing model-----------------------
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
          else if ((filePath_lower.find("segmentation") != std::string::npos)
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
        std::string filePath = subjectPath + "/CONVENTIONAL" + "/" + files[i];
        std::string extension = cbica::getFilenameExtension(filePath, false);

        if ((files[i].find("t1ce") != std::string::npos || files[i].find("T1CE") != std::string::npos || files[i].find("T1ce") != std::string::npos || files[i].find("T1-gd") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          t1ceFilePath = subjectPath + "/CONVENTIONAL" + "/" + files[i];
        else if ((files[i].find("t1") != std::string::npos || files[i].find("T1") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          t1FilePath = subjectPath + "/CONVENTIONAL" + "/" + files[i];
        else if ((files[i].find("t2") != std::string::npos || files[i].find("T2") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          t2FilePath = subjectPath + "/CONVENTIONAL" + "/" + files[i];
        else if ((files[i].find("flair") != std::string::npos || files[i].find("FLAIR") != std::string::npos || files[i].find("Flair") != std::string::npos || files[i].find("T2-Flair") != std::string::npos) && (extension == HDR_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
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
      featureFilePath = subjectPath + "/features.csv";

    if (labelPath.empty() || t1FilePath.empty() || t2FilePath.empty() || t1ceFilePath.empty() || t2FlairFilePath.empty() || rcbvFilePath.empty() || axFilePath.empty() || faFilePath.empty()
      || radFilePath.empty() || trFilePath.empty() || psrFilePath.empty() || phFilePath.empty() || featureFilePath.empty())
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
int SurvivalPredictionOnExistingModel(const std::string modeldirectory,
  const std::string inputdirectory,
  const std::string outputdirectory)
{
  std::cout << "Module loaded: Survival Prediction on Existing Model:" << std::endl;
  std::vector<double> finalresult;
  std::vector<std::map<CAPTK::ImageModalityType, std::string>> QualifiedSubjects = LoadQualifiedSubjectsFromGivenDirectory(inputdirectory);
  std::cout << "Number of subjects with required input: " << QualifiedSubjects.size() << std::endl;
  /*PopulationAtlases objPopulationAtlases;
  for (unsigned int subjectID = 0; subjectID < QualifiedSubjects.size(); subjectID++)
  {
  std::map<ImageModalityType, std::string> onesubject = QualifiedSubjects[subjectID];

  std::vector<std::map<ImageModalityType, std::string>> OneQualifiedSubject;
  OneQualifiedSubject.push_back(QualifiedSubjects[subjectID]);
  VectorDouble result = objPopulationAtlases.SurvivalPredictionOnExistingModel(modeldirectory, inputdirectory, QualifiedSubjects, outputdirectory);
  std::cout << static_cast<std::string>(onesubject[IMAGE_TYPE_SUDOID])<<" " <<result[0]<<std::endl;
  }*/
  return EXIT_SUCCESS;
}
int PrepareNewSurvivalPredictionModel(const std::string inputdirectory,
  const std::string outputdirectory)
{
  //std::cout << "Module loaded: Prepare Survival Prediction Model." << std::endl;
  //std::vector<double> finalresult;
  //std::vector<std::map<ImageModalityType, std::string>> QualifiedSubjects = LoadQualifiedSubjectsFromGivenDirectory(inputdirectory);
  //PopulationAtlases objPopulationAtlases;
  //std::cout << "Number of subjects with required input: " << QualifiedSubjects.size() << std::endl;
  //objPopulationAtlases.PrepareNewSurvivalPredictionModel(inputdirectory, QualifiedSubjects, outputdirectory);
  return EXIT_SUCCESS;
}
int main(int argc, char **argv)
{
  cbica::CmdParser parser = cbica::CmdParser(argc, argv, "PopulationAtlases");
  parser.addRequiredParameter("i", "input", cbica::Parameter::STRING, "", "The input directory.");
  parser.addRequiredParameter("l", "label", cbica::Parameter::STRING, "", "The input label file in .csv format.");
  parser.addRequiredParameter("a", "atlas", cbica::Parameter::STRING, "", "The atlas template.");
  parser.addRequiredParameter("o", "output", cbica::Parameter::STRING, "", "The output directory.");
  parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");

  parser.writeCWLFile(".", parser.getExeName(), false);
  parser.exampleUsage("PopulationAtlases -i <input dir> -l atlaslabelfile.csv -a jakob_stripped_with_cere_lps_256256128.nii.gz -o <output dir>");

  // parameters to get from the command line
  cbica::Logging logger;
  std::string loggerFile;
  bool loggerRequested = false;

  int tempPosition;
  std::string inputDirectoryName, inputLabelName, inputAtlasName, outputDirectoryName, toWrite;
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

  if (parser.compareParameter("a", tempPosition))
  {
    inputAtlasName = argv[tempPosition + 1];
  }

  if (parser.compareParameter("o", tempPosition))
  {
    outputDirectoryName = argv[tempPosition + 1];
  }
  if (parser.compareParameter("l", tempPosition))
  {
    inputLabelName = argv[tempPosition + 1];
  }
  std::cout << "Input Directory:" << inputDirectoryName << std::endl;
  std::cout << "Input Atlas File:" << inputAtlasName << std::endl;
  std::cout << "Input Label File:" << inputLabelName << std::endl;
  std::cout << "Output Directory:" << outputDirectoryName << std::endl;

  if (!cbica::directoryExists(inputDirectoryName))
  {
    std::cout << "The input directory does not exist:" << inputDirectoryName << std::endl;
    return EXIT_FAILURE;
  }
  if (!cbica::directoryExists(outputDirectoryName))
  {
    if (!cbica::createDirectory(outputDirectoryName))
      std::cout << "The output directory can not be created:" << outputDirectoryName << std::endl;
    return EXIT_FAILURE;
  }
  if (!cbica::fileExists(inputLabelName))
  {
    std::cout << "The input label file does not exist:" << inputLabelName << std::endl;
    return EXIT_FAILURE;
  }
  if (!cbica::fileExists(inputAtlasName))
  {
    std::cout << "The input atlas file does not exist:" << inputAtlasName << std::endl;
    return EXIT_FAILURE;
  }
  PopulationAtlases objPopulationAtlases;
  std::vector<typename ImageType::Pointer> atlases = objPopulationAtlases.GeneratePopualtionAtlas(inputDirectoryName, inputLabelName, inputAtlasName, outputDirectoryName);
  std::cout << "Writing atlases in the specified output directory.\n";
  for (int i = 0; i < atlases.size(); i++)
  {
    typedef itk::ImageFileWriter< ImageType > WriterType;
    WriterType::Pointer writer1 = WriterType::New();
    writer1->SetFileName(outputDirectoryName + "/Atlas_" + std::to_string(i + 1) + ".nii.gz");
    writer1->SetInput(atlases[i]);
    writer1->Update();
  }
  std::cout << "Finished successfully.\n";

  return EXIT_SUCCESS;
}