#include "ImagingSubtypePredictor.h"
#include "cbicaUtilities.h"
#include "cbicaCmdParser.h"

//------------------Survival Prediction on existing model-----------------------
std::vector<std::map<ImageModalityType, std::string>> LoadQualifiedSubjectsFromGivenDirectory(const std::string directoryname)
{
  std::map<ImageModalityType, std::string> OneQualifiedSubject;
  std::vector<std::map<ImageModalityType, std::string>> QualifiedSubjects;
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
          if ((filePath_lower.find("segmentation") != std::string::npos)
            && (extension == HDR_EXT || extension
            == NII_EXT || extension == NII_GZ_EXT))
            labelPath = subjectPath + "/SEGMENTATION" + "/" + files[i];
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

    if (labelPath.empty() || t1FilePath.empty() || t2FilePath.empty() || t1ceFilePath.empty() || t2FlairFilePath.empty() || rcbvFilePath.empty() || axFilePath.empty() || faFilePath
      == "" || radFilePath.empty() || trFilePath.empty() || psrFilePath.empty() || phFilePath.empty())
      continue;

    OneQualifiedSubject[IMAGE_TYPE_T1] = t1FilePath;
    OneQualifiedSubject[IMAGE_TYPE_T2] = t2FilePath;
    OneQualifiedSubject[IMAGE_TYPE_T1CE] = t1ceFilePath;
    OneQualifiedSubject[IMAGE_TYPE_T2FLAIR] = t2FlairFilePath;
    OneQualifiedSubject[IMAGE_TYPE_AX] = axFilePath;
    OneQualifiedSubject[IMAGE_TYPE_FA] = faFilePath;
    OneQualifiedSubject[IMAGE_TYPE_RAD] = radFilePath;
    OneQualifiedSubject[IMAGE_TYPE_TR] = trFilePath;

    OneQualifiedSubject[IMAGE_TYPE_PSR] = psrFilePath;
    OneQualifiedSubject[IMAGE_TYPE_PH] = phFilePath;
    OneQualifiedSubject[IMAGE_TYPE_RCBV] = rcbvFilePath;

    OneQualifiedSubject[IMAGE_TYPE_SEG] = labelPath;
    OneQualifiedSubject[IMAGE_TYPE_SUDOID] = subjectNames[sid];

    QualifiedSubjects.push_back(OneQualifiedSubject);
  }
  return QualifiedSubjects;
}


int SubtypePrediction(const std::string modelfile,
  const std::string inputdirectory,
  const std::string outputdirectory)
{
  std::cout << "Module loaded: Subtype Prediction:" << std::endl;
  std::vector<double> finalresult;
  std::vector<std::map<ImageModalityType, std::string>> QualifiedSubjects = LoadQualifiedSubjectsFromGivenDirectory(inputdirectory);
  std::cout << "Number of subjects with required input: " << QualifiedSubjects.size() << std::endl;
  ImagingSubtypePredictor objImagingSubtypePredictor;
  for (unsigned int subjectID = 0; subjectID < QualifiedSubjects.size(); subjectID++)
  {
    std::map<ImageModalityType, std::string> onesubject = QualifiedSubjects[subjectID];

    std::vector<std::map<ImageModalityType, std::string>> OneQualifiedSubject;
    OneQualifiedSubject.push_back(QualifiedSubjects[subjectID]);
    //VectorDouble result = objImagingSubtypePredictor.SurvivalPredictionOnExistingModel(modelfile, inputdirectory, QualifiedSubjects, outputdirectory);
    //std::cout << static_cast<std::string>(onesubject[IMAGE_TYPE_SUDOID])<<" " <<result[0]<<std::endl;
  }
  return EXIT_SUCCESS;
}
//----------------------------------------------
int main(int argc, char **argv)
{
  cbica::CmdParser parser = cbica::CmdParser(argc, argv, "ImagingSubtypePredictor");
  parser.addRequiredParameter("i", "input", cbica::Parameter::STRING, "", "The input directory having test subjects");
  parser.addOptionalParameter("m", "model", cbica::Parameter::STRING, "", "The directory having model");
  parser.addRequiredParameter("o", "output", cbica::Parameter::STRING, "", "The output direcory to write output");
  parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");
  parser.exampleUsage("ImagingSubtypePredictor -i <input dir> -o <output dir>");

  // parameters to get from the command line
  cbica::Logging logger;
  std::string loggerFile;
  bool loggerRequested = false;

  int tempPosition;
  std::string inputDirectoryName, modelFileName, outputDirectoryName, toWrite;

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
    modelFileName = argv[tempPosition + 1];
  }

  if (parser.compareParameter("o", tempPosition))
  {
    outputDirectoryName = argv[tempPosition + 1];
  }
  std::cout << "Input Directory:" << inputDirectoryName << std::endl;
  std::cout << "Output Directory:" << outputDirectoryName << std::endl;
  std::cout << "Model File:" << outputDirectoryName << std::endl;


  if (!cbica::directoryExists(inputDirectoryName))
  {
    std::cout << "The input directory does not exist:" << inputDirectoryName << std::endl;
    return EXIT_FAILURE;
  }
  if (!cbica::directoryExists(outputDirectoryName))
  {
    std::cout << "The output directory does not exist:" << outputDirectoryName << std::endl;
    return EXIT_FAILURE;
  }
  if (!cbica::fileExists(modelFileName))
  {
    std::cout << "The model file does not exist:" << modelFileName << std::endl;
    return EXIT_FAILURE;
  }
  SubtypePrediction(modelFileName, inputDirectoryName, outputDirectoryName);

  std::cout << "Finished successfully\n";

  return EXIT_SUCCESS;
}