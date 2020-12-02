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
  //std::cout << "This functionality has been removed from this CaPTk release, and we are actively working on an optimized robust implementation that should enable generalization in multi-institutional data. We expect this to be released in our next patch release, in Q4 2020.\n";
  //return EXIT_FAILURE;
  cbica::CmdParser parser = cbica::CmdParser(argc, argv, "PerfusionPCA");
  parser.addRequiredParameter("i", "input", cbica::Parameter::STRING, "", "The input directory.");
  parser.addRequiredParameter("t", "type", cbica::Parameter::INTEGER, "", "The option of preparing a new model (=0), and for testing on an existing model (=1)");
  parser.addOptionalParameter("n", "number of PCAs", cbica::Parameter::FLOAT, "", "The number of principal components.");
  parser.addOptionalParameter("m", "model", cbica::Parameter::STRING, "", "The directory having PCA models");
  parser.addOptionalParameter("vt", "variance threshold", cbica::Parameter::FLOAT, "", "The variance threshold");
  parser.addOptionalParameter("dif", "dump intermediate files", cbica::Parameter::BOOLEAN, "", "Write intermediate file containing perfusion data for whole population.");
  parser.addRequiredParameter("o", "output", cbica::Parameter::STRING, "", "The output directory.");
  parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");
  parser.addExampleUsage("-t 0 -i C:/properly/formatted/inputDir -o C:/outputDir -n 5", "Trains a new model based on the samples in inputDir");
  parser.addExampleUsage("-t 1 -i C:/input -m C:/model -o C:/output -m C:/modelDir -n 5", "Tests an existing model for inputs in 'C:/input' based on 'C:/modelDir' ");
  parser.addApplicationDescription("Perfusion PCA Training and Prediction application");

  // parameters to get from the command line
  cbica::Logging logger;
  std::string loggerFile;
  bool loggerRequested = false;

  int applicationType;
  applicationType = 0;

  float inputPCs = 0;
  float varianceThreshold = 0.0;
  std::string inputFileName, inputMaskName, outputDirectoryName, modelDirectoryName;
  bool nPCsDefined = false;
  bool varianceThresholdDefined = false;
  bool dif = false;

  parser.getParameterValue("i", inputFileName);
  parser.getParameterValue("o", outputDirectoryName);
  parser.getParameterValue("t", applicationType);
  
  if(parser.isPresent("dif"))
  {
	  parser.getParameterValue("dif", dif);
  }
  else
	  dif = false;

  std::cout << " dif: " << dif << std::endl;

  if (parser.isPresent("n"))
  {
	  nPCsDefined = true;
	  parser.getParameterValue("n", inputPCs);
  }

  if (parser.isPresent("vt"))
  {
	  varianceThresholdDefined = true;
	  parser.getParameterValue("vt", varianceThreshold);
  }

  //if user provides both n and vt, then dont run app
  //put a msg to provide any one
  if (nPCsDefined && varianceThresholdDefined)
  {
	  std::cout << "Please provide either the number of principal components or the variance threshold and re-run the application." << std::endl;
	  return EXIT_FAILURE;
  }


  if (parser.isPresent("m"))
  {
    parser.getParameterValue("m", modelDirectoryName);
  }
  if (parser.isPresent("L"))
  {
    parser.getParameterValue("L", loggerFile);
    loggerRequested = true;
    logger.UseNewFile(loggerFile);
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
  std::string inValidSubject;
  if (object_pca.LoadData(QualifiedSubjects, inValidSubject) == PerfusionPCA::DifferentTimePoints)
  {
	  std::cout << "Could not load data. Please check that all input data has the same number of time points." << std::endl;
	  return EXIT_FAILURE;
  }

  if (applicationType == CAPTK::MachineLearningApplicationSubtype::TESTING)
  {
	  //TBD: during testing read the number of PCs from file 
    std::cout << "Model directory name:" << modelDirectoryName << std::endl;
    if (!cbica::directoryExists(modelDirectoryName))
    {
      std::cout << "The model directory does not exist:" << modelDirectoryName << std::endl;
      return EXIT_FAILURE;
    }
    if (!cbica::fileExists(modelDirectoryName +"/PCA_PERF.csv") || 
		!cbica::fileExists(modelDirectoryName + "/Mean_PERF.csv"))
    {
      std::cout << "The required model files PCA_PERF.csv, Mean_PERF.csv do not exist in the model directory:" <<  modelDirectoryName << std::endl;
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

	//set number of PCs on the object
	if (nPCsDefined)
		object_pca.SetNumberOfPCs(inputPCs);
	else if (varianceThresholdDefined)
		object_pca.SetVarianceThreshold(varianceThreshold);
    PerfusionPCA::ErrorCode code = object_pca.ApplyExistingPCAModel(inputPCs, inputFileName, outputDirectoryName, QualifiedSubjects,modelDirectoryName);
	if (code == PerfusionPCA::ErrorCode::DifferentTimePoints)
	{
		std::cout << "Could not load data. Please check that all input data has the same number of time points." << std::endl;
		return EXIT_FAILURE;
	}
	else if (code == PerfusionPCA::ErrorCode::NoError)
	{
		std::cout << "principal components have been saved at the specified locations.\n";
		std::cout << "Finished successfully.\n";
		return EXIT_SUCCESS;
	}
  }
  else if (applicationType == CAPTK::MachineLearningApplicationSubtype::TRAINING)
  {
		//if user doesnt provide n or vt, give a msg stating that we will provide only PCA parameters 
		// and not the PCA images and that the run will take quite some time and request for confirmation
		// press y to continue or n to exit and provide n or vt.
	  if (!nPCsDefined && !varianceThresholdDefined)
	  {
		  std::cout << "You did not provide the number of principal components or the variance threshold. \
We will provide only the PCA parameters and not any PCA images. The application will take a long time to run. \
Do you want to continue? Press 'y' to contine or 'n' (and press Enter) to exit and provide either of the parameters.";
		  char ch = cin.get();
		  if (ch == 'y')
		  {
			  nPCsDefined = true;
			  inputPCs = 0; //we don't provide any PCA images
		  }
		  else
		  {
			  std::cout << "Please provide either the number of principal components or the variance threshold and re-run the application." << std::endl;
			  return EXIT_FAILURE;
		  }
	  }

	  //spit a txt file : NumberofPCs having only a single number of PCs
	  if (nPCsDefined)
		  object_pca.SetNumberOfPCs(inputPCs);
	  else if (varianceThresholdDefined)
		  object_pca.SetVarianceThreshold(varianceThreshold);
	  //we won't have a situation here that both nPCsDefined and varianceThresholdDefined
	  //are defined. This is handled upstream.

	  object_pca.RequestPerfusionDataWholePopulation(dif);//dump intermediate files or not

	  if (object_pca.TrainNewPerfusionModel(inputPCs, inputFileName, outputDirectoryName, QualifiedSubjects))
	  {
		  std::cout << "principal components have been saved at the specified locations.\n";
		  std::cout << "Finished successfully.\n";
		  return EXIT_SUCCESS;
	  }
	  else
	  {
		  std::cout << "Something went wrong during the calculation. Please contact software@cbica.upenn.edu \n";
		  return EXIT_FAILURE;
	  }
  }
  else
  {
    parser.echoVersion();
    return EXIT_SUCCESS;
  }

}
