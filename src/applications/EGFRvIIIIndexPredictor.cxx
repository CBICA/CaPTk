#include "EGFRvIIIIndexPredictor.h"
#include "cbicaUtilities.h"
#include "cbicaCmdParser.h"
#include "CaPTkEnums.h"

//------------------EGFRvIII Prediction on existing model-----------------------
std::vector<std::map<CAPTK::ImageModalityType, std::string>> 
LoadQualifiedSubjectsFromGivenDirectory(const std::string directoryname, 
  int applicationtype)
{
	std::map<CAPTK::ImageModalityType, std::string> OneQualifiedSubject;
	std::vector<std::map<CAPTK::ImageModalityType, std::string>> QualifiedSubjects;
	std::vector<std::string> subjectNames = cbica::subdirectoriesInDirectory(directoryname);
  subjectNames.push_back("AASF_2014.02.18");
  subjectNames.push_back("AASG_2014.02.26");
  subjectNames.push_back("AASH_2014.02.27");
  subjectNames.push_back("AASK_2013.08.26");
  subjectNames.push_back("AASN_2014.04.22");
  subjectNames.push_back("AASO_2014.05.11");
  subjectNames.push_back("AASR_2014.05.19");
  subjectNames.push_back("AASU_2014.06.03");
  subjectNames.push_back("AASV_2013.08.11");
  subjectNames.push_back("AASW_2014.06.05");
  subjectNames.push_back("AASX_2014.06.18");
  subjectNames.push_back("AASY_2014.06.25");
  subjectNames.push_back("AASZ_2006.07.19");
  subjectNames.push_back("AATB_2006.08.15");
  subjectNames.push_back("AATD_2006.09.24");
  subjectNames.push_back("AATE_2006.11.26");
  subjectNames.push_back("AATF_2006.12.19");
  subjectNames.push_back("AATP_2008.11.01");
  subjectNames.push_back("AATV_2009.01.15");
  subjectNames.push_back("AATW_2006.07.23");
  subjectNames.push_back("AATX_2006.09.07");
  subjectNames.push_back("AAUC_2008.02.17");
  subjectNames.push_back("AAUE_2006.07.11");
  subjectNames.push_back("AAUF_2008.01.07");
  subjectNames.push_back("AAUM_2008.03.06");
  subjectNames.push_back("AAUQ_2007.05.13");
  subjectNames.push_back("AAUR_2007.01.28");
  subjectNames.push_back("AAVG_2006.08.30");
  subjectNames.push_back("AAVJ_2006.12.12");
  subjectNames.push_back("AAVS_2008.01.23");
  subjectNames.push_back("AAVV_2007.05.24");
  subjectNames.push_back("AAWA_2008.10.01");
  subjectNames.push_back("AAWB_2008.12.18");
  subjectNames.push_back("AAWG_2007.04.05");
  subjectNames.push_back("AAWH_2007.01.02");
  subjectNames.push_back("AAWI_2006.11.28");
  subjectNames.push_back("AAWJ_2006.07.27");
  subjectNames.push_back("AAWK_2007.04.11");
  subjectNames.push_back("AAXI_2007.01.14");
  subjectNames.push_back("AAXJ_2007.01.25");
  subjectNames.push_back("AAXL_2006.09.01");
  subjectNames.push_back("AAXM_2007.12.20");
  subjectNames.push_back("AAXO_2007.12.13");
  subjectNames.push_back("AAXQ_2007.08.30");
  subjectNames.push_back("AAXV_2007.06.28");
  subjectNames.push_back("AAXW_2007.04.03");
  subjectNames.push_back("AAYI_2008.11.26");
  subjectNames.push_back("AAYU_2008.03.02");
  subjectNames.push_back("AAYV_2009.04.30");
  subjectNames.push_back("AAYW_2008.10.09");
  subjectNames.push_back("AAZA_2009.01.02");
  subjectNames.push_back("AAZD_2008.06.22");
  subjectNames.push_back("AAZF_2008.07.17");
  subjectNames.push_back("AAZH_2008.06.13");

  std::sort(subjectNames.begin(), subjectNames.end());

	for (unsigned int sid = 0; sid < subjectNames.size(); sid++)
	{
		std::string subjectPath = directoryname + "/" + subjectNames[sid];

		std::string t1ceFilePath    = "";
		std::string t1FilePath      = "";
		std::string t2FilePath      = "";
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
			files = cbica::filesInDirectory(subjectPath + "/SEGMENTATION",false);
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
					else if ((filePath_lower.find("label") != std::string::npos)
						&& (extension == HDR_EXT || extension
						== NII_EXT || extension == NII_GZ_EXT))
						labelPath = subjectPath + "/SEGMENTATION" + "/" + files[i];
				}
			}
		}

		if (cbica::directoryExists(subjectPath + "/CONVENTIONAL"))
		{
			files = cbica::filesInDirectory(subjectPath + "/CONVENTIONAL",false);
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
			files = cbica::filesInDirectory(subjectPath + "/PERFUSION",false);
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
			files = cbica::filesInDirectory(subjectPath + "/DTI",false);
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
        || radFilePath.empty() || trFilePath.empty() || psrFilePath.empty() || phFilePath.empty())
			continue;

    if (applicationtype == CAPTK::MachineLearningApplicationSubtype::TRAINING && featureFilePath.empty())
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
		OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_FEATURES] = featureFilePath;
		OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID] = subjectNames[sid];

		QualifiedSubjects.push_back(OneQualifiedSubject);
	}
	return QualifiedSubjects;
}

int EGFRvIIIPredictionOnExistingModel(const std::string modeldirectory,
	const std::string inputdirectory,
	const std::string outputdirectory)
{
	std::cout << "Module loaded: EGFRvIII Prediction on Existing Model:" << std::endl;
	std::vector<double> finalresult;
	std::vector<std::map<CAPTK::ImageModalityType, std::string>> QualifiedSubjects = LoadQualifiedSubjectsFromGivenDirectory(inputdirectory, CAPTK::MachineLearningApplicationSubtype::TESTING);
  if (QualifiedSubjects.size() == 0)
  {
    std::cout << "No subject found with required input. Exiting...." << std::endl;
    return EXIT_FAILURE;
  }
  else
  {
    std::cout << "Number of subjects with required input: " << QualifiedSubjects.size() << std::endl;
    EGFRvIIIIndexPredictor objEGFRvIIIPredictor;
    VectorDouble result = objEGFRvIIIPredictor.EGFRvIIIPredictionOnExistingModel(modeldirectory, inputdirectory, QualifiedSubjects, outputdirectory);
    for (unsigned int subjectID = 0; subjectID < QualifiedSubjects.size(); subjectID++)
    {
      std::map<CAPTK::ImageModalityType, std::string> onesubject = QualifiedSubjects[subjectID];
      if (result[subjectID] > 0)
        std::cout << static_cast<std::string>(onesubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]) << ": Score=" << std::to_string(result[subjectID]) << " Mutant." << std::endl;
      else
        std::cout << static_cast<std::string>(onesubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]) << ": Score=" << std::to_string(result[subjectID]) << " Wildtype." << std::endl;
    }
  }
	return EXIT_SUCCESS;
}

int PrepareNewEGFRvIIIPredictionModel(const std::string inputdirectory,const std::string outputdirectory)
{
	std::cout << "Module loaded: Prepare EGFRvIII Prediction Model." << std::endl;
	std::vector<double> finalresult;
  std::vector<std::map<CAPTK::ImageModalityType, std::string>> QualifiedSubjects = LoadQualifiedSubjectsFromGivenDirectory(inputdirectory, CAPTK::MachineLearningApplicationSubtype::TRAINING);
	EGFRvIIIIndexPredictor objEGFRvIIIPredictor;
	std::cout << "Number of subjects with required input: " << QualifiedSubjects.size() << std::endl;
	if (QualifiedSubjects.size() == 0)
		std::cout << "No subject found with required input. Exiting...." << std::endl;
	else
		objEGFRvIIIPredictor.PrepareNewEGFRvIIIPredictionModel(inputdirectory, QualifiedSubjects, outputdirectory);
	return EXIT_SUCCESS;
}

int main(int argc, char **argv)
{
	cbica::CmdParser parser = cbica::CmdParser(argc, argv, "EGFRvIIIIndexPredictor");
	parser.addRequiredParameter("t", "type", cbica::Parameter::STRING, "", "The option of preparing a new model (=0), and for testing on an existing model (=1)");
	parser.addRequiredParameter("i", "input", cbica::Parameter::STRING, "", "The input directory having test subjects");
	parser.addOptionalParameter("m", "model", cbica::Parameter::STRING, "", "The directory having SVM models");
	parser.addRequiredParameter("o", "output", cbica::Parameter::STRING, "", "The output direcory to write output");
	parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");
  //parser.exampleUsage("EGFRvIIIIndexPredictor -i <input dir> -t 0 -o <output dir>");
  parser.addExampleUsage("-t 0 -i C:/properly/formatted/inputDir -o C:/outputDir", "Trains a new model based on the samples in inputDir");
  parser.addExampleUsage("-t 1 -i C:/input -m C:/model -o C:/output", "Tests an existing model for inputs in 'C:/input' based on 'C:/model' ");
  parser.addApplicationDescription("EGFRvIII Index Training and Prediction application");

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
	if (applicationType == CAPTK::MachineLearningApplicationSubtype::TESTING && modelDirectoryName=="")
	{
		std::cout << "Please specify a directory having model file."<<std::endl;
		return EXIT_FAILURE;
	}
	std::cout << "Input directory name:" << inputDirectoryName<<std::endl;
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
    if (cbica::isFile(modelDirectoryName + "/VERSION.yaml"))
    {
      if (!cbica::IsCompatible(modelDirectoryName + "/VERSION.yaml"))
      {
        std::cerr << "The version of model is incompatible with this version of CaPTk.\n";
        return EXIT_FAILURE;
      }
    }
    EGFRvIIIPredictionOnExistingModel(modelDirectoryName, inputDirectoryName, outputDirectoryName);
  }
	else if (applicationType == CAPTK::MachineLearningApplicationSubtype::TRAINING)
		PrepareNewEGFRvIIIPredictionModel(inputDirectoryName, outputDirectoryName);
	else
	{
		parser.echoVersion();
		return EXIT_SUCCESS;
	}

	std::cout << "Finished successfully\n";

	return EXIT_SUCCESS;
}