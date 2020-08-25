#include <time.h>
#include "CaPTkUtils.h"
#include "SBRT_Analysis.h"

#include "cbicaCmdParser.h"

using namespace std;

int main(int argc, char** argv)
{
  clock_t processStart = clock();

  const int imageDimension = 3;	// image dimension, now support 3D only

  // 
  cbica::CmdParser parser(argc, argv, "SBRT_Analysis");
  parser.addRequiredParameter("i", "inputImage", cbica::Parameter::FILE, "none", "Absolute path of PET image", "For example pet.nii.gz");
  parser.addRequiredParameter("m", "maskImage", cbica::Parameter::FILE, "none", "Absolute path of mask image", "For example mask.nii.gz");
  parser.addOptionalParameter("l", "label", cbica::Parameter::INTEGER, "none", "Label value of the ROI", "For example 2");
  parser.addOptionalParameter("o", "outputFile", cbica::Parameter::FILE, "none", "Absolute path and basename of output file (without extension)", "For example radiomic_feature");
  parser.addOptionalParameter("L", "logFile", cbica::Parameter::FILE, "none", "Absolute path of log file", "For example log_file.txt");
  parser.addRequiredParameter("D", "Directory", cbica::Parameter::DIRECTORY, "none", "Absolute path of model directory", "For example C:/Model");
  parser.addExampleUsage("-i C:/PET.nii.gz -m C:/mask.nii.gz -l 1",
    "This will calculate and display the predicted risks for survival and nodal failure");
  parser.addApplicationDescription("This application does automatic prediction of risks for survival and nodal failure.");

  std::string inputFileName;
  std::string maskName;
  int roiLabel = 1;
  std::string oname;
  int outputFea = 0;
  std::string logName;
  std::string modelDir;

  if (parser.isPresent("i"))
  {
    parser.getParameterValue("i", inputFileName);
  }
  if (parser.isPresent("m"))
  {
    parser.getParameterValue("m", maskName);
  }
  if (parser.isPresent("l"))
  {
    parser.getParameterValue("l", roiLabel);
  }
  if (parser.isPresent("o"))
  {
    outputFea = 1;
    parser.getParameterValue("o", oname);
    auto temp = cbica::getFilenamePath(oname, false);
    if (!cbica::isDir(temp))
    {
      cbica::createDir(temp);
    }
  }
  if (parser.isPresent("L"))
  {
    parser.getParameterValue("L", logName);
  }
  if (parser.isPresent("D"))
  {
    parser.getParameterValue("D", modelDir);
  }

  // for test
  //string inputFileName = "/Users/hmli/Google Drive/Temp/SBRT_test_data/ANON10405/suv_PET_x_WB_CTAC_Body.nii.gz";
    //string maskName = "/Users/hmli/Google Drive/Temp/SBRT_test_data/ANON10405/c_output/seg_segmentation.nii.gz";
    //string oname = "/Users/hmli/Google Drive/Temp/SBRT_test_data/ANON10405/c_output/radiomic_features";

  //string inputFileName = argv[1];
  //string maskName = argv[2];
  //string oname = argv[3];

    //unsigned int roiLabel = 2;
    //unsigned roiLabel = atoi(argv[4]);

    //string metaName = "../../data/meta_fea_proj.txt";
    //string projName = "../../data/triFac_res_cpp_kc3_kr5_pet_cox_coeff_train_all.txt";

  std::string metaName = modelDir + "/meta_fea_proj.txt";
  std::string projName = modelDir + "/triFac_res_cpp_kc3_kr5_pet_cox_coeff_train_all.txt";

  if (inputFileName.empty() || maskName.empty() || modelDir.empty())
  {
	  std::cout << "PET image, mask image and model directory need to be loaded for SBRT Analysis." << std::endl;
	  return 0;
  }
  if (guessImageType(inputFileName) != CAPTK::ImageModalityType::IMAGE_TYPE_PET)
  {
	  std::cout << "Only PET image need to be loaded for SBRT Analysis along with mask." << std::endl;
	  return 0;
  }
  if (cbica::isFile(modelDir + "/VERSION.yaml"))
  {
      if (!cbica::IsCompatible(modelDir + "/VERSION.yaml"))
      {
          std::cerr << "The version of model is incompatible with this version of CaPTk.\n";
          return EXIT_FAILURE;
      }
  }
  
  SBRT_Analysis< float, imageDimension > anaObject;

  if (!logName.empty())
  {
    anaObject.SetLogger(logName);
  }

  anaObject.SetParameters(roiLabel);
  anaObject.Initialization(inputFileName, maskName);
  anaObject.FeaExtraction();
  anaObject.GetPredictedRisk(metaName, projName);

  if (outputFea == 1)
  {
    anaObject.OutputFeature(oname);
  }

  clock_t processEnd = clock();
  double processDuration = (double)(processEnd - processStart) / CLOCKS_PER_SEC;
  std::cout << "the whole process costs: " << processDuration << "s" << std::endl;

  return 0;
}
