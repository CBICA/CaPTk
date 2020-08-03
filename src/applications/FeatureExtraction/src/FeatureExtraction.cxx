/**
\file  featureExtraction.cxx

\brief Feature Extraction comand line interface.

Dependecies: ITK (module_review, module_skullstrip enabled)

https://www.cbica.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2015 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.cbica.upenn.edu/sbia/software/license.html

*/

#include "FeatureExtraction.h"

#include "cbicaCmdParser.h"
#include "cbicaLogging.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaUtilities.h"
//#include "itkImageFileReader.h"
//#include "itkImageFileWriter.h"
//#include "itkDOMNodeXMLReader.h"
//#include "itkDOMNodeXMLWriter.h"
//#include "itkDOMNode.h"

//#include "cbicaITKImageInfo.h"

// stuff used in the program
std::string loggerFile, multipatient_file, patient_id = "DEFAULT", image_path_string, modalities_string, maskfilename, 
selected_roi_string = "all", roi_labels_string = "all", param_file, outputDir, offset_String, outputFilename;

bool debug = false, debugWrite = false, verticalConc = false, featureMaps = false;

int threads = 1;

std::vector< std::string > modality_names, image_paths, selected_roi, roi_labels;


//! The main algorithm, which is templated across the image type
template< class TImageType >
void algorithmRunner(std::vector<typename TImageType::Pointer> inputImages, typename TImageType::Pointer inputMask)
{
  FeatureExtraction<TImageType> features;

  if (inputMask.IsNull())
  {
    if (!cbica::isFile(maskfilename))
    {
      std::cerr << "Mask file is needed [use parameter '-m' on command line for single subject or the header 'ROIFile' in batch file]; SubjectID: '" << patient_id << "'\n";
      //exit(EXIT_FAILURE);
      return;
    }
  }
  if (modality_names.empty())
  {
    std::cerr << "Modality name(s) is needed [use parameter '-t' on command line for single subject or the header 'Modalities' in batch file]; SubjectID: '" << patient_id << "'\n";
    //exit(EXIT_FAILURE);
    return;
  }
  if (inputImages.empty())
  {
    if (image_paths.empty())
    {
      std::cerr << "Input images are needed [use parameter '-i' on command line for single subject or the header 'Images' in batch file]; SubjectID: '" << patient_id << "'\n";
      //exit(EXIT_FAILURE);
      return;
    }
    if (image_paths.size() != modality_names.size())
    {
      std::cerr << "Number of images and modalities should be the same; SubjectID: '" << patient_id << "'\n";
      //exit(EXIT_FAILURE);
      return;
    }
  }
  else
  {
    if (inputImages.size() != modality_names.size())
    {
      std::cerr << "Number of images and modalities should be the same; SubjectID: '" << patient_id << "'\n";
      //exit(EXIT_FAILURE);
      return;
    }
  }

  if (patient_id.empty())
  {
    std::cerr << "Patient name or ID is needed [use parameter '-n' on command line for single subject or the header 'PATIENT_ID' in batch file]; SubjectID: '" << patient_id << "'\n";
    //exit(EXIT_FAILURE);
    return;
  }
  if (selected_roi.empty())
  {
    std::cout << "No ROI values have been selected for patient_id '" << patient_id << "', computation shall be done on all ROIs present in mask.\n";
    //selected_roi_string = "all";
  }
  if (roi_labels.empty())
  {
    std::cout << "No ROI labels have been provided for patient_id '" << patient_id << "', the ROI values will be used as labels instead.\n";
    //roi_labels_string = "all";
  }
  //param_file = cbica::dos2unix(param_file, outputDir);
  std::vector< std::string > imageNames = image_paths;

  // if input images have not been defined before, read them from the paths
  if (inputImages.empty())
  {
    //check if all the input images and mask match dimension spacing and size
    for (size_t i = 0; i < imageNames.size(); i++)
    {
      if (cbica::isDir(imageNames[i]))
      {
        std::cerr << "Images cannot have directory input. Please use absolute paths; SubjectID: '" << patient_id << "'\n";
        //exit(EXIT_FAILURE);
        return;
      }

      auto currentImage = cbica::ReadImage< TImageType >(imageNames[i]);

      if (inputMask.IsNull())
      {
        if (!cbica::ImageSanityCheck(imageNames[i], maskfilename))
        {
          std::cerr << "The input images and mask are not defined in the same physical space; SubjectID: '" << patient_id << "'\n";
          //exit(EXIT_FAILURE);
          return;
        }
      }
      else
      {
        if (!cbica::ImageSanityCheck< TImageType >(currentImage, inputMask))
        {
          std::cerr << "The input images and mask are not defined in the same physical space; SubjectID: '" << patient_id << "'\n";
          //exit(EXIT_FAILURE);
          return;
        }
      }
      inputImages.push_back(currentImage);
    }

    if (inputImages.size() != imageNames.size())
    {
      std::cerr << "Not all images are in the same space as the mask, only those that pass this sanity check will be processed.\n";
    }
  }
   
  if (debug)
  {
    std::cout << "[DEBUG] Initializing FE class.\n";
    features.EnableDebugMode();
  }
  if (debugWrite)
  {
    features.EnableWritingOfIntermediateFiles();
  }

  features.SetPatientID(patient_id);
  features.SetInputImages(inputImages, modality_names);
  if (selected_roi.empty() || roi_labels.empty())
  {
    features.SetSelectedROIsAndLabels(selected_roi_string, roi_labels_string);
  }
  else
  {
    features.SetSelectedROIsAndLabels(selected_roi, roi_labels);
  }

  // check if the provided labels are present mask image. if not exit the program 
  typename TImageType::Pointer mask = cbica::ReadImage< TImageType >(maskfilename);
  typedef itk::ImageRegionIterator< TImageType > IteratorType;
  IteratorType  imageIterator(mask, mask->GetBufferedRegion());

  if (debug)
  {
    std::cout << "[DEBUG] Checking if selected ROIs are present in mask or not.\n";
  }
  auto selectedROIs = features.GetSelectedROIs();
  if (selectedROIs.size() != 0)
  {
    for (size_t x = 0; x < selectedROIs.size(); x++)
    {
      imageIterator.GoToBegin();
      while (!imageIterator.IsAtEnd())
      {
        if (imageIterator.Get() == selectedROIs[x])
          break;
        ++imageIterator;
        if (imageIterator.IsAtEnd())
        {
          std::cout << "The ROI for calculation, '" << std::to_string(selectedROIs[x]) << "' does not exist in the mask, '" << maskfilename << "'.\n";
          //exit(EXIT_FAILURE);
          return;
        }
      }
    }
  }

  if (debug)
  {
    //std::cout << "[DEBUG] All configuration tests have passed, setting up the FE class based on the following parameters:\n"; 
    //std::cout << "[DEBUG] Patient ID: " << patient_id << "\n";
    //std::cout << "[DEBUG] Images    : " << image_path_string << "\n";
    //std::cout << "[DEBUG] Modalities: " << modalities_string << "\n";
    //std::cout << "[DEBUG] Mask File : " << maskfilename << "\n";
    //std::cout << "[DEBUG] ROI Values: " << selected_roi_string << "\n";
    //std::cout << "[DEBUG] ROI Labels: " << roi_labels_string << "\n";
    std::cout << "[DEBUG] Param File: " << param_file << "\n";
  }

  features.SetValidMask();
  features.SetMaskImage(mask);
  features.SetRequestedFeatures(param_file);
  features.SetOutputFilename(outputFilename);
  features.SetVerticallyConcatenatedOutput(verticalConc);
  features.SetWriteFeatureMaps(featureMaps);
  features.SetNumberOfThreads(threads);
  features.Update();
  outputFilename = features.GetOutputFile();
}

//! wrapper for the main algorithm
template< class TImageType >
void algorithmRunner()
{
  std::vector<typename TImageType::Pointer> tempImages;
  typename TImageType::Pointer tempMask;
  algorithmRunner< TImageType >(tempImages, tempMask);
}

//! Calls cbica::stringSplit() by checking for ",", ";" and "|" as deliminators
std::vector< std::string > splitTheString(const std::string &inputString)
{
  std::vector< std::string > returnVector;
  if (inputString.find(",") != std::string::npos)
  {
    returnVector = cbica::stringSplit(inputString, ",");
  }
  else if (inputString.find("|") != std::string::npos)
  {
    returnVector = cbica::stringSplit(inputString, "|");
  }
  else if (inputString.find(";") != std::string::npos)
  {
    returnVector = cbica::stringSplit(inputString, ";");
  }
  else if(!inputString.empty()) // only a single value is present
  {
    returnVector.push_back(inputString);
  }
  return returnVector;
}

int main(int argc, char** argv)
{
  cbica::CmdParser parser(argc, argv, "FeatureExtraction");

  parser.addOptionalParameter("p", "paramFile", cbica::Parameter::FILE, ".csv", "A csv file with all features and its parameters filled", "Default: '../data/1_params_default.csv'");
  parser.addRequiredParameter("o", "outputDir", cbica::Parameter::DIRECTORY, "none", "Absolute path of directory to save results", "Result can be a CSV or Feature Maps (for lattice)");

  parser.addOptionalParameter("b", "batchFile", cbica::Parameter::FILE, ".csv", "Input file with Multi-Patient Multi-Modality details", "Example: '${CaPTk_InstallDir}/share/featureExtractionBatch/batch_featureExtraction.csv'");

  parser.addOptionalParameter("n", "name_patient", cbica::Parameter::STRING, "none", "Patient id", "Required for single subject mode");
  parser.addOptionalParameter("i", "imagePaths", cbica::Parameter::STRING, "none", "Absolute path of each coregistered modality", "Delineate by ','", "Example: -i c:/test1.nii.gz,c:/test2.nii.gz", "Required for single subject mode");
  parser.addOptionalParameter("t", "imageTypes", cbica::Parameter::STRING, "none", "Names of modalities to be processed", "Delineate by ','", "Example: -t T1,T1Gd", "Required for single subject mode");
  parser.addOptionalParameter("m", "maskPath", cbica::Parameter::STRING, "none", "Absolute path of mask coregistered with images", "Required for single subject mode");
  parser.addOptionalParameter("r", "roi", cbica::Parameter::STRING, "none", "List of roi for which feature extraction is to be done", "Delineate by ','", "Example: -r 1,2", "Required for single subject mode");
  parser.addOptionalParameter("l", "labels", cbica::Parameter::STRING, "none", "Labels variables for selected roi numbers", "Delineate by ','", "Usage: -l Edema,Necrosis", "Required for single subject mode");
  parser.addOptionalParameter("vc", "verticalConc", cbica::Parameter::BOOLEAN, "flag", "Whether vertical concatenation is needed or not", "Horizontal concatenation is useful for training", "Defaults to '0'");
  parser.addOptionalParameter("f", "featureMapsWrite", cbica::Parameter::BOOLEAN, "flag", "Whether downsampled feature maps are written or not", "For Lattice computation ONLY", "Defaults to '0'");
  parser.addOptionalParameter("th", "threads", cbica::Parameter::INTEGER, "1-64", "Number of (OpenMP) threads to run FE on", "Defaults to '1'", "This gets disabled when lattice is disabled");
  parser.addOptionalParameter("of", "offsets", cbica::Parameter::STRING, "none", "Exact offset values to pass on for GLCM & GLRLM", "Should be same as ImageDimension and in the format '<offset1>,<offset2>,<offset3>'", "This is scaled on the basis of the radius", "Example: '-of 0x0x1,0x1x0'");

  parser.addOptionalParameter("d", "debug", cbica::Parameter::BOOLEAN, "True or False", "Whether to print out additional debugging info", "Defaults to '0'");
  parser.addOptionalParameter("dw", "debugWrite", cbica::Parameter::BOOLEAN, "True or False", "Whether to write intermediate files or not", "Defaults to '0'");
  parser.addOptionalParameter("L", "Logger", cbica::Parameter::FILE, "Text file with write access", "Full path to log file to store logging information", "By default, only console output is generated");
  //parser.exampleUsage("FeatureExtraction -n AAAC -i AAAC0_flair_pp_shrunk.nii.gz -p 1_params_default.csv -m AAAC0_flair_pp_shrunk_testTumor.nii.gz -o featExParam1.csv -t FL -r 1 -l ED,NC");
  parser.addOptionalParameter("ut", "unitTest", cbica::Parameter::FILE, "Path to reference output", "Whether to run unit test or not", "Disabled for batch processing");


  parser.addApplicationDescription("This does feature calculation based on the input image(s) and mask");
  parser.addExampleUsage("-n AAAC -i AAAC0_flair_pp_shrunk.nii.gz,AAAC0_t1_pp_shrunk.nii.gz  -t FL,T1 -m AAAC0_flair_pp_shrunk_testTumor.nii.gz -r 1,2 -l ED,NC -p 1_params_default.csv -o featExParam1.csv -vc 1", 
    "This calculates features based for subject 'AAAC' using the images 'AAAC0_flair_pp_shrunk.nii.gz,AAAC0_t1_pp_shrunk.nii.gz' represented by the modalities 'FL,T1' in the region defined by 'AAAC0_flair_pp_shrunk_testTumor.nii.gz' at label '1,2' (output as 'ED,NC') based on the parameter file '1_params_default.csv' with output written into 'featExParam1.csv' with vertical concatenation");

  //bool loggerRequested = false;
  //cbica::Logging logger;

  parser.getParameterValue("o", outputDir);

  // check if the user has passed a file or a directory
  if (cbica::isFile(outputDir) || 
    (cbica::getFilenameExtension(outputDir, false) == ".csv"))
  {
    outputFilename = outputDir;
    outputDir = cbica::getFilenamePath(outputFilename, false);
  }
  else
  {
    outputFilename = cbica::normalizePath(outputDir + "/results.csv");
  }

  if (parser.isPresent("i"))
  {
    parser.getParameterValue("i", image_path_string);
    if (cbica::isDir(image_path_string))
    {
      std::cerr << "Images cannot handle a directory input, please use absolute paths.\n";
      exit(EXIT_FAILURE);
    }
    image_paths = splitTheString(image_path_string);
  }

  if (parser.isPresent("r"))
  {
    parser.getParameterValue("r", selected_roi_string);
    selected_roi = splitTheString(selected_roi_string);
  }

  if (parser.isPresent("vc"))
  {
    parser.getParameterValue("vc", verticalConc);
  }

  if (parser.isPresent("d"))
  {
    parser.getParameterValue("d", debug);
  }

  if (parser.isPresent("dw"))
  {
    parser.getParameterValue("dw", debugWrite);
  }

  if (parser.isPresent("f"))
  {
    parser.getParameterValue("f", featureMaps);
  }

  if (parser.isPresent("th"))
  {
    parser.getParameterValue("th", threads);
  }

  if (parser.isPresent("l"))
  {
    parser.getParameterValue("l", roi_labels_string);
    roi_labels = splitTheString(roi_labels_string);
  }

  if (parser.isPresent("t"))
  {
    parser.getParameterValue("t", modalities_string);
    modality_names = splitTheString(modalities_string);
  }
  else
  {
    modality_names.push_back("DEFAULT");
  }

  if (parser.isPresent("m"))
  {
    parser.getParameterValue("m", maskfilename);
    if (cbica::isDir(maskfilename))
    {
      std::cerr << "Mask File cannot handle a directory input, please use absolute paths.\n";
      exit(EXIT_FAILURE);
    }
    else if (!cbica::fileExists(maskfilename))
    {
      std::cerr << "Mask File cannot be non-existent.\n";
      exit(EXIT_FAILURE);
    }
  }

  //if (parser.isPresent("L"))
  //{
  //  parser.getParameterValue("L", loggerFile);
  //  loggerRequested = true;
  //  logger.UseNewFile(loggerFile);
  //}

  if (debug)
  {
    std::cout << "[DEBUG] Performing dos2unix using CBICA TK function; doesn't do anything in Windows machines.\n";
  }
  if (parser.isPresent("b"))
  {
    parser.getParameterValue("b", multipatient_file);
  }

  if (parser.isPresent("n"))
  {
    parser.getParameterValue("n", patient_id);
  }

  if (parser.isPresent("p"))
  {
    parser.getParameterValue("p", param_file);
    if (!cbica::fileExists(param_file))
    {
      std::cerr << "The param file provided wasn't found: " << param_file << "\n";
      return EXIT_FAILURE;
    }
  }
  else
  {
    auto baseParamFile = "1_params_default.csv";
    auto temp = cbica::normPath(cbica::getExecutablePath() + "/../data/" + baseParamFile);
    if (cbica::isFile(temp))
    {
      param_file = temp;
    }
    else
    {
      if (!patient_id.empty())
      {      
        std::string dataDir = "";
#ifndef APPLE
        dataDir = cbica::normPath(cbica::getExecutablePath() + "/../../data/");
#else
        dataDir = cbica::normPath(cbica::getExecutablePath() + "/../Resources/data/features/";
#endif
        temp = dataDir + baseParamFile;
        if (cbica::isFile(temp))
        {
          param_file = temp;
        }
        else
        {
          std::cerr << "No default param file was found. Please set it manually using '-p'.\n";
          return EXIT_FAILURE;
        }
      }
      else if (!multipatient_file.empty())
      {
        std::cout << "Will look for parameter file definition in the batch file under 'PARAM_FILE' header.\n";
      }      
    }
  }

  if (multipatient_file.empty() && patient_id.empty())
  {
    std::cerr << "NO INPUT PROVIDED" << "\n" <<
      "Please provide a patient id or path to csv file containing patient details";
    return EXIT_FAILURE;
  }
  if (parser.isPresent("of"))
  {
    parser.getParameterValue("of", offset_String);
  }

  bool unitTestRequested = false;
  std::string unitTestReferenceFile;
  if (parser.isPresent("ut"))
  {
    unitTestRequested = true;
    parser.getParameterValue("ut", unitTestReferenceFile);
  }

  if (multipatient_file.empty())
  {
    std::cout << "Single subject computation selected.\n";

    cbica::replaceString(selected_roi_string, ",", "|");
    cbica::replaceString(roi_labels_string, ",", "|");

    if (debug)
    {
      std::cout << "[DEBUG] Patient ID: " << patient_id << "\n";
      std::cout << "[DEBUG] Images: " << image_path_string << "\n";
      std::cout << "[DEBUG] Modalities: " << modalities_string << "\n";
      std::cout << "[DEBUG] Mask File: " << maskfilename << "\n";
      std::cout << "[DEBUG] ROI Values: " << selected_roi_string << "\n";
      std::cout << "[DEBUG] ROI Labels: " << roi_labels_string << "\n";
    }

    auto genericImageInfo = cbica::ImageInfo(image_paths[0]);
    switch (genericImageInfo.GetImageDimensions())
    {
    case 2:
    {
      using ImageType = itk::Image < float, 2 >;
      ImageType::Pointer maskImage;
      auto maskImageInfo = cbica::ImageInfo(maskfilename);
      if (maskImageInfo.GetImageDimensions() == 3)
      {
        if (maskImageInfo.GetImageSize()[2] == 1)
        {
          // this is actually a 2D image so re-process accordingly
          using PresumedImageType = itk::Image< float, 3 >;
          auto temp = cbica::GetExtractedImages< PresumedImageType, ImageType >(cbica::ReadImage< PresumedImageType >(maskfilename));
          maskImage = temp[0];
          maskImage->DisconnectPipeline();
        }
      }
      else
      {
        maskImage = cbica::ReadImage< ImageType >(maskfilename);
        maskImage->DisconnectPipeline();
      }

      // get the input images into a vector
      std::vector< ImageType::Pointer > inputImages;
      inputImages.resize(image_paths.size());
      for (size_t i = 0; i < image_paths.size(); i++)
      {
        // sanity check
        if (!cbica::fileExists(image_paths[i]))
        {
          std::cerr << "Image File '" << image_paths[i] << "' does not exist on the filesystem.\n";
          exit(EXIT_FAILURE);
        }
        inputImages[i] = cbica::ReadImage< ImageType >(image_paths[i]);
        inputImages[i]->DisconnectPipeline(); // ensure the new image exists as a separate memory block
        // sanity check
        if (i > 0)
        {
          if (!cbica::ImageSanityCheck< ImageType >(inputImages[i], inputImages[i - 1]))
          {
            std::cerr << "The images corresponding to the following paths do not have consistent physical dimensions:\t" << image_paths[i - 1] <<
              "\n\t" << image_paths[i] << "\n";
            exit(EXIT_FAILURE);
          }
        }
      }

      // sanity checks
      if (!cbica::ImageSanityCheck< ImageType >(inputImages[0], maskImage))
      {
        std::cerr << "The input image(s) and mask do not have consistent physical dimensions. Please resample and re-try.\n";
        exit(EXIT_FAILURE);
      }

      algorithmRunner< ImageType >(inputImages, maskImage);
      break;
    }
    case 3:
    {
      using ImageType = itk::Image < float, 3 >;
      if (genericImageInfo.GetImageSize()[2] == 1)
      {
        // this is actually a 2D image so re-process accordingly

        using ActualImageType = itk::Image< float, 2 >;
        std::vector< ActualImageType::Pointer > inputImages; // vector of input image pointers - these are used so that I/O operations can be reduced
        inputImages.resize(image_paths.size());
        ActualImageType::Pointer maskImage;

        for (size_t i = 0; i < image_paths.size(); i++)
        {
          auto temp = cbica::GetExtractedImages< ImageType, ActualImageType >(cbica::ReadImage< ImageType >(image_paths[i]));
          inputImages[i] = temp[0];
          inputImages[i]->DisconnectPipeline(); // ensure the new image exists as a separate memory block

          // sanity check
          if (i > 0)
          {
            if (!cbica::ImageSanityCheck< ActualImageType >(inputImages[i], inputImages[i - 1]))
            {
              std::cerr << "The images corresponding to the following paths do not have consistent physical dimensions:\t" << image_paths[i - 1] <<
                "\n\t" << image_paths[i] << "\n";
              exit(EXIT_FAILURE);
            }
          }
        }

        auto maskImageInfo = cbica::ImageInfo(maskfilename);
        if (maskImageInfo.GetImageDimensions() == 3)
        {
          if (maskImageInfo.GetImageSize()[2] == 1)
          {
            // this is actually a 2D image so re-process accordingly
            auto temp = cbica::GetExtractedImages< ImageType, ActualImageType >(cbica::ReadImage< ImageType >(maskfilename));
            maskImage = temp[0];
            maskImage->DisconnectPipeline();
          }
        }
        else
        {
          maskImage = cbica::ReadImage< ActualImageType >(maskfilename);
          maskImage->DisconnectPipeline();
        }

        // sanity checks
        if (!cbica::ImageSanityCheck< ActualImageType >(inputImages[0], maskImage))
        {
          std::cerr << "The input image(s) and mask do not have consistent physical dimensions. Please resample and re-try.\n";
          exit(EXIT_FAILURE);
        }

        // run the algorithm with new images
        algorithmRunner< ActualImageType >(inputImages, maskImage);
      }
      else
      {
        // otherwise, it actually is a 3D image
        algorithmRunner< ImageType >();
      }

      break;
    }
    default:
    {
      std::cerr << "Only 2 and 3 dimension images are supported right now.\n";
      return EXIT_FAILURE;
    }
    }
  }
  else // where multi-subject file is passed 
  {
    unitTestRequested = false; // we are not comparing for batch files
    // cbica::Logging(loggerFile, "Multiple subject computation selected.\n");
    // TBD: use the size of allRows to enable parallel processing, if needed
    std::vector< std::vector < std::string > > allRows; // store the entire data of the CSV file as a vector of columns and rows (vector< rows <cols> >)

    if (debug)
    {
      std::cout << "[DEBUG] Performing dos2unix using CBICA TK function; doesn't do anything in Windows machines.\n";
    }
    
    multipatient_file = cbica::dos2unix(multipatient_file, outputDir);

    std::ifstream inFile(multipatient_file.c_str());
    std::string csvPath = cbica::getFilenamePath(multipatient_file);

    while (inFile.good())
    {
      std::string line;
      std::getline(inFile, line);
      line.erase(std::remove(line.begin(), line.end(), '"'), line.end());
      if (!line.empty())
      {
        allRows.push_back(cbica::stringSplit(line, ","));
      }
    }
    inFile.close();

    auto maxThreads = omp_get_max_threads();
    if (threads == -1)
    {
      threads = maxThreads;
    }
    else if (threads > maxThreads)
    {
      threads = maxThreads;
    }

//#pragma omp parallel for num_threads(threads)
    for (int j = 1; j < static_cast<int>(allRows.size()); j++)
    {
      patient_id.clear();
      image_paths.clear();
      modality_names.clear();
      selected_roi_string.clear();
      roi_labels_string.clear();

      for (size_t k = 0; k < allRows[0].size(); k++)
      {
        auto check_wrap = allRows[0][k];
        std::transform(check_wrap.begin(), check_wrap.end(), check_wrap.begin(), ::tolower);

        if ((check_wrap == "patient_id") || (check_wrap == "patientid") || (check_wrap == "name") || (check_wrap == "patient_name") || (check_wrap == "patientname"))
        {
          patient_id = allRows[j][k];
        }
        if ((check_wrap == "images") || (check_wrap == "inputs") || (check_wrap == "inputimages"))
        {
          std::string parsedImageNames = allRows[j][k];
          image_paths = splitTheString(parsedImageNames);
        }
        if (check_wrap == "modalities")
        {
          std::stringstream modalitynames(allRows[j][k]); std::string parsed;
          while (std::getline(modalitynames, parsed, '|'))
          {
            modality_names.push_back(parsed);
          }
        }
        if ((check_wrap == "roifile") || (check_wrap == "maskfile") || (check_wrap == "roi") || (check_wrap == "mask") || (check_wrap == "segmentation"))
        {
          maskfilename = allRows[j][k];
        }
        if ((check_wrap == "selected_roi") || (check_wrap == "selectedroi"))
        {
          selected_roi_string = allRows[j][k];
          selected_roi = cbica::stringSplit(selected_roi_string, "|");
        }
        if ((check_wrap == "roi_label") || (check_wrap == "roilabel") || (check_wrap == "label") || (check_wrap == "roi_labels") || (check_wrap == "roilabels") || (check_wrap == "labels"))
        {
          roi_labels_string = allRows[j][k];
          roi_labels = cbica::stringSplit(roi_labels_string, "|");
        }
        if ((check_wrap == "paramfile") || (check_wrap == "parameterfile") || (check_wrap == "parameters")|| (check_wrap == "param_file")|| (check_wrap == "parameter_file"))
        {
          param_file = allRows[j][k];
        }
        if ((check_wrap == "outputfile") || (check_wrap == "output") || (check_wrap == "outputDir"))
        {
          outputDir = allRows[j][k];
        }
        //else if (cbica::isDir(outputDir))
        //{
        //  outputDir = outputDir + "/" + patient_id + ".csv";
        //}
      } // end of k-loop

      if (modality_names.empty())
      {
        for (size_t numberOfInputImages = 0; numberOfInputImages < image_paths.size(); numberOfInputImages++)
        {
          modality_names.push_back("DEFAULT");
        }
      }

      // sanity check
      if (image_paths.size() != modality_names.size())
      {
        std::cerr << "Number of images and modalites do not match.\n";
        exit(EXIT_FAILURE);
      }

      if (debug)
      {
        std::cout << "[DEBUG] Patient ID: " << patient_id << "\n";
        std::cout << "[DEBUG] Images: " << image_paths[0];
        std::string temp = modality_names[0];
        for (size_t imageNum = 1; imageNum < image_paths.size(); imageNum++)
        {
          std::cout << "," + image_paths[imageNum];
          temp += "," + modality_names[imageNum];
        }
        std::cout << "\n";
        std::cout << "[DEBUG] Modalities: " << temp << "\n";
        std::cout << "[DEBUG] Mask File: " << maskfilename << "\n";
        std::cout << "[DEBUG] ROI Values: " << selected_roi_string << "\n";
        std::cout << "[DEBUG] ROI Labels: " << roi_labels_string << "\n";
      }
      // cbica::Logging(loggerFile, "Starting computation for subject '" + std::string(j) + "' of '" + std::string(allRows.size() - 1) + "' total subjects.\n");
      //std::cout << "Starting computation for subject '" << j << "' of '" << allRows.size() - 1 << "' total subjects.\n";

      // ensuring different dimensions are handled properly
      auto genericImageInfo = cbica::ImageInfo(image_paths[0]);
      switch (genericImageInfo.GetImageDimensions())
      {
      case 2:
      {
        using ImageType = itk::Image < double, 2 >;
        algorithmRunner< ImageType >();
        break;
      }
      case 3:
      {
        using ImageType = itk::Image < float, 3 >;
        algorithmRunner< ImageType >();
        break;
      }
      default:
      {
        // cbica::Logging(loggerFile, "Only 2 and 3 dimension images are supported right now.\n");
        //return EXIT_FAILURE;
      }
      }
    } // end of j-loop

  }

  if (unitTestRequested)
  {
    if (debug)
    {
      std::cout << "Started comparison.\n";
    }
    // unitTestReferenceFile
    // outputFilename

    // if results don't exist, exit with meaningful message
    if (!cbica::fileExists(unitTestReferenceFile))
    {
      std::cerr << "Could not find the reference file '" << unitTestReferenceFile << "'\n";
      return EXIT_FAILURE;
    }
    if (!cbica::fileExists(outputFilename))
    {
      std::cerr << "Could not find the final output file '" << outputFilename << "'\n";
      return EXIT_FAILURE;
    }

    //std::string dirName_Wrap = dirName;

    // store number of rows in the file
    const size_t numberOfRows = cbica::numberOfRowsInFile(unitTestReferenceFile);
    if (numberOfRows != cbica::numberOfRowsInFile(outputFilename))
    {
      std::cerr << "The number of rows in the output files is inconsistent; canont compare.\n";
      return EXIT_FAILURE;
    }

    std::ifstream file_reference(unitTestReferenceFile.c_str()), file_output(outputFilename.c_str());

    // contruct 2 structs to store the entire CSV results
    for (size_t i = 0; i < numberOfRows; i++)
    {
      std::string line_reference;
      std::getline(file_reference, line_reference, '\n');
      line_reference.erase(std::remove(line_reference.begin(), line_reference.end(), '"'), line_reference.end());
      auto currentRow_reference = cbica::stringSplit(line_reference, ",");

      std::string line_output;
      std::getline(file_output, line_output, '\n');
      line_output.erase(std::remove(line_output.begin(), line_output.end(), '"'), line_output.end());
      auto currentRow_output = cbica::stringSplit(line_output, ",");

      float reference_value = std::atof(currentRow_reference[currentRow_reference.size() - 2].c_str());
      float output_value = std::atof(currentRow_output[currentRow_output.size() - 2].c_str());

      if ((reference_value - output_value) > 1e-6)
      {
        // this is an unacceptable difference
        return EXIT_FAILURE;
      }
    }

  }

  std::cout << "Finished.\n";

  return EXIT_SUCCESS;
}
