/**
\file  featureExtraction.cxx

\brief Feature Extraction comand line interface.

Dependecies: ITK (module_review, module_skullstrip enabled)

https://www.cbica.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2015 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.cbica.upenn.edu/sbia/software/license.html

*/

#include "cbicaCmdParser.h"
#include "cbicaLogging.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaUtilities.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "FeatureExtraction.h"
#include "itkDOMNodeXMLReader.h"
#include "itkDOMNodeXMLWriter.h"
#include "itkDOMNode.h"


#include "cbicaITKImageInfo.h"

template< class TImageType >
void algorithmRunner(std::string & multipatient_file, std::string& patient_id, std::string& image_path, std::string& modality_names,
  std::string& maskfilename, std::string& modality_selected, std::string& selected_roi, std::string &roi_labels,
  std::string& selected_features, std::string& feature_param, std::string& param_file, std::string &outputfile,
  std::string& application, std::vector<std::string> & image_file_name, std::vector<std::string>& modality)
{
  std::vector<typename TImageType::Pointer> inputimages;
  auto maskInfo = cbica::ImageInfo(maskfilename);
  FeatureExtraction<TImageType> features;

  if (!patient_id.empty() && !image_path.empty())
  {
    std::vector<std::string>imagenames = cbica::stringSplit(image_path, ",");

    //check if all the input images and mask match dimension sapcing and size
    for (size_t i = 0; i < inputimages.size(); i++)
    {
      auto imgInfo = cbica::ImageInfo(imagenames[i]);

      if (imgInfo.GetImageDimensions() != maskInfo.GetImageDimensions())
      {
        std::cout << "The dimensions of the inputimage" << imagenames[i] << "and the mask image doesnot match";
      }
      if (imgInfo.GetImageSize() != maskInfo.GetImageSize())
      {
        std::cout << "The size of the inputimage" << imagenames[i] << "and the mask image doesnot match";

      }
      //if (imgInfo.GetImageSpacing() != maskInfo.GetImageSpacing())
      //{
      //  std::cout << "The spacing of the inputimage" << imagenames[i] << "and the mask image doesnot match";
      //}
    }

    inputimages.push_back(cbica::ReadImage<TImageType>(image_path));
    features.SetInputImages(inputimages, modality_names);
    features.SetSelectedLabels(selected_roi, roi_labels);

    // check if the provided labels are present mask image. if not exit the program 
    typename TImageType::Pointer mask = cbica::ReadImage<TImageType>(maskfilename);
    typedef itk::ImageRegionIterator< TImageType > IteratorType;
    IteratorType  imageIterator(mask, mask->GetRequestedRegion());

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
          //++imageIterator;
          //if (imageIterator.IsAtEnd())
          //{
          //  std::cout << "The provided label does not exist in the mask";
          //  EXIT_FAILURE;
          //}

        }
      }
    }
    features.SetMaskImage(mask);
    features.SetRequestedFetures(param_file, selected_features);

    features.Setouputfilename(outputfile);
    features.Update();
  }


  if (!multipatient_file.empty())
  {

    if (param_file.empty())
    {
      std::cerr << "Please provide parameter file";
      exit(EXIT_FAILURE);
    }
    std::ifstream ifile(multipatient_file);
    if (!ifile.good())
    {
      std::cerr << "Provided file" + multipatient_file + " is not good for reading ";
      exit(EXIT_FAILURE);
    }
    else
    {
      // if CSV file doesn't exist, exit with meaningful message
      if (!cbica::fileExists(multipatient_file))
      {
        std::cerr << "Supplied file name, '" << multipatient_file << "' wasn't found.\n";
        exit(EXIT_FAILURE);
      }
    }
    //std::string dirName_Wrap = dirName;

    // store number of rows in the file - this is used to make the program parallelize-able 
    const size_t numberOfRows = cbica::numberOfRowsInFile(multipatient_file);

    // initialize return dictionary
    std::vector< CSVDict > return_CSVDict;
    return_CSVDict.resize(numberOfRows - 1);

    std::vector< std::vector < std::string > > allRows; // store the entire data of the CSV file as a vector of colums and rows (vector< rows <cols> >)

    std::ifstream inFile(multipatient_file.c_str());
    std::string csvPath = cbica::getFilenamePath(multipatient_file);

    for (size_t i = 0; i < numberOfRows; i++)
    {
      std::string line;
      std::getline(inFile, line, '\n');
      line.erase(std::remove(line.begin(), line.end(), '"'), line.end());
      allRows.push_back(cbica::stringSplit(line, ","));
    }

    inFile.close();

    auto fileChecker = [=](std::string file, std::string path) -> std::string
    {
      if (!cbica::fileExists(file))
      {
        return path + file;
      }
    };


    for (int j = 0; j < allRows.size(); j++)
    {
      for (int k = 0; k < allRows[0].size(); k++)
      {
        if (j > 0)
        {
          if (allRows[0][k] == "Images")
          {
            std::stringstream imagenames(allRows[j][k]); std::string parsed;
            while (std::getline(imagenames, parsed, '|'))
            {
              image_file_name.push_back(parsed);
            }
          }
          if (allRows[0][k] == "Modalities")
          {
            std::stringstream modalitynames(allRows[j][k]); std::string parsed;
            while (std::getline(modalitynames, parsed, '|'))
            {
              modality.push_back(parsed);
            }
          }
          if (allRows[0][k] == "ROIFile")
          {
            maskfilename = allRows[j][k];
          }
          if (allRows[0][k] == "SELECTED_ROI")
          {
            selected_roi = allRows[j][k];

          }
          if (allRows[0][k] == "ROI_LABEL")
          {
            roi_labels = allRows[j][k];
          }

        }

        features.SetSelectedLabels(selected_roi, roi_labels);

        for (int i = 0; i < image_file_name.size(); i++)
        {
          if (cbica::ImageInfo(image_file_name[i]).GetImageDimensions() != 3)
          {
            std::cerr << "Non-3D processing will be added in a future release. Please follow discussion on our NITRC page.\n";
            exit(EXIT_FAILURE);
          }

          inputimages.push_back(cbica::ReadImage<TImageType>(image_file_name[i]));
        }

        features.SetInputImages(inputimages, modality_names);
        // check if the provided labels are present mask image. if not exit the program 
        typename    TImageType::Pointer mask = cbica::ReadImage<TImageType>(maskfilename);
        typedef itk::ImageRegionIterator< TImageType > IteratorType;
        IteratorType  imageIterator(mask, mask->GetRequestedRegion());

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
                std::cout << "The provided label does not exist in the mask";
                EXIT_FAILURE;
              }

            }
          }
        }
        features.SetMaskImage(mask);
        features.Update();
      }
    }

  }
}

//! Custom function to avoid some kind of repetition going on in the main code
std::pair< std::string, std::map< std::string, double > > customFunction(cbica::CmdParser &currentCmdParser, const std::string &inputString, const std::string &paramFile)
{
  std::string tempFeatureParam, parsed;
  if (currentCmdParser.isPresent(inputString))
  {
    currentCmdParser.getParameterValue(inputString, tempFeatureParam);
  }
  std::map<std::string, double> parameters;
  if (!tempFeatureParam.empty())
  {
    std::stringstream param_stringstream(tempFeatureParam);
    while (std::getline(param_stringstream, parsed, ','))
    {
      std::stringstream ss(parsed);
      std::vector<std::string> tokens; std::string  item;
      while (std::getline(ss, item, '=')) {
        tokens.push_back(item);
      }
      parameters[tokens[0]] = stof(tokens[1]);
    }
  }
  else (!paramFile.empty());
  {
    //check for glcm features and rplace from command line values
  }

  return std::make_pair(inputString, parameters);
}


std::vector<std::string> getNextLineAndSplitIntoTokens(std::istream& str)
{
  std::vector<std::string>   result;
  std::string                line;
  std::getline(str, line);

  std::stringstream          lineStream(line);
  std::string                cell;

  while (std::getline(lineStream, cell, ','))
  {
    result.push_back(cell);
  }
  // This checks for a trailing comma with no data after it.
  if (!lineStream && cell.empty())
  {
    // If there was a trailing comma then add an empty element.
    result.push_back("");
  }
  return result;
}

int main(int argc, char** argv)
{

  cbica::CmdParser parser(argc, argv);

  parser.addOptionalParameter("b", "batchfile", cbica::Parameter::FILE, ".csv", "Input file with Multi-Patient Multi-Modality details", "Header format is as follows:", "'PATIENT_ID,IMAGES,MASK,ROI,SELECTED_ROI,ROI_LABEL,", "SELECTED_FEATURES,PARAM_FILE'");

  parser.addOptionalParameter("n", "name_patient", cbica::Parameter::STRING, "none", "Patient id");
  parser.addOptionalParameter("i", "imagepath", cbica::Parameter::STRING, "none", "Absolute path of each modality", "Delineate by ','");
  parser.addOptionalParameter("m", "maskpath", cbica::Parameter::STRING, "none", "Absolute path of mask");

  parser.addOptionalParameter("t", "imagetype", cbica::Parameter::STRING, "none", "Names of modalities to be processed", "Delineate by ','");

  parser.addOptionalParameter("r", "roi", cbica::Parameter::STRING, "none", "List of roi for which feature extraction is to be done", "Delineate by ','", "Usage: -r 1,2,5");
  parser.addOptionalParameter("l", "labels", cbica::Parameter::STRING, "none", "Labels variables for selected roi numbers", "Delineate by ','", "Usage: -l Edema, Necorsis");
  parser.addOptionalParameter("f", "features", cbica::Parameter::STRING, "none", "Select features to be calculated . If this option is not used all features listed in param file is calculated", "Delineate by ','", "An exemplary scenario is -f GLCM,GLRLM");

  // parser.addOptionalParameter("fGLCM", "GLCM_parameter", cbica::Parameter::STRING, "none", "Each individual featurename to be parameterized", "Delineate by ','");
  // parser.addOptionalParameter("fGLRLM", "GLRLM_parameter", cbica::Parameter::STRING, "none", "Each individual featurename to be parameterized", "Delineate by ','");
  // parser.addOptionalParameter("fLBP", "LBP_parameter", cbica::Parameter::STRING, "none", "Each individual featurename to be parameterized", "Delineate by ','");
  // parser.addOptionalParameter("fGLSZ", "GLSZ_parameter", cbica::Parameter::STRING, "none", "Each individual featurename to be parameterized", "Delineate by ','");
  // parser.addOptionalParameter("fNGTDM", "NGTDM_parameter", cbica::Parameter::STRING, "none", "Each individual featurename to be parameterized", "Delineate by ','");

  parser.addRequiredParameter("p", "paramfile", cbica::Parameter::FILE, ".csv", "A csv file with all features and its parameters filled");
  //parser.addOptionalParameter("a", "application", cbica::Parameter::STRING, "none", "The application name / default if not specified");

  parser.addOptionalParameter("o", "outputFile", cbica::Parameter::FILE, "none", "Absolute path of Filename ending in .csv" "For example features.csv");
  parser.addOptionalParameter("L", "Logger", cbica::Parameter::FILE, "Log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");

  std::string loggerFile, multipatient_file, patient_id, image_path, modality_names,
    maskfilename, modality_selected, selected_roi, roi_labels, selected_features, feature_param, param_file, outputfile, application;

  std::vector<std::string> roi_label; std::vector<int> roi;  std::vector<std::string> feature_list;
  bool loggerRequested = false;
  cbica::Logging logger;


  if (parser.isPresent("f"))
  {
    parser.getParameterValue("f", selected_features);
  }

  if (parser.isPresent("i"))
  {
    parser.getParameterValue("i", image_path);
  }

  if (parser.isPresent("r"))
  {
    parser.getParameterValue("r", selected_roi);
  }

  if (parser.isPresent("l"))
  {
    parser.getParameterValue("l", roi_labels);
  }

  if (parser.isPresent("t"))
  {
    parser.getParameterValue("t", modality_names);
  }

  if (parser.isPresent("m"))
  {
    parser.getParameterValue("m", maskfilename);
  }

  if (parser.isPresent("L"))
  {
    parser.getParameterValue("L", loggerFile);
    loggerRequested = true;
    logger.UseNewFile(loggerFile);
  }

  if (parser.isPresent("o"))
  {
    parser.getParameterValue("o", outputfile);
  }

  if (parser.isPresent("p"))
  {
    parser.getParameterValue("p", param_file);
    if (!cbica::fileExists(param_file))
    {
      std::cout << "The input param file name supllied is inavalid:" << param_file;
    }
  }

  if (parser.isPresent("b"))
  {
    parser.getParameterValue("b", multipatient_file);
  }

  if (parser.isPresent("i"))
  {
    parser.getParameterValue("i", patient_id);
  }

  if (multipatient_file.empty() && patient_id.empty())
  {
    std::cerr << "NO INPUT PROVIDED" << "\n" <<
      "Please provide a patient id or path to csv file containing patient details";
    return EXIT_FAILURE;
  }
  if (!multipatient_file.empty() && !patient_id.empty())
  {
    std::cerr << "MULTIPLE INPUT PROVIDED" << "\n" <<
      "Please provide either a patient id or path to csv file containing patient details";
    return EXIT_FAILURE;
  }

  std::vector<std::string> image_file_name, modality;

  if (!maskfilename.empty())
  {
    auto maskInfo = cbica::ImageInfo(maskfilename);
    switch (maskInfo.GetImageDimensions())
    {
    case 2:
    {
      using ImageType = itk::Image < float, 2 >;
      algorithmRunner< ImageType >(multipatient_file, patient_id, image_path, modality_names, maskfilename, modality_selected,
        selected_roi, roi_labels, selected_features, feature_param, param_file, outputfile, application, image_file_name, modality);
      break;
    }
    case 3:
    {
      using ImageType = itk::Image < float, 3 >;
      algorithmRunner< ImageType >(multipatient_file, patient_id, image_path, modality_names, maskfilename, modality_selected,
        selected_roi, roi_labels, selected_features, feature_param, param_file, outputfile, application, image_file_name, modality);
      break;
    }
    default:
      break;
    }
  }

  std::cout << "Finished.\n";
  return EXIT_SUCCESS;
}
