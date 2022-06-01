#include "EGFRvIIISurrogateIndex.h"
#include "cbicaCmdParser.h"
#include "cbicaUtilities.h"
#include "cbicaITKSafeImageIO.h"
#include <regex>

//returns pair of int (return code) and vector<double> (EGFRStatus)
std::pair<int, std::vector<double>> RunForSingleSubject(const std::string& inputFile, const std::string& drawingFile,
    const int userTimePoints = -1)
{
    // returns EXIT_SUCCESS (0) if all good, EXIT_FAILURE(1) otherwise

    if (!cbica::fileExists(inputFile))
    {
        std::cerr << "Input image file not found :'" << inputFile << "'\n";
        return std::make_pair(EXIT_FAILURE, std::vector<double>());
    }

    if (!cbica::fileExists(drawingFile))
    {
        std::cerr << "Input mask file not found :'" << drawingFile << "'\n";
        return std::make_pair(EXIT_FAILURE, std::vector<double>());
    }


    std::cout << "Reading inputs.\n";
    using ImageType = itk::Image< float, 3 >;
    using ImageTypePerfusion = itk::Image< float, 4 >;
    auto perfusionImage = cbica::ReadImage< ImageTypePerfusion >(inputFile);
    auto mask = cbica::ReadImage< ImageType >(drawingFile);
    int timePointsInImage = perfusionImage->GetLargestPossibleRegion().GetSize()[3];

    int timePointsToUse = 0;

    if (userTimePoints == -1)
    {
        timePointsToUse = timePointsInImage;
        std::cout << "Found " << std::to_string(timePointsInImage) << " time points in image." << std::endl;
        std::cout << "If you wish to process fewer time points, specify the -t parameter." << std::endl;
    }

    if (userTimePoints > timePointsInImage)
    {
        std::cerr << "Error: Requested number of time points is greater than those present in the image." << std::endl;
        std::cerr << "To resolve, specify a value within bounds." << std::endl;
        return std::make_pair(EXIT_FAILURE, std::vector<double>());
    }

    // Handle case where time points are specified and valid (within image bounds)
    auto finalPerfusionImage = ImageTypePerfusion::New();
    if (userTimePoints > 1 && userTimePoints <= timePointsInImage)
    {
        std::cout << "Using " + std::to_string(userTimePoints) + " time points from perfusion image." << std::endl;
        typedef itk::ExtractImageFilter<ImageTypePerfusion, ImageTypePerfusion> ExtractFilterType;
        ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
        extractFilter->SetDirectionCollapseToSubmatrix();
        ImageTypePerfusion::RegionType inputRegion = perfusionImage->GetLargestPossibleRegion();
        ImageTypePerfusion::SizeType size = inputRegion.GetSize();
        size[3] = userTimePoints;
        ImageTypePerfusion::IndexType start = inputRegion.GetIndex();
        ImageTypePerfusion::RegionType desiredRegion;
        desiredRegion.SetSize(size);
        desiredRegion.SetIndex(start);
        extractFilter->SetExtractionRegion(desiredRegion);
        extractFilter->SetInput(perfusionImage);
        extractFilter->Update();
        finalPerfusionImage = extractFilter->GetOutput();
        finalPerfusionImage->DisconnectPipeline();

    }
    else // Use all time points instead, no extraction needed
    {
        finalPerfusionImage = perfusionImage;
    }


    std::vector<ImageType::IndexType> nearIndices, farIndices;

    itk::ImageRegionIteratorWithIndex< ImageType > maskIt(mask, mask->GetLargestPossibleRegion());
    for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++maskIt)
    {
        if (maskIt.Get() == 1)
            nearIndices.push_back(maskIt.GetIndex());
        else if (maskIt.Get() == 2)
            farIndices.push_back(maskIt.GetIndex());
    }
    if (nearIndices.size() == 0)
    {
        std::cerr << "Mask file does not have near indices with label=1. \n";
        return std::make_pair(EXIT_FAILURE, std::vector<double>());
    }
    if (farIndices.size() == 0)
    {
        std::cerr << "Mask file does not have far indices with label=2. \n";
        return std::make_pair(EXIT_FAILURE, std::vector<double>());
    }
    EGFRStatusPredictor egfrEstimator;
    auto extension = cbica::getFilenameExtension(inputFile);
    std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
    auto EGFRStatusParams = egfrEstimator.PredictEGFRStatus< ImageType, ImageTypePerfusion >(finalPerfusionImage, nearIndices, farIndices);

    std::cout << "Printing results...\n\n";
    std::cout << "PHI Value = " << EGFRStatusParams[0] << "\n";
    std::cout << "Peak Height Ratio = " << EGFRStatusParams[1] / EGFRStatusParams[2] << "\n";
    std::cout << "Number of near voxels used = " << EGFRStatusParams[3] << "\n";
    std::cout << "Number of far voxels used = " << EGFRStatusParams[4] << "\n";

    return std::make_pair(EXIT_SUCCESS, EGFRStatusParams);
}

int main(int argc, char **argv)
{
  auto parser = cbica::CmdParser(argc, argv, "EGFRvIIISurrogateIndex");
  parser.addOptionalParameter("i", "image", cbica::Parameter::FILE, "NIfTI", "Input Perfusion image on which computation is done");
  parser.addOptionalParameter("m", "mask", cbica::Parameter::FILE, "NIfTI", "Mask containing near (1) and far (2) labels");
  parser.addOptionalParameter("t", "timePoints", cbica::Parameter::INTEGER, ">1", "Number of time points to use for computation.", "Useful for using a uniform number of points across batch subjects.", "(Defaults to the number of time points present in the image)");
  parser.addOptionalParameter("b", "batchDir", cbica::Parameter::DIRECTORY, "", "Directory of subjects (subdirectories) to run PHI on.",  "Each subject subdirectory must have a file with 'DSC' or 'perf' in the basename and a mask with 'mask' in the basename.");
  parser.addOptionalParameter("o", "batchOutputFile", cbica::Parameter::FILE, "", "File to place output in for batch processing, csv");
  parser.addExampleUsage("-i DSC-MRI_data.nii.gz -m Near_Far_masks.nii.gz",
    "Based on the near-far mask and input DSC MRI data, the PHI index is calculated");
  parser.addExampleUsage("-b C:/PHI_Subjects -o C:/PHI_Output.csv -t 45",
    "For each subject dir in C:/PHI_Subjects, like C:/PHI_Subjects/AAAA/, calculate PHI based on 45 time points and place in C:/PHI_Output.csv .");
  parser.addApplicationDescription("Peritumoral Heterogeneity Index calculator");
  
  bool singleSubjectMode = false;
  bool batchMode = false;
  if (parser.isPresent("i") || parser.isPresent("m"))
  {
      singleSubjectMode = true;
  }
  if (parser.isPresent("b") || parser.isPresent("o"))
  {
      batchMode = true;
  }
  if (singleSubjectMode && batchMode)
  {
      std::cerr << "Inconsistent options provided. -i and -m parameters are for single-subject, -b and -o are for batch processing." << std::endl;
      std::cerr << "Please check parameters and try again." << std::endl;
      return EXIT_FAILURE;
  }
  if (!singleSubjectMode && !batchMode)
  {
      std::cerr << "You must specify options for either single-subject mode (-i, -m) or batch mode (-b, -o)." << std::endl;
      std::cerr << "Please check parameters and try again." << std::endl;
      return EXIT_FAILURE;
  }

  bool bRequestedTimePoints = false;
  int userSpecifiedTimePoints = 0;
  if (parser.isPresent("t"))
  {
      bRequestedTimePoints = true;
      std::string timePointsString;
      parser.getParameterValue("t", userSpecifiedTimePoints);

      if(userSpecifiedTimePoints <= 1)
      {
          std::cerr << "The number of provided time points is too low. Specify a value greater than 1." << std::endl;
          return EXIT_FAILURE;
      }
  }

  if (singleSubjectMode)
  {
      std::string inputFile, drawingFile;
      parser.getParameterValue("i", inputFile);
      parser.getParameterValue("m", drawingFile);

      bool singleSubjectRunFailed;
      if (bRequestedTimePoints)
      {
          singleSubjectRunFailed = RunForSingleSubject(inputFile, drawingFile, userSpecifiedTimePoints).first;
      }
      else
      {
          singleSubjectRunFailed = RunForSingleSubject(inputFile, drawingFile).first;
      }
      
      if (singleSubjectRunFailed)
      {
          return EXIT_FAILURE;
      }

  }
  else if (batchMode)
  {
      std::string inputDirStr, outputFileStr;
      parser.getParameterValue("b", inputDirStr);
      parser.getParameterValue("o", outputFileStr);
      std::cout << "Searching in dir " + inputDirStr + " for subject subdirectories." << std::endl;
      if (!cbica::directoryExists(inputDirStr))
      {
          std::cerr << "Could not read directory " + inputDirStr + ". Check existence and permissions and try again." << std::endl;
          return EXIT_FAILURE;
      }

      // Write header
      std::ofstream outFile;
      outFile.open(outputFileStr, ios::trunc);
      if (outFile.is_open())
      {
          outFile << "Subject Name, PHI Value, Peak Height Ratio, Num Near Voxels Used, Num Far Voxels Used\n";
          outFile.close();
      }
      else
      {
          std::cerr << "Couldn't write to output file. Check permissions." << std::endl;
          return EXIT_FAILURE;
      }
      
      auto subjectDirs = cbica::subdirectoriesInDirectory(inputDirStr, false, true);
      //QStringList subjectDirs = inputDir.entryList();
      bool anySubjectsFailed = false;
      for (int i = 0; i < subjectDirs.size(); i++)
      {
          // Subject-specific loop -- in here, don't outright return failure, but flag any failures.
          std::string currentSubjectDir = subjectDirs[i];
          std::string dirpath, dirbase, dirext; // 
          cbica::splitFileName(currentSubjectDir, dirpath, dirbase, dirext);
          std::string currentSubjectName = dirbase;
          std::cout << "Processing subject " + currentSubjectName + " ..." << std::endl;

          std::string possibleInputFile, possibleDrawingFile;
          auto filesInDir = cbica::filesInDirectory(currentSubjectDir, true);
          std::regex perfRegex("(perf|DSC|PERF|dsc)");
          std::regex drawingRegex("(mask|MASK|label|LABEL|near|NEAR|far|FAR)");
          for (int file = 0; file < filesInDir.size(); file++)
          {
              std::string path, base, ext;
              cbica::splitFileName(filesInDir[file], path, base, ext);
              if (std::regex_search(base, perfRegex))
              {
                  possibleInputFile = filesInDir[file];
                  continue;
              }
              if (std::regex_search(base, drawingRegex))
              {
                  possibleDrawingFile = filesInDir[file];
                  continue;
              }
          }
          if (possibleInputFile.empty())
          {
              std::cerr << "No input DSC/perfusion NIFTIs found for subject " + currentSubjectName + " , skipping..." << std::endl;
              anySubjectsFailed = true;
              continue;
          }
          if (possibleDrawingFile.empty())
          {
              std::cerr << "No near/far label NIFTIs found for subject " + currentSubjectName + " , skipping..." << std::endl;
              anySubjectsFailed = true;
              continue;
          }

          bool singleSubjectRunFailed;
          std::pair<int, VectorDouble> resultFromRun;
          if (bRequestedTimePoints)
          {
              resultFromRun = RunForSingleSubject(possibleInputFile, possibleDrawingFile, userSpecifiedTimePoints);
          }
          else
          {
              resultFromRun = RunForSingleSubject(possibleInputFile, possibleDrawingFile);
          }
          singleSubjectRunFailed = resultFromRun.first > 0;
          if (singleSubjectRunFailed)
          {
              anySubjectsFailed = true;
              std::cerr << "Algorithm failed for subject " + currentSubjectName + " . Check logs for details. Skipping..." << std::endl;
              continue;
          }
          else
          {
              std::vector<double> egfrStatus = resultFromRun.second;
              outFile.open(outputFileStr, ios::app);
              if (outFile.is_open())
              {
                  outFile << std::fixed << std::setprecision(4) << currentSubjectName << "," << egfrStatus[0] << ","
                      << egfrStatus[1] / egfrStatus[2] << "," << egfrStatus[3] << "," << egfrStatus[4] << "\n";
                  outFile.close();
              }
              else
              {
                  std::cerr << "Couldn't write to file. Check permissions." << std::endl;
                  return EXIT_FAILURE;
              }
              
          }


          std::cout << "Finished processing subject " + currentSubjectName + " ." << std::endl;
      }

      if (anySubjectsFailed)
      {
          std::cerr << "Failed to process one or more subjects. Check the console logs for more details." << std::endl;
          return EXIT_FAILURE;
      }
      else
      {
          std::cout << "Finished processing all subjects successfully." << std::endl;
      }

  }

  // If we get here without failing...
  return EXIT_SUCCESS;
}

