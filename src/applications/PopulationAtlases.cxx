#include "PopulationAtlases.h"
#include "cbicaUtilities.h"
#include "cbicaCmdParser.h"
#include "CaPTkEnums.h"
#include "CaPTkUtils.h"

int main(int argc, char **argv)
{
  cbica::CmdParser parser = cbica::CmdParser(argc, argv, "PopulationAtlases");
  parser.addRequiredParameter("i", "input", cbica::Parameter::STRING, "", "The input batch file.");
  parser.addRequiredParameter("a", "atlas", cbica::Parameter::STRING, "", "The atlas template.");
  parser.addRequiredParameter("o", "output", cbica::Parameter::STRING, "", "The output directory.");
  parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");
  //parser.exampleUsage("PopulationAtlases -i <input dir> -l atlaslabelfile.csv -a jakob_stripped_with_cere_lps_256256128.nii.gz -o <output dir>");
  parser.addExampleUsage("-i C:/properly/formatted/inputfile -a jakob_stripped_with_cere_lps_256256128.nii.gz -o C:/outputDir",
    "Calculates the population atlas based on the input images, atlas labels and the jakob");
  parser.addApplicationDescription("Population Atlas calculator");

  // parameters to get from the command line
  cbica::Logging logger;
  std::string loggerFile;
  bool loggerRequested = false;

  int tempPosition;
  std::string inputFileName, inputAtlasName, outputDirectoryName, toWrite;
  if (parser.compareParameter("L", tempPosition))
  {
    loggerFile = argv[tempPosition + 1];
    loggerRequested = true;
    logger.UseNewFile(loggerFile);
  }
  if (parser.compareParameter("i", tempPosition))
    inputFileName = argv[tempPosition + 1];

  if (parser.compareParameter("a", tempPosition))
    inputAtlasName = argv[tempPosition + 1];

  if (parser.compareParameter("o", tempPosition))
    outputDirectoryName = argv[tempPosition + 1];

  //read and store the entire data of csv file
  std::vector< std::vector < std::string > > allRows; // store the entire data of the CSV file as a vector of columns and rows (vector< rows <cols> >)
  inputFileName = cbica::dos2unix(inputFileName, outputDirectoryName);
  std::ifstream inFile(inputFileName.c_str());
  std::string csvPath = cbica::getFilenamePath(inputFileName);
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
  std::cout << "input file parsed. Number of entries: " << allRows.size()<< std::endl;

  // sanity check to make sure that the file is not empty
  if (allRows.size() ==0)
  {
    std::cerr << "There is no data in the given file: " << inputFileName << std::endl;
    exit(EXIT_FAILURE);
  }

  //put the data in respective vectors
  std::vector< std::string > patient_ids, image_paths, atlas_labels;
  for (int j = 1; j < allRows.size(); j++)
  {
    for (size_t k = 0; k < allRows[0].size(); k++)
    {
      auto check_wrap = allRows[0][k];
      std::transform(check_wrap.begin(), check_wrap.end(), check_wrap.begin(), ::tolower);

      if (check_wrap == "patient_ids")
        patient_ids.push_back(allRows[j][k]);
      else if (check_wrap == "images")
        image_paths.push_back(allRows[j][k]);
      else if (check_wrap == "atlas_labels")
        atlas_labels.push_back(allRows[j][k]);
    }
  }

  // sanity check to make sure that all patient ids have corresponding atlas numbers and paths
  if (image_paths.size() != patient_ids.size() || image_paths.size() != atlas_labels.size())
  {
    std::cerr << "There is a mismatch in the number of patinet ids, images, and atlas identifiers.\n";
    exit(EXIT_FAILURE);
  }

  //for (int j = 0; j < patient_ids.size(); j++)
  //  std::cout << patient_ids[j] << image_paths[j] << atlas_labels[j] << std::endl;

  //convert atlas labels from string to numbers
  std::vector<int> atlas_labels_numbers;
  for (int i = 0; i < atlas_labels.size(); i++)
    atlas_labels_numbers.push_back(std::stoi(atlas_labels[i]));


  //find number of atlas in the input file. 
  //atlas numbers should in ascending order like, 1,2,3,....,n
  int no_of_atlases = 0;
  for (int i = 0; i < atlas_labels.size(); i++)
  {
    if (atlas_labels_numbers[i] > no_of_atlases)
      no_of_atlases = atlas_labels_numbers[i];
  }
  if (no_of_atlases == 0)
  {
    std::cerr << "Please specify atleast one label for the atlases.";
    exit(EXIT_FAILURE);
  }
  std::cout << "Number of identified atlases: " << no_of_atlases << std::endl;

  //find unique number of regions in the template image
  //region numbers should be in ascending order like, 1,2,3,...,n
  ImageType::Pointer AtlasImagePointer = cbica::ReadImage<ImageType>(inputAtlasName);
  typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
  IteratorType atlasIt(AtlasImagePointer, AtlasImagePointer->GetLargestPossibleRegion());
  atlasIt.GoToBegin();
  int numberofregions = 0;
  while (!atlasIt.IsAtEnd())
  {
    if (atlasIt.Get() > numberofregions)
      numberofregions = atlasIt.Get();
    ++atlasIt;
  }
  if (numberofregions < 2)
  {
    std::cerr << "There should be atleast 2 regions in the atlas file.";
    exit(EXIT_FAILURE);
  }
  std::vector < std::string> region_names;
  for (int index = 0; index < numberofregions; index++)
    region_names.push_back("Location_" + std::to_string(index+1));

  //code to calculate atlases
  PopulationAtlases objPopulationAtlases;
  std::vector<typename ImageType::Pointer> atlases = objPopulationAtlases.GeneratePopualtionAtlas(image_paths, atlas_labels_numbers, inputAtlasName, no_of_atlases, outputDirectoryName);
  if (atlases.size() > 0)
  {
    std::cout << "Writing atlases in the specified output directory.\n";
    for (int i = 0; i < atlases.size(); i++)
    {
      typedef itk::ImageFileWriter< ImageType > WriterType;
      WriterType::Pointer writer1 = WriterType::New();
      writer1->SetFileName(outputDirectoryName + "/Atlas_" + std::to_string(i + 1) + ".nii.gz");
      writer1->SetInput(atlases[i]);
      writer1->Update();
    }
    std::cout << "Atlas writing finished.\n";
  }
  else
  {
    std::cerr << "System encountered an error in calculating atlases.\n";
    exit(EXIT_FAILURE);
  }

  //code to calculate spatial location features
  VariableSizeMatrixType LocationFeaturesAll;
  if (objPopulationAtlases.CalculateSpatialLocationFeatures(image_paths, inputAtlasName, numberofregions, LocationFeaturesAll, outputDirectoryName) == true)
  {
    WriteCSVFilesWithHorizontalAndVerticalHeaders(LocationFeaturesAll, patient_ids,region_names, outputDirectoryName + "/locationfeatures.csv");
    std::cout << "Spatial location features written in the output directory.\n";
  }
  else
  {
    std::cerr <<"System encountered an error in calculating location features.\n";
    exit(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}