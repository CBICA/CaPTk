#include <sstream>
#include "PopulationAtlases.h"
#include "CaPTkEnums.h"


typedef itk::Image< float, 3 > ImageType;
//PopulationAtlases::~PopulationAtlases()
//{
//}


std::vector<typename ImageType::Pointer> PopulationAtlases::GeneratePopualtionAtlas(const std::vector<std::string> image_paths, const std::vector<int> atlas_labels, const std::string inputatlas, const int no_of_atlases, const std::string outputdirectory)
{
	std::vector<typename ImageType::Pointer> PopulationAtlases;

	//read the template image
	ImageType::Pointer AtlasImagePointer;
	ImageType::Pointer InputImagePointer;
	try
	{
		AtlasImagePointer = cbica::ReadImage<ImageType>(inputatlas);
	}
	catch (const std::exception& e1)
	{
		mLastErrorMessage = "Cannot read the atlas file " + inputatlas + ".Error code : " + std::string(e1.what());
		std::cout << mLastErrorMessage << std::endl;

		logger.WriteError("Cannot read the atlas file " +inputatlas + ".Error code : " + std::string(e1.what()));
		return PopulationAtlases;
	}

	for (int i = 0; i < no_of_atlases; i++)
	{
		int atlas_no = i + 1;
		int counter = 0;
		std::cout << "Atlas # = " << atlas_no << std::endl;
		ImageType::Pointer atlasimage = ImageType::New();
		atlasimage->CopyInformation(AtlasImagePointer);
		atlasimage->SetRequestedRegion(AtlasImagePointer->GetLargestPossibleRegion());
		atlasimage->SetBufferedRegion(AtlasImagePointer->GetBufferedRegion());
		atlasimage->Allocate();
		atlasimage->FillBuffer(0);
    /* This loop iterates through all the ATLAS_LABELS given in the input .csv file and calculates 
    atlas for each atlas label 
    */
		for (int j = 0; j < atlas_labels.size(); j++)
		{
			if (atlas_labels[j] == atlas_no)
			{
				counter++;
				std::cout << "Subject # = " << counter << std::endl;
				ImageType::Pointer currentImagePointer = cbica::ReadImage<ImageType>(image_paths[j]);
				typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
				IteratorType imIt(currentImagePointer, currentImagePointer->GetLargestPossibleRegion());
				IteratorType atlasIt(atlasimage, atlasimage->GetLargestPossibleRegion());
				imIt.GoToBegin();
				atlasIt.GoToBegin();

				while (!imIt.IsAtEnd())
				{
					if (imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::TUMOR || imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::NONENHANCING)
						atlasIt.Set(atlasIt.Get()+1);

					++atlasIt;
					++imIt;
				}
			}
		}
		PopulationAtlases.push_back(atlasimage);
	}
	return PopulationAtlases;
}


bool PopulationAtlases::CalculateSpatialLocationFeatures(const std::vector<std::string> image_paths,
  const std::string inputatlas, 
  const int number_of_regions,
  VariableSizeMatrixType & AllLocationFeatures,
  const std::string outputdirectory)
{
  std::vector<typename ImageType::Pointer> PopulationAtlases;

  //read the template image
  ImageType::Pointer AtlasImagePointer;
  ImageType::Pointer InputImagePointer;
  try
  {
    AtlasImagePointer = cbica::ReadImage<ImageType>(inputatlas);
  }
  catch (const std::exception& e1)
  {
    mLastErrorMessage = "Cannot read the atlas file " + inputatlas + ".Error code : " + std::string(e1.what());
    std::cout << mLastErrorMessage << std::endl;
    logger.WriteError("Cannot read the atlas file " + inputatlas + ".Error code : " + std::string(e1.what()));
    return false;
  }
  /*This segment of code reads input atlas space images and passes those to 
  GetSpatialLocationFeaturesForFixedNumberOfRegions function for calculation of location features
  Final results are stored in AllLocationFeatures
  */
  AllLocationFeatures.SetSize(image_paths.size(), number_of_regions);
  for (int index1 = 0; index1 < image_paths.size(); index1++)
  {
    try
    {
      InputImagePointer = cbica::ReadImage<ImageType>(image_paths[index1]);
    }
    catch (const std::exception& e1)
    {
      mLastErrorMessage = "Cannot read the segmentation file " + image_paths[index1] + ".Error code : " + std::string(e1.what());
      std::cout << mLastErrorMessage << std::endl;
      logger.WriteError(mLastErrorMessage);
      return false;
    }
    std::cout << "Calculating location features for image: " << image_paths[index1] << std::endl;
    VectorDouble locationfeatures = GetSpatialLocationFeaturesForFixedNumberOfRegions<ImageType>(InputImagePointer,AtlasImagePointer, number_of_regions);
    for (int index2 = 0; index2 < locationfeatures.size(); index2++)
      AllLocationFeatures(index1, index2) = locationfeatures[index2];
  }
  return true;
}