#include <sstream>
#include "PopulationAtlases.h"
#include "CaPTkEnums.h"


typedef itk::Image< float, 3 > ImageType;
//PopulationAtlases::~PopulationAtlases()
//{
//}


std::vector<typename ImageType::Pointer> PopulationAtlases::GeneratePopualtionAtlas(const std::string inputdirectory, const std::string inputlabel, const std::string inputatlas, const std::string outputdirectory)
{
	std::vector<typename ImageType::Pointer> PopulationAtlases;
	//read label file

	ifstream file(inputlabel); // declare file stream: http://www.cplusplus.com/reference/iostream/ifstream/

	std::vector<double> labels;
	std::vector<std::string> subjects;

	try
	{
		while (!file.eof())
		{
			std::string value;
			getline(file, value, '\n'); // read a string until next comma: http://www.cplusplus.com/reference/string/getline/
			if (value == "")
				break;
			std::vector<std::string> oneline = cbica::stringSplit(value, ",");
			labels.push_back(std::stoi(oneline[1]));
			subjects.push_back(oneline[0]);
	}
	}
	catch (const std::exception& e1)
	{
		mLastErrorMessage = "Cannot open the file " + inputlabel + ". Error code : " + std::string(e1.what());
		std::cout << mLastErrorMessage << std::endl;
		logger.WriteError(mLastErrorMessage);
		return PopulationAtlases;
	}
	//CSVFileReaderType::Pointer reader = CSVFileReaderType::New();
	//MatrixType dataMatrix;
	//try
	//{
	//	reader->SetFileName(inputlabel);
	//	reader->SetFieldDelimiterCharacter(',');
	//	reader->HasColumnHeadersOff();
	//	reader->HasRowHeadersOff();
	//	reader->Parse();
	//	dataMatrix = reader->GetArray2DDataObject()->GetMatrix();

	//	for (unsigned int i = 0; i < dataMatrix.rows(); i++)
	//		for (unsigned int j = 0; j < dataMatrix.cols(); j++)
	//			labels.push_back(dataMatrix(i, j));
	//}

	std::vector<std::string> files = cbica::filesInDirectory(inputdirectory,false);
	if (labels.size() != files.size())
	{
		mLastErrorMessage = "The number of images in the input direcotry are not equal to the number of labels in the label file. ";
		std::cout << mLastErrorMessage << std::endl; 
		logger.WriteError(mLastErrorMessage);
		return PopulationAtlases;
	}

	for (int images = 0; images < subjects.size(); images++)
	{
		bool found = false;
		for (int j = 0; j < files.size(); j++)
		{
			if (subjects[images] == files[j])
				found = true;
		}
		if (found == false)
		{
			mLastErrorMessage = "The subject: " + subjects[images] + "is not present in the specified input directory.";
			std::cout << mLastErrorMessage << std::endl;
			logger.WriteError(mLastErrorMessage);
			return PopulationAtlases;
		}
	}

	//read the template image
	ImageType::Pointer AtlasImagePointer;
	ImageType::Pointer InputImagePointer;
	try
	{
		AtlasImagePointer = ReadNiftiImage<ImageType>(inputatlas);
	}
	catch (const std::exception& e1)
	{
		mLastErrorMessage = "Cannot read the atlas file " + inputatlas + ".Error code : " + std::string(e1.what());
		std::cout << mLastErrorMessage << std::endl;

		logger.WriteError("Cannot read the atlas file " +inputatlas + ".Error code : " + std::string(e1.what()));
		return PopulationAtlases;
	}
	try
	{
		InputImagePointer = ReadNiftiImage<ImageType>(inputdirectory +"/"+ files[0]);
	}
	catch (const std::exception& e1)
	{
		mLastErrorMessage = "Cannot read the segmentation file " + inputdirectory + "/" + files[0] + ".Error code : " + std::string(e1.what());
		std::cout << mLastErrorMessage << std::endl;
		logger.WriteError(mLastErrorMessage);
		return PopulationAtlases;
	}


	int no_of_atlases = 0;
	for (int i = 0; i < labels.size(); i++)
	{
		if (labels[i] > no_of_atlases)
			no_of_atlases = labels[i];
	}
	if (no_of_atlases==0)
	{
		mLastErrorMessage = "Please specify atleast one label for the atlases.";
		logger.WriteError(mLastErrorMessage);
		std::cout << mLastErrorMessage << std::endl;
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

		for (int j = 0; j < labels.size(); j++)
		{
			if (labels[j] == atlas_no)
			{
				counter++;
				std::cout << "Subject # = " << counter << std::endl;
				ImageType::Pointer currentImagePointer = ReadNiftiImage<ImageType>(inputdirectory + "/" + subjects[j]);
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
