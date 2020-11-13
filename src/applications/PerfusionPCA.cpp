#include "PerfusionPCA.h"

PerfusionMapType PerfusionPCA::CombineAndCalculatePerfusionPCA(PerfusionMapType PerfusionDataMap, VariableSizeMatrixType &TransformationMatrix, VariableLengthVectorType &MeanVector, vtkSmartPointer<vtkTable> &ReducedPCAs, vtkSmartPointer<vtkDoubleArray> &variance)
{
  PerfusionMapType RevisedPerfusionMap;

  std::vector<int> sizes;
  VectorVectorDouble CombinedPerfusionFeaturesMap;
  std::cout << " perfusion data map size: " << PerfusionDataMap.size() << std::endl;
  
  int nc = 0;
  for (auto const &mapiterator : PerfusionDataMap)
  {
    VariableSizeMatrixType Features = std::get<1>(mapiterator.second);
	std::cout << " rows in perfusiondatamap # " << nc << " = " << Features.Rows() << std::endl;
    sizes.push_back(Features.Rows());
    for (unsigned int i = 0; i < Features.Rows(); i++) //features.rows = # number of indices in mask
    {
      VectorDouble oneVector;
	  //Features number of columns  = total time points
      for (unsigned int j = 0; j < this->m_TotalTimePoints; j++)
        oneVector.push_back(Features(i, j));
      CombinedPerfusionFeaturesMap.push_back(oneVector);
    }
	nc++;
  }

  FeatureReductionClass m_featureReduction;
  //ReducedPCAs = m_featureReduction.GetDiscerningPerfusionTimePoints(CombinedPerfusionFeaturesMap, TransformationMatrix, MeanVector);
  ReducedPCAs = m_featureReduction.GetDiscerningPerfusionTimePointsDynamic(CombinedPerfusionFeaturesMap, TransformationMatrix, MeanVector, variance);

  std::cout << " written files " << std::endl;

  int start = 0;
  for (unsigned int index = 0; index<sizes.size(); index++)// for (auto const &mapiterator : PerfusionDataMap) 
  {
    VariableSizeMatrixType OnePatietnPerfusionData;
    OnePatietnPerfusionData.SetSize(sizes[index], 10);

    if (index != 0)
      start = start + sizes[index - 1];

    for (int i = start; i < start + sizes[index]; i++)
      for (unsigned int j = 0; j < 10; j++)
        OnePatietnPerfusionData(i - start, j) = ReducedPCAs->GetValue(i, j).ToDouble();

    OnePatietnPerfusionData = ColumnWiseScaling(OnePatietnPerfusionData);
    PerfusionTupleType new_tuple(std::get<0>(PerfusionDataMap[index]), OnePatietnPerfusionData);
    RevisedPerfusionMap[index] = new_tuple;
  }
  return RevisedPerfusionMap;
}

VariableSizeMatrixType PerfusionPCA::ColumnWiseScaling(VariableSizeMatrixType inputdata)
{
  //data(:, i) = (data(:, i) - min(data(:, i))). / (max(data(:, i) - min(data(:, i))));

  int NumberOfSamples = inputdata.Rows(); //this is number of voxels in each mask
  int NumberOfFeatures = inputdata.Cols(); //this is same as 'n' number of PC

  std::cout << " # samples = " << NumberOfSamples << std::endl;
  std::cout << " # features = " << NumberOfFeatures << std::endl;

  VariableSizeMatrixType outputdata;
  outputdata.SetSize(NumberOfSamples, NumberOfFeatures);

  //---------calculate mean and variance for each feature----------------
  VariableLengthVectorType minVector;
  VariableLengthVectorType maxVector;
  minVector.SetSize(NumberOfFeatures);
  maxVector.SetSize(NumberOfFeatures);

  //for each column
  //go over each value in each row, find min and max and rescale between 0 and 255
  for (int featureNo = 0; featureNo < NumberOfFeatures; featureNo++)
  {
    double min = inputdata(0, featureNo);
    double max = inputdata(0, featureNo);
    for (int sampleNo = 0; sampleNo < NumberOfSamples; sampleNo++)
    {
		//if( sampleNo < 100)
		//std::cout << inputdata(sampleNo, featureNo) << " ";
      if (inputdata(sampleNo, featureNo) < min)
        min = inputdata(sampleNo, featureNo);
      if (inputdata(sampleNo, featureNo) > max)
        max = inputdata(sampleNo, featureNo);
    }
	std::cout << std::endl;

	//std::cout << " max: " << max << std::endl;
	//std::cout << " min: " << min << std::endl;

	//std::cout << " scaled data " << std::endl;

	for (int sampleNo = 0; sampleNo < NumberOfSamples; sampleNo++)
	{
		//rescaling to 0 - 255
		outputdata(sampleNo, featureNo) = ((inputdata(sampleNo, featureNo) - min) * 255) / (max - min);
		//std::cout << outputdata(sampleNo, featureNo) << " ";
	}
	std::cout << std::endl;

	//exit(1);
  }

  return outputdata;
}

PerfusionMapType PerfusionPCA::CombineAndCalculatePerfusionPCAForTestData(PerfusionMapType PerfusionDataMap, VariableSizeMatrixType &TransformationMatrix, VariableLengthVectorType &MeanVector)
{
	std::cout << " Entering PerfusionPCA::CombineAndCalculatePerfusionPCAForTestData " << std::endl;
  PerfusionMapType RevisedPerfusionMap;

  std::vector<int> sizes; // number of subjects
  VectorVectorDouble CombinedPerfusionFeaturesMap;
  //size of perfusion data map is equal to the number of subjects
  for (auto const &mapiterator : PerfusionDataMap)
  {
    VariableSizeMatrixType Features = std::get<1>(mapiterator.second);
    sizes.push_back(Features.Rows());
	std::cout << " Features rows: " << Features.Rows() << std::endl;

	//Features.Rows() = total number of voxels in the mask
	//we are filling the onevector with intensities of each pixel of the perfusion data
	//see definition of perfusion data map
    for (unsigned int i = 0; i < Features.Rows(); i++)
    {
      VectorDouble oneVector;
      for (unsigned int j = 0; j < 45; j++) //45 time points
        oneVector.push_back(Features(i, j));

	  //combined Perfusion Features map contains intensities of all pixels of all masks
	  //provided by the user
      CombinedPerfusionFeaturesMap.push_back(oneVector);
    }
  }

  
  std::cout << " combined perfusion features map size: " << CombinedPerfusionFeaturesMap.size() << std::endl;
  std::cout << " sizes: " << sizes.size() << std::endl;

  FeatureReductionClass m_featureReduction;
  VectorVectorDouble ReducedPCAs = m_featureReduction.ApplyPCAOnTestDataWithGivenTransformations(CombinedPerfusionFeaturesMap, TransformationMatrix, MeanVector);

  int start = 0;

  //iterate over all subjects  (sizes.size = number of subjects)
  //sizes[0] = number of voxels in first mask and so on
  for (unsigned int index = 0; index<sizes.size(); index++)// for (auto const &mapiterator : PerfusionDataMap) 
  {
    VariableSizeMatrixType OnePatietnPerfusionData;
    OnePatietnPerfusionData.SetSize(sizes[index], 10); //extracting 10 PCAs?

    if (index != 0)
      start = start + sizes[index - 1];

    for (int i = start; i < start + sizes[index]; i++) //for all voxels in mask
      for (unsigned int j = 0; j < 10; j++)
        OnePatietnPerfusionData(i - start, j) = ReducedPCAs[i][j]; //fill matrix with projected PCA data

    OnePatietnPerfusionData = ColumnWiseScaling(OnePatietnPerfusionData);//scale the columns between 0 and 255
    PerfusionTupleType new_tuple(std::get<0>(PerfusionDataMap[index]), OnePatietnPerfusionData);
    RevisedPerfusionMap[index] = new_tuple;
  }
  return RevisedPerfusionMap;
}

PerfusionPCA::ErrorCode PerfusionPCA::ApplyExistingPCAModel(const int number, const std::string inputdirectory, const std::string outputdirectory, std::vector<std::map<CAPTK::ImageModalityType, std::string>> trainingsubjects, const std::string modelDirectoryName)
{
	//TBD: testing part in PsP + chiharu to provide matlab & python codes
	std::cout << " Entering PerfusionPCA::ApplyExistingPCAModel " << std::endl;
  //PerfusionMapType PerfusionDataMap;

  ////Extracting perfusion data of all the patients and putting in PerfusionDataMap
  //for (unsigned int sid = 0; sid < trainingsubjects.size(); sid++)
  //{
  //  std::cout << "Loading Perfusion Image: " << sid << std::endl;
  //  std::map<CAPTK::ImageModalityType, std::string> currentsubject = trainingsubjects[sid];
  //  ImageType::Pointer LabelImagePointer = cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG]));
  //  auto perfImagePointerNifti = cbica::ReadImage<PerfusionImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION]));
  //  std::vector<ImageType::IndexType> indices;

  //  VariableSizeMatrixType perfusionData = LoadPerfusionData<PerfusionImageType, ImageType>(LabelImagePointer, perfImagePointerNifti, indices);
  //  PerfusionTupleType new_tuple(indices, perfusionData);
  //  PerfusionDataMap[sid] = new_tuple;
  //}


	int nPCs = this->ReadNumberOfPCsFromModel(modelDirectoryName + "/NumberOfPCs.csv");
	std::cout << " number of PCs in model: " << nPCs << std::endl;

  //read all the model parameters
  VariableSizeMatrixType PCA_PERF;
  VariableLengthVectorType Mean_PERF;
  CSVFileReaderType::Pointer reader = CSVFileReaderType::New();
  MatrixType dataMatrix;
  VariableLengthVectorType meanMatrix;
  reader->SetFieldDelimiterCharacter(',');
  reader->HasColumnHeadersOff();
  reader->HasRowHeadersOff();

  reader->SetFileName(modelDirectoryName + "/PCA_PERF.csv");
  reader->Parse();
  dataMatrix = reader->GetArray2DDataObject()->GetMatrix();
  PCA_PERF.SetSize(dataMatrix.rows(), dataMatrix.cols());
  for (unsigned int i = 0; i < dataMatrix.rows(); i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_PERF(i, j) = dataMatrix(i, j);

  //TBD: check if the timepoints in model are same as the input data
  //timepoints in model = # columns in PCA_PERF
  if (PCA_PERF.Rows() != this->m_TotalTimePoints)
  {
	  std::cout << " timepoints in model: " << PCA_PERF.Rows() << std::endl;
	  std::cout << " time points don't match." << std::endl;
	  return ErrorCode::DifferentTimePoints;
  }

  reader->SetFileName(modelDirectoryName + "/Mean_PERF.csv");
  reader->Parse();
  dataMatrix = reader->GetArray2DDataObject()->GetMatrix();
  Mean_PERF.SetSize(dataMatrix.size());
  for (unsigned int i = 0; i < dataMatrix.size(); i++)
    Mean_PERF[i] = dataMatrix(0, i);


  //Apply existing PCA model to the test patient
  PerfusionMapType perfFeatures = CombineAndCalculatePerfusionPCAForTestData(this->m_PerfusionDataMap, PCA_PERF, Mean_PERF);
  std::vector<std::vector<ImageType::Pointer>> RevisedPerfusionImagesOfAllPatients;

  for (unsigned int sid = 0; sid < trainingsubjects.size(); sid++)
  {
    std::cout << "Revising Perfusion Image number: " << sid << std::endl;
	std::cout << " size: " << trainingsubjects.size() << std::endl;
    std::map<CAPTK::ImageModalityType, std::string> currentsubject = trainingsubjects[sid];

    auto perfImagePointerNifti = cbica::ReadImage<PerfusionImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION]));
    std::vector<ImageType::Pointer> PerfusionImageVector = cbica::GetExtractedImages<PerfusionImageType, ImageType>(perfImagePointerNifti);

    std::vector<ImageType::Pointer> OnePatientperfusionImages;

	//blah.....
    for (int i = 0; i < number; i++)
    {
      ImageType::Pointer CurrentTimePoint = PerfusionImageVector[i];
      itk::ImageRegionIteratorWithIndex <ImageType> imageIt(CurrentTimePoint, CurrentTimePoint->GetLargestPossibleRegion());
      imageIt.GoToBegin();
      while (!imageIt.IsAtEnd())
      {
        imageIt.Set(0);
        ++imageIt;
      }
      std::vector<ImageType::IndexType> indices = std::get<0>(perfFeatures[sid]);
      VariableSizeMatrixType revisedPerfData = std::get<1>(perfFeatures[sid]);
      for (int j = 0; j < indices.size(); j++)
        CurrentTimePoint.GetPointer()->SetPixel(indices[j], revisedPerfData(j, i));

      OnePatientperfusionImages.push_back(CurrentTimePoint);
    }
    RevisedPerfusionImagesOfAllPatients.push_back(OnePatientperfusionImages);
  }
  for (int index = 0; index < RevisedPerfusionImagesOfAllPatients.size(); index++)
  {
    std::cout << "Writing Perfusion Image number: " << index << std::endl;
    std::vector<ImageType::Pointer> PCAsOfOnePatient = RevisedPerfusionImagesOfAllPatients[index];
    std::map<CAPTK::ImageModalityType, std::string> currentsubject = trainingsubjects[index];

    for (int index2 = 0; index2 < PCAsOfOnePatient.size(); index2++)
    {
      std::string filename = outputdirectory + "/" + static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]) + "_PCA_" + std::to_string(index2) + ".nii.gz";
	  std::cout << "file written at path: " << filename << std::endl;
      cbica::WriteImage<ImageType>(PCAsOfOnePatient[index2], filename);
    }
  }
  return ErrorCode::NoError;
}

PerfusionPCA::ErrorCode PerfusionPCA::LoadData(std::vector<std::map<CAPTK::ImageModalityType, std::string>> trainingsubjects, std::string &inValidSubject)
{
	//PerfusionMapType PerfusionDataMap;
	//Extracting perfusion data of all the patients and putting in PerfusionDataMap
	for (unsigned int sid = 0; sid < trainingsubjects.size(); sid++)
	{
		std::cout << "Loading Perfusion Image: " << sid << std::endl;
		std::map<CAPTK::ImageModalityType, std::string> currentsubject = trainingsubjects[sid];
		ImageType::Pointer LabelImagePointer = cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG]));
		auto perfImagePointerNifti = cbica::ReadImage<PerfusionImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION]));
		std::vector<ImageType::IndexType> indices;

		//current time points
		int timepoints = perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[3];

		if (m_TotalTimePoints == 0) //true on first instantiation
			m_TotalTimePoints = timepoints; //save the timepoints from first subject read

		//if time points don't match, quit with error message.
		if (timepoints != m_TotalTimePoints)
		{
			//TBD: throw name of data with incorrect time points
			//std::cout << "incorrect subject: " << currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION] << std::endl;
			inValidSubject = currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION];
			std::cout << " Number of time points for all subjects are not equal. Cannot Proceed. Please make sure all subjects have the same number of time points. " << std::endl;
			return ErrorCode::DifferentTimePoints;
		}

		VariableSizeMatrixType perfusionData = LoadPerfusionData<PerfusionImageType, ImageType>(LabelImagePointer, perfImagePointerNifti, indices);
		PerfusionTupleType new_tuple(indices, perfusionData);
		this->m_PerfusionDataMap[sid] = new_tuple;
	}
	return ErrorCode::NoError;
}

bool PerfusionPCA::TrainNewPerfusionModel(const int number, const std::string inputdirectory, const std::string outputdirectory, std::vector<std::map<CAPTK::ImageModalityType, std::string>> trainingsubjects)
{
  //PerfusionMapType PerfusionDataMap;

  //Extracting perfusion data of all the patients and putting in PerfusionDataMap
 // for (unsigned int sid = 0; sid < trainingsubjects.size(); sid++)
 // {
 //   std::cout << "Loading Perfusion Image: " << sid << std::endl;
 //   std::map<CAPTK::ImageModalityType, std::string> currentsubject = trainingsubjects[sid];
 //   ImageType::Pointer LabelImagePointer = cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG]));
 //   auto perfImagePointerNifti = cbica::ReadImage<PerfusionImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION]));
 //   std::vector<ImageType::IndexType> indices;
	//int timepoints = perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[3];
	//if(m_TotalTimePoints == 0)
	//	m_TotalTimePoints = timepoints;
	//if (timepoints != m_TotalTimePoints)
	//{
	//	std::cout << " Number of time points for all subjects are not equal. Cannot Proceed. Please make sure all subjects have the same number of time points. " << std::endl;
	//	return EXIT_FAILURE;
	//}

 //   VariableSizeMatrixType perfusionData = LoadPerfusionData<PerfusionImageType, ImageType>(LabelImagePointer, perfImagePointerNifti, indices);
 //   PerfusionTupleType new_tuple(indices, perfusionData);
 //   PerfusionDataMap[sid] = new_tuple;
 // }

  //combining perfusion data, calcualting PCA
  VariableSizeMatrixType TransformationMatrix;
  VariableLengthVectorType MeanVector;
  vtkSmartPointer<vtkTable> TransformedData;
  vtkSmartPointer<vtkDoubleArray> variance;
  PerfusionMapType perfFeatures = CombineAndCalculatePerfusionPCA(this->m_PerfusionDataMap, TransformationMatrix, MeanVector, TransformedData, variance);

  VariableSizeMatrixType TransformedDataMatrix;
  TransformedDataMatrix.SetSize(TransformedData->GetNumberOfRows(), TransformedData->GetNumberOfColumns());
  for (unsigned int index1 = 0; index1 < TransformedData.GetPointer()->GetNumberOfRows(); index1++)
    for (unsigned int index2 = 0; index2 < TransformedData.GetPointer()->GetNumberOfColumns(); index2++)
      TransformedDataMatrix(index1, index2) = TransformedData->GetValue(index1, index2).ToDouble();

  //create vector with 1 item representing total timepoints in input data
  //vector is needed for csv writer
  //std::vector<int> nPCs = { this->DetermineNumberOfPCsFromVariance(variance)};

  int nPCs = this->DetermineNumberOfPCsFromVariance(variance);
  std::cout << " number of PCs under given threshold: " << nPCs << std::endl;

  //write model files

  //only need 3 files: mean, pca_perf and variance
  WriteCSVFiles(TransformationMatrix, outputdirectory + "/PCA_PERF.csv"); //t x t square matrix
  WriteCSVFiles(MeanVector, outputdirectory + "/Mean_PERF.csv"); //
  if(this->m_PerfusionDataForWholePopulationRequested)
	WriteCSVFiles(TransformedDataMatrix, outputdirectory + "/PCA_Data.csv"); //projected perf data in the reduced dimensionality space for whole population, should be extracted if user asks
  this->WriteNumberOfPCs(nPCs, outputdirectory + "/NumberOfPCs.csv"); //# PCs saved as part of the model
  this->WritevtkArray(variance, outputdirectory + "/PCCumulativeVariance.csv");

  //Putting back in images of respective patients
  std::vector<std::vector<ImageType::Pointer>> RevisedPerfusionImagesOfAllPatients;
  for (unsigned int sid = 0; sid < trainingsubjects.size(); sid++)
  {
    std::cout << "Revising current perfusion image: " << sid << std::endl;
    std::map<CAPTK::ImageModalityType, std::string> currentsubject = trainingsubjects[sid];
    auto perfImagePointerNifti = cbica::ReadImage<PerfusionImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION]));
    std::vector<ImageType::Pointer> PerfusionImageVector = cbica::GetExtractedImages<PerfusionImageType, ImageType>(perfImagePointerNifti);

    ImageTypeFloat4D::RegionType region = perfImagePointerNifti.GetPointer()->GetLargestPossibleRegion();
    ImageTypeFloat4D::IndexType regionIndex;
    ImageTypeFloat4D::SizeType regionSize;
    regionSize[0] = region.GetSize()[0];
    regionSize[1] = region.GetSize()[1];
    regionSize[2] = region.GetSize()[2];
    regionSize[3] = 0;
    regionIndex[0] = 0;
    regionIndex[1] = 0;
    regionIndex[2] = 0;
    regionIndex[3] = 0;

    std::vector<ImageType::Pointer> OnePatientperfusionImages;
	
    for (int i = 0; i < nPCs; i++)
    {
      ImageType::Pointer CurrentTimePoint = PerfusionImageVector[i];
      itk::ImageRegionIteratorWithIndex <ImageType> imageIt(CurrentTimePoint, CurrentTimePoint->GetLargestPossibleRegion());
      imageIt.GoToBegin();
      while (!imageIt.IsAtEnd())
      {
        imageIt.Set(0);
        ++imageIt;
      }
      std::vector<ImageType::IndexType> indices = std::get<0>(perfFeatures[sid]);
      VariableSizeMatrixType revisedPerfData = std::get<1>(perfFeatures[sid]);
      for (int j = 0; j < indices.size(); j++)
        CurrentTimePoint.GetPointer()->SetPixel(indices[j], revisedPerfData(j, i));

      OnePatientperfusionImages.push_back(CurrentTimePoint);
    }
    RevisedPerfusionImagesOfAllPatients.push_back(OnePatientperfusionImages);
  }
  for (int index = 0; index < RevisedPerfusionImagesOfAllPatients.size(); index++)
  {
    std::cout << "Writing Perfusion Image: " << index << std::endl;
    std::vector<ImageType::Pointer> PCAsOfOnePatient = RevisedPerfusionImagesOfAllPatients[index];
    std::map<CAPTK::ImageModalityType, std::string> currentsubject = trainingsubjects[index];

    for (int index2 = 0; index2 < PCAsOfOnePatient.size(); index2++)
    {
      std::string filename = outputdirectory + "/" + static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]) + "_PCA_" + std::to_string(index2) + ".nii.gz";
      cbica::WriteImage<ImageType>(PCAsOfOnePatient[index2], filename);
    }
  }
  return true;
}


void PerfusionPCA::WritevtkArray(vtkDoubleArray* inputdata, std::string filepath)
{
	std::ofstream myfile;
	myfile.open(filepath);
	for (vtkIdType i = 0; i < inputdata->GetNumberOfValues(); i++)
		myfile << std::to_string(inputdata->GetValue(i)) << std::endl;
	myfile.close();
}

int PerfusionPCA::DetermineNumberOfPCsFromVariance(vtkSmartPointer<vtkDoubleArray> variance)
{
	int numberOfPCs = 0;

	if (this->m_NumberOfPCsDefined) //number of PCs provided by user
		numberOfPCs = this->m_NumberOfPCs;
	else if (this->m_VarianceThresholdDefined)//variance threshold provided by the user
	{
		int counter = 0;
		for (vtkIdType i = 0; i < variance->GetNumberOfValues(); i++)
		{
			counter++; //we need to include the last component. 
			if (variance->GetValue(i) > (this->m_VarianceThreshold/100.0))
				break;
		}
		numberOfPCs = counter;
	}
	return numberOfPCs;
}

void PerfusionPCA::SetVarianceThreshold(float threshold)
{
	this->m_VarianceThreshold = threshold;
	this->m_VarianceThresholdDefined = true;
}

void PerfusionPCA::SetNumberOfPCs(int pcs)
{
	this->m_NumberOfPCs = pcs;
	this->m_NumberOfPCsDefined = true;
}

void PerfusionPCA::RequestPerfusionDataWholePopulation(bool request)
{
	this->m_PerfusionDataForWholePopulationRequested = request;
}

void PerfusionPCA::WriteNumberOfPCs(int n, std::string filepath)
{
	std::ofstream myfile;
	myfile.open(filepath);
	myfile << std::to_string(n);
	myfile.close();
}

int PerfusionPCA::ReadNumberOfPCsFromModel(std::string filepath)
{
	int n;
	std::ifstream myfile;
	myfile.open(filepath);
	myfile >> n;
	myfile.close();
	return n;
}
