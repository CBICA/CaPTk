#include "PerfusionPCA.h"

PerfusionMapType PerfusionPCA::CombineAndCalculatePerfusionPCA(PerfusionMapType PerfusionDataMap, VariableSizeMatrixType &TransformationMatrix, VariableLengthVectorType &MeanVector, vtkSmartPointer<vtkTable> &ReducedPCAs)
{
  PerfusionMapType RevisedPerfusionMap;

  std::vector<int> sizes;
  VectorVectorDouble CombinedPerfusionFeaturesMap;
  for (auto const &mapiterator : PerfusionDataMap)
  {
    VariableSizeMatrixType Features = std::get<1>(mapiterator.second);
    sizes.push_back(Features.Rows());
    for (unsigned int i = 0; i < Features.Rows(); i++)
    {
      VectorDouble oneVector;
      for (unsigned int j = 0; j < 45; j++)
        oneVector.push_back(Features(i, j));
      CombinedPerfusionFeaturesMap.push_back(oneVector);
    }
  }
  FeatureReductionClass m_featureReduction;
  ReducedPCAs = m_featureReduction.GetDiscerningPerfusionTimePoints(CombinedPerfusionFeaturesMap, TransformationMatrix, MeanVector);

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

  int NumberOfSamples = inputdata.Rows();
  int NumberOfFeatures = inputdata.Cols();
  VariableSizeMatrixType outputdata;
  outputdata.SetSize(NumberOfSamples, NumberOfFeatures);

  //---------calculate mean and variance for each feature----------------
  VariableLengthVectorType minVector;
  VariableLengthVectorType maxVector;
  minVector.SetSize(NumberOfFeatures);
  maxVector.SetSize(NumberOfFeatures);

  for (int featureNo = 0; featureNo < NumberOfFeatures; featureNo++)
  {
    double min = inputdata(0, featureNo);
    double max = inputdata(0, featureNo);
    for (int sampleNo = 0; sampleNo < NumberOfSamples; sampleNo++)
    {
      if (inputdata(sampleNo, featureNo) < min)
        min = inputdata(sampleNo, featureNo);
      if (inputdata(sampleNo, featureNo) > max)
        max = inputdata(sampleNo, featureNo);
    }
    for (int sampleNo = 0; sampleNo < NumberOfSamples; sampleNo++)
      outputdata(sampleNo, featureNo) = ((inputdata(sampleNo, featureNo) - min) * 255) / (max - min);
  }

  return outputdata;
}

PerfusionMapType PerfusionPCA::CombineAndCalculatePerfusionPCAForTestData(PerfusionMapType PerfusionDataMap, VariableSizeMatrixType &TransformationMatrix, VariableLengthVectorType &MeanVector)
{
  PerfusionMapType RevisedPerfusionMap;

  std::vector<int> sizes;
  VectorVectorDouble CombinedPerfusionFeaturesMap;
  for (auto const &mapiterator : PerfusionDataMap)
  {
    VariableSizeMatrixType Features = std::get<1>(mapiterator.second);
    sizes.push_back(Features.Rows());
    for (unsigned int i = 0; i < Features.Rows(); i++)
    {
      VectorDouble oneVector;
      for (unsigned int j = 0; j < 45; j++)
        oneVector.push_back(Features(i, j));
      CombinedPerfusionFeaturesMap.push_back(oneVector);
    }
  }
  FeatureReductionClass m_featureReduction;
  VectorVectorDouble ReducedPCAs = m_featureReduction.ApplyPCAOnTestDataWithGivenTransformations(CombinedPerfusionFeaturesMap, TransformationMatrix, MeanVector);

  int start = 0;
  for (unsigned int index = 0; index<sizes.size(); index++)// for (auto const &mapiterator : PerfusionDataMap) 
  {
    VariableSizeMatrixType OnePatietnPerfusionData;
    OnePatietnPerfusionData.SetSize(sizes[index], 10);

    if (index != 0)
      start = start + sizes[index - 1];

    for (int i = start; i < start + sizes[index]; i++)
      for (unsigned int j = 0; j < 10; j++)
        OnePatietnPerfusionData(i - start, j) = ReducedPCAs[i][j];

    OnePatietnPerfusionData = ColumnWiseScaling(OnePatietnPerfusionData);
    PerfusionTupleType new_tuple(std::get<0>(PerfusionDataMap[index]), OnePatietnPerfusionData);
    RevisedPerfusionMap[index] = new_tuple;
  }
  return RevisedPerfusionMap;
}

bool PerfusionPCA::ApplyExistingPCAModel(const int number, const std::string inputdirectory, const std::string outputdirectory, std::vector<std::map<CAPTK::ImageModalityType, std::string>> trainingsubjects, const std::string modelDirectoryName)
{
  PerfusionMapType PerfusionDataMap;

  //Extracting perfusion data of all the patients and putting in PerfusionDataMap
  for (unsigned int sid = 0; sid < trainingsubjects.size(); sid++)
  {
    std::cout << "Loading Perfusion Image: " << sid << std::endl;
    std::map<CAPTK::ImageModalityType, std::string> currentsubject = trainingsubjects[sid];
    ImageType::Pointer LabelImagePointer = cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG]));
    auto perfImagePointerNifti = cbica::ReadImage<PerfusionImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION]));
    std::vector<ImageType::IndexType> indices;

    VariableSizeMatrixType perfusionData = LoadPerfusionData<PerfusionImageType, ImageType>(LabelImagePointer, perfImagePointerNifti, indices);
    PerfusionTupleType new_tuple(indices, perfusionData);
    PerfusionDataMap[sid] = new_tuple;
  }

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

  reader->SetFileName(modelDirectoryName + "/Mean_PERF.csv");
  reader->Parse();
  dataMatrix = reader->GetArray2DDataObject()->GetMatrix();
  Mean_PERF.SetSize(dataMatrix.size());
  for (unsigned int i = 0; i < dataMatrix.size(); i++)
    Mean_PERF[i] = dataMatrix(0, i);

  //Apply existing PCA model to the test patient
  PerfusionMapType perfFeatures = CombineAndCalculatePerfusionPCAForTestData(PerfusionDataMap, PCA_PERF, Mean_PERF);
  std::vector<std::vector<ImageType::Pointer>> RevisedPerfusionImagesOfAllPatients;

  for (unsigned int sid = 0; sid < trainingsubjects.size(); sid++)
  {
    std::cout << "Revising Perfusion Image: " << sid << std::endl;
    std::map<CAPTK::ImageModalityType, std::string> currentsubject = trainingsubjects[sid];

    auto perfImagePointerNifti = cbica::ReadImage<PerfusionImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION]));
    std::vector<ImageType::Pointer> PerfusionImageVector = cbica::GetExtractedImages<PerfusionImageType, ImageType>(perfImagePointerNifti);

    std::vector<ImageType::Pointer> OnePatientperfusionImages;
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
bool PerfusionPCA::TrainNewPerfusionModel(const int number, const std::string inputdirectory, const std::string outputdirectory, std::vector<std::map<CAPTK::ImageModalityType, std::string>> trainingsubjects)
{
  PerfusionMapType PerfusionDataMap;

  //Extracting perfusion data of all the patients and putting in PerfusionDataMap
  for (unsigned int sid = 0; sid < trainingsubjects.size(); sid++)
  {
    std::cout << "Loading Perfusion Image: " << sid << std::endl;
    std::map<CAPTK::ImageModalityType, std::string> currentsubject = trainingsubjects[sid];
    ImageType::Pointer LabelImagePointer = cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG]));
    auto perfImagePointerNifti = cbica::ReadImage<PerfusionImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION]));
    std::vector<ImageType::IndexType> indices;

    VariableSizeMatrixType perfusionData = LoadPerfusionData<PerfusionImageType, ImageType>(LabelImagePointer, perfImagePointerNifti, indices);
    PerfusionTupleType new_tuple(indices, perfusionData);
    PerfusionDataMap[sid] = new_tuple;
  }

  //combining perfusion data, calcualting PCA
  VariableSizeMatrixType TransformationMatrix;
  VariableLengthVectorType MeanVector;
  vtkSmartPointer<vtkTable> TransformedData;
  PerfusionMapType perfFeatures = CombineAndCalculatePerfusionPCA(PerfusionDataMap, TransformationMatrix, MeanVector, TransformedData);

  VariableSizeMatrixType TransformedDataMatrix;
  TransformedDataMatrix.SetSize(TransformedData->GetNumberOfRows(), TransformedData->GetNumberOfColumns());
  for (unsigned int index1 = 0; index1 < TransformedData.GetPointer()->GetNumberOfRows(); index1++)
    for (unsigned int index2 = 0; index2 < TransformedData.GetPointer()->GetNumberOfColumns(); index2++)
      TransformedDataMatrix(index1, index2) = TransformedData->GetValue(index1, index2).ToDouble();

  WriteCSVFiles(TransformationMatrix, outputdirectory + "/PCA_PERF.csv");
  WriteCSVFiles(MeanVector, outputdirectory + "/Mean_PERF.csv");
  WriteCSVFiles(TransformedDataMatrix, outputdirectory + "/PCA_Data.csv");

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
