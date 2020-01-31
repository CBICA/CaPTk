#include "PseudoProgressionEstimator.h"
#include "fMainWindow.h"
#include "cbicaStatistics.h"
#include "CaPTkEnums.h"
#include "vtkDoubleArray.h"

typedef itk::Image< float, 3 > ImageType;

PseudoProgressionEstimator::~PseudoProgressionEstimator()
{
  //delete mNiftiLocalPtr;
  //delete mOutputLocalPtr;
  //delete mFeatureReductionLocalPtr;
  //delete mFeatureScalingLocalPtr;
  //delete mFeatureExtractionLocalPtr;
}

ImageTypeFloat3D::Pointer PseudoProgressionEstimator::RescaleImageIntensity(ImageTypeFloat3D::Pointer image)
{
  typedef ImageTypeFloat3D ImageType;
  typedef itk::RescaleIntensityImageFilter< ImageType, ImageType > RescaleFilterType;
  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput(image);
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(255);
  rescaleFilter->Update();
  ImageType::Pointer outputimage = rescaleFilter->GetOutput();
  return outputimage;
}

bool PseudoProgressionEstimator::TrainNewModelOnGivenData(const std::vector<std::map<CAPTK::ImageModalityType, std::string>> qualifiedsubjects, const std::string &outputdirectory, bool useConventionalData, bool useDTIData, bool usePerfData, bool useDistData)
{
  //extraction of features and target labels
  std::vector<double> traininglabels;
  VariableSizeMatrixType TrainingData = LoadPseudoProgressionTrainingData(qualifiedsubjects, traininglabels, outputdirectory);

  WriteCSVFiles(TrainingData, outputdirectory + "/combinedfeatures-captk-afterfixed.csv");
  WriteCSVFiles(traininglabels, outputdirectory + "/labels.csv");

  std::cout << std::endl << "Building model....." << std::endl;

  //scaling the input feature set
  VariableSizeMatrixType scaledFeatureSet;
  scaledFeatureSet.SetSize(qualifiedsubjects.size(), TrainingData.Cols());
  VariableLengthVectorType meanVector;
  VariableLengthVectorType stdVector;
  mFeatureScalingLocalPtr.ScaleGivenTrainingFeatures(TrainingData, scaledFeatureSet, meanVector, stdVector);

  //remove the nan values
  for (unsigned int index1 = 0; index1 < scaledFeatureSet.Rows(); index1++)
  {
    for (unsigned int index2 = 0; index2 < scaledFeatureSet.Cols(); index2++)
    {
      if (std::isnan(scaledFeatureSet[index1][index2]))
        scaledFeatureSet[index1][index2] = 0;
    }
  }
  for (unsigned int index1 = 0; index1 < meanVector.Size(); index1++)
  {
    if (std::isnan(meanVector[index1]))
      meanVector[index1] = 0;
    if (std::isnan(stdVector[index1]))
      stdVector[index1] = 0;
  }
  //writing the parameters of scaling in the output directory
  try
  {
    WriteCSVFiles(meanVector, outputdirectory + "/PSU_ZScore_Mean.csv");
    WriteCSVFiles(stdVector, outputdirectory + "/PSU_ZScore_Std.csv");
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Error in writing output files to the output directory = " + outputdirectory + "Error code : " + std::string(e1.what()));
    return false;
  }

  VariableSizeMatrixType PseudoModelFeatures;
  VariableSizeMatrixType RecurrenceModelFeatures;
  mFeatureExtractionLocalPtr.FormulatePseudoprogressionTrainingData(scaledFeatureSet, traininglabels, PseudoModelFeatures, RecurrenceModelFeatures);


  //  //---------------------------------------------------------------------------
  VariableSizeMatrixType PseudoModelSelectedFeatures = SelectModelFeatures(PseudoModelFeatures);
  VariableSizeMatrixType RecurrenceModelSelectedFeatures = SelectModelFeatures(RecurrenceModelFeatures);

  //myfile.open(outputdirectory + "/PseudoFeatures.csv");
  //for (unsigned int index1 = 0; index1 < PseudoModelSelectedFeatures.Rows(); index1++)
  //{
  //  for (unsigned int index2 = 0; index2 < PseudoModelSelectedFeatures.Cols(); index2++)
  //  {
  //    if (index2 == 0)
  //      myfile << std::to_string(PseudoModelSelectedFeatures[index1][index2]);
  //    else
  //      myfile << "," << std::to_string(PseudoModelSelectedFeatures[index1][index2]);
  //  }
  //  myfile << "\n";
  //}
  //myfile.close();

  //myfile.open(outputdirectory + "/RecurrenceFetures.csv");
  //for (unsigned int index1 = 0; index1 < RecurrenceModelSelectedFeatures.Rows(); index1++)
  //{
  //  for (unsigned int index2 = 0; index2 < RecurrenceModelSelectedFeatures.Cols(); index2++)
  //  {
  //    if (index2 == 0)
  //      myfile << std::to_string(RecurrenceModelSelectedFeatures[index1][index2]);
  //    else
  //      myfile << "," << std::to_string(RecurrenceModelSelectedFeatures[index1][index2]);
  //  }
  //  myfile << "\n";
  //}
  //myfile.close();
  try
  {
    trainOpenCVSVM(PseudoModelSelectedFeatures, outputdirectory + "/" + mPseudoTrainedFile, false, CAPTK::ApplicationCallingSVM::Survival);
    trainOpenCVSVM(RecurrenceModelSelectedFeatures, outputdirectory + "/" + mRecurrenceTrainedFile, false, CAPTK::ApplicationCallingSVM::Survival);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Training on the given subjects failed. Error code : " + std::string(e1.what()));
    return false;
  }
  std::cout << std::endl << "Model saved to the output directory." << std::endl;

  return true;



  //MatrixType sdata;
  //sdata.set_size(806, 15);
  //for (unsigned int i = 0; i < ScaledTrainingData.Rows(); i++)
  //   for (unsigned int j = 0; j < ScaledTrainingData.Cols(); j++)
  //		sdata(i, j) = ScaledTrainingData[i][j];
  //typedef itk::CSVNumericObjectFileWriter<double, 806, 15> WriterTypeS;
  //WriterTypeS::Pointer writerS = WriterTypeS::New();
  //writerS->SetFileName("sData.csv");
  //writerS->SetInput(&sdata);
  //writerS->Write();

  //typedef vnl_matrix<double> MatrixType;
  //MatrixType data;
  //data.set_size(316, 14);
  //for (unsigned int i = 0; i < data.rows(); i++)
  //{
  //  for (unsigned int j = 0; j < data.cols(); j++)
  //  {
  //    if ((i < ResampledTrainingData.Rows()) && (j < ResampledTrainingData.Cols()))
  //    {
  //      data(i, j) = ResampledTrainingData[i][j];
  //    }
  //    else
  //    {
  //      data(i, j) = 0;
  //    }
  //  }
  //}

  //typedef itk::CSVNumericObjectFileWriter<double, 316, 6> WriterTypeR;
  //WriterTypeR::Pointer writerR = WriterTypeR::New();
  //writerR->SetFileName(outputdirectory + "rData.csv");
  //writerR->SetInput(&data);
  //writerR->Write();

  //FILE *t;
  //t = fopen("TrainingData.txt", "w");

  //for (int i = 0; i < ResampledTrainingData.Rows(); i++)
  //{
  //	fprintf(t, "%f ", ResampledTrainingData[i][ResampledTrainingData.Cols() - 1]);
  //	for (int j = 0; j < ResampledTrainingData.Cols() - 1; j++)
  //		fprintf(t, "%d:%lf ", j + 1, ResampledTrainingData[i][j]);
  //	fprintf(t, "\n");
  //}
  //fclose(t);
  //mOutputLocalPtr.SetOutputDirectoryPath(outputdirectory);
  //int size = GetFeatureVectorSize(useConventionalData, useDTIData, usePerfData, useDistData);
  //try
  //{
  //  mOutputLocalPtr.SaveModelResults(ScaledTrainingData, mFeatureScalingLocalPtr.GetMeanVector(), mFeatureScalingLocalPtr.GetStdVector(), perfMeanVector, mFeatureReductionLocalPtr.GetPCATransformationMatrix(), useConventionalData, useDTIData, usePerfData, useDistData, size);
  //}
  //catch (const std::exception& e1)
  //{
  //  logger.WriteError("Error in writing model files to the output directory: " + outputdirectory + ". Error code : " + std::string(e1.what()));
  //  return false;
  //}

  //try
  //{
  //  std::cout << "Building SVM model." << std::endl;
  //  trainOpenCVSVM(ResampledTrainingData, outputdirectory + "/" + mTrainedModelNameXML, true, Recurrence);
  //}
  //catch (const std::exception& e1)
  //{
  //  logger.WriteError("Training on the subjects failed. Error code : " + std::string(e1.what()));
  //  return false;
  //}
  //mFeatureReductionLocalPtr.ResetParameters();
  //mFeatureScalingLocalPtr.ResetParameters();
}

int PseudoProgressionEstimator::GetFeatureVectorSize(bool &useConventionalData, bool &useDTIData, bool &usePerfData, bool &useDistData)
{
  int size = 0;
  if (useConventionalData)
    size = size + 4;
  if (useDTIData)
    size = size + 4;
  if (usePerfData)
    size = size + 5;
  if (useDistData)
    size = size + 1;
  return size;
}

bool PseudoProgressionEstimator::PseudoProgressionEstimateOnExistingModel(std::vector<std::map<CAPTK::ImageModalityType, std::string>> qualifiedsubjects,
  const std::string &modeldirectory,
  const std::string &inputdirectory,
  const std::string &outputdirectory,
  bool useConventionalrData,
  bool useDTIData, bool usePerfData, bool useDistData)
{
  typedef itk::CSVArray2DFileReader<double> ReaderType;
  VectorDouble results;
  CSVFileReaderType::Pointer reader = CSVFileReaderType::New();

  std::vector<double> traininglabels;
  VariableSizeMatrixType TrainingData = LoadPseudoProgressionTestingData(qualifiedsubjects, traininglabels, outputdirectory, modeldirectory);
  WriteCSVFiles(TrainingData, outputdirectory + "/testingfeatures.csv");

  MatrixType meanMatrix;
  VariableLengthVectorType mean;
  VariableLengthVectorType stddevition;
  try
  {
    reader->SetFileName(modeldirectory + "/PSU_ZScore_Mean.csv");
    reader->SetFieldDelimiterCharacter(',');
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();
    reader->Parse();
    meanMatrix = reader->GetArray2DDataObject()->GetMatrix();

    mean.SetSize(meanMatrix.size());
    for (unsigned int i = 0; i < meanMatrix.size(); i++)
      mean[i] = meanMatrix(0, i);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Error in reading the file: " + modeldirectory + "/PSU_ZScore_Mean.csv. Error code : " + std::string(e1.what()));
    //return results;
  }
  MatrixType stdMatrix;
  try
  {
    reader->SetFileName(modeldirectory + "/PSU_ZScore_Std.csv");
    reader->SetFieldDelimiterCharacter(',');
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();
    reader->Parse();
    stdMatrix = reader->GetArray2DDataObject()->GetMatrix();

    stddevition.SetSize(stdMatrix.size());
    for (unsigned int i = 0; i < stdMatrix.size(); i++)
      stddevition[i] = stdMatrix(0, i);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Error in reading the file: " + modeldirectory + "/PSU_ZScore_Std.csv. Error code : " + std::string(e1.what()));
    //return results;
  }
  std::cout << "parameters read." << std::endl;
  VariableSizeMatrixType ScaledTestingData = mFeatureScalingLocalPtr.ScaleGivenTestingFeatures(TrainingData, mean, stddevition);

  //remove the nan values
  for (unsigned int index1 = 0; index1 < ScaledTestingData.Rows(); index1++)
  {
    for (unsigned int index2 = 0; index2 < ScaledTestingData.Cols(); index2++)
    {
      if (std::isnan(ScaledTestingData[index1][index2]))
        ScaledTestingData[index1][index2] = 0;
    }
  }

  WriteCSVFiles(ScaledTestingData, outputdirectory + "/scaledtestingfeatures.csv");

  std::cout << "scaling done." << std::endl;
  VariableSizeMatrixType ScaledFeatureSetAfterAddingLabel;
  ScaledFeatureSetAfterAddingLabel.SetSize(ScaledTestingData.Rows(), ScaledTestingData.Cols() + 1);
  for (unsigned int i = 0; i < ScaledTestingData.Rows(); i++)
  {
    unsigned int j = 0;
    for (j = 0; j < ScaledTestingData.Cols(); j++)
      ScaledFeatureSetAfterAddingLabel(i, j) = ScaledTestingData(i, j);
    ScaledFeatureSetAfterAddingLabel(i, j) = 0;
  }
  //feature selection process for test data
  VariableLengthVectorType psuSelectedFeatures;
  VariableLengthVectorType recSelectedFeatures;
  MatrixType dataMatrix;
  try
  {
    reader->SetFileName(modeldirectory + "/PSU_SelectedFeatures.csv");
    reader->SetFieldDelimiterCharacter(',');
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();
    reader->Parse();
    dataMatrix = reader->GetArray2DDataObject()->GetMatrix();

    psuSelectedFeatures.SetSize(dataMatrix.rows());
    for (unsigned int i = 0; i < dataMatrix.rows(); i++)
      psuSelectedFeatures[i] = dataMatrix(i, 1);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Error in reading the file: " + modeldirectory + "/PSU_SlectedFeatures.csv. Error code : " + std::string(e1.what()));
    //return results;
  }

  std::cout << "PSU features all loaded." << std::endl;
  try
  {
    reader->SetFileName(modeldirectory + "/REC_SelectedFeatures.csv");
    reader->SetFieldDelimiterCharacter(',');
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();
    reader->Parse();
    dataMatrix = reader->GetArray2DDataObject()->GetMatrix();

    recSelectedFeatures.SetSize(dataMatrix.rows());
    for (unsigned int i = 0; i < dataMatrix.rows(); i++)
      recSelectedFeatures[i] = dataMatrix(i, 1);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Error in reading the file: " + modeldirectory + "/REC_SelectedFeatures.csv. Error code : " + std::string(e1.what()));
    //return results;
  }

  std::cout << "Selected features loaded." << std::endl;

  VariableSizeMatrixType PseudoModelSelectedFeatures = GetModelSelectedFeatures(ScaledFeatureSetAfterAddingLabel, psuSelectedFeatures);
  VariableSizeMatrixType RecurrenceModelSelectedFeatures = GetModelSelectedFeatures(ScaledFeatureSetAfterAddingLabel, recSelectedFeatures);

  WriteCSVFiles(PseudoModelSelectedFeatures, outputdirectory + "/PSU_SelectedTestFeatures.csv");
  WriteCSVFiles(RecurrenceModelSelectedFeatures, outputdirectory + "/REC_SelectedTestFeatures.csv");
  //  std::cout << "selected features done: size:" << PseudoModelSelectedFeatures.Rows() << " columns: " << PseudoModelSelectedFeatures.Cols() << std::endl;
  try
  {
    std::ofstream myfile;
    myfile.open(outputdirectory + "/results.csv");
    myfile << "SubjectName,Score (Pseudo), Score (Recurrence)\n";
    if (cbica::fileExists(modeldirectory + "/PSU_SVM_Model.xml") == true && cbica::fileExists(modeldirectory + "/REC_SVM_Model.xml") == true)
    {
      VectorDouble result_6;
      VectorDouble result_18;
      result_6 = testOpenCVSVM(PseudoModelSelectedFeatures, modeldirectory + "/PSU_SVM_Model.xml");
      result_18 = testOpenCVSVM(RecurrenceModelSelectedFeatures, modeldirectory + "/REC_SVM_Model.xml");
      for (size_t i = 0; i < result_6.size(); i++)
      {
        std::map<CAPTK::ImageModalityType, std::string> currentsubject = qualifiedsubjects[i];
        myfile << static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]) + "," + std::to_string(result_6[i]) + "," + std::to_string(result_18[i]) + "\n";
      }
    }
    myfile.close();
  }
  catch (itk::ExceptionObject & excp)
  {
    logger.WriteError("Error caught during testing: " + std::string(excp.GetDescription()));
    //return results;
  }
  //  return results;


  //check for the presence of model file
  //if (!cbica::fileExists(modeldirectory + "/" + mTrainedModelNameXML))
  //{
  //  if (!cbica::fileExists(modeldirectory + "/" + mTrainedModelNameCSV))
  //  {
  //    mLastEncounteredError = "SVM model file is not present in the directory: " + modeldirectory;
  //    return;
  //  }
  //}
  ////check for the file having number of modalities
  //if (!cbica::fileExists(modeldirectory + "/Recurrence_SVM_Modalities.csv"))
  //{
  //  mLastEncounteredError = "File having record of modalities is not present in the directory: " + modeldirectory;
  //  return;
  //}
  //else
  //{
  //  vnl_matrix<double> modalities = mOutputLocalPtr.ReadNumberOfModalities(modeldirectory + "/Recurrence_SVM_Modalities.csv");
  //  std::string message = "Model supports following modalities:\n";
  //  if (modalities[0][0] == 1)
  //    message += "Conventional\n";
  //  if (modalities[0][1] == 1)
  //    message += "DTI\n";
  //  if (modalities[0][2] == 1)
  //    message += "Perfusion\n";
  //  if (modalities[0][3] == 1)
  //    message += "Distance trnaform";
  //  if ((useConventionalrData == true && modalities[0][0] == 0) || (useDTIData == true && modalities[0][1] == 0) || (usePerfData == true && modalities[0][2] == 0) ||
  //    (useDistData == true && modalities[0][3] == 0))
  //  {
  //    mLastEncounteredError = message;
  //    return;
  //  }
  //}
  //check for the presence of z-score record
  //if (!cbica::fileExists(modeldirectory + "/Recurrence_ZScore_Mean.csv") || !cbica::fileExists(modeldirectory + "/Recurrence_ZScore_Std.csv"))
  //{
  //  mLastEncounteredError = "Z-score record is not present in the directory: " + modeldirectory;
  //  return;
  //}

  ////check for the presence of PCA related record
  //if (usePerfData)
  //{
  //  if (!cbica::fileExists(modeldirectory + "/Recurrence_COEF.csv") || !cbica::fileExists(modeldirectory + "/Recurrence_MR.csv"))
  //  {
  //    mLastEncounteredError = "PCA parameters are not present in the directory: " + modeldirectory;
  //    return;
  //  }
  //  else
  //  {
  //try
  //{
  //  mOutputLocalPtr.ReadModelParameters(modeldirectory + "/Recurrence_ZScore_Mean.csv", modeldirectory + "/Recurrence_ZScore_Std.csv", modeldirectory + "/Recurrence_COEF.csv", modeldirectory + "/Recurrence_MR.csv", mean, stds, pca_coefficients, pca_mean);
  //  mFeatureReductionLocalPtr.SetParameters(pca_coefficients, pca_mean);
  //  mFeatureScalingLocalPtr.SetParameters(mean, stds);
  //  //  }
  //  //}
  //  //else
  //  //{
  //  //  mOutputLocalPtr.ReadModelParameters(modeldirectory + "/Recurrence_ZScore_Mean.csv", modeldirectory + "/Recurrence_ZScore_Std.csv", mean, stds);
  //  //  mFeatureScalingLocalPtr.SetParameters(mean, stds);
  //  //}
  //}
  //catch (const std::exception& e1)
  //{
  //  logger.WriteError("Error in reading ZScore and PCA parameters from directory: " + modeldirectory + ". Error code : " + std::string(e1.what()));
  //  return false;
  //}

  //std::vector<std::string> subjectNames = cbica::subdirectoryInDirectory(modeldirectory + "/Subjects3");
  //for (int sid = 0; sid < subjectNames.size(); sid++)
  //{
  //	std::string perfPath;
  //	std::string radPath;
  //	std::string axPath;
  //	std::string trPath;
  //	std::string faPath;
  //	std::string t1Path;
  //	std::string t1cePath;
  //	std::string t2Path;
  //	std::string t2flairPath;
  //	std::string glistrPath;

  //	//std::string t2flairPath = root_directory + "/Subjects/" + subjectNames[i] + "/Flair/";
  //	//std::string t1cePath = root_directory + "/Subjects/" + subjectNames[i] + "/T1ce/";
  //	//std::string t1Path = root_directory + "/Subjects/" + subjectNames[i] + "/T1/";
  //	//std::string t2Path = root_directory + "/Subjects/" + subjectNames[i] + "/T2/";
  //	//std::string perfPath = root_directory + "/Subjects/" + subjectNames[i] + "/Perfusion/";
  //	//std::string dtiPath = root_directory + "/Subjects/" + subjectNames[i] + "/DTI/";
  //	//std::string glistrPath = root_directory + "/Subjects/" + subjectNames[i] + "/Glistr/";

  //	//std::vector<std::string> filenames = cbica::filesInDirectory(t2flairPath);
  //	std::vector<std::string> filenames;
  //	filenames.push_back(t2flairPath);
  //	std::string extension = ".dcm";
  //	if (cbica::getFilenameExtension(filenames[0]) == extension)
  //		imagetype = DICOM;
  //	else
  //		imagetype = NIfTI;

  //for (unsigned int sid = 0; sid < qualifiedsubjects.size(); sid++)
  //{
  //  std::map< ImageModalityType, std::string > currentsubject = qualifiedsubjects[sid];
  //  VectorVectorDouble perfusionIntensities;
  //  VectorVectorDouble otherIntensities;
  //  VectorDouble distanceIntensities;
  //  std::vector<ImageType::IndexType> testindices;


  //  typedef ImageTypeFloat3D ImageType;
  //  typedef ImageTypeFloat4D PerfusionImageType;
  //  ImageType::Pointer T1CEImagePointer;
  //  ImageType::Pointer T2FlairImagePointer;
  //  ImageType::Pointer T1ImagePointer;
  //  ImageType::Pointer T2ImagePointer;
  //  ImageType::Pointer AXImagePointer;
  //  ImageType::Pointer RADImagePointer;
  //  ImageType::Pointer FAImagePointer;
  //  ImageType::Pointer TRImagePointer;
  //  ImageType::Pointer LabelImagePointer;
  //  PerfusionImageType::Pointer perfImagePointer;
  //  ImageType::Pointer dilatedEdema;
  //  try
  //  {
  //    LabelImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_SEG]));
  //    if (usePerfData)
  //      perfImagePointer = mNiftiLocalPtr.Read4DNiftiImage(static_cast<std::string>(currentsubject[IMAGE_TYPE_PERFUSION]));
  //    if (useConventionalrData)
  //    {
  //      T1CEImagePointer = RescaleImageIntensity(mNiftiLocalPtr.ReadNiftiImage(static_cast<std::string>(currentsubject[IMAGE_TYPE_T1CE])));
  //      T2FlairImagePointer = RescaleImageIntensity(mNiftiLocalPtr.ReadNiftiImage(static_cast<std::string>(currentsubject[IMAGE_TYPE_T2FLAIR])));
  //      T1ImagePointer = RescaleImageIntensity(mNiftiLocalPtr.ReadNiftiImage(static_cast<std::string>(currentsubject[IMAGE_TYPE_T1])));
  //      T2ImagePointer = RescaleImageIntensity(mNiftiLocalPtr.ReadNiftiImage(static_cast<std::string>(currentsubject[IMAGE_TYPE_T2])));
  //    }
  //    if (useDTIData)
  //    {
  //      AXImagePointer = RescaleImageIntensity(mNiftiLocalPtr.ReadNiftiImage(static_cast<std::string>(currentsubject[IMAGE_TYPE_AX])));
  //      RADImagePointer = RescaleImageIntensity(mNiftiLocalPtr.ReadNiftiImage(static_cast<std::string>(currentsubject[IMAGE_TYPE_RAD])));
  //      FAImagePointer = RescaleImageIntensity(mNiftiLocalPtr.ReadNiftiImage(static_cast<std::string>(currentsubject[IMAGE_TYPE_FA])));
  //      TRImagePointer = RescaleImageIntensity(mNiftiLocalPtr.ReadNiftiImage(static_cast<std::string>(currentsubject[IMAGE_TYPE_TR])));
  //    }
  //    //cbica::Logging(loggerFile, "Image loading finished.");
  //    //typedef ImageTypeFloat3D OutputImageType;
  //    //OutputImageType::Pointer dilatedEdema;
  //    //typedef itk::BinaryBallStructuringElement<ImageType::PixelType, 3> StructuringElementType;
  //    //StructuringElementType structuringElement;
  //    //structuringElement.SetRadius(1);
  //    //structuringElement.CreateStructuringElement();
  //    //typedef itk::BinaryDilateImageFilter <ImageType, ImageType, StructuringElementType> BinaryDilateImageFilterType;
  //    //BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
  //    //dilateFilter->SetInput(LabelImagePointer);
  //    //dilateFilter->SetKernel(structuringElement);
  //    //dilateFilter->SetDilateValue(GLISTR_OUTPUT_LABELS::EDEMA);
  //    //dilateFilter->Update();
  //    //dilatedEdema = dilateFilter->GetOutput();


  //    //auto img = convertVtkToItk<float, 3>(mSlicerManagers[0]->mMask);
  //    std::vector<double> labels;
  //    labels.push_back(GLISTR_OUTPUT_LABELS::EDEMA);
  //    dilatedEdema = GetImageWithLabels<ImageType>(labels, LabelImagePointer);

  //    imagetype = NIfTI;
  //    testindices = mNiftiLocalPtr.LoadTestData(T1CEImagePointer, T2FlairImagePointer, T1ImagePointer, T2ImagePointer, perfImagePointer, AXImagePointer, FAImagePointer, RADImagePointer, TRImagePointer, LabelImagePointer, dilatedEdema, perfusionIntensities, otherIntensities, distanceIntensities, imagetype, useConventionalrData, useDTIData, usePerfData, useDistData);
  //    cbica::Logging(loggerFile, "Test data  loading finished.");
  //  }
  //  catch (const std::exception& e1)
  //  {
  //    logger.WriteError("Error in extracting features for patient ID: " + static_cast<std::string>(currentsubject[IMAGE_TYPE_SUDOID]) + ". Error code : " + std::string(e1.what()));
  //    return false;
  //  }

  //  VectorVectorDouble reducedPerfusionFeatures;
  //  if (usePerfData)
  //    reducedPerfusionFeatures = mFeatureReductionLocalPtr.ApplyPCAOnTestData(perfusionIntensities);

  //  int NumberOfPCs = 5;
  //  VectorVectorDouble globaltestintensities;

  //  for (unsigned int k = 0; k < testindices.size(); k++)
  //  {
  //    VectorDouble inten;

  //    if (usePerfData)
  //      for (int j = 0; j < NumberOfPCs; j++)
  //        inten.push_back(reducedPerfusionFeatures[k][j]);

  //    if (useOtherModalities)
  //      for (unsigned int j = 0; j < otherIntensities[0].size(); j++)
  //        inten.push_back(otherIntensities[k][j]);

  //    if (useDistData)
  //      inten.push_back(distanceIntensities[k]);

  //    if (inten.size() > 0)
  //      globaltestintensities.push_back(inten);
  //  }
  //  VariableSizeMatrixType TestingData = mFeatureExtractionLocalPtr.FormulateTestData(globaltestintensities);

  //  //typedef vnl_matrix<double> MatrixType;
  //  //MatrixType data;
  //  //data.set_size(98721, 15);

  //  //for (unsigned int i = 0; i < TestingData.Rows(); i++)
  //  //	for (unsigned int j = 0; j < TestingData.Cols(); j++)
  //  //		data(i, j) = TestingData[i][j];
  //  //typedef itk::CSVNumericObjectFileWriter<double, 98721, 15> WriterType;
  //  //WriterType::Pointer writer = WriterType::New();
  //  //writer->SetFileName("testdata.csv");
  //  //writer->SetInput(&data);
  //  //try
  //  //{
  //  //	writer->Write();
  //  //}
  //  //catch (itk::ExceptionObject & excp)
  //  //{
  //  //    std::cerr << "Error: " << excp.what();
  //  //}


  //  VariableSizeMatrixType ScaledTestingData = mFeatureScalingLocalPtr.ScaleGivenTestingFeatures(TestingData);
  //  //MatrixType sdata;
  //  //sdata.set_size(98721, 15);
  //  //for (int i = 0; i < TestingData.Rows(); i++)
  //  //	for (int j = 0; j < TestingData.Cols(); j++)
  //  //		sdata(i, j) = ScaledTestingData[i][j];
  //  //writer->SetFileName("scaledtestdata.csv");
  //  //writer->SetInput(&sdata);
  //  //try
  //  //{
  //  //	writer->Write();
  //  //}
  //  //catch (itk::ExceptionObject & excp)
  //  //{
  //  //}



  //  ImageType::RegionType region = T1CEImagePointer->GetLargestPossibleRegion();
  //  ImageType::Pointer RecProbabilityMap = ImageType::New();
  //  RecProbabilityMap->SetRegions(region);
  //  RecProbabilityMap->Allocate();
  //  RecProbabilityMap->SetSpacing(T1CEImagePointer->GetSpacing());
  //  RecProbabilityMap->SetOrigin(T1CEImagePointer->GetOrigin());
  //  RecProbabilityMap->SetDirection(T1CEImagePointer->GetDirection());

  //  VectorDouble result_modified;
  //  typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
  //  IteratorType imIt(LabelImagePointer, LabelImagePointer->GetLargestPossibleRegion());
  //  IteratorType RecIt(RecProbabilityMap, RecProbabilityMap->GetLargestPossibleRegion());
  //  imIt.GoToBegin(); RecIt.GoToBegin();

  //  while (!imIt.IsAtEnd())
  //  {
  //    if (!(imIt.Get() == GLISTR_OUTPUT_LABELS::EDEMA))
  //      RecIt.Set(VOXEL_STATUS::OFF);

  //    ++imIt;
  //    ++RecIt;
  //  }
  //  cbica::Logging(loggerFile, "Before testing.");
  //  //try
  //  //{
  //  //  if (cbica::fileExists(modeldirectory + "/" + mTrainedModelNameCSV))
  //  //  {
  //  //    cbica::Logging(loggerFile, "Before testing 1.");
  //  //    VariableLengthVectorType result;
  //  //    result = DistanceFunction(ScaledTestingData, modeldirectory + "/" + mTrainedModelNameCSV, RECURRENCE_MODEL_RHO, RECURRENCE_MODEL_G);
  //  //    for (unsigned int index = 0; index < result.Size(); index++)
  //  //    {
  //  //      RecProbabilityMap->SetPixel(testindices[index], result[index] * 1);
  //  //      result_modified.push_back(result[index]);
  //  //    }

  //  //  }
  //  //  else
  //  //  {
  //  //    cbica::Logging(loggerFile, "Before testing 2.");
  //  //    VectorDouble result;
  //  //    result = testOpenCVSVM(ScaledTestingData, modeldirectory + "/" + mTrainedModelNameXML);
  //  //    for (unsigned int index = 0; index < result.size(); index++)
  //  //      RecProbabilityMap->SetPixel(testindices[index], result[index] * 1);
  //  //    result_modified = result;
  //  //  }
  //  //}
  //  //catch (itk::ExceptionObject & excp)
  //  //{
  //  //  cbica::Logging(loggerFile, "Error caught during testing: " + std::string(excp.GetDescription()));
  //  //  return false;
  //  //}


  //  //---------------------------post-processing------------------------------------------------------------------------------------
  //  try
  //  {
  //    VectorDouble result_revised = RecurrenceMapPostprocessing<ImageType>(result_modified, testindices, RecProbabilityMap, dilatedEdema);
  //    for (unsigned int index = 0; index < result_modified.size(); index++)
  //      RecProbabilityMap->SetPixel(testindices[index], result_revised[index] * 1);

  //    //averaging filter
  //    typedef itk::MeanImageFilter<ImageType, ImageType > FilterType;
  //    FilterType::Pointer meanFilter = FilterType::New();
  //    FilterType::InputSizeType radius;
  //    radius.Fill(1);
  //    meanFilter->SetRadius(radius);
  //    meanFilter->SetInput(RecProbabilityMap);
  //    ImageType::Pointer RevisedRecurrenceMap = meanFilter->GetOutput();

  //    if (imagetype == NIfTI)
  //      mOutputLocalPtr.WriteRecurrenceOutputInNifti<ImageType>(RevisedRecurrenceMap, outputdirectory + "/" + static_cast<std::string>(currentsubject[IMAGE_TYPE_SUDOID]));
  //    else
  //      mOutputLocalPtr.WriteRecurrenceOutputInNifti<ImageType>(RevisedRecurrenceMap, outputdirectory + "/" + static_cast<std::string>(currentsubject[IMAGE_TYPE_SUDOID]));
  //  }
  //  catch (itk::ExceptionObject & excp)
  //  {
  //    cbica::Logging(loggerFile, "Error caught during post-processing of recurrence map: " + std::string(excp.GetDescription()));
  //    return false;
  //  }
  //  //----------------------------------------------------------------------------------------------------------------------------------

  //}
  //mFeatureReductionLocalPtr.ResetParameters();
  //mFeatureScalingLocalPtr.ResetParameters();
  return true;
}

VariableLengthVectorType PseudoProgressionEstimator::DistanceFunction(const VariableSizeMatrixType &testData, const std::string &filename, const double &rho, const double &bestg)
{
  CSVFileReaderType::Pointer readerMean = CSVFileReaderType::New();
  readerMean->SetFileName(filename);
  readerMean->SetFieldDelimiterCharacter(',');
  readerMean->HasColumnHeadersOff();
  readerMean->HasRowHeadersOff();
  readerMean->Parse();
  MatrixType dataMatrix = readerMean->GetArray2DDataObject()->GetMatrix();

  VariableSizeMatrixType SupportVectors;
  VariableLengthVectorType Coefficients;
  VariableLengthVectorType Distances;

  SupportVectors.SetSize(dataMatrix.rows(), dataMatrix.cols() - 1);
  Coefficients.SetSize(dataMatrix.rows(), 1);
  Distances.SetSize(testData.Rows(), 1);

  for (unsigned int i = 0; i < dataMatrix.rows(); i++)
  {
    unsigned int j = 0;
    for (j = 0; j < dataMatrix.cols() - 1; j++)
      SupportVectors(i, j) = dataMatrix(i, j);
    Coefficients[i] = dataMatrix(i, j);
  }
  VariableLengthVectorType yyy;
  yyy.SetSize(testData.Rows(), 1);


  for (unsigned int patID = 0; patID < testData.Rows(); patID++)
  {
    double distance = 0;
    for (unsigned int svID = 0; svID < SupportVectors.Rows(); svID++)
    {
      double euclideanDistance = 0;
      for (unsigned int iterator = 0; iterator < SupportVectors.Cols(); iterator++)
        euclideanDistance = euclideanDistance + (SupportVectors(svID, iterator) - testData(patID, iterator))*(SupportVectors(svID, iterator) - testData(patID, iterator));
      double result = std::exp(-1 * bestg*euclideanDistance);
      distance = distance + result * Coefficients[svID];
    }
    Distances[patID] = distance - rho;
  }
  return Distances;
}

VariableSizeMatrixType PseudoProgressionEstimator::LoadPseudoProgressionTestingData(const std::vector<std::map<CAPTK::ImageModalityType, std::string>> &testingsubjects, std::vector<double> &testinglabels, std::string outputdirectory, std::string modeldirectory)
{
  VariableSizeMatrixType FeaturesOfAllSubjects;
  FeaturesOfAllSubjects.SetSize(testingsubjects.size(), 1040);

  VariableSizeMatrixType otherFeatures;
  otherFeatures.SetSize(testingsubjects.size(), 810);

  VectorVectorDouble perfusionFeatures;

  VectorVectorDouble T1IntensityHistogram;
  VectorVectorDouble T2IntensityHistogram;
  VectorVectorDouble TCIntensityHistogram;
  VectorVectorDouble T1TCIntensityHistogram;
  VectorVectorDouble T2FLIntensityHistogram;
  VectorVectorDouble FLIntensityHistogram;
  VectorVectorDouble AXIntensityHistogram;
  VectorVectorDouble FAIntensityHistogram;
  VectorVectorDouble RDIntensityHistogram;
  VectorVectorDouble TRIntensityHistogram;
  VectorVectorDouble PHIntensityHistogram;
  VectorVectorDouble PSIntensityHistogram;
  VectorVectorDouble RCIntensityHistogram;
  VectorVectorDouble PCA1IntensityHistogram;
  VectorVectorDouble PCA2IntensityHistogram;
  VectorVectorDouble PCA3IntensityHistogram;
  VectorVectorDouble PCA4IntensityHistogram;
  VectorVectorDouble PCA5IntensityHistogram;
  VectorVectorDouble PCA6IntensityHistogram;
  VectorVectorDouble PCA7IntensityHistogram;
  VectorVectorDouble PCA8IntensityHistogram;
  VectorVectorDouble PCA9IntensityHistogram;
  VectorVectorDouble PCA10IntensityHistogram;

  PerfusionMapType PerfusionDataMap;


  //Extracting perfusion data of all the patients and putting in PerfusionDataMap
  for (unsigned int sid = 0; sid < testingsubjects.size(); sid++)
  {
    std::cout << "Loading Perfusion Image: " << sid << std::endl;
    std::map<CAPTK::ImageModalityType, std::string> currentsubject = testingsubjects[sid];
    ImageType::Pointer LabelImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG]));
    //NiftiDataManager m_obj;
    //auto perfImagePointerNifti = m_obj.Read4DNiftiImage(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION]));
    using ImageTypePerfusion = itk::Image< float, 4 >;
    auto perfImagePointerNifti = cbica::ReadImage< ImageTypePerfusion >(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION]));
    
    std::vector<ImageType::IndexType> indices;

    VariableSizeMatrixType perfusionData = LoadPerfusionData<PerfusionImageType, ImageType>(LabelImagePointer, perfImagePointerNifti, indices);
    PerfusionTupleType new_tuple(indices, perfusionData);
    PerfusionDataMap[sid] = new_tuple;
  }
  VariableSizeMatrixType PCA_PERF;
  VariableSizeMatrixType PCA_T1;
  VariableSizeMatrixType PCA_T1CE;
  VariableSizeMatrixType PCA_T2;
  VariableSizeMatrixType PCA_FL;
  VariableSizeMatrixType PCA_T1T1CE;
  VariableSizeMatrixType PCA_T2FL;
  VariableSizeMatrixType PCA_AX;
  VariableSizeMatrixType PCA_FA;
  VariableSizeMatrixType PCA_RAD;
  VariableSizeMatrixType PCA_TR;
  VariableSizeMatrixType PCA_PH;
  VariableSizeMatrixType PCA_PSR;
  VariableSizeMatrixType PCA_RCBV;
  VariableSizeMatrixType PCA_PC1;
  VariableSizeMatrixType PCA_PC2;
  VariableSizeMatrixType PCA_PC3;
  VariableSizeMatrixType PCA_PC4;
  VariableSizeMatrixType PCA_PC5;
  VariableSizeMatrixType PCA_PC6;
  VariableSizeMatrixType PCA_PC7;
  VariableSizeMatrixType PCA_PC8;
  VariableSizeMatrixType PCA_PC9;
  VariableSizeMatrixType PCA_PC10;
  VariableLengthVectorType Mean_PERF;
  VariableLengthVectorType Mean_T1;
  VariableLengthVectorType Mean_T1CE;
  VariableLengthVectorType Mean_T2;
  VariableLengthVectorType Mean_FL;
  VariableLengthVectorType Mean_T1T1CE;
  VariableLengthVectorType Mean_T2FL;
  VariableLengthVectorType Mean_AX;
  VariableLengthVectorType Mean_FA;
  VariableLengthVectorType Mean_RAD;
  VariableLengthVectorType Mean_TR;
  VariableLengthVectorType Mean_PH;
  VariableLengthVectorType Mean_PSR;
  VariableLengthVectorType Mean_RCBV;
  VariableLengthVectorType Mean_PC1;
  VariableLengthVectorType Mean_PC2;
  VariableLengthVectorType Mean_PC3;
  VariableLengthVectorType Mean_PC4;
  VariableLengthVectorType Mean_PC5;
  VariableLengthVectorType Mean_PC6;
  VariableLengthVectorType Mean_PC7;
  VariableLengthVectorType Mean_PC8;
  VariableLengthVectorType Mean_PC9;
  VariableLengthVectorType Mean_PC10;

  ReadAllTheModelParameters(modeldirectory, PCA_PERF, PCA_T1, PCA_T1CE, PCA_T2, PCA_FL, PCA_T1T1CE, PCA_T2FL,
    PCA_AX, PCA_FA, PCA_RAD, PCA_TR, PCA_PH, PCA_PSR, PCA_RCBV, PCA_PC1, PCA_PC2, PCA_PC3, PCA_PC4,
    PCA_PC5, PCA_PC6, PCA_PC7, PCA_PC8, PCA_PC9, PCA_PC10,
    Mean_PERF, Mean_T1, Mean_T1CE, Mean_T2, Mean_FL, Mean_T1T1CE, Mean_T2FL, Mean_AX, Mean_FA,
    Mean_RAD, Mean_TR, Mean_PH, Mean_PSR, Mean_RCBV, Mean_PC1, Mean_PC2,
    Mean_PC3, Mean_PC4, Mean_PC5, Mean_PC6, Mean_PC7, Mean_PC8, Mean_PC9, Mean_PC10);

  //Apply existing PCA model on the test patient
  //--------------------------------------------------------------------------------------------
  std::cout << "Combining and calculating perfusion PCA on test data: "<< std::endl;
  PerfusionMapType perfFeatures = CombineAndCalculatePerfusionPCAForTestData(PerfusionDataMap, PCA_PERF, Mean_PERF);
  std::vector<std::vector<ImageType::Pointer>> RevisedPerfusionImagesOfAllPatients;

  for (unsigned int sid = 0; sid < testingsubjects.size(); sid++)
  {
    std::cout << "Revising Perfusion Image: " << sid << std::endl;
    std::map<CAPTK::ImageModalityType, std::string> currentsubject = testingsubjects[sid];
    NiftiDataManager m_obj;
    auto perfImagePointerNifti = m_obj.Read4DNiftiImage(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION]));
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
    for (int i = 0; i < 10; i++)
    {
      regionIndex[3] = i;
      ImageTypeFloat4D::RegionType desiredRegion(regionIndex, regionSize);
      auto filter = itk::ExtractImageFilter< ImageTypeFloat4D, ImageTypeFloat3D >::New();
      filter->SetExtractionRegion(desiredRegion);
      filter->SetInput(perfImagePointerNifti);
      filter->SetDirectionCollapseToIdentity();
      filter->Update();
      ImageType::Pointer CurrentTimePoint = filter->GetOutput();

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
      //cbica::WriteImage<ImageType>(CurrentTimePoint, "E:/Projects/PSU/Data/" + std::to_string(sid) + "_" + std::to_string(i) + ".nii.gz");
    }
    RevisedPerfusionImagesOfAllPatients.push_back(OnePatientperfusionImages);
  }

  ////Load remaining data of all the patietns
  ////----------------------------------------

  for (unsigned int sid = 0; sid < testingsubjects.size(); sid++)
  {
    std::cout << "Loading Remianing Features: " << sid << std::endl;
    VectorDouble neuroScores;
    std::map<CAPTK::ImageModalityType, std::string> currentsubject = testingsubjects[sid];

    CSVFileReaderType::Pointer reader = CSVFileReaderType::New();
    MatrixType dataMatrix;
    reader->SetFileName(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_FEATURES]));
    reader->SetFieldDelimiterCharacter(',');
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();
    reader->Parse();
    dataMatrix = reader->GetArray2DDataObject()->GetMatrix();

    //for (unsigned int i = 0; i < dataMatrix.rows(); i++)
    //{
    //  neuroScores.push_back(dataMatrix(i, 0));
    //  neuroScores.push_back(dataMatrix(i, 1));
    //  neuroScores.push_back(dataMatrix(i, 2));
    //  neuroScores.push_back(dataMatrix(i, 3));
    //  neuroScores.push_back(dataMatrix(i, 4));
    //  neuroScores.push_back(dataMatrix(i, 5));
    //  testinglabels.push_back(dataMatrix(i, 6));
    //}

    testinglabels.push_back(0);
    ImageType::Pointer LabelImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG]));

    ImageType::Pointer OriginalT1CEImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1CE]));
    ImageType::Pointer OriginalT2FlairImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR]));
    ImageType::Pointer OriginalT1ImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1]));
    ImageType::Pointer OriginalT2ImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2]));

    ImageType::Pointer OriginalT1T1CEImagePointer = MakeAdditionalModality<ImageType>(OriginalT1ImagePointer, OriginalT1CEImagePointer);
    ImageType::Pointer OriginalT2FLImagePointer = MakeAdditionalModality<ImageType>(OriginalT2ImagePointer, OriginalT2FlairImagePointer);

    ImageType::Pointer T1ImagePointer = RescaleImageIntensity<ImageType>(OriginalT1ImagePointer);
    ImageType::Pointer T1CEImagePointer = RescaleImageIntensity<ImageType>(OriginalT1CEImagePointer);
    ImageType::Pointer T2ImagePointer = RescaleImageIntensity<ImageType>(OriginalT2ImagePointer);
    ImageType::Pointer T2FlairImagePointer = RescaleImageIntensity<ImageType>(OriginalT2FlairImagePointer);
    ImageType::Pointer T1T1CEImagePointer = RescaleImageIntensity<ImageType>(OriginalT1T1CEImagePointer);
    ImageType::Pointer T2FLImagePointer = RescaleImageIntensity<ImageType>(OriginalT2FLImagePointer);

    ImageType::Pointer AXImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_AX])));
    ImageType::Pointer RADImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_RAD])));
    ImageType::Pointer FAImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_FA])));
    ImageType::Pointer TRImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_TR])));

    ImageType::Pointer RCBVImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_RCBV])));
    ImageType::Pointer PHImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PH])));
    ImageType::Pointer PSRImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PSR])));



    typedef std::tuple< VectorDouble, VectorDouble, VectorDouble, VectorDouble, VectorDouble> TupleType;
    typedef std::map<std::string, TupleType> MapType;
    MapType OtherFeaturesInMap;


    //calculate perfusion based features from all the pre-calculated perfusion images
    for (int i = 0; i < 10; i++)
      OtherFeaturesInMap["Z" + std::to_string(i)] = GetAllFeaturesPerImagePerROI<ImageType>(RevisedPerfusionImagesOfAllPatients[sid][i], LabelImagePointer, "PCA_" + std::to_string(i));

    //5 shape features per patient
    VectorDouble ShapeFeatures = GetShapeFeatures<ImageType>(LabelImagePointer);

    OtherFeaturesInMap["C0"] = GetAllFeaturesPerImagePerROI<ImageType>(T1ImagePointer, LabelImagePointer, "T1");
    OtherFeaturesInMap["C1"] = GetAllFeaturesPerImagePerROI<ImageType>(T1CEImagePointer, LabelImagePointer, "TC");
    OtherFeaturesInMap["C2"] = GetAllFeaturesPerImagePerROI<ImageType>(T2ImagePointer, LabelImagePointer, "T2");
    OtherFeaturesInMap["C3"] = GetAllFeaturesPerImagePerROI<ImageType>(T2FlairImagePointer, LabelImagePointer, "FL");
    OtherFeaturesInMap["C4"] = GetAllFeaturesPerImagePerROI<ImageType>(T1T1CEImagePointer, LabelImagePointer, "T1TC");
    OtherFeaturesInMap["C5"] = GetAllFeaturesPerImagePerROI<ImageType>(T2FLImagePointer, LabelImagePointer, "T2FL");


    OtherFeaturesInMap["D0"] = GetAllFeaturesPerImagePerROI<ImageType>(AXImagePointer, LabelImagePointer, "AX");
    OtherFeaturesInMap["D1"] = GetAllFeaturesPerImagePerROI<ImageType>(FAImagePointer, LabelImagePointer, "FA");
    OtherFeaturesInMap["D2"] = GetAllFeaturesPerImagePerROI<ImageType>(RADImagePointer, LabelImagePointer, "RD");
    OtherFeaturesInMap["D3"] = GetAllFeaturesPerImagePerROI<ImageType>(TRImagePointer, LabelImagePointer, "TR");
    OtherFeaturesInMap["P0"] = GetAllFeaturesPerImagePerROI<ImageType>(PHImagePointer, LabelImagePointer, "PH");
    OtherFeaturesInMap["P1"] = GetAllFeaturesPerImagePerROI<ImageType>(PSRImagePointer, LabelImagePointer, "PS");
    OtherFeaturesInMap["P2"] = GetAllFeaturesPerImagePerROI<ImageType>(RCBVImagePointer, LabelImagePointer, "RC");


    int counter = 0;
    //for (int i = 0; i < neuroScores.size(); i++)
    //{
    //  otherFeatures[sid][counter] = neuroScores[i];
    //  counter++;
    //}
    for (int i = 0; i < ShapeFeatures.size(); i++)
    {
      otherFeatures[sid][counter] = ShapeFeatures[i];
      counter++;
    }


    std::cout << "Shape and neuro features calculated." << std::endl;
    //10 histogram, 7 intensity, 8 GLCM, 10 GLRLM
    //22 modalities * (10+7+18) = 770
    //770+5 shape features+6 neuro features = 781
    for (auto const &mapiterator : OtherFeaturesInMap)
    {
      VectorDouble Features = std::get<0>(mapiterator.second);
      for (int j = 0; j < Features.size(); j++)
      {
        otherFeatures[sid][counter] = Features[j];
        counter++;
      }

      Features = std::get<1>(mapiterator.second);
      for (int j = 0; j < Features.size(); j++)
      {
        otherFeatures[sid][counter] = Features[j];
        counter++;
      }

      Features = std::get<2>(mapiterator.second);
      for (int j = 0; j < Features.size(); j++)
      {
        otherFeatures[sid][counter] = Features[j];
        counter++;
      }

      Features = std::get<3>(mapiterator.second);
      for (int j = 0; j < Features.size(); j++)
      {
        otherFeatures[sid][counter] = Features[j];
        counter++;
      }
      std::cout << "Counter Size" << counter << std::endl;
    }
    std::cout << "Basic features copied in the OtherFeatures map." << std::endl;


    T1IntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["C0"]));
    TCIntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["C1"]));
    T2IntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["C2"]));
    FLIntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["C3"]));
    T1TCIntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["C4"]));
    T2FLIntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["C5"]));

    AXIntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["D0"]));
    FAIntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["D1"]));
    RDIntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["D2"]));
    TRIntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["D3"]));

    PHIntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["P0"]));
    PSIntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["P1"]));
    RCIntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["P2"]));

    PCA1IntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["Z0"]));
    PCA2IntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["Z1"]));
    PCA3IntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["Z2"]));
    PCA4IntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["Z3"]));
    PCA5IntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["Z4"]));
    PCA6IntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["Z5"]));
    PCA7IntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["Z6"]));
    PCA8IntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["Z7"]));
    PCA9IntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["Z8"]));
    PCA10IntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["Z9"]));
  }
  std::cout << "PCA features of all the modalities copied in the OtherFeatures map." << std::endl;
  //-----------------------------------------------------------------
  FeatureReductionClass m_featureReduction;
  //22 modalities*10 principal components = 220
  VectorVectorDouble T1ReducedIntensityHistogram = mFeatureReductionLocalPtr.ApplyPCAOnTestDataWithGivenTransformations(T1IntensityHistogram, PCA_T1, Mean_T1);
  VectorVectorDouble TCReducedIntensityHistogram = mFeatureReductionLocalPtr.ApplyPCAOnTestDataWithGivenTransformations(TCIntensityHistogram, PCA_T1CE, Mean_T1CE);
  VectorVectorDouble T1TCReducedIntensityHistogram = mFeatureReductionLocalPtr.ApplyPCAOnTestDataWithGivenTransformations(T1TCIntensityHistogram, PCA_T1T1CE, Mean_T1T1CE);
  VectorVectorDouble T2ReducedIntensityHistogram = mFeatureReductionLocalPtr.ApplyPCAOnTestDataWithGivenTransformations(T2IntensityHistogram, PCA_T2, Mean_T2);
  VectorVectorDouble FLReducedIntensityHistogram = mFeatureReductionLocalPtr.ApplyPCAOnTestDataWithGivenTransformations(FLIntensityHistogram, PCA_FL, Mean_FL);
  VectorVectorDouble T2FLReducedIntensityHistogram = mFeatureReductionLocalPtr.ApplyPCAOnTestDataWithGivenTransformations(T2FLIntensityHistogram, PCA_T2FL, Mean_T2FL);

  VectorVectorDouble AXReducedIntensityHistogram = mFeatureReductionLocalPtr.ApplyPCAOnTestDataWithGivenTransformations(AXIntensityHistogram, PCA_AX, Mean_AX);
  VectorVectorDouble FAReducedIntensityHistogram = mFeatureReductionLocalPtr.ApplyPCAOnTestDataWithGivenTransformations(FAIntensityHistogram, PCA_FA, Mean_FA);
  VectorVectorDouble RADReducedIntensityHistogram = mFeatureReductionLocalPtr.ApplyPCAOnTestDataWithGivenTransformations(RDIntensityHistogram, PCA_RAD, Mean_RAD);
  VectorVectorDouble TRReducedIntensityHistogram = mFeatureReductionLocalPtr.ApplyPCAOnTestDataWithGivenTransformations(TRIntensityHistogram, PCA_TR, Mean_TR);

  VectorVectorDouble PHReducedIntensityHistogram = mFeatureReductionLocalPtr.ApplyPCAOnTestDataWithGivenTransformations(PHIntensityHistogram, PCA_PH, Mean_PH);
  VectorVectorDouble PSRReducedIntensityHistogram = mFeatureReductionLocalPtr.ApplyPCAOnTestDataWithGivenTransformations(PSIntensityHistogram, PCA_PSR, Mean_PSR);
  VectorVectorDouble RCReducedIntensityHistogram = mFeatureReductionLocalPtr.ApplyPCAOnTestDataWithGivenTransformations(RCIntensityHistogram, PCA_RCBV, Mean_RCBV);

  std::cout << PCA1IntensityHistogram.size() << " " << PCA1IntensityHistogram[0].size() << std::endl;
  std::cout << PCA_PC1.Rows() << " " << PCA_PC1.Cols() << std::endl;
  std::cout << Mean_PC1.Size() << std::endl;


  VectorVectorDouble PC1ReducedIntensityHistogram = mFeatureReductionLocalPtr.ApplyPCAOnTestDataWithGivenTransformations(PCA1IntensityHistogram, PCA_PC1, Mean_PC1);
  VectorVectorDouble PC2ReducedIntensityHistogram = mFeatureReductionLocalPtr.ApplyPCAOnTestDataWithGivenTransformations(PCA2IntensityHistogram, PCA_PC2, Mean_PC2);
  VectorVectorDouble PC3ReducedIntensityHistogram = mFeatureReductionLocalPtr.ApplyPCAOnTestDataWithGivenTransformations(PCA3IntensityHistogram, PCA_PC3, Mean_PC3);
  VectorVectorDouble PC4ReducedIntensityHistogram = mFeatureReductionLocalPtr.ApplyPCAOnTestDataWithGivenTransformations(PCA4IntensityHistogram, PCA_PC4, Mean_PC4);
  VectorVectorDouble PC5ReducedIntensityHistogram = mFeatureReductionLocalPtr.ApplyPCAOnTestDataWithGivenTransformations(PCA5IntensityHistogram, PCA_PC5, Mean_PC5);
  VectorVectorDouble PC6ReducedIntensityHistogram = mFeatureReductionLocalPtr.ApplyPCAOnTestDataWithGivenTransformations(PCA6IntensityHistogram, PCA_PC6, Mean_PC6);
  VectorVectorDouble PC7ReducedIntensityHistogram = mFeatureReductionLocalPtr.ApplyPCAOnTestDataWithGivenTransformations(PCA7IntensityHistogram, PCA_PC7, Mean_PC7);
  VectorVectorDouble PC8ReducedIntensityHistogram = mFeatureReductionLocalPtr.ApplyPCAOnTestDataWithGivenTransformations(PCA8IntensityHistogram, PCA_PC8, Mean_PC8);
  VectorVectorDouble PC9ReducedIntensityHistogram = mFeatureReductionLocalPtr.ApplyPCAOnTestDataWithGivenTransformations(PCA9IntensityHistogram, PCA_PC9, Mean_PC9);
  VectorVectorDouble PC10ReducedIntensityHistogram = mFeatureReductionLocalPtr.ApplyPCAOnTestDataWithGivenTransformations(PCA10IntensityHistogram, PCA_PC10, Mean_PC10);
  //
  //
  VectorVectorDouble PC_Features;
  for (int i = 0; i < T1ReducedIntensityHistogram.size(); i++)
  {
    std::cout << "patient number" << i << std::endl;
    VectorDouble OnePatient;
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(T1ReducedIntensityHistogram[i][j]);
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(TCReducedIntensityHistogram[i][j]);
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(T1TCReducedIntensityHistogram[i][j]);
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(T2ReducedIntensityHistogram[i][j]);
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(FLReducedIntensityHistogram[i][j]);
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(T2FLReducedIntensityHistogram[i][j]);

    for (int j = 0; j < 10; j++)
      OnePatient.push_back(AXReducedIntensityHistogram[i][j]);
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(FAReducedIntensityHistogram[i][j]);
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(RADReducedIntensityHistogram[i][j]);
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(TRReducedIntensityHistogram[i][j]);
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(PHReducedIntensityHistogram[i][j]);
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(PSRReducedIntensityHistogram[i][j]);
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(RCReducedIntensityHistogram[i][j]);

    std::cout<<PC1ReducedIntensityHistogram.size() << " " << PC1ReducedIntensityHistogram[0].size() << std::endl;
    std::cout << "One patient size" << OnePatient.size() << std::endl;
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(PC1ReducedIntensityHistogram[i][j]);
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(PC2ReducedIntensityHistogram[i][j]);
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(PC3ReducedIntensityHistogram[i][j]);
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(PC4ReducedIntensityHistogram[i][j]);
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(PC5ReducedIntensityHistogram[i][j]);
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(PC6ReducedIntensityHistogram[i][j]);
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(PC7ReducedIntensityHistogram[i][j]);
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(PC8ReducedIntensityHistogram[i][j]);
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(PC9ReducedIntensityHistogram[i][j]);
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(PC10ReducedIntensityHistogram[i][j]);

    std::cout << "One patient size" << OnePatient.size() << std::endl;

    PC_Features.push_back(OnePatient);
  }




  ////VectorVectorDouble PC_Features = CombineAllThePerfusionFeaures(T1IntensityHistogram,
  ////  TCIntensityHistogram, T1TCIntensityHistogram, T2IntensityHistogram, FLIntensityHistogram, T2FLIntensityHistogram,
  ////  AXIntensityHistogram,
  ////  FAIntensityHistogram, RDIntensityHistogram, TRIntensityHistogram,
  ////  PHIntensityHistogram, PSIntensityHistogram, RCIntensityHistogram,
  ////  PCA1IntensityHistogram, PCA2IntensityHistogram, PCA3IntensityHistogram,
  ////  PCA4IntensityHistogram, PCA5IntensityHistogram, PCA6IntensityHistogram,
  ////  PCA7IntensityHistogram, PCA8IntensityHistogram, PCA9IntensityHistogram,
  ////  PCA10IntensityHistogram);

  //std::cout << "Final PCA calculation features finished." << std::endl;

  for (uint i = 0; i < FeaturesOfAllSubjects.Rows(); i++)
  {
    for (uint j = 0; j < otherFeatures.Cols(); j++)
      FeaturesOfAllSubjects(i, j) = otherFeatures[i][j];
    for (uint j = 0; j < PC_Features[i].size(); j++)
      FeaturesOfAllSubjects(i, j + otherFeatures.Cols()) = PC_Features[i][j];
  }
  std::cout << "FeaturesOfAllSubjects populated." << std::endl;
  return FeaturesOfAllSubjects;
}

VariableSizeMatrixType PseudoProgressionEstimator::LoadPseudoProgressionTrainingData(const std::vector<std::map<CAPTK::ImageModalityType, std::string>> &trainingsubjects, std::vector<double> &traininglabels, std::string outputdirectory)
{
  VariableSizeMatrixType FeaturesOfAllSubjects;
  FeaturesOfAllSubjects.SetSize(trainingsubjects.size(), 1040);

  VariableSizeMatrixType otherFeatures;
  otherFeatures.SetSize(trainingsubjects.size(), 810);

  VectorVectorDouble perfusionFeatures;

  VectorVectorDouble T1IntensityHistogram;
  VectorVectorDouble T2IntensityHistogram;
  VectorVectorDouble TCIntensityHistogram;
  VectorVectorDouble T1TCIntensityHistogram;
  VectorVectorDouble T2FLIntensityHistogram;
  VectorVectorDouble FLIntensityHistogram;
  VectorVectorDouble AXIntensityHistogram;
  VectorVectorDouble FAIntensityHistogram;
  VectorVectorDouble RDIntensityHistogram;
  VectorVectorDouble TRIntensityHistogram;
  VectorVectorDouble PHIntensityHistogram;
  VectorVectorDouble PSIntensityHistogram;
  VectorVectorDouble RCIntensityHistogram;
  VectorVectorDouble PCA1IntensityHistogram;
  VectorVectorDouble PCA2IntensityHistogram;
  VectorVectorDouble PCA3IntensityHistogram;
  VectorVectorDouble PCA4IntensityHistogram;
  VectorVectorDouble PCA5IntensityHistogram;
  VectorVectorDouble PCA6IntensityHistogram;
  VectorVectorDouble PCA7IntensityHistogram;
  VectorVectorDouble PCA8IntensityHistogram;
  VectorVectorDouble PCA9IntensityHistogram;
  VectorVectorDouble PCA10IntensityHistogram;

  PerfusionMapType PerfusionDataMap;


  //Extracting perfusion data of all the patients and putting in PerfusionDataMap
  for (unsigned int sid = 0; sid < trainingsubjects.size(); sid++)
  {
    std::cout << "Loading Perfusion Image: " << sid << std::endl;
    std::map<CAPTK::ImageModalityType, std::string> currentsubject = trainingsubjects[sid];
    ImageType::Pointer LabelImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG]));
    NiftiDataManager m_obj;
    auto perfImagePointerNifti = m_obj.Read4DNiftiImage(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION]));
    std::vector<ImageType::IndexType> indices;

    VariableSizeMatrixType perfusionData = LoadPerfusionData<PerfusionImageType, ImageType>(LabelImagePointer, perfImagePointerNifti, indices);
    PerfusionTupleType new_tuple(indices, perfusionData);
    PerfusionDataMap[sid] = new_tuple;
  }

  //combining perfusion data, calcualting PCA
  //-----------------------------------------
  VariableSizeMatrixType TransformationMatrix;
  VariableLengthVectorType MeanVector;
  PerfusionMapType perfFeatures = CombineAndCalculatePerfusionPCA(PerfusionDataMap, TransformationMatrix, MeanVector);

  std::ofstream myfile;
  myfile.open(outputdirectory + "/PCA_PERF.csv");
  for (unsigned int index1 = 0; index1 < TransformationMatrix.Rows(); index1++)
  {
    for (unsigned int index2 = 0; index2 < TransformationMatrix.Cols(); index2++)
    {
      if (index2 == 0)
        myfile << std::to_string(TransformationMatrix[index1][index2]);
      else
        myfile << "," << std::to_string(TransformationMatrix[index1][index2]);
    }
    myfile << "\n";
  }
  myfile.close();
  myfile.open(outputdirectory + "/Mean_PERF.csv");
  for (unsigned int index1 = 0; index1 < MeanVector.Size(); index1++)
    myfile << std::to_string(MeanVector[index1]) << ",";
  myfile << "\n";
  myfile.close();

  //Putting back in images of respective patients
  //---------------------------------------------

  std::vector<std::vector<ImageType::Pointer>> RevisedPerfusionImagesOfAllPatients;

  for (unsigned int sid = 0; sid < trainingsubjects.size(); sid++)
  {
    std::cout << "Revising Perfusion Image: " << sid << std::endl;
    std::map<CAPTK::ImageModalityType, std::string> currentsubject = trainingsubjects[sid];
    NiftiDataManager m_obj;
    auto perfImagePointerNifti = m_obj.Read4DNiftiImage(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION]));
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
    for (int i = 0; i < 10; i++)
    {
      regionIndex[3] = i;
      ImageTypeFloat4D::RegionType desiredRegion(regionIndex, regionSize);
      auto filter = itk::ExtractImageFilter< ImageTypeFloat4D, ImageTypeFloat3D >::New();
      filter->SetExtractionRegion(desiredRegion);
      filter->SetInput(perfImagePointerNifti);
      filter->SetDirectionCollapseToIdentity();
      filter->Update();
      ImageType::Pointer CurrentTimePoint = filter->GetOutput();

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
      cbica::WriteImage<ImageType>(CurrentTimePoint, outputdirectory + std::to_string(sid) + "_" + std::to_string(i) + ".nii.gz");
    }
    RevisedPerfusionImagesOfAllPatients.push_back(OnePatientperfusionImages);
  }


  //either do all the above mentioned steps or read from the already written files
  //---------------------------------------------------
  //std::vector<std::vector<ImageType::Pointer>> RevisedPerfusionImagesOfAllPatients;
  //for (unsigned int sid = 0; sid < trainingsubjects.size(); sid++)
  //{
  //  std::cout << "Revising Perfusion Image: " << sid << std::endl;
  //  std::map<CAPTK::ImageModalityType, std::string> currentsubject = trainingsubjects[sid];
  //  NiftiDataManager m_obj;
  //  std::vector<ImageType::Pointer> OnePatientperfusionImages;
  //  for (int i = 0; i < 10; i++)
  //  {
  //    ImageTypeFloat3D::Pointer perfImagePointerNifti = m_obj.ReadNiftiImage("E:/SoftwareDevelopmentProjects/PseudoprogressionRelatedMaterial/output" + std::to_string(sid) + "_" + std::to_string(i) + ".nii.gz");
  //    OnePatientperfusionImages.push_back(perfImagePointerNifti);
  //  }
  //  RevisedPerfusionImagesOfAllPatients.push_back(OnePatientperfusionImages);
  //}

  //Load remaining data of all the patietns
  //----------------------------------------

  for (unsigned int sid = 0; sid < trainingsubjects.size(); sid++)
  {
    std::cout << "Loading Remaining Features: " << sid << std::endl;
    std::map<CAPTK::ImageModalityType, std::string> currentsubject = trainingsubjects[sid];

    CSVFileReaderType::Pointer reader = CSVFileReaderType::New();
    MatrixType dataMatrix;
    reader->SetFileName(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_FEATURES]));
    reader->SetFieldDelimiterCharacter(',');
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();
    reader->Parse();
    dataMatrix = reader->GetArray2DDataObject()->GetMatrix();
    traininglabels.push_back(dataMatrix(0, 0));

    ImageType::Pointer LabelImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG]));

    ImageType::Pointer OriginalT1CEImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1CE]));
    ImageType::Pointer OriginalT2FlairImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR]));
    ImageType::Pointer OriginalT1ImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1]));
    ImageType::Pointer OriginalT2ImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2]));

    ImageType::Pointer OriginalT1T1CEImagePointer = MakeAdditionalModality<ImageType>( OriginalT1CEImagePointer, OriginalT1ImagePointer);
    ImageType::Pointer OriginalT2FLImagePointer = MakeAdditionalModality<ImageType>(OriginalT2ImagePointer, OriginalT2FlairImagePointer);

    ImageType::Pointer T1ImagePointer = RescaleImageIntensity<ImageType>(OriginalT1ImagePointer);
    ImageType::Pointer T1CEImagePointer = RescaleImageIntensity<ImageType>(OriginalT1CEImagePointer);
    ImageType::Pointer T2ImagePointer = RescaleImageIntensity<ImageType>(OriginalT2ImagePointer);
    ImageType::Pointer T2FlairImagePointer = RescaleImageIntensity<ImageType>(OriginalT2FlairImagePointer);
    ImageType::Pointer T1T1CEImagePointer = RescaleImageIntensity<ImageType>(OriginalT1T1CEImagePointer);
    ImageType::Pointer T2FLImagePointer = RescaleImageIntensity<ImageType>(OriginalT2FLImagePointer);

    ImageType::Pointer AXImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_AX])));
    ImageType::Pointer RADImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_RAD])));
    ImageType::Pointer FAImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_FA])));
    ImageType::Pointer TRImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_TR])));

    ImageType::Pointer RCBVImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_RCBV])));
    ImageType::Pointer PHImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PH])));
    ImageType::Pointer PSRImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PSR])));

    typedef std::tuple< VectorDouble, VectorDouble, VectorDouble, VectorDouble, VectorDouble> TupleType;
    typedef std::map<std::string, TupleType> MapType;
    MapType OtherFeaturesInMap;


    //calculate perfusion based features from all the pre-calculated perfusion images
    for (int i = 0; i < 10; i++)
      OtherFeaturesInMap["Z" + std::to_string(i)] = GetAllFeaturesPerImagePerROI<ImageType>(RevisedPerfusionImagesOfAllPatients[sid][i], LabelImagePointer, "PCA_" + std::to_string(i));

    //5 shape features per patient
    VectorDouble ShapeFeatures = GetShapeFeatures<ImageType>(LabelImagePointer);

    OtherFeaturesInMap["C0"] = GetAllFeaturesPerImagePerROI<ImageType>(T1ImagePointer, LabelImagePointer, "T1");
    OtherFeaturesInMap["C1"] = GetAllFeaturesPerImagePerROI<ImageType>(T1CEImagePointer, LabelImagePointer, "TC");
    OtherFeaturesInMap["C2"] = GetAllFeaturesPerImagePerROI<ImageType>(T2ImagePointer, LabelImagePointer, "T2");
    OtherFeaturesInMap["C3"] = GetAllFeaturesPerImagePerROI<ImageType>(T2FlairImagePointer, LabelImagePointer, "FL");
    OtherFeaturesInMap["C4"] = GetAllFeaturesPerImagePerROI<ImageType>(T1T1CEImagePointer, LabelImagePointer, "T1TC");
    OtherFeaturesInMap["C5"] = GetAllFeaturesPerImagePerROI<ImageType>(T2FLImagePointer, LabelImagePointer, "T2FL");


    OtherFeaturesInMap["D0"] = GetAllFeaturesPerImagePerROI<ImageType>(AXImagePointer, LabelImagePointer, "AX");
    OtherFeaturesInMap["D1"] = GetAllFeaturesPerImagePerROI<ImageType>(FAImagePointer, LabelImagePointer, "FA");
    OtherFeaturesInMap["D2"] = GetAllFeaturesPerImagePerROI<ImageType>(RADImagePointer, LabelImagePointer, "RD");
    OtherFeaturesInMap["D3"] = GetAllFeaturesPerImagePerROI<ImageType>(TRImagePointer, LabelImagePointer, "TR");
    OtherFeaturesInMap["P0"] = GetAllFeaturesPerImagePerROI<ImageType>(PHImagePointer, LabelImagePointer, "PH");
    OtherFeaturesInMap["P1"] = GetAllFeaturesPerImagePerROI<ImageType>(PSRImagePointer, LabelImagePointer, "PS");
    OtherFeaturesInMap["P2"] = GetAllFeaturesPerImagePerROI<ImageType>(RCBVImagePointer, LabelImagePointer, "RC");


    int counter = 0;
    //for (int i = 0; i < neuroScores.size(); i++)
    //{
    //  otherFeatures[sid][counter] = neuroScores[i];
    //  counter++;
    //}
    for (int i = 0; i < ShapeFeatures.size(); i++)
    {
      otherFeatures[sid][counter] = ShapeFeatures[i];
      counter++;
    }


    std::cout << "Shape and neuro features calculated." << std::endl;
    //10 histogram, 7 intensity, 8 GLCM, 10 GLRLM
    //22 modalities * (10+7+18) = 770
    //770+5 shape features+6 neuro features = 781
    for (auto const &mapiterator : OtherFeaturesInMap)
    {
      VectorDouble Features = std::get<0>(mapiterator.second);
      for (int j = 0; j < Features.size(); j++)
      {
        otherFeatures[sid][counter] = Features[j];
        counter++;
      }

      Features = std::get<1>(mapiterator.second);
      for (int j = 0; j < Features.size(); j++)
      {
        otherFeatures[sid][counter] = Features[j];
        counter++;
      }

      Features = std::get<2>(mapiterator.second);
      for (int j = 0; j < Features.size(); j++)
      {
        otherFeatures[sid][counter] = Features[j];
        counter++;
      }

      Features = std::get<3>(mapiterator.second);
      for (int j = 0; j < Features.size(); j++)
      {
        otherFeatures[sid][counter] = Features[j];
        counter++;
      }
      std::cout << "Counter Size" << counter << std::endl;
    }
    std::cout << "Basic features copied in the OtherFeatures map." << std::endl;


    T1IntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["C0"]));
    TCIntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["C1"]));
    T2IntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["C2"]));
    FLIntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["C3"]));
    T1TCIntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["C4"]));
    T2FLIntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["C5"]));

    AXIntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["D0"]));
    FAIntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["D1"]));
    RDIntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["D2"]));
    TRIntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["D3"]));

    PHIntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["P0"]));
    PSIntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["P1"]));
    RCIntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["P2"]));

    PCA1IntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["Z0"]));
    PCA2IntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["Z1"]));
    PCA3IntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["Z2"]));
    PCA4IntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["Z3"]));
    PCA5IntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["Z4"]));
    PCA6IntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["Z5"]));
    PCA7IntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["Z6"]));
    PCA8IntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["Z7"]));
    PCA9IntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["Z8"]));
    PCA10IntensityHistogram.push_back(std::get<4>(OtherFeaturesInMap["Z9"]));
  }
  std::cout << "PCA features of all the modalities copied in the OtherFeatures map." << std::endl;
  //-----------------------------------------------------------------
  FeatureReductionClass m_featureReduction;
  //22 modalities*10 principal components = 220
  VectorVectorDouble PC_Features = CombineAllThePerfusionFeaures(T1IntensityHistogram,
    TCIntensityHistogram, T1TCIntensityHistogram, T2IntensityHistogram, FLIntensityHistogram, T2FLIntensityHistogram,
    AXIntensityHistogram,
    FAIntensityHistogram, RDIntensityHistogram, TRIntensityHistogram,
    PHIntensityHistogram, PSIntensityHistogram, RCIntensityHistogram,
    PCA1IntensityHistogram, PCA2IntensityHistogram, PCA3IntensityHistogram,
    PCA4IntensityHistogram, PCA5IntensityHistogram, PCA6IntensityHistogram,
    PCA7IntensityHistogram, PCA8IntensityHistogram, PCA9IntensityHistogram,
    PCA10IntensityHistogram, outputdirectory);

  std::cout << "Final PCA calculation features finished." << std::endl;

  for (uint i = 0; i < FeaturesOfAllSubjects.Rows(); i++)
  {
    for (uint j = 0; j < otherFeatures.Cols(); j++)
      FeaturesOfAllSubjects(i, j) = otherFeatures[i][j];
    for (uint j = 0; j < PC_Features[i].size(); j++)
      FeaturesOfAllSubjects(i, j + otherFeatures.Cols()) = PC_Features[i][j];
  }
  std::cout << "FeaturesOfAllSubjects populated." << std::endl;
  return FeaturesOfAllSubjects;
}


PerfusionMapType PseudoProgressionEstimator::CombineAndCalculatePerfusionPCA(PerfusionMapType PerfusionDataMap, VariableSizeMatrixType &TransformationMatrix, VariableLengthVectorType &MeanVector)
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
      for (unsigned int j = 0; j < Features.Cols(); j++)
        oneVector.push_back(Features(i, j));
      CombinedPerfusionFeaturesMap.push_back(oneVector);
    }
  }
  FeatureReductionClass m_featureReduction;
  vtkSmartPointer<vtkTable> ReducedPCAs = m_featureReduction.GetDiscerningPerfusionTimePointsForPSU(CombinedPerfusionFeaturesMap, TransformationMatrix, MeanVector);

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


PerfusionMapType PseudoProgressionEstimator::CombineAndCalculatePerfusionPCAForTestData(PerfusionMapType PerfusionDataMap, VariableSizeMatrixType &TransformationMatrix, VariableLengthVectorType &MeanVector)
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
      for (unsigned int j = 0; j < Features.Cols(); j++)
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


PerfusionMapType PseudoProgressionEstimator::CombinePerfusionDataAndApplyExistingPerfusionModel(PerfusionMapType PerfusionDataMap, VariableSizeMatrixType TransformationMatrix, VariableLengthVectorType MeanVector)
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


VectorDouble PseudoProgressionEstimator::GetIntensityFeatures(std::vector<float> m_nonZeroPixels)
{
  std::vector<double> features;

  cbica::Statistics< typename ImageType::PixelType > m_statistics;
  m_statistics.SetInput(m_nonZeroPixels);

  features.push_back(m_statistics.GetMinimum());
  features.push_back(m_statistics.GetMaximum());
  features.push_back(m_statistics.GetMean());
  features.push_back(m_statistics.GetVariance());
  features.push_back(m_statistics.GetStandardDeviation());
  features.push_back(m_statistics.GetSkewness());
  features.push_back(m_statistics.GetKurtosis());

  return features;
}

VectorDouble PseudoProgressionEstimator::GetHistogramFeatures(std::vector<float> intensities, int m_Bins)
{
  std::vector<double> features;
  const float maxRescaleVal = 255;

  double interval = maxRescaleVal / m_Bins;
  double final_interval = (int)(interval * 100);
  final_interval = (double)final_interval / 100;

  std::vector<double> finalBins;
  std::vector<std::vector<double>> Ranges;
  double current_index = 0;
  for (double i = 0; i < m_Bins; i++)
  {
    std::vector<double> onerange;
    onerange.push_back(current_index);
    current_index = current_index + final_interval;
    onerange.push_back(current_index);
    Ranges.push_back(onerange);

    if (Ranges.size() == m_Bins)
      Ranges[Ranges.size() - 1][1] = maxRescaleVal;
  }
  //toadd the last one
  for (unsigned int j = 0; j < Ranges.size(); j++)
  {
    std::vector<double> onerange = Ranges[j];
    int counter = 0;
    for (unsigned int i = 0; i < intensities.size(); i++)
    {
      if (onerange[0] == 0)
      {
        if (intensities[i] >= onerange[0] && intensities[i] <= onerange[1])
          counter = counter + 1;
      }
      else
      {
        if (intensities[i] > onerange[0] && intensities[i] <= onerange[1])
          counter = counter + 1;
      }
    }
    finalBins.push_back(counter);
  }
  for (unsigned int j = 0; j < finalBins.size(); j++)
  {
    finalBins[j] = (finalBins[j] * 100) / intensities.size();
    features.push_back(finalBins[j]);
    //featurevec["Bin_" + std::to_string(j)] = finalBins[j];
    //featurevec["BinEndIntensity_" + std::to_string(j)] = Ranges[j][1];
  }
  return features;
}

VectorDouble PseudoProgressionEstimator::GetHistogramFeaturesWhole(std::vector<float> intensities)
{
  VectorDouble features;
  for (unsigned int i = 0; i < 255; i++)
  {
    int counter = 0;
    for (unsigned int j = 0; j < intensities.size(); j++)
    {
      if (intensities[j] > i && intensities[j] <= i + 1)
        counter++;
    }
    features.push_back(counter);
  }
  return features;
}

VectorVectorDouble PseudoProgressionEstimator::CombineAllThePerfusionFeaures(VectorVectorDouble T1IntensityHistogram,
  VectorVectorDouble TCIntensityHistogram,
  VectorVectorDouble T1TCIntensityHistogram,
  VectorVectorDouble T2IntensityHistogram,
  VectorVectorDouble FLIntensityHistogram,
  VectorVectorDouble T2FLIntensityHistogram,
  VectorVectorDouble AXIntensityHistogram,
  VectorVectorDouble FAIntensityHistogram,
  VectorVectorDouble RDIntensityHistogram,
  VectorVectorDouble TRIntensityHistogram,
  VectorVectorDouble PHIntensityHistogram,
  VectorVectorDouble PSIntensityHistogram,
  VectorVectorDouble RCIntensityHistogram,
  VectorVectorDouble PCA1IntensityHistogram,
  VectorVectorDouble PCA2IntensityHistogram,
  VectorVectorDouble PCA3IntensityHistogram,
  VectorVectorDouble PCA4IntensityHistogram,
  VectorVectorDouble PCA5IntensityHistogram,
  VectorVectorDouble PCA6IntensityHistogram,
  VectorVectorDouble PCA7IntensityHistogram,
  VectorVectorDouble PCA8IntensityHistogram,
  VectorVectorDouble PCA9IntensityHistogram,
  VectorVectorDouble PCA10IntensityHistogram,
  std::string outputdirectory)
{
  //writing of all the modalities perfusion data finished
  //WriteCSVFiles(T1IntensityHistogram, outputdirectory+ "/t1.csv");
  //WriteCSVFiles(TCIntensityHistogram, outputdirectory+ "/t1ce.csv");
  //WriteCSVFiles(T2IntensityHistogram, outputdirectory + "/t2.csv");
  //WriteCSVFiles(FLIntensityHistogram, outputdirectory + "/flair.csv");
  //WriteCSVFiles(T1TCIntensityHistogram, outputdirectory + "/t1t1ce.csv");
  //WriteCSVFiles(T2FLIntensityHistogram, outputdirectory + "/t2flair.csv");
  //WriteCSVFiles(AXIntensityHistogram, outputdirectory + "/AX.csv");
  //WriteCSVFiles(FAIntensityHistogram, outputdirectory + "/FA.csv");
  //WriteCSVFiles(RDIntensityHistogram, outputdirectory + "/RAD.csv");
  //WriteCSVFiles(TRIntensityHistogram, outputdirectory + "/TR.csv");
  //WriteCSVFiles(PHIntensityHistogram, outputdirectory + "/PH.csv");
  //WriteCSVFiles(PSIntensityHistogram, outputdirectory + "/PSR.csv");
  //WriteCSVFiles(RCIntensityHistogram, outputdirectory + "/RCBV.csv");
  //WriteCSVFiles(PCA1IntensityHistogram, outputdirectory + "/PCA1.csv");
  //WriteCSVFiles(PCA2IntensityHistogram, outputdirectory +"/PCA2.csv");
  //WriteCSVFiles(PCA3IntensityHistogram, outputdirectory +"/PCA3.csv");
  //WriteCSVFiles(PCA4IntensityHistogram, outputdirectory +"/PCA4.csv");
  //WriteCSVFiles(PCA5IntensityHistogram, outputdirectory +"/PCA5.csv");
  //WriteCSVFiles(PCA6IntensityHistogram, outputdirectory +"/PCA6.csv");
  //WriteCSVFiles(PCA7IntensityHistogram, outputdirectory +"/PCA7.csv");
  //WriteCSVFiles(PCA8IntensityHistogram, outputdirectory +"/PCA8.csv");
  //WriteCSVFiles(PCA9IntensityHistogram, outputdirectory +"/PCA9.csv");
  //WriteCSVFiles(PCA10IntensityHistogram, outputdirectory +"/PCA10.csv");


  FeatureReductionClass m_featureReduction;
  VariableSizeMatrixType PCA_T1, PCA_T1CE, PCA_T2, PCA_FL, PCA_T2FL, PCA_T1T1CE, PCA_AX, PCA_FA, PCA_RAD, PCA_TR, PCA_PC1, PCA_PC2, PCA_PC3, PCA_PC4, PCA_PC5, PCA_PC6, PCA_PC7, PCA_PC8, PCA_PC9, PCA_PC10, PCA_PH, PCA_PSR, PCA_RCBV;
  VariableLengthVectorType Mean_T1, Mean_T1CE, Mean_T2, Mean_FL, Mean_T2FL, Mean_T1T1CE, Mean_AX, Mean_FA, Mean_RAD, Mean_TR, Mean_PC1, Mean_PC2, Mean_PC3, Mean_PC4, Mean_PC5, Mean_PC6, Mean_PC7, Mean_PC8, Mean_PC9, Mean_PC10, Mean_PH, Mean_PSR, Mean_RCBV;

  int NumberOfFeatures = T1IntensityHistogram[0].size();
  int NumberOfSamples = T1IntensityHistogram.size();
  vtkSmartPointer<vtkTable> T1ReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> TCReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> T2ReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> FLReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> T1TCReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> T2FLReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> AXReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> FAReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> RDReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> TRReducedIntensityHistogram;
  
  vtkSmartPointer<vtkTable> PHReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> PSReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> RCReducedIntensityHistogram;

  vtkSmartPointer<vtkTable> PC1ReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> PC2ReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> PC3ReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> PC4ReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> PC5ReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> PC6ReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> PC7ReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> PC8ReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> PC9ReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> PC10ReducedIntensityHistogram;


  try
  {
    T1ReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(T1IntensityHistogram, PCA_T1, Mean_T1);
    std::cout << "T1" << std::endl;
  }
  catch (const std::exception& e1)
  {
    T1ReducedIntensityHistogram = MakePCAMatrix(NumberOfFeatures, NumberOfSamples); 
    logger.WriteError("Error in writing T1 reduced intensity histogram. Error code : " + std::string(e1.what()));
  }

  try
  {
    TCReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(TCIntensityHistogram, PCA_T1CE, Mean_T1CE);
    std::cout << "TC" << std::endl;
  }
  catch (const std::exception& e1)
  {
    TCReducedIntensityHistogram = MakePCAMatrix(NumberOfFeatures, NumberOfSamples);
    logger.WriteError("Error in writing TC reduced intensity histogram. Error code : " + std::string(e1.what()));
  }

  try
  {
    T2ReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(T2IntensityHistogram, PCA_T2, Mean_T2);
  }
  catch (const std::exception& e1)
  {
    T2ReducedIntensityHistogram = MakePCAMatrix(NumberOfFeatures, NumberOfSamples);
    logger.WriteError("Error in writing T2 reduced intensity histogram. Error code : " + std::string(e1.what()));
  }
  
  try
  {
    T1TCReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(T1TCIntensityHistogram, PCA_T1T1CE, Mean_T1T1CE);
    std::cout << "T1TC" << std::endl;
  }
  catch (const std::exception& e1)
  {
    T1TCReducedIntensityHistogram = MakePCAMatrix(NumberOfFeatures, NumberOfSamples);
    logger.WriteError("Error in writing T1TC reduced intensity histogram. Error code : " + std::string(e1.what()));
  }

  try
  {
  FLReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(FLIntensityHistogram, PCA_FL, Mean_FL);
  std::cout << "FL" << std::endl;
  }
  catch (const std::exception& e1)
  {
    FLReducedIntensityHistogram = MakePCAMatrix(NumberOfFeatures, NumberOfSamples);
    logger.WriteError("Error in writing FL reduced intensity histogram. Error code : " + std::string(e1.what()));
  }


  try
  {
  T2FLReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(T2FLIntensityHistogram, PCA_T2FL, Mean_T2FL);
  std::cout << "T2FL" << std::endl;
  }
  catch (const std::exception& e1)
  {
    T2FLReducedIntensityHistogram = MakePCAMatrix(NumberOfFeatures, NumberOfSamples);
    logger.WriteError("Error in writing T2FL reduced intensity histogram. Error code : " + std::string(e1.what()));
  }

  try
  {
    AXReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(AXIntensityHistogram, PCA_AX, Mean_AX);
    std::cout << "AX" << std::endl;
  }
  catch (const std::exception& e1)
  {
    AXReducedIntensityHistogram = MakePCAMatrix(NumberOfFeatures, NumberOfSamples);
    logger.WriteError("Error in writing AX reduced intensity histogram. Error code : " + std::string(e1.what()));
  }

  try
  {
    FAReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(FAIntensityHistogram, PCA_FA, Mean_FA);
    std::cout << "FA" << std::endl;
  }
  catch (const std::exception& e1)
  {
    FAReducedIntensityHistogram = MakePCAMatrix(NumberOfFeatures, NumberOfSamples);
    logger.WriteError("Error in writing FA reduced intensity histogram. Error code : " + std::string(e1.what()));
  }

  try
  {
  RDReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(RDIntensityHistogram, PCA_RAD, Mean_RAD);
  std::cout << "RD" << std::endl;
  }
  catch (const std::exception& e1)
  {
    RDReducedIntensityHistogram = MakePCAMatrix(NumberOfFeatures, NumberOfSamples);
    logger.WriteError("Error in writing RD reduced intensity histogram. Error code : " + std::string(e1.what()));
  }

  try
  {
    TRReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(TRIntensityHistogram, PCA_TR, Mean_TR);
    std::cout << "TR" << std::endl;
  }
    catch (const std::exception& e1)
    {
      TRReducedIntensityHistogram = MakePCAMatrix(NumberOfFeatures, NumberOfSamples);
      logger.WriteError("Error in writing TR reduced intensity histogram. Error code : " + std::string(e1.what()));
    }


    try
    {
  PHReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(PHIntensityHistogram, PCA_PH, Mean_PH);
  std::cout << "PH" << std::endl;
    }
    catch (const std::exception& e1)
    {
      PHReducedIntensityHistogram = MakePCAMatrix(NumberOfFeatures, NumberOfSamples);
      logger.WriteError("Error in writing PH reduced intensity histogram. Error code : " + std::string(e1.what()));
    }

    try
    {
  PSReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(PSIntensityHistogram, PCA_PSR, Mean_PSR);
  std::cout << "PS" << std::endl;
    }
    catch (const std::exception& e1)
    {
      PSReducedIntensityHistogram = MakePCAMatrix(NumberOfFeatures, NumberOfSamples);
      logger.WriteError("Error in writing PS reduced intensity histogram. Error code : " + std::string(e1.what()));
    }


    try
    {
  RCReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(RCIntensityHistogram, PCA_RCBV, Mean_RCBV);
  std::cout << "RC" << std::endl;
    }
    catch (const std::exception& e1)
    {
      RCReducedIntensityHistogram = MakePCAMatrix(NumberOfFeatures, NumberOfSamples);
      logger.WriteError("Error in writing RC reduced intensity histogram. Error code : " + std::string(e1.what()));
    }

  std::cout << "basic modalities perfusion components extracted" << std::endl;

  try
  {
    PC1ReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(PCA1IntensityHistogram, PCA_PC1, Mean_PC1);
  std::cout << "PC1" << std::endl;
  }
  catch (const std::exception& e1)
  {
    PC1ReducedIntensityHistogram = MakePCAMatrix(NumberOfFeatures, NumberOfSamples);
    logger.WriteError("Error in writing PC1 reduced intensity histogram. Error code : " + std::string(e1.what()));
  }


  try
  {
    PC2ReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(PCA2IntensityHistogram, PCA_PC2, Mean_PC2);
  std::cout << "PC2" << std::endl;
  }
  catch (const std::exception& e1)
  {
    PC2ReducedIntensityHistogram = MakePCAMatrix(NumberOfFeatures, NumberOfSamples);
    logger.WriteError("Error in writing PC2 reduced intensity histogram. Error code : " + std::string(e1.what()));
  }

  try
  {
  PC3ReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(PCA3IntensityHistogram, PCA_PC3, Mean_PC3);
  std::cout << "PC3" << std::endl;
  }
  catch (const std::exception& e1)
  {
    PC3ReducedIntensityHistogram = MakePCAMatrix(NumberOfFeatures, NumberOfSamples);
    logger.WriteError("Error in writing PC3 reduced intensity histogram. Error code : " + std::string(e1.what()));
  }
  try
  {
    PC4ReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(PCA4IntensityHistogram, PCA_PC4, Mean_PC4);
  std::cout << "PC4" << std::endl;
  }
  catch (const std::exception& e1)
  {
    PC4ReducedIntensityHistogram = MakePCAMatrix(NumberOfFeatures, NumberOfSamples);
    logger.WriteError("Error in writing PC4 reduced intensity histogram. Error code : " + std::string(e1.what()));
  }

  try
  {
    PC5ReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(PCA5IntensityHistogram, PCA_PC5, Mean_PC5);
  std::cout << "PC5" << std::endl;
  }
  catch (const std::exception& e1)
  {
    PC5ReducedIntensityHistogram = MakePCAMatrix(NumberOfFeatures, NumberOfSamples);
    logger.WriteError("Error in writing PC5 reduced intensity histogram. Error code : " + std::string(e1.what()));
  }

  try
  {
    PC6ReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(PCA6IntensityHistogram, PCA_PC6, Mean_PC6);
  std::cout << "PC6" << std::endl;
  }
  catch (const std::exception& e1)
  {
    PC6ReducedIntensityHistogram = MakePCAMatrix(NumberOfFeatures, NumberOfSamples);
    logger.WriteError("Error in writing PC6 reduced intensity histogram. Error code : " + std::string(e1.what()));
  }

  try
  {
    PC7ReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(PCA7IntensityHistogram, PCA_PC7, Mean_PC7);
  std::cout << "PC7" << std::endl;
  }
  catch (const std::exception& e1)
  {
    PC7ReducedIntensityHistogram = MakePCAMatrix(NumberOfFeatures, NumberOfSamples);
    logger.WriteError("Error in writing PC7 reduced intensity histogram. Error code : " + std::string(e1.what()));
  }

  try
  {
    PC8ReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(PCA8IntensityHistogram, PCA_PC8, Mean_PC8);
  std::cout << "PC8" << std::endl;
  }
  catch (const std::exception& e1)
  {
    PC8ReducedIntensityHistogram = MakePCAMatrix(NumberOfFeatures, NumberOfSamples);
    logger.WriteError("Error in writing PC8 reduced intensity histogram. Error code : " + std::string(e1.what()));
  }

  try
  {
    PC9ReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(PCA9IntensityHistogram, PCA_PC9, Mean_PC9);
    std::cout << "PC9" << std::endl;
  }
  catch (const std::exception& e1)
  {
    PC9ReducedIntensityHistogram = MakePCAMatrix(NumberOfFeatures, NumberOfSamples);
    logger.WriteError("Error in writing PC9 reduced intensity histogram. Error code : " + std::string(e1.what()));
  }

  try
  {
    PC10ReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(PCA10IntensityHistogram, PCA_PC10, Mean_PC10);
    std::cout << "PC10" << std::endl;
  }
  catch (const std::exception& e1)
  {
    PC10ReducedIntensityHistogram = MakePCAMatrix(NumberOfFeatures, NumberOfSamples);
    logger.WriteError("Error in writing PC10 reduced intensity histogram. Error code : " + std::string(e1.what()));
  }

  VariableSizeMatrixType AllPCAs;
  VariableSizeMatrixType AllMeans;
  AllPCAs.SetSize(23 * 20, 20);
  AllMeans.SetSize(23, 20);

  int start_counter = 0;
  for (unsigned int i = 0; i <PCA_T1.Rows(); i++)
    for (unsigned int j = 0; j < PCA_T1.Cols(); j++)
      AllPCAs(i+start_counter, j) = PCA_T1(i, j);

  start_counter = start_counter + 20;
  for (unsigned int i = 0; i <PCA_T1CE.Rows(); i++)
    for (unsigned int j = 0; j < PCA_T1CE.Cols(); j++)
      AllPCAs(i + start_counter, j) = PCA_T1CE(i, j);

  start_counter = start_counter + 20;
  for (unsigned int i = 0; i <PCA_T2.Rows(); i++)
    for (unsigned int j = 0; j < PCA_T2.Cols(); j++)
      AllPCAs(i + start_counter, j) = PCA_T2(i, j);

  start_counter = start_counter + 20;
  for (unsigned int i = 0; i <PCA_FL.Rows(); i++)
    for (unsigned int j = 0; j < PCA_FL.Cols(); j++)
      AllPCAs(i + start_counter, j) = PCA_FL(i, j);

  start_counter = start_counter + 20;
  for (unsigned int i = 0; i <PCA_T1T1CE.Rows(); i++)
    for (unsigned int j = 0; j < PCA_T1T1CE.Cols(); j++)
      AllPCAs(i + start_counter, j) = PCA_T1T1CE(i, j);

  start_counter = start_counter + 20;
  for (unsigned int i = 0; i <PCA_T2FL.Rows(); i++)
    for (unsigned int j = 0; j < PCA_T2FL.Cols(); j++)
      AllPCAs(i + start_counter, j) = PCA_T2FL(i, j);

  start_counter = start_counter + 20;
  for (unsigned int i = 0; i <PCA_AX.Rows(); i++)
    for (unsigned int j = 0; j < PCA_AX.Cols(); j++)
      AllPCAs(i + start_counter, j) = PCA_AX(i, j);

  start_counter = start_counter + 20;
  for (unsigned int i = 0; i <PCA_FA.Rows(); i++)
    for (unsigned int j = 0; j < PCA_FA.Cols(); j++)
      AllPCAs(i + start_counter, j) = PCA_FA(i, j);

  start_counter = start_counter + 20;
  for (unsigned int i = 0; i <PCA_RAD.Rows(); i++)
    for (unsigned int j = 0; j < PCA_RAD.Cols(); j++)
      AllPCAs(i + start_counter, j) = PCA_RAD(i, j);

  start_counter = start_counter + 20;
  for (unsigned int i = 0; i <PCA_TR.Rows(); i++)
    for (unsigned int j = 0; j < PCA_TR.Cols(); j++)
      AllPCAs(i + start_counter, j) = PCA_TR(i, j);

  start_counter = start_counter + 20;
  for (unsigned int i = 0; i <PCA_PH.Rows(); i++)
    for (unsigned int j = 0; j < PCA_PH.Cols(); j++)
      AllPCAs(i + start_counter, j) = PCA_PH(i, j);

  start_counter = start_counter + 20;
  for (unsigned int i = 0; i <PCA_PSR.Rows(); i++)
    for (unsigned int j = 0; j < PCA_PSR.Cols(); j++)
      AllPCAs(i + start_counter, j) = PCA_PSR(i, j);

  start_counter = start_counter + 20;
  for (unsigned int i = 0; i <PCA_RCBV.Rows(); i++)
    for (unsigned int j = 0; j < PCA_RCBV.Cols(); j++)
      AllPCAs(i + start_counter, j) = PCA_RCBV(i, j);

  start_counter = start_counter + 20;
  for (unsigned int i = 0; i <PCA_PC1.Rows(); i++)
    for (unsigned int j = 0; j < PCA_PC1.Cols(); j++)
      AllPCAs(i + start_counter, j) = PCA_PC1(i, j);

  start_counter = start_counter + 20;
  for (unsigned int i = 0; i <PCA_PC2.Rows(); i++)
    for (unsigned int j = 0; j < PCA_PC2.Cols(); j++)
      AllPCAs(i + start_counter, j) = PCA_PC2(i, j);

  start_counter = start_counter + 20;
  for (unsigned int i = 0; i <PCA_PC3.Rows(); i++)
    for (unsigned int j = 0; j < PCA_PC3.Cols(); j++)
      AllPCAs(i + start_counter, j) = PCA_PC3(i, j);

  start_counter = start_counter + 20;
  for (unsigned int i = 0; i <PCA_PC4.Rows(); i++)
    for (unsigned int j = 0; j < PCA_PC4.Cols(); j++)
      AllPCAs(i + start_counter, j) = PCA_PC4(i, j);

  start_counter = start_counter + 20;
  for (unsigned int i = 0; i <PCA_PC5.Rows(); i++)
    for (unsigned int j = 0; j < PCA_PC5.Cols(); j++)
      AllPCAs(i + start_counter, j) = PCA_PC5(i, j);

  start_counter = start_counter + 20;
  for (unsigned int i = 0; i <PCA_PC6.Rows(); i++)
    for (unsigned int j = 0; j < PCA_PC6.Cols(); j++)
      AllPCAs(i + start_counter, j) = PCA_PC6(i, j);

  start_counter = start_counter + 20;
  for (unsigned int i = 0; i <PCA_PC7.Rows(); i++)
    for (unsigned int j = 0; j < PCA_PC7.Cols(); j++)
      AllPCAs(i + start_counter, j) = PCA_PC7(i, j);

  start_counter = start_counter + 20;
  for (unsigned int i = 0; i <PCA_PC8.Rows(); i++)
    for (unsigned int j = 0; j < PCA_PC8.Cols(); j++)
      AllPCAs(i + start_counter, j) = PCA_PC8(i, j);

  start_counter = start_counter + 20;
  for (unsigned int i = 0; i <PCA_PC9.Rows(); i++)
    for (unsigned int j = 0; j < PCA_PC9.Cols(); j++)
      AllPCAs(i + start_counter, j) = PCA_PC9(i, j);

  start_counter = start_counter + 20;
  for (unsigned int i = 0; i <PCA_PC10.Rows(); i++)
    for (unsigned int j = 0; j < PCA_PC10.Cols(); j++)
      AllPCAs(i + start_counter, j) = PCA_PC10(i, j);

  for (unsigned int i = 0; i < Mean_T1.Size(); i++)
  {
    AllMeans(0,i) = Mean_T1[i];
    AllMeans(1, i) = Mean_T1CE[i];
    AllMeans(2, i) = Mean_T2[i];
    AllMeans(3, i) = Mean_FL[i];
    AllMeans(4, i) = Mean_T1T1CE[i];
    AllMeans(5, i) = Mean_T2FL[i];
    AllMeans(6, i) = Mean_AX[i];
    AllMeans(7, i) = Mean_FA[i];
    AllMeans(8, i) = Mean_RAD[i];
    AllMeans(9, i) = Mean_TR[i];
    AllMeans(10, i) = Mean_PH[i];
    AllMeans(11, i) = Mean_PSR[i];
    AllMeans(12, i) = Mean_RCBV[i];
    AllMeans(13, i) = Mean_PC1[i];
    AllMeans(14, i) = Mean_PC2[i];
    AllMeans(15, i) = Mean_PC3[i];
    AllMeans(16, i) = Mean_PC4[i];
    AllMeans(17, i) = Mean_PC5[i];
    AllMeans(18, i) = Mean_PC6[i];
    AllMeans(19, i) = Mean_PC7[i];
    AllMeans(20, i) = Mean_PC8[i];
    AllMeans(21, i) = Mean_PC9[i];
    AllMeans(22, i) = Mean_PC10[i];
  }
  WriteCSVFiles(AllPCAs, outputdirectory + "/PCA_Others.csv");
  WriteCSVFiles(AllMeans, outputdirectory + "/Mean_Others.csv");

  std::cout << "pca modalities perfusion components extracted" << std::endl;

  VectorVectorDouble Features;
  for (int i = 0; i < T1ReducedIntensityHistogram.GetPointer()->GetNumberOfRows(); i++)
  {
    std::cout << "patient number" << i << std::endl;
    VectorDouble OnePatient;
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(T1ReducedIntensityHistogram->GetValue(i, j).ToDouble());
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(TCReducedIntensityHistogram->GetValue(i, j).ToDouble());
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(T1TCReducedIntensityHistogram->GetValue(i, j).ToDouble());
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(T2ReducedIntensityHistogram->GetValue(i, j).ToDouble());
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(FLReducedIntensityHistogram->GetValue(i, j).ToDouble());
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(T2FLReducedIntensityHistogram->GetValue(i, j).ToDouble());
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(AXReducedIntensityHistogram->GetValue(i, j).ToDouble());
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(FAReducedIntensityHistogram->GetValue(i, j).ToDouble());
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(RDReducedIntensityHistogram->GetValue(i, j).ToDouble());
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(TRReducedIntensityHistogram->GetValue(i, j).ToDouble());
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(PHReducedIntensityHistogram->GetValue(i, j).ToDouble());
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(PSReducedIntensityHistogram->GetValue(i, j).ToDouble());
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(RCReducedIntensityHistogram->GetValue(i, j).ToDouble());

    std::cout << "One patient size" << OnePatient.size() << std::endl;
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(PC1ReducedIntensityHistogram->GetValue(i, j).ToDouble());
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(PC2ReducedIntensityHistogram->GetValue(i, j).ToDouble());
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(PC3ReducedIntensityHistogram->GetValue(i, j).ToDouble());
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(PC4ReducedIntensityHistogram->GetValue(i, j).ToDouble());
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(PC5ReducedIntensityHistogram->GetValue(i, j).ToDouble());
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(PC6ReducedIntensityHistogram->GetValue(i, j).ToDouble());
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(PC7ReducedIntensityHistogram->GetValue(i, j).ToDouble());
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(PC8ReducedIntensityHistogram->GetValue(i, j).ToDouble());
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(PC9ReducedIntensityHistogram->GetValue(i, j).ToDouble());
    for (int j = 0; j < 10; j++)
      OnePatient.push_back(PC10ReducedIntensityHistogram->GetValue(i, j).ToDouble());

    std::cout << "One patient size" << OnePatient.size() << std::endl;

    Features.push_back(OnePatient);
  }
  return Features;
}


VariableSizeMatrixType PseudoProgressionEstimator::ColumnWiseScaling(VariableSizeMatrixType inputdata)
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

VariableSizeMatrixType PseudoProgressionEstimator::SelectModelFeatures(const VariableSizeMatrixType &ModelFeatures)
{
  VariableSizeMatrixType Data;
  Data.SetSize(ModelFeatures.Rows(), ModelFeatures.Cols() - 1);
  std::vector<double> Labels;

  for (unsigned int i = 0; i < ModelFeatures.Rows(); i++)
    for (unsigned int j = 0; j < ModelFeatures.Cols() - 1; j++)
      Data(i, j) = ModelFeatures(i, j);

  for (unsigned int i = 0; i < ModelFeatures.Rows(); i++)
    Labels.push_back(ModelFeatures(i, ModelFeatures.Cols() - 1));

  double numberOfSelectedFeatures = 0;
  VectorDouble EffectSize = EffectSizeFeatureSelection(Data, Labels);
  for (int eSizeCounter = 0; eSizeCounter < EffectSize.size(); eSizeCounter++)
  {
    if (EffectSize[eSizeCounter] < 0)
      EffectSize[eSizeCounter] = EffectSize[eSizeCounter] * -1;
  }

  std::vector<size_t> indices = sort_indexes(EffectSize);

  //copy the selected features
  VariableSizeMatrixType reducedFeatureSet;
  reducedFeatureSet.SetSize(ModelFeatures.Rows(), 51);
  for (int featureNo = 0; featureNo < 50; featureNo++)
  {
    for (unsigned int sampleNo = 0; sampleNo < Data.Rows(); sampleNo++)
      reducedFeatureSet(sampleNo, featureNo) = Data(sampleNo, indices[featureNo]);
  }
  //copy the label
  for (unsigned int sampleNo = 0; sampleNo < Data.Rows(); sampleNo++)
    reducedFeatureSet(sampleNo, 50) = Labels[sampleNo];

  return reducedFeatureSet;
}
VectorDouble PseudoProgressionEstimator::EffectSizeFeatureSelection(const VariableSizeMatrixType training_features, std::vector<double> target)
{
  //make set 1and set2
  int NoOfSamplesC1 = 0;
  int NoOfSamplesC2 = 0;
  std::vector<double> indices_set1;
  std::vector<double> indices_set2;
  VariableSizeMatrixType features_set1;
  VariableSizeMatrixType features_set2;
  VariableLengthVectorType mean_set1;
  VariableLengthVectorType mean_set2;

  for (int index = 0; index < target.size(); index++)
  {
    if (target[index] == -1)
      NoOfSamplesC1++;
    else if (target[index] == 1)
      NoOfSamplesC2++;
  }
  features_set1.SetSize(NoOfSamplesC1, training_features.Cols());
  features_set2.SetSize(NoOfSamplesC2, training_features.Cols());
  mean_set1.SetSize(training_features.Cols());
  mean_set2.SetSize(training_features.Cols());

  NoOfSamplesC1 = 0;
  NoOfSamplesC2 = 0;
  for (int index = 0; index < target.size(); index++)
  {
    if (target[index] == -1)
    {
      for (unsigned int featureNo = 0; featureNo < training_features.Cols(); featureNo++)
        features_set1(NoOfSamplesC1, featureNo) = training_features(index, featureNo);
      NoOfSamplesC1++;
    }
    else if (target[index] == 1)
    {
      for (unsigned int featureNo = 0; featureNo < training_features.Cols(); featureNo++)
        features_set2(NoOfSamplesC2, featureNo) = training_features(index, featureNo);
      NoOfSamplesC2++;
    }
  }
  std::vector<double> EffectSize;
  for (unsigned int featureNo = 0; featureNo < training_features.Cols(); featureNo++)
  {
    double temp = 0.0;
    for (int sampleNo = 0; sampleNo < NoOfSamplesC1; sampleNo++)
      temp = temp + features_set1(sampleNo, featureNo);
    mean_set1[featureNo] = temp / NoOfSamplesC1;

    temp = 0.0;
    for (int sampleNo = 0; sampleNo < NoOfSamplesC2; sampleNo++)
      temp = temp + features_set2(sampleNo, featureNo);
    mean_set2[featureNo] = temp / NoOfSamplesC2;


    double sum1 = 0;
    double sum2 = 0;
    for (int sampleNo = 0; sampleNo < NoOfSamplesC1; sampleNo++)
      sum1 = sum1 + (features_set1(sampleNo, featureNo) - mean_set1[featureNo])*(features_set1(sampleNo, featureNo) - mean_set1[featureNo]);

    for (int sampleNo = 0; sampleNo < NoOfSamplesC2; sampleNo++)
      sum2 = sum2 + (features_set2(sampleNo, featureNo) - mean_set2[featureNo])*(features_set2(sampleNo, featureNo) - mean_set2[featureNo]);

    double SC1 = sum1 / (NoOfSamplesC1 - 1);
    double SC2 = sum2 / (NoOfSamplesC2 - 1);
    double SP = ((NoOfSamplesC1 - 1)*SC1 + (NoOfSamplesC2 - 1)*SC2) / (NoOfSamplesC1 + NoOfSamplesC2 - 2);
    EffectSize.push_back((mean_set1[featureNo] - mean_set2[featureNo]) / sqrt(SP));
  }
  //std::vector<size_t> indices = sort_indexes(EffectSize);
  //VariableSizeMatrixType selected_feature_set;

  //for (int index1 = 0; index1 < training_features.Rows(); index1++)
  //	for (int index = 0; index < no_of_features; index++)
  //		selected_feature_set(index1, index) = training_features(index1, indices[index]);

  //return selected_feature_set;
  //EffectSize(find(std::isnan(EffectSize))) = 0.0001;
  return EffectSize;
}

template <typename T>
std::vector<size_t> PseudoProgressionEstimator::sort_indexes(const std::vector<T> &v)
{
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] > v[i2]; });

  return idx;
}

VectorDouble PseudoProgressionEstimator::CombineEstimates(const VariableLengthVectorType &estimates1, const VariableLengthVectorType &estimates2)
{
  VectorDouble returnVec;
  returnVec.resize(estimates1.Size());
  for (size_t i = 0; i < estimates1.Size(); i++)
  {
    float temp_abs, temp_pos1, temp_neg1, temp_1, temp_2;
    // estimate for 1st vector
    if (std::abs(estimates1[i]) < 2)
    {
      temp_abs = estimates1[i];
    }
    else
    {
      temp_abs = 0;
    }

    if (estimates1[i] > 1)
    {
      temp_pos1 = 1;
    }
    else
    {
      temp_pos1 = 0;
    }

    if (estimates1[i] < -1)
    {
      temp_neg1 = 1;
    }
    else
    {
      temp_neg1 = 0;
    }
    temp_1 = temp_abs + (temp_pos1 - temp_neg1);

    // estimate for 2nd vector, all temp values are getting overwritten
    if (std::abs(estimates2[i]) < 2)
    {
      temp_abs = estimates2[i];
    }
    else
    {
      temp_abs = 0;
    }

    if (estimates2[i] > 1)
    {
      temp_pos1 = 1;
    }
    else
    {
      temp_pos1 = 0;
    }

    if (estimates2[i] < -1)
    {
      temp_neg1 = 1;
    }
    else
    {
      temp_neg1 = 0;
    }
    temp_2 = temp_abs + (temp_pos1 - temp_neg1);

    // combine the two
    returnVec[i] = temp_1 + temp_2;
  }
  return returnVec;
}
void PseudoProgressionEstimator::WritePCAOutputs(std::string suffix, std::string outputdirectory, const VariableLengthVectorType mean, const VariableSizeMatrixType coefs)
{
  std::ofstream myfile;
  myfile.open(outputdirectory + "/mean_" + suffix + ".csv");
  for (unsigned int index1 = 0; index1 < mFeatureReductionLocalPtr.GetPerfusionMeanVector().Size(); index1++)
    myfile << std::to_string(mFeatureReductionLocalPtr.GetPerfusionMeanVector()[index1]) + "\n";
  myfile.close();

  myfile.open(outputdirectory + "/pca_" + suffix + ".csv");
  for (unsigned int index1 = 0; index1 < mFeatureReductionLocalPtr.GetPCATransformationMatrix().Rows(); index1++)
  {
    std::string onerow;
    for (unsigned int index2 = 0; index2 < mFeatureReductionLocalPtr.GetPCATransformationMatrix().Cols(); index2++)
    {
      if (index2 == 0)
        myfile << std::to_string(mFeatureReductionLocalPtr.GetPCATransformationMatrix()[index1][index2]);
      else
        myfile << "," + std::to_string(mFeatureReductionLocalPtr.GetPCATransformationMatrix()[index1][index2]);
    }
    myfile << "\n";
  }
  myfile.close();
}

void PseudoProgressionEstimator::ReadAllTheModelParameters(std::string modeldirectory,
  VariableSizeMatrixType &PCA_PERF,
  VariableSizeMatrixType &PCA_T1,
  VariableSizeMatrixType &PCA_T1CE,
  VariableSizeMatrixType &PCA_T2,
  VariableSizeMatrixType &PCA_FL,
  VariableSizeMatrixType &PCA_T1T1CE,
  VariableSizeMatrixType &PCA_T2FL,
  VariableSizeMatrixType &PCA_AX,
  VariableSizeMatrixType &PCA_FA,
  VariableSizeMatrixType &PCA_RAD,
  VariableSizeMatrixType &PCA_TR,
  VariableSizeMatrixType &PCA_PH,
  VariableSizeMatrixType &PCA_PSR,
  VariableSizeMatrixType &PCA_RCBV,
  VariableSizeMatrixType &PCA_PC1,
  VariableSizeMatrixType &PCA_PC2,
  VariableSizeMatrixType &PCA_PC3,
  VariableSizeMatrixType &PCA_PC4,
  VariableSizeMatrixType &PCA_PC5,
  VariableSizeMatrixType &PCA_PC6,
  VariableSizeMatrixType &PCA_PC7,
  VariableSizeMatrixType &PCA_PC8,
  VariableSizeMatrixType &PCA_PC9,
  VariableSizeMatrixType &PCA_PC10,
  VariableLengthVectorType &Mean_PERF,
  VariableLengthVectorType &Mean_T1,
  VariableLengthVectorType &Mean_T1CE,
  VariableLengthVectorType &Mean_T2,
  VariableLengthVectorType &Mean_FL,
  VariableLengthVectorType &Mean_T1T1CE,
  VariableLengthVectorType &Mean_T2FL,
  VariableLengthVectorType &Mean_AX,
  VariableLengthVectorType &Mean_FA,
  VariableLengthVectorType &Mean_RAD,
  VariableLengthVectorType &Mean_TR,
  VariableLengthVectorType &Mean_PH,
  VariableLengthVectorType &Mean_PSR,
  VariableLengthVectorType &Mean_RCBV,
  VariableLengthVectorType &Mean_PC1,
  VariableLengthVectorType &Mean_PC2,
  VariableLengthVectorType &Mean_PC3,
  VariableLengthVectorType &Mean_PC4,
  VariableLengthVectorType &Mean_PC5,
  VariableLengthVectorType &Mean_PC6,
  VariableLengthVectorType &Mean_PC7,
  VariableLengthVectorType &Mean_PC8,
  VariableLengthVectorType &Mean_PC9,
  VariableLengthVectorType &Mean_PC10)
{
  CSVFileReaderType::Pointer reader = CSVFileReaderType::New();
  MatrixType dataMatrix;
  VariableLengthVectorType meanMatrix;
  reader->SetFieldDelimiterCharacter(',');
  reader->HasColumnHeadersOff();
  reader->HasRowHeadersOff();

  std::cout << "Reading model files" << std::endl;
  //-------------perfusion related data reading------------------
  reader->SetFileName(modeldirectory + "/PCA_PERF.csv");
  reader->Parse();
  dataMatrix = reader->GetArray2DDataObject()->GetMatrix();
  PCA_PERF.SetSize(dataMatrix.rows(), dataMatrix.cols());
  for (unsigned int i = 0; i < dataMatrix.rows(); i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_PERF(i, j) = dataMatrix(i, j);

  reader->SetFileName(modeldirectory + "/Mean_PERF.csv");
  reader->Parse();
  dataMatrix = reader->GetArray2DDataObject()->GetMatrix();
  Mean_PERF.SetSize(dataMatrix.size());
  for (unsigned int i = 0; i < dataMatrix.size(); i++)
    Mean_PERF[i] = dataMatrix(0, i);

  //-------------others related data reading------------------
  int PCA_Others_Size = 20;
  reader->SetFileName(modeldirectory + "/PCA_Others.csv");
  reader->Parse();
  dataMatrix = reader->GetArray2DDataObject()->GetMatrix();

  PCA_T1.SetSize(PCA_Others_Size, PCA_Others_Size);
  PCA_T1CE.SetSize(PCA_Others_Size, PCA_Others_Size);
  PCA_T2.SetSize(PCA_Others_Size, PCA_Others_Size);
  PCA_FL.SetSize(PCA_Others_Size, PCA_Others_Size);
  PCA_T1T1CE.SetSize(PCA_Others_Size, PCA_Others_Size);
  PCA_T2FL.SetSize(PCA_Others_Size, PCA_Others_Size);
  PCA_AX.SetSize(PCA_Others_Size, PCA_Others_Size);
  PCA_FA.SetSize(PCA_Others_Size, PCA_Others_Size);
  PCA_RAD.SetSize(PCA_Others_Size, PCA_Others_Size);
  PCA_TR.SetSize(PCA_Others_Size, PCA_Others_Size);
  PCA_PH.SetSize(PCA_Others_Size, PCA_Others_Size);
  PCA_PSR.SetSize(PCA_Others_Size, PCA_Others_Size);
  PCA_RCBV.SetSize(PCA_Others_Size, PCA_Others_Size);
  PCA_PC1.SetSize(PCA_Others_Size, PCA_Others_Size);
  PCA_PC2.SetSize(PCA_Others_Size, PCA_Others_Size);
  PCA_PC3.SetSize(PCA_Others_Size, PCA_Others_Size);
  PCA_PC4.SetSize(PCA_Others_Size, PCA_Others_Size);
  PCA_PC5.SetSize(PCA_Others_Size, PCA_Others_Size);
  PCA_PC6.SetSize(PCA_Others_Size, PCA_Others_Size);
  PCA_PC7.SetSize(PCA_Others_Size, PCA_Others_Size);
  PCA_PC8.SetSize(PCA_Others_Size, PCA_Others_Size);
  PCA_PC9.SetSize(PCA_Others_Size, PCA_Others_Size);
  PCA_PC10.SetSize(PCA_Others_Size, PCA_Others_Size);

  Mean_T1.SetSize(PCA_Others_Size);
  Mean_T1CE.SetSize(PCA_Others_Size);
  Mean_T2.SetSize(PCA_Others_Size);
  Mean_FL.SetSize(PCA_Others_Size);
  Mean_T1T1CE.SetSize(PCA_Others_Size);
  Mean_T2FL.SetSize(PCA_Others_Size);
  Mean_AX.SetSize(PCA_Others_Size);
  Mean_FA.SetSize(PCA_Others_Size);
  Mean_RAD.SetSize(PCA_Others_Size);
  Mean_TR.SetSize(PCA_Others_Size);
  Mean_PH.SetSize(PCA_Others_Size);
  Mean_PSR.SetSize(PCA_Others_Size);
  Mean_RCBV.SetSize(PCA_Others_Size);
  Mean_PC1.SetSize(PCA_Others_Size);
  Mean_PC2.SetSize(PCA_Others_Size);
  Mean_PC3.SetSize(PCA_Others_Size);
  Mean_PC4.SetSize(PCA_Others_Size);
  Mean_PC5.SetSize(PCA_Others_Size);
  Mean_PC6.SetSize(PCA_Others_Size);
  Mean_PC7.SetSize(PCA_Others_Size);
  Mean_PC8.SetSize(PCA_Others_Size);
  Mean_PC9.SetSize(PCA_Others_Size);
  Mean_PC10.SetSize(PCA_Others_Size);

  int start_counter = 0;
  int end_counter = 19;
  for (int i = start_counter; i <= end_counter; i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_T1(i, j) = dataMatrix(i, j);

  std::cout << "PCA_T1 read.\n";

  start_counter = start_counter + PCA_Others_Size;
  end_counter = end_counter + PCA_Others_Size;
  for (int i = start_counter; i <= end_counter; i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_T1CE(i - start_counter, j) = dataMatrix(i, j);

  std::cout << "PCA_T1CE read.\n";

  start_counter = start_counter + PCA_Others_Size;
  end_counter = end_counter + PCA_Others_Size;
  for (int i = start_counter; i <= end_counter; i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_T2(i - start_counter, j) = dataMatrix(i, j);

  std::cout << "PCA_T2 read.\n";

  start_counter = start_counter + PCA_Others_Size;
  end_counter = end_counter + PCA_Others_Size;
  for (int i = start_counter; i <= end_counter; i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_FL(i - start_counter, j) = dataMatrix(i, j);

  std::cout << "PCA_FL read.\n";

  start_counter = start_counter + PCA_Others_Size;
  end_counter = end_counter + PCA_Others_Size;
  for (int i = start_counter; i <= end_counter; i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_T1T1CE(i - start_counter, j) = dataMatrix(i, j);

  std::cout << "PCA_T1T1CE read.\n";

  start_counter = start_counter + PCA_Others_Size;
  end_counter = end_counter + PCA_Others_Size;
  for (int i = start_counter; i <= end_counter; i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_T2FL(i - start_counter, j) = dataMatrix(i, j);

  std::cout << "PCA_T2FL read.\n";

  start_counter = start_counter + PCA_Others_Size;
  end_counter = end_counter + PCA_Others_Size;
  for (int i = start_counter; i <= end_counter; i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_AX(i - start_counter, j) = dataMatrix(i, j);

  std::cout << "PCA_AX read.\n";

  start_counter = start_counter + PCA_Others_Size;
  end_counter = end_counter + PCA_Others_Size;
  for (int i = start_counter; i <= end_counter; i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_FA(i - start_counter, j) = dataMatrix(i, j);

  std::cout << "PCA_FA read.\n";

  start_counter = start_counter + PCA_Others_Size;
  end_counter = end_counter + PCA_Others_Size;
  for (int i = start_counter; i <= end_counter; i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_RAD(i - start_counter, j) = dataMatrix(i, j);

  std::cout << "PCA_RAD read.\n";

  start_counter = start_counter + PCA_Others_Size;
  end_counter = end_counter + PCA_Others_Size;
  for (int i = start_counter; i <= end_counter; i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_TR(i - start_counter, j) = dataMatrix(i, j);

  std::cout << "PCA_TR read.\n";

  start_counter = start_counter + PCA_Others_Size;
  end_counter = end_counter + PCA_Others_Size;
  for (int i = start_counter; i <= end_counter; i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_PH(i - start_counter, j) = dataMatrix(i, j);

  std::cout << "PCA_PH read.\n";

  start_counter = start_counter + PCA_Others_Size;
  end_counter = end_counter + PCA_Others_Size;
  for (int i = start_counter; i <= end_counter; i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_PSR(i - start_counter, j) = dataMatrix(i, j);

  std::cout << "PCA_PSR read.\n";

  start_counter = start_counter + PCA_Others_Size;
  end_counter = end_counter + PCA_Others_Size;
  for (int i = start_counter; i <= end_counter; i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_RCBV(i - start_counter, j) = dataMatrix(i, j);

  std::cout << "PCA_RCBV read.\n";

  start_counter = start_counter + PCA_Others_Size;
  end_counter = end_counter + PCA_Others_Size;
  for (int i = start_counter; i <= end_counter; i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_PC1(i - start_counter, j) = dataMatrix(i, j);

  std::cout << "PCA_PC1 read.\n";

  start_counter = start_counter + PCA_Others_Size;
  end_counter = end_counter + PCA_Others_Size;
  for (int i = start_counter; i <= end_counter; i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_PC2(i - start_counter, j) = dataMatrix(i, j);

  std::cout << "PCA_PC2 read.\n";

  start_counter = start_counter + PCA_Others_Size;
  end_counter = end_counter + PCA_Others_Size;
  for (int i = start_counter; i <= end_counter; i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_PC3(i - start_counter, j) = dataMatrix(i, j);

  std::cout << "PCA_PC3 read.\n";

  start_counter = start_counter + PCA_Others_Size;
  end_counter = end_counter + PCA_Others_Size;
  for (int i = start_counter; i <= end_counter; i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_PC4(i - start_counter, j) = dataMatrix(i, j);

  std::cout << "PCA_PC4 read.\n";

  start_counter = start_counter + PCA_Others_Size;
  end_counter = end_counter + PCA_Others_Size;
  for (int i = start_counter; i <= end_counter; i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_PC5(i - start_counter, j) = dataMatrix(i, j);

  std::cout << "PCA_PC5 read.\n";

  start_counter = start_counter + PCA_Others_Size;
  end_counter = end_counter + PCA_Others_Size;
  for (int i = start_counter; i <= end_counter; i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_PC6(i - start_counter, j) = dataMatrix(i, j);

  std::cout << "PCA_PC6 read.\n";

  start_counter = start_counter + PCA_Others_Size;
  end_counter = end_counter + PCA_Others_Size;
  for (int i = start_counter; i <= end_counter; i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_PC7(i - start_counter, j) = dataMatrix(i, j);

  std::cout << "PCA_PC7 read.\n";

  start_counter = start_counter + PCA_Others_Size;
  end_counter = end_counter + PCA_Others_Size;
  for (int i = start_counter; i <= end_counter; i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_PC8(i - start_counter, j) = dataMatrix(i, j);

  std::cout << "PCA_PC8 read.\n";

  start_counter = start_counter + PCA_Others_Size;
  end_counter = end_counter + PCA_Others_Size;
  for (int i = start_counter; i <= end_counter; i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_PC9(i - start_counter, j) = dataMatrix(i, j);

  std::cout << "PCA_PC9 read.\n";

  start_counter = start_counter + PCA_Others_Size;
  end_counter = end_counter + PCA_Others_Size;
  for (int i = start_counter; i <= end_counter; i++)
    for (unsigned int j = 0; j < dataMatrix.cols(); j++)
      PCA_PC10(i - start_counter, j) = dataMatrix(i, j);

  std::cout << "PCA_PC10 read.\n";

  reader->SetFileName(modeldirectory + "/Mean_Others.csv");
  reader->Parse();
  dataMatrix = reader->GetArray2DDataObject()->GetMatrix();

  for (unsigned int i = 0; i < dataMatrix.cols(); i++)
  {
    Mean_T1[i] = dataMatrix(0, i);
    Mean_T1CE[i] = dataMatrix(1, i);
    Mean_T2[i] = dataMatrix(2, i);
    Mean_FL[i] = dataMatrix(3, i);
    Mean_T1T1CE[i] = dataMatrix(4, i);
    Mean_T2FL[i] = dataMatrix(5, i);
    Mean_AX[i] = dataMatrix(6, i);
    Mean_FA[i] = dataMatrix(7, i);
    Mean_RAD[i] = dataMatrix(8, i);
    Mean_TR[i] = dataMatrix(9, i);
    Mean_PH[i] = dataMatrix(10, i);
    Mean_PSR[i] = dataMatrix(11, i);
    Mean_RCBV[i] = dataMatrix(12, i);
    Mean_PC1[i] = dataMatrix(13, i);
    Mean_PC2[i] = dataMatrix(14, i);
    Mean_PC3[i] = dataMatrix(15, i);
    Mean_PC4[i] = dataMatrix(16, i);
    Mean_PC5[i] = dataMatrix(17, i);
    Mean_PC6[i] = dataMatrix(18, i);
    Mean_PC7[i] = dataMatrix(19, i);
    Mean_PC8[i] = dataMatrix(20, i);
    Mean_PC9[i] = dataMatrix(21, i);
    Mean_PC10[i] = dataMatrix(22, i);
  }
  std::cout << "Finished reading all model files" << std::endl;
  int a = 0;
}


void PseudoProgressionEstimator::WriteCSVFiles(VariableSizeMatrixType inputdata, std::string filepath)
{
  std::ofstream myfile;
  myfile.open(filepath);
  for (unsigned int index1 = 0; index1 < inputdata.Rows(); index1++)
  {
    for (unsigned int index2 = 0; index2 < inputdata.Cols(); index2++)
    {
      if (index2 == 0)
        myfile << std::to_string(inputdata[index1][index2]);
      else
        myfile << "," << std::to_string(inputdata[index1][index2]);
    }
    myfile << "\n";
  }
  myfile.close();
}
void PseudoProgressionEstimator::WriteCSVFiles(VectorVectorDouble inputdata, std::string filepath)
{
  std::ofstream myfile;
  myfile.open(filepath);
  for (unsigned int index1 = 0; index1 < inputdata.size(); index1++)
  {
    for (unsigned int index2 = 0; index2 < inputdata[0].size(); index2++)
    {
      if (index2 == 0)
        myfile << std::to_string(inputdata[index1][index2]);
      else
        myfile << "," << std::to_string(inputdata[index1][index2]);
    }
    myfile << "\n";
  }
  myfile.close();
}
void PseudoProgressionEstimator::WriteCSVFiles(vtkSmartPointer<vtkTable> inputdata, std::string filepath)
{
  std::ofstream myfile;
  myfile.open(filepath);
  for (unsigned int index1 = 0; index1 < inputdata.GetPointer()->GetNumberOfRows(); index1++)
  {
    for (unsigned int index2 = 0; index2 < inputdata.GetPointer()->GetNumberOfColumns(); index2++)
    {
      if (index2 == 0)
        myfile << std::to_string(inputdata->GetValue(index1, index2).ToDouble());
      else
        myfile << "," << std::to_string(inputdata->GetValue(index1, index2).ToDouble());
    }
    myfile << "\n";
  }
  myfile.close();
}
void PseudoProgressionEstimator::WriteCSVFiles(VariableLengthVectorType inputdata, std::string filepath)
{
  std::ofstream myfile;
  myfile.open(filepath);
  for (unsigned int index1 = 0; index1 < inputdata.Size(); index1++)
    myfile << std::to_string(inputdata[index1]) << ",";

  myfile << "\n";
  myfile.close();
}

void PseudoProgressionEstimator::WriteCSVFiles(std::vector<double> inputdata, std::string filepath)
{
  std::ofstream myfile;
  myfile.open(filepath);
  for (unsigned int index1 = 0; index1 < inputdata.size(); index1++)
    myfile << std::to_string(inputdata[index1]) << ",";

  myfile << "\n";
  myfile.close();
}

VariableSizeMatrixType PseudoProgressionEstimator::GetModelSelectedFeatures(VariableSizeMatrixType & ScaledFeatureSetAfterAddingLabel, VariableLengthVectorType & SelectedFeatures)
{
  VariableSizeMatrixType ModelSelectedFeatures;
  ModelSelectedFeatures.SetSize(ScaledFeatureSetAfterAddingLabel.Rows(), SelectedFeatures.Size() + 1);
  int counter = 0;
  for (unsigned int i = 0; i < SelectedFeatures.Size(); i++)
  {
    for (unsigned int j = 0; j < ScaledFeatureSetAfterAddingLabel.Rows(); j++)
      ModelSelectedFeatures(j, counter) = ScaledFeatureSetAfterAddingLabel(j, SelectedFeatures[i]);
    counter++;
  }
  for (unsigned int j = 0; j < ScaledFeatureSetAfterAddingLabel.Rows(); j++)
    ModelSelectedFeatures(j, counter) = ScaledFeatureSetAfterAddingLabel(j, ScaledFeatureSetAfterAddingLabel.Cols() - 1);

  return ModelSelectedFeatures;
}


vtkSmartPointer<vtkTable> PseudoProgressionEstimator::MakePCAMatrix(int NumberOfFeatures, int NumberOfSamples)
{
  vtkSmartPointer<vtkTable> T1ReducedIntensityHistogram = vtkSmartPointer<vtkTable>::New();
  for (vtkIdType r = 0; r < static_cast<vtkIdType>(NumberOfFeatures); r++)
  {
    vtkSmartPointer<vtkDoubleArray> col = vtkSmartPointer<vtkDoubleArray>::New();
    for (vtkIdType c = 0; c < static_cast<vtkIdType>(NumberOfSamples); c++)
      col->InsertNextValue(0);
    T1ReducedIntensityHistogram->AddColumn(col);
  }
  return T1ReducedIntensityHistogram;
}