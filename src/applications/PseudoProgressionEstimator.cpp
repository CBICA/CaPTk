#include "PseudoProgressionEstimator.h"
#include "fMainWindow.h"
#include "cbicaStatistics.h"
#include "CapTkEnums.h"


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
  VariableSizeMatrixType TrainingData = LoadPseudoProgressionFeaturesData(qualifiedsubjects, traininglabels,outputdirectory);

  //std::ofstream myfile;
  //myfile.open(outputdirectory + "/AllFeatures.csv");
  //for (unsigned int index1 = 0; index1 < TrainingData.Rows(); index1++)
  //{
  //  for (unsigned int index2 = 0; index2 < TrainingData.Cols(); index2++)
  //  {
  //    if (index2 == 0)
  //      myfile << std::to_string(TrainingData[index1][index2]);
  //    else
  //      myfile << "," <<std::to_string(TrainingData[index1][index2]);
  //  }
  //  myfile << "\n";
  //}
  //myfile.close();

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
    std::ofstream myfile;
    myfile.open(outputdirectory + "/PSU_ZScore_Mean.csv");
    for (unsigned int index1 = 0; index1 < meanVector.Size(); index1++)
      myfile << std::to_string(meanVector[index1]) + "\n";
    myfile.close();
    myfile.open(outputdirectory + "/PSU_ZScore_Std.csv");
    for (unsigned int index1 = 0; index1 < stdVector.Size(); index1++)
      myfile << std::to_string(stdVector[index1]) + "\n";
    myfile.close();
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Error in writing output files to the output directory = " + outputdirectory + "Error code : " + std::string(e1.what()));
    return false;
  }

    typedef vnl_matrix<double> MatrixType;
  MatrixType data;
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
  VariableSizeMatrixType TrainingData = LoadPseudoProgressionFeaturesData(qualifiedsubjects, traininglabels,outputdirectory);

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
      mean[i] = meanMatrix(i, 0);
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
      stddevition[i] = stdMatrix(i, 0);
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
  VariableSizeMatrixType PseudoModelSelectedFeatures = SelectModelFeatures(ScaledFeatureSetAfterAddingLabel);
  VariableSizeMatrixType RecurrenceModelSelectedFeatures = SelectModelFeatures(ScaledFeatureSetAfterAddingLabel);

//  std::cout << "selected features done: size:" << PseudoModelSelectedFeatures.Rows() << " columns: " << PseudoModelSelectedFeatures.Cols() << std::endl;
  try
  {
    std::ofstream myfile;
    myfile.open(outputdirectory + "/results.csv");
    myfile << "SubjectName,Score (Pseudo), Score (Recurrence)\n";
    if (cbica::fileExists(modeldirectory + "/PSeudo_SVM_Model.csv") == true && cbica::fileExists(modeldirectory + "/Recurrence_SVM_Model.csv") == true)
    {
      VariableLengthVectorType result_6;
      VariableLengthVectorType result_18;
      result_6 = DistanceFunction(PseudoModelSelectedFeatures, modeldirectory + "/Pseudo_SVM_Model.csv", -1.0927, 0.0313);
      result_18 = DistanceFunction(RecurrenceModelSelectedFeatures, modeldirectory + "/Recurrence_SVM_Model.csv", -0.2854, 0.5);
      results = CombineEstimates(result_6, result_18);
      for (size_t i = 0; i < results.size(); i++)
      {
        std::map<CAPTK::ImageModalityType, std::string> currentsubject = qualifiedsubjects[i];
        myfile << static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]) + "," + std::to_string(result_6[i]) + "," + std::to_string(result_18[i]) + "," + std::to_string(results[i]) + "\n";
      }
    }
    else if (cbica::fileExists(modeldirectory + "/Pseudo_SVM_Model.xml") == true && cbica::fileExists(modeldirectory + "/Recurrence_SVM_Model.xml") == true)
    {
      VectorDouble result_6;
      VectorDouble result_18;
      result_6 = testOpenCVSVM(PseudoModelSelectedFeatures, modeldirectory + "/Pseudo_SVM_Model.xml");
      result_18 = testOpenCVSVM(RecurrenceModelSelectedFeatures, modeldirectory + "/Recurrence_SVM_Model.xml");
      //results = CombineEstimates(result_6, result_18);
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

VariableSizeMatrixType PseudoProgressionEstimator::LoadPseudoProgressionFeaturesData(const std::vector<std::map<CAPTK::ImageModalityType, std::string>> &trainingsubjects, std::vector<double> &traininglabels,std::string outputdirectory)
{
  VariableSizeMatrixType FeaturesOfAllSubjects;
  FeaturesOfAllSubjects.SetSize(trainingsubjects.size(), 810);

  VariableSizeMatrixType otherFeatures;
  otherFeatures.SetSize(trainingsubjects.size(), 810);

  VectorVectorDouble perfusionFeatures;

  VectorVectorDouble T1IntensityHistogram;
  VectorVectorDouble T2IntensityHistogram;
  VectorVectorDouble TCIntensityHistogram;
  VectorVectorDouble T1TCIntensityHistogram;
  VectorVectorDouble FLIntensityHistogram;
  VectorVectorDouble T2FLIntensityHistogram;
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

  //combining perfusion data, calcualting PCA, and putting back in images of respective patients
  //--------------------------------------------------------------------------------------------
  PerfusionMapType perfFeatures = CombineAndCalculatePerfusionPCA(PerfusionDataMap,outputdirectory);
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
      //cbica::WriteImage<ImageType>(CurrentTimePoint, "E:/Projects/PSU/Data/" + std::to_string(sid) + "_" + std::to_string(i) + ".nii.gz");
    }
    RevisedPerfusionImagesOfAllPatients.push_back(OnePatientperfusionImages);
  }

  //Load remaining data of all the patietns
  //----------------------------------------

  for (unsigned int sid = 0; sid < trainingsubjects.size(); sid++)
  {
    std::cout << "Loading Remianing Features: " << sid << std::endl;
    VectorDouble neuroScores;
    std::map<CAPTK::ImageModalityType, std::string> currentsubject = trainingsubjects[sid];

    CSVFileReaderType::Pointer reader = CSVFileReaderType::New();
    MatrixType dataMatrix;
    reader->SetFileName(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_FEATURES]));
    reader->SetFieldDelimiterCharacter(',');
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();
    reader->Parse();
    dataMatrix = reader->GetArray2DDataObject()->GetMatrix();

    for (unsigned int i = 0; i < dataMatrix.rows(); i++)
    {
      neuroScores.push_back(dataMatrix(i, 0));
      neuroScores.push_back(dataMatrix(i, 1));
      neuroScores.push_back(dataMatrix(i, 2));
      neuroScores.push_back(dataMatrix(i, 3));
      neuroScores.push_back(dataMatrix(i, 4));
      neuroScores.push_back(dataMatrix(i, 5));
      traininglabels.push_back(dataMatrix(i, 6));
    }
    ImageType::Pointer LabelImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG]));

    ImageType::Pointer T1CEImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1CE])));
    ImageType::Pointer T2FlairImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR])));
    ImageType::Pointer T1ImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1])));
    ImageType::Pointer T2ImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2])));
    ImageType::Pointer T1T1CEImagePointer = RescaleImageIntensity<ImageType>(MakeAdditionalModality<ImageType>(T1ImagePointer,T1CEImagePointer));
    ImageType::Pointer T2FLImagePointer = RescaleImageIntensity<ImageType>(MakeAdditionalModality<ImageType>(T2ImagePointer, T2FlairImagePointer));

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
    //OtherFeaturesInMap["C4"] = GetAllFeaturesPerImagePerROI<ImageType>(T1T1CEImagePointer, LabelImagePointer, "T1TC");
    //OtherFeaturesInMap["C5"] = GetAllFeaturesPerImagePerROI<ImageType>(T2FLImagePointer, LabelImagePointer, "T2FL");

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


    std::cout << "Shape and neuro features calcualted." << std::endl;
    //10 histogram, 7 intensity, 8 GLCM, 10 GLRLM
    //23 modalities * (10+7+18) = 805
    //805+5 shape features+6 neuro features = 816
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
      //std::cout << "Counter Size" << counter << std::endl;
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
 // VectorVectorDouble PC_Features = CombineAllThePerfusionFeaures(T1IntensityHistogram,TCIntensityHistogram, T1TCIntensityHistogram,
   //                                                              T2IntensityHistogram, FLIntensityHistogram, T2FLIntensityHistogram, AXIntensityHistogram,
     //                                                            FAIntensityHistogram, RDIntensityHistogram, TRIntensityHistogram,
       //                                                          PHIntensityHistogram, PSIntensityHistogram, RCIntensityHistogram,
         //                                                        PCA1IntensityHistogram, PCA2IntensityHistogram, PCA3IntensityHistogram,
           //                                                      PCA4IntensityHistogram, PCA5IntensityHistogram, PCA6IntensityHistogram,
             //                                                    PCA7IntensityHistogram, PCA8IntensityHistogram, PCA9IntensityHistogram,
               //                                                  PCA10IntensityHistogram);

  std::cout << "Final PCA calculation features finished." << std::endl;

  for (unsigned int i = 0; i < FeaturesOfAllSubjects.Rows(); i++)
  {
    for (unsigned int j = 0; j < otherFeatures.Cols(); j++)
      FeaturesOfAllSubjects(i, j) = otherFeatures[i][j];
  /*  for (size_t j = 0; j < PC_Features[i].size(); j++)
      FeaturesOfAllSubjects(i, j + otherFeatures.Cols()) = PC_Features[i][j];*/
  }

  std::cout << "FeaturesOfAllSubjects populated." << std::endl;


  return FeaturesOfAllSubjects;
}


PerfusionMapType PseudoProgressionEstimator::CombineAndCalculatePerfusionPCA(PerfusionMapType PerfusionDataMap,std::string outputdirectory)
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
  vtkSmartPointer<vtkTable> ReducedPCAs = m_featureReduction.GetDiscerningPerfusionTimePoints(CombinedPerfusionFeaturesMap);
  mFeatureReductionLocalPtr.GetPerfusionMeanVector();
  mFeatureReductionLocalPtr.GetPCATransformationMatrix();

  std::ofstream myfile;
  myfile.open(outputdirectory + "/PSU_FullPCA_Mean.csv");
  for (unsigned int index1 = 0; index1 < mFeatureReductionLocalPtr.GetPerfusionMeanVector().Size(); index1++)
    myfile << std::to_string(mFeatureReductionLocalPtr.GetPerfusionMeanVector()[index1]) + "\n";
  myfile.close();

  myfile.open(outputdirectory + "/PSU_FullPCA_Transformation.csv");
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
  VectorVectorDouble PCA10IntensityHistogram)
{
  FeatureReductionClass m_featureReduction;
  vtkSmartPointer<vtkTable> T1ReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> T2ReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> TCReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> T1TCReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> FLReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> T2FLReducedIntensityHistogram;

  vtkSmartPointer<vtkTable> AXReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> FAReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> RDReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> TRReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> PHReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> PSReducedIntensityHistogram;
  vtkSmartPointer<vtkTable> RCReducedIntensityHistogram;
  typedef vnl_matrix<double> MatrixType;
  MatrixType data;
  data.set_size(62, 255);
  typedef itk::CSVNumericObjectFileWriter<double, 62, 255> WriterTypeVector;
  WriterTypeVector::Pointer writerv = WriterTypeVector::New();


  //for (unsigned int i = 0; i < TCIntensityHistogram.size(); i++)
  //  for (unsigned int j = 0; j < TCIntensityHistogram[0].size(); j++)
  //    data(i, j) = TCIntensityHistogram[i][j];

  //writerv->SetFileName("E:/Projects/PSU/Data/TC.csv");
  //writerv->SetInput(&data);
  //writerv->Write();

  //for (unsigned int i = 0; i < TRIntensityHistogram.size(); i++)
  //  for (unsigned int j = 0; j < TRIntensityHistogram[0].size(); j++)
  //    data(i, j) = TRIntensityHistogram[i][j];

  //writerv->SetFileName("E:/Projects/PSU/Data/tr.csv");
  //writerv->SetInput(&data);
  //writerv->Write();
  ////---------------------
  //for (unsigned int i = 0; i < T1IntensityHistogram.size(); i++)
  //  for (unsigned int j = 0; j < T1IntensityHistogram[0].size(); j++)
  //    data(i, j) = T1IntensityHistogram[i][j];

  //writerv->SetFileName("E:/Projects/PSU/Data/T1.csv");
  //writerv->SetInput(&data);
  //writerv->Write();
  //for (unsigned int i = 0; i < T2IntensityHistogram.size(); i++)
  //  for (unsigned int j = 0; j < T2IntensityHistogram[0].size(); j++)
  //    data(i, j) = T2IntensityHistogram[i][j];

  //writerv->SetFileName("E:/Projects/PSU/Data/T2.csv");
  //writerv->SetInput(&data);
  //writerv->Write();
  //for (unsigned int i = 0; i < T1TCIntensityHistogram.size(); i++)
  //  for (unsigned int j = 0; j < T1TCIntensityHistogram[0].size(); j++)
  //    data(i, j) = T1TCIntensityHistogram[i][j];

  //writerv->SetFileName("E:/Projects/PSU/Data/t1tc.csv");
  //writerv->SetInput(&data);
  //writerv->Write();

  //for (unsigned int i = 0; i < FLIntensityHistogram.size(); i++)
  //  for (unsigned int j = 0; j < FLIntensityHistogram[0].size(); j++)
  //    data(i, j) = FLIntensityHistogram[i][j];

  //writerv->SetFileName("E:/Projects/PSU/Data/flair.csv");
  //writerv->SetInput(&data);
  //writerv->Write();

  //for (unsigned int i = 0; i < AXIntensityHistogram.size(); i++)
  //  for (unsigned int j = 0; j < AXIntensityHistogram[0].size(); j++)
  //    data(i, j) = AXIntensityHistogram[i][j];

  //writerv->SetFileName("E:/Projects/PSU/Data/ax.csv");
  //writerv->SetInput(&data);
  //writerv->Write();
  //for (unsigned int i = 0; i < FAIntensityHistogram.size(); i++)
  //  for (unsigned int j = 0; j < FAIntensityHistogram[0].size(); j++)
  //    data(i, j) = FAIntensityHistogram[i][j];

  //writerv->SetFileName("E:/Projects/PSU/Data/fa.csv");
  //writerv->SetInput(&data);
  //writerv->Write();
  //for (unsigned int i = 0; i < RDIntensityHistogram.size(); i++)
  //  for (unsigned int j = 0; j < RDIntensityHistogram[0].size(); j++)
  //    data(i, j) = RDIntensityHistogram[i][j];

  //writerv->SetFileName("E:/Projects/PSU/Data/rad.csv");
  //writerv->SetInput(&data);
  //writerv->Write();

  //for (unsigned int i = 0; i < PCA1IntensityHistogram.size(); i++)
  //  for (unsigned int j = 0; j < PCA1IntensityHistogram[0].size(); j++)
  //    data(i, j) = PCA1IntensityHistogram[i][j];

  //writerv->SetFileName("E:/Projects/PSU/Data/pca1.csv");
  //writerv->SetInput(&data);
  //writerv->Write();


  //for (unsigned int i = 0; i < PCA2IntensityHistogram.size(); i++)
  //  for (unsigned int j = 0; j < PCA2IntensityHistogram[0].size(); j++)
  //    data(i, j) = PCA2IntensityHistogram[i][j];

  //writerv->SetFileName("E:/Projects/PSU/Data/pca2.csv");
  //writerv->SetInput(&data);
  //writerv->Write();

  //for (unsigned int i = 0; i < PCA3IntensityHistogram.size(); i++)
  //  for (unsigned int j = 0; j < PCA3IntensityHistogram[0].size(); j++)
  //    data(i, j) = PCA3IntensityHistogram[i][j];

  //writerv->SetFileName("E:/Projects/PSU/Data/pca3.csv");
  //writerv->SetInput(&data);
  //writerv->Write();


  //for (unsigned int i = 0; i < PCA4IntensityHistogram.size(); i++)
  //  for (unsigned int j = 0; j < PCA4IntensityHistogram[0].size(); j++)
  //    data(i, j) = PCA4IntensityHistogram[i][j];

  //writerv->SetFileName("E:/Projects/PSU/Data/pca4.csv");
  //writerv->SetInput(&data);
  //writerv->Write();

  //for (unsigned int i = 0; i < PCA5IntensityHistogram.size(); i++)
  //  for (unsigned int j = 0; j < PCA5IntensityHistogram[0].size(); j++)
  //    data(i, j) = PCA5IntensityHistogram[i][j];

  //writerv->SetFileName("E:/Projects/PSU/Data/pca5.csv");
  //writerv->SetInput(&data);
  //writerv->Write();

  //for (unsigned int i = 0; i < PCA6IntensityHistogram.size(); i++)
  //  for (unsigned int j = 0; j < PCA6IntensityHistogram[0].size(); j++)
  //    data(i, j) = PCA6IntensityHistogram[i][j];

  //writerv->SetFileName("E:/Projects/PSU/Data/pca6.csv");
  //writerv->SetInput(&data);
  //writerv->Write();

  //for (unsigned int i = 0; i < PCA7IntensityHistogram.size(); i++)
  //  for (unsigned int j = 0; j < PCA7IntensityHistogram[0].size(); j++)
  //    data(i, j) = PCA7IntensityHistogram[i][j];

  //writerv->SetFileName("E:/Projects/PSU/Data/pca7.csv");
  //writerv->SetInput(&data);
  //writerv->Write();
  //for (unsigned int i = 0; i < PCA8IntensityHistogram.size(); i++)
  //  for (unsigned int j = 0; j < PCA8IntensityHistogram[0].size(); j++)
  //    data(i, j) = PCA8IntensityHistogram[i][j];

  //writerv->SetFileName("E:/Projects/PSU/Data/pca8.csv");
  //writerv->SetInput(&data);
  //writerv->Write();


  //for (unsigned int i = 0; i < PCA9IntensityHistogram.size(); i++)
  //  for (unsigned int j = 0; j < PCA9IntensityHistogram[0].size(); j++)
  //    data(i, j) = PCA9IntensityHistogram[i][j];

  //writerv->SetFileName("E:/Projects/PSU/Data/pca9.csv");
  //writerv->SetInput(&data);
  //writerv->Write();


  //for (unsigned int i = 0; i < PCA10IntensityHistogram.size(); i++)
  //  for (unsigned int j = 0; j < PCA10IntensityHistogram[0].size(); j++)
  //    data(i, j) = PCA10IntensityHistogram[i][j];

  //writerv->SetFileName("E:/Projects/PSU/Data/pca10.csv");
  //writerv->SetInput(&data);
  //writerv->Write();







  try
  {
    TCReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(TCIntensityHistogram);
    std::cout << "TC" << std::endl;
  }
  catch (itk::ExceptionObject & err)
  {
    std::cout << "ExceptionObject caught !" << err << std::endl;
  }
  try
  {
    T1ReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(T1IntensityHistogram);
    std::cout << "T1" << std::endl;
  }
  catch (itk::ExceptionObject & err)
  {
    std::cout << "ExceptionObject caught !" << err << std::endl;
  }
  try
  {
    T2ReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(T2IntensityHistogram);
    std::cout << "T2" << std::endl;
  }
  catch (itk::ExceptionObject & err)
  {
    std::cout << "ExceptionObject caught !" << err << std::endl;
  }




  try
  {

    T1TCReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(T1TCIntensityHistogram);
    std::cout << "T1TC" << std::endl;
  }
  catch (itk::ExceptionObject & err)
  {
    std::cout << "ExceptionObject caught !" << err << std::endl;
  }

  try
  {

    T2FLReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(T2FLIntensityHistogram);
    std::cout << "T2FL" << std::endl;
  }
  catch (itk::ExceptionObject & err)
  {
    std::cout << "ExceptionObject caught !" << err << std::endl;
  }



  try
  {


    FLReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(FLIntensityHistogram);
    std::cout << "FL" << std::endl;
  }
  catch (itk::ExceptionObject & err)
  {
    std::cout << "ExceptionObject caught !" << err << std::endl;
  }

  try
  {


    AXReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(AXIntensityHistogram);
    std::cout << "AX" << std::endl;
  }
  catch (itk::ExceptionObject & err)
  {
    std::cout << "ExceptionObject caught !" << err << std::endl;
  }

  try
  {


    FAReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(FAIntensityHistogram);
    std::cout << "FA" << std::endl;
  }
  catch (itk::ExceptionObject & err)
  {
    std::cout << "ExceptionObject caught !" << err << std::endl;
  }


  try
  {


    RDReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(RDIntensityHistogram);
    std::cout << "RD" << std::endl;
  }
  catch (itk::ExceptionObject & err)
  {
    std::cout << "ExceptionObject caught !" << err << std::endl;
  }


  try
  {
    TRReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(TRIntensityHistogram);
    std::cout << "TR" << std::endl;
  }
  catch (itk::ExceptionObject & err)
  {
    std::cout << "ExceptionObject caught !" << err << std::endl;
  }


  try
  {
    PHReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(PHIntensityHistogram);
    std::cout << "PH" << std::endl;
  }
  catch (itk::ExceptionObject & err)
  {
    std::cout << "ExceptionObject caught !" << err << std::endl;
  }


  try
  {
    PSReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(PSIntensityHistogram);
    std::cout << "PS" << std::endl;
  }
  catch (itk::ExceptionObject & err)
  {
    std::cout << "ExceptionObject caught !" << err << std::endl;
  }

  try
  {
    RCReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(RCIntensityHistogram);
    std::cout << "RC" << std::endl;
  }
  catch (itk::ExceptionObject & err)
  {
    std::cout << "ExceptionObject caught !" << err << std::endl;
  }


  std::cout << "basic modalities perfusion components extracted" << std::endl;
  vtkSmartPointer<vtkTable> PC1ReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(PCA1IntensityHistogram);
  std::cout << "PC1" << std::endl;

  vtkSmartPointer<vtkTable> PC2ReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(PCA2IntensityHistogram);
  std::cout << "PC2" << std::endl;

  vtkSmartPointer<vtkTable> PC3ReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(PCA3IntensityHistogram);
  std::cout << "PC3" << std::endl;

  vtkSmartPointer<vtkTable> PC4ReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(PCA4IntensityHistogram);
  std::cout << "PC4" << std::endl;

  vtkSmartPointer<vtkTable> PC5ReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(PCA5IntensityHistogram);
  std::cout << "PC5" << std::endl;

  vtkSmartPointer<vtkTable> PC6ReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(PCA6IntensityHistogram);
  std::cout << "PC6" << std::endl;

  vtkSmartPointer<vtkTable> PC7ReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(PCA7IntensityHistogram);
  std::cout << "PC7" << std::endl;

  vtkSmartPointer<vtkTable> PC8ReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(PCA8IntensityHistogram);
  std::cout << "PC8" << std::endl;

  vtkSmartPointer<vtkTable> PC9ReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(PCA9IntensityHistogram);
  std::cout << "PC9" << std::endl;

  vtkSmartPointer<vtkTable> PC10ReducedIntensityHistogram = m_featureReduction.GetDiscerningPerfusionTimePointsFullPCA(PCA10IntensityHistogram);
  std::cout << "PC10" << std::endl;

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