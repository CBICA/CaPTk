#include "MolecularSubtypePredictor.h"
#include "itkCSVNumericObjectFileWriter.h"
#include "CaPTkGUIUtils.h"
#include "CaPTkClassifierUtils.h"


typedef itk::Image< float, 3 > ImageType;
//MolecularSubtypePredictor::~MolecularSubtypePredictor()
//{
//}

VectorDouble MolecularSubtypePredictor::GetStatisticalFeatures(const VectorDouble &intensities)
{
  VectorDouble StatisticalFeatures;

  double temp = 0.0;
  double mean = 0.0;
  double std = 0.0;

  for (unsigned int i = 0; i < intensities.size(); i++)
    temp = temp + intensities[i];
  mean = temp / intensities.size();

  temp = 0;
  for (unsigned int i = 0; i < intensities.size(); i++)
    temp = temp + (intensities[i] - mean)*(intensities[i] - mean);
  std = std::sqrt(temp / (intensities.size() - 1));

  StatisticalFeatures.push_back(mean);
  StatisticalFeatures.push_back(std);

  return StatisticalFeatures;
}

VectorDouble MolecularSubtypePredictor::GetHistogramFeatures(const VectorDouble &intensities, const double &start, const double &interval, const double &end)
{
  VariableLengthVectorType BinCount;
  VectorDouble finalBins;
  VectorVectorDouble Ranges;
  for (double i = start; i <= end; i = i + interval)
  {
    VectorDouble onerange;
    onerange.push_back(i - (interval / 2));
    onerange.push_back(i + (interval / 2));
    Ranges.push_back(onerange);
  }
  Ranges[Ranges.size() - 1][1] = 255;
  Ranges[0][0] = 0;


  BinCount.SetSize(Ranges.size());
  for (unsigned int j = 0; j < Ranges.size(); j++)
  {
    VectorDouble onerange = Ranges[j];
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
  return finalBins;
}

VectorDouble MolecularSubtypePredictor::GetVolumetricFeatures(const double &edemaSize, const double &tuSize, const double &neSize, const double &totalSize)
{
  VectorDouble VolumetricFeatures;
  VolumetricFeatures.push_back(tuSize);
  VolumetricFeatures.push_back(neSize);
  VolumetricFeatures.push_back(edemaSize);
  VolumetricFeatures.push_back(totalSize);

  VolumetricFeatures.push_back(tuSize + neSize);
  VolumetricFeatures.push_back((tuSize + neSize) / totalSize);
  VolumetricFeatures.push_back(edemaSize / totalSize);
  return VolumetricFeatures;
}

int MolecularSubtypePredictor::PrepareNewMolecularPredictionModel(const std::string &inputdirectory, const std::vector< std::map< CAPTK::ImageModalityType, std::string > > &qualifiedsubjects, const std::string &outputdirectory)
{
  VectorDouble AllSurvival;
  VariableSizeMatrixType FeaturesOfAllSubjects;
  VariableSizeMatrixType HistogramFeaturesConfigurations;
  HistogramFeaturesConfigurations.SetSize(33, 3); //11 modalities*3 regions = 33 configurations*3 histogram features for each configuration
  CSVFileReaderType::Pointer reader = CSVFileReaderType::New();
  MatrixType dataMatrix;
  try
  {
    reader->SetFileName("../data/survival/Survival_HMFeatures_Configuration.csv");
    reader->SetFieldDelimiterCharacter(',');
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();
    reader->Parse();
    dataMatrix = reader->GetArray2DDataObject()->GetMatrix();

    for (unsigned int i = 0; i < dataMatrix.rows(); i++)
      for (unsigned int j = 0; j < dataMatrix.cols(); j++)
        HistogramFeaturesConfigurations(i, j) = dataMatrix(i, j);
  }
  catch (const std::exception& e1)
  {
	  logger.WriteError("Cannot find the file 'Survival_HMFeatures_Configuration.csv' in the ../data/survival directory. Error code : " + std::string(e1.what()));
	  return false;
  }
  //---------------------------------------------------------------------------
  FeaturesOfAllSubjects.SetSize(qualifiedsubjects.size(), MOLECULAR_NO_OF_FEATURES);
  std::vector<std::string> patient_ids;
  for (unsigned int sid = 0; sid < qualifiedsubjects.size(); sid++)
  {
    std::cout << "Patient's data loaded:" << sid + 1 << std::endl;
    std::map< CAPTK::ImageModalityType, std::string > currentsubject = qualifiedsubjects[sid];
    patient_ids.push_back(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]));
    try
    {
      ImageType::Pointer LabelImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG]));
      ImageType::Pointer AtlasImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_ATLAS]));
      ImageType::Pointer TemplateImagePointer = ReadNiftiImage<ImageType>(getCaPTkDataDir() + "/survival/Template.nii.gz");
      ImageType::Pointer RCBVImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_RCBV])));
      ImageType::Pointer PHImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PH])));
      ImageType::Pointer T1CEImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1CE])));
      ImageType::Pointer T2FlairImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR])));
      ImageType::Pointer T1ImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1])));
      ImageType::Pointer T2ImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2])));
      ImageType::Pointer AXImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_AX])));
      ImageType::Pointer RADImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_RAD])));
      ImageType::Pointer FAImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_FA])));
      ImageType::Pointer TRImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_TR])));
      ImageType::Pointer PSRImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PSR])));

      VectorDouble ages;
      VectorDouble survival;

      reader->SetFileName(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_FEATURES]));
      reader->SetFieldDelimiterCharacter(',');
      reader->HasColumnHeadersOff();
      reader->HasRowHeadersOff();
      reader->Parse();
      dataMatrix = reader->GetArray2DDataObject()->GetMatrix();

      for (unsigned int i = 0; i < dataMatrix.rows(); i++)
      {
        ages.push_back(dataMatrix(i, 0));
        survival.push_back(dataMatrix(i, 1));
        AllSurvival.push_back(dataMatrix(i, 1));
      }

      VectorDouble TestFeatures = LoadTestData<ImageType>(T1CEImagePointer, T2FlairImagePointer, T1ImagePointer, T2ImagePointer,
        RCBVImagePointer, PSRImagePointer, PHImagePointer, AXImagePointer, FAImagePointer, RADImagePointer, TRImagePointer, LabelImagePointer, AtlasImagePointer, TemplateImagePointer, HistogramFeaturesConfigurations);

      FeaturesOfAllSubjects(sid, 0) = ages[0];
      for (unsigned int i = 1; i <= TestFeatures.size(); i++)
        FeaturesOfAllSubjects(sid, i) = TestFeatures[i - 1];
    }
    catch (const std::exception& e1)
    {
		logger.WriteError("Error in calculating the features for patient ID = " + static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PSR]) + "Error code : " + std::string(e1.what()));
		return false;
    }
  }
  std::string FeatureLabels[MOLECULAR_NO_OF_FEATURES] = { "Age","ET","NCR","ED","brain size","ET+NCR","(ET+NCR)/brain size","ED/brain size",
    "Vent_Tumor_Dist","Vent_ED_Dist","Mean_ET_T1CE","STD_ET_T1CE","Mean_ET_T1","STD_ET_T1","Mean_ET_T2","STD_ET_T2","Mean_ET_Flair","STD_ET_Flair","Mean_ET_PH","STD_ET_PH","Mean_ET_PSR","STD_ET_PSR","Mean_ET_RCBV","STD_ET_RCBV","Mean_ET_FA","STD_ET_FA","Mean_ET_AX","STD_ET_AX","Mean_ET_RAD","STD_ET_RAD","Mean_ET_ADC","STD_ET_ADC",
    "ET_T1CE_Bin1","ET_T1CE_Bin2","ED_T1CE_Bin1","ED_T1CE_Bin2","ED_T1CE_Bin3","ED_T1CE_Bin4","NCR_T1CE_Bin1","NCR_T1CE_Bin2",
    "ET_T1_Bin1","ET_T1_Bin2","ED_T1_Bin1","ED_T1_Bin2","ED_T1_Bin3","ED_T1_Bin4","ED_T1_Bin5","ED_T1_Bin6","ED_T1_Bin7","ED_T1_Bin8","ED_T1_Bin9","NCR_T1_Bin1","NCR_T1_Bin2",
    "ET_T2_Bin1","ET_T2_Bin2","ED_T2_Bin1","ED_T2_Bin2","ED_T2_Bin3","ED_T2_Bin4","ED_T2_Bin5","ED_T2_Bin6","ED_T2_Bin7","ED_T2_Bin8","ED_T2_Bin9","NCR_T2_Bin1","NCR_T2_Bin2",
    "ET_Flair_Bin1","ET_Flair_Bin2","ET_Flair_Bin3","ED_Flair_Bin1","ED_Flair_Bin2","ED_Flair_Bin3","ED_Flair_Bin4","ED_Flair_Bin5","ED_Flair_Bin6","ED_Flair_Bin7","ED_Flair_Bin8","ED_Flair_Bin9","NCR_Flair_Bin1","NCR_Flair_Bin2","NCR_Flair_Bin3","NCR_Flair_Bin4","NCR_Flair_Bin5","NCR_Flair_Bin6","NCR_Flair_Bin7",
    "ET_PH_Bin1","ET_PH_Bin2","ET_PH_Bin3","ED_PH_Bin1","ED_PH_Bin2","ED_PH_Bin3","ED_PH_Bin4","ED_PH_Bin5","ED_PH_Bin6","ED_PH_Bin7","ED_PH_Bin8","NCR_PH_Bin1","NCR_PH_Bin2","NCR_PH_Bin3",
    "ET_PSR_Bin1","ET_PSR_Bin2","ET_PSR_Bin3","ET_PSR_Bin4","ED_PSR_Bin1","ED_PSR_Bin2","ED_PSR_Bin3","ED_PSR_Bin4","ED_PSR_Bin5","ED_PSR_Bin6","ED_PSR_Bin7","ED_PSR_Bin8","ED_PSR_Bin9","ED_PSR_Bin10","NCR_PSR_Bin1","NCR_PSR_Bin2","NCR_PSR_Bin3",
    "ET_RCBV_Bin1","ET_RCBV_Bin2","ET_RCBV_Bin3","ED_RCBV_Bin1","ED_RCBV_Bin2","ED_RCBV_Bin3","NCR_RCBV_Bin1","NCR_RCBV_Bin2","NCR_RCBV_Bin3",
    "ET_FA_Bin1","ET_FA_Bin2","ED_FA_Bin1","ED_FA_Bin2","NCR_FA_Bin1","NCR_FA_Bin2",
    "ET_AX_Bin1","ET_AX_Bin2","ED_AX_Bin1","ED_AX_Bin2","NCR_AX_Bin1","NCR_AX_Bin2","NCR_AX_Bin3",
    "ET_RAD_Bin1","ET_RAD_Bin2","ED_RAD_Bin1","ED_RAD_Bin2","NCR_RAD_Bin1","NCR_RAD_Bin2","NCR_RAD_Bin3",
    "ET_ADC_Bin1","ET_ADC_Bin2","ED_ADC_Bin1","ED_ADC_Bin2","NCR_ADC_Bin1","NCR_ADC_Bin2","NCR_ADC_Bin3",
    "Frontal","Temporal","Parietal","Basal G","Insula","CC_Fornix","Occipital","Cere","Brain stem" };

  //write raw extracted features to a .csv file
  std::vector<std::string> StringFeatureLabels;
  for (int index = 0; index < MOLECULAR_NO_OF_FEATURES; index++)
    StringFeatureLabels.push_back(FeatureLabels[index]);

  WriteCSVFilesWithHorizontalAndVerticalHeaders(FeaturesOfAllSubjects, patient_ids, StringFeatureLabels, outputdirectory + "/RawFeatures.csv");

  std::cout << std::endl << "Building model....." << std::endl;
  VariableSizeMatrixType scaledFeatureSet;
  scaledFeatureSet.SetSize(qualifiedsubjects.size(), 161);
  VariableLengthVectorType meanVector;
  VariableLengthVectorType stdVector;
  mFeatureScalingLocalPtr.ScaleGivenTrainingFeatures(FeaturesOfAllSubjects, scaledFeatureSet, meanVector, stdVector);

  for (unsigned int i = 0; i < scaledFeatureSet.Rows(); i++)
    for (unsigned int j = 0; j < scaledFeatureSet.Cols(); j++)
      if (std::isnan(scaledFeatureSet(i, j)))
        scaledFeatureSet(i, j) = 0;

  WriteCSVFilesWithHorizontalAndVerticalHeaders(scaledFeatureSet, patient_ids, StringFeatureLabels, outputdirectory + "/ScaledFeatures.csv");

  //VariableSizeMatrixType ScaledFeatureSetAfterAddingAge;
  //ScaledFeatureSetAfterAddingAge.SetSize(scaledFeatureSet.Rows(), scaledFeatureSet.Cols() + 1);
  //for (unsigned int i = 0; i < scaledFeatureSet.Rows(); i++)
  //{
  //  ScaledFeatureSetAfterAddingAge(i, 0) = ages[i]; 
  //  for (unsigned int j = 0; j < scaledFeatureSet.Cols(); j++)
  //  {
  //    ScaledFeatureSetAfterAddingAge(i, j+1) = scaledFeatureSet(i, j);
  //  }
  //}
  //
  //  //readerMean->SetFileName("scaledfeatures.csv");
  //  //readerMean->SetFieldDelimiterCharacter(',');
  //  //readerMean->HasColumnHeadersOff();
  //  //readerMean->HasRowHeadersOff();
  //  //readerMean->Parse();
  //  ////typedef vnl_matrix<double> MatrixType;
  //  //dataMatrix = readerMean->GetArray2DDataObject()->GetMatrix();
  //
  //  //for (unsigned int i = 0; i < dataMatrix.rows(); i++)
  //  //  for (unsigned int j = 0; j < dataMatrix.cols(); j++)
  //  //    scaledFeatureSet(i, j) = dataMatrix(i, j);
  //  //{
  //  //  ages.push_back(dataMatrix(i, 0));
  //  //  survival.push_back(dataMatrix(i, 1));
  //  //}
  //
  ////-----------------------writing in files to compare results------------------------------
  //typedef vnl_matrix<double> MatrixType;
  //MatrixType data;
  //  //data.set_size(105, 169);
  //  //for (unsigned int i = 0; i < scaledFeatureSet.Rows(); i++)
  //  //{
  //  //  for (unsigned int j = 0; j < scaledFeatureSet.Cols(); j++)
  //  //  {
  //  //    data(i, j) = scaledFeatureSet(i, j);
  //  //  }
  //  //}
  //  //typedef itk::CSVNumericObjectFileWriter<double, 105, 169> WriterTypeMatrix;
  //  //WriterTypeMatrix::Pointer writermatrix = WriterTypeMatrix::New();
  //  //writermatrix->SetFileName(outputdirectory + "/scaledfeatures.csv");
  //  //writermatrix->SetInput(&data);
  //  //writermatrix->Write();  
  //  
  VariableSizeMatrixType proneuralModelFeatures;
  VariableSizeMatrixType neuralModelFeatures;
  VariableSizeMatrixType classicalModelFeatures;
  VariableSizeMatrixType messenchymalModelFeatures;
  mFeatureExtractionLocalPtr.FormulateMolecularTrainingData(scaledFeatureSet, AllSurvival, proneuralModelFeatures, neuralModelFeatures, messenchymalModelFeatures, classicalModelFeatures);

  try
  {
    WriteCSVFiles(meanVector, outputdirectory + "/Molecular_ZScore_Mean.csv");
    WriteCSVFiles(stdVector, outputdirectory + "/Molecular_ZScore_Std.csv");
  }
  catch (const std::exception& e1)
  {
	  logger.WriteError("Error in writing output files to the output directory = " + outputdirectory + "Error code : " + std::string(e1.what()));
	  return false;
  }

  //  //---------------------------------------------------------------------------
  // VariableSizeMatrixType SixModelSelectedFeatures = SelectSixMonthsModelFeatures(SixModelFeatures);
  //  VariableSizeMatrixType EighteenModelSelectedFeatures = SelectEighteenMonthsModelFeatures(EighteenModelFeatures);
  //MatrixType data;
  // data.set_size(30, 21);
  // for (unsigned int i = 0; i < SixModelSelectedFeatures.Rows(); i++)
  // {
  //   for (unsigned int j = 0; j < SixModelSelectedFeatures.Cols(); j++)
  //   {
  //     data(i, j) = SixModelSelectedFeatures(i, j);
  //   }
  // }
  // typedef itk::CSVNumericObjectFileWriter<double, 30,21> WriterTypeMatrix;
  // WriterTypeMatrix::Pointer writermatrix = WriterTypeMatrix::New();
  // writermatrix->SetFileName(outputdirectory + "/sixmodel.csv");
  // writermatrix->SetInput(&data);
  // writermatrix->Write();


  // for (unsigned int i = 0; i < EighteenModelSelectedFeatures.Rows(); i++)
  // {
  //   for (unsigned int j = 0; j < EighteenModelSelectedFeatures.Cols(); j++)
  //   {
  //     data(i, j) = EighteenModelSelectedFeatures(i, j);
  //   }
  // }
  // typedef itk::CSVNumericObjectFileWriter<double, 30,21> WriterTypeMatrix2;
  // WriterTypeMatrix2::Pointer writermatrix2 = WriterTypeMatrix2::New();
  // writermatrix2->SetFileName(outputdirectory + "/eighteenmodel.csv");
  // writermatrix2->SetInput(&data);
  // writermatrix2->Write();
  //
  //   //--------------------------------------read whole new data for training---------------------------------
  //VariableSizeMatrixType SixModelDataFromMatlab;
  //VariableSizeMatrixType EighteenModelDataFromMatlab;

  //readerMean->SetFileName("E:/CapTKApplications/Survival/newmatlabwork/SixModelData.csv");
  //readerMean->SetFieldDelimiterCharacter(',');
  //readerMean->HasColumnHeadersOff();
  //readerMean->HasRowHeadersOff();
  //readerMean->Parse();
  ////typedef vnl_matrix<double> MatrixType;
  //dataMatrix = readerMean->GetArray2DDataObject()->GetMatrix();
  //SixModelDataFromMatlab.SetSize(dataMatrix.rows(), dataMatrix.cols());

  //for (unsigned int i = 0; i < dataMatrix.rows(); i++)
  //{
  //  for (int j = 0; j < dataMatrix.cols(); j++)
  //    SixModelDataFromMatlab(i, j) = dataMatrix(i, j);
  //}

  //readerMean->SetFileName("E:/CapTKApplications/Survival/newmatlabwork/EighteenModelData.csv");
  //readerMean->SetFieldDelimiterCharacter(',');
  //readerMean->HasColumnHeadersOff();
  //readerMean->HasRowHeadersOff();
  //readerMean->Parse();
  ////typedef vnl_matrix<double> MatrixType;
  //dataMatrix = readerMean->GetArray2DDataObject()->GetMatrix();
  //EighteenModelDataFromMatlab.SetSize(dataMatrix.rows(), dataMatrix.cols());
  //for (unsigned int i = 0; i < dataMatrix.rows(); i++)
  //  for (int j = 0; j < dataMatrix.cols(); j++)
  //    EighteenModelDataFromMatlab(i, j) = dataMatrix(i, j);

  try
  {
    trainOpenCVSVM(neuralModelFeatures, outputdirectory + "/" + mNeuralTrainedFile, false, CAPTK::ApplicationCallingSVM::Survival);
    trainOpenCVSVM(proneuralModelFeatures, outputdirectory + "/" + mProneuralTrainedFile, false, CAPTK::ApplicationCallingSVM::Survival);
    trainOpenCVSVM(classicalModelFeatures, outputdirectory + "/" + mClassicalTrainedFile, false, CAPTK::ApplicationCallingSVM::Survival);
    trainOpenCVSVM(messenchymalModelFeatures, outputdirectory + "/" + mMessenchymalTrainedFile, false, CAPTK::ApplicationCallingSVM::Survival);
  }
  catch (const std::exception& e1)
  {
	  logger.WriteError("Training on the subjects failed. Error code : " + std::string(e1.what()));
	  return false;
  }
  std::cout << std::endl << "Model saved to the output directory." << std::endl;
  return true;
}

VariableLengthVectorType MolecularSubtypePredictor::DistanceFunction(const VariableSizeMatrixType &testData, const std::string &filename, const double &rho, const double &bestg)
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

  for (unsigned int patID = 0; patID < testData.Rows(); patID++)
  {
    double distance = 0;
    for (unsigned int svID = 0; svID < SupportVectors.Rows(); svID++)
    {
      double euclideanDistance = 0;
      for (unsigned int iterator = 0; iterator < SupportVectors.Cols(); iterator++)
        euclideanDistance = euclideanDistance + (SupportVectors(svID, iterator) - testData(patID, iterator))*(SupportVectors(svID, iterator) - testData(patID, iterator));
      double result = std::exp(-1 * bestg*euclideanDistance);
      distance = distance + result*Coefficients[svID];
    }
    Distances[patID] = distance - rho;
  }
  return Distances;
}

VectorDouble MolecularSubtypePredictor::CombineEstimates(const VectorDouble &result_proneural,
  const VectorDouble &result_neural,
  const VectorDouble &result_messenchymal,
  const VectorDouble &result_classical)
{
  VectorDouble returnVec;
  for (size_t i = 0; i < result_proneural.size(); i++)
  {
    if (result_proneural[i] > result_neural[i] && result_proneural[i] > result_messenchymal[i] && result_proneural[i] > result_classical[i])
      returnVec.push_back(CAPTK::MOLECULAR_SUBTYPES::PRONEURAL);
    else if (result_neural[i] > result_proneural[i] && result_neural[i] > result_messenchymal[i] && result_neural[i] > result_classical[i])
      returnVec.push_back(CAPTK::MOLECULAR_SUBTYPES::NEURAL);
    else if (result_messenchymal[i] > result_proneural[i] && result_messenchymal[i] > result_neural[i] && result_messenchymal[i] > result_classical[i])
      returnVec.push_back(CAPTK::MOLECULAR_SUBTYPES::MESSENCHYMAL);
    else
      returnVec.push_back(CAPTK::MOLECULAR_SUBTYPES::CLASSICAL);
  }
  return returnVec;
}


VectorDouble MolecularSubtypePredictor::MolecularSubtypePredictionOnExistingModel(const std::string &modeldirectory, const std::string &inputdirectory, const std::vector < std::map < CAPTK::ImageModalityType, std::string>> &qualifiedsubjects, const std::string &outputdirectory)
{
  VectorDouble results;
  typedef itk::CSVArray2DFileReader<double> ReaderType;
  VariableSizeMatrixType HistogramFeaturesConfigurations;
  HistogramFeaturesConfigurations.SetSize(33, 3); //11 modalities*3 regions = 33 configurations*3 histogram features for each configuration

  CSVFileReaderType::Pointer reader = CSVFileReaderType::New();
  VectorDouble ages;
  MatrixType dataMatrix;
  try
  {
    reader->SetFileName(getCaPTkDataDir() + "/survival/Survival_HMFeatures_Configuration.csv");
    reader->SetFieldDelimiterCharacter(',');
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();
    reader->Parse();
    dataMatrix = reader->GetArray2DDataObject()->GetMatrix();

    for (unsigned int i = 0; i < dataMatrix.rows(); i++)
      for (unsigned int j = 0; j < dataMatrix.cols(); j++)
        HistogramFeaturesConfigurations(i, j) = dataMatrix(i, j);
  }
  catch (const std::exception& e1)
  {
	  logger.WriteError("Cannot find the file 'Survival_HMFeatures_Configuration.csv' in the ../data/survival directory. Error code : " + std::string(e1.what()));
	  return results;
  }

  MatrixType meanMatrix;
  VariableLengthVectorType mean;
  VariableLengthVectorType stddevition;
  try
  {
    reader->SetFileName(modeldirectory + "/Molecular_ZScore_Mean.csv");
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
	  logger.WriteError("Error in reading the file: " + modeldirectory + "/Molecular_ZScore_Mean.csv. Error code : " + std::string(e1.what()));
	  return results;
  }


  MatrixType stdMatrix;
  try
  {
    reader->SetFileName(modeldirectory + "/Molecular_ZScore_Std.csv");
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
	  logger.WriteError("Error in reading the file: " + modeldirectory + "/Molecular_ZScore_Std.csv. Error code : " + std::string(e1.what()));
	  return results;
  }

  //----------------------------------------------------
  VariableSizeMatrixType FeaturesOfAllSubjects;
  FeaturesOfAllSubjects.SetSize(qualifiedsubjects.size(), MOLECULAR_NO_OF_FEATURES);
  std::vector<std::string> patient_ids;
  for (unsigned int sid = 0; sid < qualifiedsubjects.size(); sid++)
  {

    std::map<CAPTK::ImageModalityType, std::string> currentsubject = qualifiedsubjects[sid];
    patient_ids.push_back(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]));
    cout << "Subject No:" << sid << " Name: "<<patient_ids[sid]<< endl;
    try
	{
    ImageType::Pointer LabelImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG]));
    ImageType::Pointer AtlasImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_ATLAS]));
    ImageType::Pointer TemplateImagePointer = ReadNiftiImage<ImageType>(getCaPTkDataDir() + "/survival/Template.nii.gz");
    ImageType::Pointer RCBVImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_RCBV])));
    ImageType::Pointer PHImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PH])));
    ImageType::Pointer T1CEImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1CE])));
    ImageType::Pointer T2FlairImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR])));
    ImageType::Pointer T1ImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1])));
    ImageType::Pointer T2ImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2])));
    ImageType::Pointer AXImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_AX])));
    ImageType::Pointer RADImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_RAD])));
    ImageType::Pointer FAImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_FA])));
    ImageType::Pointer TRImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_TR])));
    ImageType::Pointer PSRImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PSR])));

    VectorDouble TestFeatures = LoadTestData<ImageType>(T1CEImagePointer, T2FlairImagePointer, T1ImagePointer, T2ImagePointer,
      RCBVImagePointer, PSRImagePointer, PHImagePointer, AXImagePointer, FAImagePointer, RADImagePointer, TRImagePointer, LabelImagePointer, AtlasImagePointer, TemplateImagePointer, HistogramFeaturesConfigurations);

    double age;

      reader->SetFileName(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_FEATURES]));
      reader->SetFieldDelimiterCharacter(',');
      reader->HasColumnHeadersOff();
      reader->HasRowHeadersOff();
      reader->Parse();
      dataMatrix = reader->GetArray2DDataObject()->GetMatrix();

      for (unsigned int i = 0; i < dataMatrix.rows(); i++)
        age = dataMatrix(i, 0);
    FeaturesOfAllSubjects(sid, 0) = age;
    for (unsigned int i = 1; i <= TestFeatures.size(); i++)
      FeaturesOfAllSubjects(sid, i) = TestFeatures[i - 1];

	}
	catch (const std::exception& e1)
	{
		logger.WriteError("Error in calculating the features for patient ID = " + static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PSR]) + "Error code : " + std::string(e1.what()));
		return results;
	}

  }

  std::string FeatureLabels[MOLECULAR_NO_OF_FEATURES] = { "Age","ET","NCR","ED","brain size","ET+NCR","(ET+NCR)/brain size","ED/brain size",
    "Vent_Tumor_Dist","Vent_ED_Dist","Mean_ET_T1CE","STD_ET_T1CE","Mean_ET_T1","STD_ET_T1","Mean_ET_T2","STD_ET_T2","Mean_ET_Flair","STD_ET_Flair","Mean_ET_PH","STD_ET_PH","Mean_ET_PSR","STD_ET_PSR","Mean_ET_RCBV","STD_ET_RCBV","Mean_ET_FA","STD_ET_FA","Mean_ET_AX","STD_ET_AX","Mean_ET_RAD","STD_ET_RAD","Mean_ET_ADC","STD_ET_ADC",
    "ET_T1CE_Bin1","ET_T1CE_Bin2","ED_T1CE_Bin1","ED_T1CE_Bin2","ED_T1CE_Bin3","ED_T1CE_Bin4","NCR_T1CE_Bin1","NCR_T1CE_Bin2",
    "ET_T1_Bin1","ET_T1_Bin2","ED_T1_Bin1","ED_T1_Bin2","ED_T1_Bin3","ED_T1_Bin4","ED_T1_Bin5","ED_T1_Bin6","ED_T1_Bin7","ED_T1_Bin8","ED_T1_Bin9","NCR_T1_Bin1","NCR_T1_Bin2",
    "ET_T2_Bin1","ET_T2_Bin2","ED_T2_Bin1","ED_T2_Bin2","ED_T2_Bin3","ED_T2_Bin4","ED_T2_Bin5","ED_T2_Bin6","ED_T2_Bin7","ED_T2_Bin8","ED_T2_Bin9","NCR_T2_Bin1","NCR_T2_Bin2",
    "ET_Flair_Bin1","ET_Flair_Bin2","ET_Flair_Bin3","ED_Flair_Bin1","ED_Flair_Bin2","ED_Flair_Bin3","ED_Flair_Bin4","ED_Flair_Bin5","ED_Flair_Bin6","ED_Flair_Bin7","ED_Flair_Bin8","ED_Flair_Bin9","NCR_Flair_Bin1","NCR_Flair_Bin2","NCR_Flair_Bin3","NCR_Flair_Bin4","NCR_Flair_Bin5","NCR_Flair_Bin6","NCR_Flair_Bin7",
    "ET_PH_Bin1","ET_PH_Bin2","ET_PH_Bin3","ED_PH_Bin1","ED_PH_Bin2","ED_PH_Bin3","ED_PH_Bin4","ED_PH_Bin5","ED_PH_Bin6","ED_PH_Bin7","ED_PH_Bin8","NCR_PH_Bin1","NCR_PH_Bin2","NCR_PH_Bin3",
    "ET_PSR_Bin1","ET_PSR_Bin2","ET_PSR_Bin3","ET_PSR_Bin4","ED_PSR_Bin1","ED_PSR_Bin2","ED_PSR_Bin3","ED_PSR_Bin4","ED_PSR_Bin5","ED_PSR_Bin6","ED_PSR_Bin7","ED_PSR_Bin8","ED_PSR_Bin9","ED_PSR_Bin10","NCR_PSR_Bin1","NCR_PSR_Bin2","NCR_PSR_Bin3",
    "ET_RCBV_Bin1","ET_RCBV_Bin2","ET_RCBV_Bin3","ED_RCBV_Bin1","ED_RCBV_Bin2","ED_RCBV_Bin3","NCR_RCBV_Bin1","NCR_RCBV_Bin2","NCR_RCBV_Bin3",
    "ET_FA_Bin1","ET_FA_Bin2","ED_FA_Bin1","ED_FA_Bin2","NCR_FA_Bin1","NCR_FA_Bin2",
    "ET_AX_Bin1","ET_AX_Bin2","ED_AX_Bin1","ED_AX_Bin2","NCR_AX_Bin1","NCR_AX_Bin2","NCR_AX_Bin3",
    "ET_RAD_Bin1","ET_RAD_Bin2","ED_RAD_Bin1","ED_RAD_Bin2","NCR_RAD_Bin1","NCR_RAD_Bin2","NCR_RAD_Bin3",
    "ET_ADC_Bin1","ET_ADC_Bin2","ED_ADC_Bin1","ED_ADC_Bin2","NCR_ADC_Bin1","NCR_ADC_Bin2","NCR_ADC_Bin3",
    "Frontal","Temporal","Parietal","Basal G","Insula","CC_Fornix","Occipital","Cere","Brain stem" };

  //write raw extracted features to a .csv file
  std::vector<std::string> StringFeatureLabels;
  for (int index = 0; index < MOLECULAR_NO_OF_FEATURES; index++)
    StringFeatureLabels.push_back(FeatureLabels[index]);

  WriteCSVFilesWithHorizontalAndVerticalHeaders(FeaturesOfAllSubjects, patient_ids, StringFeatureLabels, outputdirectory + "/RawFeatures.csv");
  
  VariableSizeMatrixType ScaledTestingData = mFeatureScalingLocalPtr.ScaleGivenTestingFeatures(FeaturesOfAllSubjects, mean, stddevition);
  for (unsigned int i = 0; i < ScaledTestingData.Rows(); i++)
    for (unsigned int j = 0; j < ScaledTestingData.Cols(); j++)
      if (std::isnan(ScaledTestingData(i, j)))
        ScaledTestingData(i, j) = 0;

  WriteCSVFilesWithHorizontalAndVerticalHeaders(ScaledTestingData, patient_ids, StringFeatureLabels, outputdirectory + "/ScaledFeatures.csv");

  //VariableSizeMatrixType ScaledFeatureSetAfterAddingLabel;
  //ScaledFeatureSetAfterAddingLabel.SetSize(ScaledTestingData.Rows(), ScaledTestingData.Cols() + 1);
  //for (unsigned int i = 0; i < ScaledTestingData.Rows(); i++)
  //{
  //  unsigned int j = 0;
  //  for (j = 0; j < ScaledTestingData.Cols(); j++)
  //    ScaledFeatureSetAfterAddingLabel(i, j) = ScaledTestingData(i, j);
  //  ScaledFeatureSetAfterAddingLabel(i, j) = 0;
  //}

  //VariableSizeMatrixType SixModelSelectedFeatures = SelectSixMonthsModelFeatures(ScaledFeatureSetAfterAddingLabel);
  //VariableSizeMatrixType EighteenModelSelectedFeatures = SelectEighteenMonthsModelFeatures(ScaledFeatureSetAfterAddingLabel);

  try
  {
    std::ofstream myfile;
    myfile.open(outputdirectory + "/results.csv");
    myfile << "Subject Name,Predicted Subtype \n";

   if (cbica::fileExists(modeldirectory + "/" + mProneuralTrainedFile) == true && cbica::fileExists(modeldirectory + "/" + mNeuralTrainedFile) == true && cbica::fileExists(modeldirectory + "/" + mMessenchymalTrainedFile) == true && cbica::fileExists(modeldirectory + "/" + mClassicalTrainedFile) == true)
   {
      VectorDouble result_proneural;
      VectorDouble result_neural;
      VectorDouble result_classical;
      VectorDouble result_messenchymal;
      result_proneural = testOpenCVSVM(ScaledTestingData, modeldirectory + "/" + mProneuralTrainedFile);
      result_neural = testOpenCVSVM(ScaledTestingData, modeldirectory + "/" + mNeuralTrainedFile);
      result_messenchymal = testOpenCVSVM(ScaledTestingData, modeldirectory + "/" + mMessenchymalTrainedFile);
      result_classical = testOpenCVSVM(ScaledTestingData, modeldirectory + "/" + mClassicalTrainedFile);

      results = CombineEstimates(result_proneural, result_neural, result_messenchymal, result_classical);
      for (size_t i = 0; i < results.size(); i++)
      {
        std::map<CAPTK::ImageModalityType, std::string> currentsubject = qualifiedsubjects[i];
        std::string subtype;
        if (results[i] == 1)
          subtype = "Proneural";
        else if (results[i] == 2)
          subtype = "Neural";
        else if (results[i] == 3)
          subtype = "Mesenchymal";
        else if (results[i] == 4)
          subtype = "Classical";
        else
          subtype = "Unknown";
        myfile << static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]) + "," + subtype + "\n";
      }
    }

    myfile.close();
  }
  catch (itk::ExceptionObject & excp)
  {
	  logger.WriteError("Error caught during testing: " + std::string(excp.GetDescription()));
	  return results;
  }
  return results;

}

VariableSizeMatrixType MolecularSubtypePredictor::SelectSixMonthsModelFeatures(const VariableSizeMatrixType &SixModelFeatures)
{
  int selectedFeatures[20] = { 1, 5, 9, 10, 20, 23, 24, 37, 38, 43, 44, 48, 49, 50, 51, 56, 57, 61, 62, 63 };
  for (unsigned int i = 0; i <20; i++)
    selectedFeatures[i] = selectedFeatures[i] - 1;
  VariableSizeMatrixType SixModelSelectedFeatures;
  SixModelSelectedFeatures.SetSize(SixModelFeatures.Rows(), 21);
  int counter = 0;
  for (unsigned int i = 0; i < 20; i++)
  {
    for (unsigned int j = 0; j < SixModelFeatures.Rows(); j++)
      SixModelSelectedFeatures(j, counter) = SixModelFeatures(j, selectedFeatures[i]);
    counter++;
  }
  for (unsigned int j = 0; j < SixModelFeatures.Rows(); j++)
    SixModelSelectedFeatures(j, 20) = SixModelFeatures(j, 161);

  return SixModelSelectedFeatures;
}

VariableSizeMatrixType MolecularSubtypePredictor::SelectEighteenMonthsModelFeatures(const VariableSizeMatrixType &EighteenModelFeatures)
{
  int selectedFeatures[20] = { 1, 5, 10, 15, 24, 27, 37, 38, 50, 51, 53, 62, 63, 64, 67, 70, 71, 85, 158, 159 };
  for (unsigned int i = 0; i<20; i++)
    selectedFeatures[i] = selectedFeatures[i] - 1;

  VariableSizeMatrixType EighteenModelSelectedFeatures;
  EighteenModelSelectedFeatures.SetSize(EighteenModelFeatures.Rows(), 21);
  int counter = 0;
  for (unsigned int i = 0; i < 20; i++)
  {
    for (unsigned int j = 0; j < EighteenModelFeatures.Rows(); j++)
      EighteenModelSelectedFeatures(j, counter) = EighteenModelFeatures(j, selectedFeatures[i]);
    counter++;
  }
  for (unsigned int j = 0; j < EighteenModelFeatures.Rows(); j++)
    EighteenModelSelectedFeatures(j, 20) = EighteenModelFeatures(j, 161);
  return EighteenModelSelectedFeatures;
}
