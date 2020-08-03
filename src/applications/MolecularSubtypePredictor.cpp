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

  VariableSizeMatrixType scaledFeatureSet;
  scaledFeatureSet.SetSize(qualifiedsubjects.size(), MOLECULAR_NO_OF_FEATURES);
  VariableLengthVectorType meanVector;
  VariableLengthVectorType stdVector;
  mFeatureScalingLocalPtr.ScaleGivenTrainingFeatures(FeaturesOfAllSubjects, scaledFeatureSet, meanVector, stdVector);

  for (unsigned int i = 0; i < scaledFeatureSet.Rows(); i++)
    for (unsigned int j = 0; j < scaledFeatureSet.Cols(); j++)
      if (std::isnan(scaledFeatureSet(i, j)))
        scaledFeatureSet(i, j) = 0;

  try
  {
    WriteCSVFilesWithHorizontalAndVerticalHeaders(scaledFeatureSet, patient_ids, StringFeatureLabels, outputdirectory + "/ScaledFeatures.csv");
    WriteCSVFiles(meanVector, outputdirectory + "/Molecular_ZScore_Mean.csv",true);
    WriteCSVFiles(stdVector, outputdirectory + "/Molecular_ZScore_Std.csv",true);
  }
  catch (const std::exception& e1)
  {
	  logger.WriteError("Error in writing output files to the output directory = " + outputdirectory + "Error code : " + std::string(e1.what()));
	  return false;
  }

  VectorDouble proneuralModelLabels, neuralModelLabels, messModelLabels, classicalModelLabels;
  mFeatureExtractionLocalPtr.FormulateMolecularTrainingData(AllSurvival, proneuralModelLabels, neuralModelLabels, messModelLabels, classicalModelLabels);

  //select model features using routines of training module
  TrainingModule mTrainingModule;
  VectorDouble selectedfeatures_Proneural, selectedfeatures_Neural, selectedfeatures_Mess, selectedfeatures_Classical;
  
  VectorDouble EffectSize = mTrainingModule.EffectSizeFeatureSelection(scaledFeatureSet, proneuralModelLabels);
  for (unsigned int eSizeCounter = 0; eSizeCounter < EffectSize.size(); eSizeCounter++)
    EffectSize[eSizeCounter] = std::abs(EffectSize[eSizeCounter]);
  std::vector<size_t> indices = mTrainingModule.sort_indexes(EffectSize);
  for (int index = 0; index < indices.size()*0.2; index++)
    selectedfeatures_Proneural.push_back(indices[index]);

  EffectSize = mTrainingModule.EffectSizeFeatureSelection(scaledFeatureSet, neuralModelLabels);
  for (unsigned int eSizeCounter = 0; eSizeCounter < EffectSize.size(); eSizeCounter++)
    EffectSize[eSizeCounter] = std::abs(EffectSize[eSizeCounter]);
  indices = mTrainingModule.sort_indexes(EffectSize);
  for (int index = 0; index < indices.size()*0.2; index++)
    selectedfeatures_Neural.push_back(indices[index]);

  EffectSize = mTrainingModule.EffectSizeFeatureSelection(scaledFeatureSet, messModelLabels);
  for (unsigned int eSizeCounter = 0; eSizeCounter < EffectSize.size(); eSizeCounter++)
    EffectSize[eSizeCounter] = std::abs(EffectSize[eSizeCounter]);
  indices = mTrainingModule.sort_indexes(EffectSize);
  for (int index = 0; index < indices.size()*0.2; index++)
    selectedfeatures_Mess.push_back(indices[index]);

  EffectSize = mTrainingModule.EffectSizeFeatureSelection(scaledFeatureSet, classicalModelLabels);
  for (unsigned int eSizeCounter = 0; eSizeCounter < EffectSize.size(); eSizeCounter++)
    EffectSize[eSizeCounter] = std::abs(EffectSize[eSizeCounter]);
  indices = mTrainingModule.sort_indexes(EffectSize);
  for (int index = 0; index < indices.size()*0.2; index++)
    selectedfeatures_Classical.push_back(indices[index]);

  WriteCSVFiles(selectedfeatures_Proneural, outputdirectory + "/Molecular_SelectedFeatures_Proneural.csv", true);
  WriteCSVFiles(selectedfeatures_Neural, outputdirectory + "/Molecular_SelectedFeatures_Neural.csv", true);
  WriteCSVFiles(selectedfeatures_Mess, outputdirectory + "/Molecular_SelectedFeatures_Messenchymal.csv", true);  
  WriteCSVFiles(selectedfeatures_Classical, outputdirectory + "/Molecular_SelectedFeatures_Classical.csv", true);

  VariableSizeMatrixType ProneuralSelectedFeatures = SelectModelFeatures(scaledFeatureSet, selectedfeatures_Proneural);
  VariableSizeMatrixType NeuralSelectedFeatures = SelectModelFeatures(scaledFeatureSet, selectedfeatures_Neural);
  VariableSizeMatrixType MessSelectedFeatures = SelectModelFeatures(scaledFeatureSet, selectedfeatures_Mess);
  VariableSizeMatrixType ClassicalSelectedFeatures = SelectModelFeatures(scaledFeatureSet, selectedfeatures_Classical);

  //writing selected model features
  std::vector<std::string> SelectedFeatureLabels_Proneural, SelectedFeatureLabels_Neural, SelectedFeatureLabels_Mess, SelectedFeatureLabels_Classical;
  for (int index = 0; index < selectedfeatures_Proneural.size(); index++)
  {
    int currentindex = selectedfeatures_Proneural[index];
    SelectedFeatureLabels_Proneural.push_back(FeatureLabels[currentindex]);
  }
  for (int index = 0; index < selectedfeatures_Neural.size(); index++)
  {
    int currentindex = selectedfeatures_Neural[index];
    SelectedFeatureLabels_Neural.push_back(FeatureLabels[currentindex]);
  }
  for (int index = 0; index < selectedfeatures_Mess.size(); index++)
  {
    int currentindex = selectedfeatures_Mess[index];
    SelectedFeatureLabels_Mess.push_back(FeatureLabels[currentindex]);
  }
  for (int index = 0; index < selectedfeatures_Classical.size(); index++)
  {
    int currentindex = selectedfeatures_Classical[index];
    SelectedFeatureLabels_Classical.push_back(FeatureLabels[currentindex]);
  }
  WriteCSVFilesWithHorizontalAndVerticalHeaders(ProneuralSelectedFeatures, patient_ids, SelectedFeatureLabels_Proneural, outputdirectory + "/SelectedFeatures_Proneural.csv");
  WriteCSVFilesWithHorizontalAndVerticalHeaders(NeuralSelectedFeatures, patient_ids, SelectedFeatureLabels_Neural, outputdirectory + "/SelectedFeatures_Neural.csv");
  WriteCSVFilesWithHorizontalAndVerticalHeaders(MessSelectedFeatures, patient_ids, SelectedFeatureLabels_Mess, outputdirectory + "/SelectedFeatures_Messenchymal.csv");
  WriteCSVFilesWithHorizontalAndVerticalHeaders(ClassicalSelectedFeatures, patient_ids, SelectedFeatureLabels_Classical, outputdirectory + "/SelectedFeatures_Classical.csv");

  //append labels to model features as the training function expects labels in the last column
  VariableSizeMatrixType finaldatamatrix_Proneural, finaldatamatrix_Neural, finaldatamatrix_Mess, finaldatamatrix_Classical;
  finaldatamatrix_Proneural.SetSize(ProneuralSelectedFeatures.Rows(), ProneuralSelectedFeatures.Cols() + 1);
  finaldatamatrix_Neural.SetSize(NeuralSelectedFeatures.Rows(), NeuralSelectedFeatures.Cols() + 1);
  finaldatamatrix_Mess.SetSize(MessSelectedFeatures.Rows(), MessSelectedFeatures.Cols() + 1);
  finaldatamatrix_Classical.SetSize(ClassicalSelectedFeatures.Rows(), ClassicalSelectedFeatures.Cols() + 1);

  for (unsigned int i = 0; i < finaldatamatrix_Proneural.Rows(); i++)
  {
    for (unsigned int j = 0; j < finaldatamatrix_Proneural.Cols() - 1; j++)
      finaldatamatrix_Proneural(i, j) = ProneuralSelectedFeatures(i, j);
    finaldatamatrix_Proneural(i, finaldatamatrix_Proneural.Cols() - 1) = proneuralModelLabels[i];
  }
  for (unsigned int i = 0; i < finaldatamatrix_Neural.Rows(); i++)
  {
    for (unsigned int j = 0; j < finaldatamatrix_Neural.Cols() - 1; j++)
      finaldatamatrix_Neural(i, j) = NeuralSelectedFeatures(i, j);
    finaldatamatrix_Neural(i, finaldatamatrix_Neural.Cols() - 1) = neuralModelLabels[i];
  }
  for (unsigned int i = 0; i < finaldatamatrix_Mess.Rows(); i++)
  {
    for (unsigned int j = 0; j < finaldatamatrix_Mess.Cols() - 1; j++)
      finaldatamatrix_Mess(i, j) = MessSelectedFeatures(i, j);
    finaldatamatrix_Mess(i, finaldatamatrix_Mess.Cols() - 1) = messModelLabels[i];
  }
  for (unsigned int i = 0; i < finaldatamatrix_Classical.Rows(); i++)
  {
    for (unsigned int j = 0; j < finaldatamatrix_Classical.Cols() - 1; j++)
      finaldatamatrix_Classical(i, j) = ClassicalSelectedFeatures(i, j);
    finaldatamatrix_Classical(i, finaldatamatrix_Classical.Cols() - 1) = classicalModelLabels[i];
  }
  std::cout << std::endl << "Building model....." << std::endl;
  try
  {
    trainOpenCVSVM(finaldatamatrix_Proneural, outputdirectory + "/" + mNeuralTrainedFile, false, CAPTK::ApplicationCallingSVM::Survival);
    trainOpenCVSVM(finaldatamatrix_Neural, outputdirectory + "/" + mProneuralTrainedFile, false, CAPTK::ApplicationCallingSVM::Survival);
    trainOpenCVSVM(finaldatamatrix_Mess, outputdirectory + "/" + mClassicalTrainedFile, false, CAPTK::ApplicationCallingSVM::Survival);
    trainOpenCVSVM(finaldatamatrix_Classical, outputdirectory + "/" + mMessenchymalTrainedFile, false, CAPTK::ApplicationCallingSVM::Survival);
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
  VariableLengthVectorType features_Proneural, features_Neural, features_Mess, features_Classical;
  try
  {
    reader->SetFileName(modeldirectory + "/Molecular_SelectedFeatures_Proneural.csv");
    reader->SetFieldDelimiterCharacter(',');
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();
    reader->Parse();
    dataMatrix = reader->GetArray2DDataObject()->GetMatrix();
    features_Proneural.SetSize(dataMatrix.size());
    for (unsigned int i = 0; i < dataMatrix.size(); i++)
      features_Proneural[i] = dataMatrix(i, 0);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Error in reading the file: " + modeldirectory + "/Molecular_SelectedFeatures_Proneural.csv. Error code : " + std::string(e1.what()));
    return results;
  }
  try
  {
    reader->SetFileName(modeldirectory + "/Molecular_SelectedFeatures_Neural.csv");
    reader->SetFieldDelimiterCharacter(',');
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();
    reader->Parse();
    dataMatrix = reader->GetArray2DDataObject()->GetMatrix();
    features_Neural.SetSize(dataMatrix.size());
    for (unsigned int i = 0; i < dataMatrix.size(); i++)
      features_Neural[i] = dataMatrix(i, 0);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Error in reading the file: " + modeldirectory + "/Molecular_SelectedFeatures_Neural.csv. Error code : " + std::string(e1.what()));
    return results;
  }
  try
  {
    reader->SetFileName(modeldirectory + "/Molecular_SelectedFeatures_Messenchymal.csv");
    reader->SetFieldDelimiterCharacter(',');
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();
    reader->Parse();
    dataMatrix = reader->GetArray2DDataObject()->GetMatrix();
    features_Mess.SetSize(dataMatrix.size());
    for (unsigned int i = 0; i < dataMatrix.size(); i++)
      features_Mess[i] = dataMatrix(i, 0);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Error in reading the file: " + modeldirectory + "/Molecular_SelectedFeatures_Messenchymal.csv. Error code : " + std::string(e1.what()));
    return results;
  }
  try
  {
    reader->SetFileName(modeldirectory + "/Molecular_SelectedFeatures_Classical.csv");
    reader->SetFieldDelimiterCharacter(',');
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();
    reader->Parse();
    dataMatrix = reader->GetArray2DDataObject()->GetMatrix();
    features_Classical.SetSize(dataMatrix.size());
    for (unsigned int i = 0; i < dataMatrix.size(); i++)
      features_Classical[i] = dataMatrix(i, 0);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Error in reading the file: " + modeldirectory + "/Molecular_SelectedFeatures_Classical.csv. Error code : " + std::string(e1.what()));
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

  VariableSizeMatrixType ProneuralSelectedFeatures = SelectModelFeatures(ScaledTestingData, features_Proneural);
  VariableSizeMatrixType NeuralSelectedFeatures = SelectModelFeatures(ScaledTestingData, features_Neural);
  VariableSizeMatrixType MessSelectedFeatures = SelectModelFeatures(ScaledTestingData, features_Mess);
  VariableSizeMatrixType ClassicalSelectedFeatures = SelectModelFeatures(ScaledTestingData, features_Classical);

  //writing selected model features
  std::vector<std::string> SelectedFeatureLabels_Proneural, SelectedFeatureLabels_Neural, SelectedFeatureLabels_Mess, SelectedFeatureLabels_Classical;
  for (int index = 0; index < features_Proneural.Size(); index++)
  {
    int currentindex = features_Proneural[index];
    SelectedFeatureLabels_Proneural.push_back(FeatureLabels[currentindex]);
  }
  for (int index = 0; index < features_Neural.Size(); index++)
  {
    int currentindex = features_Neural[index];
    SelectedFeatureLabels_Neural.push_back(FeatureLabels[currentindex]);
  }
  for (int index = 0; index < features_Mess.Size(); index++)
  {
    int currentindex = features_Mess[index];
    SelectedFeatureLabels_Mess.push_back(FeatureLabels[currentindex]);
  }
  for (int index = 0; index < features_Classical.Size(); index++)
  {
    int currentindex = features_Classical[index];
    SelectedFeatureLabels_Classical.push_back(FeatureLabels[currentindex]);
  }
  WriteCSVFilesWithHorizontalAndVerticalHeaders(ProneuralSelectedFeatures, patient_ids, SelectedFeatureLabels_Proneural, outputdirectory + "/SelectedFeatures_Proneural.csv");
  WriteCSVFilesWithHorizontalAndVerticalHeaders(NeuralSelectedFeatures, patient_ids, SelectedFeatureLabels_Neural, outputdirectory + "/SelectedFeatures_Neural.csv");
  WriteCSVFilesWithHorizontalAndVerticalHeaders(MessSelectedFeatures, patient_ids, SelectedFeatureLabels_Mess, outputdirectory + "/SelectedFeatures_Messenchymal.csv");
  WriteCSVFilesWithHorizontalAndVerticalHeaders(ClassicalSelectedFeatures, patient_ids, SelectedFeatureLabels_Classical, outputdirectory + "/SelectedFeatures_Classical.csv");

  try
  {
    std::ofstream myfile;
    myfile.open(outputdirectory + "/results.csv");
    myfile << "Subject Name,Predicted Subtype \n";

   if (cbica::fileExists(modeldirectory + "/" + mProneuralTrainedFile) == true && cbica::fileExists(modeldirectory + "/" + mNeuralTrainedFile) == true && cbica::fileExists(modeldirectory + "/" + mMessenchymalTrainedFile) == true && cbica::fileExists(modeldirectory + "/" + mClassicalTrainedFile) == true)
   {
      VectorDouble result_proneural    = testOpenCVSVM(ProneuralSelectedFeatures, modeldirectory + "/" + mProneuralTrainedFile);
      VectorDouble result_neural       = testOpenCVSVM(NeuralSelectedFeatures, modeldirectory + "/" + mNeuralTrainedFile);
      VectorDouble result_messenchymal = testOpenCVSVM(MessSelectedFeatures, modeldirectory + "/" + mMessenchymalTrainedFile);
      VectorDouble result_classical    = testOpenCVSVM(ClassicalSelectedFeatures, modeldirectory + "/" + mClassicalTrainedFile);

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


VariableSizeMatrixType MolecularSubtypePredictor::SelectModelFeatures(const VariableSizeMatrixType &ModelFeatures, const VariableLengthVectorType &selectedFeatures)
{
  VariableSizeMatrixType ModelSelectedFeatures;
  //make a feature matrix to store selected features. rows= rows of input features, columns=number of selected features 
  ModelSelectedFeatures.SetSize(ModelFeatures.Rows(), selectedFeatures.Size());
  int counter = 0;
  //copy selected features
  for (unsigned int i = 0; i < selectedFeatures.Size(); i++)
  {
    for (unsigned int j = 0; j < ModelFeatures.Rows(); j++)
      ModelSelectedFeatures(j, counter) = ModelFeatures(j, selectedFeatures[i]);
    counter++;
  }
  return ModelSelectedFeatures;
}

VariableSizeMatrixType MolecularSubtypePredictor::SelectModelFeatures(const VariableSizeMatrixType &ModelFeatures, const VectorDouble &selectedFeatures)
{
  VariableSizeMatrixType ModelSelectedFeatures;
  //make a feature matrix to store selected features. rows= rows of input features, columns=number of selected features 
  ModelSelectedFeatures.SetSize(ModelFeatures.Rows(), selectedFeatures.size());
  int counter = 0;
  //copy selected features
  for (unsigned int i = 0; i < selectedFeatures.size(); i++)
  {
    for (unsigned int j = 0; j < ModelFeatures.Rows(); j++)
      ModelSelectedFeatures(j, counter) = ModelFeatures(j, selectedFeatures[i]);
    counter++;
  }
  return ModelSelectedFeatures;
}
