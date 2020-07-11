#include "SurvivalPredictor.h"
#include "itkCSVNumericObjectFileWriter.h"
#include "CaPTkGUIUtils.h"
#include "CaPTkClassifierUtils.h"
#include "fMainWindow.h"
#include "cbicaStatistics.h"

VectorDouble SurvivalPredictor::GetStatisticalFeatures(const VectorDouble &intensities)
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

VectorDouble SurvivalPredictor::GetHistogramFeatures(const VectorDouble &intensities, const double &start, const double &interval, const double &end)
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

VectorDouble SurvivalPredictor::GetVolumetricFeatures(const double &edemaSize,const double &tuSize, const double &neSize, const double &totalSize)
{
  //calculate volume of enhancing tumor, non-enhancing tumor, edema, whole tumor, and tumor core
  VectorDouble VolumetricFeatures;
  VolumetricFeatures.push_back(tuSize);
  VolumetricFeatures.push_back(neSize);
  VolumetricFeatures.push_back(edemaSize);
  VolumetricFeatures.push_back(totalSize);
  VolumetricFeatures.push_back(tuSize + neSize);
  
  //calculate ratios of different tumor sub-regions
  VolumetricFeatures.push_back( (tuSize + neSize)*100 / totalSize);
  VolumetricFeatures.push_back(edemaSize*100 / totalSize);
  VolumetricFeatures.push_back(tuSize*100 / (tuSize + neSize));
  VolumetricFeatures.push_back(neSize*100 / (tuSize + neSize));
  VolumetricFeatures.push_back(neSize*100 / (tuSize + neSize));

  return VolumetricFeatures;
}

int SurvivalPredictor::TrainNewSurvivalPredictionModel(const std::string &inputdirectory,const std::vector< std::map< CAPTK::ImageModalityType, std::string > > &qualifiedsubjects, const std::string &outputdirectory)
{
  VectorDouble AllSurvival; 
  VariableSizeMatrixType FeaturesOfAllSubjects;
  VariableSizeMatrixType HistogramFeaturesConfigurations;
  HistogramFeaturesConfigurations.SetSize(33, 3); //11 modalities*3 regions = 33 configurations*3 histogram features for each configuration
  CSVFileReaderType::Pointer reader = CSVFileReaderType::New();
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
	  return false;
  }
  FeaturesOfAllSubjects.SetSize(qualifiedsubjects.size(), SURVIVAL_NO_OF_FEATURES);
  std::vector<std::string> patient_ids;
  for (unsigned int sid = 0; sid < qualifiedsubjects.size(); sid++)
  {
	  std::cout << "Patient's data loaded:" << sid + 1 << std::endl;
	  std::map< CAPTK::ImageModalityType, std::string > currentsubject = qualifiedsubjects[sid];
    patient_ids.push_back(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]));
	  try
	  {
		  ImageType::Pointer LabelImagePointer = cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG]));
		  ImageType::Pointer AtlasImagePointer = cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_ATLAS]));
      ImageType::Pointer TemplateImagePointer = cbica::ReadImage<ImageType>(getCaPTkDataDir() + "/survival/Template.nii.gz");
		  ImageType::Pointer RCBVImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_RCBV])));
		  ImageType::Pointer PHImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PH])));
		  ImageType::Pointer T1CEImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1CE])));
		  ImageType::Pointer T2FlairImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR])));
		  ImageType::Pointer T1ImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1])));
		  ImageType::Pointer T2ImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2])));
		  ImageType::Pointer AXImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_AX])));
		  ImageType::Pointer RADImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_RAD])));
		  ImageType::Pointer FAImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_FA])));
		  ImageType::Pointer TRImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_TR])));
		  ImageType::Pointer PSRImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PSR])));

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
		  VectorDouble TestFeatures = CalculateFeatures<ImageType>(T1CEImagePointer, T2FlairImagePointer, T1ImagePointer, T2ImagePointer,
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
  std::string FeatureLabels[SURVIVAL_NO_OF_FEATURES] = { "Age","ET","NCR","ED","brain size","ET+NCR","(ET+NCR)x100/brain size","EDx100/brain size","ETx100/(ET+NCR)","NCRx100/(ET+NCR)","NCRx100/(ET+NCR)",
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
  for (int index = 0; index < SURVIVAL_NO_OF_FEATURES; index++)
    StringFeatureLabels.push_back(FeatureLabels[index]);

  WriteCSVFilesWithHorizontalAndVerticalHeaders(FeaturesOfAllSubjects, patient_ids, StringFeatureLabels, outputdirectory + "/RawFeatures.csv");

  VariableSizeMatrixType scaledFeatureSet;
  scaledFeatureSet.SetSize(qualifiedsubjects.size(), SURVIVAL_NO_OF_FEATURES);
  VariableLengthVectorType meanVector;
  VariableLengthVectorType stdVector;
  mFeatureScalingLocalPtr.ScaleGivenTrainingFeatures(FeaturesOfAllSubjects, scaledFeatureSet, meanVector, stdVector);

  for (unsigned int i = 0; i < scaledFeatureSet.Rows(); i++)
    for (unsigned int j = 0; j < scaledFeatureSet.Cols(); j++)
      if (std::isnan(scaledFeatureSet(i, j)))
        scaledFeatureSet(i, j) = 0;

  WriteCSVFilesWithHorizontalAndVerticalHeaders(scaledFeatureSet, patient_ids, StringFeatureLabels, outputdirectory + "/ScaledFeatures.csv");
  WriteCSVFiles(meanVector, outputdirectory + "/Survival_ZScore_Mean.csv",true);
  WriteCSVFiles(stdVector, outputdirectory + "/Survival_ZScore_Std.csv",true);

  typedef vnl_matrix<double> MatrixType;
  MatrixType data;
  VectorDouble SixModelLabels, EighteenModelLabels;
  mFeatureExtractionLocalPtr.FormulateSurvivalTrainingData(AllSurvival, SixModelLabels, EighteenModelLabels);

  //select 6-months and 18-months model features using routines of training module
  TrainingModule mTrainingModule; 
  VectorDouble cvaccuracies;
  VectorDouble EffectSize = mTrainingModule.EffectSizeFeatureSelection(scaledFeatureSet, SixModelLabels);
  for (unsigned int eSizeCounter = 0; eSizeCounter < EffectSize.size(); eSizeCounter++)
  {
    if (EffectSize[eSizeCounter] < 0)
      EffectSize[eSizeCounter] = EffectSize[eSizeCounter] * -1;
  }
  std::vector<size_t> indices = mTrainingModule.sort_indexes(EffectSize);
  VectorDouble selectedfeatures_6Months;
  for (int index = 0; index < indices.size()*0.2; index++)
    selectedfeatures_6Months.push_back(indices[index]);
  
  //VectorDouble selectedfeatures_6Months;
  //for (int index = 0; index < EffectSize.size(); index++)
  //  if(EffectSize[index]>=0.5)
  //    selectedfeatures_6Months.push_back(index);

  EffectSize = mTrainingModule.EffectSizeFeatureSelection(scaledFeatureSet, EighteenModelLabels);
  for (unsigned int eSizeCounter = 0; eSizeCounter < EffectSize.size(); eSizeCounter++)
  {
    if (EffectSize[eSizeCounter] < 0)
      EffectSize[eSizeCounter] = EffectSize[eSizeCounter] * -1;
  }
  VectorDouble selectedfeatures_18Months;
  indices = mTrainingModule.sort_indexes(EffectSize);
  for (int index = 0; index < indices.size()*0.2; index++)
    selectedfeatures_18Months.push_back(indices[index]);

  WriteCSVFiles(selectedfeatures_6Months, outputdirectory + "/Survival_SelectedFeatures_6Months.csv",true);
  WriteCSVFiles(selectedfeatures_18Months, outputdirectory + "/Survival_SelectedFeatures_18Months.csv",true);

  //extract 6-months and 18-months model features as specified by the training module
  VariableSizeMatrixType SixModelSelectedFeatures = SelectModelFeatures(scaledFeatureSet,selectedfeatures_6Months);
  VariableSizeMatrixType EighteenModelSelectedFeatures = SelectModelFeatures(scaledFeatureSet,selectedfeatures_18Months);

  //writing selected 6-months and 18-months model features
  std::vector<std::string> SelectedFeatureLabels_6months, SelectedFeatureLabels_18months;
  for (int index = 0; index < selectedfeatures_6Months.size(); index++)
  {
    int currentindex = selectedfeatures_6Months[index];
    SelectedFeatureLabels_6months.push_back(FeatureLabels[currentindex]);
  }
  for (int index = 0; index < selectedfeatures_18Months.size(); index++)
  {
    int currentindex = selectedfeatures_18Months[index];
    SelectedFeatureLabels_18months.push_back(FeatureLabels[currentindex]);
  }  
  WriteCSVFilesWithHorizontalAndVerticalHeaders(SixModelSelectedFeatures, patient_ids, SelectedFeatureLabels_6months, outputdirectory + "/SelectedFeatures_6months.csv");
  WriteCSVFilesWithHorizontalAndVerticalHeaders(EighteenModelSelectedFeatures, patient_ids, SelectedFeatureLabels_18months, outputdirectory + "/SelectedFeatures_18months.csv");

  //append labels to 6-months and 18-months model features as the training function expects labels in the last column
  VariableSizeMatrixType finaldatamatrix_6Months, finaldatamatrix_18Months;
  finaldatamatrix_6Months.SetSize(SixModelSelectedFeatures.Rows(), SixModelSelectedFeatures.Cols() + 1);
  finaldatamatrix_18Months.SetSize(EighteenModelSelectedFeatures.Rows(), EighteenModelSelectedFeatures.Cols() + 1);
  for (unsigned int i = 0; i < finaldatamatrix_6Months.Rows(); i++)
  {
    for (unsigned int j = 0; j < finaldatamatrix_6Months.Cols() - 1; j++)
      finaldatamatrix_6Months(i, j) = SixModelSelectedFeatures(i, j);
    finaldatamatrix_6Months(i, finaldatamatrix_6Months.Cols() - 1) = SixModelLabels[i];
  }
  for (unsigned int i = 0; i < finaldatamatrix_18Months.Rows(); i++)
  {
    for (unsigned int j = 0; j < finaldatamatrix_18Months.Cols() - 1; j++)
      finaldatamatrix_18Months(i, j) = EighteenModelSelectedFeatures(i, j);
    finaldatamatrix_18Months(i, finaldatamatrix_18Months.Cols() - 1) = EighteenModelLabels[i];
  }
  std::cout << std::endl << "Building model....." << std::endl;
  try
   {
	   trainOpenCVSVM(finaldatamatrix_6Months, outputdirectory + "/" + mSixTrainedFile, false, CAPTK::ApplicationCallingSVM::Survival);
	   trainOpenCVSVM(finaldatamatrix_18Months, outputdirectory + "/" + mEighteenTrainedFile, false, CAPTK::ApplicationCallingSVM::Survival);
   }
   catch (const std::exception& e1)
   {
     logger.WriteError("Training on the given subjects failed. Error code : " + std::string(e1.what()));
     return false;
   }
   std::cout << std::endl << "Model saved to the output directory." << std::endl;
   return true;
}

VariableLengthVectorType SurvivalPredictor::DistanceFunction(const VariableSizeMatrixType &testData, const std::string &filename)
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
  double rho;
  double bestg;
	SupportVectors.SetSize(dataMatrix.rows(), dataMatrix.cols() - 1);
	Coefficients.SetSize(dataMatrix.rows(), 1);
	Distances.SetSize(testData.Rows(), 1);

  //copy the support vectors, coefficients, rho, and bestg values from the model file
  for (unsigned int i = 0; i < dataMatrix.rows(); i++)
  {
    unsigned int j = 0;
    for (j = 0; j < dataMatrix.cols() - 2; j++)
      SupportVectors(i, j) = dataMatrix(i, j);
    Coefficients[i] = dataMatrix(i, j);
    if (i == 0)
      rho = dataMatrix(i, j + 1);
    if(i==1)
      bestg= dataMatrix(i, j + 1);
  }

  //calculate distance value for each given sample in the test data
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

VectorDouble SurvivalPredictor::CombineEstimates(const VariableLengthVectorType &estimates1, const VariableLengthVectorType &estimates2)
{
	VectorDouble returnVec;
	returnVec.resize(estimates1.Size());
  // +1/-1 are the class labels, and 2 is the capping limit that we apply on the predicted distances.
	for (size_t i = 0; i < estimates1.Size(); i++)
	{
		float temp_abs, temp_pos1, temp_neg1, temp_1, temp_2;
		// estimate for 1st vector
		if (std::abs(estimates1[i]) < 2)
			temp_abs = estimates1[i];
		else
			temp_abs = 0;

		if (estimates1[i] > 1)
			temp_pos1 = 1;
		else
			temp_pos1 = 0;

		if (estimates1[i] < -1)
			temp_neg1 = 1;
		else
			temp_neg1 = 0;
		temp_1 = temp_abs + (temp_pos1 - temp_neg1);

		// estimate for 2nd vector, all temp values are getting overwritten
		if (std::abs(estimates2[i]) < 2)
			temp_abs = estimates2[i];
		else
			temp_abs = 0;

		if (estimates2[i] > 1)
			temp_pos1 = 1;
		else
			temp_pos1 = 0;

		if (estimates2[i] < -1)
			temp_neg1 = 1;
		else
			temp_neg1 = 0;

    temp_2 = temp_abs + (temp_pos1 - temp_neg1);

		// combine the two
		returnVec[i] = temp_1 + temp_2;
	}
	return returnVec;
}

VectorDouble SurvivalPredictor::CombineEstimates(const VectorDouble &estimates1, const VectorDouble &estimates2)
{
	VectorDouble returnVec;
	returnVec.resize(estimates1.size());
  // +1/-1 are the class labels, and 2 is the capping limit that we apply on the predicted distances.
	for (size_t i = 0; i < estimates1.size(); i++)
	{
		float temp_abs, temp_pos1, temp_neg1, temp_1, temp_2;
		// estimate for 1st vector
		if (std::abs(estimates1[i]) < 2)
			temp_abs = estimates1[i];
    else
      temp_abs = 0;

		if (estimates1[i] > 1)
			temp_pos1 = 1;
    else
      temp_pos1 = 0;

		if (estimates1[i] < -1)
			temp_neg1 = 1;
		else
			temp_neg1 = 0;
		temp_1 = temp_abs + (temp_pos1 - temp_neg1);

		// estimate for 2nd vector, all temp values are getting overwritten
		if (std::abs(estimates2[i]) < 2)
			temp_abs = estimates2[i];
		else
			temp_abs = 0;

		if (estimates2[i] > 1)
			temp_pos1 = 1;
    else
      temp_pos1 = 0;

		if (estimates2[i] < -1)
			temp_neg1 = 1;
    else
      temp_neg1 = 0;
		temp_2 = temp_abs + (temp_pos1 - temp_neg1);

		// combine the two
		returnVec[i] = temp_1 + temp_2;
	}
	return returnVec;
}

VariableLengthVectorType SurvivalPredictor::DistanceFunctionLinear(const VariableSizeMatrixType &testData, const std::string &filename)
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
  double Rho = 0;
  SupportVectors.SetSize(dataMatrix.rows(), dataMatrix.cols() - 2);
  Coefficients.SetSize(dataMatrix.rows(), 1);
  Distances.SetSize(testData.Rows(), 1);
  
  //copy the support vectors, coefficients, and rho values from the model file
  for (unsigned int i = 0; i < dataMatrix.rows(); i++)
  {
    unsigned int j = 0;
    for (j = 0; j < dataMatrix.cols() - 2; j++)
      SupportVectors(i, j) = dataMatrix(i, j);
    Coefficients[i] = dataMatrix(i, j);
    if (i == 0)
      Rho = dataMatrix(i, j + 1);
  }
  VariableSizeMatrixType TransposedSupportVectors = MatrixTranspose(SupportVectors);
  VariableSizeMatrixType w;
  w.SetSize(TransposedSupportVectors.Rows(), 1);
  for (unsigned int svID = 0; svID < TransposedSupportVectors.Rows(); svID++)
  {
    double currentSum = 0;
    for (unsigned int iterator = 0; iterator < TransposedSupportVectors.Cols(); iterator++)
      currentSum = currentSum + TransposedSupportVectors(svID, iterator)*Coefficients[iterator];
    w(svID, 0) = currentSum;
  }
  VariableSizeMatrixType wTranspose = MatrixTranspose(w);  //1x7   1x7

  //calculate distance value for each given sample in the test data
  for (unsigned int patID = 0; patID < testData.Rows(); patID++)
  {
    double distance = 0;
    for (unsigned int svID = 0; svID < wTranspose.Cols(); svID++)
      distance = distance + wTranspose(0, svID)*testData(patID, svID);
    Distances[patID] = distance - Rho;
  }
  return Distances;
}

VariableSizeMatrixType SurvivalPredictor::MatrixTranspose(const VariableSizeMatrixType &inputmatrix)
{
  //calculate transpose of a matrix
  VariableSizeMatrixType output;
  output.SetSize(inputmatrix.Cols(), inputmatrix.Rows());

  for (unsigned int i = 0; i < output.Rows(); i++)
    for (unsigned int j = 0; j < output.Cols(); j++)
      output(i, j) = inputmatrix(j, i);
  return output;
}

VectorDouble SurvivalPredictor::SurvivalPredictionOnExistingModel(const std::string &modeldirectory, const std::string &inputdirectory, const std::vector < std::map < CAPTK::ImageModalityType, std::string>> &qualifiedsubjects, const std::string &outputdirectory)
{
	typedef itk::CSVArray2DFileReader<double> ReaderType;
	VariableSizeMatrixType HistogramFeaturesConfigurations;
	HistogramFeaturesConfigurations.SetSize(33, 3); //11 modalities*3 regions = 33 configurations*3 histogram features for each configuration
	VectorDouble results;

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
  VariableLengthVectorType selectedfeatures_6months;
  VariableLengthVectorType selectedfeatures_18months;
  
  try
	{
		reader->SetFileName(modeldirectory + "/Survival_ZScore_Mean.csv");
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
		logger.WriteError("Error in reading the file: " + modeldirectory + "/Survival_ZScore_Mean.csv. Error code : " + std::string(e1.what()));
		return results;
	}
  std::cout << "done" << std::endl;
  MatrixType stdMatrix;
	try
	{
		reader->SetFileName(modeldirectory + "/Survival_ZScore_Std.csv");
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
		logger.WriteError("Error in reading the file: " + modeldirectory + "/Survival_ZScore_Std.csv. Error code : " + std::string(e1.what()));
		return results;
	}
  //read selected features of 6-months model from .csv file
  MatrixType features6Matrix;
  try
  {
    reader->SetFileName(modeldirectory + "/Survival_SelectedFeatures_6Months.csv");
    reader->SetFieldDelimiterCharacter(',');
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();
    reader->Parse();
    features6Matrix = reader->GetArray2DDataObject()->GetMatrix();
    selectedfeatures_6months.SetSize(features6Matrix.size());
    for (unsigned int i = 0; i < features6Matrix.size(); i++)
      selectedfeatures_6months[i] = features6Matrix(i, 0);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Error in reading the file: " + modeldirectory + "/Survival_SelectedFeatures_6Months.csv. Error code : " + std::string(e1.what()));
    return results;
  }
  //read selected features of 18-months model from .csv file
  MatrixType features18Matrix;
  try
  {
    reader->SetFileName(modeldirectory + "/Survival_SelectedFeatures_18Months.csv");
    reader->SetFieldDelimiterCharacter(',');
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();
    reader->Parse();
    features18Matrix = reader->GetArray2DDataObject()->GetMatrix();
    selectedfeatures_18months.SetSize(features18Matrix.size());
    for (unsigned int i = 0; i < features18Matrix.size(); i++)
      selectedfeatures_18months[i] = features18Matrix(i, 0);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Error in reading the file: " + modeldirectory + "/Survival_SelectedFeatures_18Months.csv. Error code : " + std::string(e1.what()));
    return results;
  }

  VariableSizeMatrixType FeaturesOfAllSubjects;
	FeaturesOfAllSubjects.SetSize(qualifiedsubjects.size(), SURVIVAL_NO_OF_FEATURES);
  std::vector<std::string> patient_ids;
	for (unsigned int sid = 0; sid < qualifiedsubjects.size(); sid++)
	{
		std::cout << "Subject No:" << sid << std::endl;
		std::map<CAPTK::ImageModalityType, std::string> currentsubject = qualifiedsubjects[sid];
    patient_ids.push_back(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]));
		try
		{
			ImageType::Pointer LabelImagePointer = cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG]));
			ImageType::Pointer AtlasImagePointer = cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_ATLAS]));
      ImageType::Pointer TemplateImagePointer = cbica::ReadImage<ImageType>(getCaPTkDataDir() + "/survival/Template.nii.gz");
			ImageType::Pointer RCBVImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_RCBV])));
			ImageType::Pointer PHImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PH])));
			ImageType::Pointer T1CEImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1CE])));
			ImageType::Pointer T2FlairImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR])));
			ImageType::Pointer T1ImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1])));
			ImageType::Pointer T2ImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2])));
			ImageType::Pointer AXImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_AX])));
			ImageType::Pointer RADImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_RAD])));
			ImageType::Pointer FAImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_FA])));
			ImageType::Pointer TRImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_TR])));
			ImageType::Pointer PSRImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PSR])));

			VectorDouble TestFeatures = CalculateFeatures<ImageType>(T1CEImagePointer, T2FlairImagePointer, T1ImagePointer, T2ImagePointer,
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
  
  std::string FeatureLabels[SURVIVAL_NO_OF_FEATURES] = { "Age","ET","NCR","ED","brain size","ET+NCR","(ET+NCR)x100/brain size","EDx100/brain size","ETx100/(ET+NCR)","NCRx100/(ET+NCR)","NCRx100/(ET+NCR)",
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
  for (int index = 0; index < SURVIVAL_NO_OF_FEATURES; index++)
    StringFeatureLabels.push_back(FeatureLabels[index]);

  WriteCSVFilesWithHorizontalAndVerticalHeaders(FeaturesOfAllSubjects, patient_ids, StringFeatureLabels, outputdirectory + "/RawFeatures.csv");

	VariableSizeMatrixType ScaledTestingData = mFeatureScalingLocalPtr.ScaleGivenTestingFeatures(FeaturesOfAllSubjects, mean, stddevition);
  for (unsigned int i = 0; i < ScaledTestingData.Rows(); i++)
    for (unsigned int j = 0; j < ScaledTestingData.Cols(); j++)
      if (std::isnan(ScaledTestingData(i, j)))
        ScaledTestingData(i, j) = 0;

  //write scaled features in a .csv file
  WriteCSVFilesWithHorizontalAndVerticalHeaders(ScaledTestingData, patient_ids, StringFeatureLabels, outputdirectory + "/ScaledFeatures.csv");
  
  VariableSizeMatrixType SixModelSelectedFeatures = SelectModelFeatures(ScaledTestingData,selectedfeatures_6months);
	VariableSizeMatrixType EighteenModelSelectedFeatures = SelectModelFeatures(ScaledTestingData,selectedfeatures_18months);
  std::vector<std::string> SelectedFeatureLabels_6months, SelectedFeatureLabels_18months;
  for (int index = 0; index < selectedfeatures_6months.Size(); index++)
  {
    int currentindex = selectedfeatures_6months[index];
    SelectedFeatureLabels_6months.push_back(FeatureLabels[currentindex]);
  }
  for (int index = 0; index < selectedfeatures_18months.Size(); index++)
  {
    int currentindex = selectedfeatures_18months[index];
    SelectedFeatureLabels_18months.push_back(FeatureLabels[currentindex]);
  }

  //write selected features in a .csv file
  WriteCSVFilesWithHorizontalAndVerticalHeaders(SixModelSelectedFeatures, patient_ids, SelectedFeatureLabels_6months, outputdirectory + "/SelectedFeatures_6Months.csv");
  WriteCSVFilesWithHorizontalAndVerticalHeaders(EighteenModelSelectedFeatures, patient_ids, SelectedFeatureLabels_18months, outputdirectory + "/SelectedFeatures_18Months.csv");

  try
	{
		std::ofstream myfile;
		myfile.open(outputdirectory + "/results.csv");
		myfile << "SubjectName,SPI (6 months), SPI (18 months), Composite SPI\n";
		if (cbica::fileExists(modeldirectory + "/Survival_SVM_Model6.csv") == true && cbica::fileExists(modeldirectory + "/Survival_SVM_Model18.csv") == true)
		{
			VariableLengthVectorType result_6;
			VariableLengthVectorType result_18;
			result_6 = DistanceFunctionLinear(SixModelSelectedFeatures, modeldirectory + "/Survival_SVM_Model6.csv");
			result_18 = DistanceFunctionLinear(EighteenModelSelectedFeatures, modeldirectory + "/Survival_SVM_Model18.csv");
			results = CombineEstimates(result_6, result_18);
			for (size_t i = 0; i < results.size(); i++)
			{
				std::map<CAPTK::ImageModalityType, std::string> currentsubject = qualifiedsubjects[i];
				myfile << static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]) + "," + std::to_string(result_6[i]) + "," + std::to_string(result_18[i]) + "," + std::to_string(results[i]) + "\n";
			}
		}
		else if (cbica::fileExists(modeldirectory + "/Survival_SVM_Model6.xml") == true && cbica::fileExists(modeldirectory + "/Survival_SVM_Model18.xml") == true)
		{
			VectorDouble result_6;
			VectorDouble result_18;
			result_6 = testOpenCVSVM(SixModelSelectedFeatures, modeldirectory + "/Survival_SVM_Model6.xml");
			result_18 = testOpenCVSVM(EighteenModelSelectedFeatures, modeldirectory + "/Survival_SVM_Model18.xml");
			results = CombineEstimates(result_6, result_18);
			for (size_t i = 0; i < results.size(); i++)
			{
				std::map<CAPTK::ImageModalityType, std::string> currentsubject = qualifiedsubjects[i];
				myfile << static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]) + "," + std::to_string(result_6[i]) + "," + std::to_string(result_18[i]) + "," + std::to_string(results[i]) + "\n";
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

VariableSizeMatrixType SurvivalPredictor::SelectModelFeatures(const VariableSizeMatrixType &ModelFeatures, const VariableLengthVectorType &selectedFeatures)
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

VariableSizeMatrixType SurvivalPredictor::SelectModelFeatures(const VariableSizeMatrixType &ModelFeatures, const VectorDouble &selectedFeatures)
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
