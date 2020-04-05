#include "SurvivalPredictor.h"
#include "itkCSVNumericObjectFileWriter.h"
#include "CaPTkGUIUtils.h"
#include "CaPTkClassifierUtils.h"


typedef itk::Image< float, 3 > ImageType;
//SurvivalPredictor::~SurvivalPredictor()
//{
//}

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

int SurvivalPredictor::PrepareNewSurvivalPredictionModel(const std::string &inputdirectory,const std::vector< std::map< CAPTK::ImageModalityType, std::string > > &qualifiedsubjects, const std::string &outputdirectory)
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
 //---------------------------------------------------------------------------
  FeaturesOfAllSubjects.SetSize(qualifiedsubjects.size(), 161);
  for (unsigned int sid = 0; sid < qualifiedsubjects.size(); sid++)
  {
	  std::cout << "Patient's data loaded:" << sid + 1 << std::endl;
	  std::map< CAPTK::ImageModalityType, std::string > currentsubject = qualifiedsubjects[sid];
	  try
	  {
		  ImageType::Pointer LabelImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG]));
		  ImageType::Pointer AtlasImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_ATLAS]));
		  //ImageType::Pointer TemplateImagePointer = ReadNiftiImage<ImageType>("../data/survival/Template.nii.gz");
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
    //
//  for (int i = 0; i < 105; i++)
//    for (int j = 0; j < 166; j++)
//      data(i, j) = scaledFeatureSet(i, j);
//  writer->SetFileName("scaled_train_features.csv");
//  writer->SetInput(&data);
//  writer->Write();
//
//
//
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
  typedef vnl_matrix<double> MatrixType;
  MatrixType data;
  

  VariableSizeMatrixType SixModelFeatures;
  VariableSizeMatrixType EighteenModelFeatures;
  mFeatureExtractionLocalPtr.FormulateSurvivalTrainingData(scaledFeatureSet, AllSurvival, SixModelFeatures, EighteenModelFeatures);
  try
  {
	  data.set_size(161, 1); // TOCHECK - are these hard coded sizes fine?
	  for (unsigned int i = 0; i < meanVector.Size(); i++)
		  data(i, 0) = meanVector[i];
	  typedef itk::CSVNumericObjectFileWriter<double, 161, 1> WriterTypeVector;
	  WriterTypeVector::Pointer writerv = WriterTypeVector::New();
	  writerv->SetFileName(outputdirectory + "/Survival_ZScore_Mean.csv");
	  writerv->SetInput(&data);
	  writerv->Write();

	  for (unsigned int i = 0; i < stdVector.Size(); i++)
		  data(i, 0) = stdVector[i];
	  writerv->SetFileName(outputdirectory + "/Survival_ZScore_Std.csv");
	  writerv->SetInput(&data);
	  writerv->Write();
  }
  catch (const std::exception& e1)
  {
	  logger.WriteError("Error in writing output files to the output directory = " + outputdirectory + "Error code : " + std::string(e1.what()));
	  return false;
  }
  //---------------------------------------------------------------------------
  VariableLengthVectorType selectedfeatures_6Months;
  VariableLengthVectorType selectedfeatures_18Months;
  VariableSizeMatrixType SixModelSelectedFeatures = SelectModelFeatures(SixModelFeatures,selectedfeatures_6Months);
  VariableSizeMatrixType EighteenModelSelectedFeatures = SelectModelFeatures(EighteenModelFeatures,selectedfeatures_18Months);

  //WriteCSVFiles(FeaturesOfAllSubjects, outputdirectory + "/FeaturesOfAllSubjects.csv");
  //WriteCSVFiles(scaledFeatureSet, outputdirectory + "/scaledFeatureSet.csv");
  //WriteCSVFiles(SixModelFeatures, outputdirectory + "/SixModelFeatures.csv");
  //WriteCSVFiles(EighteenModelFeatures, outputdirectory + "/EighteenModelFeatures.csv");
  //WriteCSVFiles(SixModelSelectedFeatures, outputdirectory + "/SixModelSelectedFeatures.csv");
  //WriteCSVFiles(EighteenModelSelectedFeatures, outputdirectory + "/EighteenModelSelectedFeatures.csv");

   try
   {
	   trainOpenCVSVM(SixModelSelectedFeatures, outputdirectory + "/" + mSixTrainedFile, false, CAPTK::ApplicationCallingSVM::Survival);
	   trainOpenCVSVM(EighteenModelSelectedFeatures, outputdirectory + "/" + mEighteenTrainedFile, false, CAPTK::ApplicationCallingSVM::Survival);
   }
   catch (const std::exception& e1)
   {
     logger.WriteError("Training on the given subjects failed. Error code : " + std::string(e1.what()));
     return false;
   }
   std::cout << std::endl << "Model saved to the output directory." << std::endl;
   return true;
}
void SurvivalPredictor::WriteCSVFiles(VariableSizeMatrixType inputdata, std::string filepath)
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


VectorDouble SurvivalPredictor::CombineEstimates(const VectorDouble &estimates1, const VectorDouble &estimates2)
{
	VectorDouble returnVec;
	returnVec.resize(estimates1.size());
	for (size_t i = 0; i < estimates1.size(); i++)
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
  
  //copy the support vectors, coefficients, rho, and bestg values from the model file
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
  //----------------------------------------------------
	VariableSizeMatrixType FeaturesOfAllSubjects;
	FeaturesOfAllSubjects.SetSize(qualifiedsubjects.size(), 164);

	for (unsigned int sid = 0; sid < qualifiedsubjects.size(); sid++)
	{
		std::cout << "Subject No:" << sid << std::endl;
		std::map<CAPTK::ImageModalityType, std::string> currentsubject = qualifiedsubjects[sid];
		try
		{
			ImageType::Pointer LabelImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG]));
			ImageType::Pointer AtlasImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_ATLAS]));
			//ImageType::Pointer TemplateImagePointer = ReadNiftiImage<ImageType>("../data/survival/Template.nii.gz");
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
	VariableSizeMatrixType ScaledTestingData = mFeatureScalingLocalPtr.ScaleGivenTestingFeatures(FeaturesOfAllSubjects, mean, stddevition);
  for (unsigned int i = 0; i < ScaledTestingData.Rows(); i++)
    for (unsigned int j = 0; j < ScaledTestingData.Cols(); j++)
      if (std::isnan(ScaledTestingData(i, j)))
        ScaledTestingData(i, j) = 0;

	VariableSizeMatrixType ScaledFeatureSetAfterAddingLabel;
	ScaledFeatureSetAfterAddingLabel.SetSize(ScaledTestingData.Rows(), ScaledTestingData.Cols() + 1);
	for (unsigned int i = 0; i < ScaledTestingData.Rows(); i++)
	{
		unsigned int j = 0;
		for (j = 0; j < ScaledTestingData.Cols(); j++)
			ScaledFeatureSetAfterAddingLabel(i, j) = ScaledTestingData(i, j);
		ScaledFeatureSetAfterAddingLabel(i, j) = 0;
	}


  VariableSizeMatrixType SixModelSelectedFeatures = SelectModelFeatures(ScaledFeatureSetAfterAddingLabel,selectedfeatures_6months);
	VariableSizeMatrixType EighteenModelSelectedFeatures = SelectModelFeatures(ScaledFeatureSetAfterAddingLabel,selectedfeatures_18months);

  WriteCSVFiles(FeaturesOfAllSubjects, outputdirectory + "/raw_features.csv");
  WriteCSVFiles(ScaledTestingData, outputdirectory + "/scaled_features.csv");
  WriteCSVFiles(ScaledFeatureSetAfterAddingLabel, outputdirectory + "/scaled_features_with_label.csv");
  WriteCSVFiles(SixModelSelectedFeatures, outputdirectory + "/selectedfeatures_6months.csv");
  WriteCSVFiles(EighteenModelSelectedFeatures, outputdirectory + "/selectedfeatures_18months.csv");
  
	try
	{
		std::ofstream myfile;
		myfile.open(outputdirectory + "/results.csv");
		myfile << "SubjectName,SPI (6 months), SPI (18 months), Composite SPI\n";
		if (cbica::fileExists(modeldirectory + "/Survival_SVM_Model6.csv") == true && cbica::fileExists(modeldirectory + "/Survival_SVM_Model18.csv") == true)
		{
			VariableLengthVectorType result_6;
			VariableLengthVectorType result_18;
			result_6 = DistanceFunction(SixModelSelectedFeatures, modeldirectory + "/Survival_SVM_Model6.csv");
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
			result_6 = testOpenCVSVM(ScaledTestingData, modeldirectory + "/Survival_SVM_Model6.xml");
			result_18 = testOpenCVSVM(ScaledTestingData, modeldirectory + "/Survival_SVM_Model18.xml");
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

VariableSizeMatrixType SurvivalPredictor::SelectModelFeatures(const VariableSizeMatrixType &SixModelFeatures,const VariableLengthVectorType selectedfeatures)
{
  VariableSizeMatrixType SixModelSelectedFeatures; 
  SixModelSelectedFeatures.SetSize(SixModelFeatures.Rows(),selectedfeatures.Size()+1);
  int counter = 0;
  for (unsigned int i = 0; i < selectedfeatures.Size(); i++)
  {
    for (unsigned int j = 0; j < SixModelFeatures.Rows(); j++)
      SixModelSelectedFeatures(j, counter) = SixModelFeatures(j, selectedfeatures[i]);
    counter++;
  }
 for (unsigned int j = 0; j < SixModelFeatures.Rows(); j++)
      SixModelSelectedFeatures(j, selectedfeatures.Size()) = SixModelFeatures(j, SixModelFeatures.Cols()-1);

  return SixModelSelectedFeatures;
}
