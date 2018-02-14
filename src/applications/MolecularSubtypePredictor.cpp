#include "MolecularSubtypePredictor.h"


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

int MolecularSubtypePredictor::PrepareNewMolecularPredictionModel(const std::string &inputdirectory, const std::vector< std::map< ImageModalityType, std::string > > &qualifiedsubjects, const std::string &outputdirectory)
{
	VectorDouble AllSurvival;
	VariableSizeMatrixType FeaturesOfAllSubjects;
	VariableSizeMatrixType HistogramFeaturesConfigurations;
	HistogramFeaturesConfigurations.SetSize(33, 3); //11 modalities*3 regions = 33 configurations*3 histogram features for each configuration
	CSVFileReaderType::Pointer reader = CSVFileReaderType::New();
	MatrixType dataMatrix;
	try
	{
		reader->SetFileName("../data/molecular/Survival_HMFeatures_Configuration.csv");
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
		mLastErrorMessage = "Cannot find the file: ../data/molecular/Survival_HMFeatures_Configuration.csv. Error code : " + std::string(e1.what());
		std::cout << mLastErrorMessage << std::endl;
		logger.WriteError(mLastErrorMessage);
		return EXIT_FAILURE;
	}
	//---------------------------------------------------------------------------
	FeaturesOfAllSubjects.SetSize(qualifiedsubjects.size(), 161);
	for (unsigned int sid = 0; sid < qualifiedsubjects.size(); sid++)
	{
		std::cout << "Patient's data loaded:" << sid + 1 << std::endl;
		std::map< ImageModalityType, std::string > currentsubject = qualifiedsubjects[sid];
		try
		{
			ImageType::Pointer LabelImagePointer	= ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_SEG]));
			ImageType::Pointer AtlasImagePointer	= ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_ATLAS]));
			ImageType::Pointer TemplateImagePointer = ReadNiftiImage<ImageType>("../data/survival/Template.nii.gz");
			ImageType::Pointer RCBVImagePointer		= RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_RCBV])));
			ImageType::Pointer PHImagePointer		= RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_PH])));
			ImageType::Pointer T1CEImagePointer		= RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_T1CE])));
			ImageType::Pointer T2FlairImagePointer	= RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_T2FLAIR])));
			ImageType::Pointer T1ImagePointer		= RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_T1])));
			ImageType::Pointer T2ImagePointer		= RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_T2])));
			ImageType::Pointer AXImagePointer		= RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_AX])));
			ImageType::Pointer RADImagePointer		= RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_RAD])));
			ImageType::Pointer FAImagePointer		= RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_FA])));
			ImageType::Pointer TRImagePointer		= RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_TR])));
			ImageType::Pointer PSRImagePointer		= RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_PSR])));

			VectorDouble ages;
			VectorDouble survival;

			reader->SetFileName(static_cast<std::string>(currentsubject[IMAGE_TYPE_FEATURES]));
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
			std::string name = static_cast<std::string>(currentsubject[IMAGE_TYPE_SUDOID]);
			logger.WriteError("The data for the subject " + name + "was not loaded: " + e1.what());
		}
	}
	std::cout << std::endl << "Building model....." << std::endl;
	VariableSizeMatrixType scaledFeatureSet;
	scaledFeatureSet.SetSize(qualifiedsubjects.size(), 161);
	VariableLengthVectorType meanVector;
	VariableLengthVectorType stdVector;
	mFeatureScalingLocalPtr.ScaleGivenTrainingFeatures(FeaturesOfAllSubjects, scaledFeatureSet, meanVector, stdVector);

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
	
	  data.set_size(166, 1); // TOCHECK - are these hard coded sizes fine?
	  for (unsigned int i = 0; i < meanVector.Size(); i++)
	    data(i, 0) = meanVector[i];
	  typedef itk::CSVNumericObjectFileWriter<double, 166, 1> WriterTypeVector;
	  WriterTypeVector::Pointer writerv = WriterTypeVector::New();
	  writerv->SetFileName(outputdirectory + "/Molecular_ZScore_Mean.csv");
	  writerv->SetInput(&data);
	  writerv->Write();
	  
	  for (unsigned int i = 0; i < stdVector.Size(); i++)
	    data(i, 0) = stdVector[i];
	  writerv->SetFileName(outputdirectory + "/Molecular_ZScore_Std.csv");
	  writerv->SetInput(&data);
	  writerv->Write();
	
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
		trainOpenCVSVM(neuralModelFeatures, outputdirectory + "/" + mNeuralTrainedFile, false, Survival);
		trainOpenCVSVM(proneuralModelFeatures, outputdirectory + "/" + mProneuralTrainedFile, false, Survival);
		trainOpenCVSVM(classicalModelFeatures, outputdirectory + "/" + mClassicalTrainedFile, false, Survival);
		trainOpenCVSVM(messenchymalModelFeatures, outputdirectory + "/" + mMessenchymalTrainedFile, false, Survival);
	}
	catch (const std::exception& e1)
	{
		mLastErrorMessage = "Training on the given subjects failed. Error code : " + std::string(e1.what());
		logger.WriteError(mLastErrorMessage);
		std::cout << mLastErrorMessage << std::endl;
		return EXIT_FAILURE;
	}
	std::cout << std::endl << "Model saved to the output directory." << std::endl;
	int c;
	std::cin >> c;
	return EXIT_SUCCESS;
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
			returnVec.push_back(MOLECULAR_SUBTYPES::PRONEURAL);
		else if (result_neural[i] > result_proneural[i] && result_neural[i] > result_messenchymal[i] && result_neural[i] > result_classical[i])
			returnVec.push_back(MOLECULAR_SUBTYPES::NEURAL);
		else if (result_messenchymal[i] > result_proneural[i] && result_messenchymal[i] > result_neural[i] && result_messenchymal[i] > result_classical[i])
			returnVec.push_back(MOLECULAR_SUBTYPES::MESSENCHYMAL);
		else
			returnVec.push_back(MOLECULAR_SUBTYPES::CLASSICAL);
	}
	return returnVec;
}


VectorDouble MolecularSubtypePredictor::MolecularSubtypePredictionOnExistingModel(const std::string &modeldirectory, const std::string &inputdirectory, const std::vector < std::map < ImageModalityType, std::string>> &qualifiedsubjects, const std::string &outputdirectory)
{
	typedef itk::CSVArray2DFileReader<double> ReaderType;
	VariableSizeMatrixType HistogramFeaturesConfigurations;
	HistogramFeaturesConfigurations.SetSize(33, 3); //11 modalities*3 regions = 33 configurations*3 histogram features for each configuration

	CSVFileReaderType::Pointer reader = CSVFileReaderType::New();
	VectorDouble ages;
	MatrixType dataMatrix;
	try
	{
		reader->SetFileName("../data/molecular/Survival_HMFeatures_Configuration.csv");
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
		logger.WriteError("Cannot find the file 'Survival_HMFeatures_Configuration.csv' in the input directory. Error code : " + std::string(e1.what()));
		exit(EXIT_FAILURE);
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
		logger.WriteError("Cannot find the file 'mean.csv' in the model directory. Error code : " + std::string(e1.what()));
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
		logger.WriteError("Cannot find the file 'std.csv' in the model directory. Error code : " + std::string(e1.what()));
	}
	//----------------------------------------------------
	VariableSizeMatrixType FeaturesOfAllSubjects;
	FeaturesOfAllSubjects.SetSize(qualifiedsubjects.size(), 161);

	for (unsigned int sid = 0; sid < qualifiedsubjects.size(); sid++)
	{
		cout << "Subject No:" << sid << endl;
		std::map<ImageModalityType, std::string> currentsubject = qualifiedsubjects[sid];
		ImageType::Pointer LabelImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_SEG]));
		ImageType::Pointer AtlasImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_ATLAS]));
		ImageType::Pointer TemplateImagePointer = ReadNiftiImage<ImageType>("../data/survival/Template.nii.gz");
		ImageType::Pointer RCBVImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_RCBV])));
		ImageType::Pointer PHImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_PH])));
		ImageType::Pointer T1CEImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_T1CE])));
		ImageType::Pointer T2FlairImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_T2FLAIR])));
		ImageType::Pointer T1ImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_T1])));
		ImageType::Pointer T2ImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_T2])));
		ImageType::Pointer AXImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_AX])));
		ImageType::Pointer RADImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_RAD])));
		ImageType::Pointer FAImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_FA])));
		ImageType::Pointer TRImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_TR])));
		ImageType::Pointer PSRImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[IMAGE_TYPE_PSR])));

		VectorDouble TestFeatures = LoadTestData<ImageType>(T1CEImagePointer, T2FlairImagePointer, T1ImagePointer, T2ImagePointer,
			RCBVImagePointer, PSRImagePointer, PHImagePointer, AXImagePointer, FAImagePointer, RADImagePointer, TRImagePointer, LabelImagePointer, AtlasImagePointer, TemplateImagePointer, HistogramFeaturesConfigurations);

		double age;
		try
		{
			reader->SetFileName(static_cast<std::string>(currentsubject[IMAGE_TYPE_FEATURES]));
			reader->SetFieldDelimiterCharacter(',');
			reader->HasColumnHeadersOff();
			reader->HasRowHeadersOff();
			reader->Parse();
			dataMatrix = reader->GetArray2DDataObject()->GetMatrix();

			for (unsigned int i = 0; i < dataMatrix.rows(); i++)
				age = dataMatrix(i, 0);
		}
		catch (const std::exception& e1)
		{
			logger.WriteError("Cannot find the file 'features.csv' in the input directory. Error code : " + std::string(e1.what()));
		}
		FeaturesOfAllSubjects(sid, 0) = age;
		for (unsigned int i = 1; i <= TestFeatures.size(); i++)
			FeaturesOfAllSubjects(sid, i) = TestFeatures[i - 1];
	}
	//typedef vnl_matrix<double> MatrixType;
	//MatrixType data;
	//  data.set_size(105, 169);
	//  for (unsigned int i = 0; i < FeaturesOfAllSubjects.Rows(); i++)
	//  {
	//	  for (unsigned int j = 0; j < FeaturesOfAllSubjects.Cols(); j++)
	//    {
	//		data(i, j) = FeaturesOfAllSubjects(i, j);
	//    }
	//  }
	//  typedef itk::CSVNumericObjectFileWriter<double, 105, 169> WriterTypeMatrix;
	//  WriterTypeMatrix::Pointer writermatrix = WriterTypeMatrix::New();
	//  writermatrix->SetFileName(outputdirectory + "/plainfeatures.csv");
	//  writermatrix->SetInput(&data);
	//  writermatrix->Write();  
	  


	VariableSizeMatrixType ScaledTestingData = mFeatureScalingLocalPtr.ScaleGivenTestingFeatures(FeaturesOfAllSubjects, mean, stddevition);
	VariableSizeMatrixType ScaledFeatureSetAfterAddingLabel;
	ScaledFeatureSetAfterAddingLabel.SetSize(ScaledTestingData.Rows(), ScaledTestingData.Cols() + 1);
	for (unsigned int i = 0; i < ScaledTestingData.Rows(); i++)
	{
		unsigned int j = 0;
		for (j = 0; j < ScaledTestingData.Cols(); j++)
			ScaledFeatureSetAfterAddingLabel(i, j) = ScaledTestingData(i, j);
		ScaledFeatureSetAfterAddingLabel(i, j) = 0;
	}

	//for (unsigned int i = 0; i < ScaledFeatureSetAfterAddingLabel.Rows(); i++)
	//{
	//	for (unsigned int j = 0; j < ScaledFeatureSetAfterAddingLabel.Cols(); j++)
	//	{
	//		data(i, j) = ScaledFeatureSetAfterAddingLabel(i, j);
	//	}
	//}
	//typedef itk::CSVNumericObjectFileWriter<double, 105, 169> WriterTypeMatrix;
	//writermatrix->SetFileName(outputdirectory + "/scaledfeatures.csv");
	//writermatrix->SetInput(&data);
	//writermatrix->Write();





	//VariableSizeMatrixType SixModelSelectedFeatures = SelectSixMonthsModelFeatures(ScaledFeatureSetAfterAddingLabel);
	//VariableSizeMatrixType EighteenModelSelectedFeatures = SelectEighteenMonthsModelFeatures(ScaledFeatureSetAfterAddingLabel);

	VectorDouble results;
	try
	{
		std::ofstream myfile;
		myfile.open(outputdirectory + "/results.csv");
		myfile << "SubjectName,Result\n";

		if (cbica::fileExists(modeldirectory + "/" + mProneuralTrainedFile) == true && cbica::fileExists(modeldirectory + "/" + mNeuralTrainedFile) == true && cbica::fileExists(modeldirectory + "/" + mMessenchymalTrainedFile) == true && cbica::fileExists(modeldirectory + "/" + mClassicalTrainedFile) == true)
		{
			VectorDouble result_proneural;
			VectorDouble result_neural;
			VectorDouble result_classical;
			VectorDouble result_messenchymal;
			result_proneural = testOpenCVSVM(ScaledFeatureSetAfterAddingLabel, modeldirectory + "/" + mProneuralTrainedFile);
			result_neural = testOpenCVSVM(ScaledFeatureSetAfterAddingLabel, modeldirectory + "/" + mNeuralTrainedFile);
			result_messenchymal = testOpenCVSVM(ScaledFeatureSetAfterAddingLabel, modeldirectory + "/" + mMessenchymalTrainedFile);
			result_classical = testOpenCVSVM(ScaledFeatureSetAfterAddingLabel, modeldirectory + "/" + mClassicalTrainedFile);

			results = CombineEstimates(result_proneural, result_neural, result_messenchymal, result_classical);
			for (size_t i = 0; i < results.size(); i++)
			{
				std::map<ImageModalityType, std::string> currentsubject = qualifiedsubjects[i];
				myfile << static_cast<std::string>(currentsubject[IMAGE_TYPE_SUDOID]) + "," + std::to_string(result_proneural[i]) + "," + std::to_string(result_neural[i]) + "," + std::to_string(result_messenchymal[i]) + "," + std::to_string(result_classical[i]) + +"," + std::to_string(results[i]) + "\n";
			}
		}
		else
		{
			logger.WriteError("Error caught during testing: There is no exisitg model file in the model directory: " + modeldirectory);
			exit(EXIT_FAILURE);
		}
		myfile.close();
	}
	catch (itk::ExceptionObject & excp)
	{
		logger.WriteError("Error caught during testing: " + std::string(excp.GetDescription()));
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
