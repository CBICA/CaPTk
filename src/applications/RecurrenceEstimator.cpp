#include "RecurrenceEstimator.h"
#include "fMainWindow.h"


typedef itk::Image< float, 3 > ImageType;

RecurrenceEstimator::~RecurrenceEstimator()
{
  //delete mNiftiLocalPtr;
  //delete mOutputLocalPtr;
  //delete mFeatureReductionLocalPtr;
  //delete mFeatureScalingLocalPtr;
  //delete mFeatureExtractionLocalPtr;
}

ImageTypeFloat3D::Pointer RecurrenceEstimator::RescaleImageIntensity(ImageTypeFloat3D::Pointer image)
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


bool RecurrenceEstimator::TrainNewModelOnGivenData(const std::vector<std::map<CAPTK::ImageModalityType, std::string>> qualifiedsubjects, const std::string &outputdirectory, bool useConventionalData, bool useDTIData, bool usePerfData, bool useDistData)
{

	bool useOtherModalities = false;
	if (useConventionalData || useDTIData)
		useOtherModalities = true;

	std::vector<std::string> vector;
	VectorVectorDouble pNearIntensitiesTotal;
	VectorVectorDouble mNearIntensitiesTotal;
	VectorDouble tNearIntensitiesTotal;

	VectorVectorDouble pFarIntensitiesTotal;
	VectorVectorDouble mFarIntensitiesTotal;
	VectorDouble tFarIntensitiesTotal;

	VectorVectorDouble pIntensitiesTotal;
	VectorVectorDouble mIntensitiesTotal;
	VectorDouble tIntensitiesTotal;


	for (unsigned int sid = 0; sid < qualifiedsubjects.size(); sid++)
	{
		std::cout << "Patient's data laoding:" << sid + 1 << std::endl;
		VectorVectorDouble mNearIntensities;
		VectorVectorDouble mFarIntensities;
		VectorDouble tNearIntensities;
		VectorDouble tFarIntensities;
		VectorVectorDouble pNearIntensities;
		VectorVectorDouble pFarIntensities;

		typedef ImageTypeFloat3D ImageType;
		typedef ImageTypeFloat4D PerfusionImageType;
		ImageType::Pointer T1CEImagePointer;
		ImageType::Pointer T2FlairImagePointer;
		ImageType::Pointer T1ImagePointer;
		ImageType::Pointer T2ImagePointer;
		ImageType::Pointer AXImagePointer;
		ImageType::Pointer RADImagePointer;
		ImageType::Pointer FAImagePointer;
		ImageType::Pointer TRImagePointer;

		PerfusionImageType::Pointer perfImagePointer;

		std::map< CAPTK::ImageModalityType, std::string > currentsubject = qualifiedsubjects[sid];
		try
		{
			ImageType::Pointer LabelImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG]));
			ImageType::Pointer recurImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_NEAR]));
			ImageType::Pointer nonrecurImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_FAR]));

			if (useConventionalData)
			{
				T1CEImagePointer = RescaleImageIntensity(mNiftiLocalPtr.ReadNiftiImage(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1CE])));
				T2FlairImagePointer = RescaleImageIntensity(mNiftiLocalPtr.ReadNiftiImage(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR])));
				T1ImagePointer = RescaleImageIntensity(mNiftiLocalPtr.ReadNiftiImage(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1])));
				T2ImagePointer = RescaleImageIntensity(mNiftiLocalPtr.ReadNiftiImage(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2])));
			}
			if (useDTIData)
			{
				AXImagePointer = RescaleImageIntensity(mNiftiLocalPtr.ReadNiftiImage(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_AX])));
				RADImagePointer = RescaleImageIntensity(mNiftiLocalPtr.ReadNiftiImage(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_RAD])));
				FAImagePointer = RescaleImageIntensity(mNiftiLocalPtr.ReadNiftiImage(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_FA])));
				TRImagePointer = RescaleImageIntensity(mNiftiLocalPtr.ReadNiftiImage(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_TR])));
			}
			if (usePerfData)
				perfImagePointer = mNiftiLocalPtr.Read4DNiftiImage(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION]));

			int imagetype = CAPTK::ImageExtension::NIfTI;
			VectorVectorDouble tumorIndices;
			typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
			IteratorType LabelImageIt(LabelImagePointer, LabelImagePointer->GetLargestPossibleRegion());
			LabelImageIt.GoToBegin();
			while (!LabelImageIt.IsAtEnd())
			{
				if (LabelImageIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::TUMOR || LabelImageIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::NONENHANCING)
				{
					std::vector<double> local;
					local.push_back(LabelImageIt.GetIndex()[0]);
					local.push_back(LabelImageIt.GetIndex()[1]);
					local.push_back(LabelImageIt.GetIndex()[2]);
					tumorIndices.push_back(local);
				}
				++LabelImageIt;
			}

			PreprocessingPipelineClass mPreprocessingObj;
			ImageType::Pointer tumorMaskFinal = mPreprocessingObj.PrepareTumroImageFromPoints<ImageType>(LabelImagePointer, tumorIndices);

			mNiftiLocalPtr.LoadTrainingData(tumorMaskFinal, recurImagePointer, nonrecurImagePointer,T1CEImagePointer, T2FlairImagePointer, T1ImagePointer,T2ImagePointer, perfImagePointer, AXImagePointer,FAImagePointer, RADImagePointer, TRImagePointer,pNearIntensities, pFarIntensities, mNearIntensities,mFarIntensities, tNearIntensities, tFarIntensities, imagetype,useConventionalData, useDTIData, usePerfData, useDistData);
			for (size_t i = 0; i < mNearIntensities.size(); i++)
			{
				if (usePerfData)
					pNearIntensitiesTotal.push_back(pNearIntensities[i]);
				if (useOtherModalities)
					mNearIntensitiesTotal.push_back(mNearIntensities[i]);
				if (useDistData)
					tNearIntensitiesTotal.push_back(tNearIntensities[i]);
			}
			for (size_t i = 0; i < mFarIntensities.size(); i++)
			{
				if (usePerfData)
					pFarIntensitiesTotal.push_back(pFarIntensities[i]);
				if (useOtherModalities)
					mFarIntensitiesTotal.push_back(mFarIntensities[i]);
				if (useDistData)
					tFarIntensitiesTotal.push_back(tFarIntensities[i]);
			}
		}
		catch (const std::exception& e1)
		{
			logger.WriteError("Error in calculating the features for patient ID = " + static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PSR]) + ". Error code : " + std::string(e1.what()));
			return false;
		}
	}
	int totalNearSize = pNearIntensitiesTotal.size();
	int totalFarSize = pFarIntensitiesTotal.size();
	int totalSize = totalNearSize + totalFarSize;

	for (size_t i = 0; i < pNearIntensitiesTotal.size(); i++)
	{
		pIntensitiesTotal.push_back(pNearIntensitiesTotal[i]);
		mIntensitiesTotal.push_back(mNearIntensitiesTotal[i]);
		tIntensitiesTotal.push_back(tNearIntensitiesTotal[i]);
	}
	for (size_t i = 0; i < pFarIntensitiesTotal.size(); i++)
	{
		pIntensitiesTotal.push_back(pFarIntensitiesTotal[i]);
		mIntensitiesTotal.push_back(mFarIntensitiesTotal[i]);
		tIntensitiesTotal.push_back(tFarIntensitiesTotal[i]);
	}

	VariableLengthVectorType perfMeanVector;
	vtkSmartPointer<vtkTable> reducedPerfusionFeatures;
	//------------------------------------------reduce perfusion intensities------------------------------------------
	std::cout << "Applying pCA." << std::endl;
	if (usePerfData)
	{
		perfMeanVector = mFeatureReductionLocalPtr.ComputeMeanOfGivenFeatureVectors(pIntensitiesTotal);
		reducedPerfusionFeatures = mFeatureReductionLocalPtr.GetDiscerningPerfusionTimePoints(pIntensitiesTotal);
	}

	//-----------------------------------develope final near and far vectors------------------------------------------

	VectorVectorDouble fNearIntensities;
	VectorVectorDouble fFarIntensities;
	for (int i = 0; i < totalNearSize; i++)
	{
		VectorDouble cIntensityVectorPerSub;
		if (usePerfData)
			for (int j = 0; j < NO_OF_PCS; j++)
				cIntensityVectorPerSub.push_back(reducedPerfusionFeatures->GetValue(i, j).ToDouble());
		if (useOtherModalities)
			for (unsigned int j = 0; j < mIntensitiesTotal[i].size(); j++)
				cIntensityVectorPerSub.push_back(mIntensitiesTotal[i][j]);
		if (useDistData)
			cIntensityVectorPerSub.push_back(tIntensitiesTotal[i]);
		fNearIntensities.push_back(cIntensityVectorPerSub);
	}
	for (int i = totalNearSize; i < totalSize; i++)
	{
		VectorDouble cIntensityVectorPerSub;
		if (usePerfData)
			for (int j = 0; j < NO_OF_PCS; j++)
				cIntensityVectorPerSub.push_back(reducedPerfusionFeatures->GetValue(i, j).ToDouble());
		if (useOtherModalities)
			for (unsigned int j = 0; j < mIntensitiesTotal[i].size(); j++)
				cIntensityVectorPerSub.push_back(mIntensitiesTotal[i][j]);
		if (useDistData)
			cIntensityVectorPerSub.push_back(tIntensitiesTotal[i]);
		fFarIntensities.push_back(cIntensityVectorPerSub);
	}

	//---------------------training data formulation-----------------------------------
	std::cout << "Training data formulation. Assigning class labels." << std::endl;
	mFeatureExtractionLocalPtr.FormulateTrainingData(fNearIntensities, fFarIntensities);
	VariableSizeMatrixType TrainingData = mFeatureExtractionLocalPtr.GetTrainingData();

	//typedef vnl_matrix<double> MatrixType;
	//MatrixType data;

	//data.set_size(1404,6);
	//for (int i = 0; i < TrainingData.Rows(); i++)
	//	for (int j = 0; j < TrainingData.Cols(); j++)
	//		data(i, j) = TrainingData[i][j];
	//typedef itk::CSVNumericObjectFileWriter<double, 1404,6> WriterType;
	//WriterType::Pointer writer = WriterType::New();
	//writer->SetFileName("tData.csv");
	//writer->SetInput(&data);
	//writer->Write();
	std::cout << "Training data scaling (Z-Score)." << std::endl;
	VariableSizeMatrixType ScaledTrainingData = mFeatureScalingLocalPtr.ScaleGivenTrainingFeatures(TrainingData);
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

	std::cout << "Resampling training data. Reducing the size of majority class." << std::endl;
	VariableSizeMatrixType ResampledTrainingData = mFeatureExtractionLocalPtr.ResampleTrainingData(ScaledTrainingData, totalNearSize, totalFarSize);
	//------------------------------------------saving model related features---------------------------------------------------

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
	mOutputLocalPtr.SetOutputDirectoryPath(outputdirectory);
	int size = GetFeatureVectorSize(useConventionalData, useDTIData, usePerfData, useDistData);
	try
	{
		mOutputLocalPtr.SaveModelResults(ScaledTrainingData, mFeatureScalingLocalPtr.GetMeanVector(), mFeatureScalingLocalPtr.GetStdVector(), perfMeanVector, mFeatureReductionLocalPtr.GetPCATransformationMatrix(),useConventionalData, useDTIData, usePerfData, useDistData, size);
	}
	catch (const std::exception& e1)
	{
		logger.WriteError("Error in writing model files to the output directory: " + outputdirectory + ". Error code : " + std::string(e1.what()));
		return false;
	}

	try
	{
		std::cout << "Building SVM model." << std::endl;
		trainOpenCVSVM(ResampledTrainingData, outputdirectory + "/" + mTrainedModelNameXML, true, CAPTK::ApplicationCallingSVM::Recurrence);
	}
	catch (const std::exception& e1)
	{
		logger.WriteError("Training on the subjects failed. Error code : " + std::string(e1.what()));
		return false;
	}
	mFeatureReductionLocalPtr.ResetParameters();
  mFeatureScalingLocalPtr.ResetParameters(); 
  return true;
}

int RecurrenceEstimator::GetFeatureVectorSize(bool &useConventionalData, bool &useDTIData, bool &usePerfData, bool &useDistData)
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

bool RecurrenceEstimator::RecurrenceEstimateOnExistingModel(std::vector<std::map<CAPTK::ImageModalityType, std::string>> qualifiedsubjects,
  const std::string &modeldirectory,
  const std::string &inputdirectory,
  const std::string &outputdirectory,
  bool useConventionalrData,
  bool useDTIData, bool usePerfData, bool useDistData)
{

  bool useOtherModalities = false;
  if (useConventionalrData || useDTIData)
    useOtherModalities = true;

  // int size = GetFeatureVectorSize(useT1Data, useT1CEData, useT2FlairData, useT2Data, useDTIData, usePerfData, useDistData);
  int imagetype;

  VariableSizeMatrixType pca_coefficients;
  VariableLengthVectorType pca_mean;
  VariableLengthVectorType mean;
  VariableLengthVectorType stds;

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
  try
  {
	  mOutputLocalPtr.ReadModelParameters(modeldirectory + "/Recurrence_ZScore_Mean.csv", modeldirectory + "/Recurrence_ZScore_Std.csv", modeldirectory + "/Recurrence_COEF.csv", modeldirectory + "/Recurrence_MR.csv", mean, stds, pca_coefficients, pca_mean);
	  mFeatureReductionLocalPtr.SetParameters(pca_coefficients, pca_mean);
	  mFeatureScalingLocalPtr.SetParameters(mean, stds);
	  //  }
	  //}
	  //else
	  //{
	  //  mOutputLocalPtr.ReadModelParameters(modeldirectory + "/Recurrence_ZScore_Mean.csv", modeldirectory + "/Recurrence_ZScore_Std.csv", mean, stds);
	  //  mFeatureScalingLocalPtr.SetParameters(mean, stds);
	  //}
  }
  catch (const std::exception& e1)
  {
	  logger.WriteError("Error in reading ZScore and PCA parameters from directory: " + modeldirectory + ". Error code : " + std::string(e1.what()));
	  return false;
  }

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

  for (unsigned int sid = 0; sid < qualifiedsubjects.size(); sid++)
  {
    std::map<CAPTK::ImageModalityType, std::string > currentsubject = qualifiedsubjects[sid];
    VectorVectorDouble perfusionIntensities;
    VectorVectorDouble otherIntensities;
    VectorDouble distanceIntensities;
    std::vector<ImageType::IndexType> testindices;


    typedef ImageTypeFloat3D ImageType;
    typedef ImageTypeFloat4D PerfusionImageType;
    ImageType::Pointer T1CEImagePointer;
    ImageType::Pointer T2FlairImagePointer;
    ImageType::Pointer T1ImagePointer;
    ImageType::Pointer T2ImagePointer;
    ImageType::Pointer AXImagePointer;
    ImageType::Pointer RADImagePointer;
    ImageType::Pointer FAImagePointer;
    ImageType::Pointer TRImagePointer;
	ImageType::Pointer LabelImagePointer;
    PerfusionImageType::Pointer perfImagePointer;
	ImageType::Pointer dilatedEdema;
	try
	{
		LabelImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG]));
		if (usePerfData)
			perfImagePointer = mNiftiLocalPtr.Read4DNiftiImage(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION]));
		if (useConventionalrData)
		{
			T1CEImagePointer = RescaleImageIntensity(mNiftiLocalPtr.ReadNiftiImage(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1CE])));
			T2FlairImagePointer = RescaleImageIntensity(mNiftiLocalPtr.ReadNiftiImage(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR])));
			T1ImagePointer = RescaleImageIntensity(mNiftiLocalPtr.ReadNiftiImage(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1])));
			T2ImagePointer = RescaleImageIntensity(mNiftiLocalPtr.ReadNiftiImage(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2])));
		}
		if (useDTIData)
		{
			AXImagePointer = RescaleImageIntensity(mNiftiLocalPtr.ReadNiftiImage(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_AX])));
			RADImagePointer = RescaleImageIntensity(mNiftiLocalPtr.ReadNiftiImage(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_RAD])));
			FAImagePointer = RescaleImageIntensity(mNiftiLocalPtr.ReadNiftiImage(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_FA])));
			TRImagePointer = RescaleImageIntensity(mNiftiLocalPtr.ReadNiftiImage(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_TR])));
		}
		//cbica::Logging(loggerFile, "Image loading finished.");
		//typedef ImageTypeFloat3D OutputImageType;
		//OutputImageType::Pointer dilatedEdema;
		//typedef itk::BinaryBallStructuringElement<ImageType::PixelType, 3> StructuringElementType;
		//StructuringElementType structuringElement;
		//structuringElement.SetRadius(1);
		//structuringElement.CreateStructuringElement();
		//typedef itk::BinaryDilateImageFilter <ImageType, ImageType, StructuringElementType> BinaryDilateImageFilterType;
		//BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
		//dilateFilter->SetInput(LabelImagePointer);
		//dilateFilter->SetKernel(structuringElement);
		//dilateFilter->SetDilateValue(GLISTR_OUTPUT_LABELS::EDEMA);
		//dilateFilter->Update();
		//dilatedEdema = dilateFilter->GetOutput();


		//auto img = convertVtkToItk<float, 3>(mSlicerManagers[0]->mMask);
		std::vector<double> labels;
		labels.push_back(CAPTK::GLISTR_OUTPUT_LABELS::EDEMA);
		dilatedEdema = GetImageWithLabels<ImageType>(labels, LabelImagePointer);

		imagetype = CAPTK::ImageExtension::NIfTI;
		testindices = mNiftiLocalPtr.LoadTestData(T1CEImagePointer, T2FlairImagePointer, T1ImagePointer, T2ImagePointer, perfImagePointer, AXImagePointer, FAImagePointer, RADImagePointer, TRImagePointer, LabelImagePointer, dilatedEdema, perfusionIntensities, otherIntensities, distanceIntensities, imagetype,useConventionalrData, useDTIData, usePerfData, useDistData);
		cbica::Logging(loggerFile, "Test data  loading finished.");
	}
	catch (const std::exception& e1)
	{
		logger.WriteError("Error in extracting features for patient ID: " + static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]) + ". Error code : " + std::string(e1.what()));
		return false;
	}

    VectorVectorDouble reducedPerfusionFeatures;
    if (usePerfData)
      reducedPerfusionFeatures = mFeatureReductionLocalPtr.ApplyPCAOnTestData(perfusionIntensities);

    int NumberOfPCs = 5;
    VectorVectorDouble globaltestintensities;

    for (unsigned int k = 0; k < testindices.size(); k++)
    {
      VectorDouble inten;

      if (usePerfData)
        for (int j = 0; j < NumberOfPCs; j++)
          inten.push_back(reducedPerfusionFeatures[k][j]);

      if (useOtherModalities)
        for (unsigned int j = 0; j < otherIntensities[0].size(); j++)
          inten.push_back(otherIntensities[k][j]);

      if (useDistData)
        inten.push_back(distanceIntensities[k]);

      if (inten.size() > 0)
        globaltestintensities.push_back(inten);
    }
    VariableSizeMatrixType TestingData = mFeatureExtractionLocalPtr.FormulateTestData(globaltestintensities);

    //typedef vnl_matrix<double> MatrixType;
    //MatrixType data;
    //data.set_size(98721, 15);

    //for (unsigned int i = 0; i < TestingData.Rows(); i++)
    //	for (unsigned int j = 0; j < TestingData.Cols(); j++)
    //		data(i, j) = TestingData[i][j];
    //typedef itk::CSVNumericObjectFileWriter<double, 98721, 15> WriterType;
    //WriterType::Pointer writer = WriterType::New();
    //writer->SetFileName("testdata.csv");
    //writer->SetInput(&data);
    //try
    //{
    //	writer->Write();
    //}
    //catch (itk::ExceptionObject & excp)
    //{
    //    std::cerr << "Error: " << excp.what();
    //}


    VariableSizeMatrixType ScaledTestingData = mFeatureScalingLocalPtr.ScaleGivenTestingFeatures(TestingData);
    //MatrixType sdata;
    //sdata.set_size(98721, 15);
    //for (int i = 0; i < TestingData.Rows(); i++)
    //	for (int j = 0; j < TestingData.Cols(); j++)
    //		sdata(i, j) = ScaledTestingData[i][j];
    //writer->SetFileName("scaledtestdata.csv");
    //writer->SetInput(&sdata);
    //try
    //{
    //	writer->Write();
    //}
    //catch (itk::ExceptionObject & excp)
    //{
    //}



    ImageType::RegionType region = T1CEImagePointer->GetLargestPossibleRegion();
    ImageType::Pointer RecProbabilityMap = ImageType::New();
    RecProbabilityMap->SetRegions(region);
    RecProbabilityMap->Allocate();
    RecProbabilityMap->SetSpacing(T1CEImagePointer->GetSpacing());
    RecProbabilityMap->SetOrigin(T1CEImagePointer->GetOrigin());
    RecProbabilityMap->SetDirection(T1CEImagePointer->GetDirection());

    VectorDouble result_modified;
    typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
    IteratorType imIt(LabelImagePointer, LabelImagePointer->GetLargestPossibleRegion());
    IteratorType RecIt(RecProbabilityMap, RecProbabilityMap->GetLargestPossibleRegion());
    imIt.GoToBegin(); RecIt.GoToBegin();

    while (!imIt.IsAtEnd())
    {
      if (!(imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::EDEMA))
        RecIt.Set(CAPTK::VOXEL_STATUS::OFF);

      ++imIt;
      ++RecIt;
    }
    cbica::Logging(loggerFile, "Before testing.");
    try
    {
      if (cbica::fileExists(modeldirectory + "/" + mTrainedModelNameCSV))
      {
        cbica::Logging(loggerFile, "Before testing 1.");
        VariableLengthVectorType result;
        result = DistanceFunction(ScaledTestingData, modeldirectory + "/" + mTrainedModelNameCSV, RECURRENCE_MODEL_RHO, RECURRENCE_MODEL_G);
        for (unsigned int index = 0; index < result.Size(); index++)
        {
          RecProbabilityMap->SetPixel(testindices[index], result[index] * 1);
          result_modified.push_back(result[index]);
        }

      }
      else
      {
        cbica::Logging(loggerFile, "Before testing 2.");
        VectorDouble result;
        result = testOpenCVSVM(ScaledTestingData, modeldirectory + "/" + mTrainedModelNameXML);
        for (unsigned int index = 0; index < result.size(); index++)
          RecProbabilityMap->SetPixel(testindices[index], result[index] * 1);
        result_modified = result;
      }
    }
    catch (itk::ExceptionObject & excp)
    {
      cbica::Logging(loggerFile, "Error caught during testing: " + std::string(excp.GetDescription()));
	  return false;
    }


    //---------------------------post-processing------------------------------------------------------------------------------------
	try
	{
		VectorDouble result_revised = RecurrenceMapPostprocessing<ImageType>(result_modified, testindices, RecProbabilityMap, dilatedEdema);
		for (unsigned int index = 0; index < result_modified.size(); index++)
			RecProbabilityMap->SetPixel(testindices[index], result_revised[index] * 1);

		////averaging filter
		//typedef itk::MeanImageFilter<ImageType, ImageType > FilterType;
		//FilterType::Pointer meanFilter = FilterType::New();
		//FilterType::InputSizeType radius;
		//radius.Fill(1);
		//meanFilter->SetRadius(radius);
		//meanFilter->SetInput(RecProbabilityMap);
		//ImageType::Pointer RevisedRecurrenceMap = meanFilter->GetOutput();
	
    ImageType::Pointer RevisedRecurrenceMap = RecurrenceMapPostprocessingForBackground<ImageType>(RecProbabilityMap, dilatedEdema);
		if (imagetype == CAPTK::ImageExtension::NIfTI)
			mOutputLocalPtr.WriteRecurrenceOutputInNifti<ImageType>(RevisedRecurrenceMap, outputdirectory + "/" + static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]));
		else
			mOutputLocalPtr.WriteRecurrenceOutputInNifti<ImageType>(RevisedRecurrenceMap, outputdirectory + "/" + static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]));
	}
	catch (itk::ExceptionObject & excp)
	{
		cbica::Logging(loggerFile, "Error caught during post-processing of recurrence map: " + std::string(excp.GetDescription()));
		return false;
	}
    //----------------------------------------------------------------------------------------------------------------------------------

  }
  mFeatureReductionLocalPtr.ResetParameters();
  mFeatureScalingLocalPtr.ResetParameters();
  return true;
}



VariableLengthVectorType RecurrenceEstimator::DistanceFunction(const VariableSizeMatrixType &testData, const std::string &filename, const double &rho, const double &bestg)
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
      distance = distance + result*Coefficients[svID];
    }
    Distances[patID] = distance - rho;
  }
  return Distances;
}



