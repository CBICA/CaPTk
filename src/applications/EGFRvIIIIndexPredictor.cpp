#include "EGFRvIIIIndexPredictor.h"
#include "itkCSVNumericObjectFileWriter.h"
#include "CaPTkGUIUtils.h"
#include "CaPTkClassifierUtils.h"


typedef itk::Image< float, 3 > ImageType;
//EGFRvIIIIndexPredictor::~EGFRvIIIIndexPredictor()
//{
//}

VectorDouble EGFRvIIIIndexPredictor::GetStatisticalFeatures(const VectorDouble &intensities)
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

VectorDouble EGFRvIIIIndexPredictor::GetHistogramFeatures(const VectorDouble &intensities, const double &start, const double &interval, const double &end)
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

VectorDouble EGFRvIIIIndexPredictor::GetVolumetricFeatures(const double &edemaSize, const double &tuSize, const double &neSize, const double &totalSize)
{
  VectorDouble VolumetricFeatures;
  VolumetricFeatures.push_back(tuSize);
  VolumetricFeatures.push_back(neSize);
  VolumetricFeatures.push_back(edemaSize);
  VolumetricFeatures.push_back(totalSize);

  VolumetricFeatures.push_back(tuSize + neSize);
  VolumetricFeatures.push_back((tuSize + neSize)*100 / totalSize);
  VolumetricFeatures.push_back(edemaSize*100 / totalSize);

  VolumetricFeatures.push_back(tuSize*100 / (tuSize+neSize));
  VolumetricFeatures.push_back(neSize * 100 / (tuSize + neSize));
  VolumetricFeatures.push_back(edemaSize* 100 / (tuSize + neSize));

  return VolumetricFeatures;
}

int EGFRvIIIIndexPredictor::PrepareNewEGFRvIIIPredictionModel(const std::string &inputdirectory, const std::vector< std::map< CAPTK::ImageModalityType, std::string > > &qualifiedsubjects, const std::string &outputdirectory)
{
  VectorDouble AllSurvival;
  VariableSizeMatrixType FeaturesOfAllSubjects;
  VariableSizeMatrixType HistogramFeaturesConfigurations;
  HistogramFeaturesConfigurations.SetSize(33, 3); //11 modalities*3 regions = 33 configurations*3 histogram features for each configuration
  CSVFileReaderType::Pointer reader = CSVFileReaderType::New();
  MatrixType dataMatrix;
  try
  {
    reader->SetFileName(getCaPTkDataDir() + "/egfrv3/EGFRvIII_HMFeatures_Configuration.csv");
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
    logger.WriteError("Cannot find the file 'EGFRvIII_HMFeatures_Configuration.csv' in the ../data/egfrv3 directory. Error code : " + std::string(e1.what()));
    return false;
  }
  //---------------------------------------------------------------------------
  ImageType::Pointer NEGAtlasImagePointer = ReadNiftiImage<ImageType>(getCaPTkDataDir() + "/egfrv3/EGFRneg.nii.gz");
  ImageType::Pointer POSAtlasImagePointer = ReadNiftiImage<ImageType>(getCaPTkDataDir() + "/egfrv3/EGFRpos.nii.gz");

  FeaturesOfAllSubjects.SetSize(qualifiedsubjects.size(), 448);
  for (unsigned int sid = 0; sid < qualifiedsubjects.size(); sid++)
  {
    std::cout << "Patient's data loaded:" << sid + 1 << std::endl;
    std::map< CAPTK::ImageModalityType, std::string > currentsubject = qualifiedsubjects[sid];
    try
    {
      ImageType::Pointer LabelImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG]));
      ImageType::Pointer AtlasImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_ATLAS]));
      ImageType::Pointer TemplateImagePointer9Regions = ReadNiftiImage<ImageType>(getCaPTkDataDir() + "/egfrv3/template9regions.nii.gz");
      ImageType::Pointer TemplateImagePointer21Regions = ReadNiftiImage<ImageType>(getCaPTkDataDir() + "/egfrv3/template21regions.nii.gz");
      ImageType::Pointer TemplateImagePointerAllRegions = ReadNiftiImage<ImageType>(getCaPTkDataDir() + "/egfrv3/templateallregions.nii.gz");

      ImageType::Pointer RCBVImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_RCBV])));
      ImageType::Pointer PHImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PH])));
      ImageType::Pointer PSRImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PSR])));

      ImageType::Pointer T1CEImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1CE])));
      ImageType::Pointer T2FlairImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR])));
      ImageType::Pointer T1ImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1])));
      ImageType::Pointer T2ImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2])));

      ImageType::Pointer AXImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_AX])));
      ImageType::Pointer RADImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_RAD])));
      ImageType::Pointer FAImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_FA])));
      ImageType::Pointer TRImagePointer = RescaleImageIntensity<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_TR])));

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
        RCBVImagePointer, PSRImagePointer, PHImagePointer, AXImagePointer, FAImagePointer, RADImagePointer, TRImagePointer, LabelImagePointer, AtlasImagePointer, TemplateImagePointer9Regions, TemplateImagePointer21Regions, TemplateImagePointerAllRegions, POSAtlasImagePointer, NEGAtlasImagePointer, HistogramFeaturesConfigurations);

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
  std::cout << std::endl << "Features loaded....." << std::endl;
  VariableSizeMatrixType scaledFeatureSet;
  scaledFeatureSet.SetSize(qualifiedsubjects.size(), 161);
  VariableLengthVectorType meanVector;
  VariableLengthVectorType stdVector;
  mFeatureScalingLocalPtr.ScaleGivenTrainingFeatures(FeaturesOfAllSubjects, scaledFeatureSet, meanVector, stdVector);

  std::cout << std::endl << "Scaling done....." << std::endl;
  //typedef itk::CSVNumericObjectFileWriter<double, 20, 448> WriterTypeMatrix;
  //WriterTypeMatrix::Pointer writermatrix = WriterTypeMatrix::New();
  MatrixType data;
  /*data.set_size(20, 448);
  for (int i = 0; i < 20; i++)
    for (int j = 0; j < 448; j++)
      data(i, j) = FeaturesOfAllSubjects(i, j);
  writermatrix->SetFileName(outputdirectory + "/raw_training_features.csv");
  writermatrix->SetInput(&data);
  writermatrix->Write();

  for (int i = 0; i < 20; i++)
    for (int j = 0; j < 448; j++)
      data(i, j) = scaledFeatureSet(i, j);
  writermatrix->SetFileName(outputdirectory + "/scaled_training_features.csv");
  writermatrix->SetInput(&data);
  writermatrix->Write();
*/

  VariableSizeMatrixType ModelFeatures;
  mFeatureExtractionLocalPtr.FormulateEGFRTrainingData(scaledFeatureSet, AllSurvival, ModelFeatures);
  
  std::cout << std::endl << "Data formulation finished....." << std::endl;

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
  /*VariableSizeMatrixType SixModelFeatures;
  VariableSizeMatrixType EighteenModelFeatures;
  mFeatureExtractionLocalPtr.FormulateSurvivalTrainingData(scaledFeatureSet, AllSurvival, SixModelFeatures, EighteenModelFeatures);*/
  //MatrixType data;
  try
  {
    data.set_size(448, 1); // TOCHECK - are these hard coded sizes fine?
    for (unsigned int i = 0; i < meanVector.Size(); i++)
      data(i, 0) = meanVector[i];
    typedef itk::CSVNumericObjectFileWriter<double, 448, 1> WriterTypeVector;
    WriterTypeVector::Pointer writerv = WriterTypeVector::New();
    writerv->SetFileName(outputdirectory + "/EGFRvIII_ZScore_Mean.csv");
    writerv->SetInput(&data);
    writerv->Write();

    for (unsigned int i = 0; i < stdVector.Size(); i++)
      data(i, 0) = stdVector[i];
    writerv->SetFileName(outputdirectory + "/EGFRvIII_ZScore_Std.csv");
    writerv->SetInput(&data);
    writerv->Write();
  }
  catch (const std::exception& e1)
  {
    std::cout << "scaling writing error" << std::endl;
    logger.WriteError("Error in writing output files to the output directory = " + outputdirectory + "Error code : " + std::string(e1.what()));
    return false;
  }
  std::cout << std::endl << "Building model....." << std::endl;
  //  //---------------------------------------------------------------------------
  VariableSizeMatrixType ModelSelectedFeatures = SelectModelFeatures(ModelFeatures);
  std::cout << std::endl << "Building model....." << std::endl;

  //typedef itk::CSVNumericObjectFileWriter<double, 20, 8> WriterTypeMatrix2;
  //WriterTypeMatrix2::Pointer writermatrix2 = WriterTypeMatrix2::New();
  //data.set_size(20, 8);
  //for (int i = 0; i < 20; i++)
  //  for (int j = 0; j < 8; j++)
  //    data(i, j) = ModelSelectedFeatures(i, j);
  //writermatrix2->SetFileName(outputdirectory + "/selected_training_features.csv");
  //writermatrix2->SetInput(&data);
  //writermatrix2->Write();


  std::cout << std::endl << "Building model....." << std::endl;

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
    trainOpenCVSVM(ModelSelectedFeatures, outputdirectory + "/" + mTrainedFile, false, CAPTK::ApplicationCallingSVM::Survival);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Training on the given subjects failed. Error code : " + std::string(e1.what()));
    return false;
  }
  std::cout << std::endl << "Model saved to the output directory." << std::endl;

  return true;
}

VariableLengthVectorType EGFRvIIIIndexPredictor::DistanceFunctionLinear(const VariableSizeMatrixType &testData, const std::string &filename, const double &rho, const double &bestg)
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
  VariableSizeMatrixType TransposedSupportVectors= MatrixTranspose(SupportVectors);
  VariableSizeMatrixType w;
  w.SetSize(TransposedSupportVectors.Rows(),1);
  for (unsigned int svID = 0; svID < TransposedSupportVectors.Rows(); svID++)
  {
    double currentSum = 0;
    for (unsigned int iterator = 0; iterator < TransposedSupportVectors.Cols(); iterator++)
      currentSum = currentSum + TransposedSupportVectors(svID, iterator)*Coefficients[iterator];
    w(svID,0) = currentSum;
  }
  VariableSizeMatrixType wTranspose = MatrixTranspose(w);  //1x7   1x7


  for (unsigned int patID = 0; patID < testData.Rows(); patID++)
  {
    double distance = 0;
    for (unsigned int svID = 0; svID < wTranspose.Cols(); svID++)
      distance = distance+wTranspose(0,svID)*testData(patID, svID);

    Distances[patID] = distance - rho;
  }
  return Distances;
}

VectorDouble EGFRvIIIIndexPredictor::CombineEstimates(const VariableLengthVectorType &estimates1, const VariableLengthVectorType &estimates2)
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


VectorDouble EGFRvIIIIndexPredictor::CombineEstimates(const VectorDouble &estimates1, const VectorDouble &estimates2)
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

VectorDouble EGFRvIIIIndexPredictor::EGFRvIIIPredictionOnExistingModel(const std::string &modeldirectory, const std::string &inputdirectory, const std::vector < std::map < CAPTK::ImageModalityType, std::string>> &qualifiedsubjects, const std::string &outputdirectory)
{
  std::cout << "Started reading model parameters" << std::endl;
  //std::string local_captk_dataDir = "E:/SoftwareDevelopmentProjects/CaPTk-August2018/data";
  typedef itk::CSVArray2DFileReader<double> ReaderType;
  VariableSizeMatrixType HistogramFeaturesConfigurations;
  HistogramFeaturesConfigurations.SetSize(33, 3); //11 modalities*3 regions = 33 configurations*3 histogram features for each configuration
  VectorDouble results;

  CSVFileReaderType::Pointer reader = CSVFileReaderType::New();
  VectorDouble ages;
  MatrixType dataMatrix;
  try
  {
    reader->SetFileName(getCaPTkDataDir() +"egfrv3/EGFRvIII_HMFeatures_Configuration.csv");
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
    logger.WriteError("Cannot find the file 'EGFRvIII_HMFeatures_Configuration.csv' in the ../data/egfrv3 directory. Error code : " + std::string(e1.what()));
    return results;
  }

  MatrixType meanMatrix;
  VariableLengthVectorType mean;
  VariableLengthVectorType stddevition;
  try
  {
    reader->SetFileName(modeldirectory + "/EGFRvIII_ZScore_Mean.csv");
    reader->SetFieldDelimiterCharacter(',');
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();
    reader->Parse();
    meanMatrix = reader->GetArray2DDataObject()->GetMatrix();

    mean.SetSize(meanMatrix.size());
    for (unsigned int i = 0; i < meanMatrix.size(); i++)
      mean[i] = meanMatrix(0,i);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Error in reading the file: " + modeldirectory + "/EGFRvIII_ZScore_Mean.csv. Error code : " + std::string(e1.what()));
    return results;
  }
  MatrixType stdMatrix;
  try
  {
    reader->SetFileName(modeldirectory + "/EGFRvIII_ZScore_Std.csv");
    reader->SetFieldDelimiterCharacter(',');
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();
    reader->Parse();
    stdMatrix = reader->GetArray2DDataObject()->GetMatrix();

    stddevition.SetSize(stdMatrix.size());
    for (unsigned int i = 0; i < stdMatrix.size(); i++)
      stddevition[i] = stdMatrix(0,i);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Error in reading the file: " + modeldirectory + "/EGFRvIII_ZScore_Std.csv. Error code : " + std::string(e1.what()));
    return results;
  }
  //----------------------------------------------------
  ImageType::Pointer NEGAtlasImagePointer = ReadNiftiImage<ImageType>(getCaPTkDataDir() + "egfrv3/EGFRneg.nii.gz");
  ImageType::Pointer POSAtlasImagePointer = ReadNiftiImage<ImageType>(getCaPTkDataDir() + "egfrv3/EGFRpos.nii.gz");

  std::cout << "Finished readig model parameters" << std::endl;
  VariableSizeMatrixType FeaturesOfAllSubjects;
  FeaturesOfAllSubjects.SetSize(qualifiedsubjects.size(), 448);

  for (unsigned int sid = 0; sid < qualifiedsubjects.size(); sid++)
  {
    std::cout << "Subject No:" << sid << std::endl;
    std::map<CAPTK::ImageModalityType, std::string> currentsubject = qualifiedsubjects[sid];
    try
    {
      ImageType::Pointer LabelImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG]));
      ImageType::Pointer AtlasImagePointer = ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_ATLAS]));
      ImageType::Pointer TemplateImagePointer9Regions = ReadNiftiImage<ImageType>(getCaPTkDataDir() + "egfrv3/template9regions.nii.gz");
      ImageType::Pointer TemplateImagePointer21Regions = ReadNiftiImage<ImageType>(getCaPTkDataDir() + "egfrv3/template21regions.nii.gz");
      ImageType::Pointer TemplateImagePointerAllRegions = ReadNiftiImage<ImageType>(getCaPTkDataDir() + "egfrv3/templateallregions.nii.gz");
      ImageType::Pointer T1CEImagePointer = RescaleImageIntensity<ImageType>(CapImageIntensityWithPercentile<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1CE]))));
      ImageType::Pointer RCBVImagePointer = RescaleImageIntensity<ImageType>(CapImageIntensityWithPercentile<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_RCBV]))));
      ImageType::Pointer PHImagePointer = RescaleImageIntensity<ImageType>(CapImageIntensityWithPercentile<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PH]))));
      ImageType::Pointer T2FlairImagePointer = RescaleImageIntensity<ImageType>(CapImageIntensityWithPercentile<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR]))));
      ImageType::Pointer T1ImagePointer = RescaleImageIntensity<ImageType>(CapImageIntensityWithPercentile<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1]))));
      ImageType::Pointer T2ImagePointer = RescaleImageIntensity<ImageType>(CapImageIntensityWithPercentile<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2]))));
      ImageType::Pointer AXImagePointer = RescaleImageIntensity<ImageType>(CapImageIntensityWithPercentile<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_AX]))));
      ImageType::Pointer RADImagePointer = RescaleImageIntensity<ImageType>(CapImageIntensityWithPercentile<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_RAD]))));
      ImageType::Pointer FAImagePointer = RescaleImageIntensity<ImageType>(CapImageIntensityWithPercentile<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_FA]))));
      ImageType::Pointer TRImagePointer = RescaleImageIntensity<ImageType>(CapImageIntensityWithPercentile<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_TR]))));
      ImageType::Pointer PSRImagePointer = RescaleImageIntensity<ImageType>(CapImageIntensityWithPercentile<ImageType>(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PSR]))));

      VectorDouble TestFeatures = LoadTestData<ImageType>(T1CEImagePointer, T2FlairImagePointer, T1ImagePointer, T2ImagePointer,
        RCBVImagePointer, PSRImagePointer, PHImagePointer, AXImagePointer, FAImagePointer, RADImagePointer, TRImagePointer, LabelImagePointer, AtlasImagePointer, TemplateImagePointer9Regions, TemplateImagePointer21Regions, TemplateImagePointerAllRegions,POSAtlasImagePointer, NEGAtlasImagePointer, HistogramFeaturesConfigurations);


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
  //std::cout << "Feature writing started:" << std::endl;
  //typedef itk::CSVNumericObjectFileWriter<double, 144, 448> WriterTypeMatrix;
  //WriterTypeMatrix::Pointer writermatrix = WriterTypeMatrix::New();
  //MatrixType data;
  //data.set_size(144,448);
  //for (int i = 0; i < 144; i++)
  //	for (int j = 0; j < 448; j++)
  //		data(i, j) = FeaturesOfAllSubjects(i, j);
  //writermatrix->SetFileName(outputdirectory+ "/plain_test_features.csv");
  //writermatrix->SetInput(&data);
  //writermatrix->Write();
  

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

  //typedef itk::CSVNumericObjectFileWriter<double, 144, 449> WriterTypeMatrix1;
  //WriterTypeMatrix1::Pointer writermatrix1 = WriterTypeMatrix1::New();
  //data.set_size(144, 449);
  //for (int i = 0; i < 144; i++)
  //  for (int j = 0; j < 449; j++)
  //    data(i, j) = ScaledFeatureSetAfterAddingLabel(i, j);
  //writermatrix->SetFileName(outputdirectory + "/scaled_test_features.csv");
  //writermatrix->SetInput(&data);
  //writermatrix->Write();

  VariableSizeMatrixType ModelSelectedFeatures = SelectModelFeatures(ScaledFeatureSetAfterAddingLabel);

  //typedef itk::CSVNumericObjectFileWriter<double, 144, 8> WriterTypeMatrix2;
  //WriterTypeMatrix2::Pointer writermatrix2 = WriterTypeMatrix2::New();
  //data.set_size(144, 8);
  //for (int i = 0; i < 144; i++)
  //  for (int j = 0; j < 8; j++)
  //    data(i, j) = ModelSelectedFeatures(i, j);
  //writermatrix2->SetFileName(outputdirectory + "/selected_test_features.csv");
  //writermatrix2->SetInput(&data);
  //writermatrix2->Write();

  //---------------------------------------------------------------------------------------------------------------	
  try
  {
    std::ofstream myfile;
    myfile.open(outputdirectory + "/results.csv");
    myfile << "SubjectName,Score, Result \n";
    if (cbica::fileExists(modeldirectory + "/EGFRvIII_SVM_Model.csv") == true)
    {
      VariableLengthVectorType result;
      result = DistanceFunctionLinear(ModelSelectedFeatures, modeldirectory + "/EGFRvIII_SVM_Model.csv", -0.9925,2);
      for (size_t i = 0; i < result.Size(); i++)
      {
        std::map<CAPTK::ImageModalityType, std::string> currentsubject = qualifiedsubjects[i];
        results.push_back(-1*result[i]);
        if(results[i]<0)
            myfile << static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]) + "," + std::to_string(results[i]) + ", Wildtype \n";
        else
          myfile << static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]) + "," + std::to_string(results[i]) + ", Mutant \n";
      }
    }
    else if (cbica::fileExists(modeldirectory + "/EGFRvIII_SVM_Model.xml") == true)
    {
      VectorDouble result;
      result = testOpenCVSVM(ModelSelectedFeatures, modeldirectory + "/EGFRvIII_SVM_Model.xml");
      results = result;
      for (size_t i = 0; i < result.size(); i++)
      {
        std::map<CAPTK::ImageModalityType, std::string> currentsubject = qualifiedsubjects[i];
        if (result[i]<0)
          myfile << static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]) + "," + std::to_string(result[i]) + ", Wildtype \n";
        else
          myfile << static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]) + "," + std::to_string(result[i]) + ", Mutant \n";
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

VariableSizeMatrixType EGFRvIIIIndexPredictor::SelectModelFeatures(const VariableSizeMatrixType &ModelFeatures)
{
  int selectedFeatures[7] = { 1,2,3,4,54,420,421};

  VariableSizeMatrixType ModelSelectedFeatures;
  ModelSelectedFeatures.SetSize(ModelFeatures.Rows(), 8);
  int counter = 0;
  for (unsigned int i = 0; i < 7; i++)
  {
    for (unsigned int j = 0; j < ModelFeatures.Rows(); j++)
      ModelSelectedFeatures(j, counter) = ModelFeatures(j, selectedFeatures[i]);
    counter++;
  }
  for (unsigned int j = 0; j < ModelFeatures.Rows(); j++)
    ModelSelectedFeatures(j, 7) = ModelFeatures(j, 448);

  return ModelSelectedFeatures;
}

VariableSizeMatrixType EGFRvIIIIndexPredictor::MatrixTranspose(const VariableSizeMatrixType &inputmatrix)
{
  VariableSizeMatrixType output;
  output.SetSize(inputmatrix.Cols(), inputmatrix.Rows());

  for (unsigned int i = 0; i < output.Rows(); i++)
    for (unsigned int j = 0; j < output.Cols(); j++)
      output(i, j) = inputmatrix(j, i);
  return output;
}

VectorDouble EGFRvIIIIndexPredictor::NormalizeHistogramFeatures(VectorDouble inputHistogramFeatures,const int size)
{
  VectorDouble result;
  for (int i = 0; i < inputHistogramFeatures.size(); i++)
    result.push_back(inputHistogramFeatures[i]*100 / size);
  return result;
}

