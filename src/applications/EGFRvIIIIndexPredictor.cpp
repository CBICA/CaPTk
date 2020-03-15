#include "EGFRvIIIIndexPredictor.h"
#include "itkCSVNumericObjectFileWriter.h"
#include "CaPTkGUIUtils.h"
#include "CaPTkClassifierUtils.h"
#include "CaPTkUtils.h"
#include "cbicaUtilities.h"
#include "cbicaITKSafeImageIO.h"
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
  VectorDouble AllLabels;
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
  std::vector<int> GroundtruthLabelsAllPatients;
  std::vector<ImageType::Pointer> AtlasSegmentationsAllPatients;
  for (unsigned int sid = 0; sid < qualifiedsubjects.size(); sid++)
  {
    std::cout << "Patient's data loaded:" << sid + 1 << std::endl;
    std::map< CAPTK::ImageModalityType, std::string > currentsubject = qualifiedsubjects[sid];
    try
    {
      AtlasSegmentationsAllPatients.push_back(ReadNiftiImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_ATLAS])));
    }
    catch (const std::exception& e1)
    {
      logger.WriteError("Error in calculating the features for patient ID = " + static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID]) + "Error code : " + std::string(e1.what()));
      return false;
    }
    reader->SetFileName(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_FEATURES]));
    reader->SetFieldDelimiterCharacter(',');
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();
    reader->Parse();
    dataMatrix = reader->GetArray2DDataObject()->GetMatrix();

    for (unsigned int i = 0; i < dataMatrix.rows(); i++)
      GroundtruthLabelsAllPatients.push_back(dataMatrix(i, 1));
  }
  //---------------------------------------------------------------------
  ImageType::Pointer NEGAtlasImagePointer = MakeAtlases(AtlasSegmentationsAllPatients, GroundtruthLabelsAllPatients, 0);
  ImageType::Pointer POSAtlasImagePointer = MakeAtlases(AtlasSegmentationsAllPatients, GroundtruthLabelsAllPatients, 1);
  cbica::WriteImage<ImageType>(NEGAtlasImagePointer,outputdirectory+ "/EGFRneg.nii.gz");
  cbica::WriteImage<ImageType>(POSAtlasImagePointer, outputdirectory + "/EGFRpos.nii.gz");

  //ImageType::Pointer NEGAtlasImagePointer = ReadNiftiImage<ImageType>(getCaPTkDataDir() + "/egfrv3/EGFRneg.nii.gz");
  //ImageType::Pointer POSAtlasImagePointer = ReadNiftiImage<ImageType>(getCaPTkDataDir() + "/egfrv3/EGFRpos.nii.gz");
    
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
      VectorDouble label;

      reader->SetFileName(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_FEATURES]));
      reader->SetFieldDelimiterCharacter(',');
      reader->HasColumnHeadersOff();
      reader->HasRowHeadersOff();
      reader->Parse();
      dataMatrix = reader->GetArray2DDataObject()->GetMatrix();

      for (unsigned int i = 0; i < dataMatrix.rows(); i++)
      {
        ages.push_back(dataMatrix(i, 0));
        if (dataMatrix(i, 1) == 0)
        {
          label.push_back(-1);
          AllLabels.push_back(-1);
        }
        else
        {
          label.push_back(1);
          AllLabels.push_back(-1);
        }
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
  std::cout << "Feature writing started:" << std::endl;
  WriteCSVFiles(FeaturesOfAllSubjects, outputdirectory + "/RawFeatures.csv");
  //------------------------------------------------------------------------------------
  VariableSizeMatrixType scaledFeatureSet;
  scaledFeatureSet.SetSize(FeaturesOfAllSubjects.Rows(), FeaturesOfAllSubjects.Cols());
  VariableLengthVectorType meanVector;
  VariableLengthVectorType stdVector;
  mFeatureScalingLocalPtr.ScaleGivenTrainingFeatures(FeaturesOfAllSubjects, scaledFeatureSet, meanVector, stdVector);
  
  try
  {
    MatrixType data;
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

  std::cout << std::endl << "Scaling done....." << std::endl;
  std::cout << "Feature writing started:" << std::endl;
  WriteCSVFiles(scaledFeatureSet, outputdirectory + "/ScaledFeatures.csv");

  ////---------------------------------------------------------------------------
  //bool ModelSelectedFeatures = SelectModelFeaturesForTrainingSVMFFSel(scaledFeatureSet, AllLabels,outputdirectory);
  ////= SelectModelFeatures(ModelFeatures);


  //VariableSizeMatrixType ModelFeatures;
  //mFeatureExtractionLocalPtr.FormulateEGFRTrainingData(scaledFeatureSet, AllLabels, ModelFeatures);
  //
  //std::cout << std::endl << "Data formulation finished....." << std::endl;




  //---------------------------------------------------------------------------
  std::cout << std::endl << "Building model....." << std::endl;
  try
  {
    bool ModelSelectedFeatures = SelectModelFeaturesForTrainingSVMFFSel(scaledFeatureSet, AllLabels, outputdirectory);
//    trainOpenCVSVM(ModelSelectedFeatures, outputdirectory + "/" + mTrainedFile, false, CAPTK::ApplicationCallingSVM::Survival);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Training on the given subjects failed. Error code : " + std::string(e1.what()));
    return false;
  }
  std::cout << std::endl << "Model saved to the output directory." << std::endl;

  return true;
}

VariableLengthVectorType EGFRvIIIIndexPredictor::DistanceFunctionLinear(const VariableSizeMatrixType &testData, const std::string &filename)
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

  for (unsigned int i = 0; i < dataMatrix.rows(); i++)
  {
    unsigned int j = 0;
    for (j = 0; j < dataMatrix.cols() - 2; j++)
      SupportVectors(i, j) = dataMatrix(i, j);
    Coefficients[i] = dataMatrix(i, j);
    if (i == 0)
      Rho = dataMatrix(i, j + 1);
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

    Distances[patID] = distance - Rho;
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
  VariableLengthVectorType selectedfeatures;
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
  try
  {
    reader->SetFileName(modeldirectory + "/EGFRvIII_SelectedFeatures.csv");
    reader->SetFieldDelimiterCharacter(',');
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();
    reader->Parse();
    meanMatrix = reader->GetArray2DDataObject()->GetMatrix();

    selectedfeatures.SetSize(meanMatrix.size());
    for (unsigned int i = 0; i < meanMatrix.size(); i++)
      selectedfeatures[i] = meanMatrix(0, i);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Error in reading the file: " + modeldirectory + "/EGFRvIII_SelectedFeatures.csv. Error code : " + std::string(e1.what()));
    return results;
  }
/*
  std::vector<int> grroundtruth;
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(1);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(1);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(1);
  grroundtruth.push_back(1);
  grroundtruth.push_back(1);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(1);
  grroundtruth.push_back(1);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(1);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(1);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(1);
  grroundtruth.push_back(1);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(1);
  grroundtruth.push_back(1);
  grroundtruth.push_back(1);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
  grroundtruth.push_back(0);
*/
  //----------------------------------------------------
  ImageType::Pointer NEGAtlasImagePointer = ReadNiftiImage<ImageType>(modeldirectory+ "/EGFRneg.nii.gz");
  ImageType::Pointer POSAtlasImagePointer = ReadNiftiImage<ImageType>(modeldirectory + "/EGFRpos.nii.gz");

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
  std::cout << "Feature writing started:" << std::endl;
  WriteCSVFiles(FeaturesOfAllSubjects, outputdirectory + "/RawFeatures.csv");
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
  WriteCSVFiles(ScaledFeatureSetAfterAddingLabel, outputdirectory + "/ScaledFeatures.csv");
  //typedef itk::CSVNumericObjectFileWriter<double, 144, 449> WriterTypeMatrix1;
  //WriterTypeMatrix1::Pointer writermatrix1 = WriterTypeMatrix1::New();
  //data.set_size(144, 449);
  //for (int i = 0; i < 144; i++)
  //  for (int j = 0; j < 449; j++)
  //    data(i, j) = ScaledFeatureSetAfterAddingLabel(i, j);
  //writermatrix->SetFileName(outputdirectory + "/scaled_test_features.csv");
  //writermatrix->SetInput(&data);
  //writermatrix->Write();

  VariableSizeMatrixType ModelSelectedFeatures = SelectModelFeatures(ScaledFeatureSetAfterAddingLabel,selectedfeatures);
  WriteCSVFiles(ModelSelectedFeatures, outputdirectory + "/ScaledSelectedFeatures.csv");
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
      result = DistanceFunctionLinear(ModelSelectedFeatures, modeldirectory + "/EGFRvIII_SVM_Model.csv");
//      result = DistanceFunctionLinear(ModelSelectedFeatures, modeldirectory + "/EGFRvIII_SVM_Model.csv", -0.9092, 2);
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

VariableSizeMatrixType EGFRvIIIIndexPredictor::SelectModelFeatures(const VariableSizeMatrixType &ModelFeatures,const VariableLengthVectorType &selectedFeatures)
{
  VariableSizeMatrixType ModelSelectedFeatures;
  ModelSelectedFeatures.SetSize(ModelFeatures.Rows(), selectedFeatures.Size()+1);
  int counter = 0;
  for (unsigned int i = 0; i < selectedFeatures.Size(); i++)
  {
    for (unsigned int j = 0; j < ModelFeatures.Rows(); j++)
      ModelSelectedFeatures(j, counter) = ModelFeatures(j, selectedFeatures[i]);
    counter++;
  }
  for (unsigned int j = 0; j < ModelFeatures.Rows(); j++)
    ModelSelectedFeatures(j, selectedFeatures.Size()) = ModelFeatures(j, 448);

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


ImageType::Pointer EGFRvIIIIndexPredictor::MakeAtlases(const std::vector<ImageType::Pointer> &SegmentationaAllPatients, const std::vector<int> &LabelsAllPatientsconst, const int atlas_no)
{
    ImageType::Pointer atlasimage = ImageType::New();
    atlasimage->CopyInformation(SegmentationaAllPatients[0]);
    atlasimage->SetRequestedRegion(SegmentationaAllPatients[0]->GetLargestPossibleRegion());
    atlasimage->SetBufferedRegion(SegmentationaAllPatients[0]->GetBufferedRegion());
    atlasimage->Allocate();
    atlasimage->FillBuffer(0);

    for (int j = 0; j < LabelsAllPatientsconst.size(); j++)
    {
      if (LabelsAllPatientsconst[j] == atlas_no)
      {
        ImageType::Pointer currentImagePointer = SegmentationaAllPatients[j];
        typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
        IteratorType imIt(currentImagePointer, currentImagePointer->GetLargestPossibleRegion());
        IteratorType atlasIt(atlasimage, atlasimage->GetLargestPossibleRegion());
        imIt.GoToBegin();
        atlasIt.GoToBegin();

        while (!imIt.IsAtEnd())
        {
          if (imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::TUMOR || imIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::NONENHANCING)
            atlasIt.Set(atlasIt.Get() + 1);

          ++atlasIt;
          ++imIt;
        }
      }
    }
    return atlasimage;
}

bool EGFRvIIIIndexPredictor::SelectModelFeaturesForTrainingSVMFFSel(const VariableSizeMatrixType &scaledFeatureSet, const VectorDouble &AllLabels, const std::string outputdirectory)
{
  std::vector<int> SelectedFeatures;
  std::vector<double> CrossValidatedBalancedAccuraciesFinal;
  std::vector<int> UnselectedFeatures = UpdateUnselectedFeatures(SelectedFeatures, scaledFeatureSet.Cols());

  //feature selection mechanism
  //---------------------------
  while (UnselectedFeatures.size() > 0)
  {
    std::vector<double> CrossValidatedBalancedAccuraciesAfterAddingUnselectedFeatures;
    for (unsigned int featureNo = 0; featureNo < UnselectedFeatures.size(); featureNo++)
    {
      VariableSizeMatrixType reducedFeatureSet;
      reducedFeatureSet.SetSize(scaledFeatureSet.Rows(), SelectedFeatures.size() + 1);

      //copy the already selected features to the current feature set
      for (unsigned int j = 0; j < reducedFeatureSet.Rows(); j++)
      {
        for (unsigned int k = 0; k < SelectedFeatures.size(); k++)
        {
          double value = (int)(scaledFeatureSet(j, SelectedFeatures[k]) * 10000 + .5);
          reducedFeatureSet(j, k) = (double)value / 10000;
        }
      }

      //copy the new feature to the current feature set
      for (unsigned int j = 0; j < reducedFeatureSet.Rows(); j++)
      {
        double value = (int)(scaledFeatureSet(j, UnselectedFeatures[featureNo]) * 10000 + .5);
        reducedFeatureSet(j, reducedFeatureSet.Cols() - 1) = (double)value / 10000;
      }


      //check crossvalidated performance after adding the current feature
      double bestCV = 0;
      double bestC = 1;

      for (double cValue = -5; cValue <= 5; cValue = cValue + 1)
      {
        VectorDouble result = InternalCrossValidationResubstitution(reducedFeatureSet, AllLabels, pow(2, cValue), 0.01, 1); //classifier type 1 for linear kernel
        if (result[3] > bestCV)
        {
          bestC = pow(2, cValue);
          bestCV = result[3];
        }
      }

      CrossValidatedBalancedAccuraciesAfterAddingUnselectedFeatures.push_back(bestCV);
    }
    int index = std::distance(CrossValidatedBalancedAccuraciesAfterAddingUnselectedFeatures.begin(), std::max_element(CrossValidatedBalancedAccuraciesAfterAddingUnselectedFeatures.begin(), CrossValidatedBalancedAccuraciesAfterAddingUnselectedFeatures.end()));
    CrossValidatedBalancedAccuraciesFinal.push_back(CrossValidatedBalancedAccuraciesAfterAddingUnselectedFeatures[index]);
    SelectedFeatures.push_back(UnselectedFeatures[index]);
    std::cout << "CurrentSize=" << SelectedFeatures.size() << " Performance=" << CrossValidatedBalancedAccuraciesAfterAddingUnselectedFeatures[index] << std::endl;
    UnselectedFeatures = UpdateUnselectedFeatures(SelectedFeatures, scaledFeatureSet.Cols());
  }

  std::cout << "Feature Selection Done!!!" << std::endl;
  VariableLengthVectorType MovingAverageOnCrossValidatedPerformance;
  MovingAverageOnCrossValidatedPerformance.SetSize(CrossValidatedBalancedAccuraciesFinal.size());
  for (int index = 0; index < MovingAverageOnCrossValidatedPerformance.Size(); index++)
    MovingAverageOnCrossValidatedPerformance[index] = 0;

  for (unsigned int index = 1; index < CrossValidatedBalancedAccuraciesFinal.size() - 1; index++)
    MovingAverageOnCrossValidatedPerformance[index] = (CrossValidatedBalancedAccuraciesFinal[index] + CrossValidatedBalancedAccuraciesFinal[index - 1] + CrossValidatedBalancedAccuraciesFinal[index + 1]) / 3;

  int max_performance_counter = std::distance(CrossValidatedBalancedAccuraciesFinal.begin(), std::max_element(CrossValidatedBalancedAccuraciesFinal.begin(), CrossValidatedBalancedAccuraciesFinal.end()));
  std::vector<int> FinalSelectedFeatures;
  for (int index = 0; index <= max_performance_counter; index++)
    FinalSelectedFeatures.push_back(SelectedFeatures[index]);

  //Write the selected features
  std::ofstream myfile;
  myfile.open(outputdirectory + "/EGFRvIII_SelectedFeatures.csv");
  for (int counter = 0; counter < FinalSelectedFeatures.size(); counter++)
    myfile << std::to_string(FinalSelectedFeatures[counter]) + "\n";
  myfile.close();

  //optimize classifiers parameters on the selected feature set
  VariableSizeMatrixType FinalSelectedFeatureSet;
  FinalSelectedFeatureSet.SetSize(AllLabels.size(), FinalSelectedFeatures.size());

  //copy the already selected features to the current feature set
  for (unsigned int j = 0; j < FinalSelectedFeatureSet.Rows(); j++)
    for (unsigned int k = 0; k < FinalSelectedFeatures.size(); k++)
      FinalSelectedFeatureSet(j, k) = scaledFeatureSet(j, FinalSelectedFeatures[k]);

  double bestCV = 0;
  double bestC = 1;
  for (double cValue = -5; cValue <= 5; cValue = cValue + 1)
  {
    VectorDouble result = InternalCrossValidationResubstitution(FinalSelectedFeatureSet, AllLabels, pow(2, cValue), 0.01, 1);
    if (result[3] > bestCV)
    {
      bestC = pow(2, cValue);
      bestCV = result[3];
    }
  }
  std::cout << "Optimal C=" << bestC << std::endl;

  //model training and testing
  cv::Mat trainingData = cv::Mat::zeros(FinalSelectedFeatureSet.Rows(), FinalSelectedFeatureSet.Cols(), CV_32FC1);
  cv::Mat trainingLabels = cv::Mat::zeros(FinalSelectedFeatureSet.Rows(), 1, CV_32FC1);

  for (int copyDataCounter = 0; copyDataCounter < trainingData.rows; copyDataCounter++)
  {
    trainingLabels.ptr< float >(copyDataCounter)[0] = AllLabels[copyDataCounter];
    for (int copyDataCounter2 = 0; copyDataCounter2 < trainingData.cols; copyDataCounter2++)
      trainingData.ptr< float >(copyDataCounter)[copyDataCounter2] = FinalSelectedFeatureSet[copyDataCounter][copyDataCounter2];
  }

  trainingLabels.convertTo(trainingLabels, CV_32SC1);

  //make an SVM model
  auto svm = cv::ml::SVM::create();
  svm->setType(cv::ml::SVM::C_SVC);
  svm->setC(bestC);
  svm->setKernel(cv::ml::SVM::LINEAR);

  bool res;
  std::string msg;
  VectorDouble predictedLabels;
  VectorDouble predictedDistances;
  try
  {
    res = svm->train(trainingData, cv::ml::ROW_SAMPLE, trainingLabels);
    svm->save(outputdirectory + "/EGFRvIII_SVM_Model.xml");
  }
  catch (cv::Exception ex)
  {
    msg = ex.what();
  }
  return true;
}

bool EGFRvIIIIndexPredictor::SelectModelFeaturesForTrainingFromStudy(const VariableSizeMatrixType &scaledFeatureSet, const VectorDouble &AllLabels, const std::string outputdirectory)
{
  std::vector<int> FinalSelectedFeatures;
  FinalSelectedFeatures.push_back(1);
  FinalSelectedFeatures.push_back(2);
  FinalSelectedFeatures.push_back(3);
  FinalSelectedFeatures.push_back(4);
  FinalSelectedFeatures.push_back(54);
  FinalSelectedFeatures.push_back(420);
  FinalSelectedFeatures.push_back(421);

  //Write the selected features
  std::ofstream myfile;
  myfile.open(outputdirectory + "/EGFRvIII_SelectedFeatures.csv");
  for (int counter = 0; counter < FinalSelectedFeatures.size(); counter++)
    myfile << std::to_string(FinalSelectedFeatures[counter]) + "\n";
  myfile.close();

  //optimize classifiers parameters on the selected feature set
  VariableSizeMatrixType FinalSelectedFeatureSet;
  FinalSelectedFeatureSet.SetSize(AllLabels.size(), FinalSelectedFeatures.size());

  //copy the already selected features to the current feature set
  for (unsigned int j = 0; j < FinalSelectedFeatureSet.Rows(); j++)
    for (unsigned int k = 0; k < FinalSelectedFeatures.size(); k++)
      FinalSelectedFeatureSet(j, k) = scaledFeatureSet(j, FinalSelectedFeatures[k]);

  double bestCV = 0;
  double bestC = 1;
  for (double cValue = -5; cValue <= 5; cValue = cValue + 1)
  {
    VectorDouble result = InternalCrossValidationResubstitution(FinalSelectedFeatureSet, AllLabels, pow(2, cValue), 0.01, 1);
    if (result[3] > bestCV)
    {
      bestC = pow(2, cValue);
      bestCV = result[3];
    }
  }
  std::cout << "Optimal C=" << bestC << std::endl;

  //model training and testing
  cv::Mat trainingData = cv::Mat::zeros(FinalSelectedFeatureSet.Rows(), FinalSelectedFeatureSet.Cols(), CV_32FC1);
  cv::Mat trainingLabels = cv::Mat::zeros(FinalSelectedFeatureSet.Rows(), 1, CV_32FC1);

  for (int copyDataCounter = 0; copyDataCounter < trainingData.rows; copyDataCounter++)
  {
    trainingLabels.ptr< float >(copyDataCounter)[0] = AllLabels[copyDataCounter];
    for (int copyDataCounter2 = 0; copyDataCounter2 < trainingData.cols; copyDataCounter2++)
      trainingData.ptr< float >(copyDataCounter)[copyDataCounter2] = FinalSelectedFeatureSet[copyDataCounter][copyDataCounter2];
  }

  trainingLabels.convertTo(trainingLabels, CV_32SC1);

  //make an SVM model
  auto svm = cv::ml::SVM::create();
  svm->setType(cv::ml::SVM::C_SVC);
  svm->setC(bestC);
  svm->setKernel(cv::ml::SVM::LINEAR);

  bool res;
  std::string msg;
  VectorDouble predictedLabels;
  VectorDouble predictedDistances;
  try
  {
    res = svm->train(trainingData, cv::ml::ROW_SAMPLE, trainingLabels);
    svm->save(outputdirectory + "/EGFRvIII_SVM_Model.xml");
  }
  catch (cv::Exception ex)
  {
    msg = ex.what();
  }
//just for the testing purposes
  cv::Mat predicted(1, 1, CV_32F);
  for (int i = 0; i < trainingData.rows; i++)
  {
    cv::Mat sample = trainingData.row(i);
    svm->predict(sample, predicted, false);
    predictedLabels[i] = predicted.ptr< float >(0)[0];
  }
  for (int i = 0; i < trainingData.rows; i++)
  {
    cv::Mat sample = trainingData.row(i);
    svm->predict(sample, predicted, true);
    predictedDistances[i] = predicted.ptr< float >(0)[0];
  }
  if (predictedDistances[1] < 0 && predictedLabels[1]>0)
    std::cout << "signs flipped\n" << std::endl;



  return true;
}


//----------------------------------------------------------------------------------------------------------
std::vector<int> EGFRvIIIIndexPredictor::UpdateUnselectedFeatures(std::vector<int> SelectedFeatures, int featuresize)
{
  std::vector<int> UnselectedFeatures;
  for (unsigned int featureCounter = 0; featureCounter < featuresize; featureCounter++)
  {
    int found = 0;
    for (unsigned int index2 = 0; index2 < SelectedFeatures.size(); index2++)
    {
      if (SelectedFeatures[index2] == featureCounter)
      {
        found = 1;
        break;
      }
    }
    if (found == 0)
      UnselectedFeatures.push_back(featureCounter);
  }
  return UnselectedFeatures;
}
//-------------------------------------------------------------------------------------------------------------

VectorDouble EGFRvIIIIndexPredictor::InternalCrossValidationResubstitution(VariableSizeMatrixType inputFeatures, std::vector<double> inputLabels, double cValue, double gValue, int kerneltype)
{
  VariableLengthVectorType predictedLabels;
  predictedLabels.SetSize(inputLabels.size());
  VariableLengthVectorType predictedDistances;
  predictedDistances.SetSize(inputLabels.size());


  //copy training and test data from given dataset
  cv::Mat trainingData = cv::Mat::zeros(inputFeatures.Rows(), inputFeatures.Cols(), CV_32FC1), trainingLabels = cv::Mat::zeros(inputFeatures.Rows(), 1, CV_32FC1);
  for (int index = 0; index < trainingData.rows; index++)
  {
    trainingLabels.ptr< float >(index)[0] = inputLabels[index];
    for (int featureNo = 0; featureNo < trainingData.cols; featureNo++)
      trainingData.ptr< float >(index)[featureNo] = inputFeatures(index, featureNo);
  }

  trainingLabels.convertTo(trainingLabels, CV_32SC1);
  //make an SVM model
  auto svm = cv::ml::SVM::create();
  svm->setType(cv::ml::SVM::C_SVC);
  //svm->setTermCriteria(cv::TermCriteria(cv::TermCriteria::MAX_ITER, 100, 1e-6));
  svm->setC(cValue);

  if (kerneltype == 2)
  {
    svm->setGamma(gValue);
    svm->setKernel(cv::ml::SVM::RBF);
  }
  else
    svm->setKernel(cv::ml::SVM::LINEAR);
  bool res = true;
  std::string msg;

  try
  {
    res = svm->train(trainingData, cv::ml::ROW_SAMPLE, trainingLabels);
  }
  catch (cv::Exception ex)
  {
    msg = ex.what();
  }
  //apply SVM model on test data
  cv::Mat predicted(1, 1, CV_32F);
  for (int i = 0; i < trainingData.rows; i++)
  {
    cv::Mat sample = trainingData.row(i);
    svm->predict(sample, predicted, false);
    predictedLabels[i] = predicted.ptr< float >(0)[0];
  }
  for (int i = 0; i < trainingData.rows; i++)
  {
    cv::Mat sample = trainingData.row(i);
    svm->predict(sample, predicted, true);
    predictedDistances[i] = predicted.ptr< float >(0)[0];
  }
  if (predictedDistances[1] < 0 && predictedLabels[1]>0)
    std::cout << "signs flipped\n" << std::endl;

  VectorDouble results = CalculatePerformanceMeasures(predictedLabels, inputLabels);

  return results;
}
//-----------------------------

VectorDouble EGFRvIIIIndexPredictor::CalculatePerformanceMeasures(VariableLengthVectorType predictedLabels, VectorDouble GivenLabels)
{
  //calcualte performance measures
  double TP = 0;
  double TN = 0;
  double FP = 0;
  double FN = 0;
  VectorDouble result;

  for (int index = 0; index< predictedLabels.Size(); index++)
  {
    if (predictedLabels[index] == 1 && GivenLabels[index] == 1)
      TP++;
    else if (predictedLabels[index] == -1 && GivenLabels[index] == -1)
      TN++;
    else if (predictedLabels[index] == 1 && GivenLabels[index] == -1)
      FP++;
    else if (predictedLabels[index] == -1 && GivenLabels[index] == 1)
      FN++;
    else
    {
    }
  }
  result.push_back((TP + TN) / predictedLabels.Size());
  double sen = TP / (TP + FN);
  double spe = TN / (TN + FP);
  result.push_back(sen);
  result.push_back(spe);

  result.push_back((sen + spe) / 2);
  return result;
}

