#include "EGFRvIIIIndexPredictor.h"
#include "itkCSVNumericObjectFileWriter.h"
#include "CaPTkGUIUtils.h"
#include "CaPTkClassifierUtils.h"
#include "CaPTkUtils.h"
#include "cbicaUtilities.h"
#include "cbicaITKSafeImageIO.h"
typedef itk::Image< float, 3 > ImageType;

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
      AtlasSegmentationsAllPatients.push_back(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_ATLAS])));
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

  //ImageType::Pointer NEGAtlasImagePointer = cbica::ReadImage<ImageType>(getCaPTkDataDir() + "/egfrv3/EGFRneg.nii.gz");
  //ImageType::Pointer POSAtlasImagePointer = cbica::ReadImage<ImageType>(getCaPTkDataDir() + "/egfrv3/EGFRpos.nii.gz");
    
  FeaturesOfAllSubjects.SetSize(qualifiedsubjects.size(), EGFR_NO_OF_FEATURES);
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
      ImageType::Pointer TemplateImagePointer9Regions = cbica::ReadImage<ImageType>(getCaPTkDataDir() + "/egfrv3/template9regions.nii.gz");
      ImageType::Pointer TemplateImagePointer21Regions = cbica::ReadImage<ImageType>(getCaPTkDataDir() + "/egfrv3/template21regions.nii.gz");
      ImageType::Pointer TemplateImagePointerAllRegions = cbica::ReadImage<ImageType>(getCaPTkDataDir() + "/egfrv3/templateallregions.nii.gz");

      ImageType::Pointer RCBVImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_RCBV])));
      ImageType::Pointer PHImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PH])));
      ImageType::Pointer PSRImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PSR])));

      ImageType::Pointer T1CEImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1CE])));
      ImageType::Pointer T2FlairImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR])));
      ImageType::Pointer T1ImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1])));
      ImageType::Pointer T2ImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2])));

      ImageType::Pointer AXImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_AX])));
      ImageType::Pointer RADImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_RAD])));
      ImageType::Pointer FAImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_FA])));
      ImageType::Pointer TRImagePointer = RescaleImageIntensity<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_TR])));

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

  //WriteCSVFilesWithHorizontalAndVerticalHeaders(FeaturesOfAllSubjects, patient_ids, FeatureLabels, outputdirectory + "/RawFeatures.csv");

  VariableSizeMatrixType scaledFeatureSet;
  scaledFeatureSet.SetSize(FeaturesOfAllSubjects.Rows(), FeaturesOfAllSubjects.Cols());
  VariableLengthVectorType meanVector;
  VariableLengthVectorType stdVector;
  mFeatureScalingLocalPtr.ScaleGivenTrainingFeatures(FeaturesOfAllSubjects, scaledFeatureSet, meanVector, stdVector);
  
  try
  {
    WriteCSVFiles(meanVector, outputdirectory + "/EGFRvIII_ZScore_Mean.csv");
    WriteCSVFiles(stdVector, outputdirectory + "/EGFRvIII_ZScore_Std.csv");
  }
  catch (const std::exception& e1)
  {
    std::cout << "scaling writing error" << std::endl;
    logger.WriteError("Error in writing output files to the output directory = " + outputdirectory + "Error code : " + std::string(e1.what()));
    return false;
  }

  std::cout << std::endl << "Scaling done....." << std::endl;
  std::cout << "Feature writing started:" << std::endl;
  //WriteCSVFilesWithHorizontalAndVerticalHeaders(scaledFeatureSet, patient_ids, FeatureLabels, outputdirectory + "/ScaledFeatures.csv");

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

VectorDouble EGFRvIIIIndexPredictor::EGFRvIIIPredictionOnExistingModel(const std::string &modeldirectory, const std::string &inputdirectory, const std::vector < std::map < CAPTK::ImageModalityType, std::string>> &qualifiedsubjects, const std::string &outputdirectory)
{
  std::cout << "Started reading model parameters" << std::endl;
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
    logger.WriteError("Cannot find the file 'EGFRvIII_HMFeatures_Configuration.csv' in the ../data/egfr directory. Error code : " + std::string(e1.what()));
    return results;
  }

  MatrixType meanMatrix;
  VariableLengthVectorType mean;
  VariableLengthVectorType stddevition;
  VariableLengthVectorType selectedfeatures;

  //read z-score mean values from the model directory
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
  //read z-score std. deviation values from the model directory
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

  //read selected feature indices from the model directory
  MatrixType featuresMatrix;
  std::cout << modeldirectory + "/EGFRvIII_SelectedFeatures.csv" << std::endl;
  try
  {
    reader->SetFileName(modeldirectory + "/EGFRvIII_SelectedFeatures.csv");
    reader->SetFieldDelimiterCharacter(',');
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();
    reader->Parse();
    featuresMatrix = reader->GetArray2DDataObject()->GetMatrix();

    selectedfeatures.SetSize(featuresMatrix.size());
    for (unsigned int i = 0; i < featuresMatrix.size(); i++)
      selectedfeatures[i] = featuresMatrix(0, i);
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Error in reading the file: " + modeldirectory + "/EGFRvIII_SelectedFeatures.csv. Error code : " + std::string(e1.what()));
    return results;
  }
  //read EGFRneg and EGFRpos atlases from the model directory
  ImageType::Pointer POSAtlasImagePointer;
  ImageType::Pointer NEGAtlasImagePointer;
  try
  {
    NEGAtlasImagePointer = cbica::ReadImage<ImageType>(modeldirectory+ "/EGFRneg.nii.gz");
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Error in reading atlas files: " + modeldirectory + "/EGFRneg.csv. Error code : " + std::string(e1.what()));
    return results;
  }

  try
  {
    POSAtlasImagePointer = cbica::ReadImage<ImageType>(modeldirectory + "/EGFRpos.nii.gz");
  }
  catch (const std::exception& e1)
  {
    logger.WriteError("Error in reading atlas files: " + modeldirectory + "/EGFRpos.csv. Error code : " + std::string(e1.what()));
    return results;
  }
  std::cout << "Finished readig model parameters" << std::endl;

  //extract features from the given dataset
  VariableSizeMatrixType FeaturesOfAllSubjects;
  FeaturesOfAllSubjects.SetSize(qualifiedsubjects.size(), EGFR_NO_OF_FEATURES);
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
      ImageType::Pointer TemplateImagePointer9Regions = cbica::ReadImage<ImageType>(getCaPTkDataDir() + "egfrv3/template9regions.nii.gz");
      ImageType::Pointer TemplateImagePointer21Regions = cbica::ReadImage<ImageType>(getCaPTkDataDir() + "egfrv3/template21regions.nii.gz");
      ImageType::Pointer TemplateImagePointerAllRegions = cbica::ReadImage<ImageType>(getCaPTkDataDir() + "egfrv3/templateallregions.nii.gz");
      ImageType::Pointer T1CEImagePointer = RescaleImageIntensity<ImageType>(CapImageIntensityWithPercentile<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1CE]))));
      ImageType::Pointer RCBVImagePointer = RescaleImageIntensity<ImageType>(CapImageIntensityWithPercentile<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_RCBV]))));
      ImageType::Pointer PHImagePointer = RescaleImageIntensity<ImageType>(CapImageIntensityWithPercentile<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PH]))));
      ImageType::Pointer T2FlairImagePointer = RescaleImageIntensity<ImageType>(CapImageIntensityWithPercentile<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR]))));
      ImageType::Pointer T1ImagePointer = RescaleImageIntensity<ImageType>(CapImageIntensityWithPercentile<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1]))));
      ImageType::Pointer T2ImagePointer = RescaleImageIntensity<ImageType>(CapImageIntensityWithPercentile<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2]))));
      ImageType::Pointer AXImagePointer = RescaleImageIntensity<ImageType>(CapImageIntensityWithPercentile<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_AX]))));
      ImageType::Pointer RADImagePointer = RescaleImageIntensity<ImageType>(CapImageIntensityWithPercentile<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_RAD]))));
      ImageType::Pointer FAImagePointer = RescaleImageIntensity<ImageType>(CapImageIntensityWithPercentile<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_FA]))));
      ImageType::Pointer TRImagePointer = RescaleImageIntensity<ImageType>(CapImageIntensityWithPercentile<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_TR]))));
      ImageType::Pointer PSRImagePointer = RescaleImageIntensity<ImageType>(CapImageIntensityWithPercentile<ImageType>(cbica::ReadImage<ImageType>(static_cast<std::string>(currentsubject[CAPTK::ImageModalityType::IMAGE_TYPE_PSR]))));

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
  std::string FeatureLabels[EGFR_NO_OF_FEATURES] = { "Age","Loc_Mean_diff","Loc_Max_Diff","Loc_Mean_Ratio","Loc_Max_Ratio","ET_PatientSpace","NCR_PatientSpace","ED_PatientSpace","brain size_PatientSpace","ET+NCR_PatientSpace","ET+NCR/BS_PatientSpace","Edema/BS","ET/ET+NCR_PatientSpace","NCR/ET+NCR_PatientSpace","Edema/ET+NCR_PatientSpace",
    "ET_Atl","NCR_Atl","ED_Atl","brain size_Atl","ET+NCR_Atl","ET+NCR/BS_Atl","Edema/BS_Atl","ET/ET+NCR_Atl","NCR/ET+NCR_Atl","Edema/ET+NCR_Atl",
    "Vent_Tumor_Dist","Vent_ED_Dist","Mean_ET_T1CE","STD_ET_T1CE","Mean_ET_T1","STD_ET_T1","Mean_ET_T2","STD_ET_T2","Mean_ET_Flair","STD_ET_Flair","Mean_ET_PH","STD_ET_PH",
    "Mean_ET_PSR","STD_ET_PSR","Mean_ET_RCBV","STD_ET_RCBV","Mean_ET_FA","STD_ET_FA","Mean_ET_AX","STD_ET_AX","Mean_ET_RAD","STD_ET_RAD","Mean_ET_TR","STD_ET_TR",
    "Mean_NCR_T1CE","STD_NCR_T1CE","Mean_NCR_T1","STD_NCR_T1","Mean_NCR_T2","STD_NCR_T2","Mean_NCR_Flair","STD_NCR_Flair","Mean_NCR_PH","STD_NCR_PH",
    "Mean_NCR_PSR","STD_NCR_PSR","Mean_NCR_RCBV","STD_NCR_RCBV","Mean_NCR_FA","STD_NCR_FA","Mean_NCR_AX","STD_NCR_AX","Mean_NCR_RAD","STD_NCR_RAD","Mean_NCR_TR","STD_NCR_TR",
    "Mean_ED_T1CE","STD_ED_T1CE","Mean_ED_T1","STD_ED_T1","Mean_ED_T2","STD_ED_T2","Mean_ED_Flair","STD_ED_Flair","Mean_ED_PH","STD_ED_PH",
    "Mean_ED_PSR","STD_ED_PSR","Mean_ED_RCBV","STD_ED_RCBV","Mean_ED_FA","STD_ED_FA","Mean_ED_AX","STD_ED_AX","Mean_ED_RAD","STD_ED_RAD","Mean_ED_TR","STD_ED_TR",
    "ET_T1CE_Bin1","ET_T1CE_Bin2","ET_T1CE_Bin3","ET_T1CE_Bin4","ET_T1CE_Bin5","ET_T1CE_Bin6","ET_T1CE_Bin7","ET_T1CE_Bin8","ET_T1CE_Bin9","ET_T1CE_Bin10",
    "ED_T1CE_Bin1","ED_T1CE_Bin2","ED_T1CE_Bin3","ED_T1CE_Bin4","ED_T1CE_Bin5","ED_T1CE_Bin6","ED_T1CE_Bin7","ED_T1CE_Bin8","ED_T1CE_Bin9","ED_T1CE_Bin10",
    "NCR_T1CE_Bin1","NCR_T1CE_Bin2","NCR_T1CE_Bin3","NCR_T1CE_Bin4","NCR_T1CE_Bin5","NCR_T1CE_Bin6","NCR_T1CE_Bin7",
    "ET_T1_Bin1","ET_T1_Bin2","ET_T1_Bin3","ET_T1_Bin4","ET_T1_Bin5","ET_T1_Bin6","ET_T1_Bin7","ET_T1_Bin8","ET_T1_Bin9","ET_T1_Bin10",
    "ED_T1_Bin1","ED_T1_Bin2","ED_T1_Bin3","ED_T1_Bin4","ED_T1_Bin5","ED_T1_Bin6","ED_T1_Bin7","ED_T1_Bin8","ED_T1_Bin9","ED_T1_Bin10",
    "NCR_T1_Bin1","NCR_T1_Bin2","NCR_T1_Bin3","NCR_T1_Bin4","NCR_T1_Bin5","NCR_T1_Bin6","NCR_T1_Bin7","NCR_T1_Bin8","NCR_T1_Bin9","NCR_T1_Bin10",
    "ET_T2_Bin1","ET_T2_Bin2","ET_T2_Bin3","ET_T2_Bin4","ET_T2_Bin5","ET_T2_Bin6","ET_T2_Bin7","ET_T2_Bin8","ET_T2_Bin9","ET_T2_Bin10",
    "ED_T2_Bin1","ED_T2_Bin2","ED_T2_Bin3","ED_T2_Bin4","ED_T2_Bin5","ED_T2_Bin6","ED_T2_Bin7","ED_T2_Bin8","ED_T2_Bin9","ED_T2_Bin10",
    "NCR_T2_Bin1","NCR_T2_Bin2","NCR_T2_Bin3","NCR_T2_Bin4","NCR_T2_Bin5","NCR_T2_Bin6","NCR_T2_Bin7","NCR_T2_Bin8","NCR_T2_Bin9","NCR_T2_Bin10",
    "ET_Flair_Bin1","ET_Flair_Bin2","ET_Flair_Bin3","ET_Flair_Bin4","ET_Flair_Bin5","ET_Flair_Bin6","ET_Flair_Bin7","ET_Flair_Bin8","ET_Flair_Bin9","ET_Flair_Bin10",
    "ED_Flair_Bin1","ED_Flair_Bin2","ED_Flair_Bin3","ED_Flair_Bin4","ED_Flair_Bin5","ED_Flair_Bin6","ED_Flair_Bin7","ED_Flair_Bin8","ED_Flair_Bin9","ED_Flair_Bin10",
    "NCR_Flair_Bin1","NCR_Flair_Bin2","NCR_Flair_Bin3","NCR_Flair_Bin4","NCR_Flair_Bin5","NCR_Flair_Bin6","NCR_Flair_Bin7","NCR_Flair_Bin8","NCR_Flair_Bin9","NCR_Flair_Bin10",
    "ET_PH_Bin1","ET_PH_Bin2","ET_PH_Bin3","ET_PH_Bin4","ET_PH_Bin5","ET_PH_Bin6","ET_PH_Bin7","ET_PH_Bin8","ET_PH_Bin9","ET_PH_Bin10",
    "ED_PH_Bin1","ED_PH_Bin2","ED_PH_Bin3","ED_PH_Bin4","ED_PH_Bin5","ED_PH_Bin6","ED_PH_Bin7","ED_PH_Bin8","ED_PH_Bin9","ED_PH_Bin10",
    "NCR_PH_Bin1","NCR_PH_Bin2","NCR_PH_Bin3","NCR_PH_Bin4","NCR_PH_Bin5","NCR_PH_Bin6","NCR_PH_Bin7","NCR_PH_Bin8","NCR_PH_Bin9","NCR_PH_Bin10",
    "ET_PSR_Bin1","ET_PSR_Bin2","ET_PSR_Bin3","ET_PSR_Bin4","ET_PSR_Bin5","ET_PSR_Bin6","ET_PSR_Bin7","ET_PSR_Bin8","ET_PSR_Bin9","ET_PSR_Bin10",
    "ED_PSR_Bin1","ED_PSR_Bin2","ED_PSR_Bin3","ED_PSR_Bin4","ED_PSR_Bin5","ED_PSR_Bin6","ED_PSR_Bin7","ED_PSR_Bin8","ED_PSR_Bin9","ED_PSR_Bin10",
    "NCR_PSR_Bin1","NCR_PSR_Bin2","NCR_PSR_Bin3","NCR_PSR_Bin4","NCR_PSR_Bin5","NCR_PSR_Bin6","NCR_PSR_Bin7","NCR_PSR_Bin8","NCR_PSR_Bin9","NCR_PSR_Bin10",
    "ET_RCBV_Bin1","ET_RCBV_Bin2","ET_RCBV_Bin3","ET_RCBV_Bin4","ET_RCBV_Bin5","ET_RCBV_Bin6","ET_RCBV_Bin7","ET_RCBV_Bin8","ET_RCBV_Bin9","ET_RCBV_Bin10",
    "ED_RCBV_Bin1","ED_RCBV_Bin2","ED_RCBV_Bin3","ED_RCBV_Bin4","ED_RCBV_Bin5","ED_RCBV_Bin6","ED_RCBV_Bin7","ED_RCBV_Bin8","ED_RCBV_Bin9","ED_RCBV_Bin10",
    "NCR_RCBV_Bin1","NCR_RCBV_Bin2","NCR_RCBV_Bin3","NCR_RCBV_Bin4","NCR_RCBV_Bin5","NCR_RCBV_Bin6","NCR_RCBV_Bin7","NCR_RCBV_Bin8","NCR_RCBV_Bin9","NCR_RCBV_Bin10",
    "ET_FA_Bin1","ET_FA_Bin2","ET_FA_Bin3","ET_FA_Bin4","ET_FA_Bin5","ET_FA_Bin6","ET_FA_Bin7","ET_FA_Bin8","ET_FA_Bin9","ET_FA_Bin10",
    "ED_FA_Bin1","ED_FA_Bin2","ED_FA_Bin3","ED_FA_Bin4","ED_FA_Bin5","ED_FA_Bin6","ED_FA_Bin7","ED_FA_Bin8","ED_FA_Bin9","ED_FA_Bin10",
    "NCR_FA_Bin1","   NCR_FA_Bin2","NCR_FA_Bin3","NCR_FA_Bin4","NCR_FA_Bin5","NCR_FA_Bin6","NCR_FA_Bin7","NCR_FA_Bin8","NCR_FA_Bin9","NCR_FA_Bin10",
    "ET_AX_Bin1","ET_AX_Bin2","ET_AX_Bin3","ET_AX_Bin4","ET_AX_Bin5","ET_AX_Bin6","ET_AX_Bin7","ET_AX_Bin8","ET_AX_Bin9","ET_AX_Bin10",
    "ED_AX_Bin1","ED_AX_Bin2","ED_AX_Bin3","ED_AX_Bin4","ED_AX_Bin5","ED_AX_Bin6","ED_AX_Bin7","ED_AX_Bin8","ED_AX_Bin9","ED_AX_Bin10",
    "NCR_AX_Bin1","NCR_AX_Bin2","NCR_AX_Bin3","NCR_AX_Bin4","NCR_AX_Bin5","NCR_AX_Bin6","NCR_AX_Bin7","NCR_AX_Bin8","NCR_AX_Bin9","NCR_AX_Bin10",
    "ET_RAD_Bin1","ET_RAD_Bin2","ET_RAD_Bin3","ET_RAD_Bin4","ET_RAD_Bin5","ET_RAD_Bin6","ET_RAD_Bin7","ET_RAD_Bin8","ET_RAD_Bin9","   ET_RAD_Bin10",
    "ED_RAD_Bin1","ED_RAD_Bin2","ED_RAD_Bin3","ED_RAD_Bin4","ED_RAD_Bin5","ED_RAD_Bin6","ED_RAD_Bin7","ED_RAD_Bin8","ED_RAD_Bin9","ED_RAD_Bin10",
    "NCR_RAD_Bin1","NCR_RAD_Bin2","NCR_RAD_Bin3","NCR_RAD_Bin4","NCR_RAD_Bin5","NCR_RAD_Bin6","NCR_RAD_Bin7","NCR_RAD_Bin8","NCR_RAD_Bin9","NCR_RAD_Bin10",
    "ET_TR_Bin1","ET_TR_Bin2","ET_TR_Bin3","ET_TR_Bin4","ET_TR_Bin5","ET_TR_Bin6","ET_TR_Bin7","ET_TR_Bin8","ET_TR_Bin9","ET_TR_Bin10",
    "ED_TR_Bin1","ED_TR_Bin2","ED_TR_Bin3","ED_TR_Bin4","ED_TR_Bin5","ED_TR_Bin6","ED_TR_Bin7","ED_TR_Bin8","ED_TR_Bin9","ED_TR_Bin10",
    "NCR_TR_Bin1","NCR_TR_Bin2","NCR_TR_Bin3","NCR_TR_Bin4","NCR_TR_Bin5","NCR_TR_Bin6","NCR_TR_Bin7","NCR_TR_Bin8","NCR_TR_Bin9","NCR_TR_Bin10",
    "Frontal","Temporal","Parietal","Basal G","Insula","CC_Fornix","Occipital","Cere","Brain stem",
    "ROI1","ROI2","ROI3","ROI4","ROI5","ROI6","ROI7","ROI8","ROI9","ROI10","ROI11","ROI12","ROI13","ROI14","ROI15","ROI16","ROI18","ROI19","ROI20" };

  std::cout << "Feature writing started:" << std::endl;
  //write raw extracted features to a .csv file
  std::vector<std::string> StringFeatureLabels;
  for (int index = 0; index < EGFR_NO_OF_FEATURES; index++)
    StringFeatureLabels.push_back(FeatureLabels[index]);

  WriteCSVFilesWithHorizontalAndVerticalHeaders(FeaturesOfAllSubjects, patient_ids, StringFeatureLabels, outputdirectory + "/RawFeatures.csv");

  //scale raw features based on the z-score values
  VariableSizeMatrixType ScaledTestingData = mFeatureScalingLocalPtr.ScaleGivenTestingFeatures(FeaturesOfAllSubjects, mean, stddevition);
  
  //write scaled features in a .csv file
  WriteCSVFilesWithHorizontalAndVerticalHeaders(ScaledTestingData, patient_ids, StringFeatureLabels, outputdirectory + "/ScaledFeatures.csv");

  //select model features
  VariableSizeMatrixType ModelSelectedFeatures = SelectModelFeatures(ScaledTestingData,selectedfeatures);
  std::vector<std::string> SelectedFeatureLabels;
  for (int index = 0; index < selectedfeatures.Size(); index++)
  {
    int currentindex = selectedfeatures[index];
    SelectedFeatureLabels.push_back(FeatureLabels[currentindex]);
  }

  WriteCSVFilesWithHorizontalAndVerticalHeaders(ModelSelectedFeatures, patient_ids, SelectedFeatureLabels, outputdirectory + "/ScaledSelectedFeatures.csv");

  //apply SVM model on the test data
  try
  {
    std::ofstream myfile;
    myfile.open(outputdirectory + "/results.csv");
    myfile << "SubjectName,Score, Result \n";
    //existing model in the .csv file
    if (cbica::fileExists(modeldirectory + "/EGFRvIII_SVM_Model.csv") == true)
    {
      VariableLengthVectorType result;
      result = DistanceFunctionLinear(ModelSelectedFeatures, modeldirectory + "/EGFRvIII_SVM_Model.csv");
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
      //new models trained by the users using CaPTk
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
  ModelSelectedFeatures.SetSize(ModelFeatures.Rows(), selectedFeatures.Size());
  int counter = 0;
  for (unsigned int i = 0; i < selectedFeatures.Size(); i++)
  {
    for (unsigned int j = 0; j < ModelFeatures.Rows(); j++)
      ModelSelectedFeatures(j, counter) = ModelFeatures(j, selectedFeatures[i]);
    counter++;
  }
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

