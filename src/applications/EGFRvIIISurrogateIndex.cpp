#include "itkListSample.h"
#include "EGFRvIIISurrogateIndex.h"
#include "vnl_cholesky.h"
#include "itkCSVNumericObjectFileWriter.h"
#include "itkCovarianceSampleFilter.h"
#include "itkMeanSampleFilter.h"
#include "cbicaLogging.h"

EGFRStatusPredictor::EGFRStatusPredictor()
{
}
EGFRStatusPredictor::~EGFRStatusPredictor()
{
}
void EGFRStatusPredictor::CalculateAveragePerfusionSignal(const VectorVectorDouble &PerfusionIntensities, VectorDouble &avgSignal)
{
  for (unsigned int index = 0; index < PerfusionIntensities[0].size(); index++)
  {
    double temp = 0.0;
    for (unsigned int samples = 0; samples < PerfusionIntensities.size(); samples++)
      temp = temp + PerfusionIntensities[samples][index];
    avgSignal.push_back(temp / PerfusionIntensities.size());
  }
}

double EGFRStatusPredictor::CalculateMaximumDrop(const VectorDouble &avgSignal)
{
  double temp = 0.0;
  double mean = 0.0;
  for (int index = 0; index < 10; index++)
    temp = temp + avgSignal[index];
  mean = temp / 10;

  double min = avgSignal[0];
  for (unsigned int index = 1; index < avgSignal.size(); index++)
    if (avgSignal[index] < min)
      min = avgSignal[index];

  return (mean - min);
}

void EGFRStatusPredictor::CalculateQualifiedIndices(VectorVectorDouble &rpNearIntensities, VectorVectorDouble &rpFarIntensities, double & percentageNear, double & percentageFar)
{
  VectorVectorDouble revisedNearIntensities;
  VectorVectorDouble revisedFarIntensities;

  double originalNear = rpNearIntensities.size();
  double originalFar = rpFarIntensities.size();

  for (unsigned int index = 0; index < rpNearIntensities.size(); index++)
  {
    double mNearMean = 0.0;
    double mNearStd = 0.0;
    double temp = 0.0;
    for (unsigned int featureNo = 0; featureNo < rpNearIntensities[index].size(); featureNo++)
      temp = temp + rpNearIntensities[index][featureNo];
    mNearMean = temp / rpNearIntensities[index].size();

    temp = 0.0;
    for (unsigned int featureNo = 0; featureNo < rpNearIntensities[index].size(); featureNo++)
      temp = temp + (rpNearIntensities[index][featureNo] - mNearMean)*(rpNearIntensities[index][featureNo] - mNearMean);
    mNearStd = std::sqrt(temp / (rpNearIntensities[index].size() - 1));
    if (mNearStd != 0)
      revisedNearIntensities.push_back(rpNearIntensities[index]);
  }
  for (unsigned int index = 0; index < rpFarIntensities.size(); index++)
  {
    double mFarMean = 0.0;
    double mFarStd = 0.0;
    double temp = 0.0;
    for (unsigned int featureNo = 0; featureNo < rpFarIntensities[index].size(); featureNo++)
      temp = temp + rpFarIntensities[index][featureNo];
    mFarMean = temp / rpFarIntensities[index].size();

    temp = 0.0;
    for (unsigned int featureNo = 0; featureNo < rpFarIntensities[index].size(); featureNo++)
      temp = temp + (rpFarIntensities[index][featureNo] - mFarMean)*(rpFarIntensities[index][featureNo] - mFarMean);
    mFarStd = std::sqrt(temp / (rpFarIntensities[index].size() - 1));
    if (mFarStd != 0)
      revisedFarIntensities.push_back(rpFarIntensities[index]);
  }
  rpNearIntensities.clear();
  rpFarIntensities.clear();
  rpNearIntensities = revisedNearIntensities;
  rpFarIntensities = revisedFarIntensities;

  percentageNear = (revisedNearIntensities.size() * 100) / originalNear;
  percentageFar = (revisedFarIntensities.size() * 100) / originalFar;
}

VariableSizeMatrixType EGFRStatusPredictor::GetSumOfTwoMatrice(const VariableSizeMatrixType &matrix1, const  VariableSizeMatrixType &matrix2)
{
  VariableSizeMatrixType outputMatrix;
  outputMatrix.SetSize(matrix1.Rows(), matrix1.Cols());
  for (unsigned int i = 0; i < matrix1.Rows(); i++)
    for (unsigned int j = 0; j < matrix1.Cols(); j++)
      outputMatrix(i, j) = (matrix1(i, j) + matrix2(i, j)) / 2;
  return outputMatrix;
}

VariableSizeMatrixType EGFRStatusPredictor::MatrixTranspose(const VariableSizeMatrixType &inputmatrix)
{
  VariableSizeMatrixType output;
  output.SetSize(inputmatrix.Cols(), inputmatrix.Rows());

  for (unsigned int i = 0; i < output.Rows(); i++)
    for (unsigned int j = 0; j < output.Cols(); j++)
      output(i, j) = inputmatrix(j, i);
  return output;
}

VariableSizeMatrixType EGFRStatusPredictor::GetCovarianceMatrix(const VectorVectorDouble &inputData)
{
  int NumberOfSamples = inputData.size();
  const unsigned int MeasurementVectorLength = EGFR_PCS;
  typedef itk::Vector< double, MeasurementVectorLength > MeasurementVectorType;
  typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;
  SampleType::Pointer sample = SampleType::New();
  sample->SetMeasurementVectorSize(MeasurementVectorLength);

  for (int i = 0; i < NumberOfSamples; i++)
  {
    MeasurementVectorType mv;
    mv[0] = inputData[i][0];
    mv[1] = inputData[i][1];
    mv[2] = inputData[i][2];
    sample->PushBack(mv);
  }
  typedef itk::Statistics::CovarianceSampleFilter< SampleType >CovarianceAlgorithmType;
  CovarianceAlgorithmType::Pointer covarianceAlgorithm = CovarianceAlgorithmType::New();
  covarianceAlgorithm->SetInput(sample);
  covarianceAlgorithm->Update();
  VariableSizeMatrixType CovMatrix = covarianceAlgorithm->GetCovarianceMatrix();
  return CovMatrix;
}

VariableSizeMatrixType EGFRStatusPredictor::GetCholeskyFactorization(VariableSizeMatrixType &inputData)
{
  inputData[0][1] = 0;
  inputData[0][2] = 0;
  inputData[1][2] = 0;

  inputData[0][0] = std::sqrt(inputData[0][0]);
  inputData[1][0] = inputData[1][0] / inputData[0][0];
  inputData[2][0] = inputData[2][0] / inputData[0][0];

  inputData[1][1] = inputData[1][1] - inputData[1][0] * inputData[1][0];
  inputData[2][1] = inputData[1][2] - inputData[1][0] * inputData[2][0];
  inputData[2][2] = inputData[2][2] - inputData[2][0] * inputData[2][0];

  inputData[1][1] = std::sqrt(inputData[1][1]);
  inputData[2][1] = inputData[2][1] / inputData[1][1];

  inputData[2][2] = std::sqrt(inputData[2][2] - inputData[2][1]);

  return inputData;
}

double EGFRStatusPredictor::CalculateBhattacharyaCoefficient(const VectorVectorDouble &rpNearIntensities, const VectorVectorDouble &rpFarIntensities)
{
  double bDistance = 0;
  try
  {
    VariableLengthVectorType nMeanVector;
    VariableLengthVectorType fMeanVector;
    nMeanVector.SetSize(EGFR_PCS);
    fMeanVector.SetSize(EGFR_PCS);

    for (int featureNo = 0; featureNo < EGFR_PCS; featureNo++)
    {
      double temp = 0.0;
      for (unsigned int sampleNo = 0; sampleNo < rpNearIntensities.size(); sampleNo++)
        temp = temp + rpNearIntensities[sampleNo][featureNo];
      nMeanVector[featureNo] = temp / rpNearIntensities.size();

      temp = 0.0;
      for (unsigned int sampleNo = 0; sampleNo < rpFarIntensities.size(); sampleNo++)
        temp = temp + rpFarIntensities[sampleNo][featureNo];
      fMeanVector[featureNo] = temp / rpFarIntensities.size();
    }
    VariableSizeMatrixType nCovMatrix = GetCovarianceMatrix(rpNearIntensities);
    VariableSizeMatrixType fCovMatrix = GetCovarianceMatrix(rpFarIntensities);

    VariableSizeMatrixType sumCovariance = GetSumOfTwoMatrice(nCovMatrix, fCovMatrix);

    VariableSizeMatrixType meanDifference;
    meanDifference.SetSize(1, nMeanVector.Size());
    for (unsigned int i = 0; i < nMeanVector.Size(); i++)
      meanDifference[0][i] = nMeanVector[i] - fMeanVector[i];

    //-------------vnl matrices------------------------
    //this segment of code is applying cholesky factorization.
    //EGFR_PCS is the number of principal components that code uses to calculate covariance matrices and bhattacharyya distance afterwards
    vnl_matrix<double> covarianceMatrix;
    covarianceMatrix.set_size(EGFR_PCS, EGFR_PCS);
    for (int i = 0; i < EGFR_PCS; i++)
      for (int j = 0; j < EGFR_PCS; j++)
        covarianceMatrix(i, j) = sumCovariance(i, j);
    vnl_matrix<double> covMatrix1;
    vnl_matrix<double> covMatrix2;
    covMatrix1.set_size(EGFR_PCS, EGFR_PCS);
    covMatrix2.set_size(EGFR_PCS, EGFR_PCS);
    for (int i = 0; i < EGFR_PCS; i++)
      for (int j = 0; j < EGFR_PCS; j++)
      {
        covMatrix1(i, j) = nCovMatrix(i, j);
        covMatrix2(i, j) = fCovMatrix(i, j);
      }
    vnl_cholesky covObj(covarianceMatrix);
    vnl_matrix<double> inverse = covObj.inverse();
    vnl_cholesky covObj1(covMatrix1);
    vnl_cholesky covObj2(covMatrix2);

    if (inverse.rows() == 0)
      bDistance = 0;
    else
    {
      VectorDouble tempmatrix;
      tempmatrix.push_back(meanDifference[0][0] * inverse[0][0] + meanDifference[0][1] * inverse[1][0] + meanDifference[0][2] * inverse[2][0]);
      tempmatrix.push_back(meanDifference[0][0] * inverse[0][1] + meanDifference[0][1] * inverse[1][1] + meanDifference[0][2] * inverse[2][1]);
      tempmatrix.push_back(meanDifference[0][0] * inverse[0][2] + meanDifference[0][1] * inverse[1][2] + meanDifference[0][2] * inverse[2][2]);
      double firstfactor = (tempmatrix[0] * meanDifference[0][0] + tempmatrix[1] * meanDifference[0][1] + tempmatrix[2] * meanDifference[0][2]) / 8;
      double secondfactor = std::log(covObj.determinant() / std::sqrt(covObj1.determinant()*covObj2.determinant())) / 2;
      bDistance = firstfactor + secondfactor;
    }
  }
  catch (itk::ExceptionObject & excp)
  {
    bDistance = 0;
    cbica::Logging(loggerFile, "Exception detected while trying to Calculate Bhattacharya Coefficient: '" + std::string(excp.GetDescription()));
    exit(EXIT_FAILURE);
  }
  return bDistance;
}

