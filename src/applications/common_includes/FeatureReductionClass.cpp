/**
\file  FeatureReductionClass.cpp

\brief Implementation of the FeatureReductionClass

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

//#include <vtkVersion.h>
//#include "vtkDoubleArray.h"
//#include "vtkMultiBlockDataSet.h"
//#include "vtkPCAStatistics.h"
//#include "vtkStringArray.h"
//#include "vtkTable.h"
#include "FeatureReductionClass.h"
//#include "itkCSVNumericObjectFileWriter.h"
//#include "CAPTk.h"
#include "vtkSmartPointer.h"
#include "vtkVariant.h"
#include "vtkTable.h"
#include "vtkDoubleArray.h"
#include "vtkPCAStatistics.h"
#include "vtkStatisticsAlgorithm.h"
#include "CaPTkDefines.h"

FeatureReductionClass::FeatureReductionClass()
{
  PCATransformationMatrix.SetSize(NO_OF_PCA_FEATURES, NO_OF_PCA_FEATURES);
  mPMeanvector.SetSize(NO_OF_PCA_FEATURES);
}

FeatureReductionClass::~FeatureReductionClass()
{
  PCATransformationMatrix.SetSize(0, 0);
  mPMeanvector.SetSize(0);
}

VariableLengthVectorType FeatureReductionClass::ComputeMeanOfGivenFeatureVectors(vtkSmartPointer< vtkTable > &inputdata)
{
  vtkIdType NumberOfSamples = inputdata.GetPointer()->GetNumberOfRows();
  vtkIdType NumberOfFeatures = inputdata.GetPointer()->GetNumberOfColumns();

  VariableLengthVectorType mMeanVector;
  mMeanVector.SetSize(NumberOfFeatures);
  for (vtkIdType featureNo = 0; featureNo < NumberOfFeatures; featureNo++)
  {
    double temp = 0.0;
    for (vtkIdType sampleNo = 0; sampleNo < NumberOfSamples; sampleNo++)
      temp += inputdata.GetPointer()->GetValue(sampleNo, featureNo).ToDouble();

    mMeanVector[featureNo] = temp / NumberOfSamples;
  }
  return mMeanVector;
}

VariableLengthVectorType FeatureReductionClass::ComputeMeanOfGivenFeatureVectors(vnl_matrix<double> &inputdata)
{
  vtkIdType NumberOfSamples = inputdata.rows();
  vtkIdType NumberOfFeatures = inputdata.cols();

  VariableLengthVectorType mMeanVector;
  mMeanVector.SetSize(NumberOfFeatures);
  for (vtkIdType featureNo = 0; featureNo < NumberOfFeatures; featureNo++)
  {
    double temp = 0.0;
    for (vtkIdType sampleNo = 0; sampleNo < NumberOfSamples; sampleNo++)
      temp += inputdata(sampleNo, featureNo);

    mMeanVector[featureNo] = temp / NumberOfSamples;
  }
  return mMeanVector;
}

VariableLengthVectorType FeatureReductionClass::ComputeMeanOfGivenFeatureVectors(VectorVectorDouble &inputdata)
{
  vtkIdType NumberOfSamples = inputdata.size();
  vtkIdType NumberOfFeatures = inputdata[0].size();

  VariableLengthVectorType mMeanVector;
  mMeanVector.SetSize(NumberOfFeatures);
  for (vtkIdType featureNo = 0; featureNo < NumberOfFeatures; featureNo++)
  {
    double temp = 0.0;
    for (vtkIdType sampleNo = 0; sampleNo < NumberOfSamples; sampleNo++)
      temp += inputdata[sampleNo][featureNo];

    mMeanVector[featureNo] = temp / NumberOfSamples;
  }
  return mMeanVector;
}

VariableSizeMatrixType FeatureReductionClass::MatrixTranspose(VariableSizeMatrixType &inputmatrix)
{
  VariableSizeMatrixType output = VariableSizeMatrixType(inputmatrix.Cols(), inputmatrix.Rows());
  //vnl_matrix<double> test1 = inputmatrix.GetTranspose();// vnl_matrix cannot be directly converted into a VariableSizeMatrix 

  for (unsigned int i = 0; i < output.Rows(); i++)
  {
    for (unsigned int j = 0; j < output.Cols(); j++)
      output(i, j) = inputmatrix(j, i);
  }

  return output;
}

vtkSmartPointer< vtkTable >  FeatureReductionClass::GetDiscerningPerfusionTimePointsDynamic(VectorVectorDouble &intensities)
{
  mPMeanvector = ComputeMeanOfGivenFeatureVectors(intensities);
  size_t NumberOfFeatures = intensities[0].size();
  size_t NumberOfSamples = intensities.size();

  vtkSmartPointer<vtkTable> datasetTable = vtkSmartPointer<vtkTable>::New();
  std::vector<vtkSmartPointer<vtkDoubleArray>> VectorOfVariables;
  int counter = 0;
  int var_Counter = 0;
  std::string var_string;

  vtkSmartPointer<vtkPCAStatistics> pcaStatistics = vtkSmartPointer<vtkPCAStatistics>::New();
#if VTK_MAJOR_VERSION <= 5
  pcaStatistics->SetInput(vtkStatisticsAlgorithm::INPUT_DATA, datasetTable);
#else
  pcaStatistics->SetInputData(vtkStatisticsAlgorithm::INPUT_DATA, datasetTable);
#endif

for (int i = 0; i < NumberOfFeatures; i++)
{
    vtkSmartPointer<vtkDoubleArray> CurrentVariable = vtkSmartPointer<vtkDoubleArray>::New();
    CurrentVariable->SetNumberOfComponents(1);
    if (var_Counter == 0)
      var_string = "A";
    else if (var_Counter == 1)
      var_string = "B";
    else if (var_Counter == 2)
      var_string = "C";
    else if (var_Counter == 3)
      var_string = "D";
    else if (var_Counter == 4)
      var_string = "E";
    else if(var_Counter == 5)
      var_string = "F";
    else if (var_Counter == 6)
      var_string = "G";
    else if (var_Counter == 7)
      var_string = "H";
    else if (var_Counter == 8)
      var_string = "I";
    else if (var_Counter == 9)
      var_string = "J";

    std::string namefinal = var_string + std::to_string(counter);
    CurrentVariable->SetName(namefinal.c_str());

    for (size_t index = 0; index < NumberOfSamples; index++)
      CurrentVariable->InsertNextValue(intensities[index][i]);
    datasetTable->AddColumn(CurrentVariable);
    pcaStatistics->SetColumnStatus(namefinal.c_str(), 1);
    counter++;
    if (counter == 9)
    {
      counter = 0;
      var_Counter++;
    }
  }
  pcaStatistics->RequestSelectedColumns();
  pcaStatistics->SetDeriveOption(true);
  pcaStatistics->Update();

  VariableSizeMatrixType transposePCAMatrix;
  transposePCAMatrix.SetSize(NumberOfFeatures, NumberOfFeatures);

  vtkSmartPointer<vtkDoubleArray> eigenvectors = vtkSmartPointer<vtkDoubleArray>::New();
  pcaStatistics->GetEigenvectors(eigenvectors);

  for (vtkIdType i = 0; i < static_cast<vtkIdType>(NumberOfFeatures); i++)
  {
    double* evec = new double[eigenvectors->GetNumberOfComponents()];
    eigenvectors->GetTuple(i, evec);
    for (vtkIdType j = 0; j < eigenvectors->GetNumberOfComponents(); j++)
      transposePCAMatrix[i][j] = evec[j];
    delete evec;
  }

  PCATransformationMatrix = this->MatrixTranspose(transposePCAMatrix);

  vtkSmartPointer<vtkTable> projectedDatasetTable = vtkSmartPointer<vtkTable>::New();
  for (vtkIdType r = 0; r < static_cast<vtkIdType>(NumberOfFeatures); r++)
  {
    vtkSmartPointer<vtkDoubleArray> col = vtkSmartPointer<vtkDoubleArray>::New();
    for (vtkIdType c = 0; c < static_cast<vtkIdType>(NumberOfSamples); c++)
      col->InsertNextValue(0);
    projectedDatasetTable->AddColumn(col);
  }

  for (size_t i = 0; i < NumberOfSamples; i++)
  {
    for (size_t j = 0; j < NumberOfFeatures; j++)
    {
      double sum = 0;
      for (size_t k = 0; k < NumberOfFeatures; k++)
        sum = sum + datasetTable->GetValue(i, k).ToDouble()*PCATransformationMatrix[k][j];
      projectedDatasetTable->SetValue(i, j, vtkVariant(sum));
    }
  }
  VariableLengthVectorType mMeanVector = this->ComputeMeanOfGivenFeatureVectors(projectedDatasetTable);
  for (size_t c = 0; c < NumberOfFeatures; c++)
    for (size_t r = 0; r < NumberOfSamples; r++)
      projectedDatasetTable->SetValue(r, c, projectedDatasetTable->GetValue(r, c).ToDouble() - mMeanVector[c]);

  return projectedDatasetTable;
}

vtkSmartPointer< vtkTable >  FeatureReductionClass::GetDiscerningPerfusionTimePoints(VectorVectorDouble &intensities)
{
  mPMeanvector = ComputeMeanOfGivenFeatureVectors(intensities);
  size_t NumberOfFeatures = intensities[0].size();
  size_t NumberOfSamples = intensities.size();

  vtkSmartPointer<vtkTable> datasetTable = vtkSmartPointer<vtkTable>::New();
  vtkSmartPointer<vtkDoubleArray> A0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> A1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> B0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> B1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> C0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> C1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> D0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> D1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> E0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> E1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> F0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> F1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> G0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> G1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> H0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> H1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> I0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> I1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> J0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> J1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> K0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> K1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> L0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> L1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> M0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> M1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> N0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> N1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> O0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> O1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> P0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> P1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Q0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Q1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> R0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> R1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> S0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> S1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> T0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> T1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> U0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> U1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> V0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> V1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> W0 = vtkSmartPointer<vtkDoubleArray>::New();

  A0->SetNumberOfComponents(1);
  A1->SetNumberOfComponents(1);
  B0->SetNumberOfComponents(1);
  B1->SetNumberOfComponents(1);
  C0->SetNumberOfComponents(1);
  C1->SetNumberOfComponents(1);
  D0->SetNumberOfComponents(1);
  D1->SetNumberOfComponents(1);
  E0->SetNumberOfComponents(1);
  E1->SetNumberOfComponents(1);
  F0->SetNumberOfComponents(1);
  F1->SetNumberOfComponents(1);
  G0->SetNumberOfComponents(1);
  G1->SetNumberOfComponents(1);
  H0->SetNumberOfComponents(1);
  H1->SetNumberOfComponents(1);
  I0->SetNumberOfComponents(1);
  I1->SetNumberOfComponents(1);
  J0->SetNumberOfComponents(1);
  J1->SetNumberOfComponents(1);
  K0->SetNumberOfComponents(1);
  K1->SetNumberOfComponents(1);
  L0->SetNumberOfComponents(1);
  L1->SetNumberOfComponents(1);
  M0->SetNumberOfComponents(1);
  M1->SetNumberOfComponents(1);
  N0->SetNumberOfComponents(1);
  N1->SetNumberOfComponents(1);
  O0->SetNumberOfComponents(1);
  O1->SetNumberOfComponents(1);
  P0->SetNumberOfComponents(1);
  P1->SetNumberOfComponents(1);
  Q0->SetNumberOfComponents(1);
  Q1->SetNumberOfComponents(1);
  R0->SetNumberOfComponents(1);
  R1->SetNumberOfComponents(1);
  S0->SetNumberOfComponents(1);
  S1->SetNumberOfComponents(1);
  T0->SetNumberOfComponents(1);
  T1->SetNumberOfComponents(1);
  U0->SetNumberOfComponents(1);
  U1->SetNumberOfComponents(1);
  V0->SetNumberOfComponents(1);
  V1->SetNumberOfComponents(1);
  W0->SetNumberOfComponents(1);

  A0->SetName("a0");
  A1->SetName("a1");
  B0->SetName("b0");
  B1->SetName("b1");
  C0->SetName("c0");
  C1->SetName("c1");
  D0->SetName("d0");
  D1->SetName("d1");
  E0->SetName("e0");
  E1->SetName("e1");
  F0->SetName("f0");
  F1->SetName("f1");
  G0->SetName("g0");
  G1->SetName("g1");
  H0->SetName("h0");
  H1->SetName("h1");
  I0->SetName("i0");
  I1->SetName("i1");
  J0->SetName("j0");
  J1->SetName("j1");
  K0->SetName("k0");
  K1->SetName("k1");
  L0->SetName("l0");
  L1->SetName("l1");
  M0->SetName("m0");
  M1->SetName("m1");
  N0->SetName("n0");
  N1->SetName("n1");
  O0->SetName("o0");
  O1->SetName("o1");
  P0->SetName("p0");
  P1->SetName("p1");
  Q0->SetName("q0");
  Q1->SetName("q1");
  R0->SetName("r0");
  R1->SetName("r1");
  S0->SetName("s0");
  S1->SetName("s1");
  T0->SetName("t0");
  T1->SetName("t1");
  U0->SetName("u0");
  U1->SetName("u1");
  V0->SetName("v0");
  V1->SetName("v1");
  W0->SetName("w0");

  // this should be made parallel
  for (size_t index = 0; index < NumberOfSamples; index++)
  {
    A0->InsertNextValue(intensities[index][0]);
    A1->InsertNextValue(intensities[index][1]);
    B0->InsertNextValue(intensities[index][2]);
    B1->InsertNextValue(intensities[index][3]);
    C0->InsertNextValue(intensities[index][4]);
    C1->InsertNextValue(intensities[index][5]);
    D0->InsertNextValue(intensities[index][6]);
    D1->InsertNextValue(intensities[index][7]);
    E0->InsertNextValue(intensities[index][8]);
    E1->InsertNextValue(intensities[index][9]);
    F0->InsertNextValue(intensities[index][10]);
    F1->InsertNextValue(intensities[index][11]);
    G0->InsertNextValue(intensities[index][12]);
    G1->InsertNextValue(intensities[index][13]);
    H0->InsertNextValue(intensities[index][14]);
    H1->InsertNextValue(intensities[index][15]);
    I0->InsertNextValue(intensities[index][16]);
    I1->InsertNextValue(intensities[index][17]);
    J0->InsertNextValue(intensities[index][18]);
    J1->InsertNextValue(intensities[index][19]);
    K0->InsertNextValue(intensities[index][20]);
    K1->InsertNextValue(intensities[index][21]);
    L0->InsertNextValue(intensities[index][22]);
    L1->InsertNextValue(intensities[index][23]);
    M0->InsertNextValue(intensities[index][24]);
    M1->InsertNextValue(intensities[index][25]);
    N0->InsertNextValue(intensities[index][26]);
    N1->InsertNextValue(intensities[index][27]);
    O0->InsertNextValue(intensities[index][28]);
    O1->InsertNextValue(intensities[index][29]);
    P0->InsertNextValue(intensities[index][30]);
    P1->InsertNextValue(intensities[index][31]);
    Q0->InsertNextValue(intensities[index][32]);
    Q1->InsertNextValue(intensities[index][33]);
    R0->InsertNextValue(intensities[index][34]);
    R1->InsertNextValue(intensities[index][35]);
    S0->InsertNextValue(intensities[index][36]);
    S1->InsertNextValue(intensities[index][37]);
    T0->InsertNextValue(intensities[index][38]);
    T1->InsertNextValue(intensities[index][39]);
    U0->InsertNextValue(intensities[index][40]);
    U1->InsertNextValue(intensities[index][41]);
    V0->InsertNextValue(intensities[index][42]);
    V1->InsertNextValue(intensities[index][43]);
    W0->InsertNextValue(intensities[index][44]);
  }
  datasetTable->AddColumn(A0);
  datasetTable->AddColumn(A1);
  datasetTable->AddColumn(B0);
  datasetTable->AddColumn(B1);
  datasetTable->AddColumn(C0);
  datasetTable->AddColumn(C1);
  datasetTable->AddColumn(D0);
  datasetTable->AddColumn(D1);
  datasetTable->AddColumn(E0);
  datasetTable->AddColumn(E1);
  datasetTable->AddColumn(F0);
  datasetTable->AddColumn(F1);
  datasetTable->AddColumn(G0);
  datasetTable->AddColumn(G1);
  datasetTable->AddColumn(H0);
  datasetTable->AddColumn(H1);
  datasetTable->AddColumn(I0);
  datasetTable->AddColumn(I1);
  datasetTable->AddColumn(J0);
  datasetTable->AddColumn(J1);
  datasetTable->AddColumn(K0);
  datasetTable->AddColumn(K1);
  datasetTable->AddColumn(L0);
  datasetTable->AddColumn(L1);
  datasetTable->AddColumn(M0);
  datasetTable->AddColumn(M1);
  datasetTable->AddColumn(N0);
  datasetTable->AddColumn(N1);
  datasetTable->AddColumn(O0);
  datasetTable->AddColumn(O1);
  datasetTable->AddColumn(P0);
  datasetTable->AddColumn(P1);
  datasetTable->AddColumn(Q0);
  datasetTable->AddColumn(Q1);
  datasetTable->AddColumn(R0);
  datasetTable->AddColumn(R1);
  datasetTable->AddColumn(S0);
  datasetTable->AddColumn(S1);
  datasetTable->AddColumn(T0);
  datasetTable->AddColumn(T1);
  datasetTable->AddColumn(U0);
  datasetTable->AddColumn(U1);
  datasetTable->AddColumn(V0);
  datasetTable->AddColumn(V1);
  datasetTable->AddColumn(W0);

  vtkSmartPointer<vtkPCAStatistics> pcaStatistics = vtkSmartPointer<vtkPCAStatistics>::New();
#if VTK_MAJOR_VERSION <= 5
  pcaStatistics->SetInput(vtkStatisticsAlgorithm::INPUT_DATA, datasetTable);
#else
  pcaStatistics->SetInputData(vtkStatisticsAlgorithm::INPUT_DATA, datasetTable);
#endif

  pcaStatistics->SetColumnStatus("a0", 1);
  pcaStatistics->SetColumnStatus("a1", 1);
  pcaStatistics->SetColumnStatus("b0", 1);
  pcaStatistics->SetColumnStatus("b1", 1);
  pcaStatistics->SetColumnStatus("c0", 1);
  pcaStatistics->SetColumnStatus("c1", 1);
  pcaStatistics->SetColumnStatus("d0", 1);
  pcaStatistics->SetColumnStatus("d1", 1);
  pcaStatistics->SetColumnStatus("e0", 1);
  pcaStatistics->SetColumnStatus("e1", 1);
  pcaStatistics->SetColumnStatus("f0", 1);
  pcaStatistics->SetColumnStatus("f1", 1);
  pcaStatistics->SetColumnStatus("g0", 1);
  pcaStatistics->SetColumnStatus("g1", 1);
  pcaStatistics->SetColumnStatus("h0", 1);
  pcaStatistics->SetColumnStatus("h1", 1);
  pcaStatistics->SetColumnStatus("i0", 1);
  pcaStatistics->SetColumnStatus("i1", 1);
  pcaStatistics->SetColumnStatus("j0", 1);
  pcaStatistics->SetColumnStatus("j1", 1);
  pcaStatistics->SetColumnStatus("k0", 1);
  pcaStatistics->SetColumnStatus("k1", 1);
  pcaStatistics->SetColumnStatus("l0", 1);
  pcaStatistics->SetColumnStatus("l1", 1);
  pcaStatistics->SetColumnStatus("m0", 1);
  pcaStatistics->SetColumnStatus("m1", 1);
  pcaStatistics->SetColumnStatus("n0", 1);
  pcaStatistics->SetColumnStatus("n1", 1);
  pcaStatistics->SetColumnStatus("o0", 1);
  pcaStatistics->SetColumnStatus("o1", 1);
  pcaStatistics->SetColumnStatus("p0", 1);
  pcaStatistics->SetColumnStatus("p1", 1);
  pcaStatistics->SetColumnStatus("q0", 1);
  pcaStatistics->SetColumnStatus("q1", 1);
  pcaStatistics->SetColumnStatus("r0", 1);
  pcaStatistics->SetColumnStatus("r1", 1);
  pcaStatistics->SetColumnStatus("s0", 1);
  pcaStatistics->SetColumnStatus("s1", 1);
  pcaStatistics->SetColumnStatus("t0", 1);
  pcaStatistics->SetColumnStatus("t1", 1);
  pcaStatistics->SetColumnStatus("u0", 1);
  pcaStatistics->SetColumnStatus("u1", 1);
  pcaStatistics->SetColumnStatus("v0", 1);
  pcaStatistics->SetColumnStatus("v1", 1);
  pcaStatistics->SetColumnStatus("w0", 1);

  pcaStatistics->RequestSelectedColumns();
  pcaStatistics->SetDeriveOption(true);
  pcaStatistics->Update();

  VariableSizeMatrixType transposePCAMatrix;
  transposePCAMatrix.SetSize(NumberOfFeatures, NumberOfFeatures);

  vtkSmartPointer<vtkDoubleArray> eigenvectors = vtkSmartPointer<vtkDoubleArray>::New();
  pcaStatistics->GetEigenvectors(eigenvectors);

  for (vtkIdType i = 0; i < static_cast<vtkIdType>(NumberOfFeatures); i++)
  {
    double* evec = new double[eigenvectors->GetNumberOfComponents()];
    eigenvectors->GetTuple(i, evec);
    for (vtkIdType j = 0; j < eigenvectors->GetNumberOfComponents(); j++)
      transposePCAMatrix[i][j] = evec[j];
    delete evec;
  }

  PCATransformationMatrix = this->MatrixTranspose(transposePCAMatrix);

  vtkSmartPointer<vtkTable> projectedDatasetTable = vtkSmartPointer<vtkTable>::New();
  for (vtkIdType r = 0; r < static_cast<vtkIdType>(NumberOfFeatures); r++)
  {
    vtkSmartPointer<vtkDoubleArray> col = vtkSmartPointer<vtkDoubleArray>::New();
    for (vtkIdType c = 0; c < static_cast<vtkIdType>(NumberOfSamples); c++)
      col->InsertNextValue(0);
    projectedDatasetTable->AddColumn(col);
  }

  for (size_t i = 0; i < NumberOfSamples; i++)
  {
    for (size_t j = 0; j < NumberOfFeatures; j++)
    {
      double sum = 0;
      for (size_t k = 0; k < NumberOfFeatures; k++)
        sum = sum + datasetTable->GetValue(i, k).ToDouble()*PCATransformationMatrix[k][j];
      projectedDatasetTable->SetValue(i, j, vtkVariant(sum));
    }
  }
  VariableLengthVectorType mMeanVector = this->ComputeMeanOfGivenFeatureVectors(projectedDatasetTable);
  for (size_t c = 0; c < NumberOfFeatures; c++)
    for (size_t r = 0; r < NumberOfSamples; r++)
      projectedDatasetTable->SetValue(r, c, projectedDatasetTable->GetValue(r, c).ToDouble() - mMeanVector[c]);

  return projectedDatasetTable;
}

vtkSmartPointer< vtkTable >  FeatureReductionClass::GetDiscerningPerfusionTimePoints(vnl_matrix<double> &intensities)
{
  mPMeanvector = ComputeMeanOfGivenFeatureVectors(intensities);
  size_t NumberOfFeatures = intensities.columns();
  size_t NumberOfSamples = intensities.rows();

  vtkSmartPointer<vtkTable> datasetTable = vtkSmartPointer<vtkTable>::New();
  vtkSmartPointer<vtkDoubleArray> A0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> A1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> B0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> B1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> C0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> C1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> D0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> D1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> E0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> E1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> F0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> F1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> G0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> G1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> H0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> H1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> I0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> I1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> J0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> J1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> K0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> K1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> L0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> L1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> M0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> M1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> N0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> N1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> O0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> O1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> P0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> P1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Q0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Q1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> R0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> R1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> S0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> S1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> T0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> T1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> U0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> U1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> V0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> V1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> W0 = vtkSmartPointer<vtkDoubleArray>::New();

  A0->SetNumberOfComponents(1);
  A1->SetNumberOfComponents(1);
  B0->SetNumberOfComponents(1);
  B1->SetNumberOfComponents(1);
  C0->SetNumberOfComponents(1);
  C1->SetNumberOfComponents(1);
  D0->SetNumberOfComponents(1);
  D1->SetNumberOfComponents(1);
  E0->SetNumberOfComponents(1);
  E1->SetNumberOfComponents(1);
  F0->SetNumberOfComponents(1);
  F1->SetNumberOfComponents(1);
  G0->SetNumberOfComponents(1);
  G1->SetNumberOfComponents(1);
  H0->SetNumberOfComponents(1);
  H1->SetNumberOfComponents(1);
  I0->SetNumberOfComponents(1);
  I1->SetNumberOfComponents(1);
  J0->SetNumberOfComponents(1);
  J1->SetNumberOfComponents(1);
  K0->SetNumberOfComponents(1);
  K1->SetNumberOfComponents(1);
  L0->SetNumberOfComponents(1);
  L1->SetNumberOfComponents(1);
  M0->SetNumberOfComponents(1);
  M1->SetNumberOfComponents(1);
  N0->SetNumberOfComponents(1);
  N1->SetNumberOfComponents(1);
  O0->SetNumberOfComponents(1);
  O1->SetNumberOfComponents(1);
  P0->SetNumberOfComponents(1);
  P1->SetNumberOfComponents(1);
  Q0->SetNumberOfComponents(1);
  Q1->SetNumberOfComponents(1);
  R0->SetNumberOfComponents(1);
  R1->SetNumberOfComponents(1);
  S0->SetNumberOfComponents(1);
  S1->SetNumberOfComponents(1);
  T0->SetNumberOfComponents(1);
  T1->SetNumberOfComponents(1);
  U0->SetNumberOfComponents(1);
  U1->SetNumberOfComponents(1);
  V0->SetNumberOfComponents(1);
  V1->SetNumberOfComponents(1);
  W0->SetNumberOfComponents(1);

  A0->SetName("a0");
  A1->SetName("a1");
  B0->SetName("b0");
  B1->SetName("b1");
  C0->SetName("c0");
  C1->SetName("c1");
  D0->SetName("d0");
  D1->SetName("d1");
  E0->SetName("e0");
  E1->SetName("e1");
  F0->SetName("f0");
  F1->SetName("f1");
  G0->SetName("g0");
  G1->SetName("g1");
  H0->SetName("h0");
  H1->SetName("h1");
  I0->SetName("i0");
  I1->SetName("i1");
  J0->SetName("j0");
  J1->SetName("j1");
  K0->SetName("k0");
  K1->SetName("k1");
  L0->SetName("l0");
  L1->SetName("l1");
  M0->SetName("m0");
  M1->SetName("m1");
  N0->SetName("n0");
  N1->SetName("n1");
  O0->SetName("o0");
  O1->SetName("o1");
  P0->SetName("p0");
  P1->SetName("p1");
  Q0->SetName("q0");
  Q1->SetName("q1");
  R0->SetName("r0");
  R1->SetName("r1");
  S0->SetName("s0");
  S1->SetName("s1");
  T0->SetName("t0");
  T1->SetName("t1");
  U0->SetName("u0");
  U1->SetName("u1");
  V0->SetName("v0");
  V1->SetName("v1");
  W0->SetName("w0");

  // this should be made parallel
  for (size_t index = 0; index < NumberOfSamples; index++)
  {
    A0->InsertNextValue(intensities[index][0]);
    A1->InsertNextValue(intensities[index][1]);
    B0->InsertNextValue(intensities[index][2]);
    B1->InsertNextValue(intensities[index][3]);
    C0->InsertNextValue(intensities[index][4]);
    C1->InsertNextValue(intensities[index][5]);
    D0->InsertNextValue(intensities[index][6]);
    D1->InsertNextValue(intensities[index][7]);
    E0->InsertNextValue(intensities[index][8]);
    E1->InsertNextValue(intensities[index][9]);
    F0->InsertNextValue(intensities[index][10]);
    F1->InsertNextValue(intensities[index][11]);
    G0->InsertNextValue(intensities[index][12]);
    G1->InsertNextValue(intensities[index][13]);
    H0->InsertNextValue(intensities[index][14]);
    H1->InsertNextValue(intensities[index][15]);
    I0->InsertNextValue(intensities[index][16]);
    I1->InsertNextValue(intensities[index][17]);
    J0->InsertNextValue(intensities[index][18]);
    J1->InsertNextValue(intensities[index][19]);
    K0->InsertNextValue(intensities[index][20]);
    K1->InsertNextValue(intensities[index][21]);
    L0->InsertNextValue(intensities[index][22]);
    L1->InsertNextValue(intensities[index][23]);
    M0->InsertNextValue(intensities[index][24]);
    M1->InsertNextValue(intensities[index][25]);
    N0->InsertNextValue(intensities[index][26]);
    N1->InsertNextValue(intensities[index][27]);
    O0->InsertNextValue(intensities[index][28]);
    O1->InsertNextValue(intensities[index][29]);
    P0->InsertNextValue(intensities[index][30]);
    P1->InsertNextValue(intensities[index][31]);
    Q0->InsertNextValue(intensities[index][32]);
    Q1->InsertNextValue(intensities[index][33]);
    R0->InsertNextValue(intensities[index][34]);
    R1->InsertNextValue(intensities[index][35]);
    S0->InsertNextValue(intensities[index][36]);
    S1->InsertNextValue(intensities[index][37]);
    T0->InsertNextValue(intensities[index][38]);
    T1->InsertNextValue(intensities[index][39]);
    U0->InsertNextValue(intensities[index][40]);
    U1->InsertNextValue(intensities[index][41]);
    V0->InsertNextValue(intensities[index][42]);
    V1->InsertNextValue(intensities[index][43]);
    W0->InsertNextValue(intensities[index][44]);
  }
  datasetTable->AddColumn(A0);
  datasetTable->AddColumn(A1);
  datasetTable->AddColumn(B0);
  datasetTable->AddColumn(B1);
  datasetTable->AddColumn(C0);
  datasetTable->AddColumn(C1);
  datasetTable->AddColumn(D0);
  datasetTable->AddColumn(D1);
  datasetTable->AddColumn(E0);
  datasetTable->AddColumn(E1);
  datasetTable->AddColumn(F0);
  datasetTable->AddColumn(F1);
  datasetTable->AddColumn(G0);
  datasetTable->AddColumn(G1);
  datasetTable->AddColumn(H0);
  datasetTable->AddColumn(H1);
  datasetTable->AddColumn(I0);
  datasetTable->AddColumn(I1);
  datasetTable->AddColumn(J0);
  datasetTable->AddColumn(J1);
  datasetTable->AddColumn(K0);
  datasetTable->AddColumn(K1);
  datasetTable->AddColumn(L0);
  datasetTable->AddColumn(L1);
  datasetTable->AddColumn(M0);
  datasetTable->AddColumn(M1);
  datasetTable->AddColumn(N0);
  datasetTable->AddColumn(N1);
  datasetTable->AddColumn(O0);
  datasetTable->AddColumn(O1);
  datasetTable->AddColumn(P0);
  datasetTable->AddColumn(P1);
  datasetTable->AddColumn(Q0);
  datasetTable->AddColumn(Q1);
  datasetTable->AddColumn(R0);
  datasetTable->AddColumn(R1);
  datasetTable->AddColumn(S0);
  datasetTable->AddColumn(S1);
  datasetTable->AddColumn(T0);
  datasetTable->AddColumn(T1);
  datasetTable->AddColumn(U0);
  datasetTable->AddColumn(U1);
  datasetTable->AddColumn(V0);
  datasetTable->AddColumn(V1);
  datasetTable->AddColumn(W0);

  vtkSmartPointer<vtkPCAStatistics> pcaStatistics = vtkSmartPointer<vtkPCAStatistics>::New();
#if VTK_MAJOR_VERSION <= 5
  pcaStatistics->SetInput(vtkStatisticsAlgorithm::INPUT_DATA, datasetTable);
#else
  pcaStatistics->SetInputData(vtkStatisticsAlgorithm::INPUT_DATA, datasetTable);
#endif

  pcaStatistics->SetColumnStatus("a0", 1);
  pcaStatistics->SetColumnStatus("a1", 1);
  pcaStatistics->SetColumnStatus("b0", 1);
  pcaStatistics->SetColumnStatus("b1", 1);
  pcaStatistics->SetColumnStatus("c0", 1);
  pcaStatistics->SetColumnStatus("c1", 1);
  pcaStatistics->SetColumnStatus("d0", 1);
  pcaStatistics->SetColumnStatus("d1", 1);
  pcaStatistics->SetColumnStatus("e0", 1);
  pcaStatistics->SetColumnStatus("e1", 1);
  pcaStatistics->SetColumnStatus("f0", 1);
  pcaStatistics->SetColumnStatus("f1", 1);
  pcaStatistics->SetColumnStatus("g0", 1);
  pcaStatistics->SetColumnStatus("g1", 1);
  pcaStatistics->SetColumnStatus("h0", 1);
  pcaStatistics->SetColumnStatus("h1", 1);
  pcaStatistics->SetColumnStatus("i0", 1);
  pcaStatistics->SetColumnStatus("i1", 1);
  pcaStatistics->SetColumnStatus("j0", 1);
  pcaStatistics->SetColumnStatus("j1", 1);
  pcaStatistics->SetColumnStatus("k0", 1);
  pcaStatistics->SetColumnStatus("k1", 1);
  pcaStatistics->SetColumnStatus("l0", 1);
  pcaStatistics->SetColumnStatus("l1", 1);
  pcaStatistics->SetColumnStatus("m0", 1);
  pcaStatistics->SetColumnStatus("m1", 1);
  pcaStatistics->SetColumnStatus("n0", 1);
  pcaStatistics->SetColumnStatus("n1", 1);
  pcaStatistics->SetColumnStatus("o0", 1);
  pcaStatistics->SetColumnStatus("o1", 1);
  pcaStatistics->SetColumnStatus("p0", 1);
  pcaStatistics->SetColumnStatus("p1", 1);
  pcaStatistics->SetColumnStatus("q0", 1);
  pcaStatistics->SetColumnStatus("q1", 1);
  pcaStatistics->SetColumnStatus("r0", 1);
  pcaStatistics->SetColumnStatus("r1", 1);
  pcaStatistics->SetColumnStatus("s0", 1);
  pcaStatistics->SetColumnStatus("s1", 1);
  pcaStatistics->SetColumnStatus("t0", 1);
  pcaStatistics->SetColumnStatus("t1", 1);
  pcaStatistics->SetColumnStatus("u0", 1);
  pcaStatistics->SetColumnStatus("u1", 1);
  pcaStatistics->SetColumnStatus("v0", 1);
  pcaStatistics->SetColumnStatus("v1", 1);
  pcaStatistics->SetColumnStatus("w0", 1);

  pcaStatistics->RequestSelectedColumns();
  pcaStatistics->SetDeriveOption(true);
  pcaStatistics->Update();

  VariableSizeMatrixType transposePCAMatrix;
  transposePCAMatrix.SetSize(NumberOfFeatures, NumberOfFeatures);

  vtkSmartPointer<vtkDoubleArray> eigenvectors = vtkSmartPointer<vtkDoubleArray>::New();
  pcaStatistics->GetEigenvectors(eigenvectors);

  for (vtkIdType i = 0; i < static_cast<vtkIdType>(NumberOfFeatures); i++)
  {
    double* evec = new double[eigenvectors->GetNumberOfComponents()];
    eigenvectors->GetTuple(i, evec);
    for (vtkIdType j = 0; j < eigenvectors->GetNumberOfComponents(); j++)
      transposePCAMatrix[i][j] = evec[j];
    delete evec;
  }

  PCATransformationMatrix = this->MatrixTranspose(transposePCAMatrix);

  vtkSmartPointer<vtkTable> projectedDatasetTable = vtkSmartPointer<vtkTable>::New();
  for (vtkIdType r = 0; r < static_cast<vtkIdType>(NumberOfFeatures); r++)
  {
    vtkSmartPointer<vtkDoubleArray> col = vtkSmartPointer<vtkDoubleArray>::New();
    for (vtkIdType c = 0; c < static_cast<vtkIdType>(NumberOfSamples); c++)
      col->InsertNextValue(0);
    projectedDatasetTable->AddColumn(col);
  }

  for (size_t i = 0; i < NumberOfSamples; i++)
  {
    for (size_t j = 0; j < NumberOfFeatures; j++)
    {
      double sum = 0;
      for (size_t k = 0; k < NumberOfFeatures; k++)
        sum = sum + datasetTable->GetValue(i, k).ToDouble()*PCATransformationMatrix[k][j];
      projectedDatasetTable->SetValue(i, j, vtkVariant(sum));
    }
  }
  VariableLengthVectorType mMeanVector = this->ComputeMeanOfGivenFeatureVectors(projectedDatasetTable);
  for (size_t c = 0; c < NumberOfFeatures; c++)
    for (size_t r = 0; r < NumberOfSamples; r++)
      projectedDatasetTable->SetValue(r, c, projectedDatasetTable->GetValue(r, c).ToDouble() - mMeanVector[c]);

  return projectedDatasetTable;
}

VectorVectorDouble  FeatureReductionClass::ApplyPCAOnTestData(VectorVectorDouble &intensities)
{
  size_t NumberOfFeatures = intensities[0].size();
  size_t NumberOfSamples = intensities.size();

  for (size_t index = 0; index < NumberOfSamples; index++)
  {
    for (size_t timepoint = 0; timepoint < NumberOfFeatures; timepoint++)
      intensities[index][timepoint] = intensities[index][timepoint] - mPMeanvector[timepoint];
  }
  vtkSmartPointer<vtkTable> projectedDatasetTable = vtkSmartPointer<vtkTable>::New();
  for (vtkIdType r = 0; r < static_cast<vtkIdType>(NumberOfFeatures); r++)
  {
    vtkSmartPointer<vtkDoubleArray> col = vtkSmartPointer<vtkDoubleArray>::New();
    for (vtkIdType c = 0; c < static_cast<vtkIdType>(NumberOfSamples); c++)
      col->InsertNextValue(0);
    projectedDatasetTable->AddColumn(col);
  }
  VectorVectorDouble projectedData;

  for (size_t i = 0; i < NumberOfSamples; i++)
  {
    VectorDouble oneDataPoint;
    for (size_t j = 0; j < NumberOfFeatures; j++)
    {
      double sum = 0;
      for (size_t k = 0; k < NumberOfFeatures; k++)
        sum = sum + intensities[i][k] * PCATransformationMatrix[k][j];
      oneDataPoint.push_back(sum);
    }
    projectedData.push_back(oneDataPoint);
  }
  return projectedData;
}


vtkSmartPointer< vtkTable >  FeatureReductionClass::GetDiscerningPerfusionTimePointsFullPCA(VectorVectorDouble &intensities, VariableSizeMatrixType &TransformationMatrix, VariableLengthVectorType &MeanVector)
{
  mPMeanvector = ComputeMeanOfGivenFeatureVectors(intensities);
  size_t NumberOfFeatures = intensities[0].size();
  size_t NumberOfSamples = intensities.size();

  vtkSmartPointer<vtkTable> datasetTable = vtkSmartPointer<vtkTable>::New();
  vtkSmartPointer<vtkDoubleArray> A0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> A1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> A2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> A3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> A4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> A5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> A6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> A7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> A8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> A9 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> B0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> B1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> B2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> B3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> B4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> B5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> B6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> B7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> B8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> B9 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> C0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> C1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> C2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> C3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> C4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> C5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> C6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> C7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> C8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> C9 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> D0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> D1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> D2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> D3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> D4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> D5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> D6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> D7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> D8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> D9 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> E0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> E1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> E2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> E3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> E4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> E5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> E6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> E7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> E8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> E9 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> F0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> F1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> F2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> F3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> F4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> F5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> F6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> F7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> F8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> F9 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> G0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> G1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> G2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> G3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> G4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> G5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> G6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> G7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> G8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> G9 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> H0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> H1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> H2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> H3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> H4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> H5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> H6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> H7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> H8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> H9 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> I0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> I1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> I2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> I3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> I4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> I5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> I6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> I7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> I8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> I9 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> J0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> J1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> J2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> J3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> J4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> J5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> J6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> J7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> J8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> J9 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> K0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> K1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> K2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> K3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> K4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> K5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> K6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> K7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> K8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> K9 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> L0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> L1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> L2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> L3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> L4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> L5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> L6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> L7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> L8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> L9 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> M0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> M1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> M2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> M3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> M4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> M5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> M6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> M7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> M8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> M9 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> N0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> N1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> N2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> N3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> N4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> N5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> N6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> N7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> N8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> N9 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> O0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> O1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> O2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> O3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> O4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> O5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> O6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> O7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> O8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> O9 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> P0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> P1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> P2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> P3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> P4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> P5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> P6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> P7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> P8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> P9 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Q0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Q1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Q2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Q3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Q4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Q5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Q6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Q7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Q8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Q9 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> R0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> R1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> R2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> R3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> R4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> R5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> R6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> R7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> R8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> R9 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> S0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> S1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> S2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> S3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> S4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> S5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> S6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> S7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> S8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> S9 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> T0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> T1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> T2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> T3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> T4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> T5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> T6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> T7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> T8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> T9 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> U0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> U1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> U2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> U3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> U4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> U5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> U6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> U7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> U8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> U9 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> V0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> V1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> V2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> V3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> V4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> V5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> V6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> V7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> V8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> V9 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> W0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> W1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> W2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> W3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> W4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> W5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> W6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> W7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> W8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> W9 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> X0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> X1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> X2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> X3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> X4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> X5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> X6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> X7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> X8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> X9 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Y0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Y1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Y2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Y3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Y4 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Y5 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Y6 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Y7 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Y8 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Y9 = vtkSmartPointer<vtkDoubleArray>::New();

  vtkSmartPointer<vtkDoubleArray> Z0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Z1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Z2 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Z3 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Z4 = vtkSmartPointer<vtkDoubleArray>::New();

  A0->SetNumberOfComponents(1);
  A1->SetNumberOfComponents(1);
  A2->SetNumberOfComponents(1);
  A3->SetNumberOfComponents(1);
  A4->SetNumberOfComponents(1);
  A5->SetNumberOfComponents(1);
  A6->SetNumberOfComponents(1);
  A7->SetNumberOfComponents(1);
  A8->SetNumberOfComponents(1);
  A9->SetNumberOfComponents(1);
  B0->SetNumberOfComponents(1);
  B1->SetNumberOfComponents(1);
  B2->SetNumberOfComponents(1);
  B3->SetNumberOfComponents(1);
  B4->SetNumberOfComponents(1);
  B5->SetNumberOfComponents(1);
  B6->SetNumberOfComponents(1);
  B7->SetNumberOfComponents(1);
  B8->SetNumberOfComponents(1);
  B9->SetNumberOfComponents(1);
  C0->SetNumberOfComponents(1);
  C1->SetNumberOfComponents(1);
  C2->SetNumberOfComponents(1);
  C3->SetNumberOfComponents(1);
  C4->SetNumberOfComponents(1);
  C5->SetNumberOfComponents(1);
  C6->SetNumberOfComponents(1);
  C7->SetNumberOfComponents(1);
  C8->SetNumberOfComponents(1);
  C9->SetNumberOfComponents(1);
  D0->SetNumberOfComponents(1);
  D1->SetNumberOfComponents(1);
  D2->SetNumberOfComponents(1);
  D3->SetNumberOfComponents(1);
  D4->SetNumberOfComponents(1);
  D5->SetNumberOfComponents(1);
  D6->SetNumberOfComponents(1);
  D7->SetNumberOfComponents(1);
  D8->SetNumberOfComponents(1);
  D9->SetNumberOfComponents(1);
  E0->SetNumberOfComponents(1);
  E1->SetNumberOfComponents(1);
  E2->SetNumberOfComponents(1);
  E3->SetNumberOfComponents(1);
  E4->SetNumberOfComponents(1);
  E5->SetNumberOfComponents(1);
  E6->SetNumberOfComponents(1);
  E7->SetNumberOfComponents(1);
  E8->SetNumberOfComponents(1);
  E9->SetNumberOfComponents(1);
  F0->SetNumberOfComponents(1);
  F1->SetNumberOfComponents(1);
  F2->SetNumberOfComponents(1);
  F3->SetNumberOfComponents(1);
  F4->SetNumberOfComponents(1);
  F5->SetNumberOfComponents(1);
  F6->SetNumberOfComponents(1);
  F7->SetNumberOfComponents(1);
  F8->SetNumberOfComponents(1);
  F9->SetNumberOfComponents(1);
  G0->SetNumberOfComponents(1);
  G1->SetNumberOfComponents(1);
  G2->SetNumberOfComponents(1);
  G3->SetNumberOfComponents(1);
  G4->SetNumberOfComponents(1);
  G5->SetNumberOfComponents(1);
  G6->SetNumberOfComponents(1);
  G7->SetNumberOfComponents(1);
  G8->SetNumberOfComponents(1);
  G9->SetNumberOfComponents(1);
  H0->SetNumberOfComponents(1);
  H1->SetNumberOfComponents(1);
  H2->SetNumberOfComponents(1);
  H3->SetNumberOfComponents(1);
  H4->SetNumberOfComponents(1);
  H5->SetNumberOfComponents(1);
  H6->SetNumberOfComponents(1);
  H7->SetNumberOfComponents(1);
  H8->SetNumberOfComponents(1);
  H9->SetNumberOfComponents(1);
  I0->SetNumberOfComponents(1);
  I1->SetNumberOfComponents(1);
  I2->SetNumberOfComponents(1);
  I3->SetNumberOfComponents(1);
  I4->SetNumberOfComponents(1);
  I5->SetNumberOfComponents(1);
  I6->SetNumberOfComponents(1);
  I7->SetNumberOfComponents(1);
  I8->SetNumberOfComponents(1);
  I9->SetNumberOfComponents(1);
  J0->SetNumberOfComponents(1);
  J1->SetNumberOfComponents(1);
  J2->SetNumberOfComponents(1);
  J3->SetNumberOfComponents(1);
  J4->SetNumberOfComponents(1);
  J5->SetNumberOfComponents(1);
  J6->SetNumberOfComponents(1);
  J7->SetNumberOfComponents(1);
  J8->SetNumberOfComponents(1);
  J9->SetNumberOfComponents(1);
  K0->SetNumberOfComponents(1);
  K1->SetNumberOfComponents(1);
  K2->SetNumberOfComponents(1);
  K3->SetNumberOfComponents(1);
  K4->SetNumberOfComponents(1);
  K5->SetNumberOfComponents(1);
  K6->SetNumberOfComponents(1);
  K7->SetNumberOfComponents(1);
  K8->SetNumberOfComponents(1);
  K9->SetNumberOfComponents(1);
  L0->SetNumberOfComponents(1);
  L1->SetNumberOfComponents(1);
  L2->SetNumberOfComponents(1);
  L3->SetNumberOfComponents(1);
  L4->SetNumberOfComponents(1);
  L5->SetNumberOfComponents(1);
  L6->SetNumberOfComponents(1);
  L7->SetNumberOfComponents(1);
  L8->SetNumberOfComponents(1);
  L9->SetNumberOfComponents(1);
  M0->SetNumberOfComponents(1);
  M1->SetNumberOfComponents(1);
  M2->SetNumberOfComponents(1);
  M3->SetNumberOfComponents(1);
  M4->SetNumberOfComponents(1);
  M5->SetNumberOfComponents(1);
  M6->SetNumberOfComponents(1);
  M7->SetNumberOfComponents(1);
  M8->SetNumberOfComponents(1);
  M9->SetNumberOfComponents(1);
  N0->SetNumberOfComponents(1);
  N1->SetNumberOfComponents(1);
  N2->SetNumberOfComponents(1);
  N3->SetNumberOfComponents(1);
  N4->SetNumberOfComponents(1);
  N5->SetNumberOfComponents(1);
  N6->SetNumberOfComponents(1);
  N7->SetNumberOfComponents(1);
  N8->SetNumberOfComponents(1);
  N9->SetNumberOfComponents(1);
  O0->SetNumberOfComponents(1);
  O1->SetNumberOfComponents(1);
  O2->SetNumberOfComponents(1);
  O3->SetNumberOfComponents(1);
  O4->SetNumberOfComponents(1);
  O5->SetNumberOfComponents(1);
  O6->SetNumberOfComponents(1);
  O7->SetNumberOfComponents(1);
  O8->SetNumberOfComponents(1);
  O9->SetNumberOfComponents(1);
  P0->SetNumberOfComponents(1);
  P1->SetNumberOfComponents(1);
  P2->SetNumberOfComponents(1);
  P3->SetNumberOfComponents(1);
  P4->SetNumberOfComponents(1);
  P5->SetNumberOfComponents(1);
  P6->SetNumberOfComponents(1);
  P7->SetNumberOfComponents(1);
  P8->SetNumberOfComponents(1);
  P9->SetNumberOfComponents(1);
  Q0->SetNumberOfComponents(1);
  Q1->SetNumberOfComponents(1);
  Q2->SetNumberOfComponents(1);
  Q3->SetNumberOfComponents(1);
  Q4->SetNumberOfComponents(1);
  Q5->SetNumberOfComponents(1);
  Q6->SetNumberOfComponents(1);
  Q7->SetNumberOfComponents(1);
  Q8->SetNumberOfComponents(1);
  Q9->SetNumberOfComponents(1);
  R0->SetNumberOfComponents(1);
  R1->SetNumberOfComponents(1);
  R2->SetNumberOfComponents(1);
  R3->SetNumberOfComponents(1);
  R4->SetNumberOfComponents(1);
  R5->SetNumberOfComponents(1);
  R6->SetNumberOfComponents(1);
  R7->SetNumberOfComponents(1);
  R8->SetNumberOfComponents(1);
  R9->SetNumberOfComponents(1);
  S0->SetNumberOfComponents(1);
  S1->SetNumberOfComponents(1);
  S2->SetNumberOfComponents(1);
  S3->SetNumberOfComponents(1);
  S4->SetNumberOfComponents(1);
  S5->SetNumberOfComponents(1);
  S6->SetNumberOfComponents(1);
  S7->SetNumberOfComponents(1);
  S8->SetNumberOfComponents(1);
  S9->SetNumberOfComponents(1);
  T0->SetNumberOfComponents(1);
  T1->SetNumberOfComponents(1);
  T2->SetNumberOfComponents(1);
  T3->SetNumberOfComponents(1);
  T4->SetNumberOfComponents(1);
  T5->SetNumberOfComponents(1);
  T6->SetNumberOfComponents(1);
  T7->SetNumberOfComponents(1);
  T8->SetNumberOfComponents(1);
  T9->SetNumberOfComponents(1);
  U0->SetNumberOfComponents(1);
  U1->SetNumberOfComponents(1);
  U2->SetNumberOfComponents(1);
  U3->SetNumberOfComponents(1);
  U4->SetNumberOfComponents(1);
  U5->SetNumberOfComponents(1);
  U6->SetNumberOfComponents(1);
  U7->SetNumberOfComponents(1);
  U8->SetNumberOfComponents(1);
  U9->SetNumberOfComponents(1);
  V0->SetNumberOfComponents(1);
  V1->SetNumberOfComponents(1);
  V2->SetNumberOfComponents(1);
  V3->SetNumberOfComponents(1);
  V4->SetNumberOfComponents(1);
  V5->SetNumberOfComponents(1);
  V6->SetNumberOfComponents(1);
  V7->SetNumberOfComponents(1);
  V8->SetNumberOfComponents(1);
  V9->SetNumberOfComponents(1);
  W0->SetNumberOfComponents(1);
  W1->SetNumberOfComponents(1);
  W2->SetNumberOfComponents(1);
  W3->SetNumberOfComponents(1);
  W4->SetNumberOfComponents(1);
  W5->SetNumberOfComponents(1);
  W6->SetNumberOfComponents(1);
  W7->SetNumberOfComponents(1);
  W8->SetNumberOfComponents(1);
  W9->SetNumberOfComponents(1);
  X0->SetNumberOfComponents(1);
  X1->SetNumberOfComponents(1);
  X2->SetNumberOfComponents(1);
  X3->SetNumberOfComponents(1);
  X4->SetNumberOfComponents(1);
  X5->SetNumberOfComponents(1);
  X6->SetNumberOfComponents(1);
  X7->SetNumberOfComponents(1);
  X8->SetNumberOfComponents(1);
  X9->SetNumberOfComponents(1);
  Y0->SetNumberOfComponents(1);
  Y1->SetNumberOfComponents(1);
  Y2->SetNumberOfComponents(1);
  Y3->SetNumberOfComponents(1);
  Y4->SetNumberOfComponents(1);
  Y5->SetNumberOfComponents(1);
  Y6->SetNumberOfComponents(1);
  Y7->SetNumberOfComponents(1);
  Y8->SetNumberOfComponents(1);
  Y9->SetNumberOfComponents(1);
  Z0->SetNumberOfComponents(1);
  Z1->SetNumberOfComponents(1);
  Z2->SetNumberOfComponents(1);
  Z3->SetNumberOfComponents(1);
  Z4->SetNumberOfComponents(1);


  A0->SetName("a0");
  A1->SetName("a1");
  A2->SetName("a2");
  A3->SetName("a3");
  A4->SetName("a4");
  A5->SetName("a5");
  A6->SetName("a6");
  A7->SetName("a7");
  A8->SetName("a8");
  A9->SetName("a9");

  B0->SetName("b0");
  B1->SetName("b1");
  B2->SetName("b2");
  B3->SetName("b3");
  B4->SetName("b4");
  B5->SetName("b5");
  B6->SetName("b6");
  B7->SetName("b7");
  B8->SetName("b8");
  B9->SetName("b9");

  C0->SetName("c0");
  C1->SetName("c1");
  C2->SetName("c2");
  C3->SetName("c3");
  C4->SetName("c4");
  C5->SetName("c5");
  C6->SetName("c6");
  C7->SetName("c7");
  C8->SetName("c8");
  C9->SetName("c9");

  D0->SetName("d0");
  D1->SetName("d1");
  D2->SetName("d2");
  D3->SetName("d3");
  D4->SetName("d4");
  D5->SetName("d5");
  D6->SetName("d6");
  D7->SetName("d7");
  D8->SetName("d8");
  D9->SetName("d9");

  E0->SetName("e0");
  E1->SetName("e1");
  E2->SetName("e2");
  E3->SetName("e3");
  E4->SetName("e4");
  E5->SetName("e5");
  E6->SetName("e6");
  E7->SetName("e7");
  E8->SetName("e8");
  E9->SetName("e9");
  F0->SetName("f0");
  F1->SetName("f1");
  F2->SetName("f2");
  F3->SetName("f3");
  F4->SetName("f4");
  F5->SetName("f5");
  F6->SetName("f6");
  F7->SetName("f7");
  F8->SetName("f8");
  F9->SetName("f9");
  G0->SetName("g0");
  G1->SetName("g1");
  G2->SetName("g2");
  G3->SetName("g3");
  G4->SetName("g4");
  G5->SetName("g5");
  G6->SetName("g6");
  G7->SetName("g7");
  G8->SetName("g8");
  G9->SetName("g9");
  H0->SetName("h0");
  H1->SetName("h1");
  H2->SetName("h2");
  H3->SetName("h3");
  H4->SetName("h4");
  H5->SetName("h5");
  H6->SetName("h6");
  H7->SetName("h7");
  H8->SetName("h8");
  H9->SetName("h9");
  I0->SetName("i0");
  I1->SetName("i1");
  I2->SetName("i2");
  I3->SetName("i3");
  I4->SetName("i4");
  I5->SetName("i5");
  I6->SetName("i6");
  I7->SetName("i7");
  I8->SetName("i8");
  I9->SetName("i9");
  J0->SetName("j0");
  J1->SetName("j1");
  J2->SetName("j2");
  J3->SetName("j3");
  J4->SetName("j4");
  J5->SetName("j5");
  J6->SetName("j6");
  J7->SetName("j7");
  J8->SetName("j8");
  J9->SetName("j9");
  K0->SetName("k0");
  K1->SetName("k1");
  K2->SetName("k2");
  K3->SetName("k3");
  K4->SetName("k4");
  K5->SetName("k5");
  K6->SetName("k6");
  K7->SetName("k7");
  K8->SetName("k8");
  K9->SetName("k9");
  L0->SetName("l0");
  L1->SetName("l1");
  L2->SetName("l2");
  L3->SetName("l3");
  L4->SetName("l4");
  L5->SetName("l5");
  L6->SetName("l6");
  L7->SetName("l7");
  L8->SetName("l8");
  L9->SetName("l9");

  M0->SetName("m0");
  M1->SetName("m1");
  M2->SetName("m2");
  M3->SetName("m3");
  M4->SetName("m4");
  M5->SetName("m5");
  M6->SetName("m6");
  M7->SetName("m7");
  M8->SetName("m8");
  M9->SetName("m9");

  N0->SetName("n0");
  N1->SetName("n1");
  N2->SetName("n2");
  N3->SetName("n3");
  N4->SetName("n4");
  N5->SetName("n5");
  N6->SetName("n6");
  N7->SetName("n7");
  N8->SetName("n8");
  N9->SetName("n9");

  O0->SetName("o0");
  O1->SetName("o1");
  O2->SetName("o2");
  O3->SetName("o3");
  O4->SetName("o4");
  O5->SetName("o5");
  O6->SetName("o6");
  O7->SetName("o7");
  O8->SetName("o8");
  O9->SetName("o9");

  P0->SetName("p0");
  P1->SetName("p1");
  P2->SetName("p2");
  P3->SetName("p3");
  P4->SetName("p4");
  P5->SetName("p5");
  P6->SetName("p6");
  P7->SetName("p7");
  P8->SetName("p8");
  P9->SetName("p9");

  Q0->SetName("q0");
  Q1->SetName("q1");
  Q2->SetName("q2");
  Q3->SetName("q3");
  Q4->SetName("q4");
  Q5->SetName("q5");
  Q6->SetName("q6");
  Q7->SetName("q7");
  Q8->SetName("q8");
  Q9->SetName("q9");

  R0->SetName("r0");
  R1->SetName("r1");
  R2->SetName("r2");
  R3->SetName("r3");
  R4->SetName("r4");
  R5->SetName("r5");
  R6->SetName("r6");
  R7->SetName("r7");
  R8->SetName("r8");
  R9->SetName("r9");

  S0->SetName("s0");
  S1->SetName("s1");
  S2->SetName("s2");
  S3->SetName("s3");
  S4->SetName("s4");
  S5->SetName("s5");
  S6->SetName("s6");
  S7->SetName("s7");
  S8->SetName("s8");
  S9->SetName("s9");

  T0->SetName("t0");
  T1->SetName("t1");
  T2->SetName("t2");
  T3->SetName("t3");
  T4->SetName("t4");
  T5->SetName("t5");
  T6->SetName("t6");
  T7->SetName("t7");
  T8->SetName("t8");
  T9->SetName("t9");

  U0->SetName("u0");
  U1->SetName("u1");
  U2->SetName("u2");
  U3->SetName("u3");
  U4->SetName("u4");
  U5->SetName("u5");
  U6->SetName("u6");
  U7->SetName("u7");
  U8->SetName("u8");
  U9->SetName("u9");

  V0->SetName("v0");
  V1->SetName("v1");
  V2->SetName("v2");
  V3->SetName("v3");
  V4->SetName("v4");
  V5->SetName("v5");
  V6->SetName("v6");
  V7->SetName("v7");
  V8->SetName("v8");
  V9->SetName("v9");

  W0->SetName("w0");
  W1->SetName("w1");
  W2->SetName("w2");
  W3->SetName("w3");
  W4->SetName("w4");
  W5->SetName("w5");
  W6->SetName("w6");
  W7->SetName("w7");
  W8->SetName("w8");
  W9->SetName("w9");

  X0->SetName("x0");
  X1->SetName("x1");
  X2->SetName("x2");
  X3->SetName("x3");
  X4->SetName("x4");
  X5->SetName("x5");
  X6->SetName("x6");
  X7->SetName("x7");
  X8->SetName("x8");
  X9->SetName("x9");

  Y0->SetName("y0");
  Y1->SetName("y1");
  Y2->SetName("y2");
  Y3->SetName("y3");
  Y4->SetName("y4");
  Y5->SetName("y5");
  Y6->SetName("y6");
  Y7->SetName("y7");
  Y8->SetName("y8");
  Y9->SetName("y9");

  Z0->SetName("z0");
  Z1->SetName("z1");
  Z2->SetName("z2");
  Z3->SetName("z3");
  Z4->SetName("z4");

  // this should be made parallel
  for (size_t index = 0; index < NumberOfSamples; index++)
  {
    A0->InsertNextValue(intensities[index][0]);
    A1->InsertNextValue(intensities[index][1]);
    A2->InsertNextValue(intensities[index][2]);
    A3->InsertNextValue(intensities[index][3]);
    A4->InsertNextValue(intensities[index][4]);
    A5->InsertNextValue(intensities[index][5]);
    A6->InsertNextValue(intensities[index][6]);
    A7->InsertNextValue(intensities[index][7]);
    A8->InsertNextValue(intensities[index][8]);
    A9->InsertNextValue(intensities[index][9]);
    B0->InsertNextValue(intensities[index][10]);
    B1->InsertNextValue(intensities[index][11]);
    B2->InsertNextValue(intensities[index][12]);
    B3->InsertNextValue(intensities[index][13]);
    B4->InsertNextValue(intensities[index][14]);
    B5->InsertNextValue(intensities[index][15]);
    B6->InsertNextValue(intensities[index][16]);
    B7->InsertNextValue(intensities[index][17]);
    B8->InsertNextValue(intensities[index][18]);
    B9->InsertNextValue(intensities[index][19]);
    C0->InsertNextValue(intensities[index][20]);
    C1->InsertNextValue(intensities[index][21]);
    C2->InsertNextValue(intensities[index][22]);
    C3->InsertNextValue(intensities[index][23]);
    C4->InsertNextValue(intensities[index][24]);
    C5->InsertNextValue(intensities[index][25]);
    C6->InsertNextValue(intensities[index][26]);
    C7->InsertNextValue(intensities[index][27]);
    C8->InsertNextValue(intensities[index][28]);
    C9->InsertNextValue(intensities[index][29]);
    D0->InsertNextValue(intensities[index][30]);

    D1->InsertNextValue(intensities[index][31]);
    D2->InsertNextValue(intensities[index][32]);
    D3->InsertNextValue(intensities[index][33]);
    D4->InsertNextValue(intensities[index][34]);
    D5->InsertNextValue(intensities[index][35]);
    D6->InsertNextValue(intensities[index][36]);
    D7->InsertNextValue(intensities[index][37]);
    D8->InsertNextValue(intensities[index][38]);
    D9->InsertNextValue(intensities[index][39]);
    E0->InsertNextValue(intensities[index][40]);
    E1->InsertNextValue(intensities[index][41]);
    E2->InsertNextValue(intensities[index][42]);
    E3->InsertNextValue(intensities[index][43]);
    E4->InsertNextValue(intensities[index][44]);
    E5->InsertNextValue(intensities[index][45]);
    E6->InsertNextValue(intensities[index][46]);
    E7->InsertNextValue(intensities[index][47]);
    E8->InsertNextValue(intensities[index][48]);
    E9->InsertNextValue(intensities[index][49]);
    F0->InsertNextValue(intensities[index][50]);
    F1->InsertNextValue(intensities[index][51]);
    F2->InsertNextValue(intensities[index][52]);
    F3->InsertNextValue(intensities[index][53]);
    F4->InsertNextValue(intensities[index][54]);
    F5->InsertNextValue(intensities[index][55]);
    F6->InsertNextValue(intensities[index][56]);
    F7->InsertNextValue(intensities[index][57]);
    F8->InsertNextValue(intensities[index][58]);
    F9->InsertNextValue(intensities[index][59]);
    G0->InsertNextValue(intensities[index][60]);

    G1->InsertNextValue(intensities[index][61]);
    G2->InsertNextValue(intensities[index][62]);
    G3->InsertNextValue(intensities[index][63]);
    G4->InsertNextValue(intensities[index][64]);
    G5->InsertNextValue(intensities[index][65]);
    G6->InsertNextValue(intensities[index][66]);
    G7->InsertNextValue(intensities[index][67]);
    G8->InsertNextValue(intensities[index][68]);
    G9->InsertNextValue(intensities[index][69]);
    H0->InsertNextValue(intensities[index][70]);
    H1->InsertNextValue(intensities[index][71]);
    H2->InsertNextValue(intensities[index][72]);
    H3->InsertNextValue(intensities[index][73]);
    H4->InsertNextValue(intensities[index][74]);
    H5->InsertNextValue(intensities[index][75]);
    H6->InsertNextValue(intensities[index][76]);
    H7->InsertNextValue(intensities[index][77]);
    H8->InsertNextValue(intensities[index][78]);
    H9->InsertNextValue(intensities[index][79]);
    I0->InsertNextValue(intensities[index][80]);
    I1->InsertNextValue(intensities[index][81]);
    I2->InsertNextValue(intensities[index][82]);
    I3->InsertNextValue(intensities[index][83]);
    I4->InsertNextValue(intensities[index][84]);
    I5->InsertNextValue(intensities[index][85]);
    I6->InsertNextValue(intensities[index][86]);
    I7->InsertNextValue(intensities[index][87]);
    I8->InsertNextValue(intensities[index][88]);
    I9->InsertNextValue(intensities[index][89]);
    J0->InsertNextValue(intensities[index][90]);

    J1->InsertNextValue(intensities[index][91]);
    J2->InsertNextValue(intensities[index][92]);
    J3->InsertNextValue(intensities[index][93]);
    J4->InsertNextValue(intensities[index][94]);
    J5->InsertNextValue(intensities[index][95]);
    J6->InsertNextValue(intensities[index][96]);
    J7->InsertNextValue(intensities[index][97]);
    J8->InsertNextValue(intensities[index][98]);
    J9->InsertNextValue(intensities[index][99]);
    K0->InsertNextValue(intensities[index][100]);

    K1->InsertNextValue(intensities[index][101]);
    K2->InsertNextValue(intensities[index][102]);
    K3->InsertNextValue(intensities[index][103]);
    K4->InsertNextValue(intensities[index][104]);
    K5->InsertNextValue(intensities[index][105]);
    K6->InsertNextValue(intensities[index][106]);
    K7->InsertNextValue(intensities[index][107]);
    K8->InsertNextValue(intensities[index][108]);
    K9->InsertNextValue(intensities[index][109]);

    L0->InsertNextValue(intensities[index][110]);
    L1->InsertNextValue(intensities[index][111]);
    L2->InsertNextValue(intensities[index][112]);
    L3->InsertNextValue(intensities[index][113]);
    L4->InsertNextValue(intensities[index][114]);
    L5->InsertNextValue(intensities[index][115]);
    L6->InsertNextValue(intensities[index][116]);
    L7->InsertNextValue(intensities[index][117]);
    L8->InsertNextValue(intensities[index][118]);
    L9->InsertNextValue(intensities[index][119]);
    M0->InsertNextValue(intensities[index][120]);

    M1->InsertNextValue(intensities[index][121]);
    M2->InsertNextValue(intensities[index][122]);
    M3->InsertNextValue(intensities[index][123]);
    M4->InsertNextValue(intensities[index][124]);
    M5->InsertNextValue(intensities[index][125]);
    M6->InsertNextValue(intensities[index][126]);
    M7->InsertNextValue(intensities[index][127]);
    M8->InsertNextValue(intensities[index][128]);
    M9->InsertNextValue(intensities[index][129]);
    N0->InsertNextValue(intensities[index][130]);
    N1->InsertNextValue(intensities[index][131]);
    N2->InsertNextValue(intensities[index][132]);
    N3->InsertNextValue(intensities[index][133]);
    N4->InsertNextValue(intensities[index][134]);
    N5->InsertNextValue(intensities[index][135]);
    N6->InsertNextValue(intensities[index][136]);
    N7->InsertNextValue(intensities[index][137]);
    N8->InsertNextValue(intensities[index][138]);
    N9->InsertNextValue(intensities[index][139]);
    O0->InsertNextValue(intensities[index][140]);
    O1->InsertNextValue(intensities[index][141]);
    O2->InsertNextValue(intensities[index][142]);
    O3->InsertNextValue(intensities[index][143]);
    O4->InsertNextValue(intensities[index][144]);
    O5->InsertNextValue(intensities[index][145]);
    O6->InsertNextValue(intensities[index][146]);
    O7->InsertNextValue(intensities[index][147]);
    O8->InsertNextValue(intensities[index][148]);
    O9->InsertNextValue(intensities[index][149]);
    P0->InsertNextValue(intensities[index][150]);

    P1->InsertNextValue(intensities[index][151]);
    P2->InsertNextValue(intensities[index][152]);
    P3->InsertNextValue(intensities[index][153]);
    P4->InsertNextValue(intensities[index][154]);
    P5->InsertNextValue(intensities[index][155]);
    P6->InsertNextValue(intensities[index][156]);
    P7->InsertNextValue(intensities[index][157]);
    P8->InsertNextValue(intensities[index][158]);
    P9->InsertNextValue(intensities[index][159]);
    Q0->InsertNextValue(intensities[index][160]);

    Q1->InsertNextValue(intensities[index][161]);
    Q2->InsertNextValue(intensities[index][162]);
    Q3->InsertNextValue(intensities[index][163]);
    Q4->InsertNextValue(intensities[index][164]);
    Q5->InsertNextValue(intensities[index][165]);
    Q6->InsertNextValue(intensities[index][166]);
    Q7->InsertNextValue(intensities[index][167]);
    Q8->InsertNextValue(intensities[index][168]);
    Q9->InsertNextValue(intensities[index][169]);
    R0->InsertNextValue(intensities[index][170]);
    R1->InsertNextValue(intensities[index][171]);
    R2->InsertNextValue(intensities[index][172]);
    R3->InsertNextValue(intensities[index][173]);
    R4->InsertNextValue(intensities[index][174]);
    R5->InsertNextValue(intensities[index][175]);
    R6->InsertNextValue(intensities[index][176]);
    R7->InsertNextValue(intensities[index][177]);
    R8->InsertNextValue(intensities[index][178]);
    R9->InsertNextValue(intensities[index][179]);
    S0->InsertNextValue(intensities[index][180]);

    S1->InsertNextValue(intensities[index][181]);
    S2->InsertNextValue(intensities[index][182]);
    S3->InsertNextValue(intensities[index][183]);
    S4->InsertNextValue(intensities[index][184]);
    S5->InsertNextValue(intensities[index][185]);
    S6->InsertNextValue(intensities[index][186]);
    S7->InsertNextValue(intensities[index][187]);
    S8->InsertNextValue(intensities[index][188]);
    S9->InsertNextValue(intensities[index][189]);
    T0->InsertNextValue(intensities[index][190]);

    T1->InsertNextValue(intensities[index][191]);
    T2->InsertNextValue(intensities[index][192]);
    T3->InsertNextValue(intensities[index][193]);
    T4->InsertNextValue(intensities[index][194]);
    T5->InsertNextValue(intensities[index][195]);
    T6->InsertNextValue(intensities[index][196]);
    T7->InsertNextValue(intensities[index][197]);
    T8->InsertNextValue(intensities[index][198]);
    T9->InsertNextValue(intensities[index][199]);
    U0->InsertNextValue(intensities[index][200]);

    U1->InsertNextValue(intensities[index][201]);
    U2->InsertNextValue(intensities[index][202]);
    U3->InsertNextValue(intensities[index][203]);
    U4->InsertNextValue(intensities[index][204]);
    U5->InsertNextValue(intensities[index][205]);
    U6->InsertNextValue(intensities[index][206]);
    U7->InsertNextValue(intensities[index][207]);
    U8->InsertNextValue(intensities[index][208]);
    U9->InsertNextValue(intensities[index][209]);
    V0->InsertNextValue(intensities[index][210]);
    V1->InsertNextValue(intensities[index][211]);
    V2->InsertNextValue(intensities[index][212]);
    V3->InsertNextValue(intensities[index][213]);
    V4->InsertNextValue(intensities[index][214]);
    V5->InsertNextValue(intensities[index][215]);
    V6->InsertNextValue(intensities[index][216]);
    V7->InsertNextValue(intensities[index][217]);
    V8->InsertNextValue(intensities[index][218]);
    V9->InsertNextValue(intensities[index][219]);
    W0->InsertNextValue(intensities[index][220]);

    W1->InsertNextValue(intensities[index][221]);
    W2->InsertNextValue(intensities[index][222]);
    W3->InsertNextValue(intensities[index][223]);
    W4->InsertNextValue(intensities[index][224]);
    W5->InsertNextValue(intensities[index][225]);
    W6->InsertNextValue(intensities[index][226]);
    W7->InsertNextValue(intensities[index][227]);
    W8->InsertNextValue(intensities[index][228]);
    W9->InsertNextValue(intensities[index][229]);
    X0->InsertNextValue(intensities[index][230]);
    X1->InsertNextValue(intensities[index][231]);
    X2->InsertNextValue(intensities[index][232]);
    X3->InsertNextValue(intensities[index][233]);
    X4->InsertNextValue(intensities[index][234]);
    X5->InsertNextValue(intensities[index][235]);
    X6->InsertNextValue(intensities[index][236]);
    X7->InsertNextValue(intensities[index][237]);
    X8->InsertNextValue(intensities[index][238]);
    X9->InsertNextValue(intensities[index][239]);
    Y0->InsertNextValue(intensities[index][240]);
    Y1->InsertNextValue(intensities[index][241]);
    Y2->InsertNextValue(intensities[index][242]);
    Y3->InsertNextValue(intensities[index][243]);
    Y4->InsertNextValue(intensities[index][244]);
    Y5->InsertNextValue(intensities[index][245]);
    Y6->InsertNextValue(intensities[index][246]);
    Y7->InsertNextValue(intensities[index][247]);
    Y8->InsertNextValue(intensities[index][248]);
    Y9->InsertNextValue(intensities[index][249]);

    Z0->InsertNextValue(intensities[index][250]);
    Z1->InsertNextValue(intensities[index][251]);
    Z2->InsertNextValue(intensities[index][252]);
    Z3->InsertNextValue(intensities[index][253]);
    Z4->InsertNextValue(intensities[index][254]);

  }
  datasetTable->AddColumn(A0);
  datasetTable->AddColumn(A1);
  datasetTable->AddColumn(A2);
  datasetTable->AddColumn(A3);
  datasetTable->AddColumn(A4);
  datasetTable->AddColumn(A5);
  datasetTable->AddColumn(A6);
  datasetTable->AddColumn(A7);
  datasetTable->AddColumn(A8);
  datasetTable->AddColumn(A9);
  datasetTable->AddColumn(B0);
  datasetTable->AddColumn(B1);
  datasetTable->AddColumn(B2);
  datasetTable->AddColumn(B3);
  datasetTable->AddColumn(B4);
  datasetTable->AddColumn(B5);
  datasetTable->AddColumn(B6);
  datasetTable->AddColumn(B7);
  datasetTable->AddColumn(B8);
  datasetTable->AddColumn(B9);
  datasetTable->AddColumn(C0);
  datasetTable->AddColumn(C1);
  datasetTable->AddColumn(C2);
  datasetTable->AddColumn(C3);
  datasetTable->AddColumn(C4);
  datasetTable->AddColumn(C5);
  datasetTable->AddColumn(C6);
  datasetTable->AddColumn(C7);
  datasetTable->AddColumn(C8);
  datasetTable->AddColumn(C9);
  datasetTable->AddColumn(D0);
  datasetTable->AddColumn(D1);
  datasetTable->AddColumn(D2);
  datasetTable->AddColumn(D3);
  datasetTable->AddColumn(D4);
  datasetTable->AddColumn(D5);
  datasetTable->AddColumn(D6);
  datasetTable->AddColumn(D7);
  datasetTable->AddColumn(D8);
  datasetTable->AddColumn(D9);
  datasetTable->AddColumn(E0);
  datasetTable->AddColumn(E1);
  datasetTable->AddColumn(E2);
  datasetTable->AddColumn(E3);
  datasetTable->AddColumn(E4);
  datasetTable->AddColumn(E5);
  datasetTable->AddColumn(E6);
  datasetTable->AddColumn(E7);
  datasetTable->AddColumn(E8);
  datasetTable->AddColumn(E9);
  datasetTable->AddColumn(F0);
  datasetTable->AddColumn(F1);
  datasetTable->AddColumn(F2);
  datasetTable->AddColumn(F3);
  datasetTable->AddColumn(F4);
  datasetTable->AddColumn(F5);
  datasetTable->AddColumn(F6);
  datasetTable->AddColumn(F7);
  datasetTable->AddColumn(F8);
  datasetTable->AddColumn(F9);
  datasetTable->AddColumn(G0);
  datasetTable->AddColumn(G1);
  datasetTable->AddColumn(G2);
  datasetTable->AddColumn(G3);
  datasetTable->AddColumn(G4);
  datasetTable->AddColumn(G5);
  datasetTable->AddColumn(G6);
  datasetTable->AddColumn(G7);
  datasetTable->AddColumn(G8);
  datasetTable->AddColumn(G9);
  datasetTable->AddColumn(H0);
  datasetTable->AddColumn(H1);
  datasetTable->AddColumn(H2);
  datasetTable->AddColumn(H3);
  datasetTable->AddColumn(H4);
  datasetTable->AddColumn(H5);
  datasetTable->AddColumn(H6);
  datasetTable->AddColumn(H7);
  datasetTable->AddColumn(H8);
  datasetTable->AddColumn(H9);
  datasetTable->AddColumn(I0);
  datasetTable->AddColumn(I1);
  datasetTable->AddColumn(I2);
  datasetTable->AddColumn(I3);
  datasetTable->AddColumn(I4);
  datasetTable->AddColumn(I5);
  datasetTable->AddColumn(I6);
  datasetTable->AddColumn(I7);
  datasetTable->AddColumn(I8);
  datasetTable->AddColumn(I9);
  datasetTable->AddColumn(J0);
  datasetTable->AddColumn(J1);
  datasetTable->AddColumn(J2);
  datasetTable->AddColumn(J3);
  datasetTable->AddColumn(J4);
  datasetTable->AddColumn(J5);
  datasetTable->AddColumn(J6);
  datasetTable->AddColumn(J7);
  datasetTable->AddColumn(J8);
  datasetTable->AddColumn(J9);
  datasetTable->AddColumn(K0);
  datasetTable->AddColumn(K1);
  datasetTable->AddColumn(K2);
  datasetTable->AddColumn(K3);
  datasetTable->AddColumn(K4);
  datasetTable->AddColumn(K5);
  datasetTable->AddColumn(K6);
  datasetTable->AddColumn(K7);
  datasetTable->AddColumn(K8);
  datasetTable->AddColumn(K9);
  datasetTable->AddColumn(L0);
  datasetTable->AddColumn(L1);
  datasetTable->AddColumn(L2);
  datasetTable->AddColumn(L3);
  datasetTable->AddColumn(L4);
  datasetTable->AddColumn(L5);
  datasetTable->AddColumn(L6);
  datasetTable->AddColumn(L7);
  datasetTable->AddColumn(L8);
  datasetTable->AddColumn(L9);
  datasetTable->AddColumn(M0);
  datasetTable->AddColumn(M1);
  datasetTable->AddColumn(M2);
  datasetTable->AddColumn(M3);
  datasetTable->AddColumn(M4);
  datasetTable->AddColumn(M5);
  datasetTable->AddColumn(M6);
  datasetTable->AddColumn(M7);
  datasetTable->AddColumn(M8);
  datasetTable->AddColumn(M9);
  datasetTable->AddColumn(N0);
  datasetTable->AddColumn(N1);
  datasetTable->AddColumn(N2);
  datasetTable->AddColumn(N3);
  datasetTable->AddColumn(N4);
  datasetTable->AddColumn(N5);
  datasetTable->AddColumn(N6);
  datasetTable->AddColumn(N7);
  datasetTable->AddColumn(N8);
  datasetTable->AddColumn(N9);
  datasetTable->AddColumn(O0);
  datasetTable->AddColumn(O1);
  datasetTable->AddColumn(O2);
  datasetTable->AddColumn(O3);
  datasetTable->AddColumn(O4);
  datasetTable->AddColumn(O5);
  datasetTable->AddColumn(O6);
  datasetTable->AddColumn(O7);
  datasetTable->AddColumn(O8);
  datasetTable->AddColumn(O9);
  datasetTable->AddColumn(P0);
  datasetTable->AddColumn(P1);
  datasetTable->AddColumn(P2);
  datasetTable->AddColumn(P3);
  datasetTable->AddColumn(P4);
  datasetTable->AddColumn(P5);
  datasetTable->AddColumn(P6);
  datasetTable->AddColumn(P7);
  datasetTable->AddColumn(P8);
  datasetTable->AddColumn(P9);
  datasetTable->AddColumn(Q0);
  datasetTable->AddColumn(Q1);
  datasetTable->AddColumn(Q2);
  datasetTable->AddColumn(Q3);
  datasetTable->AddColumn(Q4);
  datasetTable->AddColumn(Q5);
  datasetTable->AddColumn(Q6);
  datasetTable->AddColumn(Q7);
  datasetTable->AddColumn(Q8);
  datasetTable->AddColumn(Q9);
  datasetTable->AddColumn(R0);
  datasetTable->AddColumn(R1);
  datasetTable->AddColumn(R2);
  datasetTable->AddColumn(R3);
  datasetTable->AddColumn(R4);
  datasetTable->AddColumn(R5);
  datasetTable->AddColumn(R6);
  datasetTable->AddColumn(R7);
  datasetTable->AddColumn(R8);
  datasetTable->AddColumn(R9);
  datasetTable->AddColumn(S0);
  datasetTable->AddColumn(S1);
  datasetTable->AddColumn(S2);
  datasetTable->AddColumn(S3);
  datasetTable->AddColumn(S4);
  datasetTable->AddColumn(S5);
  datasetTable->AddColumn(S6);
  datasetTable->AddColumn(S7);
  datasetTable->AddColumn(S8);
  datasetTable->AddColumn(S9);
  datasetTable->AddColumn(T0);
  datasetTable->AddColumn(T1);
  datasetTable->AddColumn(T2);
  datasetTable->AddColumn(T3);
  datasetTable->AddColumn(T4);
  datasetTable->AddColumn(T5);
  datasetTable->AddColumn(T6);
  datasetTable->AddColumn(T7);
  datasetTable->AddColumn(T8);
  datasetTable->AddColumn(T9);
  datasetTable->AddColumn(U0);
  datasetTable->AddColumn(U1);
  datasetTable->AddColumn(U2);
  datasetTable->AddColumn(U3);
  datasetTable->AddColumn(U4);
  datasetTable->AddColumn(U5);
  datasetTable->AddColumn(U6);
  datasetTable->AddColumn(U7);
  datasetTable->AddColumn(U8);
  datasetTable->AddColumn(U9);
  datasetTable->AddColumn(V0);
  datasetTable->AddColumn(V1);
  datasetTable->AddColumn(V2);
  datasetTable->AddColumn(V3);
  datasetTable->AddColumn(V4);
  datasetTable->AddColumn(V5);
  datasetTable->AddColumn(V6);
  datasetTable->AddColumn(V7);
  datasetTable->AddColumn(V8);
  datasetTable->AddColumn(V9);
  datasetTable->AddColumn(W0);
  datasetTable->AddColumn(W1);
  datasetTable->AddColumn(W2);
  datasetTable->AddColumn(W3);
  datasetTable->AddColumn(W4);
  datasetTable->AddColumn(W5);
  datasetTable->AddColumn(W6);
  datasetTable->AddColumn(W7);
  datasetTable->AddColumn(W8);
  datasetTable->AddColumn(W9);
  datasetTable->AddColumn(X0);
  datasetTable->AddColumn(X1);
  datasetTable->AddColumn(X2);
  datasetTable->AddColumn(X3);
  datasetTable->AddColumn(X4);
  datasetTable->AddColumn(X5);
  datasetTable->AddColumn(X6);
  datasetTable->AddColumn(X7);
  datasetTable->AddColumn(X8);
  datasetTable->AddColumn(X9);
  datasetTable->AddColumn(Y0);
  datasetTable->AddColumn(Y1);
  datasetTable->AddColumn(Y2);
  datasetTable->AddColumn(Y3);
  datasetTable->AddColumn(Y4);
  datasetTable->AddColumn(Y5);
  datasetTable->AddColumn(Y6);
  datasetTable->AddColumn(Y7);
  datasetTable->AddColumn(Y8);
  datasetTable->AddColumn(Y9);

  datasetTable->AddColumn(Z0);
  datasetTable->AddColumn(Z1);
  datasetTable->AddColumn(Z2);
  datasetTable->AddColumn(Z3);
  datasetTable->AddColumn(Z4);

  vtkSmartPointer<vtkPCAStatistics> pcaStatistics = vtkSmartPointer<vtkPCAStatistics>::New();
#if VTK_MAJOR_VERSION <= 5
  pcaStatistics->SetInput(vtkStatisticsAlgorithm::INPUT_DATA, datasetTable);
#else
  pcaStatistics->SetInputData(vtkStatisticsAlgorithm::INPUT_DATA, datasetTable);
#endif

  pcaStatistics->SetColumnStatus("a0", 1);
  pcaStatistics->SetColumnStatus("a1", 1);
  pcaStatistics->SetColumnStatus("a2", 1);
  pcaStatistics->SetColumnStatus("a3", 1);
  pcaStatistics->SetColumnStatus("a4", 1);
  pcaStatistics->SetColumnStatus("a5", 1);
  pcaStatistics->SetColumnStatus("a6", 1);
  pcaStatistics->SetColumnStatus("a7", 1);
  pcaStatistics->SetColumnStatus("a8", 1);
  pcaStatistics->SetColumnStatus("a9", 1);

  pcaStatistics->SetColumnStatus("b0", 1);
  pcaStatistics->SetColumnStatus("b1", 1);
  pcaStatistics->SetColumnStatus("b2", 1);
  pcaStatistics->SetColumnStatus("b3", 1);
  pcaStatistics->SetColumnStatus("b4", 1);
  pcaStatistics->SetColumnStatus("b5", 1);
  pcaStatistics->SetColumnStatus("b6", 1);
  pcaStatistics->SetColumnStatus("b7", 1);
  pcaStatistics->SetColumnStatus("b8", 1);
  pcaStatistics->SetColumnStatus("b9", 1);

  pcaStatistics->SetColumnStatus("c0", 1);
  pcaStatistics->SetColumnStatus("c1", 1);
  pcaStatistics->SetColumnStatus("c2", 1);
  pcaStatistics->SetColumnStatus("c3", 1);
  pcaStatistics->SetColumnStatus("c4", 1);
  pcaStatistics->SetColumnStatus("c5", 1);
  pcaStatistics->SetColumnStatus("c6", 1);
  pcaStatistics->SetColumnStatus("c7", 1);
  pcaStatistics->SetColumnStatus("c8", 1);
  pcaStatistics->SetColumnStatus("c9", 1);

  pcaStatistics->SetColumnStatus("d0", 1);
  pcaStatistics->SetColumnStatus("d1", 1);
  pcaStatistics->SetColumnStatus("d2", 1);
  pcaStatistics->SetColumnStatus("d3", 1);
  pcaStatistics->SetColumnStatus("d4", 1);
  pcaStatistics->SetColumnStatus("d5", 1);
  pcaStatistics->SetColumnStatus("d6", 1);
  pcaStatistics->SetColumnStatus("d7", 1);
  pcaStatistics->SetColumnStatus("d8", 1);
  pcaStatistics->SetColumnStatus("d9", 1);

  pcaStatistics->SetColumnStatus("e0", 1);
  pcaStatistics->SetColumnStatus("e1", 1);
  pcaStatistics->SetColumnStatus("e2", 1);
  pcaStatistics->SetColumnStatus("e3", 1);
  pcaStatistics->SetColumnStatus("e4", 1);
  pcaStatistics->SetColumnStatus("e5", 1);
  pcaStatistics->SetColumnStatus("e6", 1);
  pcaStatistics->SetColumnStatus("e7", 1);
  pcaStatistics->SetColumnStatus("e8", 1);
  pcaStatistics->SetColumnStatus("e9", 1);
  pcaStatistics->SetColumnStatus("f0", 1);
  pcaStatistics->SetColumnStatus("f1", 1);
  pcaStatistics->SetColumnStatus("f2", 1);
  pcaStatistics->SetColumnStatus("f3", 1);
  pcaStatistics->SetColumnStatus("f4", 1);
  pcaStatistics->SetColumnStatus("f5", 1);
  pcaStatistics->SetColumnStatus("f6", 1);
  pcaStatistics->SetColumnStatus("f7", 1);
  pcaStatistics->SetColumnStatus("f8", 1);
  pcaStatistics->SetColumnStatus("f9", 1);
  pcaStatistics->SetColumnStatus("g0", 1);
  pcaStatistics->SetColumnStatus("g1", 1);
  pcaStatistics->SetColumnStatus("g2", 1);
  pcaStatistics->SetColumnStatus("g3", 1);
  pcaStatistics->SetColumnStatus("g4", 1);
  pcaStatistics->SetColumnStatus("g5", 1);
  pcaStatistics->SetColumnStatus("g6", 1);
  pcaStatistics->SetColumnStatus("g7", 1);
  pcaStatistics->SetColumnStatus("g8", 1);
  pcaStatistics->SetColumnStatus("g9", 1);
  pcaStatistics->SetColumnStatus("h0", 1);
  pcaStatistics->SetColumnStatus("h1", 1);
  pcaStatistics->SetColumnStatus("h2", 1);
  pcaStatistics->SetColumnStatus("h3", 1);
  pcaStatistics->SetColumnStatus("h4", 1);
  pcaStatistics->SetColumnStatus("h5", 1);
  pcaStatistics->SetColumnStatus("h6", 1);
  pcaStatistics->SetColumnStatus("h7", 1);
  pcaStatistics->SetColumnStatus("h8", 1);
  pcaStatistics->SetColumnStatus("h9", 1);
  pcaStatistics->SetColumnStatus("i0", 1);
  pcaStatistics->SetColumnStatus("i1", 1);
  pcaStatistics->SetColumnStatus("i2", 1);
  pcaStatistics->SetColumnStatus("i3", 1);
  pcaStatistics->SetColumnStatus("i4", 1);
  pcaStatistics->SetColumnStatus("i5", 1);
  pcaStatistics->SetColumnStatus("i6", 1);
  pcaStatistics->SetColumnStatus("i7", 1);
  pcaStatistics->SetColumnStatus("i8", 1);
  pcaStatistics->SetColumnStatus("i9", 1);
  pcaStatistics->SetColumnStatus("j0", 1);
  pcaStatistics->SetColumnStatus("j1", 1);
  pcaStatistics->SetColumnStatus("j2", 1);
  pcaStatistics->SetColumnStatus("j3", 1);
  pcaStatistics->SetColumnStatus("j4", 1);
  pcaStatistics->SetColumnStatus("j5", 1);
  pcaStatistics->SetColumnStatus("j6", 1);
  pcaStatistics->SetColumnStatus("j7", 1);
  pcaStatistics->SetColumnStatus("j8", 1);
  pcaStatistics->SetColumnStatus("j9", 1);
  pcaStatistics->SetColumnStatus("k0", 1);
  pcaStatistics->SetColumnStatus("k1", 1);
  pcaStatistics->SetColumnStatus("k2", 1);
  pcaStatistics->SetColumnStatus("k3", 1);
  pcaStatistics->SetColumnStatus("k4", 1);
  pcaStatistics->SetColumnStatus("k5", 1);
  pcaStatistics->SetColumnStatus("k6", 1);
  pcaStatistics->SetColumnStatus("k7", 1);
  pcaStatistics->SetColumnStatus("k8", 1);
  pcaStatistics->SetColumnStatus("k9", 1);
  pcaStatistics->SetColumnStatus("l0", 1);
  pcaStatistics->SetColumnStatus("l1", 1);
  pcaStatistics->SetColumnStatus("l2", 1);
  pcaStatistics->SetColumnStatus("l3", 1);
  pcaStatistics->SetColumnStatus("l4", 1);
  pcaStatistics->SetColumnStatus("l5", 1);
  pcaStatistics->SetColumnStatus("l6", 1);
  pcaStatistics->SetColumnStatus("l7", 1);
  pcaStatistics->SetColumnStatus("l8", 1);
  pcaStatistics->SetColumnStatus("l9", 1);

  pcaStatistics->SetColumnStatus("m0", 1);
  pcaStatistics->SetColumnStatus("m1", 1);
  pcaStatistics->SetColumnStatus("m2", 1);
  pcaStatistics->SetColumnStatus("m3", 1);
  pcaStatistics->SetColumnStatus("m4", 1);
  pcaStatistics->SetColumnStatus("m5", 1);
  pcaStatistics->SetColumnStatus("m6", 1);
  pcaStatistics->SetColumnStatus("m7", 1);
  pcaStatistics->SetColumnStatus("m8", 1);
  pcaStatistics->SetColumnStatus("m9", 1);

  pcaStatistics->SetColumnStatus("n0", 1);
  pcaStatistics->SetColumnStatus("n1", 1);
  pcaStatistics->SetColumnStatus("n2", 1);
  pcaStatistics->SetColumnStatus("n3", 1);
  pcaStatistics->SetColumnStatus("n4", 1);
  pcaStatistics->SetColumnStatus("n5", 1);
  pcaStatistics->SetColumnStatus("n6", 1);
  pcaStatistics->SetColumnStatus("n7", 1);
  pcaStatistics->SetColumnStatus("n8", 1);
  pcaStatistics->SetColumnStatus("n9", 1);

  pcaStatistics->SetColumnStatus("o0", 1);
  pcaStatistics->SetColumnStatus("o1", 1);
  pcaStatistics->SetColumnStatus("o2", 1);
  pcaStatistics->SetColumnStatus("o3", 1);
  pcaStatistics->SetColumnStatus("o4", 1);
  pcaStatistics->SetColumnStatus("o5", 1);
  pcaStatistics->SetColumnStatus("o6", 1);
  pcaStatistics->SetColumnStatus("o7", 1);
  pcaStatistics->SetColumnStatus("o8", 1);
  pcaStatistics->SetColumnStatus("o9", 1);

  pcaStatistics->SetColumnStatus("p0", 1);
  pcaStatistics->SetColumnStatus("p1", 1);
  pcaStatistics->SetColumnStatus("p2", 1);
  pcaStatistics->SetColumnStatus("p3", 1);
  pcaStatistics->SetColumnStatus("p4", 1);
  pcaStatistics->SetColumnStatus("p5", 1);
  pcaStatistics->SetColumnStatus("p6", 1);
  pcaStatistics->SetColumnStatus("p7", 1);
  pcaStatistics->SetColumnStatus("p8", 1);
  pcaStatistics->SetColumnStatus("p9", 1);

  pcaStatistics->SetColumnStatus("q0", 1);
  pcaStatistics->SetColumnStatus("q1", 1);
  pcaStatistics->SetColumnStatus("q2", 1);
  pcaStatistics->SetColumnStatus("q3", 1);
  pcaStatistics->SetColumnStatus("q4", 1);
  pcaStatistics->SetColumnStatus("q5", 1);
  pcaStatistics->SetColumnStatus("q6", 1);
  pcaStatistics->SetColumnStatus("q7", 1);
  pcaStatistics->SetColumnStatus("q8", 1);
  pcaStatistics->SetColumnStatus("q9", 1);

  pcaStatistics->SetColumnStatus("r0", 1);
  pcaStatistics->SetColumnStatus("r1", 1);
  pcaStatistics->SetColumnStatus("r2", 1);
  pcaStatistics->SetColumnStatus("r3", 1);
  pcaStatistics->SetColumnStatus("r4", 1);
  pcaStatistics->SetColumnStatus("r5", 1);
  pcaStatistics->SetColumnStatus("r6", 1);
  pcaStatistics->SetColumnStatus("r7", 1);
  pcaStatistics->SetColumnStatus("r8", 1);
  pcaStatistics->SetColumnStatus("r9", 1);

  pcaStatistics->SetColumnStatus("s0", 1);
  pcaStatistics->SetColumnStatus("s1", 1);
  pcaStatistics->SetColumnStatus("s2", 1);
  pcaStatistics->SetColumnStatus("s3", 1);
  pcaStatistics->SetColumnStatus("s4", 1);
  pcaStatistics->SetColumnStatus("s5", 1);
  pcaStatistics->SetColumnStatus("s6", 1);
  pcaStatistics->SetColumnStatus("s7", 1);
  pcaStatistics->SetColumnStatus("s8", 1);
  pcaStatistics->SetColumnStatus("s9", 1);

  pcaStatistics->SetColumnStatus("t0", 1);
  pcaStatistics->SetColumnStatus("t1", 1);
  pcaStatistics->SetColumnStatus("t2", 1);
  pcaStatistics->SetColumnStatus("t3", 1);
  pcaStatistics->SetColumnStatus("t4", 1);
  pcaStatistics->SetColumnStatus("t5", 1);
  pcaStatistics->SetColumnStatus("t6", 1);
  pcaStatistics->SetColumnStatus("t7", 1);
  pcaStatistics->SetColumnStatus("t8", 1);
  pcaStatistics->SetColumnStatus("t9", 1);

  pcaStatistics->SetColumnStatus("u0", 1);
  pcaStatistics->SetColumnStatus("u1", 1);
  pcaStatistics->SetColumnStatus("u2", 1);
  pcaStatistics->SetColumnStatus("u3", 1);
  pcaStatistics->SetColumnStatus("u4", 1);
  pcaStatistics->SetColumnStatus("u5", 1);
  pcaStatistics->SetColumnStatus("u6", 1);
  pcaStatistics->SetColumnStatus("u7", 1);
  pcaStatistics->SetColumnStatus("u8", 1);
  pcaStatistics->SetColumnStatus("u9", 1);

  pcaStatistics->SetColumnStatus("v0", 1);
  pcaStatistics->SetColumnStatus("v1", 1);
  pcaStatistics->SetColumnStatus("v2", 1);
  pcaStatistics->SetColumnStatus("v3", 1);
  pcaStatistics->SetColumnStatus("v4", 1);
  pcaStatistics->SetColumnStatus("v5", 1);
  pcaStatistics->SetColumnStatus("v6", 1);
  pcaStatistics->SetColumnStatus("v7", 1);
  pcaStatistics->SetColumnStatus("v8", 1);
  pcaStatistics->SetColumnStatus("v9", 1);

  pcaStatistics->SetColumnStatus("w0", 1);
  pcaStatistics->SetColumnStatus("w1", 1);
  pcaStatistics->SetColumnStatus("w2", 1);
  pcaStatistics->SetColumnStatus("w3", 1);
  pcaStatistics->SetColumnStatus("w4", 1);
  pcaStatistics->SetColumnStatus("w5", 1);
  pcaStatistics->SetColumnStatus("w6", 1);
  pcaStatistics->SetColumnStatus("w7", 1);
  pcaStatistics->SetColumnStatus("w8", 1);
  pcaStatistics->SetColumnStatus("w9", 1);

  pcaStatistics->SetColumnStatus("x0", 1);
  pcaStatistics->SetColumnStatus("x1", 1);
  pcaStatistics->SetColumnStatus("x2", 1);
  pcaStatistics->SetColumnStatus("x3", 1);
  pcaStatistics->SetColumnStatus("x4", 1);
  pcaStatistics->SetColumnStatus("x5", 1);
  pcaStatistics->SetColumnStatus("x6", 1);
  pcaStatistics->SetColumnStatus("x7", 1);
  pcaStatistics->SetColumnStatus("x8", 1);
  pcaStatistics->SetColumnStatus("x9", 1);

  pcaStatistics->SetColumnStatus("y0", 1);
  pcaStatistics->SetColumnStatus("y1", 1);
  pcaStatistics->SetColumnStatus("y2", 1);
  pcaStatistics->SetColumnStatus("y3", 1);
  pcaStatistics->SetColumnStatus("y4", 1);
  pcaStatistics->SetColumnStatus("y5", 1);
  pcaStatistics->SetColumnStatus("y6", 1);
  pcaStatistics->SetColumnStatus("y7", 1);
  pcaStatistics->SetColumnStatus("y8", 1);
  pcaStatistics->SetColumnStatus("y9", 1);

  pcaStatistics->SetColumnStatus("z0", 1);
  pcaStatistics->SetColumnStatus("z1", 1);
  pcaStatistics->SetColumnStatus("z2", 1);
  pcaStatistics->SetColumnStatus("z3", 1);
  pcaStatistics->SetColumnStatus("z4", 1);

  pcaStatistics->RequestSelectedColumns();
  pcaStatistics->SetDeriveOption(true);
  pcaStatistics->Update();
  std::cout << "PCA update done" << std::endl;

  VariableSizeMatrixType transposePCAMatrix;
  transposePCAMatrix.SetSize(NumberOfFeatures, NumberOfFeatures);


  vtkSmartPointer<vtkTable> projectedDatasetTable = vtkSmartPointer<vtkTable>::New();
  for (vtkIdType r = 0; r < static_cast<vtkIdType>(NumberOfFeatures); r++)
  {
    vtkSmartPointer<vtkDoubleArray> col = vtkSmartPointer<vtkDoubleArray>::New();
    for (vtkIdType c = 0; c < static_cast<vtkIdType>(NumberOfSamples); c++)
      col->InsertNextValue(0);
    projectedDatasetTable->AddColumn(col);
  }
  std::cout << "projected table created" << std::endl;

  try
  {

    vtkSmartPointer<vtkDoubleArray> eigenvectors = vtkSmartPointer<vtkDoubleArray>::New();
    pcaStatistics->GetEigenvectors(eigenvectors);
    std::cout << "eigen vectors retrieved" << std::endl;


    for (vtkIdType i = 0; i < static_cast<vtkIdType>(NumberOfFeatures); i++)
    {
      //std::cout << eigenvectors->GetNumberOfComponents() << std::endl;
      double* evec = new double[eigenvectors->GetNumberOfComponents()];
      eigenvectors->GetTuple(i, evec);
      for (vtkIdType j = 0; j < eigenvectors->GetNumberOfComponents(); j++)
      {
        transposePCAMatrix[i][j] = evec[j];
      }
      delete evec;
    }
    std::cout << "eigen vectors retrieved" << std::endl;

    PCATransformationMatrix = this->MatrixTranspose(transposePCAMatrix);
    std::cout << "Transpose done" << std::endl;

    for (size_t i = 0; i < NumberOfSamples; i++)
    {
      for (size_t j = 0; j < NumberOfFeatures; j++)
      {
        double sum = 0;
        for (size_t k = 0; k < NumberOfFeatures; k++)
          sum = sum + datasetTable->GetValue(i, k).ToDouble()*PCATransformationMatrix[k][j];
        projectedDatasetTable->SetValue(i, j, vtkVariant(sum));
      }
    }
    std::cout << "projected table populated" << std::endl;
    VariableLengthVectorType mMeanVector = this->ComputeMeanOfGivenFeatureVectors(projectedDatasetTable);
    for (size_t c = 0; c < NumberOfFeatures; c++)
      for (size_t r = 0; r < NumberOfSamples; r++)
        projectedDatasetTable->SetValue(r, c, projectedDatasetTable->GetValue(r, c).ToDouble() - mMeanVector[c]);
    std::cout << "projected table revised" << std::endl;

    TransformationMatrix = PCATransformationMatrix;
    MeanVector = mPMeanvector;
  }
  catch (const std::exception& e1)
  {
  }

  return projectedDatasetTable;
}


//
//vtkSmartPointer< vtkTable >  FeatureReductionClass::GetDiscerningPerfusionTimePointsFullPCA(VectorVectorDouble &intensities, VariableSizeMatrixType &TransformationMatrix, VariableLengthVectorType &MeanVector)
//{
//  mPMeanvector = ComputeMeanOfGivenFeatureVectors(intensities);
//  size_t NumberOfFeatures = intensities[0].size();
//  size_t NumberOfSamples = intensities.size();
//
//  vtkSmartPointer<vtkTable> datasetTable = vtkSmartPointer<vtkTable>::New();
//  vtkSmartPointer<vtkDoubleArray> A0 = vtkSmartPointer<vtkDoubleArray>::New();
//  vtkSmartPointer<vtkDoubleArray> A1 = vtkSmartPointer<vtkDoubleArray>::New();
//  vtkSmartPointer<vtkDoubleArray> A2 = vtkSmartPointer<vtkDoubleArray>::New();
//  vtkSmartPointer<vtkDoubleArray> A3 = vtkSmartPointer<vtkDoubleArray>::New();
//  vtkSmartPointer<vtkDoubleArray> A4 = vtkSmartPointer<vtkDoubleArray>::New();
//  vtkSmartPointer<vtkDoubleArray> A5 = vtkSmartPointer<vtkDoubleArray>::New();
//  vtkSmartPointer<vtkDoubleArray> A6 = vtkSmartPointer<vtkDoubleArray>::New();
//  vtkSmartPointer<vtkDoubleArray> A7 = vtkSmartPointer<vtkDoubleArray>::New();
//  vtkSmartPointer<vtkDoubleArray> A8 = vtkSmartPointer<vtkDoubleArray>::New();
//  vtkSmartPointer<vtkDoubleArray> A9 = vtkSmartPointer<vtkDoubleArray>::New();
//  vtkSmartPointer<vtkDoubleArray> B0 = vtkSmartPointer<vtkDoubleArray>::New();
//  vtkSmartPointer<vtkDoubleArray> B1 = vtkSmartPointer<vtkDoubleArray>::New();
//  vtkSmartPointer<vtkDoubleArray> B2 = vtkSmartPointer<vtkDoubleArray>::New();
//  vtkSmartPointer<vtkDoubleArray> B3 = vtkSmartPointer<vtkDoubleArray>::New();
//  vtkSmartPointer<vtkDoubleArray> B4 = vtkSmartPointer<vtkDoubleArray>::New();
//  vtkSmartPointer<vtkDoubleArray> B5 = vtkSmartPointer<vtkDoubleArray>::New();
//  vtkSmartPointer<vtkDoubleArray> B6 = vtkSmartPointer<vtkDoubleArray>::New();
//  vtkSmartPointer<vtkDoubleArray> B7 = vtkSmartPointer<vtkDoubleArray>::New();
//  vtkSmartPointer<vtkDoubleArray> B8 = vtkSmartPointer<vtkDoubleArray>::New();
//  vtkSmartPointer<vtkDoubleArray> B9 = vtkSmartPointer<vtkDoubleArray>::New();
//
//  A0->SetNumberOfComponents(1);
//  A1->SetNumberOfComponents(1);
//  A2->SetNumberOfComponents(1);
//  A3->SetNumberOfComponents(1);
//  A4->SetNumberOfComponents(1);
//  A5->SetNumberOfComponents(1);
//  A6->SetNumberOfComponents(1);
//  A7->SetNumberOfComponents(1);
//  A8->SetNumberOfComponents(1);
//  A9->SetNumberOfComponents(1);
//  B0->SetNumberOfComponents(1);
//  B1->SetNumberOfComponents(1);
//  B2->SetNumberOfComponents(1);
//  B3->SetNumberOfComponents(1);
//  B4->SetNumberOfComponents(1);
//  B5->SetNumberOfComponents(1);
//  B6->SetNumberOfComponents(1);
//  B7->SetNumberOfComponents(1);
//  B8->SetNumberOfComponents(1);
//  B9->SetNumberOfComponents(1);
//
//  A0->SetName("a0");
//  A1->SetName("a1");
//  A2->SetName("a2");
//  A3->SetName("a3");
//  A4->SetName("a4");
//  A5->SetName("a5");
//  A6->SetName("a6");
//  A7->SetName("a7");
//  A8->SetName("a8");
//  A9->SetName("a9");
//
//  B0->SetName("b0");
//  B1->SetName("b1");
//  B2->SetName("b2");
//  B3->SetName("b3");
//  B4->SetName("b4");
//  B5->SetName("b5");
//  B6->SetName("b6");
//  B7->SetName("b7");
//  B8->SetName("b8");
//  B9->SetName("b9");
//
//
//  // this should be made parallel
//  for (size_t index = 0; index < NumberOfSamples; index++)
//  {
//    A0->InsertNextValue(intensities[index][0]);
//    A1->InsertNextValue(intensities[index][1]);
//    A2->InsertNextValue(intensities[index][2]);
//    A3->InsertNextValue(intensities[index][3]);
//    A4->InsertNextValue(intensities[index][4]);
//    A5->InsertNextValue(intensities[index][5]);
//    A6->InsertNextValue(intensities[index][6]);
//    A7->InsertNextValue(intensities[index][7]);
//    A8->InsertNextValue(intensities[index][8]);
//    A9->InsertNextValue(intensities[index][9]);
//    B0->InsertNextValue(intensities[index][10]);
//    B1->InsertNextValue(intensities[index][11]);
//    B2->InsertNextValue(intensities[index][12]);
//    B3->InsertNextValue(intensities[index][13]);
//    B4->InsertNextValue(intensities[index][14]);
//    B5->InsertNextValue(intensities[index][15]);
//    B6->InsertNextValue(intensities[index][16]);
//    B7->InsertNextValue(intensities[index][17]);
//    B8->InsertNextValue(intensities[index][18]);
//    B9->InsertNextValue(intensities[index][19]);
//  }
//  datasetTable->AddColumn(A0);
//  datasetTable->AddColumn(A1);
//  datasetTable->AddColumn(A2);
//  datasetTable->AddColumn(A3);
//  datasetTable->AddColumn(A4);
//  datasetTable->AddColumn(A5);
//  datasetTable->AddColumn(A6);
//  datasetTable->AddColumn(A7);
//  datasetTable->AddColumn(A8);
//  datasetTable->AddColumn(A9);
//  datasetTable->AddColumn(B0);
//  datasetTable->AddColumn(B1);
//  datasetTable->AddColumn(B2);
//  datasetTable->AddColumn(B3);
//  datasetTable->AddColumn(B4);
//  datasetTable->AddColumn(B5);
//  datasetTable->AddColumn(B6);
//  datasetTable->AddColumn(B7);
//  datasetTable->AddColumn(B8);
//  datasetTable->AddColumn(B9);
//
//  vtkSmartPointer<vtkPCAStatistics> pcaStatistics = vtkSmartPointer<vtkPCAStatistics>::New();
//#if VTK_MAJOR_VERSION <= 5
//  pcaStatistics->SetInput(vtkStatisticsAlgorithm::INPUT_DATA, datasetTable);
//#else
//  pcaStatistics->SetInputData(vtkStatisticsAlgorithm::INPUT_DATA, datasetTable);
//#endif
//
//  pcaStatistics->SetColumnStatus("a0", 1);
//  pcaStatistics->SetColumnStatus("a1", 1);
//  pcaStatistics->SetColumnStatus("a2", 1);
//  pcaStatistics->SetColumnStatus("a3", 1);
//  pcaStatistics->SetColumnStatus("a4", 1);
//  pcaStatistics->SetColumnStatus("a5", 1);
//  pcaStatistics->SetColumnStatus("a6", 1);
//  pcaStatistics->SetColumnStatus("a7", 1);
//  pcaStatistics->SetColumnStatus("a8", 1);
//  pcaStatistics->SetColumnStatus("a9", 1);
//
//  pcaStatistics->SetColumnStatus("b0", 1);
//  pcaStatistics->SetColumnStatus("b1", 1);
//  pcaStatistics->SetColumnStatus("b2", 1);
//  pcaStatistics->SetColumnStatus("b3", 1);
//  pcaStatistics->SetColumnStatus("b4", 1);
//  pcaStatistics->SetColumnStatus("b5", 1);
//  pcaStatistics->SetColumnStatus("b6", 1);
//  pcaStatistics->SetColumnStatus("b7", 1);
//  pcaStatistics->SetColumnStatus("b8", 1);
//  pcaStatistics->SetColumnStatus("b9", 1);
//
//  pcaStatistics->RequestSelectedColumns();
//  pcaStatistics->SetDeriveOption(true);
//  pcaStatistics->Update();
//  std::cout << "PCA update done" << std::endl;
//
//  VariableSizeMatrixType transposePCAMatrix;
//  transposePCAMatrix.SetSize(NumberOfFeatures, NumberOfFeatures);
//
//
//  vtkSmartPointer<vtkTable> projectedDatasetTable = vtkSmartPointer<vtkTable>::New();
//  for (vtkIdType r = 0; r < static_cast<vtkIdType>(NumberOfFeatures); r++)
//  {
//    vtkSmartPointer<vtkDoubleArray> col = vtkSmartPointer<vtkDoubleArray>::New();
//    for (vtkIdType c = 0; c < static_cast<vtkIdType>(NumberOfSamples); c++)
//      col->InsertNextValue(0);
//    projectedDatasetTable->AddColumn(col);
//  }
//  std::cout << "projected table created" << std::endl;
//
//  try
//  {
//
//    vtkSmartPointer<vtkDoubleArray> eigenvectors = vtkSmartPointer<vtkDoubleArray>::New();
//    pcaStatistics->GetEigenvectors(eigenvectors);
//    std::cout << "eigen vectors retrieved" << std::endl;
//
//
//    for (vtkIdType i = 0; i < static_cast<vtkIdType>(NumberOfFeatures); i++)
//    {
//      std::cout << eigenvectors->GetNumberOfComponents() << std::endl;
//      double* evec = new double[eigenvectors->GetNumberOfComponents()];
//      eigenvectors->GetTuple(i, evec);
//      std::cout << eigenvectors->GetNumberOfComponents() << std::endl;
//      for (vtkIdType j = 0; j < eigenvectors->GetNumberOfComponents(); j++)
//      {
//         transposePCAMatrix[i][j] = evec[j];
//      }
//      delete evec;
//    }
//    std::cout << "eigen vectors retrieved" << std::endl;
//
//    PCATransformationMatrix = this->MatrixTranspose(transposePCAMatrix);
//    std::cout << "Transpose done" << std::endl;
//
//    for (size_t i = 0; i < NumberOfSamples; i++)
//    {
//      for (size_t j = 0; j < NumberOfFeatures; j++)
//      {
//        double sum = 0;
//        for (size_t k = 0; k < NumberOfFeatures; k++)
//          sum = sum + datasetTable->GetValue(i, k).ToDouble()*PCATransformationMatrix[k][j];
//        projectedDatasetTable->SetValue(i, j, vtkVariant(sum));
//      }
//    }
//    std::cout << "projected table populated" << std::endl;
//    VariableLengthVectorType mMeanVector = this->ComputeMeanOfGivenFeatureVectors(projectedDatasetTable);
//    for (size_t c = 0; c < NumberOfFeatures; c++)
//      for (size_t r = 0; r < NumberOfSamples; r++)
//        projectedDatasetTable->SetValue(r, c, projectedDatasetTable->GetValue(r, c).ToDouble() - mMeanVector[c]);
//    std::cout << "projected table revised" << std::endl;
//
//    TransformationMatrix = PCATransformationMatrix;
//    MeanVector = mPMeanvector;
//  }
//  catch (const std::exception& e1)
//  {
//    std::cerr << e1.what() << "\n";
//  }
//  return projectedDatasetTable;
//}

VectorVectorDouble FeatureReductionClass::ApplyPCAOnTestDataWithGivenTransformations(VectorVectorDouble &intensities, VariableSizeMatrixType & TransformationMatrix, VariableLengthVectorType & MeanVector)
{
  size_t NumberOfFeatures = intensities[0].size();
  size_t NumberOfSamples = intensities.size();

  for (size_t index = 0; index < NumberOfSamples; index++)
  {
    for (size_t timepoint = 0; timepoint < NumberOfFeatures; timepoint++)
      intensities[index][timepoint] = intensities[index][timepoint] - MeanVector[timepoint];
  }
  vtkSmartPointer<vtkTable> projectedDatasetTable = vtkSmartPointer<vtkTable>::New();
  for (vtkIdType r = 0; r < static_cast<vtkIdType>(NumberOfFeatures); r++)
  {
    vtkSmartPointer<vtkDoubleArray> col = vtkSmartPointer<vtkDoubleArray>::New();
    for (vtkIdType c = 0; c < static_cast<vtkIdType>(NumberOfSamples); c++)
      col->InsertNextValue(0);
    projectedDatasetTable->AddColumn(col);
  }
  VectorVectorDouble projectedData;

  for (size_t i = 0; i < NumberOfSamples; i++)
  {
    VectorDouble oneDataPoint;
    for (size_t j = 0; j < NumberOfFeatures; j++)
    {
      double sum = 0;
      for (size_t k = 0; k < NumberOfFeatures; k++)
        sum = sum + intensities[i][k] * TransformationMatrix[k][j];
      oneDataPoint.push_back(sum);
    }
    projectedData.push_back(oneDataPoint);
  }
  return projectedData;
}


vtkSmartPointer< vtkTable >  FeatureReductionClass::GetDiscerningPerfusionTimePoints(VectorVectorDouble &intensities, VariableSizeMatrixType &TransformationMatrix, VariableLengthVectorType &MeanVector)
{
  mPMeanvector = ComputeMeanOfGivenFeatureVectors(intensities);
  size_t NumberOfFeatures = intensities[0].size();
  size_t NumberOfSamples = intensities.size();

  vtkSmartPointer<vtkTable> datasetTable = vtkSmartPointer<vtkTable>::New();
  vtkSmartPointer<vtkDoubleArray> A0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> A1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> B0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> B1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> C0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> C1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> D0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> D1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> E0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> E1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> F0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> F1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> G0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> G1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> H0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> H1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> I0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> I1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> J0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> J1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> K0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> K1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> L0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> L1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> M0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> M1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> N0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> N1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> O0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> O1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> P0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> P1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Q0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Q1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> R0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> R1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> S0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> S1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> T0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> T1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> U0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> U1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> V0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> V1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> W0 = vtkSmartPointer<vtkDoubleArray>::New();

  A0->SetNumberOfComponents(1);
  A1->SetNumberOfComponents(1);
  B0->SetNumberOfComponents(1);
  B1->SetNumberOfComponents(1);
  C0->SetNumberOfComponents(1);
  C1->SetNumberOfComponents(1);
  D0->SetNumberOfComponents(1);
  D1->SetNumberOfComponents(1);
  E0->SetNumberOfComponents(1);
  E1->SetNumberOfComponents(1);
  F0->SetNumberOfComponents(1);
  F1->SetNumberOfComponents(1);
  G0->SetNumberOfComponents(1);
  G1->SetNumberOfComponents(1);
  H0->SetNumberOfComponents(1);
  H1->SetNumberOfComponents(1);
  I0->SetNumberOfComponents(1);
  I1->SetNumberOfComponents(1);
  J0->SetNumberOfComponents(1);
  J1->SetNumberOfComponents(1);
  K0->SetNumberOfComponents(1);
  K1->SetNumberOfComponents(1);
  L0->SetNumberOfComponents(1);
  L1->SetNumberOfComponents(1);
  M0->SetNumberOfComponents(1);
  M1->SetNumberOfComponents(1);
  N0->SetNumberOfComponents(1);
  N1->SetNumberOfComponents(1);
  O0->SetNumberOfComponents(1);
  O1->SetNumberOfComponents(1);
  P0->SetNumberOfComponents(1);
  P1->SetNumberOfComponents(1);
  Q0->SetNumberOfComponents(1);
  Q1->SetNumberOfComponents(1);
  R0->SetNumberOfComponents(1);
  R1->SetNumberOfComponents(1);
  S0->SetNumberOfComponents(1);
  S1->SetNumberOfComponents(1);
  T0->SetNumberOfComponents(1);
  T1->SetNumberOfComponents(1);
  U0->SetNumberOfComponents(1);
  U1->SetNumberOfComponents(1);
  V0->SetNumberOfComponents(1);
  V1->SetNumberOfComponents(1);
  W0->SetNumberOfComponents(1);

  A0->SetName("a0");
  A1->SetName("a1");
  B0->SetName("b0");
  B1->SetName("b1");
  C0->SetName("c0");
  C1->SetName("c1");
  D0->SetName("d0");
  D1->SetName("d1");
  E0->SetName("e0");
  E1->SetName("e1");
  F0->SetName("f0");
  F1->SetName("f1");
  G0->SetName("g0");
  G1->SetName("g1");
  H0->SetName("h0");
  H1->SetName("h1");
  I0->SetName("i0");
  I1->SetName("i1");
  J0->SetName("j0");
  J1->SetName("j1");
  K0->SetName("k0");
  K1->SetName("k1");
  L0->SetName("l0");
  L1->SetName("l1");
  M0->SetName("m0");
  M1->SetName("m1");
  N0->SetName("n0");
  N1->SetName("n1");
  O0->SetName("o0");
  O1->SetName("o1");
  P0->SetName("p0");
  P1->SetName("p1");
  Q0->SetName("q0");
  Q1->SetName("q1");
  R0->SetName("r0");
  R1->SetName("r1");
  S0->SetName("s0");
  S1->SetName("s1");
  T0->SetName("t0");
  T1->SetName("t1");
  U0->SetName("u0");
  U1->SetName("u1");
  V0->SetName("v0");
  V1->SetName("v1");
  W0->SetName("w0");

  // this should be made parallel
  for (size_t index = 0; index < NumberOfSamples; index++)
  {
    A0->InsertNextValue(intensities[index][0]);
    A1->InsertNextValue(intensities[index][1]);
    B0->InsertNextValue(intensities[index][2]);
    B1->InsertNextValue(intensities[index][3]);
    C0->InsertNextValue(intensities[index][4]);
    C1->InsertNextValue(intensities[index][5]);
    D0->InsertNextValue(intensities[index][6]);
    D1->InsertNextValue(intensities[index][7]);
    E0->InsertNextValue(intensities[index][8]);
    E1->InsertNextValue(intensities[index][9]);
    F0->InsertNextValue(intensities[index][10]);
    F1->InsertNextValue(intensities[index][11]);
    G0->InsertNextValue(intensities[index][12]);
    G1->InsertNextValue(intensities[index][13]);
    H0->InsertNextValue(intensities[index][14]);
    H1->InsertNextValue(intensities[index][15]);
    I0->InsertNextValue(intensities[index][16]);
    I1->InsertNextValue(intensities[index][17]);
    J0->InsertNextValue(intensities[index][18]);
    J1->InsertNextValue(intensities[index][19]);
    K0->InsertNextValue(intensities[index][20]);
    K1->InsertNextValue(intensities[index][21]);
    L0->InsertNextValue(intensities[index][22]);
    L1->InsertNextValue(intensities[index][23]);
    M0->InsertNextValue(intensities[index][24]);
    M1->InsertNextValue(intensities[index][25]);
    N0->InsertNextValue(intensities[index][26]);
    N1->InsertNextValue(intensities[index][27]);
    O0->InsertNextValue(intensities[index][28]);
    O1->InsertNextValue(intensities[index][29]);
    P0->InsertNextValue(intensities[index][30]);
    P1->InsertNextValue(intensities[index][31]);
    Q0->InsertNextValue(intensities[index][32]);
    Q1->InsertNextValue(intensities[index][33]);
    R0->InsertNextValue(intensities[index][34]);
    R1->InsertNextValue(intensities[index][35]);
    S0->InsertNextValue(intensities[index][36]);
    S1->InsertNextValue(intensities[index][37]);
    T0->InsertNextValue(intensities[index][38]);
    T1->InsertNextValue(intensities[index][39]);
    U0->InsertNextValue(intensities[index][40]);
    U1->InsertNextValue(intensities[index][41]);
    V0->InsertNextValue(intensities[index][42]);
    V1->InsertNextValue(intensities[index][43]);
    W0->InsertNextValue(intensities[index][44]);
  }
  datasetTable->AddColumn(A0);
  datasetTable->AddColumn(A1);
  datasetTable->AddColumn(B0);
  datasetTable->AddColumn(B1);
  datasetTable->AddColumn(C0);
  datasetTable->AddColumn(C1);
  datasetTable->AddColumn(D0);
  datasetTable->AddColumn(D1);
  datasetTable->AddColumn(E0);
  datasetTable->AddColumn(E1);
  datasetTable->AddColumn(F0);
  datasetTable->AddColumn(F1);
  datasetTable->AddColumn(G0);
  datasetTable->AddColumn(G1);
  datasetTable->AddColumn(H0);
  datasetTable->AddColumn(H1);
  datasetTable->AddColumn(I0);
  datasetTable->AddColumn(I1);
  datasetTable->AddColumn(J0);
  datasetTable->AddColumn(J1);
  datasetTable->AddColumn(K0);
  datasetTable->AddColumn(K1);
  datasetTable->AddColumn(L0);
  datasetTable->AddColumn(L1);
  datasetTable->AddColumn(M0);
  datasetTable->AddColumn(M1);
  datasetTable->AddColumn(N0);
  datasetTable->AddColumn(N1);
  datasetTable->AddColumn(O0);
  datasetTable->AddColumn(O1);
  datasetTable->AddColumn(P0);
  datasetTable->AddColumn(P1);
  datasetTable->AddColumn(Q0);
  datasetTable->AddColumn(Q1);
  datasetTable->AddColumn(R0);
  datasetTable->AddColumn(R1);
  datasetTable->AddColumn(S0);
  datasetTable->AddColumn(S1);
  datasetTable->AddColumn(T0);
  datasetTable->AddColumn(T1);
  datasetTable->AddColumn(U0);
  datasetTable->AddColumn(U1);
  datasetTable->AddColumn(V0);
  datasetTable->AddColumn(V1);
  datasetTable->AddColumn(W0);

  vtkSmartPointer<vtkPCAStatistics> pcaStatistics = vtkSmartPointer<vtkPCAStatistics>::New();
#if VTK_MAJOR_VERSION <= 5
  pcaStatistics->SetInput(vtkStatisticsAlgorithm::INPUT_DATA, datasetTable);
#else
  pcaStatistics->SetInputData(vtkStatisticsAlgorithm::INPUT_DATA, datasetTable);
#endif

  pcaStatistics->SetColumnStatus("a0", 1);
  pcaStatistics->SetColumnStatus("a1", 1);
  pcaStatistics->SetColumnStatus("b0", 1);
  pcaStatistics->SetColumnStatus("b1", 1);
  pcaStatistics->SetColumnStatus("c0", 1);
  pcaStatistics->SetColumnStatus("c1", 1);
  pcaStatistics->SetColumnStatus("d0", 1);
  pcaStatistics->SetColumnStatus("d1", 1);
  pcaStatistics->SetColumnStatus("e0", 1);
  pcaStatistics->SetColumnStatus("e1", 1);
  pcaStatistics->SetColumnStatus("f0", 1);
  pcaStatistics->SetColumnStatus("f1", 1);
  pcaStatistics->SetColumnStatus("g0", 1);
  pcaStatistics->SetColumnStatus("g1", 1);
  pcaStatistics->SetColumnStatus("h0", 1);
  pcaStatistics->SetColumnStatus("h1", 1);
  pcaStatistics->SetColumnStatus("i0", 1);
  pcaStatistics->SetColumnStatus("i1", 1);
  pcaStatistics->SetColumnStatus("j0", 1);
  pcaStatistics->SetColumnStatus("j1", 1);
  pcaStatistics->SetColumnStatus("k0", 1);
  pcaStatistics->SetColumnStatus("k1", 1);
  pcaStatistics->SetColumnStatus("l0", 1);
  pcaStatistics->SetColumnStatus("l1", 1);
  pcaStatistics->SetColumnStatus("m0", 1);
  pcaStatistics->SetColumnStatus("m1", 1);
  pcaStatistics->SetColumnStatus("n0", 1);
  pcaStatistics->SetColumnStatus("n1", 1);
  pcaStatistics->SetColumnStatus("o0", 1);
  pcaStatistics->SetColumnStatus("o1", 1);
  pcaStatistics->SetColumnStatus("p0", 1);
  pcaStatistics->SetColumnStatus("p1", 1);
  pcaStatistics->SetColumnStatus("q0", 1);
  pcaStatistics->SetColumnStatus("q1", 1);
  pcaStatistics->SetColumnStatus("r0", 1);
  pcaStatistics->SetColumnStatus("r1", 1);
  pcaStatistics->SetColumnStatus("s0", 1);
  pcaStatistics->SetColumnStatus("s1", 1);
  pcaStatistics->SetColumnStatus("t0", 1);
  pcaStatistics->SetColumnStatus("t1", 1);
  pcaStatistics->SetColumnStatus("u0", 1);
  pcaStatistics->SetColumnStatus("u1", 1);
  pcaStatistics->SetColumnStatus("v0", 1);
  pcaStatistics->SetColumnStatus("v1", 1);
  pcaStatistics->SetColumnStatus("w0", 1);

  pcaStatistics->RequestSelectedColumns();
  pcaStatistics->SetDeriveOption(true);
  pcaStatistics->Update();

  VariableSizeMatrixType transposePCAMatrix;
  transposePCAMatrix.SetSize(NumberOfFeatures, NumberOfFeatures);

  vtkSmartPointer<vtkDoubleArray> eigenvectors = vtkSmartPointer<vtkDoubleArray>::New();
  pcaStatistics->GetEigenvectors(eigenvectors);

  vtkSmartPointer<vtkTable> projectedDatasetTable = vtkSmartPointer<vtkTable>::New();
  for (vtkIdType r = 0; r < static_cast<vtkIdType>(NumberOfFeatures); r++)
  {
    vtkSmartPointer<vtkDoubleArray> col = vtkSmartPointer<vtkDoubleArray>::New();
    for (vtkIdType c = 0; c < static_cast<vtkIdType>(NumberOfSamples); c++)
      col->InsertNextValue(0);
    projectedDatasetTable->AddColumn(col);
  }

  try
  {
    for (vtkIdType i = 0; i < static_cast<vtkIdType>(NumberOfFeatures); i++)
    {
      double* evec = new double[eigenvectors->GetNumberOfComponents()];
      eigenvectors->GetTuple(i, evec);
      for (vtkIdType j = 0; j < eigenvectors->GetNumberOfComponents(); j++)
        transposePCAMatrix[i][j] = evec[j];
      delete evec;
    }

    PCATransformationMatrix = this->MatrixTranspose(transposePCAMatrix);
    for (size_t i = 0; i < NumberOfSamples; i++)
    {
      for (size_t j = 0; j < NumberOfFeatures; j++)
      {
        double sum = 0;
        for (size_t k = 0; k < NumberOfFeatures; k++)
          sum = sum + datasetTable->GetValue(i, k).ToDouble()*PCATransformationMatrix[k][j];
        projectedDatasetTable->SetValue(i, j, vtkVariant(sum));
      }
    }
    VariableLengthVectorType mMeanVector = this->ComputeMeanOfGivenFeatureVectors(projectedDatasetTable);
    for (size_t c = 0; c < NumberOfFeatures; c++)
      for (size_t r = 0; r < NumberOfSamples; r++)
        projectedDatasetTable->SetValue(r, c, projectedDatasetTable->GetValue(r, c).ToDouble() - mMeanVector[c]);

    TransformationMatrix = PCATransformationMatrix;
    MeanVector = mPMeanvector;
  }
  catch (const std::exception& e1)
  {
    std::cerr << e1.what() << "\n";
  }

  return projectedDatasetTable;
}


vtkSmartPointer< vtkTable >  FeatureReductionClass::GetDiscerningPerfusionTimePointsForPSU(VectorVectorDouble &intensities, VariableSizeMatrixType &TransformationMatrix, VariableLengthVectorType &MeanVector)
{
  mPMeanvector = ComputeMeanOfGivenFeatureVectors(intensities);
  size_t NumberOfFeatures = intensities[0].size();
  size_t NumberOfSamples = intensities.size();

  vtkSmartPointer<vtkTable> datasetTable = vtkSmartPointer<vtkTable>::New();
  vtkSmartPointer<vtkDoubleArray> A0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> A1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> B0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> B1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> C0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> C1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> D0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> D1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> E0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> E1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> F0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> F1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> G0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> G1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> H0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> H1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> I0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> I1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> J0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> J1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> K0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> K1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> L0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> L1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> M0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> M1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> N0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> N1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> O0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> O1 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> P0 = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> P1 = vtkSmartPointer<vtkDoubleArray>::New();
  //vtkSmartPointer<vtkDoubleArray> Q0 = vtkSmartPointer<vtkDoubleArray>::New();
  //vtkSmartPointer<vtkDoubleArray> Q1 = vtkSmartPointer<vtkDoubleArray>::New();
  //vtkSmartPointer<vtkDoubleArray> R0 = vtkSmartPointer<vtkDoubleArray>::New();
  //vtkSmartPointer<vtkDoubleArray> R1 = vtkSmartPointer<vtkDoubleArray>::New();
  //vtkSmartPointer<vtkDoubleArray> S0 = vtkSmartPointer<vtkDoubleArray>::New();
  //vtkSmartPointer<vtkDoubleArray> S1 = vtkSmartPointer<vtkDoubleArray>::New();
  //vtkSmartPointer<vtkDoubleArray> T0 = vtkSmartPointer<vtkDoubleArray>::New();
  //vtkSmartPointer<vtkDoubleArray> T1 = vtkSmartPointer<vtkDoubleArray>::New();
  //vtkSmartPointer<vtkDoubleArray> U0 = vtkSmartPointer<vtkDoubleArray>::New();
  //vtkSmartPointer<vtkDoubleArray> U1 = vtkSmartPointer<vtkDoubleArray>::New();
  //vtkSmartPointer<vtkDoubleArray> V0 = vtkSmartPointer<vtkDoubleArray>::New();
  //vtkSmartPointer<vtkDoubleArray> V1 = vtkSmartPointer<vtkDoubleArray>::New();
  //vtkSmartPointer<vtkDoubleArray> W0 = vtkSmartPointer<vtkDoubleArray>::New();

  A0->SetNumberOfComponents(1);
  A1->SetNumberOfComponents(1);
  B0->SetNumberOfComponents(1);
  B1->SetNumberOfComponents(1);
  C0->SetNumberOfComponents(1);
  C1->SetNumberOfComponents(1);
  D0->SetNumberOfComponents(1);
  D1->SetNumberOfComponents(1);
  E0->SetNumberOfComponents(1);
  E1->SetNumberOfComponents(1);
  F0->SetNumberOfComponents(1);
  F1->SetNumberOfComponents(1);
  G0->SetNumberOfComponents(1);
  G1->SetNumberOfComponents(1);
  H0->SetNumberOfComponents(1);
  H1->SetNumberOfComponents(1);
  I0->SetNumberOfComponents(1);
  I1->SetNumberOfComponents(1);
  J0->SetNumberOfComponents(1);
  J1->SetNumberOfComponents(1);
  K0->SetNumberOfComponents(1);
  K1->SetNumberOfComponents(1);
  L0->SetNumberOfComponents(1);
  L1->SetNumberOfComponents(1);
  M0->SetNumberOfComponents(1);
  M1->SetNumberOfComponents(1);
  N0->SetNumberOfComponents(1);
  N1->SetNumberOfComponents(1);
  O0->SetNumberOfComponents(1);
  O1->SetNumberOfComponents(1);
  P0->SetNumberOfComponents(1);
  P1->SetNumberOfComponents(1);
  //Q0->SetNumberOfComponents(1);
  //Q1->SetNumberOfComponents(1);
  //R0->SetNumberOfComponents(1);
  //R1->SetNumberOfComponents(1);
  //S0->SetNumberOfComponents(1);
  //S1->SetNumberOfComponents(1);
  //T0->SetNumberOfComponents(1);
  //T1->SetNumberOfComponents(1);
  //U0->SetNumberOfComponents(1);
  //U1->SetNumberOfComponents(1);
  //V0->SetNumberOfComponents(1);
  //V1->SetNumberOfComponents(1);
  //W0->SetNumberOfComponents(1);

  A0->SetName("a0");
  A1->SetName("a1");
  B0->SetName("b0");
  B1->SetName("b1");
  C0->SetName("c0");
  C1->SetName("c1");
  D0->SetName("d0");
  D1->SetName("d1");
  E0->SetName("e0");
  E1->SetName("e1");
  F0->SetName("f0");
  F1->SetName("f1");
  G0->SetName("g0");
  G1->SetName("g1");
  H0->SetName("h0");
  H1->SetName("h1");
  I0->SetName("i0");
  I1->SetName("i1");
  J0->SetName("j0");
  J1->SetName("j1");
  K0->SetName("k0");
  K1->SetName("k1");
  L0->SetName("l0");
  L1->SetName("l1");
  M0->SetName("m0");
  M1->SetName("m1");
  N0->SetName("n0");
  N1->SetName("n1");
  O0->SetName("o0");
  O1->SetName("o1");
  P0->SetName("p0");
  P1->SetName("p1");
  //Q0->SetName("q0");
  //Q1->SetName("q1");
  //R0->SetName("r0");
  //R1->SetName("r1");
  //S0->SetName("s0");
  //S1->SetName("s1");
  //T0->SetName("t0");
  //T1->SetName("t1");
  //U0->SetName("u0");
  //U1->SetName("u1");
  //V0->SetName("v0");
  //V1->SetName("v1");
  //W0->SetName("w0");

  // this should be made parallel
  for (size_t index = 0; index < NumberOfSamples; index++)
  {
    A0->InsertNextValue(intensities[index][0]);
    A1->InsertNextValue(intensities[index][1]);
    B0->InsertNextValue(intensities[index][2]);
    B1->InsertNextValue(intensities[index][3]);
    C0->InsertNextValue(intensities[index][4]);
    C1->InsertNextValue(intensities[index][5]);
    D0->InsertNextValue(intensities[index][6]);
    D1->InsertNextValue(intensities[index][7]);
    E0->InsertNextValue(intensities[index][8]);
    E1->InsertNextValue(intensities[index][9]);
    F0->InsertNextValue(intensities[index][10]);
    F1->InsertNextValue(intensities[index][11]);
    G0->InsertNextValue(intensities[index][12]);
    G1->InsertNextValue(intensities[index][13]);
    H0->InsertNextValue(intensities[index][14]);
    H1->InsertNextValue(intensities[index][15]);
    I0->InsertNextValue(intensities[index][16]);
    I1->InsertNextValue(intensities[index][17]);
    J0->InsertNextValue(intensities[index][18]);
    J1->InsertNextValue(intensities[index][19]);
    K0->InsertNextValue(intensities[index][20]);
    K1->InsertNextValue(intensities[index][21]);
    L0->InsertNextValue(intensities[index][22]);
    L1->InsertNextValue(intensities[index][23]);
    M0->InsertNextValue(intensities[index][24]);
    M1->InsertNextValue(intensities[index][25]);
    N0->InsertNextValue(intensities[index][26]);
    N1->InsertNextValue(intensities[index][27]);
    O0->InsertNextValue(intensities[index][28]);
    O1->InsertNextValue(intensities[index][29]);
    P0->InsertNextValue(intensities[index][30]);
    P1->InsertNextValue(intensities[index][31]);
    //Q0->InsertNextValue(intensities[index][32]);
    //Q1->InsertNextValue(intensities[index][33]);
    //R0->InsertNextValue(intensities[index][34]);
    //R1->InsertNextValue(intensities[index][35]);
    //S0->InsertNextValue(intensities[index][36]);
    //S1->InsertNextValue(intensities[index][37]);
    //T0->InsertNextValue(intensities[index][38]);
    //T1->InsertNextValue(intensities[index][39]);
    //U0->InsertNextValue(intensities[index][40]);
    //U1->InsertNextValue(intensities[index][41]);
    //V0->InsertNextValue(intensities[index][42]);
    //V1->InsertNextValue(intensities[index][43]);
    //W0->InsertNextValue(intensities[index][44]);
  }
  datasetTable->AddColumn(A0);
  datasetTable->AddColumn(A1);
  datasetTable->AddColumn(B0);
  datasetTable->AddColumn(B1);
  datasetTable->AddColumn(C0);
  datasetTable->AddColumn(C1);
  datasetTable->AddColumn(D0);
  datasetTable->AddColumn(D1);
  datasetTable->AddColumn(E0);
  datasetTable->AddColumn(E1);
  datasetTable->AddColumn(F0);
  datasetTable->AddColumn(F1);
  datasetTable->AddColumn(G0);
  datasetTable->AddColumn(G1);
  datasetTable->AddColumn(H0);
  datasetTable->AddColumn(H1);
  datasetTable->AddColumn(I0);
  datasetTable->AddColumn(I1);
  datasetTable->AddColumn(J0);
  datasetTable->AddColumn(J1);
  datasetTable->AddColumn(K0);
  datasetTable->AddColumn(K1);
  datasetTable->AddColumn(L0);
  datasetTable->AddColumn(L1);
  datasetTable->AddColumn(M0);
  datasetTable->AddColumn(M1);
  datasetTable->AddColumn(N0);
  datasetTable->AddColumn(N1);
  datasetTable->AddColumn(O0);
  datasetTable->AddColumn(O1);
  datasetTable->AddColumn(P0);
  datasetTable->AddColumn(P1);
  //datasetTable->AddColumn(Q0);
  //datasetTable->AddColumn(Q1);
  //datasetTable->AddColumn(R0);
  //datasetTable->AddColumn(R1);
  //datasetTable->AddColumn(S0);
  //datasetTable->AddColumn(S1);
  //datasetTable->AddColumn(T0);
  //datasetTable->AddColumn(T1);
  //datasetTable->AddColumn(U0);
  //datasetTable->AddColumn(U1);
  //datasetTable->AddColumn(V0);
  //datasetTable->AddColumn(V1);
  //datasetTable->AddColumn(W0);

  vtkSmartPointer<vtkPCAStatistics> pcaStatistics = vtkSmartPointer<vtkPCAStatistics>::New();
#if VTK_MAJOR_VERSION <= 5
  pcaStatistics->SetInput(vtkStatisticsAlgorithm::INPUT_DATA, datasetTable);
#else
  pcaStatistics->SetInputData(vtkStatisticsAlgorithm::INPUT_DATA, datasetTable);
#endif

  pcaStatistics->SetColumnStatus("a0", 1);
  pcaStatistics->SetColumnStatus("a1", 1);
  pcaStatistics->SetColumnStatus("b0", 1);
  pcaStatistics->SetColumnStatus("b1", 1);
  pcaStatistics->SetColumnStatus("c0", 1);
  pcaStatistics->SetColumnStatus("c1", 1);
  pcaStatistics->SetColumnStatus("d0", 1);
  pcaStatistics->SetColumnStatus("d1", 1);
  pcaStatistics->SetColumnStatus("e0", 1);
  pcaStatistics->SetColumnStatus("e1", 1);
  pcaStatistics->SetColumnStatus("f0", 1);
  pcaStatistics->SetColumnStatus("f1", 1);
  pcaStatistics->SetColumnStatus("g0", 1);
  pcaStatistics->SetColumnStatus("g1", 1);
  pcaStatistics->SetColumnStatus("h0", 1);
  pcaStatistics->SetColumnStatus("h1", 1);
  pcaStatistics->SetColumnStatus("i0", 1);
  pcaStatistics->SetColumnStatus("i1", 1);
  pcaStatistics->SetColumnStatus("j0", 1);
  pcaStatistics->SetColumnStatus("j1", 1);
  pcaStatistics->SetColumnStatus("k0", 1);
  pcaStatistics->SetColumnStatus("k1", 1);
  pcaStatistics->SetColumnStatus("l0", 1);
  pcaStatistics->SetColumnStatus("l1", 1);
  pcaStatistics->SetColumnStatus("m0", 1);
  pcaStatistics->SetColumnStatus("m1", 1);
  pcaStatistics->SetColumnStatus("n0", 1);
  pcaStatistics->SetColumnStatus("n1", 1);
  pcaStatistics->SetColumnStatus("o0", 1);
  pcaStatistics->SetColumnStatus("o1", 1);
  pcaStatistics->SetColumnStatus("p0", 1);
  pcaStatistics->SetColumnStatus("p1", 1);
  //pcaStatistics->SetColumnStatus("q0", 1);
  //pcaStatistics->SetColumnStatus("q1", 1);
  //pcaStatistics->SetColumnStatus("r0", 1);
  //pcaStatistics->SetColumnStatus("r1", 1);
  //pcaStatistics->SetColumnStatus("s0", 1);
  //pcaStatistics->SetColumnStatus("s1", 1);
  //pcaStatistics->SetColumnStatus("t0", 1);
  //pcaStatistics->SetColumnStatus("t1", 1);
  //pcaStatistics->SetColumnStatus("u0", 1);
  //pcaStatistics->SetColumnStatus("u1", 1);
  //pcaStatistics->SetColumnStatus("v0", 1);
  //pcaStatistics->SetColumnStatus("v1", 1);
  //pcaStatistics->SetColumnStatus("w0", 1);

  pcaStatistics->RequestSelectedColumns();
  pcaStatistics->SetDeriveOption(true);
  pcaStatistics->Update();

  VariableSizeMatrixType transposePCAMatrix;
  transposePCAMatrix.SetSize(NumberOfFeatures, NumberOfFeatures);

  vtkSmartPointer<vtkDoubleArray> eigenvectors = vtkSmartPointer<vtkDoubleArray>::New();
  pcaStatistics->GetEigenvectors(eigenvectors);

  vtkSmartPointer<vtkTable> projectedDatasetTable = vtkSmartPointer<vtkTable>::New();
  for (vtkIdType r = 0; r < static_cast<vtkIdType>(NumberOfFeatures); r++)
  {
    vtkSmartPointer<vtkDoubleArray> col = vtkSmartPointer<vtkDoubleArray>::New();
    for (vtkIdType c = 0; c < static_cast<vtkIdType>(NumberOfSamples); c++)
      col->InsertNextValue(0);
    projectedDatasetTable->AddColumn(col);
  }

  try
  {
    for (vtkIdType i = 0; i < static_cast<vtkIdType>(NumberOfFeatures); i++)
    {
      double* evec = new double[eigenvectors->GetNumberOfComponents()];
      eigenvectors->GetTuple(i, evec);
      for (vtkIdType j = 0; j < eigenvectors->GetNumberOfComponents(); j++)
        transposePCAMatrix[i][j] = evec[j];
      delete evec;
    }

    PCATransformationMatrix = this->MatrixTranspose(transposePCAMatrix);
    for (size_t i = 0; i < NumberOfSamples; i++)
    {
      for (size_t j = 0; j < NumberOfFeatures; j++)
      {
        double sum = 0;
        for (size_t k = 0; k < NumberOfFeatures; k++)
          sum = sum + datasetTable->GetValue(i, k).ToDouble()*PCATransformationMatrix[k][j];
        projectedDatasetTable->SetValue(i, j, vtkVariant(sum));
      }
    }
    VariableLengthVectorType mMeanVector = this->ComputeMeanOfGivenFeatureVectors(projectedDatasetTable);
    for (size_t c = 0; c < NumberOfFeatures; c++)
      for (size_t r = 0; r < NumberOfSamples; r++)
        projectedDatasetTable->SetValue(r, c, projectedDatasetTable->GetValue(r, c).ToDouble() - mMeanVector[c]);

    TransformationMatrix = PCATransformationMatrix;
    MeanVector = mPMeanvector;
  }
  catch (const std::exception& e1)
  {
    std::cerr << e1.what() << "\n";
  }

  return projectedDatasetTable;
}