/**
\file  FeatureReductionClass.cpp

\brief Implementation of the FeatureReductionClass

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or http://www.med.upenn.edu/sbia/software/license.html

*/

#include <vtkVersion.h>
#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPCAStatistics.h"
#include "vtkStringArray.h"
#include "vtkTable.h"
#include "FeatureReductionClass.h"
#include "itkCSVNumericObjectFileWriter.h"
#include "CAPTk.h"

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
