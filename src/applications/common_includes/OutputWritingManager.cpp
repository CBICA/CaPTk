/**
\file  OutputWritingManager.cpp

\brief Implementation of the OutputWritingManager

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

#include "OutputWritingManager.h"
#include "itkCSVArray2DFileReader.h"
//#include "CAPTk.h"
#include <QFileInfo>
#include <QMessageBox>
#include <QString>
#include <QDir>
#include "cbicaLogging.h"
#include "CaPTkDefines.h"

OutputWritingManager::OutputWritingManager()
{
  mOutputDirectoryPath = ".\\";//TBD this not set correctly in reccurance app
  mLastEncounteredError = "";
}

OutputWritingManager::~OutputWritingManager()
{
}
bool OutputWritingManager::SetupOutputFolders()
{
  QFileInfo fi(QFileInfo(QString::fromStdString(mOutputDirectoryPath)).path());
  if (!fi.isWritable())
  {
    QMessageBox::warning(NULL, "Error", "Output directory is read only.", QMessageBox::Ok);
    return false;
  }

  QString t1path = QString::fromStdString(mOutputDirectoryPath + "/T1");
  QString t2path = QString::fromStdString(mOutputDirectoryPath + "/T2");
  QString t1cepath = QString::fromStdString(mOutputDirectoryPath + "/T1CE");
  QString t2flairpath = QString::fromStdString(mOutputDirectoryPath + "/T2Flair");
  QString perfusionpath = QString::fromStdString(mOutputDirectoryPath + "/Perfusion");
  QString dtipath = QString::fromStdString(mOutputDirectoryPath + "/DTI");
  QString recurrencepath = QString::fromStdString(mOutputDirectoryPath + "/RecurrenceOutput");
  QString maskspath = QString::fromStdString(mOutputDirectoryPath + "/Masks");

  QDir t1dir(t1path);
  if (!t1dir.exists())
    t1dir.mkdir(".");

  QDir t2dir(t2path);
  if (!t2dir.exists())
    t2dir.mkdir(".");

  QDir t1cedir(t1cepath);
  if (!t1cedir.exists())
    t1cedir.mkdir(".");

  QDir t2flairdir(t2flairpath);
  if (!t2flairdir.exists())
    t2flairdir.mkdir(".");

  QDir perfusiondir(perfusionpath);
  if (!perfusiondir.exists())
    perfusiondir.mkdir(".");

  QDir dtidir(dtipath);
  if (!dtidir.exists())
    dtidir.mkdir(".");

  QDir recurrencedir(recurrencepath);
  if (!recurrencedir.exists())
    recurrencedir.mkdir(".");

  QDir masksdir(maskspath);
  if (!masksdir.exists())
    masksdir.mkdir(".");

  return true;
}

void OutputWritingManager::SaveModelResults(VariableSizeMatrixType scaledFeatures, VariableLengthVectorType means, VariableLengthVectorType stds, VariableLengthVectorType pmeans, VariableSizeMatrixType pcacoeff,
  bool useConvData, bool useDTIData, bool usePerfData, bool useDistData, const int size)
{
  typedef vnl_matrix<double> MatrixType;
  MatrixType meanVector;
  MatrixType stdVector;
  MatrixType pmeanVector;
  MatrixType pcaVector;


  meanVector.set_size(1, size);
  stdVector.set_size(1, size);
  for (unsigned int j = 0; j < meanVector.size(); j++)
  {
    meanVector(0, j) = means[j];
    stdVector(0, j) = stds[j];
  }
  if (usePerfData)
  {
    pmeanVector.set_size(1, 45);
    pcaVector.set_size(45, 45);

    for (unsigned int j = 0; j < pmeanVector.size(); j++)
      pmeanVector(0, j) = pmeans[j];

    for (unsigned int i = 0; i < pcacoeff.Rows(); i++)
      for (unsigned int j = 0; j < pcacoeff.Cols(); j++)
        pcaVector(i, j) = pcacoeff[i][j];
  }
  //-----------------mean file writer------------------------

  if (size == 1)
  {
    typedef itk::CSVNumericObjectFileWriter<double, 1, 1> WriterType;
    WriterType::Pointer writermean = WriterType::New();
    WriterType::Pointer writerstd = WriterType::New();
    writermean->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Mean.csv");
    writermean->SetInput(&meanVector);
    writerstd->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Std.csv");
    writerstd->SetInput(&stdVector);
    try
    {
      writermean->Write();
      writerstd->Write();
    }
    catch (itk::ExceptionObject &excp)
    {
      cbica::Logging(loggerFile, "Exception detected while trying to write mean and std values: '" + std::string(excp.GetDescription()));
      exit(EXIT_FAILURE);
    }
  }
  else if (size == 2)
  {
    typedef itk::CSVNumericObjectFileWriter<double, 1, 2> WriterType;
    WriterType::Pointer writermean = WriterType::New();
    WriterType::Pointer writerstd = WriterType::New();
	writermean->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Mean.csv");
	writermean->SetInput(&meanVector);
	writerstd->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Std.csv");
	writerstd->SetInput(&stdVector);
    try
    {
      writermean->Write();
      writerstd->Write();
    }
    catch (itk::ExceptionObject & excp)
    {
      cbica::Logging(loggerFile, "Exception detected while trying to write mean and std values: '" + std::string(excp.GetDescription()));
      exit(EXIT_FAILURE);
    }
  }
  else if (size == 3)
  {
    typedef itk::CSVNumericObjectFileWriter<double, 1, 3> WriterType;
    WriterType::Pointer writermean = WriterType::New();
    WriterType::Pointer writerstd = WriterType::New();
	writermean->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Mean.csv");
	writermean->SetInput(&meanVector);
	writerstd->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Std.csv");
	writerstd->SetInput(&stdVector);
    try
    {
      writermean->Write();
      writerstd->Write();
    }
    catch (itk::ExceptionObject & excp)
    {
      cbica::Logging(loggerFile, "Exception detected while trying to write mean and std values: '" + std::string(excp.GetDescription()));
      exit(EXIT_FAILURE);
    }
  }
  else if (size == 4)
  {
    typedef itk::CSVNumericObjectFileWriter<double, 1, 4> WriterType;
    WriterType::Pointer writermean = WriterType::New();
    WriterType::Pointer writerstd = WriterType::New();
	writermean->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Mean.csv");
	writermean->SetInput(&meanVector);
	writerstd->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Std.csv");
	writerstd->SetInput(&stdVector);
    try
    {
      writermean->Write();
      writerstd->Write();
    }
    catch (itk::ExceptionObject & excp)
    {
      cbica::Logging(loggerFile, "Exception detected while trying to write mean and std values: '" + std::string(excp.GetDescription()));
      exit(EXIT_FAILURE);
    }
  }
  else if (size == 5)
  {
    typedef itk::CSVNumericObjectFileWriter<double, 1, 5> WriterType;
    WriterType::Pointer writermean = WriterType::New();
    WriterType::Pointer writerstd = WriterType::New();
	writermean->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Mean.csv");
	writermean->SetInput(&meanVector);
	writerstd->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Std.csv");
	writerstd->SetInput(&stdVector);
    try
    {
      writermean->Write();
      writerstd->Write();
    }
    catch (itk::ExceptionObject & excp)
    {
      cbica::Logging(loggerFile, "Exception detected while trying to write mean and std values: '" + std::string(excp.GetDescription()));
      exit(EXIT_FAILURE);
    }
  }
  else if (size == 6)
  {
    typedef itk::CSVNumericObjectFileWriter<double, 1, 6> WriterType;
    WriterType::Pointer writermean = WriterType::New();
    WriterType::Pointer writerstd = WriterType::New();
	writermean->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Mean.csv");
	writermean->SetInput(&meanVector);
	writerstd->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Std.csv");
	writerstd->SetInput(&stdVector);
    try
    {
      writermean->Write();
      writerstd->Write();
    }
    catch (itk::ExceptionObject & excp)
    {
      cbica::Logging(loggerFile, "Exception detected while trying to write mean and std values: '" + std::string(excp.GetDescription()));
      exit(EXIT_FAILURE);
    }
  }
  else if (size == 7)
  {
    typedef itk::CSVNumericObjectFileWriter<double, 1, 7> WriterType;
    WriterType::Pointer writermean = WriterType::New();
    WriterType::Pointer writerstd = WriterType::New();
	writermean->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Mean.csv");
	writermean->SetInput(&meanVector);
	writerstd->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Std.csv");
	writerstd->SetInput(&stdVector);
    try
    {
      writermean->Write();
      writerstd->Write();
    }
    catch (itk::ExceptionObject & excp)
    {
      cbica::Logging(loggerFile, "Exception detected while trying to write mean and std values: '" + std::string(excp.GetDescription()));
      exit(EXIT_FAILURE);
    }
  }
  else if (size == 8)
  {
    typedef itk::CSVNumericObjectFileWriter<double, 1, 8> WriterType;
    WriterType::Pointer writermean = WriterType::New();
    WriterType::Pointer writerstd = WriterType::New();
	writermean->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Mean.csv");
	writermean->SetInput(&meanVector);
	writerstd->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Std.csv");
	writerstd->SetInput(&stdVector);
    try
    {
      writermean->Write();
      writerstd->Write();
    }
    catch (itk::ExceptionObject & excp)
    {
      cbica::Logging(loggerFile, "Exception detected while trying to write mean and std values: '" + std::string(excp.GetDescription()));
      exit(EXIT_FAILURE);
    }
  }
  else if (size == 9)
  {
    typedef itk::CSVNumericObjectFileWriter<double, 1, 9> WriterType;
    WriterType::Pointer writermean = WriterType::New();
    WriterType::Pointer writerstd = WriterType::New();
	writermean->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Mean.csv");
	writermean->SetInput(&meanVector);
	writerstd->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Std.csv");
	writerstd->SetInput(&stdVector);
    try
    {
      writermean->Write();
      writerstd->Write();
    }
    catch (itk::ExceptionObject & excp)
    {
      cbica::Logging(loggerFile, "Exception detected while trying to write mean and std values: '" + std::string(excp.GetDescription()));
      exit(EXIT_FAILURE);
    }
  }
  else if (size == 10)
  {
    typedef itk::CSVNumericObjectFileWriter<double, 1, 10> WriterType;
    WriterType::Pointer writermean = WriterType::New();
    WriterType::Pointer writerstd = WriterType::New();
	writermean->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Mean.csv");
	writermean->SetInput(&meanVector);
	writerstd->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Std.csv");
	writerstd->SetInput(&stdVector);
    try
    {
      writermean->Write();
      writerstd->Write();
    }
    catch (itk::ExceptionObject & excp)
    {
      cbica::Logging(loggerFile, "Exception detected while trying to write mean and std values: '" + std::string(excp.GetDescription()));
      exit(EXIT_FAILURE);
    }
  }
  else if (size == 11)
  {
    typedef itk::CSVNumericObjectFileWriter<double, 1, 11> WriterType;
    WriterType::Pointer writermean = WriterType::New();
    WriterType::Pointer writerstd = WriterType::New();
	writermean->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Mean.csv");
	writermean->SetInput(&meanVector);
	writerstd->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Std.csv");
	writerstd->SetInput(&stdVector);
    try
    {
      writermean->Write();
      writerstd->Write();
    }
    catch (itk::ExceptionObject & excp)
    {
      cbica::Logging(loggerFile, "Exception detected while trying to write mean and std values: '" + std::string(excp.GetDescription()));
      exit(EXIT_FAILURE);
    }
  }
  else if (size == 12)
  {
    typedef itk::CSVNumericObjectFileWriter<double, 1, 12> WriterType;
    WriterType::Pointer writermean = WriterType::New();
    WriterType::Pointer writerstd = WriterType::New();
	writermean->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Mean.csv");
	writermean->SetInput(&meanVector);
	writerstd->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Std.csv");
	writerstd->SetInput(&stdVector);
    try
    {
      writermean->Write();
      writerstd->Write();
    }
    catch (itk::ExceptionObject & excp)
    {
      cbica::Logging(loggerFile, "Exception detected while trying to write mean and std values: '" + std::string(excp.GetDescription()));
      exit(EXIT_FAILURE);
    }
  }
  else if (size == 13)
  {
    typedef itk::CSVNumericObjectFileWriter<double, 1, 13> WriterType;
    WriterType::Pointer writermean = WriterType::New();
    WriterType::Pointer writerstd = WriterType::New();
	writermean->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Mean.csv");
	writermean->SetInput(&meanVector);
	writerstd->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Std.csv");
	writerstd->SetInput(&stdVector);
    try
    {
      writermean->Write();
      writerstd->Write();
    }
    catch (itk::ExceptionObject & excp)
    {
      cbica::Logging(loggerFile, "Exception detected while trying to write mean and std values: '" + std::string(excp.GetDescription()));
      exit(EXIT_FAILURE);
    }
  }
  else if (size == 14)
  {
    typedef itk::CSVNumericObjectFileWriter<double, 1, 14> WriterType;
    WriterType::Pointer writermean = WriterType::New();
    WriterType::Pointer writerstd = WriterType::New();
	writermean->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Mean.csv");
	writermean->SetInput(&meanVector);
	writerstd->SetFileName(this->mOutputDirectoryPath + "/Recurrence_ZScore_Std.csv");
	writerstd->SetInput(&stdVector);
    try
    {
      writermean->Write();
      writerstd->Write();
    }
    catch (itk::ExceptionObject & excp)
    {
      cbica::Logging(loggerFile, "Exception detected while trying to write mean and std values: '" + std::string(excp.GetDescription()));
      exit(EXIT_FAILURE);
    }
  }
  //-----------------pca file writer------------------------
  if (usePerfData)
  {
    typedef itk::CSVNumericObjectFileWriter<double, 45, 45> PCAWriterType;
    PCAWriterType::Pointer pcawriter = PCAWriterType::New();
    pcawriter->SetFileName(this->mOutputDirectoryPath + "/Recurrence_COEF.csv");
    pcawriter->SetInput(&pcaVector);
    typedef itk::CSVNumericObjectFileWriter<double, 1, 45> PerfWriterType;
    PerfWriterType::Pointer perfwriter = PerfWriterType::New();
    perfwriter->SetFileName(this->mOutputDirectoryPath + "/Recurrence_MR.csv");
    perfwriter->SetInput(&pmeanVector);

    try
    {
      pcawriter->Write();
      perfwriter->Write();
    }
    catch (itk::ExceptionObject & excp)
    {
      cbica::Logging(loggerFile, "Exception detected while trying to write perfusion measures: '" + std::string(excp.GetDescription()));
      exit(EXIT_FAILURE);
    }
  }
  //-----------------modalities file writer------------------------
  MatrixType modalitiesRecord;
  modalitiesRecord.set_size(1, 7);
  if (useConvData)
    modalitiesRecord(0, 0) = 1;
  else
    modalitiesRecord(0, 0) = 0;

  if (useDTIData)
    modalitiesRecord(0, 1) = 1;
  else
    modalitiesRecord(0, 1) = 0;
  if (usePerfData)
    modalitiesRecord(0, 2) = 1;
  else
    modalitiesRecord(0, 2) = 0;
  if (useDistData)
    modalitiesRecord(0, 3) = 1;
  else
    modalitiesRecord(0, 3) = 0;

  typedef itk::CSVNumericObjectFileWriter<double, 1, 7> ModWriterType;
  ModWriterType::Pointer modwriter = ModWriterType::New();
  modwriter->SetFileName(this->mOutputDirectoryPath + "/Recurrence_SVM_Modalities.csv");
  modwriter->SetInput(&modalitiesRecord);
  try
  {
    modwriter->Write();
  }
  catch (itk::ExceptionObject & excp)
  {
    cbica::Logging(loggerFile, "Exception detected while trying to write modality information: '" + std::string(excp.GetDescription()));
    exit(EXIT_FAILURE);
  }
}
vnl_matrix<double> OutputWritingManager::ReadNumberOfModalities(std::string modalitiesfile)
{
  typedef itk::CSVArray2DFileReader<double> ReaderType;
  ReaderType::Pointer readerModalities = ReaderType::New();
  readerModalities->SetFileName(modalitiesfile);
  readerModalities->SetFieldDelimiterCharacter(',');
  readerModalities->HasColumnHeadersOff();
  readerModalities->HasRowHeadersOff();
  readerModalities->Parse();
  vnl_matrix<double> dataMatrixModalities = readerModalities->GetArray2DDataObject()->GetMatrix();
  return dataMatrixModalities;
}
void OutputWritingManager::ReadModelParameters(std::string meanfile, std::string stdfile, std::string pcafile, std::string pmeanfile, VariableLengthVectorType & mean, VariableLengthVectorType & stds, VariableSizeMatrixType & pca, VariableLengthVectorType & pmean)
{
  typedef itk::CSVArray2DFileReader<double> ReaderType;
  ReaderType::Pointer readerMean = ReaderType::New();
  readerMean->SetFileName(meanfile);
  readerMean->SetFieldDelimiterCharacter(',');
  readerMean->HasColumnHeadersOff();
  readerMean->HasRowHeadersOff();
  readerMean->Parse();
  typedef vnl_matrix<double> MatrixType;
  MatrixType dataMatrixMean = readerMean->GetArray2DDataObject()->GetMatrix();

  ReaderType::Pointer readerpMean = ReaderType::New();
  readerpMean->SetFileName(pmeanfile);
  readerpMean->SetFieldDelimiterCharacter(',');
  readerpMean->HasColumnHeadersOff();
  readerpMean->HasRowHeadersOff();
  readerpMean->Parse();
  MatrixType dataMatrixpMean = readerpMean->GetArray2DDataObject()->GetMatrix();


  ReaderType::Pointer readerStd = ReaderType::New();
  readerStd->SetFileName(stdfile);
  readerStd->SetFieldDelimiterCharacter(',');
  readerStd->HasColumnHeadersOff();
  readerStd->HasRowHeadersOff();
  readerStd->Parse();
  MatrixType dataMatrixStd = readerStd->GetArray2DDataObject()->GetMatrix();

  ReaderType::Pointer readerPca = ReaderType::New();
  readerPca->SetFileName(pcafile);
  readerPca->SetFieldDelimiterCharacter(',');
  readerPca->HasColumnHeadersOff();
  readerPca->HasRowHeadersOff();
  readerPca->Parse();
  MatrixType dataMatrixPca = readerPca->GetArray2DDataObject()->GetMatrix();



  mean.SetSize(dataMatrixMean.size());
  stds.SetSize(dataMatrixStd.size());
  pca.SetSize(dataMatrixPca.rows(), dataMatrixPca.cols());
  pmean.SetSize(dataMatrixpMean.size());

  for (unsigned int i = 0; i < dataMatrixMean.cols(); i++)
  {
    mean[i] = dataMatrixMean(0, i);
    stds[i] = dataMatrixStd(0, i);
  }

  for (unsigned int i = 0; i < dataMatrixpMean.cols(); i++)
    pmean[i] = dataMatrixpMean(0, i);

  for (unsigned int i = 0; i < dataMatrixPca.rows(); i++)
    for (unsigned int j = 0; j < dataMatrixPca.cols(); j++)
      pca[i][j] = dataMatrixPca(i, j);
}
void OutputWritingManager::ReadModelParameters(std::string meanfile, std::string stdfile, VariableLengthVectorType & mean, VariableLengthVectorType & stds)
{
  typedef itk::CSVArray2DFileReader<double> ReaderType;
  ReaderType::Pointer readerMean = ReaderType::New();
  readerMean->SetFileName(meanfile);
  readerMean->SetFieldDelimiterCharacter(',');
  readerMean->HasColumnHeadersOff();
  readerMean->HasRowHeadersOff();
  readerMean->Parse();
  typedef vnl_matrix<double> MatrixType;
  MatrixType dataMatrixMean = readerMean->GetArray2DDataObject()->GetMatrix();

  ReaderType::Pointer readerStd = ReaderType::New();
  readerStd->SetFileName(stdfile);
  readerStd->SetFieldDelimiterCharacter(',');
  readerStd->HasColumnHeadersOff();
  readerStd->HasRowHeadersOff();
  readerStd->Parse();
  MatrixType dataMatrixStd = readerStd->GetArray2DDataObject()->GetMatrix();

  mean.SetSize(dataMatrixMean.size());
  stds.SetSize(dataMatrixStd.size());

  for (unsigned int i = 0; i < dataMatrixMean.cols(); i++)
  {
    mean[i] = dataMatrixMean(0, i);
    stds[i] = dataMatrixStd(0, i);
  }
}