/**
\file  OutputWritingManager.h

\brief Declaration of the OutputWritingManager

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/
#pragma once

//#include "CAPTk.h"
#include "itkImage.h"
#include "itkCSVNumericObjectFileWriter.h"
#include "itkArray2D.h"

#include "itkVariableSizeMatrix.h"
#include "itkVariableLengthVector.h"
#include "itkImageRegionIterator.h"
#include "cbicaITKSafeImageIO.h"


using VariableSizeMatrixType = itk::VariableSizeMatrix< double >;
using VariableLengthVectorType = itk::VariableLengthVector< double >;

enum WritingType
{
  SUSAN, BIAS_CORRECT, REGISTRATION, RECURRENCE_MAP
};

class OutputWritingManager
{
public:
  OutputWritingManager();
  ~OutputWritingManager();
  std::string mOutputDirectoryPath;
  std::string mLastEncounteredError;

  std::string GetLastEncounteredError()
  {
    return mLastEncounteredError;
  }
  std::string GetOutputDirectoryPath()
  {
    return mOutputDirectoryPath;
  }
  void SetOutputDirectoryPath(std::string InputPathFileName)
  {
    mOutputDirectoryPath = InputPathFileName;
  }

  template<class ImageType>
  void WriteOutput(typename ImageType::Pointer T1CEImagePointer,
    typename ImageType::Pointer T2FlairImagePointer,
    typename ImageType::Pointer T1ImagePointer,
    typename ImageType::Pointer T2ImagePointer,
    std::vector<typename ImageType::Pointer> PerfusionImagePointer,
    std::vector<typename ImageType::Pointer> DTIImagePointer,
    std::vector<std::string> filenames, int writingtype,
    bool convDataPresent, bool perfusionDataPresent, bool dtiDataPresent);

  template<class ImageType>
  void WriteRecurrenceOutput(typename ImageType::Pointer RecurrenceProbabilityMap, std::string T1CEPathName, const std::string &recurrenceMapFile);

  template<class ImageType>
  void WriteRecurrenceOutputInNifti(typename ImageType::Pointer RecProbabilityMap,
    std::string filepath);

  template<class ImageType>
  void WriteNearFarMasks(typename ImageType::Pointer image, std::string subjectname, std::vector<typename ImageType::IndexType> nearIndices, std::vector<typename ImageType::IndexType> farIndices);

  template<class ImageType>
  void WriteSeedMasks(typename ImageType::Pointer image, std::string subjectname, std::vector<typename ImageType::IndexType> seedIndices);

  template<class ImageType>
  void WriteImageWithGivenName(typename ImageType::Pointer ImagePointer, std::string FileName);

  bool SetupOutputFolders();
  void SaveModelResults(VariableSizeMatrixType scaledFeatures, VariableLengthVectorType means, VariableLengthVectorType stds, VariableLengthVectorType pmeans, VariableSizeMatrixType pcacoeff,
    bool useConvData,  bool useDTIData, bool usePerfData, bool useDistData, const int size);

  vnl_matrix<double> ReadNumberOfModalities(std::string modalitiesfile);
  void ReadModelParameters(std::string meanfile, std::string stdfile, std::string pcafile, std::string pmeanfile, VariableLengthVectorType & mean, VariableLengthVectorType & stds, VariableSizeMatrixType & pca, VariableLengthVectorType & pmean);
  void ReadModelParameters(std::string meanfile, std::string stdfile, VariableLengthVectorType & mean, VariableLengthVectorType & stds);

  template<class ImageType>
  void WriteOutputNifti(typename ImageType::Pointer image, std::string filename);
};

template<class ImageType>
void OutputWritingManager::WriteNearFarMasks(typename ImageType::Pointer image, std::string filename, std::vector<typename ImageType::IndexType> nearIndices, std::vector<typename ImageType::IndexType> farIndices)
{
  itk::ImageRegionIterator<ImageType> imageIterator(image, image->GetLargestPossibleRegion());
  while (!imageIterator.IsAtEnd())
  {
    imageIterator.Set(0);
    ++imageIterator;
  }
  for (unsigned int i = 0; i < nearIndices.size(); i++)
    image->SetPixel(nearIndices[i], /*NEAR_POINT_SAVE_LABEL*/1);
  for (unsigned int i = 0; i < farIndices.size(); i++)
    image->SetPixel(farIndices[i], /*FAR_POINT_SAVE_LABEL*/2);

  //itk::NiftiImageIO::Pointer nifti_io = itk::NiftiImageIO::New();
  typedef itk::ImageFileWriter< ImageType > WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(filename);
  //writer->SetImageIO(nifti_io);
  writer->SetInput(image);
  writer->Update();
}

template<class ImageType>
void OutputWritingManager::WriteSeedMasks(typename ImageType::Pointer image, std::string filename, std::vector<typename ImageType::IndexType> seedIndices)
{
  itk::ImageRegionIterator<ImageType> imageIterator(image, image->GetLargestPossibleRegion());
  image->FillBuffer(0);

  //filename = mOutputDirectoryPath + "C:/Projects/" + subjectname + "_InitialSeeds.nii.gz";
  for (unsigned int i = 0; i < seedIndices.size(); i++)
  {
    imageIterator.SetIndex(seedIndices[i]);
    imageIterator.Set(150/*used to be INIT_POINT_SAVE_LABEL*/);
  }

  cbica::WriteImage< ImageType >(image, filename);
}

template<class ImageType>
void OutputWritingManager::WriteOutputNifti(typename ImageType::Pointer image, std::string filename)
{
	//itk::NiftiImageIO::Pointer nifti_io = itk::NiftiImageIO::New();
	typedef itk::ImageFileWriter< ImageType > WriterType;
	typename WriterType::Pointer writer1 = WriterType::New();
	writer1->SetFileName(filename + ".nii.gz");
	//writer1->SetImageIO(nifti_io);
	writer1->SetInput(image);
	writer1->Update();
	return;

}
template<class ImageType>
void OutputWritingManager::WriteRecurrenceOutputInNifti(typename ImageType::Pointer RecProbabilityMap, std::string filename)
{
  cbica::WriteImage< ImageType >(RecProbabilityMap, filename + "_RecurrenceMap.nii.gz");

  ////itk::NiftiImageIO::Pointer nifti_io = itk::NiftiImageIO::New();
  //typedef itk::ImageFileWriter< ImageType > WriterType;
  //typename WriterType::Pointer writer1 = WriterType::New();
  //writer1->SetFileName(filename + "_RecurrenceMap.nii.gz");
  ////writer1->SetImageIO(nifti_io);
  //writer1->SetInput(RecProbabilityMap);
  //writer1->Update();
  //return;


  //typename WriterType::Pointer writer2 = WriterType::New();
  //writer2->SetFileName(filename + "_npMap.nii.gz");
  ////writer2->SetImageIO(nifti_io);
  //writer2->SetInput(NonRecProbabilityMap);
  //writer2->Update();


  //typename WriterType::Pointer writer = WriterType::New();
  //writer->SetFileName(filename + "_lMap.nii.gz");
  ////writer->SetImageIO(nifti_io);
  //writer->SetInput(LabelMap);
  //writer->Update();
}

template<class ImageType>
void OutputWritingManager::WriteRecurrenceOutput(typename ImageType::Pointer RecurrenceProbabilityMap, std::string T1CEPathName, const std::string &recurrenceMapFile)
{
  cbica::WriteImage< ImageType >(RecurrenceProbabilityMap, mOutputDirectoryPath + "/" + recurrenceMapFile);

  ////itk::NiftiImageIO::Pointer nifti_io = itk::NiftiImageIO::New();
  //typedef itk::ImageFileWriter< ImageType > WriterType;
  //typename WriterType::Pointer writer = WriterType::New();
  //std::string filename = mOutputDirectoryPath + "/" + recurrenceMapFile;
  //writer->SetFileName(filename);
  ////writer->SetImageIO(nifti_io);
  //writer->SetInput(RecurrenceProbabilityMap);
  //writer->Update();

  return;

  //std::string directoryname = mOutputDirectoryPath + "/RecurrenceOutput";
  //typedef itk::Image<unsigned short, 3> DicomImageType;
  //typedef itk::ImageSeriesReader<DicomImageType> ReaderType;
  //ReaderType::Pointer seriesreader = ReaderType::New();


  //typedef itk::GDCMImageIO ImageIOType;
  //ImageIOType::Pointer dicomIO = ImageIOType::New();
  //seriesreader->SetImageIO(dicomIO);
  //typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  //NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
  //nameGenerator->SetUseSeriesDetails(true);
  //nameGenerator->AddSeriesRestriction("0008|0021");
  //nameGenerator->SetInputDirectory(T1CEPathName);

  //try
  //{
  //  typedef std::vector< std::string > SeriesIdContainer;
  //  const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();

  //  SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
  //  SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
  //  while (seriesItr != seriesEnd)
  //  {
  //    typedef std::vector< std::string > FileNamesContainer;
  //    FileNamesContainer fileNames;
  //    fileNames = nameGenerator->GetFileNames(seriesItr->c_str());
  //    seriesreader->SetFileNames(fileNames);
  //    try
  //    {
  //      seriesreader->Update();
  //    }
  //    catch (itk::ExceptionObject & err)
  //    {
  //      std::stringstream error;
  //      error << err;
  //    }
  //    ++seriesItr;
  //  }
  //}
  //catch (itk::ExceptionObject & err)
  //{
  //  std::stringstream error;
  //  error << err;
  //}
  //itk::ImageRegionIterator<ImageType> imageIterator(RecurrenceProbabilityMap, RecurrenceProbabilityMap->GetLargestPossibleRegion());
  //while (!imageIterator.IsAtEnd())
  //{
  //  imageIterator.Set(imageIterator.Get() * 1000);
  //  ++imageIterator;
  //}

  ///*int x_size = region.GetSize()[0];
  //int y_size = region.GetSize()[1];
  //int z_size = region.GetSize()[2];

  //for (int i = 0; i < x_size; i++)
  //for (int j = 0; j < y_size; j++)
  //for (int k = 0; k < z_size; k++)
  //{
  //ImageType::IndexType index;
  //index[0] = i;
  //index[1] = j;
  //index[2] = k;
  //float val = RecurrenceProbabilityMap->GetPixel(index);
  //val = val * 1000;
  //RecurrenceProbabilityMap->SetPixel(index, val);
  //}*/
  //typedef itk::CastImageFilter< ImageType, DicomImageType > CastFilterType;
  //typename CastFilterType::Pointer castFilter = CastFilterType::New();
  //castFilter->SetInput(RecurrenceProbabilityMap);



  //typedef unsigned short OutputPixelType;
  //const unsigned int OutputDimension = 2;
  //typedef itk::Image<OutputPixelType, OutputDimension> Image2DType;
  //typedef itk::ImageSeriesWriter<DicomImageType, Image2DType> SeriesWriterType;
  //SeriesWriterType::Pointer serieswriter = SeriesWriterType::New();
  //serieswriter->SetInput(castFilter->GetOutput());
  //serieswriter->SetImageIO(dicomIO);
  //nameGenerator->SetOutputDirectory(directoryname);
  //serieswriter->SetFileNames(nameGenerator->GetOutputFileNames());
  //serieswriter->SetMetaDataDictionaryArray(seriesreader->GetMetaDataDictionaryArray());
  //try
  //{
  //  serieswriter->Update();
  //}
  //catch (itk::ExceptionObject & excp)
  //{
  //  cbica::Logging(loggerFile, "Error caught writing DICOM files: " + std::string(excp.GetDescription()));
  //  exit(EXIT_FAILURE);
  //}

}

template<class ImageType>
void OutputWritingManager::WriteOutput(typename ImageType::Pointer T1CEImagePointer,
  typename ImageType::Pointer T2FlairImagePointer,
  typename ImageType::Pointer T1ImagePointer,
  typename ImageType::Pointer T2ImagePointer,
  std::vector<typename ImageType::Pointer> PerfusionImagePointer,
  std::vector<typename ImageType::Pointer> DTIImagePointer,
  std::vector<std::string> filenames, int writingtype,
  bool convDataPresent, bool perfusionDataPresent, bool dtiDataPresent)
{
  std::string finalprefix = "";
  if (writingtype == WritingType::SUSAN)
    finalprefix = "_susan.nii.gz";
  else if (writingtype == WritingType::BIAS_CORRECT)
    finalprefix = "_n3.nii.gz";
  else if (writingtype == WritingType::REGISTRATION)
    finalprefix = "_registration.nii.gz";

  //itk::NiftiImageIO::Pointer nifti_io = itk::NiftiImageIO::New();
  //typedef itk::ImageFileWriter<ImageType > WriterType;
  //typename WriterType::Pointer writer = WriterType::New();
  std::string filename;
  if (convDataPresent)
  {
    cbica::WriteImage< ImageType >(T1CEImagePointer, mOutputDirectoryPath + "/T1CE/" + filenames[0] + finalprefix);
    cbica::WriteImage< ImageType >(T2FlairImagePointer, mOutputDirectoryPath + "/T2Flair/" + filenames[1] + finalprefix);
    cbica::WriteImage< ImageType >(T1ImagePointer, mOutputDirectoryPath + "/T1/" + filenames[2] + finalprefix);
    cbica::WriteImage< ImageType >(T2ImagePointer, mOutputDirectoryPath + "/T2/" + filenames[3] + finalprefix);

    //filename = mOutputDirectoryPath + "/T1CE/" + filenames[0] + finalprefix;
    //writer->SetFileName(filename);
    //writer->SetInput(T1CEImagePointer);
    //writer->Update();

	  //filename = mOutputDirectoryPath + "/T2Flair/" + filenames[1] + finalprefix;
    //writer->SetFileName(filename);
    //writer->SetInput(T2FlairImagePointer);
    //writer->Update();

	   //filename = mOutputDirectoryPath + "/T1/" + filenames[2] + finalprefix;
    //writer->SetFileName(filename);
    //writer->SetInput(T1ImagePointer);
    //writer->Update();

	   //filename = mOutputDirectoryPath + "/T2/" + filenames[3] + finalprefix;
    //writer->SetFileName(filename);
    //writer->SetInput(T2ImagePointer);
    //writer->Update();
  }
  if (perfusionDataPresent)
  {
    for (unsigned int i = 0; i < PerfusionImagePointer.size(); i++)
    {
      cbica::WriteImage< ImageType >(PerfusionImagePointer[i], mOutputDirectoryPath + "/Perfusion/" + filenames[4] + "_" + std::to_string(i) + finalprefix);
      //filename = mOutputDirectoryPath + "/Perfusion/" + filenames[4] + "_" + std::to_string(i) + finalprefix;
      //writer->SetFileName(filename);
      ////writer->SetImageIO(nifti_io);
      //writer->SetInput(PerfusionImagePointer[i]);
      //writer->Update();
    }
  }
  if (dtiDataPresent)
  {
    for (unsigned int i = 0; i < DTIImagePointer.size(); i++)
    {
      cbica::WriteImage< ImageType >(PerfusionImagePointer[i], mOutputDirectoryPath + "/DTI/" + filenames[5] + "_" + std::to_string(i) + finalprefix);
      //filename = mOutputDirectoryPath + "/DTI/" + filenames[5] + "_" + std::to_string(i) + finalprefix;
      //writer->SetFileName(filename);
      ////writer->SetImageIO(nifti_io);
      //writer->SetInput(DTIImagePointer[i]);
      //writer->Update();
    }
  }
}
template<class ImageType>
void OutputWritingManager::WriteImageWithGivenName(typename ImageType::Pointer ImagePointer, std::string filename)
{
  cbica::WriteImage< ImageType >(ImagePointer, mOutputDirectoryPath + filename);
  ////itk::NiftiImageIO::Pointer nifti_io = itk::NiftiImageIO::New();
  //typedef itk::ImageFileWriter< ImageType > WriterType;
  //typename WriterType::Pointer writer = WriterType::New();

  //writer->SetFileName(filename);
  ////writer->SetImageIO(nifti_io);
  //writer->SetInput(ImagePointer);
  //writer->Update();
}

