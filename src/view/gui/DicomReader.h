///////////////////////////////////////////////////////////////////////////////////////
// DicomReader.h
//
// Copyright (c) 2018. All rights reserved.
// Section of Biomedical Image Analysis
// Center for Biomedical Image Computing and Analytics
// Department of Radiology
// Perelman School of Medicine
// University of Pennsylvania
//
// Contact details: software@cbica.upenn.edu
//
// License Agreement: https://www.med.upenn.edu/sbia/software-agreement.html
///////////////////////////////////////////////////////////////////////////////////////

#ifndef DCMREADER_H
#define DCMREADER_H

#include "itkImageSeriesReader.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include <itkMapContainer.h>

class DicomReader
{
public:
  DicomReader();
  ~DicomReader();
  void SetDirectoryPath(std::string path);
  void ReadMetaData();

  std::map<std::string, std::pair<std::string, std::string>> GetMetaDataMap();
  bool GetTagValue(std::string tag, std::string &label, std::string &value);

private:
  void PrintMetaData(); // for testing purpose

  typedef signed short       PixelType;
  const unsigned int         Dimension = 3;

  typedef itk::Image< PixelType, 3 >      ImageType;

  typedef itk::ImageSeriesReader< ImageType >     ReaderType;
  typedef itk::GDCMSeriesFileNames     NamesGeneratorType;
  typedef itk::GDCMImageIO       ImageIOType;
  typedef std::vector<std::string>    FileNamesContainer;
  typedef itk::MetaDataDictionary   DictionaryType;
  typedef itk::MetaDataObject< std::string > MetaDataStringType;

  ReaderType::Pointer reader;
  NamesGeneratorType::Pointer nameGenerator;
  FileNamesContainer fileNames;
  ImageIOType::Pointer dicomIO;
  typedef itk::MapContainer<std::string, std::pair<std::string, std::string>> MapContainerType;     //map description: <tag, <description, value> >
  MapContainerType::Pointer tagvalueMap;
  DictionaryType dictionary;
};

#endif // DCMREADER_H

