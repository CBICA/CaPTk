///////////////////////////////////////////////////////////////////////////////////////
// DicomMetadataReader.h
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

#ifndef DICOMMMETADATAREADER_H
#define DICOMMMETADATAREADER_H

#include "itkImageSeriesReader.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include <itkMapContainer.h>

class DicomMetadataReader
{
public:

  //! we may need to handle pixeltype in a template if we run into issues
  typedef signed short       PixelType;
  const unsigned int         Dimension = 2;
  typedef itk::Image< PixelType, 2 >      ImageType; 

  typedef itk::ImageFileReader< ImageType >     ReaderType;
  typedef itk::GDCMSeriesFileNames     NamesGeneratorType;
  typedef itk::GDCMImageIO       ImageIOType;
  typedef std::vector<std::string>    FileNamesContainer;
  typedef itk::MetaDataDictionary   DictionaryType;
  typedef itk::MetaDataObject< std::string > MetaDataStringType;
  typedef itk::MapContainer<std::string, std::pair<std::string, std::string>> MapContainerType;     //map description: <tag, <description, value> >

  //! Constructor/Destructor
  DicomMetadataReader();
  ~DicomMetadataReader();

  //! Set the input file path to dicom image
  void SetFilePath(std::string path);

  //! Read meta data
  bool ReadMetaData();

  //! Get all metadata as a map <tag, <description, value> >
  std::map<std::string, std::pair<std::string, std::string>> GetMetaDataMap();

  //! Get individual tag description and value for a given tag
  bool GetTagValue(std::string tag, std::string &label, std::string &value);

private:
  void PrintMetaData(); // for testing purpose
   
  ReaderType::Pointer m_reader;
  NamesGeneratorType::Pointer m_nameGenerator;
  FileNamesContainer m_fileNames;
  ImageIOType::Pointer m_dicomIO;
  MapContainerType::Pointer m_tagvalueMap;
  DictionaryType m_dictionary;
  std::string m_FilePath;
};

#endif // DICOMMMETADATAREADER_H

