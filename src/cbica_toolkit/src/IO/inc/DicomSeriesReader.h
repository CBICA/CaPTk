///////////////////////////////////////////////////////////////////////////////////////
// DicomSeriesReader.h
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

#ifndef DICOMSERIESREADER_H
#define DICOMSERIESREADER_H

#include "itkImageSeriesReader.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include <itkMapContainer.h>

class DicomSeriesReader
{
public:

  typedef itk::Image< float, 3 > ImageType3DFloat;
  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  typedef itk::GDCMImageIO ImageIOType;
  typedef std::vector<std::string> FileNamesContainer;

  DicomSeriesReader();
  ~DicomSeriesReader();

  //! set the input directory containing dicom series
  void SetDirectoryPath(std::string path);
 
  //! get the read dicom data as 3D float ITK image
  DicomSeriesReader::ImageType3DFloat::Pointer GetITKImage();

  //! load dicom data
  bool LoadDicom();

  //! Read dicom series
  template <class TInputImage>
  typename TInputImage::Pointer ReadDicomSeries(bool &readStatus);

  ////! helper to write itk image
  //template <class TInputImage>
  //void WriteITKImage(typename TInputImage::Pointer image, std::string filename);

  ////! convert itk image to float 3D itk image
  //template <class TInputImage>
  //DicomSeriesReader::ImageType3DFloat::Pointer ConvertImage3DToFloatImage3D(typename TInputImage::Pointer image);

private:

  ImageType3DFloat::Pointer m_image3dfloat; //! image 3D as float
  std::string m_dir;  //! input directory path
};

template<class TInputImage>
inline typename TInputImage::Pointer DicomSeriesReader::ReadDicomSeries(bool &readStatus)
{
  readStatus = false;
  typedef itk::ImageSeriesReader< TInputImage>     DicomReaderType;
  auto reader = DicomReaderType::New();
  auto dicomIO = ImageIOType::New();
  auto nameGenerator = NamesGeneratorType::New();

  nameGenerator->SetInputDirectory(this->m_dir);
  reader->SetImageIO(dicomIO);
  std::vector<std::string> fileNames = nameGenerator->GetInputFileNames();
  reader->SetFileNames(fileNames);

  try
  {
    reader->Update();
  }
  catch (itk::ExceptionObject &ex)
  {
    readStatus = false;
    std::cout << "unsupported dicom" << std::endl;
    std::cout << ex << std::endl;
    return nullptr;
  }
  
  readStatus = true;
  return reader->GetOutput();

}
//
//template<class TInputImage>
//inline void DicomSeriesReader::WriteITKImage(typename TInputImage::Pointer image, std::string filename)
//{
//  typedef  itk::ImageFileWriter< TInputImage  > WriterType;
//  auto writer = WriterType::New();
//  writer->SetFileName(filename);
//  writer->SetInput(image);
//  writer->Update();
//}
//
//template<class TInputImage>
//inline DicomSeriesReader::ImageType3DFloat::Pointer DicomSeriesReader::ConvertImage3DToFloatImage3D(typename TInputImage::Pointer image)
//{
//  typedef itk::CastImageFilter<TInputImage, ImageType3DFloat> CastFilterType;
//  auto castFilter = CastFilterType::New();
//  castFilter->SetInput(image);
//  castFilter->Update();
//  return castFilter->GetOutput();
//}

#endif // DICOMSERIESREADER_H