#ifndef DICOMIOMANAGER_H
#define DICOMIOMANAGER_H

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include <string>

template <class T>
class DicomIOManager
{
public:

  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  typedef itk::GDCMImageIO ImageIOType;
  typedef std::vector<std::string> FileNamesContainer;

  DicomIOManager();
  ~DicomIOManager();

  //! set the input directory containing dicom series
  void SetDirectoryPath(std::string path);

  //! get the read dicom data as 3D float ITK image
  typename T::Pointer GetITKImage();

  //! load dicom data
  bool LoadDicom();

  //! helper to write itk image
  //template <class TInputImage>
  //void WriteITKImage(typename T::Pointer image, std::string filename);

  //! convert itk image to float 3D itk image
  template <class TInputImage>
  typename T::Pointer ConvertImage3DToFloatImage3D(typename TInputImage::Pointer image);

private:
  std::string m_dir;  //! input directory path
  typename T::Pointer m_image3d;
};

#include "dicomiomanager.hxx"
#endif // DICOMIOMANAGER_H