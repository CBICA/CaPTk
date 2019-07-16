#ifndef DICOMIOMANAGER_H
#define DICOMIOMANAGER_H

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include <string>

namespace cbica
{
  /**
  \brief Check if the given file is a valid DICOM image or not

  \param fileNameToCheck The input file
  */
  inline bool IsDicom(const std::string fileNameToCheck)
  {
    gdcm::Reader reader;
    reader.SetFileName(fileNameToCheck.c_str());
    return reader.CanRead();
  }
}

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

  //! check if file is dicom
  static bool IsDicom(std::string path);

  //! check if file can be read or not and return the base imageIO if readable
  static bool CanReadFile(std::string path, itk::ImageIOBase::Pointer &imageIO);

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

#include "DicomIOManager.hxx"
#endif // DICOMIOMANAGER_H