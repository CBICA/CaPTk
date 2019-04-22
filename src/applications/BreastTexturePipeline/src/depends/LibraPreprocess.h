/**
\file LibraPreprocess.h

\brief This holds the Preprocessing class for LIBRA
*/

#pragma once

#include <map>

//using TImageType = itk::Image< float, 2 >;

template< class TImageType = itk::Image< float, 2 > >
class LibraPreprocess
{
public:
  using TIteratorType = itk::ImageRegionIteratorWithIndex< TImageType >;
  LibraPreprocess();
  ~LibraPreprocess();

  //! Set the File name to read
  void SetInputFileName(const std::string &input);

  //! Set the resizing factor; defaults to none
  void SetResizingFactor(size_t resize = 100);

  //! Set the resizing factor in float; defaults to none
  void SetResizingFactor(float resize = 1.0);
  
  //! Gives more output on console
  void EnableDebugMode();

  //! Actually do the preprocessing
  void Update();
  
  //! Get the preprocessed image
  typename TImageType::Pointer GetOutputImage();

  //! Get the relevant DICOM tags
  std::map< std::string, std::string > GetRelevantDicomTags();

  //! Get the original input image
  typename TImageType::Pointer GetInputImage();

  //! this is to pass onto segmentation class
  std::string GetIntermediateSaveLocation() { return m_intermediateOutputDir; };

private:
  std::string m_fileName; //! original dicom file name
  typename TImageType::Pointer m_inputImage, //! raw input image
    m_output; //! processed output
  std::map< std::string, std::string > m_dicomTags; //! the relevant dicom tags that are to be returned
  bool m_algorithmDone = false;
  bool m_debugMode = false;
  size_t m_resizingFactor;
  std::string m_intermediateOutputDir;
};

#include "LibraPreprocess.hxx"