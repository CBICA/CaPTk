/**
\file  cbicaITKImageInfo.cpp

\brief Implementation of the ImageInfo class

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/
#include "cbicaITKImageInfo.h"

#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"

#include "gdcmReader.h"

namespace cbica
{
  // ============================================================================ //

  ImageInfo::ImageInfo(const std::string &fName)
  {
    auto fName_norm = cbica::normPath(fName);
    auto fName_ext = cbica::getFilenameExtension(fName);
    if (!fName_ext.empty())
    {
      std::transform(fName_ext.begin(), fName_ext.end(), fName_ext.begin(), ::tolower);
    }

    if (cbica::isFile(fName) /*&& (fName_ext != ".dcm")*/)
    {
      m_fileName = fName_norm;

      m_itkImageIOBase = itk::ImageIOFactory::CreateImageIO(m_fileName.c_str(), itk::ImageIOFactory::ReadMode);
      
      if (cbica::IsDicom(fName_norm))
      {
        gdcm::Reader reader;
        reader.SetFileName(m_fileName.c_str());

        if (!reader.CanRead())
        {
          itkGenericExceptionMacro("Cannot read '" << m_fileName << "'\n");
        }
      }
      
      if (!m_itkImageIOBase->CanReadFile(m_fileName.c_str()))
      {
        itkGenericExceptionMacro("Cannot read '" << m_fileName << "'\n");
      }

      // exception handling in case of NULL pointer initialization
      if (m_itkImageIOBase)
      {
        m_itkImageIOBase->SetFileName(m_fileName);
        m_itkImageIOBase->ReadImageInformation();
      }
      else
      {
        itkGenericExceptionMacro("Could not read the input image information from '" << m_fileName << "'\n");
      }

      m_itkImageIOBase->SetFileName(m_fileName);
      m_itkImageIOBase->ReadImageInformation();

      m_IOComponentType = m_itkImageIOBase->GetComponentType();
      m_pixelType = m_itkImageIOBase->GetPixelType();
      m_IOComponentType_asString = m_itkImageIOBase->GetComponentTypeAsString(m_IOComponentType);
      m_pixelType_asString = m_itkImageIOBase->GetPixelTypeAsString(m_pixelType);

      for (size_t i = 0; i<m_itkImageIOBase->GetNumberOfDimensions(); i++)
      {
        m_spacings.push_back(m_itkImageIOBase->GetSpacing(i));
        m_origins.push_back(m_itkImageIOBase->GetOrigin(i));
        m_size.push_back(m_itkImageIOBase->GetDimensions(i));
      }
    }
    //else if (fName_ext == ".dcm")
    //{
    //  m_dicomDetected = true;
    //  // nothing to do here
    //  return;
    //}
    else
    {
      std::cerr << "Please pass a supported image.\n";
      return;
    }
  }

  ImageInfo::~ImageInfo()
  {

  }

  itk::SmartPointer<itk::ImageIOBase> ImageInfo::GetImageIOBase()
  {
    return m_itkImageIOBase;
  }

  std::vector<itk::SizeValueType> ImageInfo::GetImageSize()
  {
    return m_size;
  }

  std::vector<double> ImageInfo::GetImageSpacings()
  {
    return m_spacings;
  }

  std::vector<double> ImageInfo::GetImageOrigins()
  {
    return m_origins;
  }

  std::string ImageInfo::GetComponentTypeAsString()
  {
    return m_IOComponentType_asString;
  }

  itk::ImageIOBase::IOComponentType ImageInfo::GetComponentType()
  {
    return m_IOComponentType;
  }

  std::string ImageInfo::GetPixelTypeAsString()
  {
    return m_pixelType_asString;
  }

  itk::ImageIOBase::IOPixelType ImageInfo::GetPixelType()
  {
    return m_pixelType;
  }

}