/**
\file  cbicaITKImageInfo.cpp

\brief Implementation of the ImageInfo class

https://www.med.upenn.edu/sbia/software/ <br>
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
    if (cbica::isFile(fName) /*&& (fName_ext != ".dcm")*/)
    {
      auto fName_ext = cbica::getFilenameExtension(fName);
      if (!fName_ext.empty())
      {
        std::transform(fName_ext.begin(), fName_ext.end(), fName_ext.begin(), ::tolower);
      }

      m_fileName = fName_norm;
    }
    else if (cbica::isDir(fName_norm))
    {
      auto filesInDir = cbica::filesInDirectory(fName_norm);
      if (!filesInDir.empty())
      {
        if (cbica::IsDicom(filesInDir[0]))
        {
          m_fileName = filesInDir[0];
        }
        else
        {
          std::cerr << "Please pass a single file.\n";
          return;
        }
      }
      else
      {
        std::cerr << "Empty Directory passed.\n";
        return;
      }
    }

    m_itkImageIOBase = itk::ImageIOFactory::CreateImageIO(m_fileName.c_str(), itk::ImageIOFactory::ReadMode);

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
    m_dimensions = m_itkImageIOBase->GetNumberOfDimensions();

    for (size_t i = 0; i < m_dimensions; i++)
    {
      m_spacings.push_back(m_itkImageIOBase->GetSpacing(i));
      m_origins.push_back(m_itkImageIOBase->GetOrigin(i));
      m_size.push_back(m_itkImageIOBase->GetDimensions(i));
      m_directions.push_back(m_itkImageIOBase->GetDirection(i));
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

  std::vector< std::vector< double > > ImageInfo::GetImageDirections()
  {
    return m_directions;
  }

}