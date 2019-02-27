/**
\file  cbicaITKImageInfo.cpp

\brief Implementation of the ImageInfo class

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2016 University of Pennsylvania. All rights reserved. <br>
See COPYING file or http://www.med.upenn.edu/sbia/software/license.html

*/
#include "cbicaITKImageInfo.h"

#include "cbicaUtilities.h"

namespace cbica
{
  // ============================================================================ //
  
  ImageInfo::ImageInfo( const std::string &fName )
  {
    m_fileName = fName;
    m_itkImageIOBase = itk::ImageIOFactory::CreateImageIO( 
                       fName.c_str(), itk::ImageIOFactory::ReadMode );
    
    // exception handling in case of NULL pointer initialization
    if ( m_itkImageIOBase )
    {
      m_itkImageIOBase->SetFileName(fName);
      m_itkImageIOBase->ReadImageInformation();
    }
    else
    {
      itkGenericExceptionMacro("Could not read the input image information from '" << fName << "'\n");
    }

    m_itkImageIOBase->SetFileName(fName);
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