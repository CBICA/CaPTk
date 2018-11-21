/**
\file  cbicaITKImageInfo.h

\brief Declaration of the ImageInfo class

Extracts basic information from any image (like dimensions, spacing, pixel type, etc.).

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/
#pragma once

#include <algorithm>
#include <string>
#include <vector>
//#include <tuple>

#include "itkImage.h"
#include "itkImageIOBase.h"
#include "itkImageFileReader.h"

namespace cbica
{  
  /**
  \class ImageInfo

  \brief Reads any image from file name and generates relevant data.

  It can give information like dimensions, pixel type, component type, spacing, etc. about 
  the image which mitigates the requirement of explicit declration in the program.
  */
  class ImageInfo
  {
  public:
    
    /**
    \brief The Constructor
    
    Note: There the class "itk::MetaDataDictionary" does a far better job of obtaining the metadata of
    an image, if that is what the user requires. This function just obtains some basic information, 
    namely, the spacing and dimensions of the image.
    
    \param fName The image file name for which information is required
    */
    explicit ImageInfo(const std::string &fName);
    
    /**
    \brief The Destructor
    */
    virtual ~ImageInfo();
    
    /**
    \brief Get the imageIOBase of the specified image
    
    \return An itk::ImageIOBase which is overwritten with information
    */
    itk::SmartPointer<itk::ImageIOBase> GetImageIOBase();
    
    /**
    \brief Get the Size of the specified image
    
    \return A vector of itk::SizeValueType which gets overwritten with information
    */
    std::vector<itk::SizeValueType> GetImageSize();
    
    /**
    \brief Get the Spacing of the specified image

    \return A vector of itk::SizeValueType which gets overwritten with information
    */
    std::vector<itk::SizeValueType> GetImageSpacing();
    /**

    \brief Get the dimensions of the specified image
    
    \return A const unsigned int with the number of dimensions of the image
    */
    const unsigned int GetImageDimensions()
    { 
      return m_itkImageIOBase->GetNumberOfDimensions();

    };
    
    /**
    \brief Get the Spacings of the specified image
    
    \return A vector of double which gets overwritten with information
    */
    std::vector<double> GetImageSpacings();
    
    /**
    \brief Get the Origins of the specified image

    \return A vector of double which gets overwritten with information
    */
    std::vector<double> GetImageOrigins();

    /**
    \brief Get the type of pixel in the image as a string

    \return Pixel type as a std::string
    */
    std::string GetComponentTypeAsString();
    
    /**
    \brief Get the type of pixel in the image as an itk IOComponentType

    \return Pixel type as an itk IOComponentType
    */
    itk::ImageIOBase::IOComponentType GetComponentType();
    
    /**
    \brief Get the type of pixel in the image as an itk IOComponentType

    \return Pixel type as an itk IOComponentType
    */
    std::string GetPixelTypeAsString();

    /**
    \brief Get the type of pixel in the image as an itk IOComponentType

    \return Pixel type as an itk IOComponentType
    */
    itk::ImageIOBase::IOPixelType GetPixelType();

    /**
    \brief Is the supplied image defined as a DICOM or not
    */
    bool IsDicom()
    {
      return m_dicomDetected;
    }
        
  protected:
    std::string m_fileName;
    itk::SmartPointer<itk::ImageIOBase> m_itkImageIOBase;
    std::vector<double> m_spacings;
    std::vector<double> m_origins;
    std::vector<itk::SizeValueType> m_size;
    unsigned int m_dimensions;
    std::string m_pixelType_asString, m_IOComponentType_asString;
    itk::ImageIOBase::IOComponentType m_IOComponentType;
    itk::ImageIOBase::IOPixelType m_pixelType;
    bool m_dicomDetected = false;
    std::vector< std::string > m_uids; // only for DICOM
  };
}