/**
\file  SimpleImageManager.cpp

\brief Implementation of SimpleImageManager class

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

#include "SimpleImageManager.h"
#include "NiftiDataManager.h"
#include "CaPTkEnums.h"
#include <QDirIterator>
#include <QStringList>
#include <QDir>

SimpleImageManager::SimpleImageManager()
{
  mPathFileName = "";
  mFileName = "";
  mBaseFileName = "";
  mLastEncounteredError = "";
  //mDTIObject = new cbica::DTIProcessingManager();
  //mDTIProcessingObject  = new cbica::cbicaDTI();
  mImageType = CAPTK::ImageModalityType::IMAGE_TYPE_UNDEFINED;
  mImageSubType = CAPTK::ImageModalityType::IMAGE_TYPE_UNDEFINED;
}

SimpleImageManager::~SimpleImageManager()
{
}

bool SimpleImageManager::ReadGivenNonViewingNiftiImage(std::string directoryname, int subimagetype)
{
  itk::Image<float, 3>::Pointer image;
  itk::ImageIOBase::Pointer reader = itk::ImageIOFactory::CreateImageIO(directoryname.c_str(), itk::ImageIOFactory::ReadMode);
  if (!reader)
    return false;

  reader->SetFileName(directoryname);
  reader->ReadImageInformation();
  std::string InputPixelType = reader->GetComponentTypeAsString(reader->GetComponentType());
  // unsigned int VImageDimension = reader->GetNumberOfDimensions();

  NiftiDataManager obj;
  if (subimagetype == CAPTK::ImageModalityType::IMAGE_TYPE_DTI)
    this->mNVImageSeriesReaderPointer.push_back(obj.ReadNiftiImage(directoryname));
  if (subimagetype == CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION)
    this->mPerfusionImagePointer = obj.Read4DNiftiImage(directoryname);
  return true;
}


bool SimpleImageManager::ReadGivenNonViewingDicomImage(std::string directoryname, int imagetype)
{
  std::string FileName = GetFileNameInADicomDirectory(QString::fromStdString(directoryname));
  unsigned int VImageDimension = 0;
  std::string InputPixelType = "";
  std::string OutputPixelType = "";

  itk::ImageIOBase::Pointer reader = itk::ImageIOFactory::CreateImageIO(FileName.c_str(), itk::ImageIOFactory::ReadMode);
  if (!reader)
  {
    mLastEncounteredError = "Unable to read given file.";
    return false;
  }
  reader->SetFileName(FileName);
  reader->ReadImageInformation();
  InputPixelType = reader->GetComponentTypeAsString(reader->GetComponentType());
  OutputPixelType = reader->GetComponentTypeAsString(reader->GetComponentType());
  VImageDimension = reader->GetNumberOfDimensions();

  if (VImageDimension == 2)
    ReadImageWithGivenDimensionValue<2>(directoryname, InputPixelType, imagetype);
  else
    ReadImageWithGivenDimensionValue<3>(directoryname, InputPixelType, imagetype);

  return true;
}
std::string SimpleImageManager::GetFileNameInADicomDirectory(QString directoryname)
{
  std::string FirstFileName = "";
  QDirIterator directoryIterator(directoryname, QStringList() << "*.dcm", QDir::Files, QDirIterator::Subdirectories);
  while (directoryIterator.hasNext())
  {
    FirstFileName = directoryIterator.next().toStdString();
    break;
  }
  return FirstFileName;
}


