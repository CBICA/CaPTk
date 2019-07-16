
/**
\file  SimpleImageManager.h

\brief Declaration of SimpleImageManager class

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

#ifndef _SimpleImageManager_h_
#define _SimpleImageManager_h_

//#include "CAPTk.h"
#include "cbicaDTIProcessingManager.h"
#include <QString>
#include "CaPTkDefines.h"
#include "cbicaLogging.h"
#include <vtksys/SystemTools.hxx>
#include "CaPTkEnums.h"

//#include "cbicaDTI.h"

class SimpleImageManager
{
public:
  SimpleImageManager();
  ~SimpleImageManager();
  std::string mFileName;
  std::string mPathFileName;
  std::string mBaseFileName;
  int mImageType;
  int mImageSubType;

  std::string mLastEncounteredError;
  cbica::DTIProcessingManager mDTIObject;
  //cbica::cbicaDTI *mDTIProcessingObject;

  std::string GetLastEncounteredError()
  {
    return mLastEncounteredError;
  }
  std::string GetPathFileName()
  {
    return mPathFileName;
  }
  std::string GetFileName()
  {
    return mFileName;
  }
  std::string GetBaseFileName()
  {
    return mBaseFileName;
  }
  void SetPathFileName(const std::string &InputPathFileName)
  {
    mPathFileName = InputPathFileName;
  }
  void SetFileName(const std::string &InputFileName)
  {
    mPathFileName = InputFileName;
    mFileName = vtksys::SystemTools::GetFilenameName(mPathFileName);
    mBaseFileName = vtksys::SystemTools::GetFilenameName(vtksys::SystemTools::GetFilenameWithoutLastExtension(mFileName));
  }
  void SetBaseFileName(const std::string &InputBaseFileName)
  {
    mBaseFileName = InputBaseFileName;
  }
  bool ReadGivenNonViewingDicomImage(std::string directoryname, int imagetype);
  bool ReadGivenNonViewingNiftiImage(std::string directoryname, int subimagetype);
  std::string GetFileNameInADicomDirectory(QString directoryname);

  template<unsigned int VImageDimension>
  void ReadImageWithGivenDimensionValue(std::string filename, std::string InputPixelType, int imagetype);

  template<class InputPixelType, unsigned int VImageDimension>
  void ReadImageWithDimAndInputPixelType(std::string directoryName, int imagetype);

  std::vector<ImageTypeFloat3D::Pointer> mNVImageSeriesReaderPointer;

  ImageTypeFloat4D::Pointer mPerfusionImagePointer;



};

template<unsigned int VImageDimension>
void SimpleImageManager::ReadImageWithGivenDimensionValue(std::string filename, std::string InputPixelType, int imagetype)
{
  if (InputPixelType == "short")
  {
    ReadImageWithDimAndInputPixelType<short, VImageDimension>(filename, imagetype);
  }
  else if (InputPixelType == "signed short")
  {
    ReadImageWithDimAndInputPixelType<short, VImageDimension>(filename, imagetype);
  }
  else if (InputPixelType == "unsigned_short" || InputPixelType == "unsigned short")
  {
    ReadImageWithDimAndInputPixelType<unsigned short, VImageDimension>(filename, imagetype);
  }
  else if (InputPixelType == "char")
  {
    ReadImageWithDimAndInputPixelType<char, VImageDimension>(filename, imagetype);
  }
  else if (InputPixelType == "unsigned_char")
  {
    ReadImageWithDimAndInputPixelType<unsigned char, VImageDimension>(filename, imagetype);
  }
  else if (InputPixelType == "int")
  {
    ReadImageWithDimAndInputPixelType<int, VImageDimension>(filename, imagetype);
  }
  else if (InputPixelType == "unsigned_int")
  {
    ReadImageWithDimAndInputPixelType<unsigned int, VImageDimension>(filename, imagetype);
  }
  else if (InputPixelType == "double")
  {
    ReadImageWithDimAndInputPixelType<double, VImageDimension>(filename, imagetype);
  }
  else if (InputPixelType == "float")
  {
    ReadImageWithDimAndInputPixelType<float, VImageDimension>(filename, imagetype);
  }
  else
  {
    cbica::Logging(loggerFile, "Error, input pixel type, '" + InputPixelType + "' is unknown!");
  }
}


template<class InputPixelType, unsigned int VImageDimension>
void SimpleImageManager::ReadImageWithDimAndInputPixelType(std::string directoryName, int imagetype)
{
  if (imagetype == CAPTK::ImageModalityType::IMAGE_TYPE_DTI)
  {
    //	int i;
    //	if (!system(NULL))
    //	{
    //		mLastEncounteredError = "System call is not possible for DTI data loading.";
    //			return;
    //	}
    //	std::string cwd = cbica::getCWD();
    //	cwd = cbica::replaceString(cwd, "src/", "");

    //	std::string tempDirectoryName = "";
    //	cbica::createTmpDir(tempDirectoryName);
    //	std::string systemCall = cwd +"data/dtistuff/dtistuff/dcm2nii -o " + tempDirectoryName + " " + directoryName;
    //	i = system(systemCall.c_str());
    //	if (i != 0)
    //	{
    //		mLastEncounteredError = "Error reading dicom data.";
    //		return;
    //	}
    //	std::vector<std::string> fileNames;
    //	fileNames = cbica::filesInDirectory(tempDirectoryName);
    //	std::string bvalFile	= "";
    //	std::string bvecFile	= "";
    //	std::string dwiFile		= "";
    //for (unsigned int i = 0; i < fileNames.size(); i++)
    //{
    //	if (cbica::getFilenameExtension(fileNames[i]) == "bval")
    //		bvalFile = fileNames[i];
    //	else if (cbica::getFilenameExtension(fileNames[i]) == "bvec")
    //		bvecFile = fileNames[i];
    //	else if (cbica::getFilenameExtension(fileNames[i]) == "gz")
    //		dwiFile = fileNames[i];
    //}
    //mDTIObject->FlipYOrientationInBVecFile(tempDirectoryName + bvecFile);
    //mDTIObject->FlipAndShiftNiftiFile(tempDirectoryName + dwiFile, tempDirectoryName + "lps_" + dwiFile);
    //std::string image = tempDirectoryName + "lps_" + dwiFile;
    //std::string bval = tempDirectoryName + bvalFile;
    //std::string bvec = tempDirectoryName + bvecFile;
    //std::string mask = cwd + "data/dtistuff/dtistuff/mask.nii.gz";
    //typedef itk::Image<float, 3> ScalarImageType;
    //std::vector<ScalarImageType::Pointer> vectorOfDTIScalars;

    //vectorOfDTIScalars = mDTIObject->dtiMainFunction(image, mask, bval, bvec, tempDirectoryName);
    //for (unsigned int i = 0; i < vectorOfDTIScalars.size();i++)
    //	mNVImageSeriesReaderPointer.push_back(vectorOfDTIScalars[i]);
    //mImageType = imagetype;
    //cbica::removeDirectoryRecursively(tempDirectoryName,true);

    std::string maskFileName = "";
    typedef ImageTypeFloat3D ScalarImageType;
    std::vector<ScalarImageType::Pointer> vectorOfDTIScalars = mDTIObject.ConvertDWIToScalars(directoryName, maskFileName);
    for (unsigned int i = 0; i < vectorOfDTIScalars.size(); i++)
      mNVImageSeriesReaderPointer.push_back(vectorOfDTIScalars[i]);
  }
  else
  {
    typedef itk::Image<InputPixelType, VImageDimension> InputImageType;
    typedef itk::ImageSeriesReader< InputImageType > ReaderType;
    typename ReaderType::Pointer seriesreader = ReaderType::New();


    typedef itk::GDCMImageIO ImageIOType;
    ImageIOType::Pointer dicomIO = ImageIOType::New();
    seriesreader->SetImageIO(dicomIO);
    typedef itk::GDCMSeriesFileNames NamesGeneratorType;
    NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
    nameGenerator->SetUseSeriesDetails(true);
    nameGenerator->AddSeriesRestriction("0008|0021");
    nameGenerator->AddSeriesRestriction("0020|0012");   //for perfusion data
    nameGenerator->SetInputDirectory(directoryName);
    try
    {
      typedef std::vector< std::string > SeriesIdContainer;
      const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();

      SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
      SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
      while (seriesItr != seriesEnd)
      {
        mImageType = imagetype;
        typedef std::vector< std::string > FileNamesContainer;
        FileNamesContainer fileNames;
        fileNames = nameGenerator->GetFileNames(seriesItr->c_str());
        seriesreader->SetFileNames(fileNames);
        try
        {
          seriesreader->Update();
          typedef ImageTypeFloat3D OutputImageType;
          typedef itk::CastImageFilter< InputImageType, OutputImageType > CastFilterType;
          typename CastFilterType::Pointer castFilter = CastFilterType::New();
          castFilter->SetInput(seriesreader->GetOutput());
          typename OutputImageType::Pointer outputImage = castFilter->GetOutput();
          outputImage->Update();
          mNVImageSeriesReaderPointer.push_back(outputImage);
        }
        catch (itk::ExceptionObject & err)
        {
          std::stringstream error;
          error << err;
          mLastEncounteredError = error.str();
        }
        ++seriesItr;
      }
    }
    catch (itk::ExceptionObject & err)
    {
      std::stringstream error;
      error << err;
      mLastEncounteredError = error.str();

      return;
    }
  }
}
#endif