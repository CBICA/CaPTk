/**
\file  cbicaITKSafeImageIO.h

\brief Defines safe input and output of itk::Images

Read and Write itk::Image data in a safe manner. Header-only

https://www.cbica.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.cbica.upenn.edu/sbia/software/license.html

*/
#pragma once

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "itkCastImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkImageIOBase.h"
#include "itkImageIOFactory.h"
#include "itkNiftiImageIO.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
//#include "itkDCMTKImageIO.h"
//#include "itkDCMTKSeriesFileNames.h"
#include "itkNumericSeriesFileNames.h"
#include "itkOrientImageFilter.h"
#include "itkChangeInformationImageFilter.h"

#if ITK_VERSION_MAJOR >= 4
#include "gdcmUIDGenerator.h"
#else
#include "gdcm/src/gdcmFile.h"
#include "gdcm/src/gdcmUtil.h"
#endif

#include "cbicaUtilities.h"
#include "cbicaITKImageInfo.h"
#include "cbicaITKUtilities.h"
#include "DicomIOManager.h"

using ImageTypeFloat3D = itk::Image< float, 3 >;
using TImageType = ImageTypeFloat3D;
using MaskType = itk::Image<unsigned int, 3>;

namespace cbica
{ 
  /**
  \brief Get the itk::ImageFileReader from input file name. This is useful for scenarios where reader meta information is needed for later writing step(s).

  Usage:
  \verbatim
  typedef itk::Image< float, 3 > ExpectedImageType;
  std::string inputFileName = parser.getParameterValue("inputImage");
  ExpectedImageType::Pointer inputImage_1 = GetImageReader< ExpectedImageType >(inputFileName)->GetOutput();
  ExpectedImageType::Pointer inputImage_2 = GetImageReader< ExpectedImageType >(inputFileName, ".nii.gz,.img")->GetOutput();
  DoAwesomeStuffWithImage( inputImage );
  \endverbatim

  \param fName name of the image
  \param supportedExtensions Supported extensions, defaults to ".nii.gz,.nii"
  \return itk::ImageFileReader::Pointer templated over the same as requested by user
  */
  template <class TImageType = ImageTypeFloat3D >
  typename itk::ImageFileReader< TImageType >::Pointer GetImageReader(const std::string &fName, const std::string &supportedExtensions = ".nii.gz,.nii,.dcm", const std::string &delimitor = ",")
  {
    //// check read access
    //if (((_access(fName.c_str(), 4)) == -1) || ((_access(fName.c_str(), 6)) == -1))
    //{
    //  ShowErrorMessage("You don't have read access in selected location. Please check.");
    //  exit(EXIT_FAILURE);
    //}

    std::string fName_wrap = cbica::normPath(fName);

    std::string fileExtension = cbica::getFilenameExtension(fName_wrap);
    std::transform(fileExtension.begin(), fileExtension.end(), fileExtension.begin(), ::tolower);

    auto reader = /*typename*/ itk::ImageFileReader< TImageType >::New();

    if (fileExtension == ".dcm")
    {
      auto filesInDir = cbica::filesInDirectory(cbica::getFilenamePath(fName_wrap));
      if (filesInDir.size() > 1)
      {
        std::cerr << "Trying to read DICOM file. Please use DICOMImageReader.\n";
        return reader;
      }
    }

    if (supportedExtensions != "")
    {
      std::vector< std::string > extensions = cbica::stringSplit(supportedExtensions, delimitor);

      bool supportedExtensionFound = false;
      for (size_t i = 0; i < extensions.size(); i++)
      {
        if (extensions[i] == fileExtension)
        {
          supportedExtensionFound = true;
        }
      }

      if (!supportedExtensionFound)
      {
        std::cerr << "Supplied file name '" << fName_wrap << "' doesn't have a supported extension. \nSupported Extensions: " << supportedExtensions << "\n";
        return reader;
      }
    }

    // ensure that the requested image dimensions and read image dimensions match up
    auto imageInfo = cbica::ImageInfo(fName_wrap);

    // perform basic sanity check
    if ((imageInfo.GetImageDimensions() != TImageType::ImageDimension) &&
      !((TImageType::ImageDimension == 2) && (imageInfo.GetImageSize()[2] == 1))) // this to check for a 2D DICOM
    {
      std::cerr << "Image Dimension mismatch. Return image is expected to be '" << TImageType::ImageDimension <<
        "'D and doesn't match the image dimension read from the input file, which is '" << imageInfo.GetImageDimensions() << "'.\n";
      return reader;
    }

    reader->SetFileName(fName_wrap);

    auto supportedExtsVector = cbica::stringSplit(supportedExtensions, ",");

    if (std::find(supportedExtsVector.begin(), supportedExtsVector.end(), fileExtension) == supportedExtsVector.end())
    {
      std::cerr << "Extension of file doesn't match the supported extensions; can't read.\n";
      return reader;
    }

    // set image IO type
    if ((fileExtension == ".dcm") || (fileExtension == ".dicom"))
    {
      //auto ioType = itk::GDCMKImageIO::New();
      //ioType->SetFileName(fName_wrap);
      //reader->SetImageIO(ioType);
    }
    else if ((fileExtension == ".nii") || (fileExtension == ".nii.gz"))
    {
      auto ioType = itk::NiftiImageIO::New();
      ioType->SetFileName(fName_wrap);
      reader->SetImageIO(ioType);
    }

    try
    {
      reader->Update();
    }
    catch (itk::ExceptionObject& e)
    {
      std::cerr << "Exception caught while reading the image '" << fName_wrap << "': " << e.what() << "\n";
      return reader;
    }

    return reader;
  }

  /**
  \brief Returns the unique series IDs in the specified directory

  The check is only done on the DICOM tag provided, so if there are series with the same UID information (but are indeed different images),
  this function will not able to handle it.

  \param dirName The directory in question
  \param tagToCheck The tag on the basis of which the test is done; defaults to "0x0020|0x00E"
  \return Vector of Series UIDs and fileName collection pairs, with each fileName collection corresponding to a UID
  */
  //inline std::vector< std::pair< std::string , std::vector< std::string > > > GetDICOMSeriesAndFilesInDir(const std::string &dirName,
  //  const std::string tagToCheck = "0x0020|0x00E")
  //{
  //  std::vector< 
  //    std::pair< 
  //    std::string, // this is the series UID information
  //    std::vector< std::string > > // these are the fileNames corresponding to each UID
  //  > returnVector;

  //  auto dirName_wrap = cbica::normPath(dirName);
  //  auto allFilesInDir = cbica::filesInDirectory(dirName_wrap);
  //  
  //  // initialize the returnVector with the first series UID and fileName
  //  returnVector.push_back(
  //    std::make_pair(cbica::GetDICOMTagValue(allFilesInDir[0], tagToCheck), // get the first series UID 
  //    std::vector< std::string >({ allFilesInDir[0] }) // construct a initial vector
  //    ));

  //  std::vector< std::string > volumeSeries;
  //  const std::string volumeSeriesTag = "0x0018|0x1030";
  //  volumeSeries.push_back(cbica::GetDICOMTagValue(allFilesInDir[0], volumeSeriesTag));

  //  // looping through all the found files
  //  for (size_t i = 1; i < allFilesInDir.size(); i++)
  //  {
  //    auto temp = cbica::GetDICOMTagValue(allFilesInDir[i], tagToCheck);
  //    auto temp_volSeries = cbica::GetDICOMTagValue(allFilesInDir[i], volumeSeriesTag);

  //    bool newUIDFound = true;
  //    for (size_t j = 0; j < returnVector.size(); j++)
  //    {
  //      if (returnVector[j].first == temp)
  //      {
  //        bool newVolSeriesFound = true;
  //        for (size_t k = 0; k < volumeSeries.size(); k++)
  //        {
  //          if (volumeSeries[k] == temp_volSeries)
  //          {
  //            newVolSeriesFound = false;
  //          }
  //        }
  //        if (!newVolSeriesFound)
  //        {
  //          returnVector[j].second.push_back(allFilesInDir[i]);
  //          newUIDFound = false;
  //          break;
  //        }
  //        else
  //        {
  //          volumeSeries.push_back(temp_volSeries); // the new volume has same series UID information so nothing changes there
  //        }
  //      }
  //    }
  //    if (newUIDFound)
  //    {
  //      // add a new seriesUID-fileNames pair
  //      returnVector.push_back(
  //        std::make_pair(temp, // this is the UID
  //        std::vector< std::string >({ allFilesInDir[i] }) // first filename corresponding to the UID
  //        ));
  //    }
  //  }

  //  return returnVector;

  //  //// this implementation takes a *lot* of time
  //  //auto dicomIO = itk::DCMTKImageIO::New();
  //  //auto inputNames = itk::DCMTKSeriesFileNames::New();
  //  //inputNames->SetInputDirectory(dirName_wrap);
  //  //inputNames->SetLoadPrivateTags(true);
  //  //auto UIDs = inputNames->GetSeriesUIDs(); // this is the primary bottle-neck, I think because it does checks on multiple different things

  //  //return cbica::GetUniqueElements< std::string >(UIDs);
  //}

  /**
  \brief Get the Dicom image reader (not the image, the READER). This is useful for scenarios where reader meta information is needed for later writing step(s).

  Usage:
  \verbatim
  typedef itk::Image< float, 3 > ExpectedImageType;
  std::string inputDirName = parser.getParameterValue("inputDirName");
  auto inputImageReader = GetDicomImageReader< ExpectedImageType >(inputDirName); // reads *all* DICOM images
  auto inputImage = inputImageReader->GetOutput();
  DoAwesomeStuffWithImage( inputImage );
  \endverbatim

  \param dirName This is the directory name of the DICOM image which needs to be loaded - if this is an image, the underlying path of the image is considered
  */
  template <class TImageType = ImageTypeFloat3D >
  typename itk::ImageSeriesReader< TImageType >::Pointer GetDicomImageReader(const std::string &dirName)
  {
    std::string dirName_wrap = cbica::normPath(dirName);
    if (!cbica::isDir(dirName_wrap))
    {
      dirName_wrap = cbica::getFilenamePath(dirName);
    }
    if (dirName_wrap[dirName_wrap.length() - 1] == '/')
      dirName_wrap.pop_back(); // this is done to ensure the last "/" isn't taken into account for file name generation

    //// check read access
    //if (((_access(dirName_wrap.c_str(), 4)) == -1) || ((_access(dirName_wrap.c_str(), 6)) == -1))
    //{
    //  ShowErrorMessage("You don't have read access in selected location. Please check.");
    //  exit(EXIT_FAILURE);
    //}

    auto dicomIO = itk::GDCMImageIO::New();
    auto inputNames = itk::GDCMSeriesFileNames::New();
    inputNames->SetInputDirectory(dirName_wrap);
    inputNames->SetLoadPrivateTags(true);
    auto UIDs = inputNames->GetSeriesUIDs();

    auto UIDs_unique = cbica::GetUniqueElements(UIDs);

    if (UIDs_unique.size() > 1)
    {
      std::cout << "Multiple DICOM series detected.\n";
    }

    inputNames->SetInputDirectory(dirName_wrap);
    //inputNames->SetLoadPrivateTags(true);

    auto filenames = inputNames->GetInputFileNames();

    auto seriesReader = /*typename*/ itk::ImageSeriesReader< TImageType >::New();
    seriesReader->SetImageIO(dicomIO);
    seriesReader->SetFileNames(filenames);

    try
    {
      seriesReader->Update();
    }
    catch (itk::ExceptionObject & err)
    {
      std::cerr << "Error while loading DICOM images: " << err.what() << "\n";
    }

    return seriesReader;
  }

  /**
  \brief Get the itk::Image from input dir name

  Usage:
  \verbatim
  typedef itk::Image< float, 3 > ExpectedImageType;
  std::string inputDirName = parser.getParameterValue("inputDirName");
  ExpectedImageType::Pointer inputImage_1 = ReadDicomImage< ExpectedImageType >(inputFileName); // reads MRI and perfusion data by default tags "0008|0021,0020|0012"
  ExpectedImageType::Pointer inputImage_2 = ReadDicomImage< ExpectedImageType >(inputDirName, "0008|0021")->GetOutput(); // only reads images with tag "0008|0021"
  DoAwesomeStuffWithImage( inputImage );
  \endverbatim

  \param fName name of the image
  \param supportedExtensions Supported extensions
  \return itk::ImageFileReader::Pointer templated over the same as requested by user
  */
  //template <class TImageType = ImageTypeFloat3D >
  //typename TImageType::Pointer ReadDicomImage(const std::string &dirName)
  //{
  //  auto dicomReader = /*typename*/ itk::ImageSeriesReader< TImageType >::New();
  //  if (cbica::isFile(dirName))
  //  {
  //    auto reader =GetImageReader< TImageType >(dirName);
  //    return reader->GetOutput();
  //  }
  //  else
  //  {
  //    //dicomReader = GetDicomImageReader< TImageType >(dirName);
  //  }

  //  // the code below is to ensure that whatever ITK reads aligns with DICOM header information
  //  auto inputDict = (*(dicomReader->GetMetaDataDictionaryArray()))[0];
  //  std::string origin, pixelSpacing, direction, sliceSpacing_1, sliceSpacing_2;
  //  typename TImageType::SpacingType outputSpacing;
  //  typename TImageType::PointType outputOrigin, outputDirection;

  //  itk::ExposeMetaData<std::string>(*inputDict, "0020|0032", origin);
  //  itk::ExposeMetaData<std::string>(*inputDict, "0028|0030", pixelSpacing);
  //  //itk::ExposeMetaData<std::string>(*inputDict, "0020|0037", direction);

  //  if (TImageType::ImageDimension > 2)
  //  {
  //    itk::ExposeMetaData<std::string>(*inputDict, "0018|0050", sliceSpacing_1);
  //    itk::ExposeMetaData<std::string>(*inputDict, "0018|0088", sliceSpacing_2);
  //    if (sliceSpacing_1 == sliceSpacing_2)
  //    {
  //      outputSpacing[2] = static_cast< typename TImageType::PixelType >(std::atof(sliceSpacing_1.c_str()));
  //    }
  //  }

  //  if (!pixelSpacing.empty())
  //  {
  //    auto temp = cbica::stringSplit(pixelSpacing, "\\");
  //    outputSpacing[0] = std::atof(temp[0].c_str());
  //    outputSpacing[1] = std::atof(temp[1].c_str());
  //  }

  //  if (!origin.empty())
  //  {
  //    auto temp = cbica::stringSplit(origin, "\\");
  //    for (unsigned int i = 0; i < TImageType::ImageDimension; i++)
  //    {
  //      outputOrigin[i] = std::atof(temp[i].c_str());
  //    }
  //  }

  //  //if (!direction.empty())
  //  //{
  //  //  auto temp = cbica::stringSplit(direction, "\\");
  //  //  for (auto i = 0; i < TImageType::ImageDimension; i++)
  //  //  {
  //  //    outputOrigin[i] = std::atof(temp[i].c_str());
  //  //  }
  //  //}

  //  auto infoChangeFilter = itk::ChangeInformationImageFilter< TImageType >::New();
  //  infoChangeFilter->SetInput(dicomReader->GetOutput());
  //  infoChangeFilter->SetChangeOrigin(true);
  //  infoChangeFilter->SetChangeSpacing(true);
  //  infoChangeFilter->SetOutputOrigin(outputOrigin);
  //  infoChangeFilter->SetOutputSpacing(outputSpacing);
  //  infoChangeFilter->Update();

  //  return infoChangeFilter->GetOutput();
  //}

  /**
  \brief Get the itk::Image from input dir name

  Usage:
  \verbatim
  typedef itk::Image< float, 3 > ExpectedImageType;
  std::string inputDirName = parser.getParameterValue("inputDirName");
  ExpectedImageType::Pointer inputImage_1 = ReadDicomImage< ExpectedImageType >(inputFileName); // reads MRI and perfusion data by default tags "0008|0021,0020|0012"
  ExpectedImageType::Pointer inputImage_2 = ReadDicomImage< ExpectedImageType >(inputDirName, "0008|0021")->GetOutput(); // only reads images with tag "0008|0021"
  DoAwesomeStuffWithImage( inputImage );
  \endverbatim

  This function calls ReadDicomImage<> internally

  \param fName name of the image
  \param supportedExtensions Supported extensions
  \return itk::ImageFileReader::Pointer templated over the same as requested by user
  */
  //template <class TImageType = ImageTypeFloat3D >
  //typename TImageType::Pointer GetDicomImage(const std::string &dirName)
  //{
  //  return ReadDicomImage< TImageType >(dirName);
  //}


  /**
  \brief Write the itk::Image to the file name

  Usage:
  \verbatim
  typedef itk::Image< float, 3 > ComputedImageType;
  typedef itk::Image< unsigned char, 3 > WrittenImageType;
  ComputedImageType::Pointer imageToWrite = ComputedImageType::New();
  imageToWrite = GetImageSomehow();
  WriteImage< ComputedImageType >(imageToWrite, fileNameToWriteImage); // casts imageToWrite to WrittenImageType
  WriteImage< ComputedImageType, WrittenImageType >(imageToWrite, fileNameToWriteImage);  // writes imageToWrite as ComputedImageType
  // at this point, the image has already been written
  \endverbatim

  \param inputImage Pointer to processed image data which is to be written
  \param fileName File containing the image
  \return itk::Image of specified pixel and dimension type
  */
  template <typename ComputedImageType = ImageTypeFloat3D, typename ExpectedImageType = ComputedImageType>
  void WriteImage(typename ComputedImageType::Pointer imageToWrite, const std::string &fileName)
  {
    //// check write access
    //if (((_access(fileName.c_str(), 2)) == -1) || ((_access(fileName.c_str(), 6)) == -1))
    //{
    //  ShowErrorMessage("You don't have write access in selected location. Please check.");
    //  return;
    //}

    auto filter = /*typename*/ itk::CastImageFilter<ComputedImageType, ExpectedImageType>::New();
    filter->SetInput(imageToWrite);
    filter->Update();

    auto writer = /*typename*/ itk::ImageFileWriter< ExpectedImageType >::New();

    auto ext = cbica::getFilenameExtension(fileName, false);
    if ((ext == ".nii") || (ext == ".nii.gz"))
    {
      writer->SetImageIO(itk::NiftiImageIO::New());
    }

    writer->SetInput(filter->GetOutput());
    writer->SetFileName(fileName);

    try
    {
      writer->Write();
    }
    catch (itk::ExceptionObject &e)
    {
      std::cerr << "Error occurred while trying to write the image '" << fileName << "': " << e.what() << "\n";
      //exit(EXIT_FAILURE);//TBD all exit(EXIT_FAILURE) should be removed 
    }

    return;
  }

  /*
  \brief Write itk::ImageReader as DICOM to specified directory

  This uses default dictionary created by GDCM::ImageIO

  Usage:
  \verbatim
  typedef itk::Image< float, 3 > ComputedImageType;
  typedef itk::Image< unsigned char, 3 > WrittenImageType;
  itk::ImageSeriesReader< ComputedImageType >::Pointer inputImageReader = GetDicomImageReader< ComputedImageType >(inputDirName);
  ComputedImageType::Pointer imageToWrite = GetImageAfterProcessing( inputImageReader->GetOutput() );
  WriteImage< ComputedImageType, WrittenImageType >(imageToWrite, dirNameToWriteImage); // casts imageToWrite to WrittenImageType
  WriteImage< ComputedImageType >(imageToWrite, dirNameToWriteImage); // writes imageToWrite as ComputedImageType
  // at this point, the image has already been written
  \endverbatim

  \param imageToWrite Pointer to processed image data which is to be written
  \param dirName File containing the image
  \return itk::Image of specified pixel and dimension type
  */
  template <typename ComputedImageType = ImageTypeFloat3D>
  void WriteDicomImage(const typename ComputedImageType::Pointer imageToWrite, const std::string &dirName)
  {
    auto reader = typename itk::ImageSeriesReader< ComputedImageType >::New();
    WriteDicomImage< ComputedImageType >(reader, imageToWrite, dirName);
  }

  /**
  \brief Write the itk::Image as a DICOM to the specified directory

  Usage:
  \verbatim
  typedef itk::Image< float, 3 > ComputedImageType;
  typedef itk::Image< unsigned char, 3 > WrittenImageType;
  itk::ImageSeriesReader< ComputedImageType >::Pointer inputImageReader = GetDicomImageReader< ComputedImageType >(inputDirName);
  ComputedImageType::Pointer imageToWrite = GetImageAfterProcessing( inputImageReader->GetOutput() );
  WriteImage< ComputedImageType, WrittenImageType >(inputImageReader, imageToWrite, dirNameToWriteImage); // casts imageToWrite to WrittenImageType
  WriteImage< ComputedImageType >(inputImageReader, imageToWrite, dirNameToWriteImage); // writes imageToWrite as ComputedImageType
  // at this point, the image has already been written
  \endverbatim

  \param inputImageReader The image reader for DICOM - this is necessary to populate the DICOM dictionary properly
  \param imageToWrite Pointer to processed image data which is to be written
  \param dirName File containing the image
  \return itk::Image of specified pixel and dimension type
  */
  template <typename ComputedImageType = ImageTypeFloat3D>
  void WriteDicomImage(const typename itk::ImageSeriesReader< ComputedImageType >::Pointer inputImageReader, const typename ComputedImageType::Pointer imageToWrite, const std::string &dirName)
  {
    if (!cbica::isDir(dirName))
    {
      std::cout << "Specified directory wasn't found, creating...\n";
      cbica::createDir(dirName);
    }

    // check write access
    //if (((_access(dirName.c_str(), 2)) == -1) || ((_access(dirName.c_str(), 6)) == -1))
    //{
    //  ShowErrorMessage("You don't have write access in selected location. Please check.");
    //  return;
    //}

    using ExpectedImageType = itk::Image< short, ComputedImageType::ImageDimension >; // this is needed because DICOM currently only supports short/int
    typedef itk::CastImageFilter<ComputedImageType, ExpectedImageType> CastFilterType;
    typename CastFilterType::Pointer castFilter = CastFilterType::New();
    castFilter->SetInput(imageToWrite);
    castFilter->Update();

    //  typedef typename ExpectedImageType::PixelType DicomPixelType;

    auto dicomIO = itk::GDCMImageIO::New();
    //auto dicomIO = MyGDCMImageIO::New();
    dicomIO->SetComponentType(itk::ImageIOBase::IOComponentType::SHORT);

    auto seriesWriter = itk::ImageSeriesWriter< ExpectedImageType, itk::Image<typename ExpectedImageType::PixelType, 2> >::New();

    auto namesGenerator = itk::NumericSeriesFileNames::New();
    //namesGenerator->SetUseSeriesDetails(false);
    auto start = imageToWrite->GetLargestPossibleRegion().GetIndex();
    auto size = imageToWrite->GetLargestPossibleRegion().GetSize();
    namesGenerator->SetSeriesFormat((dirName + "/image%03d.dcm").c_str());
    namesGenerator->SetStartIndex(start[2]);
    namesGenerator->SetEndIndex(start[2] + size[2] - 1);
    namesGenerator->SetIncrementIndex(1);

    seriesWriter->SetInput(castFilter->GetOutput());
    seriesWriter->SetImageIO(dicomIO);
    seriesWriter->SetFileNames(namesGenerator->GetFileNames());

    typename itk::ImageSeriesReader< ComputedImageType >::DictionaryArrayType outputArray;
    if (inputImageReader.IsNull() || (inputImageReader->GetImageIO() == NULL))
    {
      //dicomIO->SetOrigin(0, imageToWrite->GetOrigin()[0]);
      //dicomIO->SetOrigin(1, imageToWrite->GetOrigin()[1]);
      //dicomIO->SetOrigin(2, imageToWrite->GetOrigin()[2]);
      //dicomIO->SetSpacing(0, imageToWrite->GetSpacing()[0]);
      //dicomIO->SetSpacing(1, imageToWrite->GetSpacing()[1]);
      //dicomIO->SetSpacing(2, imageToWrite->GetSpacing()[2]);
      //dicomIO->SetDimensions(0, imageToWrite->GetLargestPossibleRegion().GetSize()[0]);
      //dicomIO->SetDimensions(1, imageToWrite->GetLargestPossibleRegion().GetSize()[1]);
      //dicomIO->SetDimensions(2, imageToWrite->GetLargestPossibleRegion().GetSize()[2]);

      typename ExpectedImageType::IndexType index;
      index[0] = 0;
      index[1] = 0;
      for (size_t i = 0; i < imageToWrite->GetLargestPossibleRegion().GetSize()[2]; i++)
      {
        auto dict = new typename itk::ImageSeriesReader< ComputedImageType >::DictionaryType;
        typename ExpectedImageType::PointType position;
        index[2] = i;
        imageToWrite->TransformIndexToPhysicalPoint(index, position);
        itk::EncapsulateMetaData<std::string>(*dict, "0020|0032", std::to_string(position[0]) + "\\" + std::to_string(position[1]) + "\\" + std::to_string(position[2])); // patient position
        itk::EncapsulateMetaData<std::string>(*dict, "0020|1041", std::to_string(position[0]) + "\\" + std::to_string(position[1]) + "\\" + std::to_string(position[2])); // slice location
        //itk::EncapsulateMetaData<std::string>(*dict, "0020|0011", std::to_string(1)); 
        //itk::EncapsulateMetaData<std::string>(*dict, "0020|0013", std::to_string(i)); 
        //itk::EncapsulateMetaData<std::string>(*dict, "0018|5100", std::to_string(position[0]) + "\\" + std::to_string(position[1]) + "\\" + std::to_string(position[2]));
        //itk::EncapsulateMetaData<std::string>(*dict, "2020|0010", std::to_string(position[0]) + "\\" + std::to_string(position[1]) + "\\" + std::to_string(position[2]));
        //itk::EncapsulateMetaData<std::string>(*dict, "0018|5101", std::to_string(position[0]) + "\\" + std::to_string(position[1]) + "\\" + std::to_string(position[2]));
        // direction
        //if (ComputedImageType::ImageDimension == 2)
        //{
        //  itk::EncapsulateMetaData<std::string>(*dict, "0020|0037", std::to_string(*imageToWrite->GetDirection()[0]) + "\\" + std::to_string(*imageToWrite->GetDirection()[1]) + "\\0\\" + std::to_string(*imageToWrite->GetDirection()[2]) + "\\" + std::to_string(*imageToWrite->GetDirection()[3]) + "\\0"); // orientation
        //}
        //else if (ComputedImageType::ImageDimension == 3)
        //{
        //  itk::EncapsulateMetaData<std::string>(*dict, "0020|0037", 
        //    std::to_string(*imageToWrite->GetDirection()[0]) + "\\" + std::to_string(*imageToWrite->GetDirection()[1]) + "\\" + std::to_string(*imageToWrite->GetDirection()[2]) + "\\" + 
        //    std::to_string(*imageToWrite->GetDirection()[3]) + "\\" + std::to_string(*imageToWrite->GetDirection()[4]) + "\\" + std::to_string(*imageToWrite->GetDirection()[5]) + "\\" +
        //    std::to_string(*imageToWrite->GetDirection()[6]) + "\\" + std::to_string(*imageToWrite->GetDirection()[7]) + "\\" + std::to_string(*imageToWrite->GetDirection()[8])
        //    ); // orientation
        //}
        itk::EncapsulateMetaData<std::string>(*dict, "0018|0050", std::to_string(imageToWrite->GetSpacing()[2])); // Slice Thickness
        itk::EncapsulateMetaData<std::string>(*dict, "0018|0088", std::to_string(imageToWrite->GetSpacing()[2])); // Spacing Between Slices
        itk::EncapsulateMetaData<std::string>(*dict, "0028|0030", std::to_string(imageToWrite->GetSpacing()[0]) + "\\" + std::to_string(imageToWrite->GetSpacing()[1]));
        //itk::EncapsulateMetaData<std::string>(*dict, "0008|0008", "DERIVED\\SECONDARY"); // Image Type
        //itk::EncapsulateMetaData<std::string>(*dict, "0008|0064", "DV"); // Conversion Type
        //itk::EncapsulateMetaData<std::string>(*dict, "0008|0060", "MR"); // Modality - can never gurantee MR
        //itk::EncapsulateMetaData<std::string>(*dict, "0018|0088", std::to_string(imageToWrite->GetSpacing()[2]));

        outputArray.push_back(dict);
      }

      seriesWriter->SetMetaDataDictionaryArray(&outputArray);
    }
    else
    {
      dicomIO->SetMetaDataDictionary(inputImageReader->GetMetaDataDictionary());
      seriesWriter->SetMetaDataDictionaryArray(inputImageReader->GetMetaDataDictionaryArray()); // no dictionary information present without seriesReader
    }

    try
    {
      seriesWriter->Write();
    }
    catch (itk::ExceptionObject &e)
    {
      std::cerr << "Error occurred while trying to write the image '" << dirName << "': " << e.what() << "\n";
      exit(EXIT_FAILURE);
    }

  }


  /**
  \brief Get the itk::Image from input file name

  Usage:
  \verbatim
  typedef itk::Image< float, 3 > ExpectedImageType;
  std::string inputFileName = parser.getParameterValue("inputImage");
  ExpectedImageType::Pointer inputImage_1 = ReadImage< ExpectedImageType >(inputFileName);
  ExpectedImageType::Pointer inputImage_2 = ReadImage< ExpectedImageType >(inputFileName, ".nii.gz,.img");
  DoAwesomeStuffWithImage( inputImage );
  \endverbatim

  \param fName name of the image
  \param supportedExtensions Supported extensions, defaults to ".nii.gz,.nii"
  \return itk::ImageFileReader::Pointer templated over the same as requested by user
  */
  template <class TImageType = ImageTypeFloat3D >
  typename TImageType::Pointer ReadImage(const std::string &fName, const std::string &supportedExtensions = ".nii.gz,.nii,.dcm", const std::string &delimitor = ",")
  {
    if (cbica::IsDicom(fName))
    {
      DicomIOManager< TImageType > dcmSeriesReader;
      dcmSeriesReader.SetDirectoryPath(fName);
      bool loadstatus = dcmSeriesReader.LoadDicom();
      if (!loadstatus)
      {
        //QMessageBox::critical(this, "Dicom Loading", "Dicom Load Failed");
        return nullptr;
      }
      return dcmSeriesReader.GetITKImage();
    }
    else
    {
      return GetImageReader< TImageType >(fName, supportedExtensions, delimitor)->GetOutput();
    }
  }

  /**
  \brief This is an inline function used to correct the orientation for correct visualization

  \param inputImage The input image
  \return TImageType::Pointer templated over the same as requested by user
  */
  template< class TImageType >
  inline typename TImageType::Pointer GetImageWithOrientFix(const typename TImageType::Pointer inputImage)
  {
    auto orienter = itk::OrientImageFilter<TImageType, TImageType>::New();
    orienter->UseImageDirectionOn();
    orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
    orienter->SetInput(inputImage);
    orienter->Update();

    return orienter->GetOutput();
  }

  /**
  \brief The reads the image according to the appropriate extension and outputs the result in ITK's RAI orientation for visualization

  Usage:
  \verbatim
  using ExpectedImageType = itk::Image< float, 3 >;
  std::string inputFileName = parser.getParameterValue("inputImage");
  auto inputImage_1 = ReadImageWithOrientFix< ExpectedImageType >(inputFileName);
  auto inputImage_2 = ReadImageWithOrientFix< ExpectedImageType >(inputFileName, ".nii.gz,.img");
  DoAwesomeStuffWithImage( inputImage );
  \endverbatim

  \param fName File name of the image
  \param supportedExtensions Supported extensions, defaults to ".nii.gz,.nii"
  \return TImageType::Pointer templated over the same as requested by user
  */
  template< class TImageType >
  typename TImageType::Pointer ReadImageWithOrientFix(const std::string &fName, const std::string &supportedExtensions = ".nii.gz,.nii", const std::string &delimitor = ",")
  {
    std::string extension = cbica::getFilenameExtension(fName);
    std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
    /*   if (cbica::isDir(fName) || (extension == ".dcm") || (extension == ".dicom"))
    {
    return GetImageWithOrientFix<TImageType>(GetDicomImage< TImageType >(fName));
    }
    else
    {*/
    return GetImageWithOrientFix<TImageType>(GetImageReader< TImageType >(fName, supportedExtensions, delimitor)->GetOutput());
    //   }
  }

  /**
  \brief Get the itk::Image from input file name

  Usage:
  \verbatim
  using ExpectedImageType = itk::Image< float, 3 >;
  std::string inputFileName = parser.getParameterValue("inputImage");
  auto inputImage_1 = cbica::ReadImage< ExpectedImageType >(inputFileName);
  auto inputImage_2 = cbica::ReadImage< ExpectedImageType >(inputFileName, ".nii.gz,.img");
  DoAwesomeStuffWithImage( inputImage );
  \endverbatim

  \param fName File name of the image
  \param supportedExtensions Supported extensions, defaults to ".nii.gz,.nii"
  \return itk::ImageFileReader::Pointer templated over the same as requested by user
  */
  template <class TImageType = ImageTypeFloat3D >
  typename TImageType::Pointer GetImage(const std::string &fName, const std::string &supportedExtensions = ".nii.gz,.nii", const std::string &delimitor = ",")
  {
    return ReadImage< TImageType >(fName, supportedExtensions, delimitor);
  }

}
