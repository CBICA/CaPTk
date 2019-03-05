#include "DicomImageReader.h"
//#include <QString>
//#include <QDir>
//#include <QDirIterator>
//#include <QTime>

DicomImageReader::DicomImageReader()
{
}

DicomImageReader::~DicomImageReader()
{
}

void DicomImageReader::SetDirectoryPath(std::string path)
{
  this->m_dir = path;
}

//bool DicomImageReader::LoadDicom()
//{
//  bool loadStatus = false;

//  NamesGeneratorType::Pointer nG = NamesGeneratorType::New();
//  nG->SetDirectory(this->m_dir);
//  std::vector<std::string> files = nG->GetInputFileNames();

//  if (files.empty())
//  {
//    //! no files in the given directory
//    loadStatus = false;
//  }
//  else
//  {
//    //! get the first file in the given directory and check its pixel type
//    //! and component type so we can invoke the load with appropriate
//    //! template type
//    //std::string fname = files[0].absoluteFilePath().toStdString();
//    std::string fname = files[0].c_str();

//    itk::ImageIOBase::Pointer imageIO;
//    itk::ImageIOBase::IOPixelType       pixelType;
//    itk::ImageIOBase::IOComponentType   componentType;
//    unsigned int dimensions;
//    imageIO = itk::ImageIOFactory::CreateImageIO(fname.c_str(), itk::ImageIOFactory::ReadMode);
//    bool canRead = imageIO->CanReadFile(fname.c_str());
//    if (canRead)
//    {
//      imageIO->SetFileName(fname);
//      imageIO->ReadImageInformation();
//      pixelType = imageIO->GetPixelType();
//      componentType = imageIO->GetComponentType();
//      dimensions = imageIO->GetNumberOfDimensions();

//      if (pixelType == itk::ImageIOBase::SCALAR)
//      {
//        if (dimensions == 3)
//        {
//          switch (componentType)
//          {
//          case itk::ImageIOBase::SHORT:
//          {
//            using PixelType = short;
//            using ImageType = itk::Image<PixelType, 3>;
//            typename ImageType::Pointer img = ImageType::New();
//            bool readStatus = false;
//            img = ReadDicomSeries<ImageType>(readStatus);
//            if (readStatus)
//            {
//              m_image3dfloat = ConvertImage3DToFloatImage3D<ImageType>(img);
//              loadStatus = true;
//            }
//            else
//              loadStatus = false;
//            break;
//          }
//          case itk::ImageIOBase::USHORT:
//          {
//            using PixelType = unsigned short;
//            using ImageType = itk::Image<PixelType, 3>;
//            typename ImageType::Pointer img = ImageType::New();
//            bool readStatus = false;
//            img = ReadDicomSeries<ImageType>(readStatus);
//            if (readStatus)
//            {
//              m_image3dfloat = ConvertImage3DToFloatImage3D<ImageType>(img);
//              loadStatus = true;
//            }
//            else
//              loadStatus = false;
//            break;
//          }
//          case itk::ImageIOBase::INT:
//          {
//            using PType = int;
//            using ImageType = itk::Image<PType, 3>;
//            typename ImageType::Pointer img = ImageType::New();
//            bool readStatus = false;
//            img = ReadDicomSeries<ImageType>(readStatus);
//            if (readStatus)
//            {
//              m_image3dfloat = ConvertImage3DToFloatImage3D<ImageType>(img);
//              loadStatus = true;
//            }
//            else
//              loadStatus = false;
//            break;
//          }
//          case itk::ImageIOBase::CHAR:
//          {
//            //! need this type of data
//            //! this needs to be handled when we get the data
//            loadStatus = false;
//            break;
//          }
//          case itk::ImageIOBase::UCHAR:
//          {
//            //! need this type of data
//            //! this needs to be handled when we get the data
//            loadStatus = false;
//            break;
//          }
//          case itk::ImageIOBase::DOUBLE:
//          {
//            //! need this type of data
//            //! this needs to be handled when we get the data
//            loadStatus = false;
//            break;
//          }
//          case itk::ImageIOBase::FLOAT:
//          {
//            //! need this type of data
//            //! this needs to be handled when we get the data
//            loadStatus = false;
//            break;
//          }
//          default:
//            loadStatus = false;
//            break;
//          }
//        }
//        else
//          loadStatus = false;
//      }
//      else
//        loadStatus = false;
//    }
//    else
//      loadStatus = false;
//  }
//  return loadStatus;
//}

//DicomImageReader::ImageType3DFloat::Pointer DicomImageReader::GetITKImage()
//{
//  return m_image3dfloat;
//}
