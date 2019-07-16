#include "DicomIOManager.h"
#include "itkImageSeriesReader.h"
#include "DicomSeriesReader.h"
#include "DicomImageReader.h"
#include "gdcmReader.h"

template <class T>
DicomIOManager<T>::DicomIOManager()
{

}
template <class T>
DicomIOManager<T>::~DicomIOManager()
{

}
template <class T>
void DicomIOManager<T>::SetDirectoryPath(std::string path)
{
  this->m_dir = path;
}
template <class T>
bool DicomIOManager<T>::LoadDicom()
{
  bool loadStatus = false;

  NamesGeneratorType::Pointer nG = NamesGeneratorType::New();
  nG->SetDirectory(this->m_dir);
  std::vector<std::string> files = nG->GetInputFileNames();

  if (files.empty())
  {
    //! no files in the given directory
    loadStatus = false;
  }
  else
  {
    //! get the first file in the given directory and check its pixel type
    //! and component type so we can invoke the load with appropriate
    //! template type
    //std::string fname = files[0].absoluteFilePath().toStdString();
    std::string fname = files[0].c_str();

    bool isDicom = this->IsDicom(fname);

    itk::ImageIOBase::Pointer imageIO;
    itk::ImageIOBase::IOPixelType       pixelType;
    itk::ImageIOBase::IOComponentType   componentType;
    unsigned int dimensions;
	bool canRead = this->CanReadFile(fname.c_str(),imageIO);
    if (canRead && isDicom)
    {
      imageIO->SetFileName(fname);
      imageIO->ReadImageInformation();
      pixelType = imageIO->GetPixelType();
      componentType = imageIO->GetComponentType();
      dimensions = imageIO->GetNumberOfDimensions();

      if (pixelType == itk::ImageIOBase::SCALAR)
      {
        if (dimensions == 3)
        {
          switch (componentType)
          {
          case itk::ImageIOBase::SHORT:
          {
            using PixelType = short;
            using ImageType = itk::Image<PixelType, T::ImageDimension>;
            typename ImageType::Pointer img = ImageType::New();
            bool readStatus = false;
            DicomSeriesReader *serReader = new DicomSeriesReader();
            serReader->SetDirectoryPath(this->m_dir);
            img = serReader->ReadDicomSeries<ImageType>(readStatus);
            if (readStatus)
            {
              m_image3d = ConvertImage3DToFloatImage3D<ImageType>(img);
              loadStatus = true;
            }
            else
              loadStatus = false;
            delete serReader;
            break;
          }
          case itk::ImageIOBase::USHORT:
          {
            using PixelType = unsigned short;
            using ImageType = itk::Image<PixelType, T::ImageDimension>;
            typename ImageType::Pointer img = ImageType::New();
            bool readStatus = false;
            DicomSeriesReader *serReader = new DicomSeriesReader();
            serReader->SetDirectoryPath(this->m_dir);
            img = serReader->ReadDicomSeries<ImageType>(readStatus);
            if (readStatus)
            {
              m_image3d = ConvertImage3DToFloatImage3D<ImageType>(img);
              loadStatus = true;
            }
            else
              loadStatus = false;
            delete serReader;
            break;
          }
          case itk::ImageIOBase::INT:
          {
            using PType = int;
            using ImageType = itk::Image<PType, T::ImageDimension>;
            typename ImageType::Pointer img = ImageType::New();
            bool readStatus = false;
            DicomSeriesReader *serReader = new DicomSeriesReader();
            serReader->SetDirectoryPath(this->m_dir);
            img = serReader->ReadDicomSeries<ImageType>(readStatus);
            if (readStatus)
            {
              m_image3d = ConvertImage3DToFloatImage3D<ImageType>(img);
              loadStatus = true;
            }
            else
              loadStatus = false;
            delete serReader;
            break;
          }
          case itk::ImageIOBase::UINT:
          case itk::ImageIOBase::LONG:
          case itk::ImageIOBase::LONGLONG:
          case itk::ImageIOBase::ULONG:
          case itk::ImageIOBase::ULONGLONG:
          {
            //! need this type of data
            //! this needs to be handled when we get the data
            loadStatus = false;
            break;
          }
          case itk::ImageIOBase::CHAR:
          {
            //! need this type of data
            //! this needs to be handled when we get the data
            loadStatus = false;
            break;
          }
          case itk::ImageIOBase::UCHAR:
          {
            //! need this type of data
            //! this needs to be handled when we get the data
            loadStatus = false;
            break;
          }
          case itk::ImageIOBase::DOUBLE:
          {
            //! need this type of data
            //! this needs to be handled when we get the data
            loadStatus = false;
            break;
          }
          case itk::ImageIOBase::FLOAT:
          {
            //! need this type of data
            //! this needs to be handled when we get the data
            loadStatus = false;
            break;
          }
          default:
            loadStatus = false;
            break;
          }
        }
        else if (dimensions == 2)
        {
          switch (componentType)
          {
          case itk::ImageIOBase::SHORT:
          {
            using PixelType = short;
            using ImageType = itk::Image<PixelType, T::ImageDimension>;
            typename ImageType::Pointer img = ImageType::New();
            bool readStatus = false;
            DicomImageReader *imgReader = new DicomImageReader();
            imgReader->SetDirectoryPath(this->m_dir);
            img = imgReader->ReadDicomImage<ImageType>(readStatus);
            if (readStatus)
            {
              m_image3d = ConvertImage3DToFloatImage3D<ImageType>(img);
              loadStatus = true;
            }
            else
              loadStatus = false;
            delete imgReader;
            break;
          }
          case itk::ImageIOBase::USHORT:
          {
            using PixelType = unsigned short;
            using ImageType = itk::Image<PixelType, T::ImageDimension>;
            typename ImageType::Pointer img = ImageType::New();
            bool readStatus = false;
            DicomImageReader *imgReader = new DicomImageReader();
            imgReader->SetDirectoryPath(this->m_dir);
            img = imgReader->ReadDicomImage<ImageType>(readStatus);
            if (readStatus)
            {
              m_image3d = ConvertImage3DToFloatImage3D<ImageType>(img);
              loadStatus = true;
            }
            else
              loadStatus = false;
            delete imgReader;
            break;
          }
          case itk::ImageIOBase::INT:
          {
            using PType = int;
            using ImageType = itk::Image<PType, T::ImageDimension>;
            typename ImageType::Pointer img = ImageType::New();
            bool readStatus = false;
            DicomImageReader *imgReader = new DicomImageReader();
            imgReader->SetDirectoryPath(this->m_dir);
            img = imgReader->ReadDicomImage<ImageType>(readStatus);
            if (readStatus)
            {
              m_image3d = ConvertImage3DToFloatImage3D<ImageType>(img);
              loadStatus = true;
            }
            else
              loadStatus = false;
            delete imgReader;
            break;
          }
          case itk::ImageIOBase::UINT:
          case itk::ImageIOBase::LONG:
          case itk::ImageIOBase::LONGLONG:
          case itk::ImageIOBase::ULONG:
          case itk::ImageIOBase::ULONGLONG:
          {
            //! need this type of data
            //! this needs to be handled when we get the data
            loadStatus = false;
            break;
          }
          case itk::ImageIOBase::CHAR:
          {
            //! need this type of data
            //! this needs to be handled when we get the data
            loadStatus = false;
            break;
          }
          case itk::ImageIOBase::UCHAR:
          {
            //! need this type of data
            //! this needs to be handled when we get the data
            loadStatus = false;
            break;
          }
          case itk::ImageIOBase::DOUBLE:
          {
            //! need this type of data
            //! this needs to be handled when we get the data
            loadStatus = false;
            break;
          }
          case itk::ImageIOBase::FLOAT:
          {
            //! need this type of data
            //! this needs to be handled when we get the data
            loadStatus = false;
            break;
          }
          default:
            loadStatus = false;
            break;
          }
        }
        else
          loadStatus = false;
      }
      else
        loadStatus = false;
    }
    else
      loadStatus = false;
  }
  return loadStatus;
}

template<class T>
inline bool DicomIOManager<T>::IsDicom(std::string path)
{
  return cbica::IsDicom(path);
}

template<class T>
inline bool DicomIOManager<T>::CanReadFile(std::string path, itk::ImageIOBase::Pointer &imageIO)
{
	imageIO = itk::ImageIOFactory::CreateImageIO(path.c_str(), itk::ImageIOFactory::ReadMode);
	return imageIO->CanReadFile(path.c_str());
}

template <class T>
template<class TInputImage>
typename T::Pointer DicomIOManager<T>::ConvertImage3DToFloatImage3D(typename TInputImage::Pointer image)
{
  typedef itk::CastImageFilter<TInputImage, T> CastFilterType;
  auto castFilter = CastFilterType::New();
  castFilter->SetInput(image);
  castFilter->Update();
  return castFilter->GetOutput();
}

template <class T>
typename T::Pointer DicomIOManager<T>::GetITKImage()
{
  return m_image3d;
}
