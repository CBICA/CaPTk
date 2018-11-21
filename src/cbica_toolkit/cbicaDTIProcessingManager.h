/**
\file  cbicaDTIProcessingManager.h

\brief File that holds the DTIProcessingManager class

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

#pragma once

#include "string.h"
//#include <type_traits>
#include <algorithm>
#include "itkImageIOBase.h"

#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageFileReader.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkExtractImageFilter.h"
#include "itkVectorImage.h"
#include "itkDiffusionTensor3DReconstructionImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkDiffusionTensor3D.h"
#include "gdcmDictEntry.h"
#include "gdcmDict.h"           // access to dictionary
#include "gdcmFile.h"           // access to dictionary
#include "gdcmDictEntry.h"      // access to dictionary
#include "gdcmGlobal.h"         // access to dictionary
#include "gdcmElement.h"
#include "gdcmPrivateTag.h"
#include "gdcmStringFilter.h"
#include "cbicaUtilities.h"
#include "itkImageMaskSpatialObject.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkNumericTraits.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkDiffusionTensor3D.h"
#include "itkTensorFractionalAnisotropyImageFilter.h"
#include "itkTensorRelativeAnisotropyImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"

#include "itkNrrdImageIO.h"

//#include "CAPTk.h"

const int Dimensions = 3;
typedef short PixelValueType;
typedef itk::Image<PixelValueType, Dimensions> VolumeType;
typedef itk::ImageSeriesReader<VolumeType> ReaderType;
typedef itk::GDCMImageIO ImageIOType;
typedef itk::GDCMSeriesFileNames InputNamesGeneratorType;
typedef itk::MetaDataDictionary DictionaryType;
typedef itk::VectorImage<PixelValueType, Dimensions> VectorImageType;
using ImageTypeScalar3D = itk::Image< float, 3 >;

// relevant GE private tags
const gdcm::DictEntry GEDictBValue("0x0043", "0x1039", gdcm::VR::IS, gdcm::VM::VM1, "B Value of diffusion weighting");
const gdcm::DictEntry GEDictXGradient("0x0019", "0x10bb", gdcm::VR::DS, gdcm::VM::VM1, "X component of gradient direction");
const gdcm::DictEntry GEDictYGradient("0x0019", "0x10bc", gdcm::VR::DS, gdcm::VM::VM1, "Y component of gradient direction");
const gdcm::DictEntry GEDictZGradient("0x0019", "0x10bd", gdcm::VR::DS, gdcm::VM::VM1, "Z component of gradient direction");

// relevant Siemens private tags
const gdcm::DictEntry SiemensMosiacParameters("0x0051", "0x100b", gdcm::VR::IS, gdcm::VM::VM1, "Mosiac Matrix Size");
const gdcm::DictEntry SiemensDictNMosiac("0x0019", "0x100a", gdcm::VR::US, gdcm::VM::VM1, "Number of Images In Mosaic");
const gdcm::DictEntry SiemensDictBValue("0x0019", "0x100c", gdcm::VR::IS, gdcm::VM::VM1, "B Value of diffusion weighting");
const gdcm::DictEntry SiemensDictDiffusionDirection("0x0019", "0x100e", gdcm::VR::FD, gdcm::VM::VM3, "Diffusion Gradient Direction");
const gdcm::DictEntry SiemensDictDiffusionMatrix("0x0019", "0x1027", gdcm::VR::FD, gdcm::VM::VM6, "Diffusion Matrix");


namespace cbica
{
  /**
  \class DTIProcessingManager

  \brief A small description of the class

  A detailed description.
  }
  */
  class DTIProcessingManager
  {
  public:
    DTIProcessingManager();
    ~DTIProcessingManager();
    template<typename ImageType, typename TensorImageType>
    typename ImageType::Pointer allocateImage(typename TensorImageType::Pointer tenIm);

    void GetImageInfo(std::string fName, itk::ImageIOBase::IOPixelType *pixelType, itk::ImageIOBase::IOComponentType *componentType);

    int GetNumberOfVolumes(VolumeType::Pointer rawVol, int nVolume, int nSliceInVolume);

    template < typename TInputPixelType, typename TMaskPixelType, typename TOutputTensorCompType>
    std::vector<ImageTypeScalar3D::Pointer> dtiRecon(VectorImageType::Pointer inputImage, std::string maskFile, int verbose, int inputIsVectorImage, std::vector< vnl_vector_fixed<double, 3> > diffuionVector, double bValue);

    std::vector<ImageTypeScalar3D::Pointer> ConvertDWIToScalars(std::string inputDirName, std::string maskFileName);


  };


  template<typename ImageType, typename TensorImageType>
  typename ImageType::Pointer DTIProcessingManager::allocateImage(typename TensorImageType::Pointer tenIm)
  {
    typename ImageType::Pointer outputImage = ImageType::New();

    outputImage->SetOrigin(tenIm->GetOrigin());
    outputImage->SetSpacing(tenIm->GetSpacing());
    outputImage->SetDirection(tenIm->GetDirection());
    outputImage->SetLargestPossibleRegion(tenIm->GetLargestPossibleRegion());
    outputImage->SetRequestedRegion(tenIm->GetRequestedRegion());
    outputImage->SetBufferedRegion(tenIm->GetBufferedRegion());
    outputImage->Allocate();

    return outputImage;
  }
  //--------------------------------------------------------------
  //----------------------------------------------------------
  template < typename TInputPixelType, typename TMaskPixelType, typename TOutputTensorCompType>
  std::vector<itk::Image<ImageTypeScalar3D::PixelType, Dimensions>::Pointer> DTIProcessingManager::dtiRecon(VectorImageType::Pointer inputImage, std::string maskFile, int verbose, int inputIsVectorImage, std::vector< vnl_vector_fixed<double, 3> > DiffusionVectorWrite, double bValue)
  {
    //bool readb0 = false;
    //double b0 = 0;
    typedef itk::VectorImage<TInputPixelType, Dimensions> GradientImageType;
    //typedef itk::Image<TInputPixelType, 4>   InputImageType;
    std::cout << "Converting DWI to vector image\n";

    //Set up the gradient image size
    typename GradientImageType::Pointer gradIm = inputImage;
    typedef itk::DiffusionTensor3DReconstructionImageFilter< TInputPixelType, TInputPixelType, TOutputTensorCompType > TensorReconstructionImageFilterType;
    // -------------------------------------------------------------------------

    //int NumberOfGradients = DiffusionVectorWrite.size();
    //double bValue = 0;
    //std::ifstream bvalIn(bvalFile.c_str());
    //std::string line;
    //while (!bvalIn.eof())
    //{
    //  getline(bvalIn, line);
    //  double val;
    //  std::stringstream ss(line);

    //  while (ss >> val)
    //  {
    //    ++NumberOfGradients;
    //    if (val != 0 && bValue == 0)
    //      bValue = val;
    //    else if (val != 0 && bValue != val)
    //    {
    //      std::cerr << "multiple bvalues not allowed" << std::endl;
    //      //return EXIT_FAILURE;
    //    }
    //  }
    //}
    //std::ifstream bvecIn(gradFile.c_str());
    //getline(bvecIn, line);
    //std::stringstream Xss(line);
    //getline(bvecIn, line);
    //std::stringstream Yss(line);
    //getline(bvecIn, line);
    //std::stringstream Zss(line);

    typename TensorReconstructionImageFilterType::GradientDirectionType vect3d;
    typename TensorReconstructionImageFilterType::GradientDirectionContainerType::Pointer DiffusionVectors = TensorReconstructionImageFilterType::GradientDirectionContainerType::New();

    for (unsigned int index1 = 0; index1 < DiffusionVectorWrite.size(); index1++)
    {
      vnl_vector_fixed<double, 3> val = DiffusionVectorWrite[index1];
      vect3d[0] = val[0];
      vect3d[1] = val[1];
      vect3d[2] = val[2];
      DiffusionVectors->InsertElement(index1, vect3d);
    }

    //int counter = 0;
    //double x, y, z;
    //while (Xss >> x)
    //{
    //  Yss >> y;
    //  Zss >> z;
    //  vect3d[0] = x; vect3d[1] = y; vect3d[2] = z;
    //  DiffusionVectors->InsertElement(counter, vect3d);
    //  ++counter;
    //}
    // -------------------------------------------------------------------------
    typename TensorReconstructionImageFilterType::Pointer tensorReconstructionFilter = TensorReconstructionImageFilterType::New();

    typedef itk::Image<TMaskPixelType, Dimensions >   ImageMaskType;
    typedef itk::ImageMaskSpatialObject< Dimensions >   MaskType;
    typename MaskType::Pointer  spatialObjectMask = MaskType::New();
    typedef itk::ImageFileReader<ImageMaskType>    MaskReaderType;
    typename MaskReaderType::Pointer  maskReader = MaskReaderType::New();
    typedef itk::BinaryThresholdImageFilter<ImageMaskType, typename MaskType::ImageType> ThresholderType;
    typename ThresholderType::Pointer thresholder = ThresholderType::New();
    thresholder->SetOutsideValue(itk::NumericTraits< MaskType::ImageType::PixelType>::Zero);
    thresholder->SetInsideValue(itk::NumericTraits< MaskType::ImageType::PixelType>::One);
    thresholder->SetLowerThreshold(itk::NumericTraits<unsigned short>::One);
    thresholder->InPlaceOn();

    maskReader->SetFileName(maskFile);
    try
    {
      maskReader->Update();
    }
    catch (itk::ExceptionObject & err)
    {
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      //return EXIT_FAILURE;
    }
    thresholder->SetInput(maskReader->GetOutput());
    thresholder->Update();
    spatialObjectMask->SetImage(thresholder->GetOutput());
    //tensorReconstructionFilter->SetMaskSpatialObject(spatialObjectMask);

    //---------------------------------------------------------------------------------
    tensorReconstructionFilter->SetGradientImage(DiffusionVectors, gradIm);
    //tensorReconstructionFilter->SetNumberOfThreads(1);
    tensorReconstructionFilter->SetBValue(bValue);
    tensorReconstructionFilter->Update();

    //--------------------------computing scalars----------------------------
    std::string default_ext = ".nii.gz";
    // perhaps the following should be replaced with bools
    static int writeFA = 1;
    static int writeTR = 1;
    static int writeEign = 0;
    static int writeGeo = 0;
    static int writeGordR = 0;
    static int writeGordK = 0;
    static int writeRadAx = 1;
    static int writeSkew = 0;
    static int writeKurt = 0;

    typedef float ScalarPixelType;
    typedef itk::DiffusionTensor3D< float > TensorPixelType;
    typedef itk::Image< float, Dimensions > ScalarImageType;
    typedef itk::VectorImage<float, Dimensions > VectorImageType;
    typedef itk::Image< TensorPixelType, Dimensions > TensorImageType;

    std::vector<ScalarImageType::Pointer> vectorOfDTIScalars;

    typedef itk::ImageRegionConstIteratorWithIndex < TensorImageType > ConstIterType;

    try
    {
      TensorImageType::Pointer tensorIm = TensorImageType::New();
      tensorIm = tensorReconstructionFilter->GetOutput();
      // Allocate each image that we will want to make..
      //FA
      ScalarImageType::Pointer faIm = ScalarImageType::New();

      //TR
      ScalarImageType::Pointer trIm = ScalarImageType::New();

      //Eigensys
      ScalarImageType::Pointer l1Im = ScalarImageType::New();
      ScalarImageType::Pointer l2Im = ScalarImageType::New();
      ScalarImageType::Pointer l3Im = ScalarImageType::New();
      VectorImageType::Pointer v1Im = VectorImageType::New();
      VectorImageType::Pointer v2Im = VectorImageType::New();
      VectorImageType::Pointer v3Im = VectorImageType::New();

      //Skewness & Kurtosis
      ScalarImageType::Pointer skIm = ScalarImageType::New();
      ScalarImageType::Pointer kuIm = ScalarImageType::New();

      //Geometric Features
      ScalarImageType::Pointer clIm = ScalarImageType::New();
      ScalarImageType::Pointer cpIm = ScalarImageType::New();
      ScalarImageType::Pointer csIm = ScalarImageType::New();

      //Radial Axial
      ScalarImageType::Pointer rdIm = ScalarImageType::New();
      ScalarImageType::Pointer adIm = ScalarImageType::New();

      //Gordons R features
      ScalarImageType::Pointer r1Im = ScalarImageType::New();
      ScalarImageType::Pointer r2Im = ScalarImageType::New();
      ScalarImageType::Pointer r3Im = ScalarImageType::New();

      //Gordons K features
      ScalarImageType::Pointer k1Im = ScalarImageType::New();
      ScalarImageType::Pointer k2Im = ScalarImageType::New();
      ScalarImageType::Pointer k3Im = ScalarImageType::New();

      //Allocate all the images...
      if (writeFA)
        faIm = allocateImage<ScalarImageType, TensorImageType>(tensorIm);

      if (writeTR)
        trIm = allocateImage<ScalarImageType, TensorImageType>(tensorIm);

      if (writeEign)
      {
        l1Im = allocateImage<ScalarImageType, TensorImageType>(tensorIm);
        l2Im = allocateImage<ScalarImageType, TensorImageType>(tensorIm);
        l3Im = allocateImage<ScalarImageType, TensorImageType>(tensorIm);
        v1Im->SetVectorLength(Dimensions);
        v2Im->SetVectorLength(Dimensions);
        v3Im->SetVectorLength(Dimensions);
        v1Im = allocateImage<VectorImageType, TensorImageType>(tensorIm);
        v2Im = allocateImage<VectorImageType, TensorImageType>(tensorIm);
        v3Im = allocateImage<VectorImageType, TensorImageType>(tensorIm);
      }

      if (writeSkew)
        skIm = allocateImage<ScalarImageType, TensorImageType>(tensorIm);

      if (writeKurt)
        kuIm = allocateImage<ScalarImageType, TensorImageType>(tensorIm);

      if (writeGeo)
      {
        clIm = allocateImage<ScalarImageType, TensorImageType>(tensorIm);
        cpIm = allocateImage<ScalarImageType, TensorImageType>(tensorIm);
        csIm = allocateImage<ScalarImageType, TensorImageType>(tensorIm);
      }
      if (writeRadAx)
      {
        rdIm = allocateImage<ScalarImageType, TensorImageType>(tensorIm);
        adIm = allocateImage<ScalarImageType, TensorImageType>(tensorIm);
      }
      if (writeGordR)
      {
        r1Im = allocateImage<ScalarImageType, TensorImageType>(tensorIm);
        r2Im = allocateImage<ScalarImageType, TensorImageType>(tensorIm);
        r3Im = allocateImage<ScalarImageType, TensorImageType>(tensorIm);
      }
      if (writeGordK)
      {
        k1Im = allocateImage<ScalarImageType, TensorImageType>(tensorIm);
        k2Im = allocateImage<ScalarImageType, TensorImageType>(tensorIm);
        k3Im = allocateImage<ScalarImageType, TensorImageType>(tensorIm);
      }
      //Loop though all the voxels and if compute the needed measures!
      ConstIterType iter(tensorIm, tensorIm->GetLargestPossibleRegion());
      for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter)
      {
        TensorPixelType tmp = iter.Get();
        typename TensorPixelType::EigenValuesArrayType     lambda;
        typename TensorPixelType::EigenVectorsMatrixType   vMat;
        tmp.ComputeEigenAnalysis(lambda, vMat);

        typename ScalarImageType::IndexType index = iter.GetIndex();

        if (writeTR)
          trIm->SetPixel(index, lambda[0] + lambda[1] + lambda[2]);

        if (writeFA)
        {
          faIm->SetPixel(index, tmp.GetFractionalAnisotropy());
        }

        if (writeEign)
        {
          l1Im->SetPixel(index, lambda[2]);
          l2Im->SetPixel(index, lambda[1]);
          l3Im->SetPixel(index, lambda[0]);

          typename VectorImageType::PixelType vec1(3);
          typename VectorImageType::PixelType vec2(3);
          typename VectorImageType::PixelType vec3(3);

          vec1[0] = vMat[2][0];
          vec1[1] = vMat[2][1];
          vec1[2] = vMat[2][2];

          vec2[0] = vMat[1][0];
          vec2[1] = vMat[1][1];
          vec2[2] = vMat[1][2];

          vec3[0] = vMat[0][0];
          vec3[1] = vMat[0][1];
          vec3[2] = vMat[0][2];

          v1Im->SetPixel(index, vec1);
          v2Im->SetPixel(index, vec2);
          v3Im->SetPixel(index, vec3);

        }

        if (writeSkew)
        {
          ScalarPixelType m1, m3, l1, l2, l3;
          l1 = abs(lambda[0]);
          l2 = abs(lambda[1]);
          l3 = abs(lambda[2]);
          m1 = (l1 + l2 + l3) / 3.0;

          if (m1 > 0)
          {
            m3 = (std::pow(l1 - m1, 3) + std::pow(l2 - m1, 3) + std::pow(l3 - m1, 3)) / (std::pow(l1, 3) + std::pow(l2, 3) + std::pow(l3, 3));
            if (m3 > 0)
            {
              skIm->SetPixel(index, std::pow(m3, static_cast<ScalarPixelType>(1.0 / 3.0)));
            }
            else
            {
              skIm->SetPixel(index, -1 * std::pow((-1 * m3), static_cast<ScalarPixelType>(1.0 / 3.0)));
            }
          }
          else
          {
            skIm->SetPixel(index, static_cast<ScalarPixelType>(0));
          }
        }

        if (writeKurt)
        {
          ScalarPixelType m1, m4, l1, l2, l3;
          l1 = abs(lambda[0]);
          l2 = abs(lambda[1]);
          l3 = abs(lambda[2]);
          m1 = (l1 + l2 + l3) / 3.0;
          if (m1 > 0)
          {
            m4 = (std::pow(l1 - m1, 4) + std::pow(l2 - m1, 4) + std::pow(l3 - m1, 4)) / (std::pow(l1, 4) + std::pow(l2, 4) + std::pow(l3, 4));
            kuIm->SetPixel(index, std::pow(m4, static_cast<ScalarPixelType>(1.0 / 4.0)));
          }
          else
          {
            kuIm->SetPixel(index, static_cast<ScalarPixelType>(0));
          }
        }

        if (writeGeo)
        {
          if (lambda[2] > 0)
          {
            clIm->SetPixel(index, (lambda[2] - lambda[1]) / lambda[2]);
            cpIm->SetPixel(index, (lambda[1] - lambda[0]) / lambda[2]);
            csIm->SetPixel(index, lambda[0] / lambda[2]);
          }
          else
          {
            clIm->SetPixel(index, 0);
            cpIm->SetPixel(index, 0);
            csIm->SetPixel(index, 0);
          }
        }

        if (writeRadAx)
        {
          rdIm->SetPixel(index, (lambda[1] + lambda[0]) / 2);
          adIm->SetPixel(index, lambda[2]);
        }

        if (writeGordR)
        {
          //Compute the moments...
          ScalarPixelType m1, m2, m3;
          ScalarPixelType r1, r2, r3;
          m1 = (lambda[0] + lambda[1] + lambda[2]) / 3.0;
          m2 = (std::pow(lambda[0] - m1, 2) + std::pow(lambda[1] - m1, 2)
            + std::pow(lambda[2] - m1, 2)) / 3.0;
          m3 = (std::pow(lambda[0] - m1, 3) + std::pow(lambda[1] - m1, 3)
            + std::pow(lambda[2] - m1, 3)) / 3.0;

          r1 = sqrt(3 * (std::pow(m1, 2) + m2));
          r2 = sqrt(3 * m2 / 2 / (std::pow(m1, 2) + m2));
          //        r3 = sqrt(2) * m3 / std::pow(static_cast<double>(m2), static_cast<double>(1.5));
          r3 = (lambda[0] * lambda[1] * lambda[2]) / std::pow(static_cast<double>(sqrt(3 * m2)), 3);

          r1Im->SetPixel(index, r1);
          r2Im->SetPixel(index, r2);
          r3Im->SetPixel(index, r3);
        }

        if (writeGordK)
        {
          //Compute the moments...
          ScalarPixelType m1, m2, m3;
          ScalarPixelType k1, k2, k3;

          m1 = (lambda[0] + lambda[1] + lambda[2]) / 3.0;
          m2 = (std::pow(lambda[0] - m1, 2) + std::pow(lambda[1] - m1, 2)
            + std::pow(lambda[2] - m1, 2)) / 3.0;
          m3 = (std::pow(lambda[0] - m1, 3) + std::pow(lambda[1] - m1, 3)
            + std::pow(lambda[2] - m1, 3)) / 3.0;

          k1 = 3 * m1;
          k2 = sqrt(3 * m2);
          //        k3 = sqrt(2) * m3 / std::pow(static_cast<double>(m2), static_cast<double>(1.5));
          k3 = (lambda[0] * lambda[1] * lambda[2]) / std::pow(static_cast<double>(sqrt(3 * m2)), 3);


          k1Im->SetPixel(index, k1);
          k2Im->SetPixel(index, k2);
          k3Im->SetPixel(index, k3);
        }
      }

      std::cout << "Done Computing Scalars\n";

      //typedef itk::ImageFileWriter< ScalarImageType >  ScalarWriterType;
      //typedef itk::ImageFileWriter< VectorImageType >  VectorWriterType;

      if (writeFA)
        vectorOfDTIScalars.push_back(faIm);

      if (writeTR)
        vectorOfDTIScalars.push_back(trIm);

      if (writeRadAx)
      {
        vectorOfDTIScalars.push_back(rdIm);
        vectorOfDTIScalars.push_back(adIm);
      }
    }
    catch (itk::ExceptionObject & excp)
    {
      std::cerr << "Exception caught - " << excp.what() << "\n";
    }
    return vectorOfDTIScalars;
  }
  //------------------------------------------------------------------



}
