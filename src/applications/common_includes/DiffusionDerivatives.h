/**
\file  DiffusionDerivatives.h

\brief The header file containing the DiffusionDerivatives class, used to calculate PSR, RCBV, and PH
Author: Saima Rathore
Library Dependecies: ITK 4.7+ <br>

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or http://www.med.upenn.edu/sbia/software/license.html

*/


#ifndef _DiffusionDerivatives_h_
#define _DiffusionDerivatives_h_

#include "qmessagebox.h"
#include "qstring.h"
#include "cbicaUtilities.h"
#include "ApplicationBase.h"
#include "FeatureReductionClass.h"
#include "CAPTk.h"
#include <itkExtractImageFilter.h>

/**
\class DiffusionDerivatives

\brief Calculates Diffusion Derivatives

Reference:

@article{paulson2008comparison,
title={Comparison of dynamic susceptibility-weighted contrast-enhanced MR methods: recommendations for measuring relative cerebral blood volume in brain tumors},
author={Paulson, Eric S and Schmainda, Kathleen M},
journal={Radiology},
volume={249},
number={2},
pages={601--613},
year={2008},
publisher={Radiological Society of North America}
}

}

*/


class DiffusionDerivatives
{


private:
  inline void SetLongRunning(bool longRunning);

public:
  //! Default constructor
  DiffusionDerivatives()
  {
    //mLastErrorMessage = "";
  }

  //! Default destructor
  ~DiffusionDerivatives() {};

  std::vector<itk::Image<float, 3>::Pointer> Run(std::string dwiFile, std::string maskFile, std::string bvalFile, std::string gradFile, std::string outputDir);

  template < typename TInputPixelType, typename TMaskPixelType, typename TOutputTensorCompType >
  std::vector<itk::Image<float, 3>::Pointer> dtiRecon(std::string dwiFile, std::string gradFile, std::string bvalFile, std::string outputBasename, std::string maskFile, int verbose, int inputIsVectorImage);

  template <typename TMaskPixelType>
  std::vector<itk::Image<float, 3>::Pointer> runOverMaskType(std::string dwiFile, std::string gradFile, std::string bvalFile, std::string outputBasename, std::string maskFile, int verbose);

  template<typename ScalarImageType, typename TensorImageType>
  void allocateScalarIm(typename ScalarImageType::Pointer scIm, typename TensorImageType::Pointer tenIm);

  template<typename VectorImageType, typename TensorImageType>
  void allocateVectorIm(typename VectorImageType::Pointer scIm, typename TensorImageType::Pointer tenIm);
};


template<typename ScalarImageType, typename TensorImageType>
void DiffusionDerivatives::allocateScalarIm(typename ScalarImageType::Pointer scIm, typename TensorImageType::Pointer tenIm)
{
  scIm->SetOrigin(tenIm->GetOrigin());
  scIm->SetSpacing(tenIm->GetSpacing());
  scIm->SetDirection(tenIm->GetDirection());
  scIm->SetLargestPossibleRegion(tenIm->GetLargestPossibleRegion());
  scIm->SetRequestedRegion(tenIm->GetRequestedRegion());
  scIm->SetBufferedRegion(tenIm->GetBufferedRegion());
  scIm->Allocate();
}

template< typename VectorImageType, typename TensorImageType >
void DiffusionDerivatives::allocateVectorIm(typename VectorImageType::Pointer scIm, typename TensorImageType::Pointer tenIm)
{
  scIm->SetOrigin(tenIm->GetOrigin());
  scIm->SetSpacing(tenIm->GetSpacing());
  scIm->SetDirection(tenIm->GetDirection());
  scIm->SetLargestPossibleRegion(tenIm->GetLargestPossibleRegion());
  scIm->SetRequestedRegion(tenIm->GetRequestedRegion());
  scIm->SetBufferedRegion(tenIm->GetBufferedRegion());
  scIm->Allocate();
}

void GetImageInfo(std::string fName, itk::ImageIOBase::IOPixelType *pixelType, itk::ImageIOBase::IOComponentType *componentType)
{
  itk::ImageIOBase::Pointer imageIO;
  //~ try
  //~ {
  imageIO = itk::ImageIOFactory::CreateImageIO(fName.c_str(), itk::ImageIOFactory::ReadMode);
  if (imageIO)
  {
    imageIO->SetFileName(fName);
    imageIO->ReadImageInformation();
    *pixelType = imageIO->GetPixelType();
    *componentType = imageIO->GetComponentType();
  }
  else
  {
    std::cout << "Could not read the input image information from " <<
      fName << std::endl;
    //TODO should throw exception
    exit(EXIT_FAILURE);
  }

}

std::vector<itk::Image<float, 3>::Pointer> DiffusionDerivatives::Run(std::string dwiFile, std::string maskFile, std::string bvalFile, std::string gradFile, std::string outputDir)
{
  typedef itk::Image<float, 3> ScalarImageType;
  std::vector<ScalarImageType::Pointer> vectorOfDTIScalars;

  std::string prefix = "tensor";
  bool prefix_flag = 1;
  bool outputDir_flag = 1;

  std::string outputBasename = "";
  if (outputDir_flag)
    outputBasename += outputDir;

  if (prefix_flag)
    outputBasename += prefix;
  else
    outputBasename += "";

  itk::ImageIOBase::IOPixelType       maskPixelType;
  itk::ImageIOBase::IOComponentType   maskCompType;
  GetImageInfo(maskFile, &maskPixelType, &maskCompType);

  //Die if it isn't scalar.
  if (maskPixelType != itk::ImageIOBase::SCALAR)
  {
    std::cerr << "Mask is expected to be a Scalar image" << std::endl;
  }

  switch (maskCompType)
  {
  case itk::ImageIOBase::UCHAR:
  {
    vectorOfDTIScalars = runOverMaskType<unsigned char>(dwiFile, gradFile, bvalFile, outputBasename, maskFile, true);
    break;
  }
  case itk::ImageIOBase::SHORT:
  {
    vectorOfDTIScalars = runOverMaskType<short>(dwiFile, gradFile, bvalFile, outputBasename, maskFile, true);
    break;
  }
  default:
  {
    std::cerr << "Mask is expected to be either a unsigned char or a signed short Image" << std::endl;
    std::cerr << "Please invectigate the image type of your supplied mask" << std::endl;
    break;
  }
  }
  return vectorOfDTIScalars;
}
template <typename TMaskPixelType>

std::vector<itk::Image<float, 3>::Pointer>  DiffusionDerivatives::runOverMaskType(std::string dwiFile, std::string gradFile, std::string bvalFile, std::string outputBasename, std::string maskFile, int verbose)
{
  itk::ImageIOBase::IOPixelType       dwiPixelType;
  itk::ImageIOBase::IOComponentType   dwiCompType;
  GetImageInfo(dwiFile, &dwiPixelType, &dwiCompType);

  int inputIsVectorImage = 1;
  if (dwiPixelType == itk::ImageIOBase::SCALAR)
  {
    inputIsVectorImage = 0;
  }
  else
  {
    std::cerr << "Vector IMAGE\n";
  }

  switch (dwiCompType)
  {
  case itk::ImageIOBase::CHAR:
  {
    return dtiRecon<char, TMaskPixelType, float>(dwiFile, gradFile, bvalFile, outputBasename, maskFile, verbose, inputIsVectorImage);
    break;
  }
  case itk::ImageIOBase::UCHAR:
  {
    return dtiRecon<unsigned char, TMaskPixelType, float>(dwiFile, gradFile, bvalFile, outputBasename, maskFile, verbose, inputIsVectorImage);
    break;
  }
  case itk::ImageIOBase::USHORT:
  {
    return dtiRecon<unsigned short, TMaskPixelType, float>(dwiFile, gradFile, bvalFile, outputBasename, maskFile, verbose, inputIsVectorImage);
    break;
  }
  case itk::ImageIOBase::SHORT:
  {
    return dtiRecon<short, TMaskPixelType, float >(dwiFile, gradFile, bvalFile, outputBasename, maskFile, verbose, inputIsVectorImage);
    break;
  }
  case itk::ImageIOBase::UINT:
  {
    return dtiRecon<unsigned int, TMaskPixelType, float>(dwiFile, gradFile, bvalFile, outputBasename, maskFile, verbose, inputIsVectorImage);
    break;
  }
  case itk::ImageIOBase::INT:
  {
    return dtiRecon<int, TMaskPixelType, float>(dwiFile, gradFile, bvalFile, outputBasename, maskFile, verbose, inputIsVectorImage);
    break;
  }
  case itk::ImageIOBase::ULONG:
  {
    return dtiRecon<unsigned long, TMaskPixelType, float>(dwiFile, gradFile, bvalFile, outputBasename, maskFile, verbose, inputIsVectorImage);
    break;
  }
  case itk::ImageIOBase::LONG:
  {
    return dtiRecon<long, TMaskPixelType, float>(dwiFile, gradFile, bvalFile, outputBasename, maskFile, verbose, inputIsVectorImage);
    break;
  }
  case itk::ImageIOBase::FLOAT:
  {
    return dtiRecon<float, TMaskPixelType, float>(dwiFile, gradFile, bvalFile, outputBasename, maskFile, verbose, inputIsVectorImage);
    break;
  }
  case itk::ImageIOBase::DOUBLE:
  {
    return dtiRecon<double, TMaskPixelType, float>(dwiFile, gradFile, bvalFile, outputBasename, maskFile, verbose, inputIsVectorImage);
    break;
  }
  default:
  {
    std::cerr << "No match found.\n";
    std::vector<itk::Image<float, 3>::Pointer> emptyVec;
    return emptyVec;
    break;
  }
  }
}
template < typename TInputPixelType, typename TMaskPixelType, typename TOutputTensorCompType>

std::vector<itk::Image<float, 3>::Pointer>  DiffusionDerivatives::dtiRecon(std::string dwiFile, std::string gradFile, std::string bvalFile, std::string outputBasename, std::string maskFile, int verbose, int inputIsVectorImage)
{
  const unsigned int Dimension = 3;
  bool readb0 = false;
  double b0 = 0;
  typedef itk::VectorImage<TInputPixelType, 3> GradientImageType;
  typedef itk::Image<TInputPixelType, 4>   InputImageType;


  std::cout << "Converting DWI to vector image\n";

  ///TODO this should be a filter. cause we'll need to do it alot I think.
  typedef itk::ImageFileReader< InputImageType >  ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(dwiFile);
  typename InputImageType::Pointer img4D = reader->GetOutput();
  reader->Update();

  //Set up the gradient image size
  typename GradientImageType::Pointer gradIm = GradientImageType::New();


  typedef typename GradientImageType::RegionType   GradRegionType;
  typedef typename GradientImageType::SizeType     GradSizeType;
  typedef typename GradientImageType::IndexType    GradIndexType;
  typedef typename InputImageType::IndexType       Img4dIndexType;
  GradSizeType  sizeGradImage;
  typename InputImageType::SizeType size4D = img4D->GetLargestPossibleRegion().GetSize();
  sizeGradImage[0] = size4D[0];
  sizeGradImage[1] = size4D[1];
  sizeGradImage[2] = size4D[2];
  gradIm->SetVectorLength(size4D[3]);

  GradIndexType   indexGradImage = { { 0, 0, 0 } };
  GradRegionType  regionGradImage;
  regionGradImage.SetSize(sizeGradImage);
  regionGradImage.SetIndex(indexGradImage);
  gradIm->SetRegions(regionGradImage);
  typename InputImageType::SpacingType img4Dspacing = img4D->GetSpacing();
  typename InputImageType::PointType img4Dorigin = img4D->GetOrigin();
  typename InputImageType::DirectionType img4Ddir = img4D->GetDirection();

 typename GradientImageType::SpacingType gradSpacing;
  typename GradientImageType::PointType gradOrigin;
  typename GradientImageType::DirectionType gradDirs;

  gradSpacing[0] = img4Dspacing[0];  gradSpacing[1] = img4Dspacing[1];   gradSpacing[2] = img4Dspacing[2];
  gradOrigin[0] = img4Dorigin[0];   gradOrigin[1] = img4Dorigin[1];    gradOrigin[2] = img4Dorigin[2];

  for (unsigned int i = 0; i<3; ++i)
  {
    for (unsigned int j = 0; j<3; ++j)
    {
      gradDirs[i][j] = img4Ddir[i][j];
    }
  }

  gradIm->SetSpacing(gradSpacing);
  gradIm->SetOrigin(gradOrigin);
  gradIm->SetDirection(gradDirs);

  gradIm->Allocate();

  ///Copy data from img4d to gradim
  typedef itk::ImageRegionIteratorWithIndex< GradientImageType > IteratorType;
  typedef typename GradientImageType::PixelType GradPixelType;

  IteratorType it(gradIm, gradIm->GetRequestedRegion());

  ///Probably a better way to do this but I don't really know what it is.
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    GradIndexType   gradIndex = it.GetIndex();
    GradPixelType   gradPix = it.Get();
    Img4dIndexType  img4dIndex;
    img4dIndex[0] = gradIndex[0];
    img4dIndex[1] = gradIndex[1];
    img4dIndex[2] = gradIndex[2];

    for (unsigned int i = 0; i<size4D[3]; ++i)
    {
      img4dIndex[3] = i;
      gradPix.SetElement(i, img4D->GetPixel(img4dIndex));
    }
    it.Set(gradPix);
  }
  typedef itk::DiffusionTensor3DReconstructionImageFilter< TInputPixelType, TInputPixelType, TOutputTensorCompType > TensorReconstructionImageFilterType;
  // -------------------------------------------------------------------------

  int NumberOfGradients = 0;
  double bValue = 0;
  std::ifstream bvalIn(bvalFile.c_str());
  std::string line;
  while (!bvalIn.eof())
  {
    getline(bvalIn, line);
    double val;
    std::stringstream ss(line);

    while (ss >> val)
    {
      ++NumberOfGradients;
      if (val != 0 && bValue == 0)
        bValue = val;
      else if (val != 0 && bValue != val)
      {
        std::cerr << "multiple bvalues not allowed" << std::endl;
        //return EXIT_FAILURE;
      }
    }
  }
  std::ifstream bvecIn(gradFile.c_str());
  getline(bvecIn, line);
  std::stringstream Xss(line);
  getline(bvecIn, line);
  std::stringstream Yss(line);
  getline(bvecIn, line);
  std::stringstream Zss(line);

  typename TensorReconstructionImageFilterType::GradientDirectionType vect3d;
  typename TensorReconstructionImageFilterType::GradientDirectionContainerType::Pointer DiffusionVectors = TensorReconstructionImageFilterType::GradientDirectionContainerType::New();

  int counter = 0;
  double x, y, z;
  while (Xss >> x)
  {
    Yss >> y;
    Zss >> z;
    vect3d[0] = x; vect3d[1] = y; vect3d[2] = z;
    DiffusionVectors->InsertElement(counter, vect3d);
    ++counter;
  }
  // -------------------------------------------------------------------------

  typename TensorReconstructionImageFilterType::Pointer tensorReconstructionFilter = TensorReconstructionImageFilterType::New();

  typedef itk::Image<TMaskPixelType, 3 >   ImageMaskType;
  typedef itk::ImageMaskSpatialObject< 3 >   MaskType;
  MaskType::Pointer  spatialObjectMask = MaskType::New();
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
  tensorReconstructionFilter->SetMaskSpatialObject(spatialObjectMask);

  //---------------------------------------------------------------------------------
  tensorReconstructionFilter->SetGradientImage(DiffusionVectors, gradIm);
  tensorReconstructionFilter->SetNumberOfThreads(1);
  tensorReconstructionFilter->SetBValue(bValue);
  //CommandProgressUpdate::Pointer observer = CommandProgressUpdate::New();
  //tensorReconstructionFilter->AddObserver(itk::ProgressEvent(), observer);
  tensorReconstructionFilter->Update();

  //----------------------------------------------------------------

  //typedef itk::ImageFileWriter<TensorReconstructionImageFilterType::OutputImageType > TensorWriterType;
  //TensorWriterType::Pointer tensorWriter = TensorWriterType::New();
  //tensorWriter->SetFileName("C:/Projects/SampleData/DTI/AALZOutput/tensorFile.nii.gz");
  //tensorWriter->SetInput(tensorReconstructionFilter->GetOutput());
  //tensorWriter->Update();

  //--------------------------computing scalars----------------------------
  std::string default_ext = ".nii.gz";
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
  typedef itk::Image< float, 3 > ScalarImageType;
  typedef itk::VectorImage<float, 3 > VectorImageType;
  typedef itk::Image< TensorPixelType, 3 > TensorImageType;

  std::vector<ScalarImageType::Pointer> vectorOfDTIScalars;

  typedef itk::ImageRegionConstIteratorWithIndex < TensorImageType > ConstIterType;

  try
  {
    typename TensorImageType::Pointer tensorIm = TensorImageType::New();
    tensorIm = tensorReconstructionFilter->GetOutput();
    // Allocate each image that we will want to make..
    //FA
    typename ScalarImageType::Pointer faIm = ScalarImageType::New();

    //TR
    typename ScalarImageType::Pointer trIm = ScalarImageType::New();

    //Eigensys
    typename ScalarImageType::Pointer l1Im = ScalarImageType::New();
    typename ScalarImageType::Pointer l2Im = ScalarImageType::New();
    typename ScalarImageType::Pointer l3Im = ScalarImageType::New();
    typename VectorImageType::Pointer v1Im = VectorImageType::New();
    typename VectorImageType::Pointer v2Im = VectorImageType::New();
    typename VectorImageType::Pointer v3Im = VectorImageType::New();

    //Skewness & Kurtosis
    typename ScalarImageType::Pointer skIm = ScalarImageType::New();
    typename ScalarImageType::Pointer kuIm = ScalarImageType::New();

    //Geometric Features
    typename ScalarImageType::Pointer clIm = ScalarImageType::New();
    typename ScalarImageType::Pointer cpIm = ScalarImageType::New();
    typename ScalarImageType::Pointer csIm = ScalarImageType::New();

    //Radial Axial
    typename ScalarImageType::Pointer rdIm = ScalarImageType::New();
    typename ScalarImageType::Pointer adIm = ScalarImageType::New();

    //Gordons R features
    typename ScalarImageType::Pointer r1Im = ScalarImageType::New();
    typename ScalarImageType::Pointer r2Im = ScalarImageType::New();
    typename ScalarImageType::Pointer r3Im = ScalarImageType::New();

    //Gordons K features
    typename ScalarImageType::Pointer k1Im = ScalarImageType::New();
    typename ScalarImageType::Pointer k2Im = ScalarImageType::New();
    typename ScalarImageType::Pointer k3Im = ScalarImageType::New();

    //Allocate all the images...
    if (writeFA)
      allocateScalarIm<ScalarImageType, TensorImageType>(faIm, tensorIm);

    if (writeTR)
      allocateScalarIm<ScalarImageType, TensorImageType>(trIm, tensorIm);

    if (writeEign)
    {
      allocateScalarIm<ScalarImageType, TensorImageType>(l1Im, tensorIm);
      allocateScalarIm<ScalarImageType, TensorImageType>(l2Im, tensorIm);
      allocateScalarIm<ScalarImageType, TensorImageType>(l3Im, tensorIm);
      v1Im->SetVectorLength(3);
      v2Im->SetVectorLength(3);
      v3Im->SetVectorLength(3);
      allocateVectorIm<VectorImageType, TensorImageType>(v1Im, tensorIm);
      allocateVectorIm<VectorImageType, TensorImageType>(v2Im, tensorIm);
      allocateVectorIm<VectorImageType, TensorImageType>(v3Im, tensorIm);

    }

    if (writeSkew)
    {
      allocateScalarIm<ScalarImageType, TensorImageType>(skIm, tensorIm);
    }

    if (writeKurt)
    {
      allocateScalarIm<ScalarImageType, TensorImageType>(kuIm, tensorIm);
    }

    if (writeGeo)
    {
      allocateScalarIm<ScalarImageType, TensorImageType>(clIm, tensorIm);
      allocateScalarIm<ScalarImageType, TensorImageType>(cpIm, tensorIm);
      allocateScalarIm<ScalarImageType, TensorImageType>(csIm, tensorIm);
    }
    if (writeRadAx)
    {
      allocateScalarIm<ScalarImageType, TensorImageType>(rdIm, tensorIm);
      allocateScalarIm<ScalarImageType, TensorImageType>(adIm, tensorIm);
    }
    if (writeGordR)
    {
      allocateScalarIm<ScalarImageType, TensorImageType>(r1Im, tensorIm);
      allocateScalarIm<ScalarImageType, TensorImageType>(r2Im, tensorIm);
      allocateScalarIm<ScalarImageType, TensorImageType>(r3Im, tensorIm);
    }
    if (writeGordK)
    {
      allocateScalarIm<ScalarImageType, TensorImageType>(k1Im, tensorIm);
      allocateScalarIm<ScalarImageType, TensorImageType>(k2Im, tensorIm);
      allocateScalarIm<ScalarImageType, TensorImageType>(k3Im, tensorIm);
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
          m3 = (vcl_pow(l1 - m1, 3) + vcl_pow(l2 - m1, 3) + vcl_pow(l3 - m1, 3)) / (vcl_pow(l1, 3) + vcl_pow(l2, 3) + vcl_pow(l3, 3));
          if (m3 > 0)
          {
            skIm->SetPixel(index, vcl_pow(m3, static_cast<ScalarPixelType>(1.0 / 3.0)));
          }
          else
          {
            skIm->SetPixel(index, -1 * vcl_pow((-1 * m3), static_cast<ScalarPixelType>(1.0 / 3.0)));
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
          m4 = (vcl_pow(l1 - m1, 4) + vcl_pow(l2 - m1, 4) + vcl_pow(l3 - m1, 4)) / (vcl_pow(l1, 4) + vcl_pow(l2, 4) + vcl_pow(l3, 4));
          kuIm->SetPixel(index, vcl_pow(m4, static_cast<ScalarPixelType>(1.0 / 4.0)));
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
        m2 = (vcl_pow(lambda[0] - m1, 2) + vcl_pow(lambda[1] - m1, 2)
          + vcl_pow(lambda[2] - m1, 2)) / 3.0;
        m3 = (vcl_pow(lambda[0] - m1, 3) + vcl_pow(lambda[1] - m1, 3)
          + vcl_pow(lambda[2] - m1, 3)) / 3.0;

        r1 = sqrt(3 * (vcl_pow(m1, 2) + m2));
        r2 = sqrt(3 * m2 / 2 / (vcl_pow(m1, 2) + m2));
        //        r3 = sqrt(2) * m3 / vcl_pow(static_cast<double>(m2), static_cast<double>(1.5));
        r3 = (lambda[0] * lambda[1] * lambda[2]) / vcl_pow(static_cast<double>(sqrt(3 * m2)), 3);

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
        m2 = (vcl_pow(lambda[0] - m1, 2) + vcl_pow(lambda[1] - m1, 2)
          + vcl_pow(lambda[2] - m1, 2)) / 3.0;
        m3 = (vcl_pow(lambda[0] - m1, 3) + vcl_pow(lambda[1] - m1, 3)
          + vcl_pow(lambda[2] - m1, 3)) / 3.0;

        k1 = 3 * m1;
        k2 = sqrt(3 * m2);
        //        k3 = sqrt(2) * m3 / vcl_pow(static_cast<double>(m2), static_cast<double>(1.5));
        k3 = (lambda[0] * lambda[1] * lambda[2]) / vcl_pow(static_cast<double>(sqrt(3 * m2)), 3);


        k1Im->SetPixel(index, k1);
        k2Im->SetPixel(index, k2);
        k3Im->SetPixel(index, k3);
      }
    }

    std::cout << "Done Computing Scalars\n";

    typedef itk::ImageFileWriter< ScalarImageType >  ScalarWriterType;
    typedef itk::ImageFileWriter< VectorImageType >  VectorWriterType;

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
    cbica::Logging(loggerFile, "Something went wrong in Diffusion Derivatives L730: '" + std::string(excp.GetDescription()));
  }
  return vectorOfDTIScalars;
}

#endif