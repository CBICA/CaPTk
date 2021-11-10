/**
\file  DiffusionDerivatives.h

\brief The header file containing the DiffusionDerivatives class, used to calculate PSR, RCBV, and PH
Author: Saima Rathore
Library Dependecies: ITK 4.7+ <br>

https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/captk/license.html

*/


#ifndef _DiffusionDerivatives_h_
#define _DiffusionDerivatives_h_

#include "qstring.h"
#include "cbicaUtilities.h"
#include "FeatureReductionClass.h"
//#include "CaPTk.h"
#include "itkImageIOBase.h"
#include "itkExtractImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNiftiImageIO.h"
#include "itkDiffusionTensor3DReconstructionImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkNaryAddImageFilter.h"
#include "itkDivideImageFilter.h"
#include "cbicaLogging.h"
#include "CaPTkDefines.h"

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

#define Malloc(type, n) (type *)malloc((n)*sizeof(type))

#ifdef APP_BASE_CaPTk_H
#include "ApplicationBase.h"
#endif


class DiffusionDerivatives
#ifdef APP_BASE_CaPTk_H
	: public ApplicationBase
#endif
{
private:
  inline void SetLongRunning(bool longRunning);

public:
  std::string m_LastError;
	int max_line_len;
	char *line;

  //! Default constructor
  DiffusionDerivatives()
  {
	  max_line_len = 5000;
	  line = NULL;
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

  void FlipYOrientationInBVecFile(std::string inputfilename, std::string outputfilename);
  char* readline(FILE *input);
  void FlipAndShiftNiftiFile(std::string inputfilename, std::string outputfilename);
  void FlipAndShiftNiftiMask(std::string inputfilename, std::string outputfilename);

  void GetImageInfo(std::string fName, itk::ImageIOBase::IOPixelType *pixelType, itk::ImageIOBase::IOComponentType *componentType);
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

void DiffusionDerivatives::GetImageInfo(std::string fName, itk::ImageIOBase::IOPixelType *pixelType, itk::ImageIOBase::IOComponentType *componentType)
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
    std::cout << "Could not read the input image information from " <<fName << std::endl;
	m_LastError = "Could not read the input image information from " + fName;
  }
}
char* DiffusionDerivatives::readline(FILE *input)
{
	int len;

	if (fgets(line, max_line_len, input) == NULL)
		return NULL;

	while (strrchr(line, '\n') == NULL)
	{
		max_line_len *= 2;
		line = (char *)realloc(line, max_line_len);
		len = (int)strlen(line);
		if (fgets(line + len, max_line_len - len, input) == NULL)
			break;
	}
	return line;
}

void DiffusionDerivatives::FlipYOrientationInBVecFile(std::string inputfilename, std::string outputfilename)
{
	std::vector<long double> XOrientation;
	std::vector<long double> YOrientation;
	std::vector<long double> ZOrientation;
	//std::string::size_type sz;

	int LineNumber = 1;
	max_line_len = 5000;
	FILE *fp = fopen(inputfilename.c_str(), "r");
	if (fp == NULL)
		exit(1);
	line = Malloc(char, max_line_len);
	while (readline(fp) != NULL)
	{
		char *p = strtok(line, "\t"); // label
		char * first = strtok(line, " ");
		std::string firstElement(first);
		//char* pEnd;
		if (LineNumber == 1)
			XOrientation.push_back(std::stold(firstElement));
		else if (LineNumber == 2)
			YOrientation.push_back(std::stold(firstElement));
		else if (LineNumber == 3)
			ZOrientation.push_back(std::stold(firstElement));
		while (1)
		{
			p = strtok(NULL, " ");
			if (p == NULL || *p == '\n') // check '\n' as ' ' may be after the last feature
				break;
			std::string element(p);
			if (LineNumber == 1)
				XOrientation.push_back(std::stold(element));
			else if (LineNumber == 2)
				YOrientation.push_back(std::stold(element));
			else if (LineNumber == 3)
				ZOrientation.push_back(std::stold(element));
		}
		LineNumber++;
	}
	fclose(fp);
	for (int i = 1; i < YOrientation.size(); i++)
		YOrientation[i] = -1 * YOrientation[i];

	fp = fopen(outputfilename.c_str(), "w");
	for (int index = 0; index < XOrientation.size(); index++)
	{
		if (index<XOrientation.size() - 1)
			fprintf(fp, "%.15Lf ", XOrientation[index]);
		else
			fprintf(fp, "%.15Lf\n", XOrientation[index]);
	}
	for (int index = 0; index < YOrientation.size(); index++)
	{
		if (index<YOrientation.size() - 1)
			fprintf(fp, "%.15Lf ", YOrientation[index]);
		else
			fprintf(fp, "%.15Lf\n", YOrientation[index]);
	}
	for (int index = 0; index < ZOrientation.size(); index++)
	{
		if (index<ZOrientation.size() - 1)
			fprintf(fp, "%.15Lf ", ZOrientation[index]);
		else
			fprintf(fp, "%.15Lf", ZOrientation[index]);
	}
	fclose(fp);
}

void DiffusionDerivatives::FlipAndShiftNiftiMask(std::string inputfilename, std::string outputfilename)
{
 /* typedef short                      InputPixelType;
  typedef short                       MiddlePixelType;
  typedef short                       OutputPixelType;
  typedef itk::Image< InputPixelType, 3 >    InputImageType;
  typedef itk::Image< MiddlePixelType, 3 >    MiddleImageType;
  typedef itk::Image< OutputPixelType, 3 >    OutputImageType;

  typedef itk::FlipImageFilter< InputImageType> FlipType;
  FlipType::Pointer flip = FlipType::New();

  typedef itk::ImageFileReader< InputImageType  >  ReaderType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  reader->SetFileName(inputfilename);
  writer->SetFileName(outputfilename);

  FlipType::FlipAxesArrayType flipAxesSet;
  flipAxesSet[0] = 0;
  flipAxesSet[1] = -1;
  flipAxesSet[2] = 0;
  flip->SetFlipAxes(flipAxesSet);
  flip->FlipAboutOriginOff();
  flip->SetInput(reader->GetOutput());
  flip->Update();

  itk::NiftiImageIO::Pointer nifti_io = itk::NiftiImageIO::New();
  OutputImageType::Pointer image = OutputImageType::New();
  image = flip->GetOutput();
  image->SetOrigin(reader->GetOutput()->GetOrigin());
  writer->SetInput(image);
  writer->SetImageIO(nifti_io);
  try
  {
    writer->Update();
  }
  catch (itk::ExceptionObject & err)
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    m_LastError = err.GetDescription();
  }*/
  //itk::AnalyzeImageIO::Pointer io = itk::AnalyzeImageIO::New();
  //ImageReaderType::Pointer fileReader = ImageReaderType::New();
  //fileReader->SetImageIO(io);
  //fileReader->SetFileName(path);
  //fileReader->Update();
  //ImageType::Pointer rval = fileReader->GetOutput();
  //itk::OrientImageFilter<ImageType, ImageType>::Pointer orienter =
  //  itk::OrientImageFilter<ImageType, ImageType>::New();
  //orienter->UseImageDirectionOn();
  //orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIP);
  //orienter->SetInput(rval);
  //orienter->Update();
  //rval = orienter->GetOutput();
  //return rval;
}


void DiffusionDerivatives::FlipAndShiftNiftiFile(std::string inputfilename, std::string outputfilename)
{

	typedef short                      InputPixelType;
	typedef short                       MiddlePixelType;
	typedef short                       OutputPixelType;
	typedef itk::Image< InputPixelType, 4 >    InputImageType;
	typedef itk::Image< MiddlePixelType, 4 >    MiddleImageType;
	typedef itk::Image< OutputPixelType, 4 >    OutputImageType;

	typedef itk::FlipImageFilter< InputImageType> FlipType;
	FlipType::Pointer flip = FlipType::New();

	typedef itk::ImageFileReader< InputImageType  >  ReaderType;
	typedef itk::ImageFileWriter< OutputImageType >  WriterType;

	ReaderType::Pointer reader = ReaderType::New();
	WriterType::Pointer writer = WriterType::New();

	reader->SetFileName(inputfilename);
	writer->SetFileName(outputfilename);


	FlipType::FlipAxesArrayType flipAxesSet;

	flipAxesSet[0] = 0;
	flipAxesSet[1] = -1;
	flipAxesSet[2] = 0;
	flipAxesSet[3] = 0;

	flip->SetFlipAxes(flipAxesSet);
	flip->FlipAboutOriginOff();
	flip->SetInput(reader->GetOutput());
	flip->Update();



	itk::NiftiImageIO::Pointer nifti_io = itk::NiftiImageIO::New();
	OutputImageType::Pointer image = OutputImageType::New();
	image = flip->GetOutput();
	image->SetOrigin(reader->GetOutput()->GetOrigin());
	writer->SetInput(image);
	writer->SetImageIO(nifti_io);
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		m_LastError = err.GetDescription();
	}

}

std::vector<itk::Image<float, 3>::Pointer> DiffusionDerivatives::Run(std::string dwiFile, std::string maskFile, std::string bvalFile, std::string gradFile, std::string outputDir)
{
	std::string dwiFileNmae = cbica::getFilenameBase(dwiFile) + cbica::getFilenameExtension(dwiFile);
	std::string gradFileNmae = cbica::getFilenameBase(gradFile) + cbica::getFilenameExtension(gradFile);
  std::string maskFileName = cbica::getFilenameBase(maskFile) + cbica::getFilenameExtension(maskFile);
	//this->FlipYOrientationInBVecFile(gradFile, outputDir + "/"+gradFileNmae);
	//this->FlipAndShiftNiftiFile(dwiFile, outputDir + "/lps_" + dwiFileNmae);
 // this->FlipAndShiftNiftiMask(maskFile, outputDir + "/lps_" + maskFileName);



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
	m_LastError = "Mask is expected to be a Scalar image";
  }

  switch (maskCompType)
  {
  case itk::ImageIOBase::UCHAR:
  {
	  vectorOfDTIScalars = runOverMaskType<unsigned char>(outputDir + "/lps_" + dwiFileNmae, outputDir + "/" + gradFileNmae, bvalFile, outputBasename, outputDir + "/lps_" + maskFileName, true);
    break;
  }
  case itk::ImageIOBase::SHORT:
  {
	  //vectorOfDTIScalars = runOverMaskType<short>(outputDir + "/" + dwiFileNmae, outputDir + "/" + gradFileNmae, bvalFile, outputBasename, outputDir + "/" + maskFileName, true);
    vectorOfDTIScalars = runOverMaskType<short>(dwiFile, gradFile, bvalFile, outputBasename, maskFile, true);
    break;
  }
  default:
  {
	  m_LastError = "Mask is expected to be either a unsigned char or a signed short Image. Please investigate the image type of your supplied mask.";
      std::cerr << "Mask is expected to be either a unsigned char or a signed short Image. Please investigate the image type of your supplied mask." << std::endl;
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
  std::vector<int> indexesWhereBValZero;
  int currentTimeIndex = 0;
  while (!bvalIn.eof())
  {
    getline(bvalIn, line);
    double val;
    std::stringstream ss(line);

    while (ss >> val)
    {
      if (val == 0) //  mark these indices for b0 extraction
         indexesWhereBValZero.push_back(currentTimeIndex);

      ++NumberOfGradients;
      if (val != 0 && bValue == 0)
        bValue = val;
      else if (val != 0 && bValue != val)
      {
        std::cerr << "multiple unique bvalues are not currently allowed" << std::endl;
		m_LastError = "multiple unique bvalues are not currently allowed";
      }
      currentTimeIndex++;
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
	m_LastError = err.GetDescription();
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
  static int writeBZero = 1; 
  // The way this is currently done means the calling code can break depending on if these features ever get enabled.
  // TODO: Fix this up

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

    // b0 image
    typename ScalarImageType::Pointer bZeroIm = ScalarImageType::New();

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
    if (writeBZero)
    {
      allocateScalarIm<ScalarImageType, TensorImageType>(bZeroIm, tensorIm);
      typename InputImageType::RegionType inputRegion = img4D->GetLargestPossibleRegion();
      InputImageType::SizeType desiredSize = inputRegion.GetSize();
      desiredSize[3] = 0; // We'll always collapse along the 4th dimension

      // Add all b0 together first to average the images
      typedef itk::NaryAddImageFilter<ScalarImageType, ScalarImageType> AddImageFilterType; // NaryAdd makes sure we can loop easily
      typename AddImageFilterType::Pointer nAdder = AddImageFilterType::New();
      for (int i = 0; i < indexesWhereBValZero.size(); i++)
      {
          typedef itk::ExtractImageFilter<InputImageType, ScalarImageType> ExtractorType;
          typename ExtractorType::Pointer extractor = ExtractorType::New();
          typename InputImageType::IndexType index = inputRegion.GetIndex();
          index[3] = indexesWhereBValZero[i]; // Extract only the current time index

          InputImageType::RegionType targetRegion;
          targetRegion.SetSize(desiredSize);
          targetRegion.SetIndex(index);

          extractor->SetInput(img4D);
          extractor->SetExtractionRegion(targetRegion);
          extractor->SetDirectionCollapseToSubmatrix();
          extractor->Update();

          nAdder->SetInput(i, extractor->GetOutput());
      }
      nAdder->Update();
      // Now divide the added image by the number of images N
      typedef itk::DivideImageFilter<ScalarImageType, ScalarImageType, ScalarImageType> DividerType;
      typename DividerType::Pointer divider = DividerType::New();
      divider->SetInput(nAdder->GetOutput());
      divider->SetConstant((double)indexesWhereBValZero.size());
      divider->Update();

      // Apply mask to the now-averaged b0
      typedef itk::MaskImageFilter<ScalarImageType, ImageMaskType, ScalarImageType> MaskerType;
      typename MaskerType::Pointer masker = MaskerType::New();
      masker->SetInput(divider->GetOutput());
      masker->SetMaskImage(maskReader->GetOutput());
      masker->Update();
      bZeroIm = masker->GetOutput();
      
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
        l1 = std::abs(lambda[0]);
        l2 = std::abs(lambda[1]);
        l3 = std::abs(lambda[2]);
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
        l1 = std::abs(lambda[0]);
        l2 = std::abs(lambda[1]);
        l3 = std::abs(lambda[2]);
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
	else
		vectorOfDTIScalars.push_back(NULL);

    if (writeTR)
      vectorOfDTIScalars.push_back(trIm);
	else
		vectorOfDTIScalars.push_back(NULL);

    if (writeRadAx)
    {
      vectorOfDTIScalars.push_back(rdIm);
      vectorOfDTIScalars.push_back(adIm);
    }
	else
	{
	vectorOfDTIScalars.push_back(NULL);
	vectorOfDTIScalars.push_back(NULL);
	}

    if (writeBZero)
    {
        vectorOfDTIScalars.push_back(bZeroIm);
    }
    else
    {
        vectorOfDTIScalars.push_back(NULL);
    }
  }
  catch (itk::ExceptionObject & excp)
  {
    cbica::Logging(loggerFile, "Something went wrong in Diffusion Derivatives L730: '" + std::string(excp.GetDescription()));
	m_LastError = "Error in computing DTI scalars.";
  }
  return vectorOfDTIScalars;
}

#endif