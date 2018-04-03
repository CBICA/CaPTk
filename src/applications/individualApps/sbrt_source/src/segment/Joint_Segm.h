#pragma once
#include <iostream>
#include <iterator>
#include <stdlib.h>
#include <string>
#include <algorithm>
#include <exception>
#include <typeinfo>
#include <stdexcept>
#include <stdio.h>
#include <numeric>
#include <functional>
#include <vector> 

#include "cbicaCmdParser.h"
#include "cbicaITKSafeImageIO.h"

#include "cbicaUtilities.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMedianImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryMask3DMeshSource.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkBinaryImageToLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelSelectionLabelMapFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkLabelToRGBImageFilter.h"
#include "itkSobelEdgeDetectionImageFilter.h"
#include "itkBinaryContourImageFilter.h"
#include "itkAddImageFilter.h"


using TPixel = float;
const unsigned int  Dimension = 3;
typedef itk::RGBPixel<unsigned char>         RGBPixelType;
class Joint_Segm
{
public:
  Joint_Segm();
  ~Joint_Segm();



  //template<typename TPixel, unsigned int VImageDimension  >
  template< class TImageType >
  typename TImageType::Pointer generateseeds(typename TImageType::Pointer CT_image, const typename TImageType::Pointer PET_imgin)
  {
    typename TImageType::Pointer Pet_bound = TImageType::New();
    Pet_bound->SetRegions(PET_imgin->GetLargestPossibleRegion());
    Pet_bound->SetSpacing(PET_imgin->GetSpacing());
    Pet_bound->SetOrigin(PET_imgin->GetOrigin());
    Pet_bound->SetDirection(PET_imgin->GetDirection());
    Pet_bound->Allocate();
    Pet_bound->FillBuffer(0);
    typename TImageType::Pointer Pet_bound1 = TImageType::New();
    Pet_bound1->SetRegions(PET_imgin->GetLargestPossibleRegion());
    Pet_bound1->SetSpacing(PET_imgin->GetSpacing());
    Pet_bound1->SetOrigin(PET_imgin->GetOrigin());
    Pet_bound1->SetDirection(PET_imgin->GetDirection());
    Pet_bound1->Allocate();
    Pet_bound1->FillBuffer(0);
    typename TImageType::IndexType index;
    index[0] = 100;
    index[1] = 100;
    index[2] = 50;
    //  auto test0 = PET_imgin->GetPixel(index);
    //  auto test1 = Pet_bound->GetPixel(index);
    //   auto test2 = Pet_bound1->GetPixel(index);



    //Smooth image 
    typedef TImageType ImageType;
    typedef itk::MedianImageFilter<ImageType, ImageType > FilterType;
    typename  FilterType::Pointer medianFilter = FilterType::New();
    typedef  itk::ImageFileWriter< ImageType  > WriterType;

    typename  WriterType::Pointer writer = WriterType::New();
    typename  FilterType::InputSizeType radius1;
    typedef itk::MaskImageFilter< ImageType, ImageType > MaskFilterType;
    typename  MaskFilterType::Pointer maskFilter = MaskFilterType::New();


    radius1.Fill(2);
    medianFilter->SetRadius(radius1);
    medianFilter->SetInput(CT_image);
    medianFilter->Update();

    typedef itk::SmoothingRecursiveGaussianImageFilter<
      ImageType, ImageType >  Filter;
    typename  Filter::Pointer smoothingfilter = Filter::New();
    smoothingfilter->SetInput(medianFilter->GetOutput());
    const double sigma = 0.5;
    smoothingfilter->SetSigma(sigma);
    smoothingfilter->Update();

    typedef itk::BinaryThresholdImageFilter < ImageType, ImageType >
      BinaryThresholdImageFilterType;

    typename  BinaryThresholdImageFilterType::Pointer binarythresholdFilter
      = BinaryThresholdImageFilterType::New();
    binarythresholdFilter->SetInput(smoothingfilter->GetOutput());
    binarythresholdFilter->SetLowerThreshold(-8000);
    binarythresholdFilter->SetUpperThreshold(-3);
    binarythresholdFilter->SetInsideValue(255);
    binarythresholdFilter->SetOutsideValue(0);
    binarythresholdFilter->Update();

    typedef itk::BinaryBallStructuringElement <typename ImageType::PixelType, Dimension >
      StructuringElementType;
    StructuringElementType structuringElement;
    radius1.Fill(3);
    structuringElement.SetRadius(radius1);
    structuringElement.CreateStructuringElement();
    typedef itk::BinaryMorphologicalOpeningImageFilter < ImageType, ImageType, StructuringElementType >
      BinaryMorphologicalOpeningImageFilterType;
    typename   BinaryMorphologicalOpeningImageFilterType::Pointer openingFilter
      = BinaryMorphologicalOpeningImageFilterType::New();
    openingFilter->SetInput(binarythresholdFilter->GetOutput());
    openingFilter->SetKernel(structuringElement);
    openingFilter->Update();
    for (int i = 0; i < 4; i++)
    {
      openingFilter->SetInput(openingFilter->GetOutput());
      openingFilter->SetKernel(structuringElement);
      openingFilter->Update();
    }
    typedef itk::BinaryMorphologicalClosingImageFilter < ImageType, ImageType, StructuringElementType >
      BinaryMorphologicalClosingImageFilterType;
    typename BinaryMorphologicalClosingImageFilterType::Pointer closingFilter
      = BinaryMorphologicalClosingImageFilterType::New();
    closingFilter->SetInput(openingFilter->GetOutput());

    closingFilter->SetKernel(structuringElement);
    closingFilter->Update();
    for (int i = 0; i < 4; i++)
    {
      closingFilter->SetInput(closingFilter->GetOutput());
      closingFilter->SetKernel(structuringElement);
      closingFilter->Update();

    }

    maskFilter->SetInput(smoothingfilter->GetOutput());
    maskFilter->SetMaskImage(closingFilter->GetOutput());
    maskFilter->Update();

    //*Generate mask from PET image *//
    radius1.Fill(2);
    medianFilter->SetRadius(radius1);
    medianFilter->SetInput(PET_imgin);
    medianFilter->Update();
    maskFilter->SetInput(medianFilter->GetOutput());
    maskFilter->SetMaskImage(closingFilter->GetOutput());
    maskFilter->Update();
    //cbica::WriteImage<ImageType>(maskFilter->GetOutput(),"C:/Users/sridharp/Documents/SBRT-lung/trunk/data/training_data/ANON17670/CT_mask.nii");

    typedef itk::MinimumMaximumImageCalculator < ImageType >
      ImageCalculatorFilterType;
    typename ImageCalculatorFilterType::Pointer imageCalculatorFilter
      = ImageCalculatorFilterType::New();
    imageCalculatorFilter->SetImage(medianFilter->GetOutput());
    imageCalculatorFilter->Compute();
    int max_int = imageCalculatorFilter->GetMaximum();
    float mean_int = max_int / 2.7;

    typedef itk::ThresholdImageFilter <ImageType>  ThresholdImageFilterType;
    typename ThresholdImageFilterType::Pointer thresholdFilter = ThresholdImageFilterType::New();
    thresholdFilter->SetInput(maskFilter->GetOutput());
    thresholdFilter->ThresholdOutside(mean_int, max_int);

    //  thresholdFilter->SetOutsideValue(0);
    thresholdFilter->Update();

    //   InputImageFileReaderPointerType PET_maskreader = ImageReaderType::New();
    //  PET_maskreader->SetFileName("C:\\Users\\sridharp\\Documents\\SBRT-lung\\trunk\\data\\test_data\\iur_PET.nii");
    //  PET_maskreader->Update();

    typedef itk::Image<  unsigned int, Dimension > OutputImageType;
    //  typedef itk::LabelObject< PixelType, Dimension >  LabelObjectType;
    //	  typedef itk::LabelMap< LabelObjectType >          LabelMapType;
    typedef itk::ConnectedComponentImageFilter < ImageType, OutputImageType > ConnectedComponentImageFilterType;

    typename ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();

    connected->SetInput(thresholdFilter->GetOutput());
    connected->Update();
    std::cout << "regions" << connected->GetObjectCount();
    /* typedef itk::Image<RGBPixelType, Dimension>  RGBImageType;
    typedef itk::LabelToRGBImageFilter<OutputImageType, RGBImageType> RGBFilterType;
    RGBFilterType::Pointer rgbFilter = RGBFilterType::New();
    rgbFilter->SetInput(connected->GetOutput());*/

    //typedef  itk::ImageFileWriter< RGBImageType  > WriterType1;
    //WriterType1::Pointer writer1 = WriterType1::New();
    //writer1->SetFileName("C:\\Users\\sridharp\\Documents\\SBRT-lung\\trunk\\data\\test_data\\rgb.nii");
    //writer1->SetInput(rgbFilter->GetOutput());
    //writer1->Update();
    typedef itk::BinaryContourImageFilter <OutputImageType, ImageType >
      binaryContourImageFilterType;

    typename binaryContourImageFilterType::Pointer binaryContourFilter
      = binaryContourImageFilterType::New();
    binaryContourFilter->SetInput(connected->GetOutput());

    typedef itk::LabelImageToLabelMapFilter<OutputImageType> LabelImageToLabelMapFilterType;
    typename  LabelImageToLabelMapFilterType::Pointer labelImageToLabelMapFilter = LabelImageToLabelMapFilterType::New();

    labelImageToLabelMapFilter->SetInput(connected->GetOutput());
    labelImageToLabelMapFilter->Update();
    //std::cout << "regions" << labelImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects();

    typedef itk::LabelMapToLabelImageFilter< LabelImageToLabelMapFilterType::OutputImageType, OutputImageType> LabelMapToLabelImageFilterType;
    typename  LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
    labelMapToLabelImageFilter->SetInput(labelImageToLabelMapFilter->GetOutput());
    labelMapToLabelImageFilter->Update();

    typedef itk::SobelEdgeDetectionImageFilter < ImageType, ImageType >
      SobelEdgeDetectionImageFilterType;
    typename SobelEdgeDetectionImageFilterType::Pointer sobelFilter
      = SobelEdgeDetectionImageFilterType::New();
    sobelFilter->SetInput(binaryContourFilter->GetOutput());
    sobelFilter->Update();

    typename ThresholdImageFilterType::Pointer thresholdFilter1 = ThresholdImageFilterType::New();
    thresholdFilter1->SetInput(sobelFilter->GetOutput());
    thresholdFilter1->ThresholdAbove(0);
    thresholdFilter1->SetOutsideValue(1);

    // cbica::WriteImage<ImageType>(thresholdFilter1->GetOutput(), "C:\\Users\\sridharp\\Documents\\SBRT-lung\\trunk\\data\\test_data\\bound.nii");

    typedef itk::LabelStatisticsImageFilter< ImageType, OutputImageType > LabelStatisticsImageFilterType;
    typename LabelStatisticsImageFilterType::Pointer labelStatisticsImageFilter = LabelStatisticsImageFilterType::New();
    labelStatisticsImageFilter->SetLabelInput(labelMapToLabelImageFilter->GetOutput());
    labelStatisticsImageFilter->SetInput(thresholdFilter->GetOutput());
    labelStatisticsImageFilter->Update();

    typedef typename LabelStatisticsImageFilterType::ValidLabelValuesContainerType ValidLabelValuesType;
    typedef typename LabelStatisticsImageFilterType::LabelPixelType  LabelPixelType;
    std::map<float, int>value;
    std::multimap<float, int>  region_value;
    std::vector<typename ImageType::IndexType>   regionindex;
    std::vector<typename ImageType::IndexType>     foregroundseeds;
    std::vector<typename ImageType::IndexType>     backgroundseeds;
    int mfseeds = 200;


    typedef itk::AddImageFilter <ImageType, ImageType > AddImageFilterType;
    typename AddImageFilterType::Pointer addFilter = AddImageFilterType::New();



    for (typename ValidLabelValuesType::const_iterator vIt = labelStatisticsImageFilter->GetValidLabelValues().begin();
      vIt != labelStatisticsImageFilter->GetValidLabelValues().end();
      ++vIt)
    {
      if (labelStatisticsImageFilter->HasLabel(*vIt) && *vIt > 0 && labelStatisticsImageFilter->GetCount(*vIt) > 2)
      {
        LabelPixelType labelValue = *vIt;

        std::cout << "min: " << labelStatisticsImageFilter->GetMinimum(labelValue) << std::endl;
        std::cout << "max: " << labelStatisticsImageFilter->GetMaximum(labelValue) << std::endl;
        std::cout << "count: " << labelStatisticsImageFilter->GetCount(labelValue) << std::endl;
        typename TImageType::IndexType start = labelStatisticsImageFilter->GetRegion(labelValue).GetIndex();
        std::cout << "index: " << labelStatisticsImageFilter->GetRegion(labelValue).GetUpperIndex();
        typename TImageType::SizeType size = labelStatisticsImageFilter->GetRegion(labelValue).GetSize();
        typename TImageType::RegionType desiredRegion;
        desiredRegion.SetSize(size);
        desiredRegion.SetIndex(start);
        int counter = 0;
        itk::ImageRegionIteratorWithIndex< ImageType> imageIterator(thresholdFilter->GetOutput(), desiredRegion);
        while (!imageIterator.IsAtEnd())
        {

          region_value.insert(std::pair<float, int>(imageIterator.Get(), counter));
          regionindex.push_back(imageIterator.GetIndex());
          counter++;
          ++imageIterator;
        }
        std::map<float, int>::reverse_iterator itr = region_value.rbegin();
        for (int i = 0; i < mfseeds; i++)
        {
          if (itr != region_value.rend())
          {
            foregroundseeds.push_back(regionindex.at(itr->second));
            itr++;
          }
        }
        region_value.clear();
        regionindex.clear();

      }
    }

    for (unsigned int i = 0; i < foregroundseeds.size(); i++)
    {
      Pet_bound->SetPixel(foregroundseeds.at(i), 2);
    }

    //cbica::WriteImage<ImageType>(Pet_bound, "C:\\Users\\sridharp\\Documents\\SBRT-lung\\trunk\\data\\test_data\\foreground_seeds.nii");

    //pasteFilter->SetSourceImage(thresholdFilter->GetOutput());
    //pasteFilter->SetDestinationImage(newimg);
    //pasteFilter->SetSourceRegion(labelStatisticsImageFilter->GetRegion(labelValue));
    //pasteFilter->SetDestinationIndex(start);
    //pasteFilter->Update();

    //thresholdFilter->SetInput(pasteFilter->GetOutput());
    //thresholdFilter->ThresholdOutside(mean_int, imageCalculatorFilter->GetMaximum());
    //thresholdFilter->SetOutsideValue(0);
    //thresholdFilter->Update();


    //binaryContourFilter->SetInput(roi->GetOutput());
    //binaryContourFilter->SetForegroundValue(1);
    //binaryContourFilter->SetBackgroundValue(0);
    //binaryContourFilter->Update();
    //Pet_bound = binaryContourFilter->GetOutput();
    //Pet_bound->SetOrigin(Pet_bound1->GetOrigin());
    //Pet_bound->SetDirection(Pet_bound1->GetDirection());
    //Pet_bound->SetSpacing(Pet_bound1->GetSpacing());

    // foreground_seeds.push_back(imageCalculatorFilter->GetIndexOfMaximum());
    //   Pet_bound->SetPixel(imageCalculatorFilter->GetIndexOfMaximum(), 1);

    //  itk::ImageRegionConstIterator<ImageType> imageIterator(PET_mask->GetOutput(),
    //    labelStatisticsImageFilter->GetRegion(labelValue));
    //  while (!imageIterator.IsAtEnd())
    //  {
    //    // Get the value of the current pixel
    //    int val = imageIterator.Get();
    //       //std::cout << imageIterator.Get();
    //    if (val <= mean_int && val > 0)
    //    {
    //      backgroundseeds.push_back(imageIterator.GetIndex());
    //      Pet_bound->SetPixel(imageIterator.GetIndex(), 2);
    //    }
    //    ++imageIterator;
    //  }

    addFilter->SetInput1(thresholdFilter1->GetOutput());

    addFilter->SetInput2(Pet_bound);
    addFilter->Update();
    Pet_bound = addFilter->GetOutput();
    itk::ImageRegionIterator<ImageType> it(Pet_bound, Pet_bound->GetRequestedRegion());
    it.GoToBegin();
    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
      if (it.Get() > 1)
      {
        it.Set(2);
      }
      else
      {
        it.Set(0);
      }
    }

    return Pet_bound;
    //find foreground and background seeds
  }

};
