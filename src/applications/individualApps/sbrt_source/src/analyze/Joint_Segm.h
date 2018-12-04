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

#include "cbicaUtilities.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMedianImageFilter.h"
#include <itkMaskImageFilter.h>
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


class Joint_Segm
{
  public:
    Joint_Segm();
    ~Joint_Segm();



    template<typename TPixel, unsigned int VImageDimension  >
    void generateseeds(itk::Image<TPixel, VImageDimension> *image, itk::Image<TPixel, VImageDimension> *PET_image, itk::Image<TPixel, VImageDimension> *mask);


};


template<typename TPixel, unsigned int VImageDimension >
void Joint_Segm::generateseeds(itk::Image< TPixel, VImageDimension> *CT_image, itk::Image< TPixel, VImageDimension> *PET_image,
  itk::Image< TPixel, VImageDimension> *mask)
{
  //Smooth image 
  typedef itk::Image< TPixel, VImageDimension > ImageType;
  typedef itk::MedianImageFilter<ImageType, ImageType > FilterType;
  FilterType::Pointer medianFilter = FilterType::New();
  WriterType::Pointer writer = WriterType::New();
  FilterType::InputSizeType radius1;
  typedef itk::MaskImageFilter< ImageType, ImageType > MaskFilterType;
  MaskFilterType::Pointer maskFilter = MaskFilterType::New();


  			  radius1.Fill(2);
  			  medianFilter->SetRadius(radius1);
  			  medianFilter->SetInput(CT_image);
  			  medianFilter->Update();
  
  			  typedef itk::SmoothingRecursiveGaussianImageFilter<
  				  ImageType, ImageType >  Filter;
  			  Filter::Pointer smoothingfilter = Filter::New();
  			  smoothingfilter->SetInput(medianFilter->GetOutput());
  			 const double sigma = 0.5;
  			 smoothingfilter->SetSigma(sigma);
  			 smoothingfilter->Update();

         typedef itk::BinaryThresholdImageFilter <ImageType, ImageType>
           				  BinaryThresholdImageFilterType;
           
           			  BinaryThresholdImageFilterType::Pointer binarythresholdFilter
           				  = BinaryThresholdImageFilterType::New();
                   binarythresholdFilter->SetInput(smoothingfilter->GetOutput());
                   binarythresholdFilter->SetLowerThreshold(100);
                   binarythresholdFilter->SetUpperThreshold(500);
                   binarythresholdFilter->SetInsideValue(255);
                  binarythresholdFilter->SetOutsideValue(0);
                   binarythresholdFilter->Update();

                   	typedef itk::BinaryBallStructuringElement<ImageType::PixelType, ImageType::ImageDimension>
                   		StructuringElementType;
                   	StructuringElementType structuringElement;
                   	radius1.Fill(3);
                   	structuringElement.SetRadius(radius1);
                   	structuringElement.CreateStructuringElement();
                   	typedef itk::BinaryMorphologicalOpeningImageFilter <ImageType, ImageType, StructuringElementType>
                   		BinaryMorphologicalOpeningImageFilterType;
                   	BinaryMorphologicalOpeningImageFilterType::Pointer openingFilter
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
                   	typedef itk::BinaryMorphologicalClosingImageFilter <ImageType, ImageType, StructuringElementType>
                   		BinaryMorphologicalClosingImageFilterType;
                   	BinaryMorphologicalClosingImageFilterType::Pointer closingFilter
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

                typedef itk::MinimumMaximumImageCalculator <ImageType>
                  ImageCalculatorFilterType;
                ImageCalculatorFilterType::Pointer imageCalculatorFilter
                  = ImageCalculatorFilterType::New();
                imageCalculatorFilter->SetImage(PET_imgin);
                imageCalculatorFilter->Compute();
                int max_int = imageCalculatorFilter->GetMaximum();
                float mean_int = max_int / 4;
                InputImageFileReaderPointerType PET_mask = ImageReaderType::New();
                PET_mask->SetFileName("C:\\Users\\sridharp\\Documents\\SBRT-lung\\trunk\\data\\test_data\\mask_PET.nii");
                PET_mask->Update();
                typedef itk::ThresholdImageFilter <ImageType>
                  ThresholdImageFilterType;
                            	
                ThresholdImageFilterType::Pointer thresholdFilter
                  = ThresholdImageFilterType::New();
                thresholdFilter->SetInput(PET_mask->GetOutput());

                thresholdFilter->ThresholdOutside(mean_int,max_int);
                thresholdFilter->SetOutsideValue(0);
                  writer->SetFileName("C:\\Users\\sridharp\\Documents\\SBRT-lung\\trunk\\data\\test_data\\iur_PET.nii");
                  writer->SetInput(thresholdFilter->GetOutput());
                  writer->Update();
                  InputImageFileReaderPointerType PET_maskreader = ImageReaderType::New();
                  PET_maskreader->SetFileName("C:\\Users\\sridharp\\Documents\\SBRT-lung\\trunk\\data\\test_data\\iur_PET.nii");
                  PET_maskreader->Update();

                typedef itk::Image<  unsigned int , Dimension > OutputImageType;
              //  typedef itk::LabelObject< PixelType, Dimension >  LabelObjectType;
            //	  typedef itk::LabelMap< LabelObjectType >          LabelMapType;
                typedef itk::ConnectedComponentImageFilter <ImageType, OutputImageType >
                  ConnectedComponentImageFilterType;

                ConnectedComponentImageFilterType::Pointer connected =
                  ConnectedComponentImageFilterType::New();

                connected->SetInput(thresholdFilter->GetOutput());
                connected->Update();
                  std::cout << "regions" << connected->GetObjectCount();
                typedef itk::LabelImageToLabelMapFilter<OutputImageType> LabelImageToLabelMapFilterType;
                LabelImageToLabelMapFilterType::Pointer labelImageToLabelMapFilter = LabelImageToLabelMapFilterType::New();

                labelImageToLabelMapFilter->SetInput(connected->GetOutput());
                labelImageToLabelMapFilter->Update();
                  std::cout << "regions" << labelImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects();
                            	  
                typedef itk::LabelMapToLabelImageFilter< LabelImageToLabelMapFilterType::OutputImageType, OutputImageType> LabelMapToLabelImageFilterType;
                LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
                labelMapToLabelImageFilter->SetInput(labelImageToLabelMapFilter->GetOutput());
                labelMapToLabelImageFilter->Update();
                                
                            	 
                typedef itk::LabelStatisticsImageFilter< ImageType, OutputImageType > LabelStatisticsImageFilterType;
                LabelStatisticsImageFilterType::Pointer labelStatisticsImageFilter = LabelStatisticsImageFilterType::New();
                labelStatisticsImageFilter->SetLabelInput(labelMapToLabelImageFilter->GetOutput());
                  labelStatisticsImageFilter->SetInput(PET_maskreader->GetOutput());
                labelStatisticsImageFilter->Update();
                                  
                typedef LabelStatisticsImageFilterType::ValidLabelValuesContainerType ValidLabelValuesType;
                typedef LabelStatisticsImageFilterType::LabelPixelType                LabelPixelType;
                std::vector<ImageType::IndexType> foreground_seeds;
                  std::vector<ImageType::IndexType>     backgroundseeds;
                  ImageType::Pointer Pet_bound = ImageType::New();
                  Pet_bound->SetRegions(PET_imgin->GetLargestPossibleRegion());
                  Pet_bound->Allocate();
                  Pet_bound->FillBuffer(0);
                  ImageType::Pointer Pet_bound1 = ImageType::New();
                  Pet_bound1->SetRegions(PET_imgin->GetLargestPossibleRegion());
                  Pet_bound1->Allocate();
                  Pet_bound1->FillBuffer(0);
                  typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > regionFilterType;
                  regionFilterType::Pointer regionfilter = regionFilterType::New();

                  typedef itk::AddImageFilter <ImageType, ImageType >
                    AddImageFilterType;

                  AddImageFilterType::Pointer addFilter
                    = AddImageFilterType::New();

                  typedef itk::BinaryContourImageFilter <ImageType, ImageType >
                    binaryContourImageFilterType;
                  binaryContourImageFilterType::Pointer binaryContourFilter
                    = binaryContourImageFilterType::New();
                  typedef itk::PasteImageFilter <ImageType, ImageType >
                    PasteImageFilterType;
                  PasteImageFilterType::Pointer pasteFilter
                    = PasteImageFilterType::New();
                  InputImagePointerType newimg= 0;
                  newimg = PET_maskreader->GetOutput();


                  for (ValidLabelValuesType::const_iterator vIt = labelStatisticsImageFilter->GetValidLabelValues().begin();
                    vIt != labelStatisticsImageFilter->GetValidLabelValues().end();
                    ++vIt)
                  {
                    if (labelStatisticsImageFilter->HasLabel(*vIt) && *vIt > 0)
                    {
                      LabelPixelType labelValue = *vIt;

                      std::cout << "min: " << labelStatisticsImageFilter->GetMinimum(labelValue) << std::endl;
                      std::cout << "max: " << labelStatisticsImageFilter->GetMaximum(labelValue) << std::endl;
                      std::cout << "count: " << labelStatisticsImageFilter->GetCount(labelValue) << std::endl;
                      ImageType::IndexType start = labelStatisticsImageFilter->GetRegion(labelValue).GetIndex();
                      std::cout << "index: " << labelStatisticsImageFilter->GetRegion(labelValue).GetUpperIndex();
                      ImageType::SizeType size = labelStatisticsImageFilter->GetRegion(labelValue).GetSize();


                      imageCalculatorFilter->SetRegion(labelStatisticsImageFilter->GetRegion(labelValue));
                      imageCalculatorFilter->SetImage(PET_maskreader->GetOutput());
                      imageCalculatorFilter->Compute();
                      std::cout << "max of region" << imageCalculatorFilter->GetMaximum() << std::endl;
                      std::cout << "min of region" << imageCalculatorFilter->GetMinimum() << std::endl;
                      std::cout << imageCalculatorFilter->GetIndexOfMaximum() << std::endl;
                      std::cout << imageCalculatorFilter->GetIndexOfMinimum() << std::endl;

                      ImageType::RegionType desiredRegion;
                      desiredRegion.SetSize(size);
                      desiredRegion.SetIndex(start);
                                    
                                     
                      pasteFilter->SetSourceImage(PET_maskreader->GetOutput());
                      pasteFilter->SetDestinationImage(newimg);
                      pasteFilter->SetSourceRegion(labelStatisticsImageFilter->GetRegion(labelValue));
                      pasteFilter->SetDestinationIndex(start);


                      thresholdFilter->SetInput(pasteFilter->GetOutput());
                      thresholdFilter->ThresholdOutside(mean_int, imageCalculatorFilter->GetMaximum());
                      thresholdFilter->SetOutsideValue(0);

                      binaryContourFilter->SetInput(thresholdFilter->GetOutput());
                      binaryContourFilter->SetForegroundValue(0);
                      binaryContourFilter->SetBackgroundValue(255);
                      binaryContourFilter->Update();
                      Pet_bound = binaryContourFilter->GetOutput();
                      Pet_bound->SetOrigin(Pet_bound1->GetOrigin());
                      Pet_bound->SetDirection(Pet_bound1->GetDirection());
                      Pet_bound->SetSpacing(Pet_bound1->GetSpacing());

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

                    }

                    addFilter->SetInput1(Pet_bound1);
                                   
                    addFilter->SetInput2(Pet_bound);
                    addFilter->Update();
                    Pet_bound1 = addFilter->GetOutput();
                }
                  writer->SetFileName("C:\\Users\\sridharp\\Documents\\SBRT-lung\\trunk\\data\\test_data\\Pet_bound_new.nii");
                  writer->SetInput(Pet_bound1);
                  writer->Update();
  //find foreground and background seeds
}