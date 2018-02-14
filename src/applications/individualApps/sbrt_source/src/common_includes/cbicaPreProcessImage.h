/**
\file  cbicaPreProcessImage.h

\brief Pre-process an input image

Currently, the implementation is valid ONLY for 3D images. 

The operations being done are:

Histogram matching using windowed filter

http://www.cbica.upenn.edu/sbia/software/ <br>
sbia-software@uphs.upenn.edu

Copyright (c) 2015 University of Pennsylvania. All rights reserved. <br>
See COPYING file or http://www.cbica.upenn.edu/sbia/software/license.html

*/
#pragma once

#include "itkImage.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkImageToHistogramFilter.h"
#include "itkIntensityWindowingImageFilter.h"


namespace cbica
{

  template <class TImageType>
  class PreProcessImage
  {
  public:
    /**
    \brief The Constructor
    */
    explicit PreProcessImage()
    {
      m_binsPerDimension = 32; // hard coded
      m_output = TImageType::New();
      m_input  = TImageType::New();
    }

    /**
    \brief The Default Constructor
    */
    PreProcessImage<TImageType>(const PreProcessImage<TImageType>&);

    /**
    \brief The Destructor
    */
    virtual ~PreProcessImage(){ };


    /**
    \brief Set the Input Image

    \param inputImage Pointer to an ITK Image
    */
    void SetInputImage(const typename TImageType::Pointer inputImage)
    {
      //m_input = inputImage;
      m_input->Graft(inputImage);
      m_input->DisconnectPipeline();
    }

    /**
    \brief Run the algorithm
    */
    void Update()
    {
      typedef itk::Statistics::ImageToHistogramFilter<TImageType> HistogramFilterType;
      typedef typename HistogramFilterType::InputBooleanObjectType InputBooleanObjectType;
      typedef typename HistogramFilterType::HistogramSizeType HistogramSizeType;

      HistogramSizeType histogramSize(m_binsPerDimension); // hard coded
      histogramSize[0] = 256;

      typename InputBooleanObjectType::Pointer autoMinMaxInputObject = InputBooleanObjectType::New();
      autoMinMaxInputObject->Set(true);

      typename HistogramFilterType::Pointer histogramFilter = HistogramFilterType::New();
      histogramFilter->SetInput(m_input);
      histogramFilter->SetAutoMinimumMaximumInput(autoMinMaxInputObject);
      histogramFilter->SetHistogramSize(histogramSize);
      histogramFilter->SetMarginalScale(10.0);
      histogramFilter->Update();

      // set lower and upper limits of windowing filter
      typename TImageType::PixelType winsorizeLowerQuantile = 0.0;
      typename TImageType::PixelType winsorizeUpperQuantile = 1.0;

      // range of intensities to map windowing filter output 
      typename TImageType::PixelType lowerScaleValue = static_cast<typename TImageType::PixelType>(0);
      typename TImageType::PixelType upperScaleValue = static_cast<typename TImageType::PixelType>(1);

      // obtain the quantile of the histogram which defines the minimum and maximum intensities
      typename TImageType::PixelType lowerValue = histogramFilter->GetOutput()->Quantile(0, winsorizeLowerQuantile);
      typename TImageType::PixelType upperValue = histogramFilter->GetOutput()->Quantile(0, winsorizeUpperQuantile);

      //std::cout << "LowerValue = " << lowerValue << "\nUpper Value = " << upperValue << "\n";

      typedef itk::IntensityWindowingImageFilter<TImageType, TImageType> IntensityWindowingImageFilterType;
      typename IntensityWindowingImageFilterType::Pointer windowingFilter = IntensityWindowingImageFilterType::New();
      windowingFilter->SetInput(m_input);
      windowingFilter->SetWindowMinimum(lowerValue);
      windowingFilter->SetWindowMaximum(upperValue);
      windowingFilter->SetOutputMinimum(lowerScaleValue);
      windowingFilter->SetOutputMaximum(upperScaleValue);
      windowingFilter->Update();

      m_output = windowingFilter->GetOutput();
      m_output->Update();
      m_output->DisconnectPipeline();
    }

    /**
    \brief Get the output

    \return An Image pointers containing the output of the pipeline
    */
    typename TImageType::Pointer GetOutput()
    {
      return m_output;
    }

  private:
    typename TImageType::Pointer m_input; // input image
    unsigned int m_binsPerDimension;
    typename TImageType::Pointer m_output;
  };

}
