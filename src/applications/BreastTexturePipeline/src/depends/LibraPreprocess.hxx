#pragma once

//#include "LibraPreprocess.h"

#include "itkFlipImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkSquareImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLogImageFilter.h"
#include "itkAbsoluteValueDifferenceImageFilter.h"
#include "itkAddImageFilter.h"

#include "cbicaITKSafeImageIO.h"
#include "cbicaITKUtilities.h"

using LibraImageType = itk::Image< float, 2 >;

template< typename TImageType >
LibraPreprocess< TImageType >::LibraPreprocess()
{

}

template< typename TImageType >
LibraPreprocess< TImageType >::~LibraPreprocess()
{

}

template< typename TImageType >
void LibraPreprocess< TImageType >::SetInputFileName(const std::string &input)
{
  m_fileName = input;
  m_intermediateOutputDir = cbica::getFilenamePath(m_fileName) + "cxx_";
}

template< typename TImageType >
void LibraPreprocess< TImageType >::SetResizingFactor(size_t resize)
{
  m_resizingFactor = resize;
}

template< typename TImageType >
void LibraPreprocess< TImageType >::SetResizingFactor(float resize)
{
  if (resize <= 1)
  {
    SetResizingFactor(static_cast<size_t>(resize * 100));
  }
}

template< typename TImageType >
void LibraPreprocess< TImageType >::EnableDebugMode()
{
  m_debugMode = true;
}

template< typename TImageType >
void LibraPreprocess< TImageType >::Update()
{
  if (!m_algorithmDone)
  {
    auto dicomReader = /*typename*/ itk::ImageSeriesReader< itk::Image<float,2> >::New();
    dicomReader->SetImageIO(itk::GDCMImageIO::New());
    dicomReader->SetFileName(m_fileName);
    try
    {
      dicomReader->Update();
    }
    catch (itk::ExceptionObject & err)
    {
      std::cerr << "Error while loading DICOM image(s): " << err.what() << "\n";
    }
    m_inputImage = dicomReader->GetOutput();
    m_output = m_inputImage;
    m_output->DisconnectPipeline();

    if (m_debugMode)
    {
      cbica::WriteImage< TImageType >(m_inputImage, m_intermediateOutputDir + "inputLoaded.nii.gz");
    }

    if (m_resizingFactor < 100) // only do a decrease, not an increase
    {
      std::string interpolator = "bicubic";
      m_output = cbica::ResizeImage< TImageType >(m_output, m_resizingFactor, interpolator); // TBD: unsure if "bicubic" is the same as "bspline" or not
      m_intermediateOutputDir += "resize-" + std::to_string(m_resizingFactor) + "-" + interpolator + "_";
      cbica::WriteImage< TImageType >(m_output, m_intermediateOutputDir + "inputLoaded.nii.gz");
    }

    auto inputDict = (*(dicomReader->GetMetaDataDictionaryArray()))[0];
    std::string presentationIntent, patientOrientation, photometricInterpretation, manufacturer, viewPosition, 
      imageLaterality, laterality, seriesTime, fieldOfViewHorzFlip, side = "Special-View", origin, codeMeaning, seriesDescription;

    itk::ExposeMetaData<std::string>(*inputDict, "0008|0068", presentationIntent);
    itk::ExposeMetaData<std::string>(*inputDict, "0020|0020", patientOrientation);
    itk::ExposeMetaData<std::string>(*inputDict, "0018|7034", fieldOfViewHorzFlip);
    itk::ExposeMetaData<std::string>(*inputDict, "0028|0004", photometricInterpretation);
    itk::ExposeMetaData<std::string>(*inputDict, "0008|0070", manufacturer);
    itk::ExposeMetaData<std::string>(*inputDict, "0018|5101", viewPosition);
    itk::ExposeMetaData<std::string>(*inputDict, "0020|0060", laterality);
    itk::ExposeMetaData<std::string>(*inputDict, "0020|0062", imageLaterality);
    itk::ExposeMetaData<std::string>(*inputDict, "0008|0031", seriesTime);
    itk::ExposeMetaData<std::string>(*inputDict, "0018|7030", origin);
    itk::ExposeMetaData<std::string>(*inputDict, "0008|0104", codeMeaning);
    itk::ExposeMetaData<std::string>(*inputDict, "0008|103E", seriesDescription);

    // If the ImageLaterality field isn't there for some reason, check the Laterality field and copy over
    if (laterality.empty())
    {
      laterality = imageLaterality;
    }
    else if (imageLaterality.empty())
    {
      imageLaterality = laterality;
    }

    m_flipAxes.Fill(false);
    auto flipFilter = itk::FlipImageFilter < TImageType >::New();
    flipFilter->SetInput(m_output);

    //// standardize orientation of DICOM -- based on MATLAB code
    //if (m_debugMode)
    //{
    //  std::cout << "Orientation fix based on MATLAB code.\n";
    //}
    //auto isFlipped = fieldOfViewHorzFlip.find("YES") != std::string::npos;
    //bool sideL = true; // check if this is the left side
    //if (imageLaterality.find("R") != std::string::npos)
    //{
    //  side = "Right";
    //  sideL = false;
    //}
    //else if (imageLaterality.find("L") != std::string::npos)
    //{
    //  side = "Left";
    //}

    //// perform flip if imageLaterality is "L" and is flipped OR if imageLaterality is "R"
    //if ((isFlipped && sideL) || !sideL) 
    //{
    //  m_flipAxes[0] = true;
    //  m_flipAxes[1] = false;
    //  flipFilter->Setm_flipAxes(m_flipAxes);
    //  flipFilter->Update();
    //  m_output = flipFilter->GetOutput();
    //}
    //m_dicomTags["Side"] = side;

    /// new code based on the texture pipeline
    auto orientationAxes = cbica::stringSplit(patientOrientation, "\\");

    if (std::strcmp(orientationAxes[0].c_str(), "P") == 0)
    {
      m_flipAxes[0] = true;
      orientationAxes[0] = "A"; // no idea why the DICOM header information needs to be changed but whatever
      if (std::strcmp(fieldOfViewHorzFlip.c_str(), "YES") == 0)
      {
        fieldOfViewHorzFlip = "NO";
      }
      else
      {
        fieldOfViewHorzFlip = "YES";
      }
    }

    if (std::strcmp(orientationAxes[1].c_str(), "H") == 0)
    {
      m_flipAxes[1] = true;
      orientationAxes[1] = "F"; // no idea why the DICOM header information needs to be changed but whatever
    }
    else if (std::strcmp(orientationAxes[1].c_str(), "L") == 0)
    {
      m_flipAxes[1] = true;
      orientationAxes[1] = "R"; // no idea why the DICOM header information needs to be changed but whatever
    }

    flipFilter->Setm_flipAxes(m_flipAxes);
    flipFilter->Update();
    m_output = flipFilter->GetOutput();
    m_output->DisconnectPipeline();
    /// new code based on the texture pipeline

    m_dicomTags["PresentationIntent"] = presentationIntent;
    m_dicomTags["PatientOrientation"] = patientOrientation;
    m_dicomTags["PatientOrientation"] = fieldOfViewHorzFlip;
    m_dicomTags["PhotometricInterpretation"] = photometricInterpretation;
    m_dicomTags["Manufacturer"] = manufacturer;
    m_dicomTags["ViewPosition"] = viewPosition;
    m_dicomTags["Laterality"] = laterality;
    m_dicomTags["ImageLaterality"] = imageLaterality;
    m_dicomTags["SeriesTime"] = seriesTime;
    m_dicomTags["CodeMeaning"] = codeMeaning;
    m_dicomTags["SeriesDescription"] = seriesDescription;

    if (m_debugMode)
    {
      cbica::WriteImage< TImageType >(m_output, m_intermediateOutputDir + "orStandardized.nii.gz");
    }

    /// standardize intensities
    auto imageCalculatorFilter = itk::MinimumMaximumImageCalculator< TImageType >::New();
    imageCalculatorFilter->SetImage(m_output);
    imageCalculatorFilter->Compute();

    // Raw needs to be log transformed, inverted and squared for preprocessing, in that order
    bool isRaw = false;
    if ((presentationIntent.find("FOR PROCESSING") != std::string::npos) ||
      (seriesDescription.find("Raw") != std::string::npos))
    {
      isRaw = true;
    }
    // Step_1: raw gets log transform
    if (isRaw)
    {
      if (m_debugMode)
      {
        std::cout << "Preprocessing raw mammogram...\n";
      }
      auto min = imageCalculatorFilter->GetMinimum();
      if (min < 1)
      {
        if (m_debugMode)
        {
          std::cout << "Adding [min + 1] to entire image...\n";
        }
        if (m_debugMode)
        {
          cbica::WriteImage< TImageType >(m_output, m_intermediateOutputDir + "intStandardized_adder-input.nii.gz");
        }
        auto adder = itk::AddImageFilter< TImageType >::New();
        adder->SetInput1(m_output);
        adder->SetConstant2(std::abs(min) + 1);
        adder->Update();
        m_output = adder->GetOutput();
        m_output->DisconnectPipeline();
        if (m_debugMode)
        {
          cbica::WriteImage< TImageType >(m_output, m_intermediateOutputDir + "intStandardized_adder.nii.gz");
        }
      }

      // performing log
      if (m_debugMode)
      {
        std::cout << "Performing log...\n";
      }
      auto logImage = itk::LogImageFilter< TImageType, TImageType >::New();
      logImage->SetInput(m_output);
      logImage->Update();
      m_output = logImage->GetOutput();
      m_output->DisconnectPipeline();
      if (m_debugMode)
      {
        cbica::WriteImage< TImageType >(m_output, m_intermediateOutputDir + "intStandardized_adder-log.nii.gz");
      }
    }

    // Step_2: Do pixel inversion if Needed
    if (photometricInterpretation.find("MONOCHROME1") != std::string::npos)
    {
      if (m_debugMode)
      {
        std::cout << "Inverting Pixel Intensities...\n";
      }
      imageCalculatorFilter->SetImage(m_output);
      imageCalculatorFilter->Compute();

      auto absoluteDiff = itk::AbsoluteValueDifferenceImageFilter< TImageType, TImageType, TImageType >::New();
      absoluteDiff->SetInput1(m_output);
      absoluteDiff->SetInput2(cbica::CreateImage< TImageType >(m_output, imageCalculatorFilter->GetMaximum()));
      absoluteDiff->Update();
      m_output = absoluteDiff->GetOutput();
      m_output->DisconnectPipeline();
      if (m_debugMode)
      {
        cbica::WriteImage< TImageType >(m_output, m_intermediateOutputDir + "intStandardized_invert.nii.gz");
      }
    }

    // Step_3: square transform for density contrast
    if (isRaw)
    {
      if (m_debugMode)
      {
        std::cout << "Square-Transform the intensities...\n";
      }

      auto squareImageFilter = itk::SquareImageFilter< TImageType, TImageType >::New();
      squareImageFilter->SetInput(m_output);
      squareImageFilter->Update();
      m_output = squareImageFilter->GetOutput();
      m_output->DisconnectPipeline();
      if (m_debugMode)
      {
        cbica::WriteImage< TImageType >(m_output, m_intermediateOutputDir + "intStandardized_square.nii.gz");
      }
    }

    if (m_debugMode)
    {
      cbica::WriteImage< TImageType >(m_output, m_intermediateOutputDir + "intStandardized.nii.gz");
    }

  }
  m_algorithmDone = true;
}

template< typename TImageType >
typename TImageType::Pointer LibraPreprocess< TImageType >::GetOutputImage()
{
  if (!m_algorithmDone)
  {
    Update();
  }
  return m_output;
}

template< typename TImageType >
std::map< std::string, std::string > LibraPreprocess< TImageType >::GetRelevantDicomTags()
{
  if (!m_algorithmDone)
  {
    Update();
  }
  return m_dicomTags;
}

template< typename TImageType >
typename TImageType::Pointer LibraPreprocess< TImageType >::GetInputImage()
{
  if (!m_algorithmDone)
  {
    Update();
  }
  return m_inputImage;
}