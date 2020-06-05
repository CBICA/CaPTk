/**
\file cbicaBiasCorrection.h

This file holds the declaration of the class cbicaBiasCorrection.

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2020 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html
*/

#include "iostream"
#include "itkImageRegionIterator.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkExpImageFilter.h"
#include "itkN3MRIBiasFieldCorrectionImageFilter.h"
#include "itkN4BiasFieldCorrectionImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkShrinkImageFilter.h"

#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"
#include "cbicaLogging.h"
//#include "CaPTkDefines.h"


/**
* \class BiasCorrection
*
* \brief This class acts as a central access point to bias correction-related functionality.
* Supports N3 and N4 bias correction using a standard method of downscaling and reconstructing a bias
* field (as per the ANTs software).
*
* \note Default values are not provided by this class (at the moment). Please handle defaults at the caller end.
*/

class BiasCorrection
{
  // TODO: inherit from ApplicationBase for the sake of progress updates on the GUI

protected:
  cbica::Logging logger;

public:
  // default values
  static const int default_splineOrder = 3, default_otsuBins = 10, default_maxIterations = 100, default_fittingLevels = 4;
  static constexpr float default_filterNoise = 0.01, default_fwhm = 0.15;


  BiasCorrection() {};
  ~BiasCorrection() {};

  /**
  * \fn Run
  *
  * \brief This runs bias correction using the provided parameters.
  * \param correctionType String indicating which algorithm to use (n3 or n4)
  * \param inputImage Image to bias-correct
  * \param maskImage Mask to use for bias-correction
  * \param splineOrder, maxIterations, fittingLevels, filterNoise, fwhm parameters to use
  * for bias-correction
  * \param otsuBins number of bins to use when generating a thresholded mask
  * (only used when no mask is provided)
  */
  template<class TImageType, class TMaskImageType>
  typename TImageType::Pointer Run(std::string correctionType, // n3, N3, N4, or n4
    typename TImageType::Pointer inputImage,
    typename TMaskImageType::Pointer maskImage,
    int splineOrder = default_splineOrder,
    int maxIterations = default_maxIterations,
    int fittingLevels = default_fittingLevels,
    float filterNoise = default_filterNoise,
    float fwhm = default_fwhm,
    int otsuBins = default_otsuBins
  );

  /**
  * \fn Run(correctionType, inputImage)
  *
  * \brief This is a convenience function which internally calls the full Run function.
  * Use this when you don't know or care about mask characteristics.
  */
  template<class TImageType>
  typename TImageType::Pointer Run(std::string correctionType,
    typename TImageType::Pointer inputImage,
    int splineOrder = default_splineOrder,
    int maxIterations = default_maxIterations,
    int fittingLevels = default_fittingLevels,
    float filterNoise = default_filterNoise,
    float fwhm = default_fwhm,
    int otsuBins = default_otsuBins
  );

private:

    //! Overload function to configure N4 filter.
    template <typename TInput, typename TMask, typename TOutput>
    static void configureBiasCorrectionFilter(
        typename itk::N4BiasFieldCorrectionImageFilter<TInput, TMask, TOutput>::Pointer filter,
        int splineOrder,
        int maxIterations,
        int fittingLevels,
        float filterNoise,
        float fwhm,
        int otsuBins)   
    {
        filter->SetSplineOrder(splineOrder);
        filter->SetWienerFilterNoise(filterNoise); 
        filter->SetBiasFieldFullWidthAtHalfMaximum(fwhm);
        // need to pass an array of maxIterations integers (1 per fitting level)
        filter->SetNumberOfFittingLevels(fittingLevels);
        typename itk::N4BiasFieldCorrectionImageFilter<TInput, TMask, TOutput>::VariableSizeArrayType maximumNumberOfIterations(fittingLevels);
        if (!(maxIterations > fittingLevels))
        {
            maximumNumberOfIterations.Fill(1);
        }
        else // general case
        {
            maximumNumberOfIterations.Fill(maxIterations / fittingLevels);
        }
        filter->SetMaximumNumberOfIterations(maximumNumberOfIterations);
        filter->SetConvergenceThreshold(0.0000001);
    }

    //! Overload function to configure N3 filter.
    template <typename TInput, typename TMask, typename TOutput>
    static void configureBiasCorrectionFilter(
        typename itk::N3MRIBiasFieldCorrectionImageFilter<TInput, TMask, TOutput>::Pointer filter,
        int splineOrder,
        int maxIterations,
        int fittingLevels,
        float filterNoise,
        float fwhm,
        int otsuBins)
    {
        // Perform any steps needed to set options for this particular type of bias-correction
        filter->SetSplineOrder(splineOrder);
        filter->SetWeinerFilterNoise(filterNoise); // not there for n4
        filter->SetBiasFieldFullWidthAtHalfMaximum(fwhm);
        filter->SetMaximumNumberOfIterations(maxIterations); // not there for n4
        filter->SetConvergenceThreshold(0.0000001);
        filter->SetNumberOfFittingLevels(fittingLevels);
    }

  //! provide a single template function to handle the bias correction
  template< class TImageType, class TMaskImageType,
    class TBiasCorrectorType >
    typename TImageType::Pointer
    algorithmRunner(
      typename TImageType::Pointer inputImage,
      typename TMaskImageType::Pointer inputMask,
      int splineOrder,
      int maxIterations,
      int fittingLevels,
      float filterNoise,
      float fwhm,
      int otsuBins,
      std::string correctionType
    )
  {
    // Downsample image for fast bias correction as recommended by ANTs developers
    typedef itk::ShrinkImageFilter<TImageType, TImageType> ShrinkerType;
    typename ShrinkerType::Pointer shrinker = ShrinkerType::New();
    shrinker->SetInput(inputImage);
    shrinker->SetShrinkFactors(4); // 4 is an arbitrary, but good, value
    shrinker->Update();
    auto inputShrunkImage = shrinker->GetOutput();

    // Downsample the mask to match the downsampled image
    typedef itk::ShrinkImageFilter<TMaskImageType, TMaskImageType> MaskShrinkerType;
    typename MaskShrinkerType::Pointer maskShrinker = MaskShrinkerType::New();
    maskShrinker->SetInput(inputMask);
    maskShrinker->SetShrinkFactors(4);
    maskShrinker->Update();
    auto inputShrunkMask = maskShrinker->GetOutput();

    auto biasCorrector = TBiasCorrectorType::New();
    biasCorrector->SetInput(inputShrunkImage);
    biasCorrector->SetMaskImage(inputShrunkMask);

    // defer miscellaneous filter-configuration to overloaded function
    configureBiasCorrectionFilter<TImageType, TMaskImageType, TImageType>(biasCorrector, 
        splineOrder, maxIterations, fittingLevels, filterNoise, fwhm, otsuBins);
    try
    {
      biasCorrector->Update();
    }
    catch (std::exception &e)
    {
      logger.Write("Error caught while running " + correctionType + " bias correction: " + std::string(e.what()));
    }


    // Map the B-spline lattice from the filter back into the original input image's space
    typedef itk::BSplineControlPointImageFilter<typename TBiasCorrectorType::BiasFieldControlPointLatticeType,
      typename TBiasCorrectorType::ScalarImageType> BSplinerType;
    typename BSplinerType::Pointer bspliner = BSplinerType::New();
    bspliner->SetInput(biasCorrector->GetLogBiasFieldControlPointLattice());
    bspliner->SetSplineOrder(biasCorrector->GetSplineOrder());
    bspliner->SetSize(inputImage->GetLargestPossibleRegion().GetSize());
    bspliner->SetOrigin(inputImage->GetOrigin());
    bspliner->SetDirection(inputImage->GetDirection());
    bspliner->SetSpacing(inputImage->GetSpacing());
    bspliner->Update();

    auto logField = TImageType::New();
    // logField is the log of the bias field 
    logField->SetOrigin(bspliner->GetOutput()->GetOrigin());
    logField->SetSpacing(bspliner->GetOutput()->GetSpacing());
    logField->SetRegions(bspliner->GetOutput()->GetLargestPossibleRegion().GetSize());
    logField->SetDirection(bspliner->GetOutput()->GetDirection());
    logField->Allocate();

    // Copy the bspliner output into logField
    itk::ImageRegionIterator<typename TBiasCorrectorType::ScalarImageType> ItB(bspliner->GetOutput(), bspliner->GetOutput()->GetLargestPossibleRegion());
    itk::ImageRegionIterator<TImageType> ItF(logField, logField->GetLargestPossibleRegion());
    for (ItB.GoToBegin(), ItF.GoToBegin(); !ItB.IsAtEnd(); ++ItB, ++ItF)
    {
      ItF.Set(ItB.Get()[0]);
    }

    // reconstruct original bias field from its log
    typedef itk::ExpImageFilter<TImageType, TImageType> ExpFilterType;
    typename ExpFilterType::Pointer expFilter = ExpFilterType::New();
    expFilter->SetInput(logField);
    expFilter->Update();

    // Divide the input image by the reconstructed bias field (finally correcting the bias)
    typedef itk::DivideImageFilter<TImageType, TImageType, TImageType> DividerType;
    typename DividerType::Pointer divider = DividerType::New();
    divider->SetInput1(inputImage);
    divider->SetInput2(expFilter->GetOutput());
    divider->Update();

    auto outputImage = cbica::CreateImage<TImageType>(inputImage); // initialize outputImage

    itk::ImageRegionIterator<TImageType> outputIterator(outputImage, outputImage->GetLargestPossibleRegion());
    itk::ImageRegionIterator<TImageType> dividerIterator(divider->GetOutput(), divider->GetOutput()->GetLargestPossibleRegion());

    // output the newly bias-corrected image
    for (outputIterator.GoToBegin(), dividerIterator.GoToBegin(); !dividerIterator.IsAtEnd(); ++dividerIterator, ++outputIterator)
    {
      outputIterator.Set(dividerIterator.Get());
    }

    return outputImage;
  }
};

// Begin templated function definitions
template<class TImageType, class TMaskImageType>
typename TImageType::Pointer BiasCorrection::Run(std::string correctionType, // n3, N3, N4, or n4
  typename TImageType::Pointer inputImage,
  typename TMaskImageType::Pointer maskImage,
  int splineOrder,
  int maxIterations,
  int fittingLevels,
  float filterNoise,
  float fwhm,
  int otsuBins
)
{

  if (correctionType == "n3" || correctionType == "N3")
  {
    logger.Write("Selected N3 bias correction.");
  }
  else if (correctionType == "n4" || correctionType == "N4")
  {
    logger.Write("Selected N4 bias correction.");
  }
  else {
    logger.WriteError("An unsupported algorithm was selected for bias correction. Exiting.");

  }


  // Generate a mask using Otsu's thresholding if one isn't provided by the user
  if (!maskImage) // pointer to itk image is not initialized
  {
    typedef itk::OtsuThresholdImageFilter<TImageType, TMaskImageType> ThresholderType;
    typename ThresholderType::Pointer otsu = ThresholderType::New();
    otsu->SetInput(inputImage);
    otsu->SetNumberOfHistogramBins(otsuBins);
    otsu->SetInsideValue(0);
    otsu->SetOutsideValue(1);
    otsu->Update();
    maskImage = otsu->GetOutput();
  }

  typename TImageType::Pointer outputImage = cbica::CreateImage<TImageType>(inputImage);

  typename TImageType::Pointer logField = TImageType::New();
  if (correctionType == "n3" || correctionType == "N3")
  {
    using TBiasCorrectorType = itk::N3MRIBiasFieldCorrectionImageFilter< TImageType, TMaskImageType, TImageType >;

    // Run the correction on the shurnken image and mask
    outputImage = algorithmRunner<
    TImageType, TMaskImageType, TBiasCorrectorType >(
      inputImage, maskImage,
      splineOrder,
      maxIterations,
      fittingLevels,
      filterNoise,
      fwhm,
      otsuBins,
      correctionType
      );
    
  } // End of N3
  else if (correctionType == "n4" || correctionType == "N4")
  {
    using TBiasCorrectorType = itk::N4BiasFieldCorrectionImageFilter< TImageType, TMaskImageType, TImageType >;

    // Run the correction on the shurnken image and mask
    outputImage = algorithmRunner<
      TImageType, TMaskImageType, TBiasCorrectorType >(
        inputImage, maskImage,
        splineOrder,
        maxIterations,
        fittingLevels,
        filterNoise,
        fwhm,
        otsuBins,
        correctionType
        );

  } // End of N4

  return outputImage;
}

template<class TImageType>
typename TImageType::Pointer BiasCorrection::Run(std::string correctionType,
  typename TImageType::Pointer inputImage,
  int splineOrder,
  int maxIterations,
  int fittingLevels,
  float filterNoise,
  float fwhm,
  int otsuBins)
{
  return Run<TImageType, TImageType>(correctionType, inputImage, NULL, splineOrder, maxIterations, fittingLevels, filterNoise, fwhm, otsuBins);
}
