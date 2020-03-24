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
#include "cbicaLogging.h"
#include "CaPTkDefines.h"


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
                                int splineOrder,
                                int maxIterations,
                                int fittingLevels,
                                float filterNoise,
                                float fwhm,
                                int otsuBins
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
                                int splineOrder,
                                int maxIterations,
                                int fittingLevels,
                                float filterNoise,
                                float fwhm,
                                int otsuBins);

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

    // Ensure the future output has the same properties as the input
    const typename TImageType::SpacingType spacing = inputImage->GetSpacing();
    const typename TImageType::PointType origin = inputImage->GetOrigin();
    const typename TImageType::DirectionType direction = inputImage->GetDirection();
    const typename TImageType::RegionType region = inputImage->GetLargestPossibleRegion();
    typename TImageType::Pointer outputImage = TImageType::New();

    outputImage->SetSpacing(spacing);
    outputImage->SetLargestPossibleRegion(region);
    outputImage->SetBufferedRegion(region);
    outputImage->Allocate();
    outputImage->SetOrigin(origin);
    outputImage->SetDirection(direction);


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


    typename TImageType::Pointer logField = TImageType::New();
    if (correctionType == "n3" || correctionType == "N3")
    {
        // Downsample image for fast bias correction as recommended by ANTs developers
        typedef itk::ShrinkImageFilter<TImageType, TImageType> ShrinkerType;
        typename ShrinkerType::Pointer shrinker = ShrinkerType::New();
        shrinker->SetInput(inputImage);
        shrinker->SetShrinkFactors(4); // 4 is an arbitrary, but good, value
        shrinker->Update();
        
        // Downsample the mask to match the downsampled image
        typedef itk::ShrinkImageFilter<TMaskImageType, TMaskImageType> MaskShrinkerType;
        typename MaskShrinkerType::Pointer maskshrinker = MaskShrinkerType::New();
        maskshrinker->SetInput(maskImage);
        maskshrinker->SetShrinkFactors(4);
        maskshrinker->Update();

        // Run the correction on the shurnken image and mask
        typedef itk::N3MRIBiasFieldCorrectionImageFilter<TImageType, TMaskImageType, TImageType> CorrectorType;
        typename CorrectorType::Pointer corrector = CorrectorType::New();
        corrector->SetInput(shrinker->GetOutput());
        corrector->SetMaskImage(maskshrinker->GetOutput());
        corrector->SetSplineOrder(splineOrder);
        corrector->SetWeinerFilterNoise(filterNoise);
        corrector->SetBiasFieldFullWidthAtHalfMaximum(fwhm);
        corrector->SetMaximumNumberOfIterations(maxIterations);
        corrector->SetConvergenceThreshold(0.0000001);
        corrector->SetNumberOfFittingLevels(fittingLevels); // increasing this exponentially increases execution time
        try
        {
            corrector->Update();
        }
        catch (std::exception &e)
        {
            logger.Write("Error caught while running " + correctionType +  " bias correction: " + std::string(e.what()));
        }

        // Map the B-spline lattice from the N3 filter back into the original input image's space
        typedef itk::BSplineControlPointImageFilter<typename CorrectorType::BiasFieldControlPointLatticeType, typename CorrectorType::ScalarImageType> BSplinerType;
        typename BSplinerType::Pointer bspliner = BSplinerType::New();
        bspliner->SetInput(corrector->GetLogBiasFieldControlPointLattice());
        bspliner->SetSplineOrder(corrector->GetSplineOrder());
        bspliner->SetSize(inputImage->GetLargestPossibleRegion().GetSize());
        bspliner->SetOrigin(inputImage->GetOrigin());
        bspliner->SetDirection(inputImage->GetDirection());
        bspliner->SetSpacing(inputImage->GetSpacing());
        bspliner->Update();

        // logField is the log of the bias field
        logField->SetOrigin(bspliner->GetOutput()->GetOrigin());
        logField->SetSpacing(bspliner->GetOutput()->GetSpacing());
        logField->SetRegions(bspliner->GetOutput()->GetLargestPossibleRegion().GetSize());
        logField->SetDirection(bspliner->GetOutput()->GetDirection());
        logField->Allocate();

        // Copy the bspliner output into logField
        typedef itk::N3MRIBiasFieldCorrectionImageFilter<TImageType, TMaskImageType, TImageType> CorrectorType;
        itk::ImageRegionIterator<typename CorrectorType::ScalarImageType> ItB(bspliner->GetOutput(), bspliner->GetOutput()->GetLargestPossibleRegion());
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


        itk::ImageRegionIterator<TImageType> outputIterator(outputImage, outputImage->GetLargestPossibleRegion());
        itk::ImageRegionIterator<TImageType> dividerIterator(divider->GetOutput(), corrector->GetOutput()->GetLargestPossibleRegion());

        // output the newly bias-corrected image
       for (outputIterator.GoToBegin(), dividerIterator.GoToBegin(); !dividerIterator.IsAtEnd(); ++dividerIterator, ++outputIterator)
       {
         outputIterator.Set(dividerIterator.Get());
       }
        
       return outputImage;
        
    } // End of N3

    else if (correctionType == "n4" || correctionType == "N4")
    {
        // N4 currently doesn't use the downsamplng method used in N3.
        // TODO: Apply the standard downsample-correct-reconstruct methodology to both N3 and N4.
        // When this was last attempted, a bug arose where the BSplineControlPointImageFilter couldn't
        // accept a LogBiasFieldControlPointLattice from the N4 correction filter.  
        typedef itk::N4BiasFieldCorrectionImageFilter<TImageType, TMaskImageType, TImageType> CorrectorType;
        typename CorrectorType::Pointer corrector = CorrectorType::New();
        corrector->SetInput(inputImage);
        corrector->SetMaskImage(maskImage);
        corrector->SetSplineOrder(splineOrder);
        corrector->SetWienerFilterNoise(filterNoise);
        corrector->SetBiasFieldFullWidthAtHalfMaximum(fwhm);
        // maxIterations doesn't work currently for N4 -- need to pass an array of maxIterations integers (1 per fitting level)
        //corrector->SetMaximumNumberOfIterations(maxIterations);
        corrector->SetConvergenceThreshold(0.0000001);
        // fittingLevels are also disabled (see above comment) for the time being 
        //corrector->SetNumberOfFittingLevels(fittingLevels);
        try
        {
            corrector->Update();
        }
        catch (std::exception &e)
        {
            logger.Write("Error caught while running " + correctionType +  " bias correction: " + std::string(e.what()));
        }

        itk::ImageRegionIterator<TImageType> outputIterator(outputImage, outputImage->GetLargestPossibleRegion());
        itk::ImageRegionIterator<TImageType> correctedIterator(corrector->GetOutput(), corrector->GetOutput()->GetLargestPossibleRegion());

        // output the newly bias-corrected image
       for (outputIterator.GoToBegin(), correctedIterator.GoToBegin(); !correctedIterator.IsAtEnd(); ++correctedIterator, ++outputIterator)
       {
         outputIterator.Set(correctedIterator.Get());
       }
       
       return outputImage;
        
    } // End of N4


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
