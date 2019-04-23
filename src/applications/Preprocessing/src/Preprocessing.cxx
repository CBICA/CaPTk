#include "cbicaCmdParser.h"
#include "cbicaLogging.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"

#include "ZScoreNormalizer.h"
#include "P1P2Normalizer.h"
#include "itkN3MRIBiasFieldCorrectionImageFilter.h"
#include "itkN4BiasFieldCorrectionImageFilter.h"
#include "SusanDenoising.h"

#include "itkBoundingBox.h"
#include "itkPointSet.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"

//! Detail the available algorithms to make it easier to initialize
enum AvailableAlgorithms
{
  None,
  HistogramMatching,
  ZScoreNormalize,
  P1P2Preprocess,
  BiasCorrectionN3,
  BiasCorrectionN4,
  SusanDenoisingAlgo
};

int requestedAlgorithm = 0;

std::string inputImageFile, inputMaskFile, outputImageFile, targetImageFile;
int histoMatchQuantiles = 40, histoMatchBins = 100;
float zNormCutLow = 3, zNormCutHigh = 3, zNormQuantLow = 5, zNormQuantHigh = 95, n3Bias_fwhm = 0.15,
ssSigma = 0.5, ssIntensityThreshold = 80;
int n3Bias_iterations = 50, n3Bias_fittingLevels = 4, n3Bias_otsuBins = 200, ssRadius = 1;

bool uniqueValsSort = true, boundingBoxIsotropic = true;

template< class TImageType >
int algorithmsRunner()
{
  if (requestedAlgorithm == HistogramMatching)
  {
    cbica::WriteImage< TImageType >(
      cbica::GetHistogramMatchedImage< TImageType >(
        cbica::ReadImage< TImageType >(inputImageFile), cbica::ReadImage< TImageType >(targetImageFile), histoMatchQuantiles, histoMatchBins), outputImageFile);
    std::cout << "Histogram matching completed.\n";
    return EXIT_SUCCESS;
  }

  if (requestedAlgorithm == ZScoreNormalize)
  {
    ZScoreNormalizer< TImageType > normalizer;
    normalizer.SetInputImage(cbica::ReadImage< TImageType >(inputImageFile));
    if (!inputMaskFile.empty())
    {
      normalizer.SetInputMask(cbica::ReadImage< TImageType >(inputMaskFile));
    }
    normalizer.SetCutoffs(zNormCutLow, zNormCutHigh);
    normalizer.SetQuantiles(zNormQuantLow, zNormQuantHigh);
    normalizer.Update();
    cbica::WriteImage< TImageType >(normalizer.GetOutput(), outputImageFile);
    return EXIT_SUCCESS;
  }

  if (requestedAlgorithm == P1P2Preprocess)
  {
    P1P2Normalizer< TImageType > normalizer;
    normalizer.SetInputImage(cbica::ReadImage< TImageType >(inputImageFile));
    normalizer.Update();
    cbica::WriteImage< TImageType >(normalizer.GetOutput(), outputImageFile);
    return EXIT_SUCCESS;
  }

  if (requestedAlgorithm == BiasCorrectionN3)
  {
    auto inputImage = cbica::ReadImage< TImageType >(inputImageFile);
    auto corrector = itk::N3MRIBiasFieldCorrectionImageFilter< TImageType, TImageType, TImageType >::New();
    corrector->SetInput(inputImage);
    corrector->SetMaximumNumberOfIterations(n3Bias_iterations);
    corrector->SetNumberOfFittingLevels(n3Bias_fittingLevels);
    corrector->SetBiasFieldFullWidthAtHalfMaximum(n3Bias_fwhm);
    if (!inputMaskFile.empty())
    {
      corrector->SetMaskImage(cbica::ReadImage< TImageType >(inputMaskFile));
    }
    else
    {
      auto otsu = itk::OtsuThresholdImageFilter< TImageType, TImageType >::New();
      otsu->SetInput(inputImage);
      otsu->SetNumberOfHistogramBins(n3Bias_otsuBins);
      otsu->SetInsideValue(0);
      otsu->SetOutsideValue(1);
      otsu->Update();
      corrector->SetMaskImage(otsu->GetOutput());
    }
    corrector->Update();

    cbica::WriteImage< TImageType >(corrector->GetOutput(), outputImageFile);
    return EXIT_SUCCESS;
  }

  if (requestedAlgorithm == BiasCorrectionN4)
  {
    auto inputImage = cbica::ReadImage< TImageType >(inputImageFile);
    using TBiasCorrectorType = itk::N4BiasFieldCorrectionImageFilter< TImageType, TImageType, TImageType >;
    auto corrector = itk::N4BiasFieldCorrectionImageFilter< TImageType, TImageType, TImageType >::New();
    typename itk::N4BiasFieldCorrectionImageFilter< TImageType, TImageType, TImageType >::VariableSizeArrayType iterations;
    iterations.Fill(n3Bias_iterations);
    corrector->SetInput(inputImage);
    corrector->SetMaximumNumberOfIterations(iterations);
    corrector->SetNumberOfFittingLevels(n3Bias_fittingLevels);
    corrector->SetBiasFieldFullWidthAtHalfMaximum(n3Bias_fwhm);
    if (!inputMaskFile.empty())
    {
      corrector->SetMaskImage(cbica::ReadImage< TImageType >(inputMaskFile));
    }
    else
    {
      auto otsu = itk::OtsuThresholdImageFilter< TImageType, TImageType >::New();
      otsu->SetInput(inputImage);
      otsu->SetNumberOfHistogramBins(n3Bias_otsuBins);
      otsu->SetInsideValue(0);
      otsu->SetOutsideValue(1);
      otsu->Update();
      corrector->SetMaskImage(otsu->GetOutput());
    }
    corrector->Update();

    cbica::WriteImage< TImageType >(corrector->GetOutput(), outputImageFile);

  }

  if (requestedAlgorithm == SusanDenoisingAlgo)
  {
    SusanDenoising denoiser;
    denoiser.SetSigma(ssSigma);
    denoiser.SetIntensityVariationThreshold(ssIntensityThreshold);
    denoiser.SetRadius(ssRadius);

    cbica::WriteImage< TImageType >(
      denoiser.Run< TImageType >(cbica::ReadImage< TImageType >(inputImageFile)),
      outputImageFile
      );

    return EXIT_SUCCESS;
  }

  return EXIT_SUCCESS;
}


int main(int argc, char** argv)
{
  cbica::CmdParser parser(argc, argv);

  parser.addOptionalParameter("i", "inputImage", cbica::Parameter::FILE, "NIfTI", "Input Image for processing");
  parser.addOptionalParameter("m", "maskImage", cbica::Parameter::FILE, "NIfTI", "Input Mask for processing");
  parser.addOptionalParameter("o", "outputImage", cbica::Parameter::FILE, "NIfTI", "Output Image for processing");
  parser.addOptionalParameter("hi", "histoMatch", cbica::Parameter::FILE, "NIfTI Target", "Match inputImage with the file provided in with this parameter", "Pass the target image after '-hi'");
  parser.addOptionalParameter("hb", "hMatchBins", cbica::Parameter::INTEGER, "1-1000", "Number of histogram bins for histogram matching", "Only used for histoMatching", "Defaults to " + std::to_string(histoMatchBins));
  parser.addOptionalParameter("hq", "hMatchQnts", cbica::Parameter::INTEGER, "1-1000", "Number of quantile values to match for histogram matching", "Only used for histoMatching", "Defaults to " + std::to_string(histoMatchQuantiles));
  parser.addOptionalParameter("zn", "zScoreNorm", cbica::Parameter::BOOLEAN, "N.A.", "Z-Score normalization");
  parser.addOptionalParameter("zq", "zNormQuant", cbica::Parameter::FLOAT, "0-100", "The Lower-Upper Quantile range to remove", "Default: 5,95");
  parser.addOptionalParameter("zc", "zNormCut", cbica::Parameter::FLOAT, "0-10", "The Lower-Upper Cut-off (multiple of stdDev) to remove", "Default: 3,3");
  parser.addOptionalParameter("n3", "n3BiasCorr", cbica::Parameter::STRING, "N.A.", "Runs the N3 bias correction", "Optional parameters: mask or bins, iterations, fitting levels");
  parser.addOptionalParameter("n3I", "n3BiasIter", cbica::Parameter::INTEGER, "N.A.", "Number of iterations the algorithm needs to run for", "Defaults to " + std::to_string(n3Bias_iterations));
  parser.addOptionalParameter("n3F", "n3BiasFitLevl", cbica::Parameter::INTEGER, "N.A.", "Number of fitting levels the algorithm needs obtain", "Defaults to " + std::to_string(n3Bias_fittingLevels));
  parser.addOptionalParameter("n3B", "n3BiasBins", cbica::Parameter::INTEGER, "N.A.", "If no mask is specified, N3 Bias correction makes one using Otsu", "This parameter specifies the number of histogram bins for Otsu", "Defaults to " + std::to_string(n3Bias_otsuBins));
  parser.addOptionalParameter("n3W", "n3BiasWidth", cbica::Parameter::INTEGER, "N.A.", "The full width at half maximum", "Defaults to " + std::to_string(n3Bias_fwhm));
  parser.addOptionalParameter("n4", "n4BiasCorr", cbica::Parameter::STRING, "N.A.", "Runs the N4 bias correction", "Optional parameters: mask or bins, iterations, fitting levels");
  parser.addOptionalParameter("n4I", "n4BiasIter", cbica::Parameter::INTEGER, "N.A.", "Number of iterations the algorithm needs to run for", "Defaults to " + std::to_string(n3Bias_iterations));
  parser.addOptionalParameter("n4F", "n4BiasFitLevl", cbica::Parameter::INTEGER, "N.A.", "Number of fitting levels the algorithm needs obtain", "Defaults to " + std::to_string(n3Bias_fittingLevels));
  parser.addOptionalParameter("n4B", "n4BiasBins", cbica::Parameter::INTEGER, "N.A.", "If no mask is specified, N3 Bias correction makes one using Otsu", "This parameter specifies the number of histogram bins for Otsu", "Defaults to " + std::to_string(n3Bias_otsuBins));
  parser.addOptionalParameter("n4W", "n4BiasWidth", cbica::Parameter::INTEGER, "N.A.", "The full width at half maximum", "Defaults to " + std::to_string(n3Bias_fwhm));
  parser.addOptionalParameter("ss", "susanSmooth", cbica::Parameter::STRING, "N.A.", "Susan smoothing of an image");
  parser.addOptionalParameter("ssS", "susanSigma", cbica::Parameter::FLOAT, "N.A.", "Susan smoothing Sigma", "Defaults to " + std::to_string(ssSigma));
  parser.addOptionalParameter("ssR", "susanRadius", cbica::Parameter::INTEGER, "N.A.", "Susan smoothing Radius", "Defaults to " + std::to_string(ssRadius));
  parser.addOptionalParameter("ssT", "susanThresh", cbica::Parameter::FLOAT, "N.A.", "Susan smoothing Intensity Variation Threshold", "Defaults to " + std::to_string(ssIntensityThreshold));
  parser.addOptionalParameter("p12", "p1p2norm", cbica::Parameter::STRING, "N.A.", "P1-P2 normalization required for skull stripping");
  

  if (parser.isPresent("i"))
  {
    parser.getParameterValue("i", inputImageFile);
  }
  if (parser.isPresent("m"))
  {
    parser.getParameterValue("m", inputMaskFile);
  }
  if (parser.isPresent("o"))
  {
    parser.getParameterValue("o", outputImageFile);
  }
  if (parser.isPresent("p12"))
  {
    requestedAlgorithm = P1P2Preprocess;
  }
  if (parser.isPresent("zn"))
  {
    requestedAlgorithm = ZScoreNormalize;
    std::string tempCutOff, tempQuant;
    if (parser.isPresent("zc"))
    {
      parser.getParameterValue("zc", tempCutOff);
      auto temp = cbica::stringSplit(tempCutOff, ",");
      if (temp.size() == 2)
      {
        zNormCutLow = std::atof(temp[0].c_str());
        zNormCutHigh = std::atof(temp[1].c_str());

        if (zNormCutHigh < zNormCutLow)
        {
          std::swap(zNormCutHigh, zNormCutLow);
        }
      }
    }
    if (parser.isPresent("zq"))
    {
      parser.getParameterValue("zq", tempCutOff);
      auto temp = cbica::stringSplit(tempCutOff, ",");
      if (temp.size() == 2)
      {
        zNormQuantLow = std::atof(temp[0].c_str());
        zNormQuantHigh = std::atof(temp[1].c_str());

        if (zNormQuantHigh < zNormQuantLow)
        {
          std::swap(zNormQuantHigh, zNormQuantLow);
        }
      }
    }
  }
  if (parser.isPresent("n3"))
  {
    if (parser.isPresent("n3I"))
    {
      parser.getParameterValue("n3I", n3Bias_iterations);
    }
    if (parser.isPresent("n3F"))
    {
      parser.getParameterValue("n3F", n3Bias_fittingLevels);
    }
    if (parser.isPresent("n3B"))
    {
      parser.getParameterValue("n3B", n3Bias_otsuBins);
    }
    if (parser.isPresent("n3W"))
    {
      parser.getParameterValue("n3W", n3Bias_fwhm);
    }

    requestedAlgorithm = BiasCorrectionN3;
  }

  if (parser.isPresent("n4"))
  {
    if (parser.isPresent("n4I"))
    {
      parser.getParameterValue("n4I", n3Bias_iterations);
    }
    if (parser.isPresent("n4F"))
    {
      parser.getParameterValue("n4F", n3Bias_fittingLevels);
    }
    if (parser.isPresent("n4B"))
    {
      parser.getParameterValue("n4B", n3Bias_otsuBins);
    }
    if (parser.isPresent("n4W"))
    {
      parser.getParameterValue("n4W", n3Bias_fwhm);
    }

    requestedAlgorithm = BiasCorrectionN4;
  }

  if (parser.isPresent("ss"))
  {
    requestedAlgorithm = SusanDenoisingAlgo;

    if (parser.isPresent("ssS"))
    {
      parser.getParameterValue("ssS", ssSigma);
    }
    if (parser.isPresent("ssR"))
    {
      parser.getParameterValue("ssR", ssRadius);
    }
    if (parser.isPresent("ssT"))
    {
      parser.getParameterValue("ssT", ssIntensityThreshold);
    }
  }

  auto inputImageInfo = cbica::ImageInfo(inputImageFile);

  switch (inputImageInfo.GetImageDimensions())
  {
  case 2:
  {
    using ImageType = itk::Image< float, 2 >;
    return algorithmsRunner< ImageType >();

    break;
  }
  case 3:
  {
    using ImageType = itk::Image< float, 3 >;
    return algorithmsRunner< ImageType >();

    break;
  }
  //case 4:
  //{
  //  using ImageType = itk::Image< float, 4 >;
  //  return algorithmsRunner< ImageType >();

  //  break;
  //}
  default:
    std::cerr << "Supplied image has an unsupported dimension of '" << inputImageInfo.GetImageDimensions() << "'; only 2 and 3 D images are supported.\n";
    return EXIT_FAILURE; // exiting here because no further processing should be done on the image
  }

  return EXIT_SUCCESS;
}