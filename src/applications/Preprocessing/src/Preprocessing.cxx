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
  SusanDenoisingAlgo,
  Registration
};

// helper enum to make things smoother
enum RegistrationTypeEnum
{
  Rigid, Affine, Deformable
};

int requestedAlgorithm = 0;

std::string inputImageFile, inputMaskFile, outputImageFile, targetImageFile;
std::string registrationFixedImageFile, registrationType = "Affine", registrationMetrics = "SSD", registrationIterations = "100,50,5",
registrationAffineTransformInput, registrationDeformableTransformInput;

int histoMatchQuantiles = 40, histoMatchBins = 100,
registrationTypeInt;
bool registrationIntermediate = false, registrationSegmentationMoving = false;
float zNormCutLow = 3, zNormCutHigh = 3, zNormQuantLow = 5, zNormQuantHigh = 95, n3Bias_fwhm = 0.15,
ssSigma = 0.5, ssIntensityThreshold = 80;
int n3Bias_iterations = 50, n3Bias_fittingLevels = 4, n3Bias_otsuBins = 200, ssRadius = 1;

bool uniqueValsSort = true, boundingBoxIsotropic = true, debugMode = false;

template< class TImageType >
int algorithmsRunner()
{
  if (requestedAlgorithm == HistogramMatching)
  {
    if (debugMode)
    {
      std::cout << "Starting Histogram Matching.\n";
    }

    auto input = cbica::ReadImage< TImageType >(inputImageFile);
    auto target = cbica::ReadImage< TImageType >(targetImageFile);
    if (debugMode)
    {
      std::cout << "Finished reading input and target images.\n";
    }

    cbica::WriteImage< TImageType >(
      cbica::GetHistogramMatchedImage< TImageType >(
        input, target, histoMatchQuantiles, histoMatchBins), outputImageFile);
    std::cout << "Histogram matching completed.\n";
    return EXIT_SUCCESS;
  }

  else if (requestedAlgorithm == ZScoreNormalize)
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

  else if (requestedAlgorithm == P1P2Preprocess)
  {
    P1P2Normalizer< TImageType > normalizer;
    normalizer.SetInputImage(cbica::ReadImage< TImageType >(inputImageFile));
    normalizer.Update();
    cbica::WriteImage< TImageType >(normalizer.GetOutput(), outputImageFile);
    return EXIT_SUCCESS;
  }

  else if (requestedAlgorithm == BiasCorrectionN3)
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

  else if (requestedAlgorithm == BiasCorrectionN4)
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

  else if (requestedAlgorithm == SusanDenoisingAlgo)
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

  else if (requestedAlgorithm == Registration)
  {
    // call the greedy executable here with the proper API. 
    // see TumorGrowthModelling regarding how it is done there
    std::string greedyPathAndDim = cbica::getExecutablePath() + "greedy" + 
#if WIN32
      ".exe" +
#endif
      " -d " + std::to_string(TImageType::ImageDimension);

    std::string commonCommands; // put all the common things for the affine/deform/reslice in single place

    // add iterations to command
    std::string iterations;
    {
      auto temp = cbica::stringSplit(registrationIterations, ",");
      iterations += temp[0];
      for (size_t i = 1; i < temp.size(); i++)
      {
        iterations += "x" + temp[i];
      }
    }

    commonCommands += " -n " + iterations;

    // add mask file to command
    if (cbica::fileExists(inputMaskFile))
    {
      commonCommands += " -mm " + inputMaskFile;
    }

    // add the fixed and moving files
    commonCommands += " -i " + cbica::normPath(registrationFixedImageFile) + " " + cbica::normPath(inputImageFile);

    std::string metricsCommand = " -m ";
    if (registrationMetrics.find("NCC") != std::string::npos)
    {
      // convert 'NCC-AxBxC' to 'NCC AxBxC' for Greedy's API
      metricsCommand += cbica::stringReplace(registrationMetrics, "-", " ");
    }
    else
    {
      metricsCommand += registrationMetrics;
    }

    if (outputImageFile.empty())
    {
      std::cerr << "WARNING: Output filename is not defined; will try to save in input directory.\n";
      outputImageFile = cbica::getFilenamePath(inputImageFile) + "/registrationOutput.nii.gz";
    }
    auto outputDir = cbica::getFilenamePath(outputImageFile, false);
    auto inputFile_base = cbica::getFilenameBase(inputImageFile);
    auto fixedFile_base = cbica::getFilenameBase(registrationFixedImageFile);
    std::map< std::string, std::string > intermediateFiles;
    auto const _registrationMetrics = "_" + registrationMetrics;
    auto const _registrationMetricsNII = _registrationMetrics + ".nii.gz";
    auto const _fixedFileTOInputFileBase = fixedFile_base + "TO" + inputFile_base;
    auto const _inputFileTOFixedFileBase = inputFile_base + "TO" + fixedFile_base;

    bool defaultNamedUsed = false;
    // populate default names for intermediate files
    if (registrationAffineTransformInput.empty())
    {
      intermediateFiles["Affine"] = outputDir + "/affine_" + _fixedFileTOInputFileBase + _registrationMetrics + ".mat";
      defaultNamedUsed = true;
    }
    else
    {
      intermediateFiles["Affine"] = registrationAffineTransformInput;
    }
    if (registrationDeformableTransformInput.empty())
    {
      intermediateFiles["Deform"] = outputDir + "/deform_" + _fixedFileTOInputFileBase + _registrationMetricsNII;
      intermediateFiles["DeformInv"] = outputDir + "/deformInv_" + _inputFileTOFixedFileBase + _registrationMetricsNII;
      defaultNamedUsed = true;
    }
    else
    {
      intermediateFiles["Deform"] = registrationDeformableTransformInput;
      std::string path, base, ext;
      cbica::splitFileName(registrationDeformableTransformInput, path, base, ext);
      intermediateFiles["DeformInv"] = path + "/" + base + "-Inv.nii.gz";
    }

    std::string commandToCall;
    if (!cbica::fileExists(registrationAffineTransformInput))
    {
      // we always do affine
      commandToCall = greedyPathAndDim + " -a" +
        commonCommands +
        metricsCommand +
        " -ia-image-centers -o " + intermediateFiles["Affine"];
      ;
      if (debugMode)
      {
        std::cout << "Starting Affine registration.\n";
        std::cout << "commandToCall: \n" << commandToCall << "\n";
      }
      if (std::system(commandToCall.c_str()) != 0)
      {
        std::cerr << "Something went wrong when calling Greedy Affine.\n";
        return EXIT_FAILURE;
      }
    }
    
    switch (registrationTypeInt)
    {
    case RegistrationTypeEnum::Rigid:
    {
      // not going to be defined
      break;
    }
    case RegistrationTypeEnum::Affine:
    {
      if (registrationSegmentationMoving)
      {
        commandToCall = greedyPathAndDim +
          " -rf " + registrationFixedImageFile +
          " -rm " + inputImageFile +
          " " + outputImageFile +
          " -ri NN -r " + intermediateFiles["Affine"];
      }
      else
      {
        commandToCall = greedyPathAndDim +
          " -rf " + registrationFixedImageFile +
          " -rm " + inputImageFile +
          " " + outputImageFile +
          " -ri LABEL 0.2vox -r " + intermediateFiles["Affine"];
      }

      if (std::system(commandToCall.c_str()) != 0)
      {
        std::cerr << "Something went wrong when calling Greedy Reslice Affine.\n";
        return EXIT_FAILURE;
      }

      break;
    }
    default: // we shall always assume deformable
    {
      if (!cbica::fileExists(registrationDeformableTransformInput))
      {
        if (debugMode)
        {
          std::cout << "Starting Deformable registration.\n";
        }
        commandToCall = greedyPathAndDim +
          commonCommands +
          metricsCommand +
          " -it " + intermediateFiles["Affine"] +
          " -o " + intermediateFiles["Deform"] +
          " -oinv " + intermediateFiles["DeformInv"];

        if (std::system(commandToCall.c_str()) != 0)
        {
          std::cerr << "Something went wrong when calling Greedy Deformable.\n";
          return EXIT_FAILURE;
        }
      }

      if (registrationSegmentationMoving)
      {
        commandToCall = greedyPathAndDim +
          " -rf " + registrationFixedImageFile +
          " -rm " + inputImageFile +
          " " + outputImageFile +
          " -ri NN -r " + intermediateFiles["Deform"] + 
          " " + intermediateFiles["Affine"];
      }
      else
      {
        commandToCall = greedyPathAndDim +
          " -rf " + registrationFixedImageFile +
          " -rm " + inputImageFile +
          " " + outputImageFile +
          " -ri LABEL 0.2vox -r " + intermediateFiles["Deform"] +
          " " + intermediateFiles["Affine"];
      }

      if (std::system(commandToCall.c_str()) != 0)
      {
        std::cerr << "Something went wrong when calling Greedy Reslice Deform.\n";
        return EXIT_FAILURE;
      }

      auto outputImageFileInv = outputImageFile;
      {
        std::string path, base, ext;
        cbica::splitFileName(outputImageFileInv, path, base, ext);
        outputImageFileInv = cbica::normPath(path + "/" + base + "_inv" + ext);
      }
      if (registrationSegmentationMoving)
      {
        commandToCall = greedyPathAndDim +
          " -rf " + registrationFixedImageFile +
          " -rm " + inputImageFile +
          " " + outputImageFileInv +
          " -ri NN -r " + intermediateFiles["DeformInv"] +
          " " + intermediateFiles["Affine"] + ",-1";
      }
      else
      {
        commandToCall = greedyPathAndDim +
          " -rf " + registrationFixedImageFile +
          " -rm " + inputImageFile +
          " " + outputImageFileInv +
          " -ri LABEL 0.2vox -r " + intermediateFiles["DeformInv"] +
          " " + intermediateFiles["Affine"] + ",-1";
      }

      if (std::system(commandToCall.c_str()) != 0)
      {
        std::cerr << "Something went wrong when calling Greedy Reslice Deform-Inverse.\n";
        return EXIT_FAILURE;
      }

      break;
    }
    }

    // delete all intermediate files if the flag is not set
    if (!registrationIntermediate)
    {
      // only do the deletion if default names are used
      if (defaultNamedUsed)
      {
        for (const auto& it : intermediateFiles)
        {
          std::remove(it.second.c_str());
        }
      }
    }
  }

  return EXIT_SUCCESS;
}

//! specific implementations for the 4D images, here the images are essentially treated as series images
template< >
int algorithmsRunner< itk::Image< float, 4 > >()
{
  if (requestedAlgorithm == Registration)
  {
    using ImageTypeFloat4D = itk::Image< float, 4 >;
    using ImageTypeFloat3D = itk::Image< float, 3 >;
    // 4D series image stuff goes here
    auto movingImages = cbica::GetExtractedImages< ImageTypeFloat4D, ImageTypeFloat3D >(cbica::ReadImage< ImageTypeFloat4D >(inputImageFile));
    std::vector< ImageTypeFloat3D::Pointer > movingImage_registered;
    movingImage_registered.resize(movingImages.size());

    auto tempFolder = cbica::normPath(cbica::getFilenamePath(outputImageFile, false) + "/temp");
    cbica::createDir(tempFolder);

    auto fixedImageInfo = cbica::ImageInfo(registrationFixedImageFile);
    // if the fixed image file is also a 4D image, we shall assume the first image in the series to be the fixed image for the entire series
    if (fixedImageInfo.GetImageDimensions() == 4)
    {
      auto fixedImages = cbica::GetExtractedImages< ImageTypeFloat4D, ImageTypeFloat3D >(cbica::ReadImage< ImageTypeFloat4D >(registrationFixedImageFile));
      registrationFixedImageFile = cbica::normalizePath(tempFolder + "/fixed.nii.gz");
      cbica::WriteImage< ImageTypeFloat3D >(fixedImages[0], registrationFixedImageFile);
    }

    auto actualOutputImageFile = outputImageFile; // save the original output image file

    // perform registration for each 3D image in the 4D stack using the parameters defined by the user 
    for (size_t i = 0; i < movingImages.size(); i++)
    {
      inputImageFile = cbica::normalizePath(tempFolder + "/moving_" + std::to_string(i) + ".nii.gz");
      outputImageFile = cbica::normalizePath(tempFolder + "/output_" + std::to_string(i) + ".nii.gz");
      cbica::WriteImage< ImageTypeFloat3D >(movingImages[i], inputImageFile);
      algorithmsRunner< ImageTypeFloat3D >();
      movingImage_registered[i] = cbica::ReadImage< ImageTypeFloat3D >(outputImageFile); // store the registered image
    }
    auto finalOutput = cbica::GetJoinedImage< ImageTypeFloat3D, ImageTypeFloat4D >(movingImage_registered);
    cbica::WriteImage< ImageTypeFloat4D >(finalOutput, actualOutputImageFile);
    // delete all intermediate files if the flag is not set
    if (!registrationIntermediate)
    {
      cbica::removeDir(tempFolder); // ensure temporary folder is removed
    }
  }
  else
  {
    std::cerr << "The requested algorithm is not defined for 4D images; please retry using stacked 3D images.\n";
    return EXIT_FAILURE;
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
  parser.addOptionalParameter("reg", "registration", cbica::Parameter::STRING, "Affine | Deformable", "The kind of registration to perform", "Defaults to '" + registrationType, "Can use Mask File");
  parser.addOptionalParameter("rFI", "regFixedImg", cbica::Parameter::FILE, "NIfTI", "The Fixed Image for the registration", "Needed for registration");
  parser.addOptionalParameter("rME", "regMetrics", cbica::Parameter::STRING, "SSD | MI | NMI | NCC-AxBxC", "The kind of metris to use: SSD (Sum of Squared Differences) or MI (Mutual Information) or", "NMI (Normalized Mutual Information) or NCC-AxBxC (Normalized Cross correlation with integer radius for 3D image)", "Defaults to " + registrationMetrics);
  parser.addOptionalParameter("rNI", "regNoIters", cbica::Parameter::STRING, "N1,N2,N3", "The umber of iterations per level of multi-res", "Defaults to " + registrationIterations);
  parser.addOptionalParameter("rIS", "regInterSave", cbica::Parameter::BOOLEAN, "0 or 1", "Whether the intermediate files are to be saved or not", "Defaults to " + std::to_string(registrationIntermediate));
  parser.addOptionalParameter("rSg", "regSegMoving", cbica::Parameter::BOOLEAN, "0 or 1", "Whether the Moving Image is a segmentation file", "If 1, the 'Nearest Label' Interpolation is applied", "Defaults to " + std::to_string(registrationSegmentationMoving));
  parser.addOptionalParameter("rIA", "regInterAffn", cbica::Parameter::FILE, "mat", "The path to the affine transformation to apply to moving image", "If this is present, the Affine registration step will be skipped");
  parser.addOptionalParameter("rID", "regInterDefm", cbica::Parameter::FILE, "NIfTI", "The path to the deformable transformation to apply to moving image", "If this is present, the Deformable registration step will be skipped");

  parser.addOptionalParameter("d", "debugMode", cbica::Parameter::BOOLEAN, "0 or 1", "Enabled debug mode", "Default: 0");

  parser.addApplicationDescription("This contains all of the Preprocessing tools available in CaPTk.");
  parser.addExampleUsage("-i C:/input.nii.gz -o C:/output.nii.gz -hi C:/target.nii.gz", "Histogram Matching of 'C:/input.nii.gz' to 'C:/target.nii.gz' using the default parameters");
  parser.addExampleUsage("-i C:/input.nii.gz -o C:/output.nii.gz -zn", "Z-Score normalization of 'C:/input.nii.gz' using the default parameters");
  parser.addExampleUsage("-i C:/input.nii.gz -o C:/output.nii.gz -n3", "N3 Bias correction of 'C:/input.nii.gz' using the default parameters");
  parser.addExampleUsage("-i C:/input.nii.gz -o C:/output.nii.gz -ss", "Susan Smoothing/Denoising of 'C:/input.nii.gz' using the default parameters");
  parser.addExampleUsage("-i C:/input.nii.gz -o C:/output.nii.gz -reg Deformable -rFI C:/fixed.nii.gz", "Deformable Registration of 'C:/input.nii.gz' to 'C:/fixed.nii.gz' using the default parameters");
  parser.addExampleUsage("-i C:/input.nii.gz -o C:/output.nii.gz -reg Deformable -rFI C:/fixed.nii.gz -rIA C:/input2fixed.mat", "Deformable Registration of 'C:/input.nii.gz' to 'C:/fixed.nii.gz' using the default parameters and the affine transformation (if found) 'C:/input2fixed.mat'");
  parser.addExampleUsage("-i C:/input.nii.gz -o C:/output.nii.gz -reg Deformable -rFI C:/fixed.nii.gz -rIA C:/input2fixed.mat -rID C:/input2fixed.nii.gz", "Applies the Deformation Field (if found) in 'C:/input2fixed.nii.gz' to 'C:/input.nii.gz' after applying the affine transformation (if found) 'C:/input2fixed.mat'");
  
  if (parser.isPresent("d"))
  {
    parser.getParameterValue("d", debugMode);
  }

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
  if (parser.isPresent("hi"))
  {
    parser.getParameterValue("hi", targetImageFile);
    requestedAlgorithm = HistogramMatching;
    if (parser.isPresent("hb"))
    {
      parser.getParameterValue("hb", histoMatchBins);
    }
    if (parser.isPresent("hq"))
    {
      parser.getParameterValue("hq", histoMatchQuantiles);
    }
  }

  if (parser.isPresent("reg"))
  {
    requestedAlgorithm = Registration;

    parser.getParameterValue("reg", registrationType);
    std::transform(registrationType.begin(), registrationType.end(), registrationType.begin(), ::toupper);
    if ((registrationType.find("RIGID") != std::string::npos) || 
      (registrationType.find("AFFINE") != std::string::npos))
    {
      registrationTypeInt = RegistrationTypeEnum::Affine;
    }
    else // default case
    {
      registrationTypeInt = RegistrationTypeEnum::Deformable;
    }
    
    // detect customizations
    if (!parser.isPresent("rFI"))
    {
      std::cerr << "Registration cannot proceed without a fixed image.\n";
      return EXIT_FAILURE;
    }
    else
    {
      parser.getParameterValue("rFI", registrationFixedImageFile);
      if (!cbica::fileExists(registrationFixedImageFile))
      {
        std::cerr << "The Fixed Image was not detected, please check file name.\n";
        return EXIT_FAILURE;
      }
    }
    if (parser.isPresent("rME"))
    {
      parser.getParameterValue("rME", registrationMetrics);
      std::transform(registrationMetrics.begin(), registrationMetrics.end(), registrationMetrics.begin(), ::toupper);
    }
    if (parser.isPresent("rNI"))
    {
      parser.getParameterValue("rNI", registrationIterations);
    }
    if (parser.isPresent("rIS"))
    {
      parser.getParameterValue("rIS", registrationIntermediate);
    }
    if (parser.isPresent("rSg"))
    {
      parser.getParameterValue("rSg", registrationSegmentationMoving);
    }
    if (parser.isPresent("rIA"))
    {
      parser.getParameterValue("rIA", registrationAffineTransformInput);
    }
    if (parser.isPresent("rID"))
    {
      parser.getParameterValue("rID", registrationDeformableTransformInput);
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
  case 4:
  {
    using ImageType = itk::Image< float, 4 >;
    return algorithmsRunner< ImageType >();

    break;
  }
  default:
    std::cerr << "Supplied image has an unsupported dimension of '" << inputImageInfo.GetImageDimensions() << "'; only 2 and 3 D images are supported.\n";
    return EXIT_FAILURE; // exiting here because no further processing should be done on the image
  }

  return EXIT_SUCCESS;
}