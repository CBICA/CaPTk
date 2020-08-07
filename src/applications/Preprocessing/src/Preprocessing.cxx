#include "cbicaCmdParser.h"
#include "cbicaLogging.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"

#include "ZScoreNormalizer.h"
#include "P1P2Normalizer.h"
//#include "itkN3MRIBiasFieldCorrectionImageFilter.h"
#include "BiasCorrection.hpp"
#include "itkN4BiasFieldCorrectionImageFilter.h"
#include "SusanDenoising.h"

#include "itkBoundingBox.h"
#include "itkPointSet.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkStatisticsImageFilter.h"

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
  Registration,
  Rescaling
};

// helper enum to make things smoother
enum RegistrationTypeEnum
{
  Rigid, Affine, Deformable
};

int requestedAlgorithm = 0;

std::string inputImageFile, inputMaskFile, outputImageFile, outputDir, targetImageFile;
std::vector< std::string > inputImageFiles; // store multiple image files
std::string registrationFixedImageFile, registrationType = "Affine", registrationMetrics = "NMI", registrationIterations = "100,50,5",
registrationAffineTransformInput, registrationDeformableTransformInput;

int histoMatchQuantiles = 40, histoMatchBins = 100,
  registrationTypeInt, registrationRigidDof = 12;
bool registrationIntermediate = false, registrationSegmentationMoving = false;
float zNormCutLow = 3, zNormCutHigh = 3, zNormQuantLow = 5, zNormQuantHigh = 95,
  bias_fwhm = BiasCorrection::default_fwhm, rescaleLower = 0, rescaleUpper = 1000,
  ssSigma = 0.5, ssIntensityThreshold = 80,
  bias_filterNoise = BiasCorrection::default_filterNoise;
int bias_splineOrder = BiasCorrection::default_splineOrder, 
  bias_otsuBins = BiasCorrection::default_otsuBins, ssRadius = 1, 
  bias_maxIterations = BiasCorrection::default_maxIterations,
  bias_fittingLevels = BiasCorrection::default_fittingLevels;
// Note: increases to bias_fittingLevels cause exponential increases in execution time, be warned

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

    if (inputImageFiles.size() > 1) // multiple images passed
    {
      std::cerr << "This operation cannot currently be performed with multiple images.\n";
      return EXIT_FAILURE;
    }
    else
    {
      auto input = cbica::ReadImage< TImageType >(inputImageFile);
      auto target = cbica::ReadImage< TImageType >(targetImageFile);
      if (debugMode)
      {
        std::cout << "Finished reading input and target images.\n";
      }

      cbica::WriteImage< TImageType >(
        cbica::GetHistogramMatchedImage< TImageType >(
          input, target, histoMatchQuantiles, histoMatchBins), outputImageFile);
    }

    std::cout << "Histogram matching completed.\n";
    return EXIT_SUCCESS;
  }

  else if (requestedAlgorithm == ZScoreNormalize)
  {
    if (inputImageFiles.size() > 1) // multiple images passed
    {
      std::vector< typename TImageType::Pointer > inputImages, inputMask;
      for (size_t i = 0; i < inputImageFiles.size(); i++)
      {
        auto currentImage = cbica::ReadImage< TImageType >(inputImageFiles[i]);
        inputImages.push_back(currentImage);
        if (!inputMaskFile.empty())
        {
          inputMask.push_back(cbica::ReadImage< TImageType >(inputMaskFile));
        }
      }

      using NewImageType = itk::Image< typename TImageType::PixelType, TImageType::ImageDimension + 1 >;
      auto combinedImage = cbica::GetJoinedImage< TImageType, NewImageType >(inputImages);
      auto combinedMask = combinedImage;
      if (!inputMask.empty())
      {
        combinedMask = cbica::GetJoinedImage< TImageType, NewImageType >(inputMask);
      }

      ZScoreNormalizer< NewImageType > normalizer;
      normalizer.SetInputImage(combinedImage);
      if (!inputMask.empty())
      {
        normalizer.SetInputMask(combinedMask);
      }
      normalizer.SetCutoffs(zNormCutLow, zNormCutHigh);
      normalizer.SetQuantiles(zNormQuantLow, zNormQuantHigh);
      normalizer.Update();

      auto combinedOutput = normalizer.GetOutput();
      auto extractedOutputs = cbica::GetExtractedImages< NewImageType, TImageType >(combinedOutput);

      for (size_t i = 0; i < extractedOutputs.size(); i++)
      {
        cbica::WriteImage< TImageType >(extractedOutputs[i], outputDir + "/" + cbica::getFilenameBase(inputImageFiles[i]) + "_zscored.nii.gz");
      }
    }
    else
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
    }
    return EXIT_SUCCESS;
  }

  else if (requestedAlgorithm == P1P2Preprocess)
  {
    if (!inputImageFiles.empty()) // multiple images passed
    {
      std::vector< typename TImageType::Pointer > inputImages, inputMask;
      for (size_t i = 0; i < inputImageFiles.size(); i++)
      {
        auto currentImage = cbica::ReadImage< TImageType >(inputImageFiles[i]);
        inputImages.push_back(currentImage);
        /// commented for future expansion
        //if (!inputMaskFile.empty())
        //{
        //  inputMask.push_back(cbica::ReadImage< TImageType >(inputMaskFile));
        //}
      }

      using NewImageType = itk::Image< typename TImageType::PixelType, TImageType::ImageDimension + 1 >;
      auto combinedImage = cbica::GetJoinedImage< TImageType, NewImageType >(inputImages);
      //auto combinedMask = combinedImage;
      /// commented for future expansion
      //if (!inputMask.empty())
      //{
      //  combinedMask = cbica::GetJoinedImage< TImageType, NewImageType >(inputMask);
      //}

      P1P2Normalizer< NewImageType > normalizer;
      normalizer.SetInputImage(combinedImage);
      normalizer.Update();

      auto combinedOutput = normalizer.GetOutput();
      auto extractedOutputs = cbica::GetExtractedImages< NewImageType, TImageType >(combinedOutput);

      for (size_t i = 0; i < extractedOutputs.size(); i++)
      {
        cbica::WriteImage< TImageType >(extractedOutputs[i], outputDir + "/" + cbica::getFilenameBase(inputImageFiles[i]) + ".nii.gz");
      }
    }
    else
    {
      P1P2Normalizer< TImageType > normalizer;
      normalizer.SetInputImage(cbica::ReadImage< TImageType >(inputImageFile));
      normalizer.Update();
      cbica::WriteImage< TImageType >(normalizer.GetOutput(), outputImageFile);
      return EXIT_SUCCESS;
    }
  }

  else if (requestedAlgorithm == BiasCorrectionN3)
  {
    if (inputImageFiles.size() > 1) // multiple images passed
    {
      std::cerr << "This operation cannot currently be performed with multiple images.\n";
      return EXIT_FAILURE;
    }
    else
    {
      auto inputImage = cbica::ReadImage<TImageType>(inputImageFile);
      typedef itk::Image<unsigned char, TImageType::ImageDimension> TMaskImageType;
      typename TMaskImageType::Pointer maskImage; // mask inits to null
      if (!inputMaskFile.empty())
      {
        maskImage = cbica::ReadImage<TMaskImageType>(inputMaskFile);
      }

      BiasCorrection biasCorrector;
      auto outputImage = biasCorrector.Run<TImageType, TMaskImageType>("n3",
        inputImage,
        maskImage,
        bias_splineOrder,
        bias_maxIterations,
        bias_fittingLevels,
        bias_filterNoise,
        bias_fwhm,
        bias_otsuBins);


      cbica::WriteImage< TImageType >(outputImage, outputImageFile);
      return EXIT_SUCCESS;
    }
  }

  else if (requestedAlgorithm == BiasCorrectionN4)
  {
    if (inputImageFiles.size() > 1) // multiple images passed
    {
      std::cerr << "This operation cannot currently be performed with multiple images.\n";
      return EXIT_FAILURE;
    }
    else
    {
      auto inputImage = cbica::ReadImage<TImageType>(inputImageFile);
      typedef itk::Image<unsigned char, TImageType::ImageDimension> TMaskImageType;
      typename TMaskImageType::Pointer maskImage; // mask inits to null
      if (!inputMaskFile.empty())
      {
        maskImage = cbica::ReadImage<TMaskImageType>(inputMaskFile);
      }

      BiasCorrection biasCorrector;
      auto outputImage = biasCorrector.Run<TImageType, TMaskImageType>("n4",
        inputImage,
        maskImage,
        bias_splineOrder,
        bias_maxIterations,
        bias_fittingLevels,
        bias_filterNoise,
        bias_fwhm,
        bias_otsuBins);


      cbica::WriteImage< TImageType >(outputImage, outputImageFile);
      return EXIT_SUCCESS;
    }
  }

  else if (requestedAlgorithm == SusanDenoisingAlgo)
  {
    if (inputImageFiles.size() > 1) // multiple images passed
    {
      std::cerr << "This operation cannot currently be performed with multiple images.\n";
      return EXIT_FAILURE;
    }
    else
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
  }

  else if (requestedAlgorithm == Registration)
  {
    if (inputImageFiles.empty()) // multiple images passed
    {
      inputImageFiles.push_back(inputImageFile);
    }

    for (size_t i = 0; i < inputImageFiles.size(); i++)
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
        if (temp.size() == 1)
        {
          temp = cbica::stringSplit(registrationIterations, "x");
        }
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
      commonCommands += " -i " + cbica::normPath(registrationFixedImageFile) + " " + cbica::normPath(inputImageFiles[i]);

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
      auto inputFile_base = cbica::getFilenameBase(inputImageFile);
      auto fixedFile_base = cbica::getFilenameBase(registrationFixedImageFile);
      std::string interimFiles_affineTransform, interimFiles_deformField, interimFiles_invDeformField;
      auto const _registrationMetrics = "_" + registrationMetrics;
      auto const _registrationMetricsNII = _registrationMetrics + ".nii.gz";
      auto const _fixedFileTOInputFileBase = fixedFile_base + "TO" + inputFile_base;
      auto const _inputFileTOFixedFileBase = inputFile_base + "TO" + fixedFile_base;

      bool affine_defaultNamedUsed = false, deformable_defaultNamedUsed = false;
      // populate default names for intermediate files
      if (registrationAffineTransformInput.empty())
      {
        interimFiles_affineTransform = outputDir + "/affine_" + _fixedFileTOInputFileBase + _registrationMetrics + ".mat";
        affine_defaultNamedUsed = true;
      }
      else
      {
        interimFiles_affineTransform = registrationAffineTransformInput;
      }
      if (registrationDeformableTransformInput.empty())
      {
        interimFiles_deformField = outputDir + "/deform_" + _fixedFileTOInputFileBase + _registrationMetricsNII;
        interimFiles_invDeformField = outputDir + "/deformInv_" + _inputFileTOFixedFileBase + _registrationMetricsNII;
        deformable_defaultNamedUsed = true;
      }
      else
      {
        interimFiles_deformField = registrationDeformableTransformInput;
        std::string path, base, ext;
        cbica::splitFileName(registrationDeformableTransformInput, path, base, ext);
        interimFiles_invDeformField = path + "/" + base + "-Inv.nii.gz";
      }

      std::string commandToCall;
      if (!cbica::fileExists(interimFiles_affineTransform))
      {
        if (registrationTypeInt == RegistrationTypeEnum::Rigid)
        {
          if (TImageType::ImageDimension == 3)
          {
            registrationRigidDof = 6;
          }
          else if (TImageType::ImageDimension == 2)
          {
            registrationRigidDof = 3;
          }
          commandToCall = greedyPathAndDim + " -a" +
            commonCommands +
            metricsCommand +
            " -ia-image-centers -dof " + std::to_string(registrationRigidDof) + " -o " + interimFiles_affineTransform;
          ;
        }
        else
        {
          commandToCall = greedyPathAndDim + " -a" +
            commonCommands +
            metricsCommand +
            " -ia-image-centers -dof " + std::to_string(registrationRigidDof) + " -o " + interimFiles_affineTransform;
          ;
        }
        if (debugMode)
        {
          std::cout << "Starting Affine/Rigid registration.\n";
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
      case RegistrationTypeEnum::Deformable:
      {
        if (!cbica::fileExists(interimFiles_deformField) || cbica::fileExists(interimFiles_invDeformField))
        {
          if (debugMode)
          {
            std::cout << "Starting Deformable registration.\n";
          }
          commandToCall = greedyPathAndDim +
            commonCommands +
            metricsCommand +
            " -it " + interimFiles_affineTransform +
            " -o " + interimFiles_deformField +
            " -oinv " + interimFiles_invDeformField;

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
            " -ri LABEL 0.2vox -r " + interimFiles_deformField +
            " " + interimFiles_affineTransform +
            " -rm " + inputImageFile +
            " " + outputImageFile;
        }
        else
        {
          commandToCall = greedyPathAndDim +
            " -rf " + registrationFixedImageFile +
            " -ri LINEAR -r " + interimFiles_deformField +
            " " + interimFiles_affineTransform +
            " -rm " + inputImageFile +
            " " + outputImageFile;
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
            " -ri NN -r " + interimFiles_invDeformField +
            " " + interimFiles_affineTransform + ",-1";
        }
        else
        {
          commandToCall = greedyPathAndDim +
            " -rf " + registrationFixedImageFile +
            " -rm " + inputImageFile +
            " " + outputImageFileInv +
            " -ri LABEL 0.2vox -r " + interimFiles_invDeformField +
            " " + interimFiles_affineTransform + ",-1";
        }

        if (std::system(commandToCall.c_str()) != 0)
        {
          std::cerr << "Something went wrong when calling Greedy Reslice Deform-Inverse.\n";
          return EXIT_FAILURE;
        }
        break;
      }
      default: // we shall always assume affine/rigid
      {
        if (registrationSegmentationMoving)
        {
          commandToCall = greedyPathAndDim +
            " -rf " + registrationFixedImageFile +
            " -ri LABEL 0.2vox -r " + interimFiles_affineTransform +
            " -rm " + inputImageFile +
            " " + outputImageFile;
        }
        else
        {
          commandToCall = greedyPathAndDim +
            " -rf " + registrationFixedImageFile +
            " -ri LINEAR -r " + interimFiles_affineTransform +
            " -rm " + inputImageFile +
            " " + outputImageFile;
        }

        if (std::system(commandToCall.c_str()) != 0)
        {
          std::cerr << "Something went wrong when calling Greedy Reslice Affine.\n";
          return EXIT_FAILURE;
        }

        break;
      }
      }

      // delete all intermediate files if the flag is not set
      if (!registrationIntermediate)
      {
        // only do the deletion if default names are used
        if (affine_defaultNamedUsed)
        {
          std::remove(interimFiles_affineTransform.c_str());
        }
        if (deformable_defaultNamedUsed)
        {
          std::remove(interimFiles_deformField.c_str());
          std::remove(interimFiles_invDeformField.c_str());
        }
      }
    }
  }

  else if (requestedAlgorithm == Rescaling)
  {
    typename TImageType::PixelType minimum = std::numeric_limits< typename TImageType::PixelType >::min(),
      maximum = std::numeric_limits< typename TImageType::PixelType >::max();

    if (inputImageFiles.size() > 1) // multiple input images passed
    {
      std::vector< typename TImageType::Pointer > inputImages;
      for (size_t i = 0; i < inputImageFiles.size(); i++)
      {
        auto currentImage = cbica::ReadImage< TImageType >(inputImageFiles[i]);
        inputImages.push_back(currentImage);
      }

      using NewImageType = itk::Image< typename TImageType::PixelType, TImageType::ImageDimension + 1 >;
      auto combinedImage = cbica::GetJoinedImage< TImageType, NewImageType >(inputImages);

      auto rescaler = itk::RescaleIntensityImageFilter< NewImageType >::New();
      rescaler->SetInput(combinedImage);
      rescaler->SetOutputMaximum(rescaleUpper);
      rescaler->SetOutputMinimum(rescaleLower);
      try
      {
        rescaler->Update();
      }
      catch (const std::exception&e)
      {
        std::cerr << "Error caught during rescaling: " << e.what() << "\n";
        return EXIT_FAILURE;
      }

      auto outputImages = cbica::GetExtractedImages< NewImageType, TImageType >(rescaler->GetOutput());

      // at this point, we have found the global minimum and maximum
      auto fileEnding = "_rescaled-" + std::to_string(rescaleLower) + "-" + std::to_string(rescaleUpper) + ".nii.gz";
      for (size_t i = 0; i < outputImages.size(); i++)
      {
        cbica::WriteImage< TImageType >(outputImages[i], outputDir + "/" + cbica::getFilenameBase(inputImageFiles[i]) + fileEnding);
      }
    }
    else
    {
      auto inputImage = cbica::ReadImage< TImageType >(inputImageFile);
      auto statsCalculator = itk::StatisticsImageFilter< TImageType >::New();
      statsCalculator->SetInput(inputImage);

      auto rescaler = itk::RescaleIntensityImageFilter< TImageType >::New();
      rescaler->SetInput(inputImage);
      rescaler->SetOutputMaximum(rescaleUpper);
      rescaler->SetOutputMinimum(rescaleLower);

      cbica::WriteImage< TImageType >(rescaler->GetOutput(), outputImageFile);
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
  parser.addOptionalParameter("zq", "zNormQuant", cbica::Parameter::STRING, "0-100", "The Lower-Upper Quantile range to remove", "Default: 5,95");
  parser.addOptionalParameter("zc", "zNormCut", cbica::Parameter::STRING, "0-10", "The Lower-Upper Cut-off (multiple of stdDev) to remove", "Default: 3,3");
  parser.addOptionalParameter("n3", "n3BiasCorr", cbica::Parameter::STRING, "N.A.", "Runs the N3 bias correction",
                              "Optional parameters: mask or bins, spline order, filter noise level, fitting levels, max iterations, full-width-at-half-maximum");
  parser.addOptionalParameter("n4", "n4BiasCorr", cbica::Parameter::STRING, "N.A.", "Runs the N4 bias correction",
                              "Optional parameters: mask or bins, spline order, filter noise level, fitting levels, full-width-at-half-maximum");
  parser.addOptionalParameter("nS", "nSplineOrder", cbica::Parameter::INTEGER, "N.A.", "The spline order for the bias correction", "Defaults to " + std::to_string(bias_splineOrder));
  parser.addOptionalParameter("nF", "nFilterNoise", cbica::Parameter::FLOAT, "N.A.", "The filter noise level for the bias correction", "Defaults to " + std::to_string(bias_filterNoise));
  parser.addOptionalParameter("nB", "nBiasBins", cbica::Parameter::INTEGER, "N.A.", "If no mask is specified, N3/N4 bias correction makes one using Otsu", "This parameter specifies the number of histogram bins for Otsu", "Defaults to " + std::to_string(bias_otsuBins));
  parser.addOptionalParameter("nFL", "nFittingLevels", cbica::Parameter::INTEGER, "N.A.", "The number of fitting levels to use for bias correction", "Defaults to " + std::to_string(bias_fittingLevels));
  parser.addOptionalParameter("nMI", "nMaxIterations", cbica::Parameter::INTEGER, "N.A.", "The maximum number of iterations for bias correction (only works for N3)", "Defaults to " + std::to_string(bias_maxIterations));
  parser.addOptionalParameter("nFWHM", "nFullWidthHalfMaximum", cbica::Parameter::INTEGER, "N.A.", "Set the full-width-at-half-maximum value for bias correction", "Defaults to " + std::to_string(bias_fwhm));
  parser.addOptionalParameter("ss", "susanSmooth", cbica::Parameter::STRING, "N.A.", "Susan smoothing of an image");
  parser.addOptionalParameter("ssS", "susanSigma", cbica::Parameter::FLOAT, "N.A.", "Susan smoothing Sigma", "Defaults to " + std::to_string(ssSigma));
  parser.addOptionalParameter("ssR", "susanRadius", cbica::Parameter::INTEGER, "N.A.", "Susan smoothing Radius", "Defaults to " + std::to_string(ssRadius));
  parser.addOptionalParameter("ssT", "susanThresh", cbica::Parameter::FLOAT, "N.A.", "Susan smoothing Intensity Variation Threshold", "Defaults to " + std::to_string(ssIntensityThreshold));
  parser.addOptionalParameter("p12", "p1p2norm", cbica::Parameter::STRING, "N.A.", "P1-P2 normalization required for skull stripping");
  parser.addOptionalParameter("reg", "registration", cbica::Parameter::STRING, "Affine-DOF | Deformable | Rigid", "The kind of registration to perform", "Defaults to " + registrationType, "Can use Mask File with '-m' and multiple moving images with '-i'", "For Affine, the second number defines the degrees of freedom, eg: '-ref Affine-12'");
  parser.addOptionalParameter("rFI", "regFixedImg", cbica::Parameter::FILE, "NIfTI", "The Fixed Image for the registration", "Needed for registration");
  parser.addOptionalParameter("rME", "regMetrics", cbica::Parameter::STRING, "SSD | MI | NMI | NCC-AxBxC", "The kind of metrics to use: SSD (Sum of Squared Differences) or MI (Mutual Information) or", "NMI (Normalized Mutual Information) or NCC-AxBxC (Normalized Cross correlation with integer radius for 3D image)", "Defaults to " + registrationMetrics);
  parser.addOptionalParameter("rNI", "regNoIters", cbica::Parameter::STRING, "N1,N2,N3", "The number of iterations per level of multi-res", "Defaults to " + registrationIterations);
  parser.addOptionalParameter("rIS", "regInterSave", cbica::Parameter::BOOLEAN, "0 or 1", "Whether the intermediate files are to be saved or not", "Defaults to " + std::to_string(registrationIntermediate));
  parser.addOptionalParameter("rSg", "regSegMoving", cbica::Parameter::BOOLEAN, "0 or 1", "Whether the Moving Image(s) is a segmentation file", "If 1, the 'Nearest Label' Interpolation is applied", "Defaults to " + std::to_string(registrationSegmentationMoving));
  parser.addOptionalParameter("rIA", "regInterAffn", cbica::Parameter::FILE, "mat", "The path to the affine transformation to apply to moving image", "If this is present, the Affine registration step will be skipped", "Also used for rigid transformation");
  parser.addOptionalParameter("rID", "regInterDefm", cbica::Parameter::FILE, "NIfTI", "The path to the deformable transformation to apply to moving image", "If this is present, the Deformable registration step will be skipped");
  parser.addOptionalParameter("rsc", "rescaleImage", cbica::Parameter::STRING, "Output Intensity range", "The output intensity range after image rescaling", "Defaults to " + std::to_string(rescaleLower) + ":" + std::to_string(rescaleUpper), "If multiple inputs are passed (comma-separated), the rescaling is done in a cumulative manner,", "i.e., stats from all images are considered for the scaling");

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
    if (inputImageFile.find(",") != std::string::npos) // multiple images are getting passed;
    {
      inputImageFiles = cbica::stringSplit(inputImageFile, ",");
      for (size_t n = 0; n < inputImageFiles.size(); n++) // perform sanity check
      {
        if (!cbica::ImageSanityCheck(inputImageFiles[0], inputImageFiles[1]))
        {
          std::cerr << "The images cannot be processed together since they aren't defined in the same physical space.\n";
          return EXIT_FAILURE;
        }
      }
    }
    else
    {
      inputImageFiles.push_back(inputImageFile);
    }
  }
  if (parser.isPresent("m"))
  {
    parser.getParameterValue("m", inputMaskFile);
  }
  if (parser.isPresent("o"))
  {
    parser.getParameterValue("o", outputImageFile);
    outputDir = cbica::getFilenamePath(outputImageFile, false);
    cbica::createDir(outputDir);
  }

  // parse all options from here
  if (parser.isPresent("p12"))
  {
    requestedAlgorithm = P1P2Preprocess;
  }
  else if (parser.isPresent("zn"))
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
  else if (parser.isPresent("n3"))
  {
    if (parser.isPresent("nS"))
    {
      parser.getParameterValue("nS", bias_splineOrder);
    }
    if (parser.isPresent("nF"))
    {
      parser.getParameterValue("nF", bias_filterNoise);
    }
    if (parser.isPresent("nB"))
    {
      parser.getParameterValue("nB", bias_otsuBins);
    }
    if (parser.isPresent("nFL"))
    {
      parser.getParameterValue("nFL", bias_fittingLevels);
    }
    if (parser.isPresent("nMI"))
    {
      parser.getParameterValue("nMI", bias_maxIterations);
    }
    if (parser.isPresent("nFWHM"))
    {
      parser.getParameterValue("nFWHM", bias_fwhm);
    }

    requestedAlgorithm = BiasCorrectionN3;
  }
  else if (parser.isPresent("n4"))
  {
    if (parser.isPresent("nS"))
    {
      parser.getParameterValue("nS", bias_splineOrder);
    }
    if (parser.isPresent("nF"))
    {
      parser.getParameterValue("nF", bias_filterNoise);
    }
    if (parser.isPresent("nB"))
    {
      parser.getParameterValue("nB", bias_otsuBins);
    }
    if (parser.isPresent("nFL"))
    {
      parser.getParameterValue("nFL", bias_fittingLevels);
    }
    //if (parser.isPresent("nMI")) // This doesn't work for N4 (currently).
    //{
    //  parser.getParameterValue("nMI", bias_maxIterations);
    //}
    if (parser.isPresent("nFWHM"))
    {
      parser.getParameterValue("nFWHM", bias_fwhm);
    }

    requestedAlgorithm = BiasCorrectionN4;
  }
  else if (parser.isPresent("ss"))
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
  else if (parser.isPresent("hi"))
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
  else if (parser.isPresent("reg"))
  {
    requestedAlgorithm = Registration;

    parser.getParameterValue("reg", registrationType);
    std::transform(registrationType.begin(), registrationType.end(), registrationType.begin(), ::toupper);
    if (registrationType.find("RIGID") != std::string::npos)
    {
      registrationTypeInt = RegistrationTypeEnum::Rigid;
    }
    else if (registrationType.find("AFFINE") != std::string::npos)
    {
      registrationTypeInt = RegistrationTypeEnum::Affine;
      auto temp = cbica::stringSplit(registrationType, "-");
      // check for different delimiters
      if (temp.size() == 1)
      {
        temp = cbica::stringSplit(registrationType, "x");
        if (temp.size() == 1)
        {
          temp = cbica::stringSplit(registrationType, ":");
        }
      }
      if (temp.size() == 2)
      {
        registrationRigidDof = std::atoi(temp[1].c_str());
        if ((registrationRigidDof != 12) && (registrationRigidDof != 6) && (registrationRigidDof != 7))
        {
          std::cerr << "Greedy only accepts 6, 7 or 12 as the affine DOF.\n";
          return EXIT_FAILURE;
        }
      }
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
      //std::transform(registrationMetrics.begin(), registrationMetrics.end(), registrationMetrics.begin(), ::toupper);
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
  else if (parser.isPresent("rsc"))
  {
    requestedAlgorithm = Rescaling;
    std::string temp;
    parser.getParameterValue("rsc", temp);
    auto delimitersToCheck = { ":", ",", "x" };
    for (auto delimIter = delimitersToCheck.begin(); delimIter != delimitersToCheck.end(); ++delimIter)
    {
      if (temp.find(*delimIter) != std::string::npos)
      {
        auto bounds = cbica::stringSplit(temp, *delimIter);
        rescaleLower = std::atof(bounds[0].c_str());
        rescaleUpper = std::atof(bounds[1].c_str());

        if (rescaleUpper < rescaleLower)
        {
          std::swap(rescaleLower, rescaleUpper);
        }
        break;
      }
    }
  }

  if (!inputImageFiles.empty()) // if multiple images are passed, set up the algorithm runner using the first image
  {
    inputImageFile = inputImageFiles[0]; 
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
