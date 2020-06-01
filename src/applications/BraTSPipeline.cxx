#include "cbicaCmdParser.h"
#include "cbicaLogging.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"
#include "CaPTkGUIUtils.h"

#include "BiasCorrection.hpp"

#include <map>

std::map< std::string, std::string > inputFiles;

std::string outputDir;

bool debug = true, intermediateFiles = true;

int main(int argc, char** argv)
{
  cbica::CmdParser parser(argc, argv, "BraTSPipeline");

  parser.addRequiredParameter("t1c", "t1ceImage", cbica::Parameter::STRING, "Input Image (DICOM or NIfTI)", "Input structural T1-weighted post-contrast image");
  parser.addRequiredParameter("t1", "t1Image", cbica::Parameter::STRING, "Input Image (DICOM or NIfTI)", "Input structural T1-weighted pre-contrast image");
  parser.addRequiredParameter("t2", "t2Image", cbica::Parameter::STRING, "Input Image (DICOM or NIfTI)", "Input structural T2-weighted contrast image");
  parser.addRequiredParameter("fl", "flImage", cbica::Parameter::STRING, "Input Image (DICOM or NIfTI)", "Input structural FLAIR contrast image");
  parser.addRequiredParameter("o", "outputDir", cbica::Parameter::DIRECTORY, "Directory", "Output directory for final output");
  parser.addOptionalParameter("d", "debug", cbica::Parameter::BOOLEAN, "0 or 1", "Print debugging information", "Defaults to 1");
  parser.addOptionalParameter("i", "interFiles", cbica::Parameter::BOOLEAN, "0 or 1", "Save intermediate files", "Defaults to 1");

  parser.addExampleUsage("-t1c C:/input/t1ce/image.dcm -t1 C:/input/t1/image.dcm -t2 C:/input/t2/image.dcm -fl C:/input/flair/image.dcm -o C:/input/output", "Run full BraTS pipeline for specified DICOM images");
  parser.addExampleUsage("-t1c C:/input/t1ce.nii.gz -t1 C:/input/t1.nii.gz -t2 C:/input/t2.nii.gz -fl C:/input/flair.nii.gz -o C:/input/output", "Run full BraTS pipeline for specified (raw) NIfTI images");

  parser.addApplicationDescription("This application performs the BraTS challenge preprocessing pipeline. Please delete contents of output directory or fresh run.");

  parser.getParameterValue("t1c", inputFiles["T1CE"]);
  parser.getParameterValue("t1", inputFiles["T1"]);
  parser.getParameterValue("t2", inputFiles["T2"]);
  parser.getParameterValue("fl", inputFiles["FL"]);
  parser.getParameterValue("o", outputDir);

  cbica::createDir(outputDir);

  if (parser.isPresent("d"))
  {
    parser.getParameterValue("d", debug);
  }
  if (parser.isPresent("i"))
  {
    parser.getParameterValue("i", intermediateFiles);
  }
  
  if (debug)
  {
    std::cout << "Performing sanity checks for input images.\n";
  }
  // sanity checks
  for (auto it = inputFiles.begin(); it != inputFiles.end(); it++)
  {
    if (!cbica::exists(it->second))
    {
      std::cerr << "Couldn't find the modality '" << it->first << "', denoted by '" << it->second << "'.\n";
      return EXIT_FAILURE;
    }

    auto inputImageInfo = cbica::ImageInfo(it->second);
    if (inputImageInfo.GetImageDimensions() != 3)
    {
      std::cerr << "The BraTS pipeline is only valid for 3D images, whereas the image '" 
        << it->second << "' for modality '" << it->first << "' has " <<
        inputImageInfo.GetImageDimensions() << " dimentions.\n";
      return EXIT_FAILURE;
    }
  }

  using ImageType = itk::Image< float, 3 >; // default image type

  // variables to store various images
  std::map< std::string, ImageType::Pointer > inputImages, inputImages_processed;

  // default names
  std::map< std::string, std::string > outputNames,
    inputReorientedFiles, inputReorientedBiasFiles, // filenames for reoriented and bias-corrected files
    outputMatFiles, outputRegisteredImages; // filenames for  matrices and images

  if (debug)
  {
    std::cout << "Reading input images.\n";
  }
  for (auto it = inputFiles.begin(); it != inputFiles.end(); it++)
  {
    auto modality = it->first;
    /// [1] read image - DICOM to NIfTI conversion, if applicable
    inputImages[modality] = cbica::ReadImage< ImageType >(it->second);

    if (inputImages[modality].IsNotNull())
    {
      if (intermediateFiles)
      {
        if (debug)
        {
          std::cout << "Writing raw input (post DICOM conversion, if applicable) for modality '" << modality << "'.\n";
        }
        cbica::WriteImage< ImageType >(inputImages[modality], outputDir + "/" + modality + "_raw.nii.gz");
      }
    }
    else
    {
      if (cbica::IsDicom(it->second))
      {
        std::cerr << "Something went wrong with the DICOM to NIfTI conversion for modality '" <<
          modality << "' with filename '" << it->second << "'"
          << ", please use another package to conver to NIfTI and try again.\n";
        return EXIT_FAILURE;
      }
      else
      {
        std::cerr << "Something went wrong with reading the raw input image, please re-try or contact sofware@cbica.upenn.edu.\n";
        return EXIT_FAILURE;
      }
    }

    /// [2] LPS/RAI re-orientation
    if (debug)
    {
      std::cout << "Performing re-orientation to LPS for modality '" << modality << "'.\n";
    }
    inputImages_processed[modality] = cbica::GetImageOrientation(inputImages[modality], "RAI").second;
    if (inputImages_processed[modality].IsNull())
    {
      std::cerr << "Something went wrong with re-orienting the input image, please re-try or contact sofware@cbica.upenn.edu.\n";
      return EXIT_FAILURE;
    }
    else
    {
      inputReorientedFiles[modality] = outputDir + "/" + modality + "_rai.nii.gz";
      // the re-oriented images need to be written because these are passed on to greedy
      cbica::WriteImage< ImageType >(inputImages_processed[modality], inputReorientedFiles[modality]);
    }

    /// [3] N4 bias correction

    if (debug)
    {
      std::cout << "Starting bias correction for modality '" << modality << "'.\n";
    }

    // the bias-corrected images need to be written because these are passed on to greedy
    inputReorientedBiasFiles[modality] = outputDir + "/" + modality + "_rai_n4.nii.gz";

    if (!cbica::fileExists(inputReorientedBiasFiles[modality]))
    {
      BiasCorrection biasCorrector;
      {
        using MaskImageType = itk::Image<unsigned char, ImageType::ImageDimension>;
        typename MaskImageType::Pointer maskImage; // mask inits to null
        auto outputImage = biasCorrector.Run<TImageType, MaskImageType>("n4",
          inputImages_processed[modality],
          maskImage,
          BiasCorrection::default_splineOrder,
          BiasCorrection::default_maxIterations,
          BiasCorrection::default_fittingLevels,
          BiasCorrection::default_filterNoise,
          BiasCorrection::default_fwhm,
          BiasCorrection::default_otsuBins);
        if (outputImage.IsNotNull())
        {
          inputImages_processed[modality] = outputImage;
          inputImages_processed[modality]->DisconnectPipeline();
        }
        else
        {
          std::cerr << "Something went wrong with bias-correcting the re-oriented image, please re-try or contact sofware@cbica.upenn.edu.\n";
          return EXIT_FAILURE;
        }
      }
      cbica::WriteImage< ImageType >(inputImages_processed[modality], inputReorientedBiasFiles[modality]);
    }

    outputNames[modality] = modality + "_to_SRI"; // all output names can be controlled from here
  } // end inputFiles iterator
  
  /// [4] Registration using Greedy
  // we do T1CE to Atlas registration first because other registrations are dependent on this
  if (debug)
  {
    std::cout << "Registering T1CE to SRI atlas.\n";
  }

  auto greedyPathAndDim = getApplicationPath("greedy") + " -d 3";

  auto captkDataDir = getCaPTkDataDir();
  auto atlasImage = captkDataDir + "/sri24/atlastImage.nii.gz";
  outputMatFiles["T1CE"] = outputDir + "/" + outputNames["T1CE"] + ".mat";
  outputRegisteredImages["T1CE"] = outputDir + "/" + outputNames["T1CE"] + ".nii.gz";

  std::string fullCommand;

  if (!cbica::exists(outputMatFiles["T1CE"]))
  {
    fullCommand = " -a -m NMI -i " + atlasImage + " " + inputReorientedBiasFiles["T1CE"]
      + " -o " + outputMatFiles["T1CE"] + " -ia-image-centers -n 100x50x10 -dof 6";

    if (debug)
    {
      std::cout << "Greedy command: " << greedyPathAndDim + fullCommand << "\n";
    }
    if (std::system((greedyPathAndDim + fullCommand).c_str()) != 0)
    {
      std::cerr << "Something went wrong when registering T1CE image to SRI atlas, please re-try or contact sofware@cbica.upenn.edu.\n";
      return EXIT_FAILURE;
    }
  } // end outputMatFiles["T1CE"] check

  if (!cbica::exists(outputRegisteredImages["T1CE"]))
  {
    fullCommand = " -rf " + atlasImage + " -ri LINEAR -rm " +
      inputReorientedFiles["T1CE"] + " " + outputRegisteredImages["T1CE"] + " -r " +
      outputMatFiles["T1CE"];

    if (debug)
    {
      std::cout << "Greedy command: " << greedyPathAndDim + fullCommand << "\n";
    }

    if (std::system((greedyPathAndDim + fullCommand).c_str()) != 0)
    {
      std::cerr << "Something went wrong when applying registration matrix to generate T1CE image in SRI atlas space, please re-try or contact sofware@cbica.upenn.edu.\n";
      return EXIT_FAILURE;
    }
  } // end outputRegisteredImages["T1CE"] check

  for (auto it = inputFiles.begin(); it != inputFiles.end(); it++)
  {
    auto modality = it->first;
    if (modality != "T1CE") // T1CE registration has happened before
    {
      outputMatFiles[modality] = outputDir + "/" + outputNames[modality] + ".mat";
      outputRegisteredImages[modality] = outputDir + "/" + outputNames[modality] + ".nii.gz";

      if (!cbica::exists(outputMatFiles[modality]))
      {
        // we use the bias-corrected image for registration as it is easier localize transformations
        fullCommand = " -a -m NMI -i " + inputReorientedBiasFiles["T1CE"] + " " + inputReorientedBiasFiles[modality]
          + " -o " + outputMatFiles[modality] + " -ia-image-centers -n 100x50x10 -dof 6";
        if (debug)
        {
          std::cout << "Registering " << modality << " to T1CE.\n";
          std::cout << "Greedy command: " << greedyPathAndDim + fullCommand << "\n";
        }

        if (std::system((greedyPathAndDim + fullCommand).c_str()) != 0)
        {
          std::cerr << "Something went wrong when registering " << modality
            << "to T1CE image, please re-try or contact sofware@cbica.upenn.edu.\n";
          return EXIT_FAILURE;
        }
      } // end outputMatFiles[modality] check

      if (debug)
      {
        std::cout << "Generating image for " << modality << " registered to the atlas.\n";
      }

      if (!cbica::exists(outputRegisteredImages[modality]))
      {
        // the final registration is applied on the original image after re-orientation (not bias-corrected) to
        // ensure maximum fidelity with original image
        fullCommand = " -rf " + atlasImage + " -ri LINEAR -rm " + inputReorientedFiles[modality] + " " +
          outputRegisteredImages[modality] + " -r "
          + outputMatFiles["T1CE"] + " "
          + outputMatFiles[modality];

        if (debug)
        {
          std::cout << "Greedy command: " << greedyPathAndDim + fullCommand << "\n";
        }

        if (std::system((greedyPathAndDim + fullCommand).c_str()) != 0)
        {
          std::cerr << "Something went wrong when applying registration matrix to generate " << modality << " image in SRI atlas space, please re-try or contact sofware@cbica.upenn.edu.\n";
          return EXIT_FAILURE;
        }
      } // end outputRegisteredImages[modality] check
    } // end modality check
  } // end modality loop

  /// [5] Skull-stripping using DeepMedic
  if (debug)
  {
    std::cout << "Starting skull-stripping using DeepMedic.\n";
  }

  if (!cbica::exists(outputDir + "/dmOut/brainMask.nii.gz"))
  {
    auto deepMedicExe = getApplicationPath("DeepMedic");

    fullCommand = " -md " + captkDataDir + "/deepMedic/saved_models/skullStripping/ " +
      "-t1c " + outputRegisteredImages["T1CE"] + " -t1 " +
      outputRegisteredImages["T1"] + " -t2 " +
      outputRegisteredImages["T2"] + " -fl " +
      outputRegisteredImages["FL"] + " -o " +
      outputDir + "/dmOut/brainMask.nii.gz";

    if (debug)
    {
      std::cout << "Command for DeepMedic: " << deepMedicExe + fullCommand << "\n";
    }

    if (std::system((deepMedicExe + fullCommand).c_str()) != 0)
    {
      std::cerr << "Something went wrong when performing skull-stripping using DeepMedic, please re-try or contact sofware@cbica.upenn.edu.\n";
      return EXIT_FAILURE;
    }
  }

  std::cout << "Finished, please perform manual quality-check of generated brain mask before applying to input images.\n";

  return EXIT_SUCCESS;
}


