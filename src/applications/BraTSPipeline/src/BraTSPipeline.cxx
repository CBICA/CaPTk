#include "cbicaCmdParser.h"
#include "cbicaLogging.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"

#include "BiasCorrection.hpp"

#include <map>

std::map< std::string, std::string > inputFiles;

std::string outputDir;

bool debug = true, intermediateFiles = true;

//! didn't want to include qt for this function
std::string getDataDir()
{
  auto currentAppPath = cbica::getExecutablePath();
  auto captk_dataDir = currentAppPath + "../data/";
  if (!cbica::exists(captk_dataDir))
  {
    captk_dataDir = currentAppPath + "../../data/";
    if (!cbica::exists(captk_dataDir))
    {
      captk_dataDir = currentAppPath + "../Resources/data/";
      if (!cbica::exists(captk_dataDir))
      {
        std::cerr << "Data Directory not found. Please re-install")
        return "";
      }
    }
  }

  return captk_dataDir;
}

int main(int argc, char** argv)
{
  cbica::CmdParser parser(argc, argv, "Utilities");

  parser.addRequiredParameter("t1c", "t1ceImage", cbica::Parameter::STRING, "Input Image (DICOM or NIfTI)", "Input structural T1-weighted post-contrast image");
  parser.addRequiredParameter("t1", "t1Image", cbica::Parameter::STRING, "Input Image (DICOM or NIfTI)", "Input structural T1-weighted pre-contrast image");
  parser.addRequiredParameter("t2", "t2Image", cbica::Parameter::STRING, "Input Image (DICOM or NIfTI)", "Input structural T2-weighted contrast image");
  parser.addRequiredParameter("fl", "flImage", cbica::Parameter::STRING, "Input Image (DICOM or NIfTI)", "Input structural FLAIR contrast image");
  parser.addRequiredParameter("o", "outputDir", cbica::Parameter::DIRECTORY, "Directory", "Output directory for final output");
  parser.addOptionalParameter("d", "debug", cbica::Parameter::BOOLEAN, "0 or 1", "Print debugging information", "Defaults to 1");
  parser.addOptionalParameter("i", "interFiles", cbica::Parameter::BOOLEAN, "0 or 1", "Save intermediate files", "Defaults to 1");

  parser.addExampleUsage("-t1c C:/input/t1ce/image.dcm -t1 C:/input/t1/image.dcm -t2 C:/input/t2/image.dcm -fl C:/input/flair/image.dcm -o C:/input/output", "Run full BraTS pipeline for specified DICOM images");
  parser.addExampleUsage("-t1c C:/input/t1ce.nii.gz -t1 C:/input/t1.nii.gz -t2 C:/input/t2.nii.gz -fl C:/input/flair.nii.gz -o C:/input/output", "Run full BraTS pipeline for specified (raw) NIfTI images");

  parser.addApplicationDescription("This application performs the BraTS challenge preprocessing pipeline.");

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
      if (intermediateFiles)
      {
        if (debug)
        {
          std::cout << "Writing re-oriented image for modality '" << modality << "'.\n";
        }
        cbica::WriteImage< ImageType >(inputImages_processed[modality], outputDir + "/" + modality + "_rai.nii.gz");
      }
    }

    /// [3] N4 bias correction

    if (debug)
    {
      std::cout << "Starting bias correction for modality '" << modality << "'.\n";
    }

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
    // the bias-corrected images need to be written because these are passed on to greedy
    cbica::WriteImage< ImageType >(inputImages_processed[modality], outputDir + "/" + modality + "_rai_n4.nii.gz");
  } // end inputFiles iterator
  

  /// [4] Registration using Greedy
  if (debug)
  {
    std::cout << "Registering T1CE to SRI atlas.\n";
  }

  auto greedyPathAndDim = cbica::getExecutablePath() + "greedy" +
#if WIN32
    ".exe" +
#endif
    " -d 3";

  auto atlasImage = getDataDir() + "/sri24/atlastImage.nii.gz";
  auto image_t1ce = outputDir + "/T1CE_rai_n4.nii.gz";

  auto fullCommand = " -a -m NMI -i " + atlasImage + " " + image_t1ce + " -o " + outputDir + "/t1ceToSRI.mat -ia-image-centers -n 100x50x10 -dof 6";

  if (std::system((greedyPathAndDim + fullCommand).c_str()) != 0)
  {
    std::cerr << "Something went wrong when registering T1CE image to SRI atlas, please re-try or contact sofware@cbica.upenn.edu.\n";
    return EXIT_FAILURE;
  }

  // perform greedy to register image_t1ce to SRI24 atlas and save matrix

  ///cbica/home/sakoc/software/greedy/v1.0.1/build/greedy -d 3 -a -m NMI 
  // -i /scratch/braintumoruser/rGreedy_Rigid.xP9PqORqxJ/atlas.nii.gz 
  // /scratch/braintumoruser/rGreedy_Rigid.xP9PqORqxJ/input.nii.gz 
  // -o /scratch/braintumoruser/rGreedy_Rigid.xP9PqORqxJ/affine_matrix.mat -ia-image-centers -n 100x50x10 -dof 6

  for (auto it = inputFiles.begin(); it != inputFiles.end(); it++)
  {
    auto modality = it->first;
    if (modality != "T1CE")
    {
      auto image_current = outputDir + "/" + modality + "_rai_n4.nii.gz";

      // perform greedy to register image_current to image_t1ce and save matrix
    }
  }
  /*
  4.  Registration (Greedy)
     *   N4-biascorrected t1/t2/flair to N4-biascorrected t1ce, save matrix
     *   Registration of N4-biascorrected LPS t1ce to SRI, save matrix
     *   Registration of LPS t1/t1ce/t2/flair (output of step 2) to SRI space using transformation matrix saved from 4a 4b. (only 1 interpolation for  all modalities)
  5.  Generating brain mask (for BraTS, we QC here, and correct if needed)
  6.  Skull stripping of registered Images
  */

  return EXIT_SUCCESS;
}


