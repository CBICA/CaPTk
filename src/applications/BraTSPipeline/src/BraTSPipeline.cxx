#include "cbicaCmdParser.h"
#include "cbicaLogging.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"

#include <map>

std::map< std::string, std::string > inputFiles;

std::string outputDir;

bool debug = true, intermediateFiles = true;

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
  /*
  1.  Dicom to Nifti
  */
  std::map< std::string, ImageType::Pointer > inputImages;

  if (debug)
  {
    std::cout << "Reading input images.\n";
  }
  for (auto it = inputFiles.begin(); it != inputFiles.end(); it++)
  {
    inputImages[it->first] = cbica::ReadImage< ImageType >(it->second);

    if (intermediateFiles)
    {
      if (debug)
      {
        std::cout << "Writing raw input (post DICOM conversion, if applicable) for modality '" << it->first << "'.\n";
      }
      cbica::WriteImage< ImageType >(inputImages[it->first], outputDir + "/" + it->first + "_raw.nii.gz");
    }
  }
  // full pipeline goes here
  /*
  2.  LPS reorientation
  3.  N4 bias correction (intermediate step)
     *   No mask
     *   shrinkFactor=4
  4.  Registration (Greedy)
     *   N4-biascorrected t1/t2/flair to N4-biascorrected t1ce, save matrix
     *   Registration of N4-biascorrected LPS t1ce to SRI, save matrix
     *   Registration of LPS t1/t1ce/t2/flair (output of step 2) to SRI space using transformation matrix saved from 4a 4b. (only 1 interpolation for  all modalities)
  5.  Generating brain mask (for BraTS, we QC here, and correct if needed)
  6.  Skull stripping of registered Images
  */

  return EXIT_SUCCESS;
}


