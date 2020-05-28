#include "cbicaCmdParser.h"
#include "cbicaLogging.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"

#include <map>

std::map< std::string, std::string > inputFiles;

std::string outputDir;

template< class TImageType >
int algorithmsRunner()
{
  // full pipeline goes here
  /*
  1.  Dicom to Nifti
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
  
  if (parser.isPresent("i"))
  {
    parser.getParameterValue("i", inputImageFile);
  }
  if (parser.isPresent("o"))
  {
    parser.getParameterValue("o", outputImageFile);
  }

  if (!cbica::isFile(inputImageFile))
  {
    std::cerr << "Input file '" << inputImageFile << "' not found.\n";
    return EXIT_FAILURE;
  }
  auto inputImageInfo = cbica::ImageInfo(inputImageFile);

  switch (inputImageInfo.GetImageDimensions())
  {
  case 2:
  {
    // this shall never be used for this application
    using ImageType = itk::Image< float, 2 >;
    return algorithmsRunner< ImageType >();

    break;
  }
  case 3:
  {
    using ImageType = itk::Image< float, 3 >;

    if (requestedAlgorithm == OrientImage) // this does not work for 2 or 4-D images
    {
      auto output = cbica::GetImageOrientation< ImageType >(cbica::ReadImage< ImageType >(inputImageFile), orientationDesired);
      std::cout << "Original Image Orientation: " << output.first << "\n";
      std::string path, base, ext;
      cbica::splitFileName(outputImageFile, path, base, ext);

      if (ext.find(".nii") != std::string::npos)
      {
        std::cerr << "WARNING: NIfTI files do NOT support orientation properly [https://github.com/InsightSoftwareConsortium/ITK/issues/1042].\n";
      }
      if (ext != ".mha")
      {
        auto tempOutputFile = path + "/" + base + ".mha"; // this is done to ensure NIfTI IO issues are taken care of
        cbica::WriteImage< ImageType >(output.second, tempOutputFile);
        auto reorientedInput = cbica::ReadImage< ImageType >(tempOutputFile);
        cbica::WriteImage< ImageType >(reorientedInput, outputImageFile);
        std::remove(tempOutputFile.c_str());
      }
      else
      {
        cbica::WriteImage< ImageType >(output.second, outputImageFile);
      }
      return EXIT_SUCCESS;
    }

    return algorithmsRunner< ImageType >();

    break;
  }
  default:
    std::cerr << "Supplied image has an unsupported dimension of '" << inputImageInfo.GetImageDimensions() << "'; only 2, 3 and 4 D images are supported.\n";
    return EXIT_FAILURE; // exiting here because no further processing should be done on the image
  }

  return EXIT_SUCCESS;
}


