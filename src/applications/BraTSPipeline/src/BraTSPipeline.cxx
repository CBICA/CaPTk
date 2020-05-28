#include "cbicaCmdParser.h"
#include "cbicaLogging.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"

template< class TImageType >
int algorithmsRunner()
{
  // full pipeline goes here
  
  return EXIT_SUCCESS;
}

int main(int argc, char** argv)
{
  cbica::CmdParser parser(argc, argv, "Utilities");

  parser.addOptionalParameter("i", "inputImage", cbica::Parameter::STRING, "File or Dir", "Input Image (all CaPTk supported images) for processing", "Directory to a single series DICOM only");
  parser.addOptionalParameter("o", "outputImage", cbica::Parameter::FILE, "NIfTI", "Output Image for processing");
  

  parser.addExampleUsage("-i C:/test.nii.gz -o C:/test_int.nii.gz -c int", "Cast an image pixel-by-pixel to a signed integer");

  parser.addApplicationDescription("This application performs the BraTS challenge preprocessing pipeline.");
  
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


