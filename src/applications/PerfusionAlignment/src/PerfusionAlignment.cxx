#include "cbicaCmdParser.h"
#include "cbicaLogging.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"

#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"

std::string inputImageFile, outputImageFile;

int main(int argc, char** argv)
{
  cbica::CmdParser parser(argc, argv);

  parser.addRequiredParameter("i", "inputImage", cbica::Parameter::FILE, "NIfTI", "Input Image for processing");
  parser.addRequiredParameter("o", "outputImage", cbica::Parameter::FILE, "NIfTI", "Output Image for processing");

  if (parser.isPresent("i"))
  {
    parser.getParameterValue("i", inputImageFile);
  }
  if (parser.isPresent("o"))
  {
    parser.getParameterValue("o", outputImageFile);
  }

  auto inputImageInfo = cbica::ImageInfo(inputImageFile);

  switch (inputImageInfo.GetImageDimensions())
  {
  case 4:
  {
    using ImageType = itk::Image< float, 4 >;

    auto inputImage = cbica::ReadImage< ImageType >(inputImageFile);

    auto outputSize = inputImage->GetLargestPossibleRegion().GetSize();
    auto outputSpacing = inputImage->GetSpacing();

    outputSpacing[3] = 1;

    // update the output image size
    outputSize[3] = outputSize[3] * outputSpacing[3];

    auto resampler = itk::ResampleImageFilter< ImageType, ImageType >::New();
    resampler->SetInput(inputImage);
    resampler->SetSize(outputSize);
    resampler->SetOutputSpacing(outputSpacing);
    resampler->SetOutputOrigin(inputImage->GetOrigin());
    resampler->SetOutputDirection(inputImage->GetDirection());
    resampler->SetOutputStartIndex(inputImage->GetLargestPossibleRegion().GetIndex());
    resampler->SetTransform(itk::IdentityTransform< double, ImageType::ImageDimension >::New());
    resampler->UpdateLargestPossibleRegion();

    cbica::WriteImage< ImageType >(resampler->GetOutput(), outputImageFile);

    break;
  }
  default:
    std::cerr << "Supplied image has an unsupported dimension of '" << inputImageInfo.GetImageDimensions() << "'; only 4-D images are supported.\n";
    return EXIT_FAILURE; // exiting here because no further processing should be done on the image
  }

  return EXIT_SUCCESS;
}