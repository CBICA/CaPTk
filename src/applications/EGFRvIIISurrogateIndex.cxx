#include "EGFRvIIISurrogateIndex.h"
#include "cbicaCmdParser.h"
#include "cbicaUtilities.h"
#include "cbicaITKSafeImageIO.h"

int main(int argc, char **argv)
{
  auto parser = cbica::CmdParser(argc, argv, "EGFRvIIISurrogateIndex");
  parser.addRequiredParameter("i", "image", cbica::Parameter::FILE, "NIfTI", "Input Perfusion image on which computation is done");
  parser.addRequiredParameter("m", "mask", cbica::Parameter::FILE, "NIfTI", "Mask containing near (1) and far (2) labels");
  parser.addExampleUsage("-i DSC-MRI_data.nii.gz -m Near_Far_masks.nii.gz",
    "Based on the near-far mask and input DSC MRI data, the PHI index is calculated");
  parser.addApplicationDescription("Peritumoral Heterogeneity Index calculator");
  
  std::string inputFile, drawingFile;

  parser.getParameterValue("i", inputFile);
  parser.getParameterValue("m", drawingFile);

  if (!cbica::fileExists(inputFile))
  {
    std::cerr << "Input image file not found :'" << inputFile << "'\n";
    return EXIT_FAILURE;
  }

  if (!cbica::fileExists(drawingFile))
  {
    std::cerr << "Input mask file not found :'" << drawingFile << "'\n";
    return EXIT_FAILURE;
  }

  std::cout << "Reading inputs.\n";
  using ImageType = itk::Image< float, 3 >;
  using ImageTypePerfusion = itk::Image< float, 4 >;
  auto perfusionImage = cbica::ReadImage< ImageTypePerfusion >(inputFile);
  auto mask = cbica::ReadImage< ImageType >(drawingFile);

  std::vector<ImageType::IndexType> nearIndices, farIndices;

  itk::ImageRegionIteratorWithIndex< ImageType > maskIt(mask, mask->GetLargestPossibleRegion());
  for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++maskIt)
  {
    if (maskIt.Get() == 1)
      nearIndices.push_back(maskIt.GetIndex());
    else if (maskIt.Get() == 2)
      farIndices.push_back(maskIt.GetIndex());
  }
  if (nearIndices.size()==0)
  {
    std::cerr << "Mask file does not have near indices with label=1. \n";
    return EXIT_FAILURE;
  }
  if (farIndices.size() == 0)
  {
    std::cerr << "Mask file does not have far indices with label=2. \n";
    return EXIT_FAILURE;
  }
  EGFRStatusPredictor egfrEstimator;
  auto extension = cbica::getFilenameExtension(inputFile);
  std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
  auto EGFRStatusParams = egfrEstimator.PredictEGFRStatus< ImageType, ImageTypePerfusion >(perfusionImage, nearIndices, farIndices);

  std::cout << "Printing results...\n\n";
  std::cout << "PHI Value = " << EGFRStatusParams[0] << "\n";
  std::cout << "Peak Height Ratio = " << EGFRStatusParams[1] / EGFRStatusParams[2] << "\n";
  std::cout << "Number of near voxels used = " << EGFRStatusParams[3] << "\n";
  std::cout << "Number of far voxels used = " << EGFRStatusParams[4] << "\n";

  std::cout << "Finished successfully.\n";
  return EXIT_SUCCESS;
}