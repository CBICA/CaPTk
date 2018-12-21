#include <algorithm>

#include "WhiteStripe.h"

#include "cbicaCmdParser.h"
#include "cbicaUtilities.h"
#include "cbicaITKSafeImageIO.h"

//#include "CAPTk.h"

int main(int argc, char **argv)
{
  auto parser = cbica::CmdParser(argc, argv);
  parser.addRequiredParameter("i", "inputImages", cbica::Parameter::FILE, "3D NIfTI Image(s)", "Input Images on which WhiteStripe needs to be applied", "Delieanted by ','");
  parser.addRequiredParameter("o", "outputDir", cbica::Parameter::DIRECTORY, "Dir with write access", "Output directory where results are to be saved");
  parser.addOptionalParameter("r", "radius", cbica::Parameter::FLOAT, "0.0 to 5.0", "WhiteStripe Radius", "Default: 0.05");
  parser.addOptionalParameter("sk", "skullStrippedImg", cbica::Parameter::BOOLEAN, "1 or 0", "Whether skull stripped image is passed", "If not, 'zSliceRange' is needed", "Default: 1");
  parser.addOptionalParameter("z", "zSliceRange", cbica::Parameter::INTEGER, "50 to 150", "z-slice Range for cropping", "Delieanted by '-'", "Default (if image is not skull-stripped): 80-120");
  parser.addOptionalParameter("t", "tissuesMax", cbica::Parameter::INTEGER, "5|10|20", "Max Tissues (5 or 10 or 20)", "Default: 5");
  parser.addOptionalParameter("m", "maxSmooth", cbica::Parameter::FLOAT, "0.0 to 50.0", "Max Smoothing", "Default: 10.0");
  parser.addOptionalParameter("d", "deltaSmooth", cbica::Parameter::FLOAT, "0.0 to 10.0", "Smoothing Delta", "Default: 0.5");
  parser.addOptionalParameter("b", "binsHist", cbica::Parameter::INTEGER, "100 to 3000", "Number of Histogram bins to do processing", "Default: 2000");
  parser.addOptionalParameter("t1", "t1Image", cbica::Parameter::BOOLEAN, "0 or 1", "T1 Image being passed or not", "Default: 1");
  parser.exampleUsage("WhiteStripe -i in.nii.gz -o <output dir>");

  bool t1Image = true, skullStrippedImage = true;
  int zSliceStart = -1, zSliceEnd = -1, tissuesMax = 5, histSize = 2000;
  float radius = 0.05, maxSmooth = 10.0, deltaSmooth = 0.5;
  std::string inputImages, outputDir;

  parser.getParameterValue("i", inputImages);
  parser.getParameterValue("o", outputDir);

  if (parser.isPresent("t1"))
  {
    int temp;
    parser.getParameterValue("t1", temp);
    if (temp != 1)
    {
      t1Image = false;
    }
  }

  if (parser.isPresent("sk"))
  {
    int temp;
    parser.getParameterValue("sk", temp);
    if (temp != 1)
    {
      if (!parser.isPresent("z"))
      {
        std::cerr << "Skull stripped image disabled; axial slicing information is needed  using parameter 'zSliceRange'.\n";
        return EXIT_FAILURE;
      }
      else
      {
        std::string zSlice;
        parser.getParameterValue("z", zSlice);
        auto parts = cbica::stringSplit(zSlice, "-");
        if (parts.size() != 2)
        {
          std::cerr << "Number of zSlices not passed correctly. Needs to be in format 'START-END'.\n";
          return EXIT_FAILURE;
        }
        else
        {
          zSliceStart = std::atoi(parts[0].c_str());
          zSliceEnd = std::atoi(parts[1].c_str());
        }
      }
    }
  }

  if (parser.isPresent("r"))
  {
    parser.getParameterValue("r", radius);
  }

  if (parser.isPresent("t"))
  {
    parser.getParameterValue("t", tissuesMax);
  }

  if (parser.isPresent("m"))
  {
    parser.getParameterValue("m", maxSmooth);
  }

  if (parser.isPresent("d"))
  {
    parser.getParameterValue("d", deltaSmooth);
  }

  if (parser.isPresent("b"))
  {
    parser.getParameterValue("b", histSize);
  }

  auto allImages = cbica::stringSplit(inputImages, ",");

  auto originalInputImages = allImages.size();

  // removing duplicate files 
  std::sort(allImages.begin(), allImages.end());
  allImages.erase(std::unique(allImages.begin(), allImages.end()), allImages.end());

  if (originalInputImages != allImages.size())
  {
    std::cout << "Duplicate Images found and removed from processing.\n";
  }

  for (size_t i = 0; i < allImages.size(); i++)
  {
    auto currentImagePath = cbica::normPath(allImages[i]);
    if (!cbica::fileExists(currentImagePath))
    {
      std::cerr << "Image Path provided '" << currentImagePath << "' does not exist; continuing with rest.\n";
    }
    else
    {
      auto currentImage = cbica::ReadImage< ImageTypeFloat3D >(currentImagePath);

      WhiteStripe normalizer;
      normalizer.setParams(radius, zSliceStart, zSliceEnd, tissuesMax, maxSmooth, deltaSmooth, histSize, t1Image);

      ImageTypeFloat3D::Pointer mask;
      auto normImage = normalizer.process(currentImage, mask);

      if (normImage.IsNotNull())
      {
        std::string path, base, ext;
        cbica::splitFileName(currentImagePath, path, base, ext);
        cbica::WriteImage< ImageTypeFloat3D >(normImage, outputDir + "/" + base + "_wsNormalized" + ext);
      }
    }
  }

  return EXIT_SUCCESS;
}