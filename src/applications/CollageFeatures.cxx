#include "cbicaCmdParser.h"
#include "cbicaLogging.h"


// This application simply wraps the collageradiomics cli (included in bin as part of FeatureExtraction being built) and shows warnings.
int main(int argc, char **argv)
{
  cbica::CmdParser parser = cbica::CmdParser(argc, argv, "CollageFeatures");
  
  // Right now we add these parameters explicitly, perhaps we can just use the wrapped executable's text.
  parser.addRequiredParameter("i", "input", cbica::Parameter::STRING, "", "Path to an input image from which features will be extracted.");
  parser.addRequiredParameter("m", "mask", cbica::Parameter::STRING, "", "Path to a mask that will be considered as binary. The highest pixel value will be considered as information and all other values will be considered outside the mask.");
  parser.addRequiredParameter("o", "outputfile", cbica::Parameter::STRING, "", "Path to the output CSV file.");
  parser.addOptionalParameter("v", "verbose", cbica::Parameter::BOOLEAN, "", "Provides additional debug output.");
  parser.addOptionalParameter("d", "dimensions", cbica::Parameter::INTEGER, "2-3", "Optional number of dimensions upon which to run collage. Supported values are 2 and 3. If left out, we will default to the dimensionality of the image itself, which may not reflect expected behavior if the image has an alpha channel.");
  parser.addOptionalParameter("s", "svdradius", cbica::Parameter::INTEGER, "", "SVD radius is used for the dominant angle calculation pixel radius. DEFAULTS to 5 and is suggested to remain at the default.");
  parser.addOptionalParameter("h", "haralickwindow", cbica::Parameter::INTEGER, "", "Number of pixels around each pixel used to calculate the haralick texture. DEFAULTS to svdradius * 2 - 1.");
  parser.addOptionalParameter("b", "binsize", cbica::Parameter::INTEGER, "", "Number of bins to use while calculating the grey level cooccurence matrix. DEFAULTS to 64.");
  
  parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");
  parser.addExampleUsage("-i C:/inputImage.nii.gz -m C:/inputMask.nii.gz -o C:/output.csv", "Generates COLLAGE features with default parameters");
  parser.addApplicationDescription("CoLlAGe captures subtle anisotropic differences in disease pathologies by measuring entropy of co-occurrences of voxel-level gradient orientations onn imaging computed within a local neighborhood. The CoLlAGe features were developed by BrIC Lab and the RadxTools team (https://doi.org/10.1038/srep37241, https://github.com/radxtools/collageradiomics). Feature requests, bug reports, and any other CoLlAGe related considerations, should be directed to https://github.com/radxtools/collageradiomics/issues .");
  
  std::string outputFile;
  parser.getParameterValue("o", outputFile);
  
  std::string argv_complete;
  for (int i = 1; i < argc; i++) // 1 to ignore the current executable name
  {
    argv_complete += " \"" + std::string(argv[i]) + "\""; // double quote for command tokenization
  }
  
  auto collage_cli = cbica::getExecutablePath() + "/collageradiomics"
#ifdef _WIN32
      + ".exe"
#else

#endif
      ;

  if (!cbica::isFile(collage_cli))
  {
      // If this happens, this is on our packaging and not on the RadxTools folks.
      std::cerr << "Could not find the Collage CLI, so cannot extract collage features. Please report this bug to CBICA Software.\n";
      return EXIT_FAILURE;
  }


  std::string commandToRun = collage_cli + argv_complete;
  if(std::system(commandToRun.c_str()) > 0) 
  {
    std::cerr << "Errors occurred while running the collageradiomics application. See logs for more details." << std::endl;
	return EXIT_FAILURE;
  }
  else 
  {
     std::cout << "Finished successfully.\n";
  }
  return EXIT_SUCCESS;
}
