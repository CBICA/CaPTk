/**
\file  GreedyRegistration.cxx

\brief GreedyRegistration comand line interface.

Dependecies: ITK (module_review, module_skullstrip enabled)

https://www.cbica.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.cbica.upenn.edu/sbia/software/license.html

*/

#include <iostream>
#include <cstdlib>
#include <string>
#include <algorithm>

#include "cbicaCmdParser.h"
#include "cbicaUtilities.h"
#include "cbicaITKSafeImageIO.h"

/// original main function
int main(int argc, char** argv)
{
  std::string argv_common; // this is where we will be putting stuff for Preprocessing -reg

  // helper variables
  std::string fixedImage, noIterations = "100,50,5", metrics = "NMI", nccRadius = "2x2x2";
  std::vector< std::string > inputImageFiles, outputImageFiles, matrixImageFiles;

  cbica::CmdParser parser(argc, argv, "GreedyRegistration");

  parser.addRequiredParameter("i", "movingImage", cbica::Parameter::FILE, "NIfTI", "Input Image for processing", "Becomes moving image in Registration mode");
  parser.addRequiredParameter("f", "fixedImage", cbica::Parameter::STRING, "NIfTI", "Fixed Image for registration");
  parser.addRequiredParameter("o", "outputImage", cbica::Parameter::FILE, "NIfTI", "Output Image for processing");
  parser.addRequiredParameter("t", "matrix", cbica::Parameter::STRING, "N.A", "Registration Matrix");

  parser.addOptionalParameter("reg", "registration", cbica::Parameter::NONE, "N.A", "Switch to registration mode");
  parser.addOptionalParameter("trf", "transformation", cbica::Parameter::NONE, "N.A", "Switch to transformation mode");

  parser.addOptionalParameter("a", "affine", cbica::Parameter::NONE, "N.A", "Affine Registration(Default)");
  parser.addOptionalParameter("r", "rigid", cbica::Parameter::NONE, "N.A", "Rigid Registration");
  parser.addOptionalParameter("m", "metrics", cbica::Parameter::STRING, "none", "MI: mutual information", "NMI(Default): normalized mutual information", "NCC -r 2x2x2: normalized cross-correlation");
  parser.addOptionalParameter("ri", "radius", cbica::Parameter::STRING, "none", "Patch radius for metrics", "Eg: 2x2x2 (Default)");

  parser.addOptionalParameter("n", "greedyIterations", cbica::Parameter::STRING, "none", "Number of iterations per level of multi-res (Default: 100x50x5)", "Corresponds to low level, Mid Level and High Level resolution", "Pattern: NxNxN");
  parser.addOptionalParameter("th", "threads", cbica::Parameter::INTEGER, "none", "Number of threads for algorithm", "IGNORED");
  //parser.exampleUsage("-reg -trf -i moving.nii.gz -f fixed.nii.gz -o output.nii.gz -t matrix.mat -a -m MI -n 100x50x5 -th 4");

  parser.addApplicationDescription("This does affine registration based on Greedy");
  parser.addExampleUsage("-reg -trf -i moving.nii.gz -f fixed.nii.gz -o output.nii.gz -t matrix.mat -a -m MI -n 100x50x5",
    "This registers the moving image 'moving.nii.gz' with fixed image 'fixed.nii.gz.' with output at 'output.nii.gz'");


  if (parser.isPresent("i"))
  {
    std::string inputFileString;
    parser.getParameterValue("i", inputFileString);
    //parser.getParameterValue("i", inputImageFile);
    inputImageFiles = cbica::stringSplit(inputFileString, ",");
    for (size_t i = 0; i < inputImageFiles.size(); i++)
    {
      if (cbica::isDir(inputImageFiles[i]))
      {
        std::cerr << "Directory detected for '" << inputImageFiles[i] << "', please check to file.\n";
        return EXIT_FAILURE;
      }
    }
  }

  if (parser.isPresent("f"))
  {
    parser.getParameterValue("f", fixedImage);
  }

  if (parser.isPresent("o"))
  {
    std::string outputFileString;
    parser.getParameterValue("o", outputFileString);
    //parser.getParameterValue("o", outputImageFile);
    outputImageFiles = cbica::stringSplit(outputFileString, ",");
  }

  if (parser.isPresent("t"))
  {
    std::string matrixFileString;
    parser.getParameterValue("t", matrixFileString);
    matrixImageFiles = cbica::stringSplit(matrixFileString, ",");
  }
  
  // sanity check
  if (inputImageFiles.size() != outputImageFiles.size()
    || inputImageFiles.size() != matrixImageFiles.size() || outputImageFiles.size() != matrixImageFiles.size())
  {
    std::cerr << "ERROR: Number of moving images does not match with output images and matrices.\n";
    return EXIT_FAILURE;
  }

  std::string tempFolderLocation = cbica::normPath(cbica::getUserHomeDirectory() + "/.CaPTk");
  cbica::createDir(tempFolderLocation);
  if (cbica::IsDicom(fixedImage))
  {
    // dicom image detected and is assumed to be 3D
    cbica::WriteImage< ImageTypeFloat3D >(
      cbica::ReadImage< ImageTypeFloat3D >(fixedImage),
      tempFolderLocation + "/tempDicomConverted_fixed.nii.gz"
      );
    fixedImage = tempFolderLocation + "/tempDicomConverted_fixed.nii.gz";
  }

  if (parser.isPresent("n"))
  {
    parser.getParameterValue("n", noIterations);
    auto tempIters = cbica::stringSplit(noIterations, "x");
    noIterations = tempIters[0];
    for (size_t n = 1; n < tempIters.size(); n++)
    {
      noIterations += "," + tempIters[n];
    }
  }
  argv_common += " -rNI " + noIterations; // if "n" is not present, we use default

  // fixed image
  parser.getParameterValue("f", fixedImage);
  if (!cbica::fileExists(fixedImage))
  {
    std::cerr << "Couldn't find fixed image: " << fixedImage << ".\n";
    return EXIT_FAILURE;
  }
  argv_common += " -rFI " + fixedImage;
  auto fixedImageInfo = cbica::ImageInfo(fixedImage);

  // metrics
  if (parser.isPresent("m"))
  {
    parser.getParameterValue("m", metrics);

    if (metrics == "NCC" || metrics == "ncc")
    {
      metrics = "NCC";
      if (parser.isPresent("ri"))
      {
        parser.getParameterValue("ri", nccRadius);
      }
      metrics += "-" + nccRadius;
    }
    else if (metrics == "MI" || metrics == "mi")
    {
      metrics = "MI";
    }
    else if (metrics == "NMI" || metrics == "nmi")
    {
      metrics = "NMI";
    }
    else if (metrics == "ssd" || metrics == "SSD")
    {
      metrics = "SSD";
    }
  }
  argv_common += metrics;

  argv_common += " -reg Affine"; // TBD: switch for affine or deformable to be added later
  for (int i = 0; i < inputImageFiles.size(); i++)
  {
    std::string argv_currentMovingImage;

    // moving image
    argv_currentMovingImage += " -i ";
    if (cbica::IsDicom(inputImageFiles[i]))
    {
      // dicom image detected
      cbica::WriteImage< ImageTypeFloat3D >(
        cbica::ReadImage< ImageTypeFloat3D >(inputImageFiles[i]),
        tempFolderLocation + "/tempDicomConverted_moving.nii.gz"
        );
      argv_currentMovingImage += tempFolderLocation + "/tempDicomConverted_moving.nii.gz";
    }
    else
    {
      if (!cbica::fileExists(inputImageFiles[i]))
      {
        std::cerr << "Couldn't find moving image: " << inputImageFiles[i] << "\n";
        return EXIT_FAILURE;
      }
      argv_currentMovingImage += inputImageFiles[i];
    }

    argv_currentMovingImage += " -rIA " + matrixImageFiles[i];

    auto movingImageInfo = cbica::ImageInfo(inputImageFiles[i]);

    if (fixedImageInfo.GetImageDimensions() != movingImageInfo.GetImageDimensions())
    {
      std::cerr << "Image dimensions of fixed image '" << fixedImage << "' and moving image '" <<
        inputImageFiles[i] << "'do not match.\n";
      return EXIT_FAILURE;
    }

    // everything is okay, start registration

    auto commandToRun = cbica::getExecutablePath() + "/Preprocessing" +
#if WIN32
      ".exe" +
#endif
      argv_common + argv_currentMovingImage;

    return std::system(commandToRun.c_str());
  }

}
///


int main(int argc, char** argv)
{
  cbica::CmdParser parser(argc, argv, "GreedyRegistration");

  parser.addRequiredParameter("i", "movingImage", cbica::Parameter::FILE, "NIfTI", "Input Image for processing", "Becomes moving image in Registration mode");
  parser.addRequiredParameter("o", "outputImage", cbica::Parameter::FILE, "NIfTI", "Output Image for processing");
  parser.addOptionalParameter("m", "maskImage", cbica::Parameter::FILE, "NIfTI", "Input Mask for processing");
  
  parser.addOptionalParameter("reg", "registration", cbica::Parameter::STRING, "Affine | Deformable", "The kind of registration to perform", "Can use Mask File");
  parser.addOptionalParameter("rFI", "regFixedImg", cbica::Parameter::FILE, "NIfTI", "The Fixed Image for the registration", "Needed for registration");
  parser.addOptionalParameter("rME", "regMetrics", cbica::Parameter::STRING, "SSD | MI | NMI | NCC-AxBxC", "The kind of metris to use: SSD (Sum of Squared Differences) or MI (Mutual Information) or", "NMI (Normalized Mutual Information) or NCC-AxBxC (Normalized Cross correlation with integer radius for 3D image)");
  parser.addOptionalParameter("rNI", "regNoIters", cbica::Parameter::STRING, "N1,N2,N3", "The umber of iterations per level of multi-res");
  parser.addOptionalParameter("rIS", "regInterSave", cbica::Parameter::BOOLEAN, "0 or 1", "Whether the intermediate files are to be saved or not");
  parser.addOptionalParameter("rSg", "regSegMoving", cbica::Parameter::BOOLEAN, "0 or 1", "Whether the Moving Image is a segmentation file", "If 1, the 'Nearest Label' Interpolation is applied");
  parser.addOptionalParameter("rIA", "regInterAffn", cbica::Parameter::FILE, "mat", "The path to the affine transformation to apply to moving image", "If this is present, the Affine registration step will be skipped");
  parser.addOptionalParameter("rID", "regInterDefm", cbica::Parameter::FILE, "NIfTI", "The path to the deformable transformation to apply to moving image", "If this is present, the Deformable registration step will be skipped");

  parser.addApplicationDescription("This does affine registration based on Greedy and internally calls Preprocessing; run 'Preprocessing -h' for details");
  parser.addExampleUsage("-i moving.nii.gz -rFI fixed.nii.gz -o output.nii.gz -reg Affine -rME NCC-2x2x2 -rIS 1 -rNI 100,50,5",
    "This registers the moving image affinely 'moving.nii.gz' with fixed image 'fixed.nii.gz.' with output at 'output.nii.gz' using NCC metric with kernel size 2x2x2");

  std::string argv_common;
  if (argc > 1)
  {
    for (size_t i = 1; i < argc; i++)
    {
      argv_common += " " + std::string(argv[i]);
    }
  }
  auto commandToRun = cbica::getExecutablePath() + "/Preprocessing" + 
#if WIN32
    ".exe" +
#endif
    argv_common;

  return std::system(commandToRun.c_str());
}
