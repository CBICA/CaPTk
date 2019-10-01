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

/// original main function
int main(int argc, char** argv)
{
  std::string argv_complete; // this is where we will be putting stuff for Preprocessing -reg

  // helper variables
  std::string fixedImage;
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
  parser.addOptionalParameter("ri", "radius", cbica::Parameter::STRING, "none", "Patch radius for metrics", "Eg: 2x2x2");

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
    // dicom image detected
    cbica::WriteImage< ImageTypeFloat3D >(
      cbica::ReadImage< ImageTypeFloat3D >(fixedImage),
      tempFolderLocation + "/tempDicomConverted_fixed.nii.gz"
      );
    fixedImage = tempFolderLocation + "/tempDicomConverted_fixed.nii.gz";
  }
  for (int i = 0; i < inputImageFiles.size(); i++)
  {

    if (parser.isPresent("reg")) {
      std::cout << "--> Analyzing parameters:" << std::endl;

      InterpSpec interp_current;
      bool affineMode = false;
      float current_weight = 1.0;

      ImagePairSpec ip;
      ip.weight = current_weight;
      ip.fixed = fixedImage;
      if (cbica::IsDicom(inputImageFiles[i]))
      {
        // dicom image detected
        cbica::WriteImage< ImageTypeFloat3D >(
          cbica::ReadImage< ImageTypeFloat3D >(inputImageFiles[i]),
          tempFolderLocation + "/tempDicomConverted_moving.nii.gz"
          );
        ip.moving = tempFolderLocation + "/tempDicomConverted_moving.nii.gz";
      }
      else
      {
        ip.moving = inputImageFiles[i];
      }
      param.inputs.push_back(ip);

      if (parser.isPresent("a")) {
        std::cout << "--> Registration Mode: Affine " << std::endl;

        param.mode = GreedyParameters::AFFINE;
        param.affine_dof = GreedyParameters::DOF_AFFINE;
      }

      if (parser.isPresent("r")) {
        std::cout << "--> Registration Mode: Rigid " << std::endl;

        param.mode = GreedyParameters::AFFINE;
        param.affine_dof = GreedyParameters::DOF_RIGID;
      }

      if (parser.isPresent("n")) {
        parser.getParameterValue("n", g_iterations);
        std::cout << "--> Number of iterations: " << g_iterations << std::endl;
        param.iter_per_level = cl.read_int_vector(g_iterations);
      }
      else {
        std::cout << "--> Using default iterations: 100x100" << std::endl;
      }

      if (parser.isPresent("m"))
      {
        parser.getParameterValue("m", metrics);

        if (metrics == "NCC" || metrics == "ncc")
        {
          metrics = "NCC";
          param.metric = GreedyParameters::NCC;
        }
        else if (metrics == "MI" || metrics == "mi")
        {
          metrics = "MI";
          param.metric = GreedyParameters::MI;
        }
        else if (metrics == "NMI" || metrics == "nmi")
        {
          metrics = "NMI";
          param.metric = GreedyParameters::NMI;
        }
        else if (metrics == "ssd" || metrics == "SSD")
        {
          metrics = "SSD";
          param.metric = GreedyParameters::SSD;
        }

        std::cout << "--> Metric used: " << metrics << std::endl;

        if (parser.isPresent("ri")) {
          std::string val;
          parser.getParameterValue("ri", val);
          std::cout << "--> Patch radius used for metrics: " << val << std::endl;
          param.metric_radius = cl.read_int_vector(val);
        }
      }

      std::cout << "--> Transformation matrix output file name: " << matrixImageFiles[i] << std::endl;
      param.output = matrixImageFiles[i];

      std::cout << "--> Looking for fixed image: " << std::endl;

      if (!cbica::fileExists(fixedImage))
      {
        std::cerr << "--> Fixed image file not found :'" << fixedImage << "'\n";
        return EXIT_FAILURE;
      }
      std::cout << "--> Fixed image found: " + fixedImage << std::endl;

      std::cout << "--> Looking for moving image: " << std::endl;

      if (!cbica::fileExists(inputImageFiles[i]))
      {
        std::cerr << "--> Moving image file not found :'" << inputImageFiles[i] << "'\n";
        return EXIT_FAILURE;
      }
      std::cout << "--> Moving image found: " + inputImageFiles[i] << std::endl;

      if (param.threads > 0)
      {
        std::cout << "--> Limiting the number of threads to " << param.threads << std::endl;
        itk::MultiThreader::SetGlobalMaximumNumberOfThreads(param.threads);
      }
      else
      {
        std::cout << "--> Executing with the default number of threads: " << itk::MultiThreader::GetGlobalDefaultNumberOfThreads() << std::endl;

      }

      // Some parameters may be specified as either vector or scalar, and need to be verified
      if (param.epsilon_per_level.size() != param.iter_per_level.size())
      {
        if (param.epsilon_per_level.size() == 1)
        {
          param.epsilon_per_level =
            std::vector<double>(param.iter_per_level.size(), param.epsilon_per_level.back());
        }
        else
        {
          throw GreedyException("--> Mismatch in size of vectors supplied with -n and -e options");
        }
      }

      auto fixedImageInfo = cbica::ImageInfo(fixedImage);
      auto movingImageInfo = cbica::ImageInfo(inputImageFiles[i]);

      if (fixedImageInfo.GetImageDimensions() != movingImageInfo.GetImageDimensions())
      {
        std::cerr << "--> Image dimensions do not match." << std::endl;
        return EXIT_FAILURE;
      }
      else
      {
        param.dim = fixedImageInfo.GetImageDimensions();
      }

      if (param.flag_float_math)
      {
        switch (fixedImageInfo.GetImageDimensions())
        {
        case 2: GreedyRunner<2, float>::Run(param); break;
        case 3: GreedyRunner<3, float>::Run(param); break;
        case 4: GreedyRunner<4, float>::Run(param); break;
        default: throw GreedyException("--> Wrong number of dimensions requested: %d", param.dim);
        }
      }
      else
      {
        switch (fixedImageInfo.GetImageDimensions())
        {
        case 2: GreedyRunner<2, double>::Run(param); break;
        case 3: GreedyRunner<3, double>::Run(param); break;
        case 4: GreedyRunner<4, double>::Run(param); break;
        default: throw GreedyException("--> Wrong number of dimensions requested: %d", param.dim);
        }
      }

      std::cout << "--> Finished registration.\n";

      //GreedyParameters param;
      //GreedyParameters::SetToDefaults(param);

      if (parser.isPresent("trf"))
      {
        std::cout << "--> Starting transformation.\n";


        if (parser.isPresent("f"))
        {
          parser.getParameterValue("f", fixedImage);
        }

        param.reslice_param.ref_image = fixedImage;
        TransformSpec spec;
        ResliceSpec reslice;

        InterpSpec interp_current;
        reslice.interp = interp_current;
        reslice.moving = inputImageFiles[i];
        reslice.output = outputImageFiles[i];

        param.reslice_param.images.push_back(reslice);

        param.mode = GreedyParameters::RESLICE;

        std::cout << "--> Looking for transformation matrix: " + matrixImageFiles[i] << std::endl;

        if (!cbica::fileExists(matrixImageFiles[i])) {
          std::cerr << "--> Transformation matrix not found at: " << matrixImageFiles[i] << "'\n";
          return EXIT_FAILURE;
        }
        std::cout << "--> Transformation Matrix found: " + matrixImageFiles[i] << std::endl;
        spec = cl.read_transform_spec(matrixImageFiles[i]);

        param.reslice_param.transforms.push_back(spec);


        /*if (param.threads > 0)
        {
            std::cout << "--> Limiting the number of threads to " << param.threads << std::endl;
            itk::MultiThreader::SetGlobalMaximumNumberOfThreads(param.threads);
        }
        else
        {
            std::cout << "--> Executing with the default number of threads: " << itk::MultiThreader::GetGlobalDefaultNumberOfThreads() << std::endl;
        }*/
        // Some parameters may be specified as either vector or scalar, and need to be verified
        if (param.epsilon_per_level.size() != param.iter_per_level.size())
        {
          if (param.epsilon_per_level.size() == 1)
          {
            param.epsilon_per_level =
              std::vector<double>(param.iter_per_level.size(), param.epsilon_per_level.back());
          }
          else
          {
            throw GreedyException("--> Mismatch in size of vectors supplied with -n and -e options");
          }
        }

        auto fixedImageInfo = cbica::ImageInfo(fixedImage);
        auto movingImageInfo = cbica::ImageInfo(inputImageFiles[i]);

        if (fixedImageInfo.GetImageDimensions() != movingImageInfo.GetImageDimensions())
        {
          std::cerr << "--> Image dimensions do not match." << std::endl;
          return EXIT_FAILURE;
        }
        else
        {
          param.dim = fixedImageInfo.GetImageDimensions();
        }

        std::cout << "--> Applied transformation to moving image: " << outputImageFiles[i] << std::endl;

        if (param.flag_float_math)
        {
          switch (fixedImageInfo.GetImageDimensions())
          {
          case 2: GreedyRunner<2, float>::Run(param); continue;
          case 3: GreedyRunner<3, float>::Run(param); continue;
          case 4: GreedyRunner<4, float>::Run(param); continue;
          default: throw GreedyException("--> Wrong number of dimensions requested: %d", param.dim);
          }
        }
        else
        {
          switch (fixedImageInfo.GetImageDimensions())
          {
          case 2: GreedyRunner<2, double>::Run(param); continue;
          case 3: GreedyRunner<3, double>::Run(param); continue;
          case 4: GreedyRunner<4, double>::Run(param); continue;
          default: throw GreedyException("--> Wrong number of dimensions requested: %d", param.dim);
          }
        }

        std::cout << "--> Transformation complete " << std::endl;
      }

      if (argc > 1)
      {
        for (size_t i = 1; i < argc; i++)
        {
          argv_complete += " " + std::string(argv[i]);
        }
      }
      auto commandToRun = cbica::getExecutablePath() + "/Preprocessing" +
#if WIN32
        ".exe" +
#endif
        argv_complete;

      return std::system(commandToRun.c_str());
    }
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

  std::string argv_complete;
  if (argc > 1)
  {
    for (size_t i = 1; i < argc; i++)
    {
      argv_complete += " " + std::string(argv[i]);
    }
  }
  auto commandToRun = cbica::getExecutablePath() + "/Preprocessing" + 
#if WIN32
    ".exe" +
#endif
    argv_complete;

  return std::system(commandToRun.c_str());
}
