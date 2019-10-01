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

int main(int argc, char** argv)
{
  cbica::CmdParser parser(argc, argv, "GreedyRegistration");

  parser.addRequiredParameter("i", "movingImage", cbica::Parameter::FILE, "NIfTI", "Input Image for processing", "Becomes moving image in Registration mode");
  parser.addRequiredParameter("o", "outputImage", cbica::Parameter::FILE, "NIfTI", "Output Image for processing");
  parser.addOptionalParameter("m", "maskImage", cbica::Parameter::FILE, "NIfTI", "Input Mask for processing");
  
  parser.addOptionalParameter("reg", "registration", cbica::Parameter::STRING, "Affine | Deformable", "The kind of registration to perform", "Can use Mask File");
  parser.addOptionalParameter("rFI", "regFixedImg", cbica::Parameter::FILE, "NIfTI", "The Fixed Image for the registration", "Needed for registration");
  parser.addOptionalParameter("rME", "regMetrics", cbica::Parameter::STRING, "SSD | MI | NMI | NCC-AxBxC", "The kind of metris to use: SSD (Sum of Squared Differences) or MI (Mutual Information) or", "NMI (Normalized Mutual Information) or NCC-AxBxC (Normalized Cross correlation with integer radius for 3D image)", "Defaults to " + registrationMetrics);
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
  auto commandToRun = cbica::getExecutablePath() + "/Preprocessing" + argv_complete;

  return std::system(commandToRun.c_str());
}
