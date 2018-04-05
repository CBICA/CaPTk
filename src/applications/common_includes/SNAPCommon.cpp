/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: SNAPCommon.cxx.in,v $
  Language:  C++
  Date:      $Date: 2009/02/03 20:28:14 $
  Version:   $Revision: 1.4 $
  Copyright (c) 2003 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
#include "SNAPCommon.h"

// Non-cvs version
const char SNAPSoftVersion[] = "ITK-SNAP @SNAP_VERSION_FULL@";

// Just the number part of the SNAP version
unsigned int SNAPVersionMajor = 3;
unsigned int SNAPVersionMinor = 6;
unsigned int SNAPVersionPatch = 0;
const char SNAPVersionQualifier[] = "";

// Hardware architecture for this build
const char SNAPArch[] = "NSIS";

// A string that appears in the user interface
const char SNAPUISoftVersion[] = 
  "Version 3.6.0";

// Release date of the current version
const char SNAPCurrentVersionReleaseDate[] = "";

// String describing the current build environment
const char SNAPBuildInfo[] = "";

// Release date of the latest version whose user preference files are
// incompatible with the current version and will be erased
const char SNAPLastIncompatibleReleaseDate[] = "";

// Build date - shown to help debugging nightly builds
const char SNAPBuildDate[] = "";

// GIT signature
const char SNAPGitSignature[] = "";
