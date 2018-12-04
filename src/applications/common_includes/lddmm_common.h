/*=========================================================================

  Program:   ALFABIS fast medical image registration programs
  Language:  C++
  Website:   github.com/pyushkevich/greedy
  Copyright (c) Paul Yushkevich, University of Pennsylvania. All rights reserved.

  This program is part of ALFABIS: Adaptive Large-Scale Framework for
  Automatic Biomedical Image Segmentation.

  ALFABIS development is funded by the NIH grant R01 EB017255.

  ALFABIS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ALFABIS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ALFABIS.  If not, see <http://www.gnu.org/licenses/>.

=========================================================================*/
#ifndef _LDDMM_COMMON_H_
#define _LDDMM_COMMON_H_

typedef unsigned int uint;

typedef double myreal;

// A macro for defining named inputs and outputs to ITK filters
#define itkNamedInputMacro(name, type, key) \
  virtual void Set##name (type *_arg) \
    { \
    itk::ProcessObject::SetInput(key, _arg); \
    } \
  \
  virtual type * Get##name() \
    { \
    return dynamic_cast<type *>(itk::ProcessObject::GetInput(key)); \
    }

#define itkNamedInputGetMacro(name, type, key) \
  virtual type * Get##name() \
    { \
    return dynamic_cast<type *>(itk::ProcessObject::GetInput(key)); \
    }

// A macro for defining named inputs and outputs to ITK filters
#define itkNamedOutputMacro(name, type, key) \
  virtual void Set##name (type *_arg) \
    { \
    itk::ProcessObject::SetOutput(key, _arg); \
    } \
  \
  virtual type * Get##name() \
    { \
    return dynamic_cast<type *>(itk::ProcessObject::GetOutput(key)); \
    }

#define itkNamedOutputGetMacro(name, type, key) \
  virtual type * Get##name() \
    { \
    return dynamic_cast<type *>(itk::ProcessObject::GetOutput(key)); \
    }





#endif
