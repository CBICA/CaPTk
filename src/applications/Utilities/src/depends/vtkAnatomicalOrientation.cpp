/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkAnatomicalOrientation.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkAnatomicalOrientation.h"

#include <iostream>
#include <sstream>
#include <algorithm>

//----------------------------------------------------------------------------
const vtkAnatomicalOrientation::Axis vtkAnatomicalOrientation::ValidAxes[6] =
{
  Axis::L,
  Axis::R,
  Axis::P,
  Axis::A,
  Axis::S,
  Axis::I
};

//----------------------------------------------------------------------------
vtkAnatomicalOrientation::Axis vtkAnatomicalOrientation::AxisInverse(Axis axis)
{
  switch (axis)
  {
    case Axis::L: return Axis::R;
    case Axis::R: return Axis::L;
    case Axis::P: return Axis::A;
    case Axis::A: return Axis::P;
    case Axis::S: return Axis::I;
    case Axis::I: return Axis::S;
    case Axis::None: return Axis::None;
  }
}

//----------------------------------------------------------------------------
char vtkAnatomicalOrientation::AxisToChar(Axis dir)
{
  switch (dir)
  {
    case Axis::L: return 'L';
    case Axis::R: return 'R';
    case Axis::P: return 'P';
    case Axis::A: return 'A';
    case Axis::S: return 'S';
    case Axis::I: return 'I';
    case Axis::None: return 0;
  }
}

//----------------------------------------------------------------------------
vtkAnatomicalOrientation::Axis vtkAnatomicalOrientation::AxisFromChar(char letter)
{
  for (const Axis& axis : ValidAxes)
  {
    if (letter == AxisToChar(axis))
    {
      return axis;
    }
  }
  return Axis::None;
}

//----------------------------------------------------------------------------
std::string vtkAnatomicalOrientation::AxisToString(Axis axis)
{
  switch (axis)
  {
    case Axis::L: return "Left";
    case Axis::R: return "Right";
    case Axis::P: return "Posterior";
    case Axis::A: return "Anterior";
    case Axis::S: return "Superior";
    case Axis::I: return "Inferior";
    case Axis::None: return "";
  }
}

//----------------------------------------------------------------------------
std::string vtkAnatomicalOrientation::AxisToLowercaseString(Axis axis)
{
  std::string name = AxisToString(axis);
  std::transform(name.begin(), name.end(), name.begin(), ::tolower);
  return name;
}

//----------------------------------------------------------------------------
vtkAnatomicalOrientation::Axis vtkAnatomicalOrientation::AxisFromString(std::string name)
{
  for (const Axis& axis : ValidAxes)
  {
    if (name == AxisToString(axis) || name == AxisToLowercaseString(axis))
    {
      return axis;
    }
  }
  return Axis::None;
}

//----------------------------------------------------------------------------
vtkAnatomicalOrientation::vtkAnatomicalOrientation(Axis X, Axis Y, Axis Z)
: X(X), Y(Y), Z(Z)
{
}

//----------------------------------------------------------------------------
vtkAnatomicalOrientation::vtkAnatomicalOrientation(std::string str)
{
  this->SetForAcronym(str);
  if (X == Axis::None || Y == Axis::None)
  {
    this->SetForString(str, '-');
  }
}

//----------------------------------------------------------------------------
const vtkAnatomicalOrientation vtkAnatomicalOrientation::LPS = vtkAnatomicalOrientation("LPS");
const vtkAnatomicalOrientation vtkAnatomicalOrientation::RAS = vtkAnatomicalOrientation("RAS");
const vtkAnatomicalOrientation vtkAnatomicalOrientation::LAS = vtkAnatomicalOrientation("LAS");

//----------------------------------------------------------------------------
bool vtkAnatomicalOrientation::IsValid() const
{
  return
  this->X != Axis::None &&
  this->Y != Axis::None &&
  this->Z != Axis::None &&
  this->Z != this->Y &&
  this->Z != this->X &&
  this->Y != this->X &&
  this->Z != AxisInverse(this->Y) &&
  this->Z != AxisInverse(this->X) &&
  this->Y != AxisInverse(this->X);
}

//----------------------------------------------------------------------------
std::string vtkAnatomicalOrientation::GetAsAcronym() const
{
  std::string acronym = "";
  acronym += AxisToChar(this->X);
  acronym += AxisToChar(this->Y);
  acronym += AxisToChar(this->Z);
  return acronym;
}

//----------------------------------------------------------------------------
void vtkAnatomicalOrientation::SetForAcronym(std::string acronym)
{
  try {
    this->X = AxisFromChar(acronym.at(0));
  } catch (const std::exception&) {
    this->X = Axis::None;
  }
  try {
    this->Y = AxisFromChar(acronym.at(1));
  } catch (const std::exception&) {
    this->Y = Axis::None;
  }
  try {
    this->Z = AxisFromChar(acronym.at(2));
  } catch (const std::exception&) {
    this->Z = Axis::None;
  }
}

//----------------------------------------------------------------------------
std::string vtkAnatomicalOrientation::GetAsString(std::string separator) const
{
  std::string XStr = AxisToLowercaseString(this->X);
  std::string YStr = AxisToLowercaseString(this->Y);
  std::string ZStr = AxisToLowercaseString(this->Z);
  std::string str = XStr;
  if (YStr.length() > 0)
  {
    if (str.length() > 0)
    {
      str += separator;
    }
    str += YStr;
  }
  if (ZStr.length() > 0)
  {
    if (str.length() > 0)
    {
      str += separator;
    }
    str += ZStr;
  }
  return str;
}

//----------------------------------------------------------------------------
void vtkAnatomicalOrientation::SetForString(std::string str, char delim)
{
  std::string XYZStr[3];
  int i = 0;
  std::stringstream ss(str);
  std::string token;
  while (std::getline(ss, token, delim) && i < 3)
  {
    XYZStr[i] = token;
    i++;
  }
  this->X = AxisFromString(XYZStr[0]);
  this->Y = AxisFromString(XYZStr[1]);
  this->Z = AxisFromString(XYZStr[2]);
}

//----------------------------------------------------------------------------
void vtkAnatomicalOrientation::GetTransformTo(vtkAnatomicalOrientation finalOrientation,
                                              double transform[9]) const
{
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if ((*this)[i] == finalOrientation[j])
        transform[i*3+j] = 1;
      else if ((*this)[i] == AxisInverse(finalOrientation[j]))
        transform[i*3+j] = -1;
      else
        transform[i*3+j] = 0;
    }
  }
}

//----------------------------------------------------------------------------
vtkAnatomicalOrientation::Axis vtkAnatomicalOrientation::operator[](int i) const
{
  if (i == 0) return this->X;
  if (i == 1) return this->Y;
  if (i == 2) return this->Z;
  else
  {
    std::cerr << "Out of bounds of the coordinate system" << std::endl;
    exit(0);
  }
}

//----------------------------------------------------------------------------
vtkAnatomicalOrientation::Axis& vtkAnatomicalOrientation::operator[](int i)
{
  if (i == 0) return this->X;
  if (i == 1) return this->Y;
  if (i == 2) return this->Z;
  else
  {
    std::cerr << "Out of bounds of the coordinate system" << std::endl;
    exit(0);
  }
}

//----------------------------------------------------------------------------
bool vtkAnatomicalOrientation::operator==(const vtkAnatomicalOrientation& rhs) const
{
  return this->X == rhs.X && this->Y == rhs.X && this->Z == rhs.Z;
}

//----------------------------------------------------------------------------
bool vtkAnatomicalOrientation::operator!=(const vtkAnatomicalOrientation& rhs) const
{
  return !(*this == rhs);
}
