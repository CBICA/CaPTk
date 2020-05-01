/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkAnatomicalOrientation.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkAnatomicalOrientation
 * @brief   Define anatomical orientation axis, coordinates and transforms
 *
 * This class defines the different anatomical axes that are used in
 * orientated images, typically in the medical field: L (left), R (right),
 * P (posterior), A (anterior), S (superior) and I (inferior). This class
 * exposes static methods to handle those axes. An instance of this class
 * represents a coordinate system composed of three of those axes (XYZ),
 * which can be accessed by the X, Y, Z variables, or with the [] operator.
 *
 * In VTK, images are defined in the LPS coordinate system. The main
 * purpose of this class is to be able to convert images read from
 * third-party formats which are not in LPS coordinates, by providing a
 * transform to go from that coordinate system to LPS. That transform can
 * be applied to the data or to the direction matrix of the image.
 *
 * @sa
 * vtkImageReader2
*/

#ifndef vtkAnatomicalOrientation_h
#define vtkAnatomicalOrientation_h

#include "vtkCommonMathModule.h" // For export macro

#include <string>

class VTKCOMMONMATH_EXPORT vtkAnatomicalOrientation
{
public:
  //@{
  /**
   * Anatomical Axis enumeration. The valid ones are stored in ValidAxes.
   */
  enum class Axis { L, R, P, A, S, I, None };
  static const Axis ValidAxes[6];
  //@}

  /**
   * Get the opposite axis of an axis
   */
  static Axis AxisInverse(Axis axis);

  //@{
  /**
   * Get a string or char representation of an axis, and vice-versa
   */
  static char AxisToChar(Axis dir);
  static Axis AxisFromChar(char letter);
  static std::string AxisToString(Axis axis);
  static std::string AxisToLowercaseString(Axis axis);
  static Axis AxisFromString(std::string name);
  //@}

  /**
   * Constructor defaulting to invalid axes (None).
   * X, Y and Z can be inputed.
   */
  vtkAnatomicalOrientation(Axis X = Axis::None,
                           Axis Y = Axis::None,
                           Axis Z = Axis::None);

  /**
   * String parsing constructor.
   * @sa AxisFromString AxisFromChar
   */
  vtkAnatomicalOrientation(std::string str);

  //@{
  /**
   * Static constant instances of the class. LPS is the coordinate
   * system used in VTK. Other are defined for convenience.
   */
  static const vtkAnatomicalOrientation LPS;
  static const vtkAnatomicalOrientation RAS;
  static const vtkAnatomicalOrientation LAS;
  //@}

  /**
   * Returns true if the coordinate system is 3D (no None), and the axes
   * all differ from each other, and each other's opposite.
   */
  bool IsValid() const;

  //@{
  /**
   * Get/Set the string or char representation of the coordinate system
   */
  std::string GetAsString(std::string separator) const;
  void SetForString(std::string str, char delim);
  std::string GetAsAcronym() const;
  void SetForAcronym(std::string acronym);
  //@}

  /**
   * Get the  9 values of the 3x3 matrix to project from this
   * system of coordinate to a given system of coordinates
   */
  void GetTransformTo(vtkAnatomicalOrientation finalOrientation, double transform[9]) const;

  //@{
  /**
   * The world axis values of the coordinate system. They can be accessed
   * and modified directly, or using the [] operator (indexed 0 to 2).
   */
  Axis X;
  Axis Y;
  Axis Z;
  //@}

  Axis operator[](int i) const;
  Axis& operator[](int i);
  bool operator==(const vtkAnatomicalOrientation& rhs) const;
  bool operator!=(const vtkAnatomicalOrientation& rhs) const;
};

#endif
