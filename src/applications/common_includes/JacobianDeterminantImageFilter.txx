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
#ifndef __JacobianDeterminantImageFilter_txx_
#define __JacobianDeterminantImageFilter_txx_

#include "itkImageRegionIterator.h"

template <unsigned int VDim>
class DeformedCubeVolumeFunctor
{
public:
  static double CubeVolume(double **Y);
};

template <>
class DeformedCubeVolumeFunctor<3>
{
public:
  static double CubeVolume(double Y[8][3])
    {
    // Pointers to the vertices of the cube
    double *Y000 = Y[0], *Y001 = Y[1], *Y010 = Y[2], *Y011 = Y[3], 
           *Y100 = Y[4], *Y101 = Y[5], *Y110 = Y[6], *Y111 = Y[7];

    // Faces enumerated in the right winding
    return 
      FaceVolume(Y000, Y001, Y101, Y100) +
      FaceVolume(Y100, Y101, Y111, Y110) +
      FaceVolume(Y001, Y011, Y111, Y101) + 
      FaceVolume(Y000, Y100, Y110, Y010) +
      FaceVolume(Y000, Y010, Y011, Y001) +
      FaceVolume(Y010, Y110, Y111, Y011); 
/*
    return 
      FaceVolume(Y000,Y001,Y011,Y010) + 
      FaceVolume(Y100,Y110,Y111,Y101) + 
      FaceVolume(Y000,Y100,Y101,Y001) + 
      FaceVolume(Y010,Y011,Y111,Y110) + 
      FaceVolume(Y000,Y010,Y110,Y100) + 
      FaceVolume(Y001,Y101,Y111,Y011); */
    }

protected:
  static double FaceVolume(double *A, double *B, double *C, double *D)
    {
    // (A - C) * (B x D)
    double t1 = 
      (A[0] - C[0]) * (B[1] * D[2] - B[2] * D[1]) + 
      (A[1] - C[1]) * (B[2] * D[0] - B[0] * D[2]) + 
      (A[2] - C[2]) * (B[0] * D[1] - B[1] * D[0]);

    // (B - D) * (A x C)
    double t2 = 
      (B[0] - D[0]) * (A[1] * C[2] - A[2] * C[1]) + 
      (B[1] - D[1]) * (A[2] * C[0] - A[0] * C[2]) + 
      (B[2] - D[2]) * (A[0] * C[1] - A[1] * C[0]);

    return (t2 - t1) / 12.0;
    }
};

template <>
class DeformedCubeVolumeFunctor<2>
{
public:

  static double CubeVolume(double Y[4][2])
    {
    // Get the four corners
    double *Y00 = Y[0], *Y01 = Y[1], *Y10 = Y[2], *Y11 = Y[3];

    // Compute the area of the polygon
    return 
      EdgeArea(Y00, Y01) + 
      EdgeArea(Y01, Y11) + 
      EdgeArea(Y11, Y10) + 
      EdgeArea(Y10, Y00);
    }

protected:
  static double EdgeArea(double *A, double *B)
    {
    return A[0] * B[1] - A[1] * B[0];
    }
};


template <>
class DeformedCubeVolumeFunctor<4>
{
public:
  static double CubeVolume(double Y[16][4]) { return 0; }
};


template <class TInputImage, class TOutputImage>
void
JacobianDeterminantImageFilter<TInputImage,TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread,
                       itk::ThreadIdType threadId)
{
  // Create an interpolator for the warp field
  typedef FastLinearInterpolator<InputImageType, double, ImageDimension> FastInterpolator;
  InputImageType *warp = const_cast<InputImageType *>(this->GetInput());
  FastInterpolator it(warp);

  // Array of voxel coordinates
  const unsigned int nCorners = (1 << ImageDimension);
  double Y[nCorners][ImageDimension];

  // Index array where we are sampling
  double cix[ImageDimension], cix_ctr[ImageDimension];

  // Array of offsets for the voxel corners
  double off[nCorners][ImageDimension];
  for(int j = 0; j < nCorners; j++)
    for(int a = 0; a < ImageDimension; a++)
      off[j][a] = (j & (1 << a)) ? -0.5 : 0.5;

  // Create a canonical voxel in physical space
  double CV[nCorners][ImageDimension];
  for(int j = 0; j < nCorners; j++)
    {
    // Create an itk::ContinuousIndex for mapping
    itk::ContinuousIndex<double, ImageDimension> c_index;
    for(int a = 0; a < ImageDimension; a++)
      c_index[a] = off[j][a];

    // Map to an itk::Point
    itk::Point<double, ImageDimension> pt;
    warp->TransformContinuousIndexToPhysicalPoint(c_index, pt);
    std::cout << "CV[" << j << "]=" << std::endl; 
    for(int a = 0; a < ImageDimension; a++)
      {
      CV[j][a] = pt[a];
      std::cout << CV[j][a] << " ";
      }
    std::cout << std::endl;
    }

  // Get the volume of the canonical voxel
  double volRef = DeformedCubeVolumeFunctor<ImageDimension>::CubeVolume(CV);

  std::cout << "reference voxel volume " << volRef << std::endl;

  // Go through the voxels
  typedef itk::ImageRegionConstIteratorWithIndex<InputImageType> InputIter;
  typedef itk::ImageRegionIterator<OutputImageType> OutputIter;
  InputIter it_input(warp, outputRegionForThread);
  OutputIter it_output(this->GetOutput(), outputRegionForThread);
  typename FastInterpolator::OutputComponentType phi_x;
  for (; !it_output.IsAtEnd(); ++it_output, ++it_input)
    {
    // Get the continuos index of the voxel center
    for(int a = 0; a < ImageDimension; a++)
      cix_ctr[a] = it_input.GetIndex()[a];

    // Sample the warp at each of the voxel corners
    for(int j = 0; j < nCorners; j++)
      {
      // Set the sampling coordinate
      for(int a = 0; a < ImageDimension; a++)
        cix[a] = cix_ctr[a] + off[j][a];

      // Sample the displacement at the coordinate
      typename FastInterpolator::InOut rc = it.Interpolate(cix, &phi_x); 

      // Add the displacement to the canonical voxel
      for(int a = 0; a < ImageDimension; a++)
        Y[j][a] = CV[j][a] + phi_x[a];
      }

    // Compute the volume of the deformed voxel
    double volVoxel = DeformedCubeVolumeFunctor<ImageDimension>::CubeVolume(Y);

    // Compute the jacobian determinant
    it_output.Set(volVoxel / volRef);
    }
}

#endif
