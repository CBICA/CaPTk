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
#ifndef LINEARTRANSFORMTOWARPFILTER_TXX
#define LINEARTRANSFORMTOWARPFILTER_TXX

#include "LinearTransformToWarpFilter.h"
#include "ImageRegionConstIteratorWithIndexOverride.h"


template <class TInputImage, class TDeformationField, class TTransform>
void
LinearTransformToWarpFilter<TInputImage,TDeformationField,TTransform>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread,
                       itk::ThreadIdType threadId)
{



  // Get a pointer to the output deformation field
  DeformationVectorType *b_phi = this->GetOutput()->GetBufferPointer();

  // Affine transform matrix and vector
  vnl_matrix_fixed<double, ImageDimension, ImageDimension> M =
      this->GetTransform()->GetMatrix().GetVnlMatrix();
  vnl_vector_fixed<double, ImageDimension> off =
      this->GetTransform()->GetOffset().GetVnlVector();

  // Create an iterator over the deformation field
  typedef itk::ImageLinearIteratorWithIndex<DeformationFieldType> IterBase;
  typedef IteratorExtender<IterBase> IterType;

  // Loop over the lines in the image
  for(IterType it(this->GetOutput(), outputRegionForThread); !it.IsAtEnd(); it.NextLine())
    {
    // Get the index at the current location. For the rest of the line, the index will
    // increment by one
    IndexType idx = it.GetIndex();

    // Displacement vector for the first position in the line and a delta corresponding to a
    // step along the line
    DeformationVectorType disp, delta_disp;

    // Map to a position at which to interpolate
    for(int i = 0; i < ImageDimension; i++)
      {
      disp[i] = off[i] - idx[i];
      delta_disp[i] = M(i, 0);
      for(int j = 0; j < ImageDimension; j++)
        disp[i] += M(i,j) * idx[j];
      }
    delta_disp[0] -= 1.0;

    // Pointer to the start and end of the line
    DeformationVectorType *p_phi = const_cast<DeformationVectorType *>(it.GetPosition());
    DeformationVectorType *p_phi_end = p_phi + outputRegionForThread.GetSize(0);

    // Run a loop filling out the displacement field
    for( ; p_phi < p_phi_end; ++p_phi, disp += delta_disp)
      {
      *p_phi = disp;
      }

    /*
    for( ; p_phi < p_phi_end; ++p_phi, ++idx[0])
      {
      for(int i = 0; i < ImageDimension; i++)
        {
        (*p_phi)[i] = off[i] - idx[i];
        for(int j = 0; j < ImageDimension; j++)
          (*p_phi)[i] += M(i,j) * idx[j];
        }
      }

*/
    }
}


#endif // LINEARTRANSFORMTOWARPFILTER_TXX

