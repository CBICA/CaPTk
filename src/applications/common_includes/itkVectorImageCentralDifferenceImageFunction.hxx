/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkVectorImageCentralDifferenceImageFunction_hxx
#define itkVectorImageCentralDifferenceImageFunction_hxx

#include "itkVectorImageCentralDifferenceImageFunction.h"

namespace itk
{
/**
 * Constructor
 */
template< typename TInputImage, typename TCoordRep >
VectorImageCentralDifferenceImageFunction< TInputImage, TCoordRep >
::VectorImageCentralDifferenceImageFunction()
{
}

/**
 *
 */
template< typename TInputImage, typename TCoordRep >
void
VectorImageCentralDifferenceImageFunction< TInputImage, TCoordRep >
::PrintSelf(std::ostream & os, Indent indent) const
{
  this->Superclass::PrintSelf(os, indent);
}

/**
 *
 */
template< typename TInputImage, typename TCoordRep >
void VectorImageCentralDifferenceImageFunction<TInputImage, TCoordRep>
::EvaluateAtIndex(const IndexType & index, OutputType &derivative) const
{
  int nc = this->GetInputImage()->GetNumberOfComponentsPerPixel();
  derivative.Fill(0.0);

  IndexType neighIndex = index;

  const typename InputImageType::SizeType & size =
    this->GetInputImage()->GetBufferedRegion().GetSize();
  const typename InputImageType::IndexType & start =
    this->GetInputImage()->GetBufferedRegion().GetIndex();

  int p = 0;
  for ( unsigned int dim = 0; dim < TInputImage::ImageDimension; dim++ )
    {
    // bounds checking
    if ( index[dim] < start[dim] + 1
         || index[dim] > ( start[dim] + static_cast< OffsetValueType >( size[dim] ) - 2 ) )
      {
      p += nc;
      }
    else
      {
      // compute derivative
      const double deriv_weight = 0.5 / this->GetInputImage()->GetSpacing()[dim];

      neighIndex[dim] += 1;
      const InputPixelType pixf = this->GetInputImage()->GetPixel(neighIndex);

      neighIndex[dim] -= 2;
      const InputPixelType pixb = this->GetInputImage()->GetPixel(neighIndex);

      neighIndex[dim] += 1;

      for ( unsigned int vdim = 0; vdim < nc; ++vdim )
        {
        derivative[p++] = (pixf[vdim] - pixb[vdim]) * deriv_weight;
        }
      }
    }
}
} // end namespace itk

#endif
