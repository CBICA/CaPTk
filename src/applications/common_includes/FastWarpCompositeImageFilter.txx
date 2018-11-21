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
#ifndef FASTWARPCOMPOSITEIMAGEFILTER_TXX
#define FASTWARPCOMPOSITEIMAGEFILTER_TXX

#include "FastLinearInterpolator.h"
#include "FastWarpCompositeImageFilter.h"
#include "ImageRegionConstIteratorWithIndexOverride.h"

template <class TInputImage, class TOutputImage, class TDeformationField>
void
FastWarpCompositeImageFilter<TInputImage,TOutputImage,TDeformationField>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread,
                       itk::ThreadIdType threadId)
{
  // Our images
  const DeformationFieldType *def = this->GetDeformationField();
  InputImageType *input = this->GetMovingImage();

  int line_len = outputRegionForThread.GetSize(0);

  // Create an iterator over the deformation field
  typedef itk::ImageLinearIteratorWithIndex<OutputImageType> IterBase;
  typedef IteratorExtender<IterBase> IterType;

  // Determine the appropriate float/double type for the interpolator.
  typedef typename itk::NumericTraits<OutputComponentType>::MeasurementVectorType::ValueType FloatType;

  // Create a fast interpolator for the moving image
  typedef FastLinearInterpolator<TInputImage, FloatType, ImageDimension> FastInterpolator;
  FastInterpolator fi(input);
  fi.SetOutsideValue(m_OutsideValue);

  int ncomp = fi.GetPointerIncrement();

  // Loop over the lines in the image
  for(IterType it(this->GetOutput(), outputRegionForThread); !it.IsAtEnd(); it.NextLine())
    {
    // Get the pointer to the displacement vector line and the output vector line
    const DeformationVectorType *phi = it.GetPixelPointer(def);
    OutputComponentType *out = it.GetPixelPointer(this->GetOutput());

    // Voxel index
    IndexType idx = it.GetIndex();

    // The current sample position
    itk::ContinuousIndex<FloatType, ImageDimension> cix;
    typename InputImageType::PointType p, pd, p_step;

    if(m_UsePhysicalSpace)
      {
      // Compute starting point and point step
      this->GetDeformationField()->TransformIndexToPhysicalPoint(idx, p);
      idx[0] += 1;
      this->GetDeformationField()->TransformIndexToPhysicalPoint(idx, p_step);
      for(uint j = 0; j < ImageDimension; j++)
        p_step[j] -= p[j];
      }

    // Loop over the line
    for(int i = 0; i < line_len; i++)
      {
      if(m_UsePhysicalSpace)
        {
        for(uint j = 0; j < ImageDimension; j++)
          {
          pd[j] = p[j] + phi[i][j] * m_DeformationScaling;
          p[j] += p_step[j];
          }

        // TODO: this calls IsInside() internally, which limits efficiency
        input->TransformPhysicalPointToContinuousIndex(pd, cix);
        }
      else
        {
        for(uint j = 0; j < ImageDimension; j++)
          cix[j] = idx[j] + phi[i][j] * m_DeformationScaling;
        idx[0]++;
        }

      // Perform the interpolation
      typename  FastInterpolator::InOut status =
          m_UseNearestNeighbor
          ? fi.InterpolateNearestNeighbor(cix.GetDataPointer(), out)
          : fi.Interpolate(cix.GetDataPointer(), out);

      if(status == FastInterpolator::INSIDE ||
         (status == FastInterpolator::BORDER && m_ExtrapolateBorders))
        {
        out += ncomp;
        }
      else
        {
        for(int k = 0; k < ncomp; k++)
          *out++ = m_OutsideValue;
        }
      }
    }
}

template <class TInputImage, class TOutputImage, class TDeformationField>
void
FastWarpCompositeImageFilter<TInputImage,TOutputImage,TDeformationField>
::GenerateInputRequestedRegion()
{
  this->GetDeformationField()->SetRequestedRegion(this->GetOutput()->GetRequestedRegion());
  this->GetMovingImage()->SetRequestedRegionToLargestPossibleRegion();

}

template <class TInputImage, class TOutputImage, class TDeformationField>
void
FastWarpCompositeImageFilter<TInputImage,TOutputImage,TDeformationField>
::GenerateOutputInformation()
{
  Superclass::GenerateOutputInformation();

  this->GetOutput()->SetNumberOfComponentsPerPixel(this->GetMovingImage()->GetNumberOfComponentsPerPixel());
  this->GetOutput()->SetLargestPossibleRegion(this->GetDeformationField()->GetLargestPossibleRegion());
}



#endif // FASTWARPCOMPOSITEIMAGEFILTER_TXX
