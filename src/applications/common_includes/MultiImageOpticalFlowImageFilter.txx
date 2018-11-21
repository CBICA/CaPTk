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
#ifndef __MultiImageOpticalFlowImageFilter_txx
#define __MultiImageOpticalFlowImageFilter_txx

#include "MultiImageOpticalFlowImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNumericTraits.h"
#include "itkContinuousIndex.h"
#include "vnl/vnl_math.h"
#include "FastLinearInterpolator.h"
#include "ImageRegionConstIteratorWithIndexOverride.h"



/**
 * Compute the output for the region specified by outputRegionForThread.
 */
template <class TMetricTraits>
void
MultiImageOpticalFlowImageFilter<TMetricTraits>
::ThreadedGenerateData(
  const OutputImageRegionType& outputRegionForThread,
  itk::ThreadIdType threadId )
{
  // Get the number of components
  int ncomp = this->GetFixedImage()->GetNumberOfComponentsPerPixel();

  // Create an iterator specialized for going through metrics
  typedef MultiComponentMetricWorker<TMetricTraits, MetricImageType> InterpType;
  InterpType iter(this, this->GetMetricOutput(), outputRegionForThread);

  // Get the per-component metric array
  typename Superclass::ThreadData &td = this->m_ThreadData[threadId];

  // Iterate over the lines
  for(; !iter.IsAtEnd(); iter.NextLine())
    {
    // If we are computing deforamble gradient, get a pointer to the gradient line
    GradientPixelType *grad_line =
        (this->m_ComputeGradient && !this->m_ComputeAffine)
        ? this->GetDeformationGradientOutput()->GetBufferPointer() + iter.GetOffsetInPixels()
        : NULL;

    // Temporary gradient pixel
    GradientPixelType grad_metric;

    // Iterate over the pixels in the line
    for(; !iter.IsAtEndOfLine(); ++iter)
      {
      // Metric value at this pixel
      double metric = 0.0;

      // Gradient vector at this pixel.
      if(this->m_ComputeGradient)
        grad_metric.Fill(0.0);

      // Check the fixed mask at this location
      if(iter.CheckFixedMask())
        {
        // Interpolate the moving image at the current position. The worker knows
        // whether to interpolate the gradient or not
        typedef FastLinearInterpolator<InputImageType, RealType, ImageDimension> FastInterpolator;
        typename FastInterpolator::InOut status = iter.Interpolate();

        // Outside interpolations are ignored
        if(status != FastInterpolator::OUTSIDE)
          {
          // Iterate over the components
          for(int k = 0; k < ncomp; k++)
            {
            // Compute the weighted squared difference
            double del = iter.GetFixedLine()[k] - iter.GetMovingSample()[k];
            double delw = this->m_Weights[k] * del;
            double del2w = delw * del;

            // Add the value to the metric
            metric += del2w;

            // Add the value to each component
            td.comp_metric[k] += del2w;

            // This is currently computing negative half of the gradient
            if(this->m_ComputeGradient)
              {
              const RealType *grad_mov_k = iter.GetMovingSampleGradient(k);
              for(int i = 0; i < ImageDimension; i++)
                grad_metric[i] += delw * grad_mov_k[i];
              }
            }

          // If on the border, apply the mask
          if(status == FastInterpolator::BORDER && this->m_ComputeMovingDomainMask)
            {
            double mask = iter.GetMask();

            // Update the metric and the gradient vector to reflect multiplication with mask
            for(int i = 0; i < ImageDimension; i++)
              grad_metric[i] = grad_metric[i] * mask - 0.5 * iter.GetMaskGradient()[i] * metric;
            metric = metric * mask;

            td.mask += mask;
            }
          else
            {
            td.mask += 1.0;
            }

          // Accumulate the metric
          td.metric += metric;

          // Do the gradient computation
          if(this->m_ComputeGradient)
            {
            // If affine, use the gradient to accumulte the affine partial derivs
            if(this->m_ComputeAffine)
              {
              for(int i = 0, q = 0; i < ImageDimension; i++)
                {
                td.gradient[q++] += grad_metric[i];
                for(int j = 0; j < ImageDimension; j++)
                  td.gradient[q++] += grad_metric[i] * iter.GetIndex()[j];
                }

              if(status == FastInterpolator::BORDER)
                {
                for(int i = 0, q = 0; i < ImageDimension; i++)
                  {
                  td.grad_mask[q++] += iter.GetMaskGradient()[i];
                  for(int j = 0; j < ImageDimension; j++)
                    td.grad_mask[q++] += iter.GetMaskGradient()[i] * iter.GetIndex()[j];
                  }
                }
              }
            }
          } // not outside
        } // check fixed mask

      // Last thing - update the output voxels
      *iter.GetOutputLine() = metric;
      if(grad_line)
        grad_line[iter.GetLinePos()] = grad_metric;
      }
    }
}


#endif
