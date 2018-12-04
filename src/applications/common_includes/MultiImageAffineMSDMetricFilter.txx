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
#ifndef __MultiImageAffineMSDMetricFilter_txx
#define __MultiImageAffineMSDMetricFilter_txx
#include "MultiImageAffineMSDMetricFilter.h"

#define _USING_MASK_ 1


#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNumericTraits.h"
#include "itkContinuousIndex.h"
#include "ImageRegionConstIteratorWithIndexOverride.h"
#include "vnl/vnl_math.h"
#include "FastLinearInterpolator.h"

/**
 * Setup state of filter before multi-threading.
 * InterpolatorType::SetInputImage is not thread-safe and hence
 * has to be setup before ThreadedGenerateData
 */
template <class TMetricTraits>
void
MultiImageAffineMetricFilter<TMetricTraits>
::BeforeThreadedGenerateData()
{
  // Initialize the per thread data
  m_ThreadData.resize(this->GetNumberOfThreads(), ThreadData());
}

/**
 * Setup state of filter after multi-threading.
 */
template <class TMetricTraits>
void
MultiImageAffineMetricFilter<TMetricTraits>
::AfterThreadedGenerateData()
{
  // Add up all the thread data
  ThreadData summary;
  for(int i = 0; i < m_ThreadData.size(); i++)
    {
    summary.metric += m_ThreadData[i].metric;
    summary.mask += m_ThreadData[i].mask;
    summary.gradient += m_ThreadData[i].gradient;
    summary.grad_mask += m_ThreadData[i].grad_mask;
    }

  // Initialize the output metric gradient
  m_MetricGradient = TransformType::New();

  // Compute the objective value
#ifdef _USING_MASK_
  m_MetricValue = summary.metric / summary.mask;
#else
  m_MetricValue = summary.metric;
#endif

  // Compute the gradient
  vnl_vector<double> grad_metric(summary.gradient.size());
  for(int j = 0; j < summary.gradient.size(); j++)
    {
#ifdef _USING_MASK_
    grad_metric[j] =
        (m_GradientScalingFactor * summary.gradient[j] - m_MetricValue * summary.grad_mask[j]) / summary.mask;
#else
    grad_metric[j] = m_GradientScalingFactor * summary.gradient[j];
#endif
    }

  // Pack into the output
  unflatten_affine_transform(grad_metric.data_block(), m_MetricGradient.GetPointer());
}

template <class TMetricTraits>
void
MultiImageAffineMetricFilter<TMetricTraits>
::GenerateInputRequestedRegion()
{
  // Call the superclass's implementation
  Superclass::GenerateInputRequestedRegion();

  // Set all regions to max
  this->GetMetricImage()->SetRequestedRegionToLargestPossibleRegion();
  this->GetMovingDomainMaskImage()->SetRequestedRegionToLargestPossibleRegion();

  if(m_ComputeGradient)
    {
    this->GetGradientImage()->SetRequestedRegionToLargestPossibleRegion();
    this->GetMovingDomainMaskGradientImage()->SetRequestedRegionToLargestPossibleRegion();
    }
}

template <class TMetricTraits>
void
MultiImageAffineMetricFilter<TMetricTraits>
::EnlargeOutputRequestedRegion(itk::DataObject *data)
{
  Superclass::EnlargeOutputRequestedRegion(data);
  data->SetRequestedRegionToLargestPossibleRegion();
}


template <class TMetricTraits>
void
MultiImageAffineMetricFilter<TMetricTraits>
::AllocateOutputs()
{
  // Propagate input
  this->GraftOutput(this->GetMetricImage());
}


template <class TMetricTraits>
void
MultiImageAffineMetricFilter<TMetricTraits>
::ThreadedGenerateData(
  const OutputImageRegionType& outputRegionForThread,
  itk::ThreadIdType threadId )
{
  // Get the pointers to the input and output images
  InputImageType *metric = this->GetMetricImage();
  InputImageType *mask = this->GetMovingDomainMaskImage();

  // Get the pointers to the buffers of these four images
  const InputPixelType *b_metric = metric->GetBufferPointer();
  const InputPixelType *b_mask = mask->GetBufferPointer();

  const GradientPixelType *b_gradient =
      m_ComputeGradient ? this->GetGradientImage()->GetBufferPointer() : NULL;

  const GradientPixelType *b_mask_gradient =
      m_ComputeGradient ? this->GetMovingDomainMaskGradientImage()->GetBufferPointer() : NULL;

  // Iterate over the deformation field and the output image. In reality, we don't
  // need to waste so much time on iteration, so we use a specialized iterator here
  typedef itk::ImageLinearConstIteratorWithIndex<InputImageType> IterBase;
  typedef IteratorExtender<IterBase> Iter;

  // The thread data to accumulate
  ThreadData &td = m_ThreadData[threadId];

  // Gradient accumulator
  vnl_vector_fixed<double, ImageDimension> grad, gradM;

  itk::ImageRegion<ImageDimension> workRegion = outputRegionForThread;
  // workRegion.ShrinkByRadius(3);

  // Iterate over the fixed space region
  for(Iter it(metric, workRegion); !it.IsAtEnd(); it.NextLine())
    {
    // Process the whole line using pointer arithmetic. We have to deal with messy behavior
    // of iterators on vector images. Trying to avoid using accessors and Set/Get
    long offset_in_pixels = it.GetPosition() - b_metric;

    // Get pointers to the start of the line
    const InputPixelType *p_metric = b_metric + offset_in_pixels;
    const InputPixelType *p_mask = b_mask + offset_in_pixels;

    // Loop over the line
    const InputPixelType *p_metric_end = p_metric + outputRegionForThread.GetSize(0);

    // Loop if we are computing gradient
    if(m_ComputeGradient)
      {
      const GradientPixelType *p_gradient = b_gradient + offset_in_pixels;
      const GradientPixelType *p_mask_gradient = b_mask_gradient + offset_in_pixels;

      // Get the index at the current location
      IndexType idx = it.GetIndex();

      // Do we need the gradient?
      for(; p_metric < p_metric_end; ++p_metric, ++p_mask, ++p_gradient, ++p_mask_gradient, ++idx[0])
        {
        // Accumulators for the gradients
        double *out_grad = td.gradient.data_block();
        double *out_grad_mask = td.grad_mask.data_block();

        // For border regions, we need to explicitly deal with the mask
        const InputPixelType &metric = *p_metric;
        const InputPixelType &mask = *p_mask;
        const GradientPixelType &grad = *p_gradient;
        const GradientPixelType &gradM = *p_mask_gradient;

        // Compute the mask and metric gradient contributions
#ifdef _USING_MASK_
        for(int i = 0; i < ImageDimension; i++)
          {
          double v = grad[i] * mask - 0.5 * gradM[i] * metric;
          *(out_grad++) += v;
          *(out_grad_mask++) += gradM[i];
          for(int j = 0; j < ImageDimension; j++)
            {
            *(out_grad++) += v * idx[j];
            *(out_grad_mask++) += gradM[i] * idx[j];
            }
          }

        td.metric += metric * mask;
        td.mask += mask;
#else
        // NOT USING MASK
        for(int i = 0; i < ImageDimension; i++)
          {
          double v = grad[i];
          *(out_grad++) += v;
          for(int j = 0; j < ImageDimension; j++)
            {
            *(out_grad++) += v * idx[j];
            }
          }

        td.metric += metric;
#endif
        }
      }
    else
      {
      // Do we need the gradient?
      for(; p_metric < p_metric_end; ++p_metric, ++p_mask)
        {
        // For border regions, we need to explicitly deal with the mask
        const InputPixelType &metric = *p_metric;
        const InputPixelType &mask = *p_mask;
        td.metric += metric * mask;
        td.mask += mask;
        }
      }
    }
}


#endif
