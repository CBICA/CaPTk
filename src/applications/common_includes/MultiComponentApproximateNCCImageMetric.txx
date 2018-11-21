/*=========================================================================

  Program:   ALFABIS fast image registration
  Language:  C++

  Copyright (c) Paul Yushkevich. All rights reserved.

  This program is part of ALFABIS

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
#ifndef __MultiComponentApproximateNCCImageMetric_txx
#define __MultiComponentApproximateNCCImageMetric_txx

#include "MultiComponentApproximateNCCImageMetric.h"
#include "OneDimensionalInPlaceAccumulateFilter.h"
#include "itkImageFileWriter.h"



/* ==========================================================================
 *
 *                        PRECOMPUTE filter implementation
 *
 * ========================================================================== */

template <class TMetricTraits, class TOutputImage>
MultiImageApproximateNCCPrecomputeFilter<TMetricTraits,TOutputImage>
::MultiImageApproximateNCCPrecomputeFilter()
{
  m_Parent = NULL;
  m_Stage = FIRST;
}

/**
 * Generate output information, which will be different from the default
 */
template <class TMetricTraits, class TOutputImage>
void
MultiImageApproximateNCCPrecomputeFilter<TMetricTraits,TOutputImage>
::GenerateOutputInformation()
{
  // Call the parent method to set up all the outputs
  Superclass::GenerateOutputInformation();

  // Set the number of components in the primary output
  int ncomp = this->GetNumberOfOutputComponents();
  this->GetOutput()->SetNumberOfComponentsPerPixel(ncomp);
}

template <class TMetricTraits, class TOutputImage>
int
MultiImageApproximateNCCPrecomputeFilter<TMetricTraits,TOutputImage>
::GetNumberOfOutputComponents()
{
  // This is complex! The number of components depends on what we are computing.
  int nc = m_Parent->GetFixedImage()->GetNumberOfComponentsPerPixel();

  // If there are no gradients computed, we just need 7 output pixels per input comp.
  return 1 + nc * 7;
}

/**
 * Compute the output for the region specified by outputRegionForThread.
 *
 * Coding:
 *    1
 *    f
 *    m
 *    f
 *    m
 *    f*f
 *    m*m
 *    f*m
 *    ...
 */
template <class TMetricTraits, class TOutputImage>
void
MultiImageApproximateNCCPrecomputeFilter<TMetricTraits,TOutputImage>
::ThreadedGenerateData(
  const OutputImageRegionType& outputRegionForThread,
  itk::ThreadIdType threadId )
{
  typedef FastLinearInterpolator<InputImageType, RealType, ImageDimension> FastInterpolator;

  // Get the number of input and output components
  int ncomp_in = m_Parent->GetFixedImage()->GetNumberOfComponentsPerPixel();

  // Create an iterator specialized for going through metrics
  typedef MultiComponentMetricWorker<TMetricTraits, TOutputImage> InterpType;
  InterpType iter(m_Parent, this->GetOutput(), outputRegionForThread);

  // Iterate over the lines
  for(; !iter.IsAtEnd(); iter.NextLine())
    {
    // Iterate over the pixels in the line
    for(; !iter.IsAtEndOfLine(); ++iter)
      {
      // Get the output pointer for this voxel
      OutputComponentType *out = iter.GetOutputLine();

      // Interpolate the moving image at the current position. The worker knows
      // whether to interpolate the gradient or not
      typename FastInterpolator::InOut status = iter.Interpolate();

      // Border should be ok here because we are not sampling the gradient
      if(status == FastInterpolator::OUTSIDE)
        {
        // Place a zero in the output
        *out++ = 1.0;

        // Iterate over the components
        for(int k = 0; k < ncomp_in; k++)
          {
          InputComponentType x_fix = iter.GetFixedLine()[k];
          *out++ = x_fix;
          *out++ = 0.0;
          *out++ = x_fix;
          *out++ = 0.0;
          *out++ = x_fix * x_fix;
          *out++ = 0.0;
          *out++ = 0.0;
          }
        }

      else
        {
        if(m_Stage == FIRST)
          {
          // Place a zero in the output
          *out++ = 1.0;

          // Iterate over the components
          for(int k = 0; k < ncomp_in; k++)
            {
            InputComponentType x_mov = iter.GetMovingSample()[k];
            InputComponentType x_fix = iter.GetFixedLine()[k];

            // Write the five components
            *out++ = x_fix;
            *out++ = x_mov;
            *out++ = 0.0;
            *out++ = 0.0;
            *out++ = 0.0;
            *out++ = 0.0;
            *out++ = 0.0;
            }
          }
        else
          {
          // Get the number of neighbors from last round
          InputComponentType n_nbr = *out;

          // Place a one in the output
          out++;

          // Iterate over the components
          for(int k = 0; k < ncomp_in; k++)
            {
            InputComponentType x_mov_raw = iter.GetMovingSample()[k];
            InputComponentType x_fix_raw = iter.GetFixedLine()[k];

            // Write the five components
            InputComponentType x_fix = x_fix_raw - *out / n_nbr;
            *out++ = x_fix;
            InputComponentType x_mov = x_mov_raw - *out / n_nbr;
            *out++ = x_mov;
            *out++ = x_fix;
            *out++ = x_mov;
            *out++ = x_fix * x_fix;
            *out++ = x_mov * x_mov;
            *out++ = x_fix * x_mov;
            }
          }
        }
      }
    }
}



/* ==========================================================================
 *
 *                        main filter implementation
 *
 * ========================================================================== */



/**
 * This is a similar function to above, but uses the approximate algorithm described by Avants et al.
 * in the 2008 NeuroImage paper. This is not the exact gradient of the metric. This is useful for
 * comparing Greedy with ANTS
 */
template <class TPixel, class TWeight, class TMetric, class TGradient>
TPixel *
MultiImageApproximateNNCPostComputeFunction(
    TPixel *ptr, TPixel *ptr_end, TWeight *weights, TMetric *ptr_metric, TGradient *ptr_gradient, int ImageDimension,
    const TPixel *grad_fix, bool debug)
{
  // Get the size of the mean filter kernel
  TPixel n = *ptr++, one_over_n = 1.0 / n;

  // Loop over components
  int i_wgt = 0;

  // Initialize metric to zero
  *ptr_metric = 0;

  for(; ptr < ptr_end; ++i_wgt)
    {
    TPixel I_bar = *ptr++;
    TPixel J_bar = *ptr++;
    TPixel I_2bar = *ptr++;
    TPixel J_2bar = *ptr++;
    TPixel x_fix_sq = *ptr++;
    TPixel x_mov_sq = *ptr++;
    TPixel x_fix_mov = *ptr++;

    TPixel sff = x_fix_sq - I_2bar * I_2bar * one_over_n;
    TPixel smm = x_mov_sq - J_2bar * J_2bar * one_over_n;
    TPixel smf = x_fix_mov - I_2bar * J_2bar * one_over_n;


    if(sff == 0 || smm == 0)
      {
      if(ptr_gradient)
        ptr += 3 * ImageDimension;
      continue;
      }

    TWeight w = weights[i_wgt];
    TPixel zoop = w * (smf / (sff * smm));

    if(ptr_gradient)
      {

      // Term to multiply the gradient by...
      TPixel factor = -2.0 * zoop * (J_bar - (smf / sff) * I_bar);

      //  2.0 * sfm / (sff * smm) * ( Ji - sfm / sff * Ii )

      if(debug)
        {
        printf("Ii = %f, Ji = %f, sff = %f, sfm = %f, smm = %f\n", I_bar, J_bar, sff, smf, smm);
        printf("Metric = %f\n", zoop * smf);
        printf("GradF = %f, %f, %f\n", grad_fix[0], grad_fix[1], grad_fix[2]);
        }

      for(int i = 0; i < ImageDimension; i++)
        {
        (*ptr_gradient)[i] += factor * (*grad_fix++);
        }

      if(debug)
        {
        printf("Deriv = %g, %g, %g\n", (*ptr_gradient)[0], (*ptr_gradient)[1], (*ptr_gradient)[2]);
        }
      }


    // Accumulate the metric
    *ptr_metric += zoop * smf;
    }

  return ptr;
}


// #define DUMP_NCC 1

template <class TMetricTraits>
void
MultiComponentApproximateNCCImageMetric<TMetricTraits>
::BeforeThreadedGenerateData()
{
  // Call the parent method
  Superclass::BeforeThreadedGenerateData();

  // Pre-compute filter 1
  typedef MultiImageApproximateNCCPrecomputeFilter<TMetricTraits, InputImageType> PreFilterType;
  typename PreFilterType::Pointer preFilter1 = PreFilterType::New();

  // Configure the precompute filter
  preFilter1->SetParent(this);
  preFilter1->SetInput(this->GetFixedImage());

  // Number of components in the working image
  int ncomp = preFilter1->GetNumberOfOutputComponents();

  // If the user supplied a working image, configure it and graft it as output
  if(m_WorkingImage)
    {
    // Configure the working image
    m_WorkingImage->CopyInformation(this->GetFixedImage());
    m_WorkingImage->SetNumberOfComponentsPerPixel(ncomp);
    m_WorkingImage->SetRegions(this->GetFixedImage()->GetBufferedRegion());
    m_WorkingImage->Allocate();

    // Graft the working image onto the filter's output
    preFilter1->GraftOutput(m_WorkingImage);
    }

  // Execute the filter
  preFilter1->SetStage(PreFilterType::FIRST);
  preFilter1->Update();

#ifdef DUMP_NCC
  typename itk::ImageFileWriter<InputImageType>::Pointer pwriter = itk::ImageFileWriter<InputImageType>::New();
  pwriter->SetInput(preFilter1->GetOutput());
  pwriter->SetFileName("nccpre1.nii.gz");
  pwriter->Update();
#endif

  // First round of accumulation
  typedef OneDimensionalInPlaceAccumulateFilter<InputImageType> AccumFilterType;

  // Create a chain of separable 1-D filters
  typename itk::ImageSource<InputImageType>::Pointer pipeTail;
  for(int dir = 0; dir < ImageDimension; dir++)
    {
    typename AccumFilterType::Pointer accum = AccumFilterType::New();
    if(pipeTail.IsNull())
      accum->SetInput(preFilter1->GetOutput());
    else
      accum->SetInput(pipeTail->GetOutput());
    accum->SetComponentRange(0, 5);
    accum->SetDimension(dir);
    accum->SetRadius(m_Radius[dir]);
    pipeTail = accum;

    accum->Update();
    }

  // At this point, we want to pull out components 1-3 from the working image

  // Perform the second round of accumulation
  typename PreFilterType::Pointer preFilter2 = PreFilterType::New();

  // Configure the precompute filter
  preFilter2->SetParent(this);
  preFilter2->SetInput(this->GetFixedImage());
  preFilter2->GraftOutput(m_WorkingImage);

  // Execute the filter
  preFilter2->SetStage(PreFilterType::SECOND);
  preFilter2->Update();

#ifdef DUMP_NCC
  typename itk::ImageFileWriter<InputImageType>::Pointer pwriter1 = itk::ImageFileWriter<InputImageType>::New();
  pwriter1->SetInput(preFilter2->GetOutput());
  pwriter1->SetFileName("nccpre2.nii.gz");
  pwriter1->Update();
#endif

  // Create a chain of separable 1-D filters
  pipeTail = NULL;
  for(int dir = 0; dir < ImageDimension; dir++)
    {
    typename AccumFilterType::Pointer accum = AccumFilterType::New();
    if(pipeTail.IsNull())
      accum->SetInput(preFilter2->GetOutput());
    else
      accum->SetInput(pipeTail->GetOutput());
    accum->SetComponentRange(3, 0);
    accum->SetDimension(dir);
    accum->SetRadius(m_Radius[dir]);
    pipeTail = accum;

    accum->Update();
    }

#ifdef DUMP_NCC
  typename itk::ImageFileWriter<InputImageType>::Pointer pwriter2 = itk::ImageFileWriter<InputImageType>::New();
  pwriter2->SetInput(pipeTail->GetOutput());
  pwriter2->SetFileName("nccaccum.nii.gz");
  pwriter2->Update();
#endif
}

// TODO: this only needs to be computed once per resolution level
#include <itkVectorImageCentralDifferenceImageFunction.h>

template <class TMetricTraits>
void
MultiComponentApproximateNCCImageMetric<TMetricTraits>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
  int nc = m_WorkingImage->GetNumberOfComponentsPerPixel();
  int line_len = outputRegionForThread.GetSize()[0];

  // Our thread data
  typename Superclass::ThreadData &td = this->m_ThreadData[threadId];

  // Set up an iterator for the working image
  typedef itk::ImageLinearConstIteratorWithIndex<InputImageType> InputIteratorTypeBase;
  typedef IteratorExtender<InputIteratorTypeBase> InputIteratorType;
  InputIteratorType it(m_WorkingImage, outputRegionForThread);

  // TODO: this can be pre-computed!
  /* BEGIN: fixed gradient comp */
  // Set up the fixed gradient calculator

  typedef itk::VectorImageCentralDifferenceImageFunction<InputImageType, RealType> CDType;
  typename CDType::Pointer cdiff = CDType::New();
  cdiff->SetInputImage(this->GetFixedImage());
  IndexType index;
  int grad_fixed_length = ImageDimension * this->GetFixedImage()->GetNumberOfComponentsPerPixel();
  RealType *grad_fixed = new RealType[grad_fixed_length];
  itk::VariableLengthVector<RealType> grad_fixed_vec(grad_fixed, grad_fixed_length);

  // Loop over the lines
  for (; !it.IsAtEnd(); it.NextLine())
    {
    // Get the pointer to the input line
    long offset_in_pixels = it.GetPosition() - m_WorkingImage->GetBufferPointer();

    // Pointer to the input pixel data for this line
    InputComponentType *p_input = m_WorkingImage->GetBufferPointer() + nc * offset_in_pixels;

    // Pointer to the metric data for this line
    MetricPixelType *p_metric = this->GetMetricOutput()->GetBufferPointer() + offset_in_pixels;

    // The gradient output is optional
    GradientPixelType *p_grad_metric = (this->m_ComputeGradient && !this->m_ComputeAffine)
                                       ? this->GetDeformationGradientOutput()->GetBufferPointer() + offset_in_pixels
                                       : NULL;

    // Index: for gradient computation (remove later)
    index = it.GetIndex();

    // Case 1 - dense gradient field requested
    if(!this->m_ComputeAffine)
      {
      if(this->m_ComputeGradient)
        {
        // Loop over the pixels in the line
        for(int i = 0; i < line_len; ++i)
          {
          // Clear the metric and the gradient
          *p_metric = itk::NumericTraits<MetricPixelType>::Zero;

          // Gradient of the fixed image is passed in through the p_grad_metric!
          *p_grad_metric = itk::NumericTraits<GradientPixelType>::Zero;

          // Compute the gradient
          index[0] = i;
          cdiff->EvaluateAtIndex(index, grad_fixed_vec);

          bool debugff = (it.GetIndex()[1] == 33 && it.GetIndex()[2] == 22 && i == 44);
          p_input = MultiImageApproximateNNCPostComputeFunction(p_input, p_input + nc, this->m_Weights.data_block(),
                                                                p_metric, p_grad_metric++, ImageDimension, grad_fixed, debugff);
          // Accumulate the total metric
          td.metric += *p_metric;
          td.mask += 1.0;
          }
        }
      else
        {
        // Loop over the pixels in the line
        for(int i = 0; i < line_len; ++i)
          {
          // Clear the metric and the gradient
          *p_metric = itk::NumericTraits<MetricPixelType>::Zero;

          // Apply the post computation
          p_input = MultiImageApproximateNNCPostComputeFunction(p_input, p_input + nc, this->m_Weights.data_block(),
                                                                p_metric, (GradientPixelType *)(NULL), ImageDimension,
                                                                (const RealType *) NULL, false);

          // Accumulate the total metric
          td.metric += *p_metric;
          td.mask += 1.0;
          }
        }
      }
    }

  delete[] grad_fixed;
}


#endif // __MultiComponentApproximateNCCImageMetric_txx

