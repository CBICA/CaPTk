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
#ifndef __MultiComponentNCCImageMetric_txx
#define __MultiComponentNCCImageMetric_txx

#include "MultiComponentNCCImageMetric.h"
#include "OneDimensionalInPlaceAccumulateFilter.h"
#include "itkImageFileWriter.h"



/* ==========================================================================
 *
 *                        PRECOMPUTE filter implementation
 *
 * ========================================================================== */

template <class TMetricTraits, class TOutputImage>
MultiImageNCCPrecomputeFilter<TMetricTraits,TOutputImage>
::MultiImageNCCPrecomputeFilter()
{
  m_Parent = NULL;
  m_FlagGenerateFixedComponents = true;
}

/**
 * Generate output information, which will be different from the default
 */
template <class TMetricTraits, class TOutputImage>
void
MultiImageNCCPrecomputeFilter<TMetricTraits,TOutputImage>
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
MultiImageNCCPrecomputeFilter<TMetricTraits,TOutputImage>
::GetNumberOfOutputComponents()
{
  // This is complex! The number of components depends on what we are computing.
  int nc = m_Parent->GetFixedImage()->GetNumberOfComponentsPerPixel();

  // If there are no gradients computed, we just need 5 output pixels per input comp.
  if(!m_Parent->GetComputeGradient())
    return 1 + nc * 5;

  // If we are not computing affine transform, then the gradient requires 3 comps per dimension
  if(!m_Parent->GetComputeAffine())
    return 1 + nc * (5 + 3 * ImageDimension);

  // Otherwise we need a ton of components, because for each gradient component we need it also
  // scaled by x, by y and by z.
  return 1 + nc * (5 + 3 * ImageDimension * (1 + ImageDimension));
}

template <class TMetricTraits, class TOutputImage>
int
MultiImageNCCPrecomputeFilter<TMetricTraits,TOutputImage>
::GetNumberOfFixedOnlyOutputComponents()
{
  // This is complex! The number of components depends on what we are computing.
  int nc = m_Parent->GetFixedImage()->GetNumberOfComponentsPerPixel();
  return 1 + nc * 2;
}

/**
 * Compute the output for the region specified by outputRegionForThread.
 */
template <class TMetricTraits, class TOutputImage>
void
MultiImageNCCPrecomputeFilter<TMetricTraits,TOutputImage>
::ThreadedGenerateData(
  const OutputImageRegionType& outputRegionForThread,
  itk::ThreadIdType threadId )
{
  typedef FastLinearInterpolator<InputImageType, RealType, ImageDimension> FastInterpolator;

  // Get the number of input and output components
  int ncomp_in = m_Parent->GetFixedImage()->GetNumberOfComponentsPerPixel();
  int ncomp_out = this->GetNumberOfOutputComponents();

  // Number of invariant components (fixed-only)
  int ncomp_fixed = ncomp_in * 2 + 1;

  // Get the pointer to the output image
  OutputImageType *out = this->GetOutput();

  // Create an iterator specialized for going through metrics
  typedef MultiComponentMetricWorker<TMetricTraits, TOutputImage> InterpType;
  InterpType iter(m_Parent, this->GetOutput(), outputRegionForThread);

  // Iterate over the lines
  for(; !iter.IsAtEnd(); iter.NextLine())
    {
    // Iterate over the pixels in the line
    for(; !iter.IsAtEndOfLine(); ++iter)
      {
      // The mask is 0.5 for points in range of the user mask, and 1.0 for the user mask
      // We can safely ignore the points outside of range - they have no impact on the
      // region that we are going to measure
      if(iter.CheckFixedMask(0.0))
        {
        // Get the output pointer for this voxel
        OutputComponentType *out;

        // Is this the first time the filter is being run for this output working image. If
        // so, we store the fixed components for future accumulation.
        if(this->m_FlagGenerateFixedComponents)
          {
          out = iter.GetOutputLine();
          // Store the components that are invariant of the moving image at the beginning. This
          // helps avoid having to do repeated accumulation on these components
          *out++ = 1.0;
          for(int k = 0; k < ncomp_in; k++)
            {
            InputComponentType x_fix = iter.GetFixedLine()[k];
            *out++ = x_fix;
            *out++ = x_fix * x_fix;
            }
          }
        else
          {
          out = iter.GetOutputLine() + ncomp_fixed;
          }

        // Interpolate the moving image at the current position. The worker knows
        // whether to interpolate the gradient or not
        typename FastInterpolator::InOut status = iter.Interpolate();

        // TODO: debugging border issues: setting mask=0 on the border
        // TODO: i added this check for affine/non-affine. Seems like there were some problems previously
        // with the NCC metric at the border in deformable mode, but in affine mode you need the border
        // values to be included.
        // Outside interpolations are ignored
        if((m_Parent->GetComputeAffine() && status == FastInterpolator::OUTSIDE) ||
           (!m_Parent->GetComputeAffine() && (status == FastInterpolator::OUTSIDE || status == FastInterpolator::BORDER)))
          {
          // Iterate over the components
          for(int k = 0; k < ncomp_in; k++)
            {
            InputComponentType x_fix = iter.GetFixedLine()[k];
            *out++ = 0.0;
            *out++ = 0.0;
            *out++ = 0.0;

            int n = m_Parent->GetComputeGradient()
                    ? ( m_Parent->GetComputeAffine()
                        ? 3 * ImageDimension * (1 + ImageDimension)
                        : 3 * ImageDimension)
                    : 0;

            for(int j = 0; j < n; j++)
              *out++ = 0.0;
            }
          }
        else
          {
          // Iterate over the components
          for(int k = 0; k < ncomp_in; k++)
            {
            InputComponentType x_mov = iter.GetMovingSample()[k];
            InputComponentType x_fix = iter.GetFixedLine()[k];

            // Write the five components
            *out++ = x_mov;
            *out++ = x_mov * x_mov;
            *out++ = x_fix * x_mov;

            // If gradient, do more
            if(m_Parent->GetComputeGradient())
              {
              const InputComponentType *mov_grad = iter.GetMovingSampleGradient(k);
              for(int i = 0; i < ImageDimension; i++)
                {
                InputComponentType mov_grad_i = mov_grad[i];
                *out++ = mov_grad_i;
                *out++ = x_fix * mov_grad_i;
                *out++ = x_mov * mov_grad_i;

                // If affine, tack on additional components
                if(m_Parent->GetComputeAffine())
                  {
                  for(int j = 0; j < ImageDimension; j++)
                    {
                    double x = iter.GetIndex()[j];
                    *out++ = mov_grad_i * x;
                    *out++ = x_fix * mov_grad_i * x;
                    *out++ = x_mov * mov_grad_i * x;
                    }
                  }
                }


              }
            }
          }
        }
      else
        {
        // TODO: do we need this? Zero out the images
        OutputComponentType *out = iter.GetOutputLine();
        for(int q = 0; q < ncomp_out; q++)
          *out++ = 0;
        }
      }
    }
}




/* ==========================================================================
 *
 *                        main filter implementation
 *
 * ========================================================================== */


template <class TPixel, class TWeight, class TMetric, class TGradient>
TPixel *
MultiImageNNCPostComputeFunction(
    TPixel *ptr, TPixel *ptr_end, int n_comp, TWeight *weights, TMetric *ptr_metric, TGradient *ptr_gradient, int ImageDimension)
{
  // IMPORTANT: this code uses double precision because single precision float seems
  // to mess up and lead to unstable computations

  // Get the size of the mean filter kernel
  double n = *ptr++, one_over_n = 1.0 / n;

  // Get the pointer for the fixed only components, and for the moving components
  TPixel *ptr_fix = ptr; ptr += n_comp * 2;

  // Loop over components
  int i_wgt = 0;
  const double eps = 1e-8;

  // Initialize metric to zero
  *ptr_metric = 0;

  for(; ptr < ptr_end; ++i_wgt)
    {
    double x_fix = *ptr_fix++;
    double x_mov = *ptr++;
    double x_fix_sq = *ptr_fix++;
    double x_mov_sq = *ptr++;
    double x_fix_mov = *ptr++;

    double x_fix_over_n = x_fix * one_over_n;
    double x_mov_over_n = x_mov * one_over_n;

    double var_fix = x_fix_sq - x_fix * x_fix_over_n;
    double var_mov = x_mov_sq - x_mov * x_mov_over_n;

    if(var_fix < eps || var_mov < eps)
      {
      if(ptr_gradient)
        ptr += 3 * ImageDimension;
      continue;
      }

    double cov_fix_mov = x_fix_mov - x_fix * x_mov_over_n;
    double one_over_denom = 1.0 / (var_fix * var_mov);
    double cov_fix_mov_over_denom = cov_fix_mov * one_over_denom;
    double ncc_fix_mov = cov_fix_mov * cov_fix_mov_over_denom;

    // Weight - includes scaling of squared covariance by direction
    TWeight w = (cov_fix_mov < 0) ? -weights[i_wgt] : weights[i_wgt];

    if(ptr_gradient)
      {
      for(int i = 0; i < ImageDimension; i++)
        {
        double x_grad_mov_i = *ptr++;
        double x_fix_grad_mov_i = *ptr++;
        double x_mov_grad_mov_i = *ptr++;

        // Derivative of cov_fix_mov
        double grad_cov_fix_mov_i = x_fix_grad_mov_i - x_fix_over_n * x_grad_mov_i;

        // One half derivative of var_mov
        double half_grad_var_mov_i = x_mov_grad_mov_i - x_mov_over_n * x_grad_mov_i;

        double grad_ncc_fix_mov_i =
            2 * cov_fix_mov_over_denom * (grad_cov_fix_mov_i - var_fix * half_grad_var_mov_i * cov_fix_mov_over_denom);

        if(fabs(grad_ncc_fix_mov_i) > 1000)
          {
          printf("big value: %f\n", grad_ncc_fix_mov_i);
          }

        (*ptr_gradient)[i] += (TPixel) (w * grad_ncc_fix_mov_i);
        }
      }

    // Accumulate the metric
    *ptr_metric += (TPixel) (w * ncc_fix_mov);
    }

  return ptr;
}



template <class TPixel, class TWeight, class TMetric, class TGradient>
TPixel *
MultiImageNNCPostComputeAffineGradientFunction(
    TPixel *ptr, TPixel *ptr_end, int n_comp, TWeight *weights, TMetric *ptr_metric, TGradient *ptr_affine_gradient, int ImageDimension)
{
  // Get the size of the mean filter kernel
  TPixel n = *ptr++, one_over_n = 1.0 / n;

  // Get the pointer for the fixed only components, and for the moving components
  TPixel *ptr_fix = ptr; ptr += n_comp * 2;

  // Loop over components
  int i_wgt = 0;
  const TPixel eps = 1e-2;

  // Initialize metric to zero
  *ptr_metric = 0;

  for(; ptr < ptr_end; ++i_wgt)
    {
    TPixel x_fix = *ptr_fix++;
    TPixel x_mov = *ptr++;
    TPixel x_fix_sq = *ptr_fix++;
    TPixel x_mov_sq = *ptr++;
    TPixel x_fix_mov = *ptr++;

    TPixel x_fix_over_n = x_fix * one_over_n;
    TPixel x_mov_over_n = x_mov * one_over_n;

    TPixel var_fix = x_fix_sq - x_fix * x_fix_over_n + eps;
    TPixel var_mov = x_mov_sq - x_mov * x_mov_over_n + eps;

    TPixel cov_fix_mov = x_fix_mov - x_fix * x_mov_over_n + eps;
    TPixel one_over_denom = 1.0 / (var_fix * var_mov);
    TPixel cov_fix_mov_over_denom = cov_fix_mov * one_over_denom;
    TPixel ncc_fix_mov = cov_fix_mov * cov_fix_mov_over_denom;

    // Weight - includes scaling of squared covariance by direction
    TWeight w = (cov_fix_mov < 0) ? -weights[i_wgt] : weights[i_wgt];

    if(ptr_affine_gradient)
      {
      // There are 12 components to compute
      TGradient *paff = ptr_affine_gradient;

      for(int i = 0; i < ImageDimension; i++)
        {
        for(int j = 0; j <= ImageDimension; j++)
          {
          TPixel x_grad_mov_i = *ptr++;
          TPixel x_fix_grad_mov_i = *ptr++;
          TPixel x_mov_grad_mov_i = *ptr++;

          // Derivative of cov_fix_mov
          TPixel grad_cov_fix_mov_i = x_fix_grad_mov_i - x_fix_over_n * x_grad_mov_i;

          // One half derivative of var_mov
          TPixel half_grad_var_mov_i = x_mov_grad_mov_i - x_mov_over_n * x_grad_mov_i;

          TPixel grad_ncc_fix_mov_i =
              2 * cov_fix_mov_over_denom * (grad_cov_fix_mov_i - var_fix * half_grad_var_mov_i * cov_fix_mov_over_denom);

          *paff++ += w * grad_ncc_fix_mov_i;
          }
        }
      }

    // Accumulate the metric
    *ptr_metric += w * ncc_fix_mov;
    }

  return ptr;
}




template <class TMetricTraits>
void
MultiComponentNCCImageMetric<TMetricTraits>
::BeforeThreadedGenerateData()
{
  // Call the parent method
  Superclass::BeforeThreadedGenerateData();

  // Pre-compute filter
  typedef MultiImageNCCPrecomputeFilter<TMetricTraits, InputImageType> PreFilterType;
  typename PreFilterType::Pointer preFilter = PreFilterType::New();

  // Configure the precompute filter
  preFilter->SetParent(this);
  preFilter->SetInput(this->GetFixedImage());

  // Number of components in the working image
  unsigned ncomp = preFilter->GetNumberOfOutputComponents();
  unsigned ncomp_ignore = 0;

  // If the working image is supplied, make sure that it has sufficient size
  bool flag_reuse = m_ReuseWorkingImageFixedComponents;

  if(m_WorkingImage)
    {
    // Check if the working image needs to be allocated
    if(m_WorkingImage->GetBufferedRegion() != this->GetFixedImage()->GetBufferedRegion()
       || m_WorkingImage->GetNumberOfComponentsPerPixel() < ncomp)
      {
      // Configure the working image
      m_WorkingImage->CopyInformation(this->GetFixedImage());
      m_WorkingImage->SetNumberOfComponentsPerPixel(ncomp);
      m_WorkingImage->SetRegions(this->GetFixedImage()->GetBufferedRegion());
      m_WorkingImage->Allocate();

      // Can't reuse the fixed components
      flag_reuse = false;
      }

    // Graft the working image onto the filter's output
    preFilter->GraftOutput(m_WorkingImage);
    }
  else
    {
    // No working image, can't reuse anything
    flag_reuse = false;
    }

  // Set reuse information
  if(flag_reuse)
    {
    // Reuse the fixed components
    preFilter->SetFlagGenerateFixedComponents(false);
    ncomp_ignore = preFilter->GetNumberOfFixedOnlyOutputComponents();
    }
  else
    {
    // Tell the filter to compute the fixed components
    preFilter->SetFlagGenerateFixedComponents(true);
    }

  // Execute the filter
  preFilter->Update();

  // Get the output image
  InputImageType *img_pre = preFilter->GetOutput();

#ifdef DUMP_NCC
  typename itk::ImageFileWriter<InputImageType>::Pointer pwriter = itk::ImageFileWriter<InputImageType>::New();
  pwriter->SetInput(img_pre);
  pwriter->SetFileName("nccpre.nii.gz");
  pwriter->Update();
#endif

  // It is feasible that the working image has been allocated for more components that
  // are currently used. In this case, we skip those components at the end
  int n_overalloc_comp = m_WorkingImage->GetNumberOfComponentsPerPixel() - ncomp;

  // Currently, we have all the stuff we need to compute the metric in the working
  // image. Next, we run the fast sum computation to give us the local average of
  // intensities, products, gradients in the working image
  typename InputImageType::Pointer img_accum =
      AccumulateNeighborhoodSumsInPlace(img_pre, m_Radius, ncomp_ignore, n_overalloc_comp);

#ifdef DUMP_NCC
  typename itk::ImageFileWriter<InputImageType>::Pointer pwriter = itk::ImageFileWriter<InputImageType>::New();
  pwriter->SetInput(img_accum);
  pwriter->SetFileName("nccaccum.nii.gz");
  pwriter->Update();
#endif

  // At this point, the working image will hold the proper neighborhood sums (I, J, I^2, J^2, IJ, etc).
  // The last step is to use this information to compute the gradients. This is done by the threaded
  // portion of the filter.
}


template <class TMetricTraits>
void
MultiComponentNCCImageMetric<TMetricTraits>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
  int nc_img = this->GetFixedImage()->GetNumberOfComponentsPerPixel();
  int nc = m_WorkingImage->GetNumberOfComponentsPerPixel();
  int line_len = outputRegionForThread.GetSize()[0];

  // Our thread data
  typename Superclass::ThreadData &td = this->m_ThreadData[threadId];

  // Set up an iterator for the working image
  typedef itk::ImageLinearConstIteratorWithIndex<InputImageType> InputIteratorTypeBase;
  typedef IteratorExtender<InputIteratorTypeBase> InputIteratorType;
  InputIteratorType it(m_WorkingImage, outputRegionForThread);

  // Loop over the lines
  for (; !it.IsAtEnd(); it.NextLine())
    {
    // Get the pointer to the input line
    long offset_in_pixels = it.GetPosition() - m_WorkingImage->GetBufferPointer();

    // Pointer to the input pixel data for this line
    const InputComponentType *p_input = m_WorkingImage->GetBufferPointer() + nc * offset_in_pixels;

    // Pointer to the metric data for this line
    MetricPixelType *p_metric = this->GetMetricOutput()->GetBufferPointer() + offset_in_pixels;

    // The gradient output is optional
    GradientPixelType *p_grad_metric = (this->m_ComputeGradient && !this->m_ComputeAffine)
                                       ? this->GetDeformationGradientOutput()->GetBufferPointer() + offset_in_pixels
                                       : NULL;

    // Get the fixed mask like
    typename MaskImageType::PixelType *fixed_mask_line =
        this->GetFixedMaskImage()
        ? this->GetFixedMaskImage()->GetBufferPointer() + offset_in_pixels
        : NULL;

    // Case 1 - dense gradient field requested
    if(!this->m_ComputeAffine)
      {
      if(this->m_ComputeGradient)
        {
        // Loop over the pixels in the line
        for(int i = 0; i < line_len; ++i)
          {
          // Clear the metric and the gradient
          *p_metric = itk::NumericTraits<MetricPixelType>::ZeroValue();
          *p_grad_metric = itk::NumericTraits<GradientPixelType>::ZeroValue();

          if(!fixed_mask_line || fixed_mask_line[i] > 0.5)
            {
            p_input = MultiImageNNCPostComputeFunction(p_input, p_input + nc, nc_img, this->m_Weights.data_block(),
                                                       p_metric, p_grad_metric++, ImageDimension);
            }
          else
            {
            p_grad_metric++;
            p_input+=nc;
            }

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
          if(!fixed_mask_line || fixed_mask_line[i] > 0.5)
            {
            p_input = MultiImageNNCPostComputeFunction(p_input, p_input + nc, nc_img, this->m_Weights.data_block(),
                                                       p_metric, (GradientPixelType *)(NULL), ImageDimension);
            }
          else
            {
            p_input+=nc;
            }

          // Accumulate the total metric
          td.metric += *p_metric;
          td.mask += 1.0;
          }
        }
      }
    // Computing affine
    else
      {
      double *p_grad = this->m_ComputeGradient ? td.gradient.data_block() : NULL;

      // Is there a mask?
      typename MaskImageType::PixelType *mask_line = NULL;
      if(this->GetFixedMaskImage())
        mask_line = this->GetFixedMaskImage()->GetBufferPointer() + offset_in_pixels;

      // Loop over the pixels in the line
      for(int i = 0; i < line_len; ++i)
        {
        // Use mask
        if(!fixed_mask_line || fixed_mask_line[i] > 0.5)
          {
          // Clear the metric and the gradient
          *p_metric = itk::NumericTraits<MetricPixelType>::Zero;

          // Apply the post computation
          p_input = MultiImageNNCPostComputeAffineGradientFunction(
                      p_input, p_input + nc, nc_img, this->m_Weights.data_block(),
                      p_metric, p_grad, ImageDimension);

          // Accumulate the total metric
          td.metric += *p_metric;
          td.mask += 1.0;
          }
        else
          {
          p_input+=nc;
          }
        }
      }
    }
}


#endif // __MultiComponentNCCImageMetric_txx

