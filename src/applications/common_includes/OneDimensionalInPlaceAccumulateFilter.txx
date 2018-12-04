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
#include "OneDimensionalInPlaceAccumulateFilter.h"
#include <itkImageLinearIteratorWithIndex.h>
#include "ImageRegionConstIteratorWithIndexOverride.h"


template <class TInputImage>
OneDimensionalInPlaceAccumulateFilter<TInputImage>
::OneDimensionalInPlaceAccumulateFilter()
{
  m_Radius = 0;
  m_Dimension = 0;
  m_ComponentOffsetFront = m_ComponentOffsetBack = 0;
  m_Splitter = SplitterType::New();
  this->InPlaceOn();
}

template <class TInputImage>
const itk::ImageRegionSplitterBase *
OneDimensionalInPlaceAccumulateFilter<TInputImage>
::GetImageRegionSplitter(void) const
{
  m_Splitter->SetDirection(m_Dimension);
  return m_Splitter;
}


template <class TInputImage>
void
OneDimensionalInPlaceAccumulateFilter<TInputImage>
::SetComponentRange(int num_ignored_at_start, int num_ignored_at_end)
{
  m_ComponentOffsetFront = num_ignored_at_start;
  m_ComponentOffsetBack = num_ignored_at_end;
  this->Modified();
}

/**
 * This worker class is defined to allow partial specialization of the ThreadedGenerateData
 * based on the pixel type (float/double)
 */
template <class TPixel, class TInputImage>
class OneDimensionalInPlaceAccumulateFilterWorker
{
public:
  typedef OneDimensionalInPlaceAccumulateFilter<TInputImage> FilterType;
  typedef typename FilterType::OutputImageRegionType OutputImageRegionType;
  typedef typename FilterType::InputImageType InputImageType;
  static void ThreadedGenerateData(FilterType *filter,
                                   const OutputImageRegionType & outputRegionForThread,
                                   itk::ThreadIdType threadId);
};


template <class TPixel, class TInputImage>
void
OneDimensionalInPlaceAccumulateFilterWorker<TPixel, TInputImage>
::ThreadedGenerateData(FilterType *filter,
                       const OutputImageRegionType & outputRegionForThread,
                       itk::ThreadIdType threadId)
{
  // Get filter parameters
  int radius = filter->GetRadius();
  int dimension = filter->GetDimension();

  // Get the image
  InputImageType *image = const_cast<InputImageType *>(filter->GetInput());

  // Set up the iterator that will go through all the lines in the
  // output region. We assume that the lines span the whole length of
  // the input, i.e., the threading direction does not interfere
  typedef itk::ImageLinearIteratorWithIndex<TInputImage> IteratorBaseType;
  typedef IteratorExtenderWithOffset<IteratorBaseType> IteratorType;

  // This is the line iterator, although for even greater speed we operate
  // directly on pointers, so we only use it's NextLine functionality()
  IteratorType itLine(image, outputRegionForThread);
  itLine.SetDirection(dimension);

  // Get the number of components
  int nc = image->GetNumberOfComponentsPerPixel();

  // Get the first and last component for accumulation - these are optionally
  // specified by the user
  int c_first = filter->GetComponentOffsetFront(),
      c_last = (nc - 1) - filter->GetComponentOffsetBack();
  int n_skipped = filter->GetComponentOffsetFront() + filter->GetComponentOffsetBack();

  // Get the offset corresponding to a move along the line for this iterator
  typename IteratorType::OffsetValueType jump = itLine.GetOffset(dimension) * nc;

  // Length of the line being traversed (in whole pixels, then in components)
  int line_length = outputRegionForThread.GetSize(dimension),
      line_length_comp = line_length * nc;

  // Width of the kernel (in whole pixels, then in components)
  int kernel_width = 2 * radius + 1;

  // Allocate an array of the length of the line in components
  TPixel *line = new TPixel[line_length_comp];
  // double *line = new double[line_length_comp];
  // double *sum = new double[nc], *sum_end = sum + nc, *p_sum;

  // Allocate an array to hold the current running sum
  // OutputImageComponentType *sum = new OutputImageComponentType[nc], *sum_end = sum + nc, *p_sum;
  TPixel *sum = new TPixel[nc];

  // Pointers into the sum array for the included components
  TPixel *sum_start = sum + c_first, *sum_end = sum + c_last + 1;

  // Two versions of the code - I thought that maybe the second version (further down) would be
  // more optimized by the compiler, but if anything, I see an opposite effect (although tiny)

#ifdef _ACCUM_ITER_CODE_

  TPixel *p_sum;

  // Start iterating over lines
  for(itLine.GoToBegin(); !itLine.IsAtEnd(); itLine.NextLine())
    {
    int i;

    // Initialize the sum to zero
    for(p_sum = sum_start; p_sum < sum_end; p_sum++)
      *p_sum  = itk::NumericTraits<TPixel>::Zero;

    // Pointer to the current position in the line
    TPixel *p_line = line + c_first, *p_tail = p_line;

    // Pointer to the beginning of the scan line
    long offset_in_pixels = itLine.GetPosition() - image->GetBufferPointer();
    long offset_in_comp = offset_in_pixels * nc;
    const TPixel *p_scan_pixel = image->GetBufferPointer() + offset_in_comp + c_first, *p_scan;

    // Pointer used for writing, it will trail the scan pointer
    TPixel *p_write_pixel = const_cast<TPixel *>(p_scan_pixel), *p_write;

    // Compute the initial sum
    for(i = 0; i < radius; i++)
      {
      for(p_scan = p_scan_pixel, p_sum = sum_start;
          p_sum < sum_end;
          p_sum++, p_line++, p_scan++)
        {
        *p_sum += *p_line = *p_scan;
        }

      p_scan_pixel += jump;
      p_line += n_skipped;
      }

    // For the next Radius + 1 values, add to the sum and write
    for(; i < kernel_width; i++)
      {
      for(p_scan = p_scan_pixel, p_write = p_write_pixel, p_sum = sum_start;
          p_sum < sum_end;
          p_sum++, p_line++, p_scan++, p_write++)
        {
        *p_line = *p_scan;
        *p_sum += *p_line;
        *p_write = *p_sum;
        }

      p_scan_pixel += jump;
      p_write_pixel += jump;
      p_line += n_skipped;
      }

    // Continue until we hit the end of the scanline
    for(; i < line_length; i++)
      {
      for(p_scan = p_scan_pixel, p_write = p_write_pixel, p_sum = sum_start;
          p_sum < sum_end;
          p_sum++, p_line++, p_scan++, p_write++, p_tail++)
        {
        *p_line = *p_scan;
        *p_sum += *p_line - *p_tail;
        *p_write = *p_sum;
        }

      p_scan_pixel += jump;
      p_write_pixel += jump;
      p_line += n_skipped;
      p_tail += n_skipped;
      }

    // Fill out the last bit
    for(; i < line_length + radius; i++)
      {
      for(p_write = p_write_pixel, p_sum = sum_start;
          p_sum < sum_end;
          p_sum++, p_write++, p_tail++)
        {
        *p_sum -= *p_tail;
        *p_write = *p_sum;
        }

      p_write_pixel += jump;
      p_tail += n_skipped;
      }
    }

#else

  // Start iterating over lines
  for(itLine.GoToBegin(); !itLine.IsAtEnd(); itLine.NextLine())
    {

    int i, k, m;

    // Initialize the sum to zero
    for(int k = c_first; k <= c_last; k++)
      sum[k] = itk::NumericTraits<TPixel>::Zero;

    // Pointer to the current position in the line
    TPixel *p_line = line, *p_tail = p_line;

    // Pointer to the beginning of the scan line
    long offset_in_pixels = itLine.GetPosition() - image->GetBufferPointer();
    long offset_in_comp = offset_in_pixels * nc;

    // Where we are scanning from
    const TPixel *p_scan_pixel = image->GetBufferPointer() + offset_in_comp;

    // Pointer used for writing, it will trail the scan pointer
    TPixel *p_write_pixel = const_cast<TPixel *>(p_scan_pixel);

    // Compute the initial sum
    for(i = 0, m = 0; i < radius; i++)
      {
      for(k = c_first; k <= c_last; k++)
        {
        sum[k] += p_line[k] = p_scan_pixel[k];
        }
      p_scan_pixel += jump;
      p_line += nc;
      }

    // For the next Radius + 1 values, add to the sum and write
    for(; i < kernel_width; i++)
      {
      for(k = c_first; k <= c_last; k++)
        {
        p_write_pixel[k] = (sum[k] += p_line[k] = p_scan_pixel[k]);
        }

      p_scan_pixel += jump;
      p_write_pixel += jump;
      p_line += nc;
      }

    // Continue until we hit the end of the scanline
    for(; i < line_length; i++)
      {
      for(k = c_first; k <= c_last; k++)
        {
        p_write_pixel[k] = (sum[k] += (p_line[k] = p_scan_pixel[k]) - p_tail[k]);
        }

      p_scan_pixel += jump;
      p_write_pixel += jump;
      p_line += nc;
      p_tail += nc;
      }

    // Fill out the last bit
    for(; i < line_length + radius; i++)
      {
      for(k = c_first; k <= c_last; k++)
        {
        p_write_pixel[k] = (sum[k] -= p_tail[k]);
        }

      p_write_pixel += jump;
      p_tail += nc;
      }
    }
#endif

  delete[] sum;
  delete[] line;
}

// ### PY, 05/18/2016: checked again, this SSE code is not causing any differences in float/double
//                     processing, safe to keep as is!
#define _NCC_SSE_
#ifdef _NCC_SSE_

#include <xmmintrin.h>

#ifdef WIN32
inline void allocate_aligned(int elements, float ** pointer)
{
  *pointer = (float *)_aligned_malloc(elements * sizeof(float), 16);
  if (*pointer == NULL)
  {
    std::cerr << "_aligned_malloc returned NULL input " << elements * sizeof(float) << std::endl;
    throw std::string("_aligned_malloc allocation error");
  }
}
#else
inline void allocate_aligned(int elements, float ** pointer)
{
  int rc = posix_memalign( (void **) pointer, 16, elements * sizeof(float));
  if(rc != 0)
  {
    std::cerr << "posix_memalign return value " << rc << " input " << elements * sizeof(float) << std::endl;
    throw std::string("posix_memalign allocation error");
  }
}
#endif

/**
 * A specialization of the threaded generate data method for floating point images that uses
 * SSE intrinsics for faster computation
 */
template <class TInputImage>
class OneDimensionalInPlaceAccumulateFilterWorker<float, TInputImage>
{
public:
  typedef OneDimensionalInPlaceAccumulateFilter<TInputImage> FilterType;
  typedef typename FilterType::OutputImageRegionType OutputImageRegionType;
  typedef typename FilterType::InputImageType InputImageType;
  static void ThreadedGenerateData(FilterType *filter,
                                   const OutputImageRegionType & outputRegionForThread,
                                   itk::ThreadIdType threadId);
};

template <class TInputImage>
void
OneDimensionalInPlaceAccumulateFilterWorker<float, TInputImage>
::ThreadedGenerateData(FilterType *filter,
                       const OutputImageRegionType & outputRegionForThread,
                       itk::ThreadIdType threadId)
{
  // Get filter parameters
  int dimension = filter->GetDimension();
  int radius = filter->GetRadius();
  int skip_front = filter->GetComponentOffsetFront();
  int skip_back = filter->GetComponentOffsetBack();

  // Get the image
  InputImageType *image = const_cast<InputImageType *>(filter->GetInput());

  // Set up the iterator that will go through all the lines in the
  // output region. We assume that the lines span the whole length of
  // the input, i.e., the threading direction does not interfere
  typedef itk::ImageLinearIteratorWithIndex<InputImageType> IteratorBaseType;
  typedef IteratorExtenderWithOffset<IteratorBaseType> IteratorType;

  // This is the line iterator, although for even greater speed we operate
  // directly on pointers, so we only use it's NextLine functionality()
  IteratorType itLine(image, outputRegionForThread);
  itLine.SetDirection(dimension);

  // Get the number of components
  int nc = image->GetNumberOfComponentsPerPixel();

  // Get the first and last component for accumulation - these are optionally
  // specified by the user. The remaining components are left untouched
  int c_first = skip_front, c_last = (nc - 1) - skip_back;
  int n_skipped = skip_front + skip_back;

  // Get the offset corresponding to a move along the line for this iterator
  typename IteratorType::OffsetValueType jump = itLine.GetOffset(dimension) * nc;

  // Length of the line being traversed (in whole pixels, then in components)
  int line_length = outputRegionForThread.GetSize(dimension);

  // Width of the kernel (in whole pixels, then in components)
  int kernel_width = 2 * radius + 1;

  // We want some alignment for SIMD purposes. So we need to make a stride be a factor of 16 bytes
  int nc_used = nc - n_skipped;
  int bytes_per_pixel = sizeof(float) * nc_used;

  // Round up, so it works out to 16 bytes
  int align_stride = 4 * sizeof(float);
  int padded_bytes_per_pixel = (bytes_per_pixel % align_stride) == 0
      ? bytes_per_pixel : align_stride * (1 + bytes_per_pixel / align_stride);

  // Number of chunks of four components per pixel
  int nc_padded = padded_bytes_per_pixel / sizeof(float);

  // The following arrays are allocated temporarily
  float *scanline, *tailline, *sum_align;

  // This is a byte-aligned copy of the pixel column from the image
  allocate_aligned(line_length * nc_padded, &scanline);

  // This is a second aligned copy
  allocate_aligned(line_length * nc_padded, &tailline);

  // End of the scanline
  float *p_scanline_end = scanline + line_length * nc_padded;

  // Aligned sum array - where the sums are computed
  allocate_aligned(nc_padded, &sum_align);

  // Start iterating over lines
  for(itLine.GoToBegin(); !itLine.IsAtEnd(); itLine.NextLine())
    {
    int i, k;

    // Pointer to the beginning of the scan line
    long offset_in_pixels = itLine.GetPosition() - image->GetBufferPointer();
    long offset_in_comp = offset_in_pixels * nc;

    // Get the pointer to first component in first pixel
    const float *p_scan_pixel = image->GetBufferPointer() + offset_in_comp + c_first;

    // Registers
    __m128 m_line, m_tail, m_sum_cur, m_sum_new;

    // Copy the contents of the image into the aligned line
    float *p_copy = scanline;
    const float *p_src = p_scan_pixel;
    for(; p_copy < p_scanline_end; p_copy += nc_padded, p_src += jump)
      {
#ifndef WIN32
      __builtin_prefetch(p_src + 5 * jump, 0, 0);
#endif 

      for (i = 0; i < nc_used; i++)
        p_copy[i] = p_src[i];
      }

    // Make a copy of the scan line
    for(p_src = scanline, p_copy = tailline; p_src < p_scanline_end; p_copy+=4, p_src+=4)
      {
      m_line = _mm_load_ps(p_src);
      _mm_store_ps(p_copy, m_line);
      }

    // Clear the sum array at the beginning
    for(k = 0; k < nc_padded; k++)
      sum_align[k] = 0.0;

    // Pointer to the current position in the line
    float *p_line = scanline, *p_tail = tailline;

    // Pointer used for writing, it will trail the scan pointer
    float *p_write_pixel = scanline;

    // Pointer used for writing, it will trail the scan pointer
    float *p_sum_end = sum_align + nc_padded, *p_sum;

    // Compute the initial sum
    for(i = 0; i < radius; i++)
      {
      for(p_sum = sum_align; p_sum < p_sum_end; p_sum+=4, p_line+=4)
        {
        m_line = _mm_load_ps(p_line);
        m_sum_cur = _mm_load_ps(p_sum);
        m_sum_new = _mm_add_ps(m_sum_cur, m_line);
        _mm_store_ps(p_sum, m_sum_new);
        }
      }

    // For the next Radius + 1 values, add to the sum and write
    for(; i < kernel_width; i++)
      {
      for(p_sum = sum_align; p_sum < p_sum_end; p_sum+=4, p_line+=4, p_write_pixel+=4)
        {
        m_line = _mm_load_ps(p_line);
        m_sum_cur = _mm_load_ps(p_sum);
        m_sum_new = _mm_add_ps(m_sum_cur, m_line);
        _mm_store_ps(p_sum, m_sum_new);
        _mm_store_ps(p_write_pixel, m_sum_new);
        }
      }

    // Continue until we hit the end of the scanline
    for(; i < line_length; i++)
      {
      for(p_sum = sum_align; p_sum < p_sum_end; p_sum+=4, p_line+=4, p_tail+=4, p_write_pixel+=4)
        {
        m_line = _mm_load_ps(p_line);
        m_tail = _mm_load_ps(p_tail);
        m_sum_cur = _mm_load_ps(p_sum);
        m_sum_new = _mm_add_ps(m_sum_cur, _mm_sub_ps(m_line, m_tail));
        _mm_store_ps(p_sum, m_sum_new);
        _mm_store_ps(p_write_pixel, m_sum_new);
        }
      }

    // Fill out the last bit
    for(; i < line_length + radius; i++)
      {
      for(p_sum = sum_align; p_sum < p_sum_end; p_sum+=4, p_tail+=4, p_write_pixel+=4)
        {
        m_tail = _mm_load_ps(p_tail);
        m_sum_cur = _mm_load_ps(p_sum);
        m_sum_new = _mm_sub_ps(m_sum_cur, m_tail);
        _mm_store_ps(p_sum, m_sum_new);
        _mm_store_ps(p_write_pixel, m_sum_new);
        }
      }

    // Copy the accumulated pixels back into the main image
    float *p_copy_back = const_cast<float *>(p_scan_pixel);
    const float *p_src_back = scanline;
    for(; p_src_back < p_scanline_end; p_src_back += nc_padded, p_copy_back += jump)
      {
#ifndef WIN32
      __builtin_prefetch(p_copy_back + 5 * jump, 1, 0);
#endif
      for(i = 0; i < nc_used; i++)
        p_copy_back[i] = p_src_back[i];
      }
    }

  // Free allocated memory
  free(tailline);
  free(scanline);
  free(sum_align);
}

#endif // _NCC_SSE_


/**
 * Default implementaton of the threaded generate data method
 */
template <class TInputImage>
void
OneDimensionalInPlaceAccumulateFilter<TInputImage>
::ThreadedGenerateData(
    const OutputImageRegionType & outputRegionForThread,
    itk::ThreadIdType threadId)
{
  typedef OneDimensionalInPlaceAccumulateFilterWorker<OutputImageComponentType, InputImageType> WorkerType;
  WorkerType::ThreadedGenerateData(this, outputRegionForThread, threadId);
}




template <class TInputImage>
typename TInputImage::Pointer
AccumulateNeighborhoodSumsInPlace(TInputImage *image, const typename TInputImage::SizeType &radius,
                                  int num_ignored_at_start, int num_ignored_at_end)
{
  typedef OneDimensionalInPlaceAccumulateFilter<TInputImage> AccumFilterType;

  typename itk::ImageSource<TInputImage>::Pointer pipeTail;
  for(int dir = 0; dir < TInputImage::ImageDimension; dir++)
    {
    typename AccumFilterType::Pointer accum = AccumFilterType::New();
    accum->SetInput(pipeTail.IsNull() ? image : pipeTail->GetOutput());
    accum->SetDimension(dir);
    accum->SetRadius(radius[dir]);
    accum->SetComponentRange(num_ignored_at_start, num_ignored_at_end);
    pipeTail = accum;

    accum->Update();
    }

  return pipeTail->GetOutput();
}
