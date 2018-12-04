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
#ifndef LieBracketFilter_hxx
#define LieBracketFilter_hxx
#include "LieBracketFilter.h"

#include "itkImageLinearIteratorWithIndex.h"
#include "itkProgressReporter.h"

template< typename TInputImage, typename TOutputImage>
LieBracketFilter< TInputImage, TOutputImage>
::LieBracketFilter()
{
}

template< typename TInputImage, typename TOutputImage>
LieBracketFilter< TInputImage, TOutputImage>
::~LieBracketFilter()
{
}

template< typename TInputImage, typename TOutputImage>
void
LieBracketFilter< TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
  // Call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // For the images U and V, we need to expand their region by radius of one
  for(int i = 0; i < 2; i++)
    {
    InputImageType *input = (i == 0) ? this->GetFieldU() : this->GetFieldV(); 
    typename InputImageType::RegionType inputRR = input->GetRequestedRegion();
    inputRR.PadByRadius(1);
    
    if ( inputRR.Crop( input->GetLargestPossibleRegion() ) )
      {
      input->SetRequestedRegion(inputRR);
      }
    else
      {
      // Couldn't crop the region (requested region is outside the largest
      // possible region).  Throw an exception.

      // store what we tried to request (prior to trying to crop)
      input->SetRequestedRegion(inputRR);

      // build an exception
      itk::InvalidRequestedRegionError e(__FILE__, __LINE__);
      e.SetLocation(ITK_LOCATION);
      e.SetDescription("Requested region is (at least partially) outside the largest possible region.");
      e.SetDataObject(input);
      throw e;
      }
    }
}

template< typename TInputImage, typename TOutputImage>
void
LieBracketFilter< TInputImage, TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread,
                       itk::ThreadIdType threadId)
{
  // Get the input and output
  OutputImageType *out = this->GetOutput();
  const InputImageType *u = this->GetFieldU();
  const InputImageType *v = this->GetFieldV();
  const InputImageType *x = this->GetFieldX();

  // Define the zero vector
  InputPixelType zero_vec(0.0);

  // As the first pass, initialize the output to x if x was passed in, zero otherwise
  if(x)
    {
    itk::ImageRegionIterator<OutputImageType> it0_out(out, outputRegionForThread);
    itk::ImageRegionConstIterator<InputImageType> it0_in(x, outputRegionForThread);
    for(; !it0_out.IsAtEnd(); ++it0_in, ++it0_out)
      it0_out.Value() = it0_in.Value();
    }
  else
    {
    itk::ImageRegionIterator<OutputImageType> it0_out(out, outputRegionForThread);
    for(; !it0_out.IsAtEnd(); ++it0_out)
      it0_out.Value() = zero_vec;
    }

  // The scaling factor
  double scale = 0.5;

  // Iterate over the x, y and z derivatives
  for(unsigned int i = 0; i < InputImageDimension; ++i)
    {
    // Create an iterator over the output region that traverses in lines
    typedef itk::ImageLinearIteratorWithIndex<OutputImageType> IterBase;
    typedef IteratorExtender<IterBase> IterType;
    IterType it(out, outputRegionForThread); it.SetDirection(i);

    // Get the stride of the u and v images
    typename InputImageType::OffsetValueType stride_u = u->GetOffsetTable()[i];
    typename InputImageType::OffsetValueType stride_v = v->GetOffsetTable()[i];

    // The length of the line in this dimension
    int line_len = outputRegionForThread.GetSize(i);

    // Test whether there is data to the left and to the right (or we should use zeros)
    typename InputImageType::IndexType test_idx = it.GetIndex(); 
    test_idx[i] = it.GetIndex()[i] - 1;
    bool have_u_left = u->GetBufferedRegion().IsInside(test_idx);
    bool have_v_left = v->GetBufferedRegion().IsInside(test_idx);
    test_idx[i] = it.GetIndex()[i] + line_len;
    bool have_u_right = u->GetBufferedRegion().IsInside(test_idx);
    bool have_v_right = v->GetBufferedRegion().IsInside(test_idx);

    // Iterate over the lines of the output
    for(; !it.IsAtEnd(); it.NextLine())
      {
      // Get the pointers into the u and v images
      const InputPixelType *ptr_u = u->GetBufferPointer() + u->ComputeOffset(it.GetIndex());
      const InputPixelType *ptr_v = v->GetBufferPointer() + v->ComputeOffset(it.GetIndex());

      // Get the output pointer
      OutputPixelType *ptr_out = out->GetBufferPointer() + out->ComputeOffset(it.GetIndex());

      // These variables store the u and v vectors before and after the current point
      const InputPixelType *u_left = (have_u_left) ? (ptr_u - stride_u) : &zero_vec;
      const InputPixelType *v_left = (have_v_left) ? (ptr_v - stride_v) : &zero_vec;
      const InputPixelType *u_right, *v_right;

      // Iterate over the current line, except the last pixel
      for(int j = 0; j < line_len - 1; j++)
        {
        // Set the left and right pointers
        u_right = ptr_u + stride_u; v_right = ptr_v + stride_v;

        // Compute the directional derivative of u and v
        for(uint k = 0; k < InputImageDimension; k++)
          ptr_out[k] += scale * (
            ((*u_right)[k] - (*u_left)[k]) * (*ptr_v)[i] - 
            ((*v_right)[k] - (*v_left)[k]) * (*ptr_u)[i]);

        // Update the pointers
        u_left = ptr_u; v_left = ptr_v;
        ptr_u = u_right; ptr_v = v_right;
        }

      // Last pixel - do the right side
      u_right = have_u_right ? (ptr_u + stride_u) : &zero_vec;
      v_right = have_v_right ? (ptr_v + stride_v) : &zero_vec;

      // Compute the directional derivative of u and v
      for(uint k = 0; k < InputImageDimension; k++)
        ptr_out[k] += scale * (
          ((*u_right)[k] - (*u_left)[k]) * (*ptr_v)[i] - 
          ((*v_right)[k] - (*v_left)[k]) * (*ptr_u)[i]);
      }
    }
}

#endif
