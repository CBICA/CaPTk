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
#ifndef __MahalanobisDistanceToTargetWarpMetric_txx
#define __MahalanobisDistanceToTargetWarpMetric_txx

#include "MahalanobisDistanceToTargetWarpMetric.h"

template <class TFloat, unsigned int VDim>
class MahalanobisDistanceFunctor
{
public:
  static TFloat compute(const TFloat *mah_data, const TFloat *x, TFloat *grad);
};

template <class TFloat>
class MahalanobisDistanceFunctor<TFloat, 2>
{
public:
  static TFloat compute(const TFloat *mah_data, const TFloat *x, TFloat *grad)
    {
    // Load the elements of the matrix
    TFloat d_x = x[0] - mah_data[0];
    TFloat d_y = x[1] - mah_data[1];
    TFloat s_xx = mah_data[2];
    TFloat s_xy = mah_data[3];
    TFloat s_yy = mah_data[4];

    // Compute
    TFloat md = d_x * d_x * s_xx + d_y * d_y * s_yy + 2 * d_x * d_y * s_xy;

    // Compute gradient
    if(grad)
      {
      grad[0] = - (d_x * s_xx + d_y * s_xy);
      grad[1] = - (d_x * s_xy + d_y * s_yy);
      }

    // Done!
    return md;
    }
};

template <class TFloat>
class MahalanobisDistanceFunctor<TFloat, 3>
{
public:
  static TFloat compute(const TFloat *mah_data, const TFloat *x, TFloat *grad)
    {
    // Load the elements of the matrix
    TFloat d_x = x[0] - mah_data[0];
    TFloat d_y = x[1] - mah_data[1];
    TFloat d_z = x[2] - mah_data[2];
    TFloat s_xx = mah_data[3];
    TFloat s_xy = mah_data[4];
    TFloat s_xz = mah_data[5];
    TFloat s_yy = mah_data[6];
    TFloat s_yz = mah_data[7];
    TFloat s_zz = mah_data[8];

    // Compute
    TFloat md = d_x * d_x * s_xx + d_y * d_y * s_yy + d_z * d_z * s_zz + 
      2 * (d_x * d_y * s_xy + d_x * d_z * s_xz + d_y * d_z * s_yz);

    // Compute gradient
    if(grad)
      {
      grad[0] = -(d_x * s_xx + d_y * s_xy + d_z * s_xz);
      grad[1] = -(d_x * s_xy + d_y * s_yy + d_z * s_yz);
      grad[2] = -(d_x * s_xz + d_y * s_yz + d_z * s_zz);
      }

    // Done!
    return md;
    }
};

template <class TFloat>
class MahalanobisDistanceFunctor<TFloat, 4>
{
public:
  static TFloat compute(const TFloat *mah_data, const TFloat *x, TFloat *grad)
    {
    // Load the elements of the matrix
    TFloat d_x = x[0] - mah_data[0];
    TFloat d_y = x[1] - mah_data[1];
    TFloat d_z = x[2] - mah_data[2];
    TFloat d_w = x[3] - mah_data[3];
    TFloat s_xx = mah_data[4];
    TFloat s_xy = mah_data[5];
    TFloat s_xz = mah_data[6];
    TFloat s_xw = mah_data[7];
    TFloat s_yy = mah_data[8];
    TFloat s_yz = mah_data[9];
    TFloat s_yw = mah_data[10];
    TFloat s_zz = mah_data[11];
    TFloat s_zw = mah_data[12];
    TFloat s_ww = mah_data[13];

    // Compute
    TFloat md = d_x * d_x * s_xx + d_y * d_y * s_yy + d_z * d_z * s_zz + d_w * s_ww +
      2 * (d_x * d_y * s_xy + d_x * d_z * s_xz + d_x * d_w * s_xw + 
           d_y * d_z * s_yz + d_y * d_w * s_yw +
           d_z * d_w * s_zw);

    // Compute gradient
    if(grad)
      {
      grad[0] = -(d_x * s_xx + d_y * s_xy + d_z * s_xz + d_w * s_xw);
      grad[1] = -(d_x * s_xy + d_y * s_yy + d_z * s_yz + d_w * s_yw);
      grad[2] = -(d_x * s_xz + d_y * s_yz + d_z * s_zz + d_w * s_zw);
      grad[3] = -(d_x * s_xw + d_y * s_yw + d_z * s_zw + d_w * s_ww);
      }

    // Done!
    return md;
    }
};



template <class TMetricTraits>
void
MahalanobisDistanceToTargetWarpMetric<TMetricTraits>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
  // Get the fixed image (means and inverse covariances)
  InputImageType *stats = this->GetFixedImage();
  int nc = stats->GetNumberOfComponentsPerPixel();

  // Get the deformation field
  DeformationFieldType *phi = this->GetDeformationField();

  // Create an iterator specialized for going through metrics
  typedef MultiComponentMetricWorker<TMetricTraits, MetricImageType> InterpType;
  InterpType iter(this, this->GetMetricOutput(), outputRegionForThread);

  // The functor that does the actual computation
  typedef MahalanobisDistanceFunctor<typename TMetricTraits::RealType, ImageDimension> Functor;

  // Get the per-component metric array
  typename Superclass::ThreadData &td = this->m_ThreadData[threadId];

  // Loop over the lines
  for (; !iter.IsAtEnd(); iter.NextLine())
    {
    // If we are computing deforamble gradient, get a pointer to the gradient line

    // Temporary gradient pixel
    GradientPixelType grad_metric;

    // Are we computing affine?
    if(this->m_ComputeAffine && this->m_ComputeGradient)
      {
      GradientPixelType grad_metric;
      for(; !iter.IsAtEndOfLine(); ++iter)
        {
        if(iter.CheckFixedMask())
          {
          RealType m = Functor::compute(iter.GetFixedLine(), iter.GetSamplePos().data_block(), 
            grad_metric.GetVnlVector().data_block());
          *iter.GetOutputLine() = m;
          td.metric += m;
          td.mask += 1.0;

          for(int i = 0, q = 0; i < ImageDimension; i++)
            {
            td.gradient[q++] += grad_metric[i];
            for(int j = 0; j < ImageDimension; j++)
              td.gradient[q++] += grad_metric[i] * iter.GetIndex()[j];
            }
          }
        }
      }
    else if(this->m_ComputeGradient)
      {
      // Pointer to the gradient
      GradientPixelType *grad_line =
        this->GetDeformationGradientOutput()->GetBufferPointer() + iter.GetOffsetInPixels();

        for(; !iter.IsAtEndOfLine(); ++iter, ++grad_line)
          {
          if(iter.CheckFixedMask())
            {
            RealType m = Functor::compute(
              iter.GetFixedLine(), 
              iter.GetDisplacement()->GetDataPointer(), 
              grad_line->GetDataPointer());
            *iter.GetOutputLine() = m;
            td.metric += m;
            td.mask += 1.0;
            }
          }
      }
    else
      {
      // No affine, no gradient
      for(; !iter.IsAtEndOfLine(); ++iter)
        {
        if(iter.CheckFixedMask())
          {
          RealType m = Functor::compute(
            iter.GetFixedLine(), iter.GetDisplacement()->GetDataPointer(), NULL);
          *iter.GetOutputLine() = m;
          td.metric += m;
          td.mask += 1.0;
          }
        }
      }
    } // iteration over lines
}




#endif // __MahalanobisDistanceToTargetWarpMetric_txx

