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
#include "GreedyAPI.h"

#include <iostream>
#include <sstream>
#include <cstdio>
#include <vector>
#include <string>
#include <algorithm>

#include "lddmm_common.h"
#include "lddmm_data.h"

#include <itkImageFileReader.h>
#include <itkAffineTransform.h>
#include <itkTransformFactory.h>
#include <itkTimeProbe.h>
#include <itkImageFileWriter.h>

#include "MultiImageRegistrationHelper.h"
#include "FastWarpCompositeImageFilter.h"
#include "MultiComponentImageMetricBase.h"

#include <vnl/algo/vnl_powell.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/vnl_trace.h>
#include <vnl/vnl_numeric_traits.h>

// Little helper functions
template <unsigned int VDim> class array_caster
{
public:
  template <class T> static itk::Size<VDim> to_itkSize(const T &t)
  {
    itk::Size<VDim> sz;
    for(uint i = 0; i < VDim; i++)
      sz[i] = t[i];
    return sz;
  }
};

template <class TITKMatrix, class TVNLMatrix>
void vnl_matrix_to_itk_matrix(
    const TVNLMatrix &vmat,
    TITKMatrix &imat)
{
  for(uint r = 0; r < TITKMatrix::RowDimensions; r++)
    for(uint c = 0; c < TITKMatrix::ColumnDimensions; c++)
      imat(r,c) = static_cast<typename TITKMatrix::ValueType>(vmat(r,c));
}

template <class TITKVector, class TVNLVector>
void vnl_vector_to_itk_vector(
    const TVNLVector &vvec,
    TITKVector &ivec)
{
  for(uint r = 0; r < TITKVector::Dimension; r++)
    ivec[r] = static_cast<typename TITKVector::ValueType>(vvec(r));
}

template <class TITKMatrix, class TVNL>
void itk_matrix_to_vnl_matrix(
    const TITKMatrix &imat,
    vnl_matrix_fixed<TVNL,TITKMatrix::RowDimensions,TITKMatrix::ColumnDimensions>  &vmat)
{
  for(uint r = 0; r < TITKMatrix::RowDimensions; r++)
    for(uint c = 0; c < TITKMatrix::ColumnDimensions; c++)
      vmat(r,c) = static_cast<TVNL>(imat(r,c));
}

template <class TITKMatrix, class TVNL>
void itk_matrix_to_vnl_matrix(
    const TITKMatrix &imat,
    vnl_matrix<TVNL>  &vmat)
{
  vmat.set_size(TITKMatrix::RowDimensions,TITKMatrix::ColumnDimensions);
  for(uint r = 0; r < TITKMatrix::RowDimensions; r++)
    for(uint c = 0; c < TITKMatrix::ColumnDimensions; c++)
      vmat(r,c) = static_cast<TVNL>(imat(r,c));
}

template <class TITKVector, class TVNL>
void itk_vector_to_vnl_vector(
    const TITKVector &ivec,
    vnl_vector_fixed<TVNL,TITKVector::Dimension> &vvec)
{
  for(uint r = 0; r < TITKVector::Dimension; r++)
    vvec(r) = static_cast<TVNL>(ivec[r]);
}

template <class TITKVector, class TVNL>
void itk_vector_to_vnl_vector(
    const TITKVector &ivec,
    vnl_vector<TVNL> &vvec)
{
  vvec.set_size(TITKVector::Dimension);
  for(uint r = 0; r < TITKVector::Dimension; r++)
    vvec(r) = static_cast<TVNL>(ivec[r]);
}

// Helper function to map from ITK coordiante space to RAS space
template<unsigned int VDim, class TMat, class TVec>
void
GetVoxelSpaceToNiftiSpaceTransform(itk::ImageBase<VDim> *image,
                                   TMat &A,
                                   TVec &b)
{
  // Generate intermediate terms
  typedef typename TMat::element_type TReal;
  vnl_matrix<double> m_dir, m_ras_matrix;
  vnl_diag_matrix<double> m_scale, m_lps_to_ras;
  vnl_vector<double> v_origin, v_ras_offset;

  // Compute the matrix
  m_dir = image->GetDirection().GetVnlMatrix();
  m_scale.set(image->GetSpacing().GetVnlVector());
  m_lps_to_ras.set(vnl_vector<double>(VDim, 1.0));
  m_lps_to_ras[0] = -1;
  m_lps_to_ras[1] = -1;
  A = m_lps_to_ras * m_dir * m_scale;

  // Compute the vector
  v_origin = image->GetOrigin().GetVnlVector();
  b = m_lps_to_ras * v_origin;
}

// Helper function to get the RAS coordinate of the center of 
// an image
template <unsigned int VDim>
vnl_vector<double>
GetImageCenterinNiftiSpace(itk::ImageBase<VDim> *image)
{
  itk::ImageRegion<VDim> r = image->GetBufferedRegion();
  itk::ContinuousIndex<double, VDim> idx;
  itk::Point<double, VDim> ctr;
  for(uint d = 0; d < VDim; d++)
    idx[d] = r.GetIndex()[d] + r.GetSize()[d] * 0.5;
  image->TransformContinuousIndexToPhysicalPoint(idx, ctr);

  // Map to RAS (flip first two coordinates)
  for(uint d = 0; d < 2 && d < VDim; d++)
    ctr[d] = -ctr[d];

  return ctr.GetVnlVector();
}



template <unsigned int VDim, typename TReal>
GreedyApproach<VDim, TReal>::PureAffineCostFunction
::PureAffineCostFunction(GreedyParameters *param, ParentType *parent, int level, OFHelperType *helper)
  : AbstractAffineCostFunction(VDim * (VDim + 1))
{
  // Store the data
  m_Param = param;
  m_OFHelper = helper;
  m_Level = level;
  m_Parent = parent;

  // Allocate the working images, but do not allocate. We will allocate on demand because
  // these affine cost functions may be created without needing to do any computation
  m_Allocated = false;

  m_Phi = VectorImageType::New();
  m_Phi->CopyInformation(helper->GetReferenceSpace(level));
  m_Phi->SetRegions(helper->GetReferenceSpace(level)->GetBufferedRegion());

  m_GradMetric = VectorImageType::New();
  m_GradMetric->CopyInformation(helper->GetReferenceSpace(level));
  m_GradMetric->SetRegions(helper->GetReferenceSpace(level)->GetBufferedRegion());

  m_GradMask = VectorImageType::New();
  m_GradMask->CopyInformation(helper->GetReferenceSpace(level));
  m_GradMask->SetRegions(helper->GetReferenceSpace(level)->GetBufferedRegion());

  m_Metric = ImageType::New();
  m_Metric->CopyInformation(helper->GetReferenceSpace(level));
  m_Metric->SetRegions(helper->GetReferenceSpace(level)->GetBufferedRegion());

  m_Mask = ImageType::New();
  m_Mask->CopyInformation(helper->GetReferenceSpace(level));
  m_Mask->SetRegions(helper->GetReferenceSpace(level)->GetBufferedRegion());
}


template <unsigned int VDim, typename TReal>
void
GreedyApproach<VDim, TReal>::PureAffineCostFunction
::compute(const vnl_vector<double> &x, double *f, vnl_vector<double> *g)
{
  // Form a matrix/vector from x
  typename LinearTransformType::Pointer tran = LinearTransformType::New();

  // Set the components of the transform
  unflatten_affine_transform(x.data_block(), tran.GetPointer());

  // Allocate the memory if needed
  if(!m_Allocated)
    {
    m_Phi->Allocate();
    m_GradMetric->Allocate();
    m_GradMask->Allocate();
    m_Metric->Allocate();
    m_Mask->Allocate();
    m_Allocated = true;
    }

  // Compute the gradient
  double val = 0.0;
  if(g)
    {
    typename LinearTransformType::Pointer grad = LinearTransformType::New();

    if(m_Param->metric == GreedyParameters::SSD)
      {
      val = m_OFHelper->ComputeAffineMSDMatchAndGradient(
              m_Level, tran, m_Metric, m_Mask, m_GradMetric, m_GradMask, m_Phi, grad);

      flatten_affine_transform(grad.GetPointer(), g->data_block());
      }
    else if(m_Param->metric == GreedyParameters::NCC)
      {

      val = m_OFHelper->ComputeAffineNCCMatchAndGradient(
              m_Level, tran, array_caster<VDim>::to_itkSize(m_Param->metric_radius),
              m_Metric, m_Mask, m_GradMetric, m_GradMask, m_Phi, grad);

      flatten_affine_transform(grad.GetPointer(), g->data_block());

      // NCC should be maximized
      (*g) *= -10000.0;
      val *= -10000.0;
      }
    else if(m_Param->metric == GreedyParameters::MI || m_Param->metric == GreedyParameters::NMI)
      {
      val = m_OFHelper->ComputeAffineMIMatchAndGradient(
              m_Level, m_Param->metric == GreedyParameters::NMI,
              tran, m_Metric, m_Mask, m_GradMetric, m_GradMask, m_Phi, grad);

      flatten_affine_transform(grad.GetPointer(), g->data_block());

      val *= -10000.0;
      (*g) *= -10000.0;

      }
    }
  else
    {
    if(m_Param->metric == GreedyParameters::SSD)
      {
      val = m_OFHelper->ComputeAffineMSDMatchAndGradient(
              m_Level, tran, m_Metric, m_Mask, m_GradMetric, m_GradMask, m_Phi, NULL);
      }
    else if(m_Param->metric == GreedyParameters::NCC)
      {
      val = m_OFHelper->ComputeAffineNCCMatchAndGradient(
              m_Level, tran, array_caster<VDim>::to_itkSize(m_Param->metric_radius)
              , m_Metric, m_Mask, m_GradMetric, m_GradMask, m_Phi, NULL);

      // NCC should be maximized
      val *= -10000.0;
      }
    else if(m_Param->metric == GreedyParameters::MI || m_Param->metric == GreedyParameters::NMI)
      {
      val = m_OFHelper->ComputeAffineMIMatchAndGradient(
              m_Level, m_Param->metric == GreedyParameters::NMI,
              tran, m_Metric, m_Mask, m_GradMetric, m_GradMask, m_Phi, NULL);

      val *= -10000.0;
      }
    }

  // Has the metric improved?
  if(m_Parent->GetMetricLog().size())
    {
    const std::vector<double> &log = m_Parent->GetMetricLog().back();
    if(log.size() == 0 || log.back() > val)
      {
      // Record the metric value
      m_Parent->RecordMetricValue(val);

      // Write out the current iteration transform
      if(m_Param->output_intermediate.length())
        {
        vnl_matrix<double> Q_physical = MapAffineToPhysicalRASSpace(*m_OFHelper, m_Level, tran);
        m_Parent->WriteAffineMatrixViaCache(m_Param->output_intermediate, Q_physical);
        }
      }
    }

  if(f)
    *f = val;
}

template <unsigned int VDim, typename TReal>
vnl_vector<double>
GreedyApproach<VDim, TReal>::PureAffineCostFunction
::GetCoefficients(LinearTransformType *tran)
{
  vnl_vector<double> x_true(this->get_number_of_unknowns());
  flatten_affine_transform(tran, x_true.data_block());
  return x_true;
}

template <unsigned int VDim, typename TReal>
void
GreedyApproach<VDim, TReal>::PureAffineCostFunction
::GetTransform(const vnl_vector<double> &coeff, LinearTransformType *tran)
{
  unflatten_affine_transform(coeff.data_block(), tran);
}

template <unsigned int VDim, typename TReal>
vnl_vector<double>
GreedyApproach<VDim, TReal>::PureAffineCostFunction
::GetOptimalParameterScaling(const itk::Size<VDim> &image_dim)
{
  // Initialize the scaling vector
  vnl_vector<double> scaling(this->get_number_of_unknowns());

  // Set the scaling of the parameters based on image dimensions. This makes it
  // possible to set tolerances in units of voxels. The order of change in the
  // parameters is comparable to the displacement of any point inside the image
  typename LinearTransformType::MatrixType matrix;
  typename LinearTransformType::OffsetType offset;

  for(uint i = 0; i < VDim; i++)
    {
    offset[i] = 1.0;
    for(uint j = 0; j < VDim; j++)
      matrix(i, j) = image_dim[j];
    }

  typename LinearTransformType::Pointer transform = LinearTransformType::New();
  transform->SetMatrix(matrix);
  transform->SetOffset(offset);
  flatten_affine_transform(transform.GetPointer(), scaling.data_block());

  return scaling;
}

/**
 * PHYSICAL SPACE COST FUNCTION - WRAPS AROUND AFFINE
 */
template <unsigned int VDim, typename TReal>
GreedyApproach<VDim, TReal>::PhysicalSpaceAffineCostFunction
::PhysicalSpaceAffineCostFunction(GreedyParameters *param, ParentType *parent, int level, OFHelperType *helper)
  : AbstractAffineCostFunction(VDim * (VDim + 1)), m_PureFunction(param, parent, level, helper)
{
  // The rigid transformation must be rigid in physical space, not in voxel space
  // So in the constructor, we must compute the mappings from the two spaces
  GetVoxelSpaceToNiftiSpaceTransform(helper->GetReferenceSpace(level), Q_fix, b_fix);
  GetVoxelSpaceToNiftiSpaceTransform(helper->GetMovingReferenceSpace(level), Q_mov, b_mov);

  // Compute the inverse transformations
  Q_fix_inv = vnl_matrix_inverse<double>(Q_fix);
  b_fix_inv = - Q_fix_inv * b_fix;

  Q_mov_inv = vnl_matrix_inverse<double>(Q_mov);
  b_mov_inv = - Q_mov_inv * b_mov;

  // Take advantage of the fact that the transformation is linear in A and b to compute
  // the Jacobian of the transformation ahead of time, and "lazily", using finite differences
  int n = VDim * (VDim + 1);
  J_phys_vox.set_size(n, n);
  vnl_vector<double> x_phys(n, 0), x_vox_0(n), x_vox(n);

  // Voxel parameter vector corresponding to zero transform
  this->map_phys_to_vox(x_phys, x_vox_0);

  // Compute each column of the jacobian
  for(int i = 0; i < n; i++)
    {
    x_phys.fill(0);
    x_phys[i] = 1;
    this->map_phys_to_vox(x_phys, x_vox);
    J_phys_vox.set_column(i, x_vox - x_vox_0);
    }


}

template <unsigned int VDim, typename TReal>
void
GreedyApproach<VDim, TReal>::PhysicalSpaceAffineCostFunction
::map_phys_to_vox(const vnl_vector<double> &x_phys, vnl_vector<double> &x_vox)
{
  Mat A_phys;
  Vec b_phys;

  // unflatten the input parameters into A and b
  unflatten_affine_transform(x_phys.data_block(), A_phys, b_phys);

  // convert into voxel-space affine transform
  Mat A_vox = Q_mov_inv * A_phys * Q_fix;
  Vec b_vox = Q_mov_inv * (A_phys * b_fix + b_phys) + b_mov_inv;

  // Flatten back
  x_vox.set_size(m_PureFunction.get_number_of_unknowns());
  flatten_affine_transform(A_vox, b_vox, x_vox.data_block());
}


template <unsigned int VDim, typename TReal>
void
GreedyApproach<VDim, TReal>::PhysicalSpaceAffineCostFunction
::compute(const vnl_vector<double> &x, double *f, vnl_vector<double> *g)
{
  // Map to voxel space
  vnl_vector<double> x_vox(m_PureFunction.get_number_of_unknowns());
  this->map_phys_to_vox(x, x_vox);

  // Do we need the gradient?
  if(g)
    {
    // Compute the function and gradient wrt voxel parameters
    vnl_vector<double> g_vox(m_PureFunction.get_number_of_unknowns());
    m_PureFunction.compute(x_vox, f, &g_vox);

    // Transform voxel-space gradient into physical-space gradient
    *g = J_phys_vox.transpose() * g_vox;
    }
  else
    {
    // Just compute the function
    m_PureFunction.compute(x_vox, f, NULL);
    }
}

template <unsigned int VDim, typename TReal>
vnl_vector<double>
GreedyApproach<VDim, TReal>::PhysicalSpaceAffineCostFunction
::GetOptimalParameterScaling(const itk::Size<VDim> &image_dim)
{
  // TODO: work out scaling for this
  return m_PureFunction.GetOptimalParameterScaling(image_dim);
}

template <unsigned int VDim, typename TReal>
vnl_vector<double>
GreedyApproach<VDim, TReal>::PhysicalSpaceAffineCostFunction
::GetCoefficients(LinearTransformType *tran)
{
  // The input transform is in voxel space, we must return parameters in physical space
  Mat A_vox, A_phys;
  Vec b_vox, b_phys;

  itk_matrix_to_vnl_matrix(tran->GetMatrix(), A_vox);
  itk_vector_to_vnl_vector(tran->GetOffset(), b_vox);

  // convert into physical-space affine transform
  A_phys = Q_mov * A_vox * Q_fix_inv;
  b_phys = Q_mov * (b_vox - b_mov_inv) - A_phys * b_fix;

  // Flatten
  vnl_vector<double> x(m_PureFunction.get_number_of_unknowns());
  flatten_affine_transform(A_phys, b_phys, x.data_block());

  return x;
}

template <unsigned int VDim, typename TReal>
void
GreedyApproach<VDim, TReal>::PhysicalSpaceAffineCostFunction
::GetTransform(const vnl_vector<double> &x, LinearTransformType *tran)
{
  // Get voxel-space tranform corresponding to the parameters x
  vnl_vector<double> x_vox(m_PureFunction.get_number_of_unknowns());
  this->map_phys_to_vox(x, x_vox);

  // Unflatten into a transform
  unflatten_affine_transform(x_vox.data_block(), tran);
}


/**
 * SCALING COST FUNCTION - WRAPS AROUND AFFINE
 */
template <unsigned int VDim, typename TReal>
void
GreedyApproach<VDim, TReal>::ScalingCostFunction
::compute(const vnl_vector<double> &x, double *f, vnl_vector<double> *g)
{
  // Scale the parameters so they are in unscaled units
  vnl_vector<double> x_scaled = element_quotient(x, m_Scaling);

  // Call the wrapped method
  if(g)
    {
    vnl_vector<double> g_scaled(x_scaled.size());
    m_PureFunction->compute(x_scaled, f, &g_scaled);
    *g = element_quotient(g_scaled, m_Scaling);
    }
  else
    {
    m_PureFunction->compute(x_scaled, f, g);
    }
}

// Get the parameters for the specified initial transform
template <unsigned int VDim, typename TReal>
vnl_vector<double>
GreedyApproach<VDim, TReal>::ScalingCostFunction
::GetCoefficients(LinearTransformType *tran)
{
  vnl_vector<double> x_true = m_PureFunction->GetCoefficients(tran);
  return element_product(x_true, m_Scaling);
}

// Get the transform for the specificed coefficients
template <unsigned int VDim, typename TReal>
void
GreedyApproach<VDim, TReal>::ScalingCostFunction
::GetTransform(const vnl_vector<double> &coeff, LinearTransformType *tran)
{
  vnl_vector<double> x_true = element_quotient(coeff, m_Scaling);
  m_PureFunction->GetTransform(x_true, tran);
}



/**
 * RIGID COST FUNCTION - WRAPS AROUND AFFINE
 */
template <unsigned int VDim, typename TReal>
GreedyApproach<VDim, TReal>::RigidCostFunction
::RigidCostFunction(GreedyParameters *param, ParentType *parent, int level, OFHelperType *helper)
  : AbstractAffineCostFunction(VDim * 2), m_AffineFn(param, parent, level, helper)
{
}

template <unsigned int VDim, typename TReal>
void
GreedyApproach<VDim, TReal>::RigidCostFunction
::compute(const vnl_vector<double> &x, double *f, vnl_vector<double> *g)
{
  // Place parameters into q and b
  Vec3 q, b;
  q[0] = x[0]; q[1] = x[1]; q[2] = x[2];
  b[0] = x[3]; b[1] = x[4]; b[2] = x[5];

  // Compute theta
  double theta = q.magnitude();

  // Predefine the rotation matrix
  Mat3 R; R.set_identity();

  // Create the Q matrix
  Mat3 Qmat; Qmat.fill(0.0);
  Qmat(0,1) = -q[2]; Qmat(1,0) =  q[2];
  Qmat(0,2) =  q[1]; Qmat(2,0) = -q[1];
  Qmat(1,2) = -q[0]; Qmat(2,1) =  q[0];

  // Compute the square of the matrix
  Mat3 QQ = vnl_matrix_fixed_mat_mat_mult(Qmat, Qmat);

  // A small epsilon for which a better approximation is R = I + Q
  double eps = 1.0e-4;
  double a1, a2;

  // When theta = 0, rotation is identity
  if(theta > eps)
    {
    // Compute the constant terms in the Rodriguez formula
    a1 = sin(theta) / theta;
    a2 = (1 - cos(theta)) / (theta * theta);

    // Compute the rotation matrix
    R += a1 * Qmat + a2 * QQ;
    }
  else
    {
    R += Qmat;
    }

  // Now we have a rotation and a translation, convert to parameters for the affine function
  vnl_vector<double> x_affine(m_AffineFn.get_number_of_unknowns());
  flatten_affine_transform(R, b, x_affine.data_block());

  // Split depending on whether there is gradient to compute
  if(g)
    {
    // Create a vector to store the affine gradient
    vnl_vector<double> g_affine(m_AffineFn.get_number_of_unknowns());
    m_AffineFn.compute(x_affine, f, &g_affine);

    // Compute the matrices d_Qmat
    Mat3 d_Qmat[3], d_R[3];
    d_Qmat[0].fill(0); d_Qmat[0](1,2) = -1; d_Qmat[0](2,1) =  1;
    d_Qmat[1].fill(0); d_Qmat[1](0,2) =  1; d_Qmat[1](2,0) = -1;
    d_Qmat[2].fill(0); d_Qmat[2](0,1) = -1; d_Qmat[2](1,0) =  1;

    // Compute partial derivatives of R wrt q
    if(theta > eps)
      {
      // Compute the scaling factors in the Rodriguez formula
      double d_a1 = (theta * cos(theta) - sin(theta)) / (theta * theta * theta);
      double d_a2 = (theta * sin(theta) + 2 * cos(theta) - 2) /
                    (theta * theta * theta * theta);

      // Loop over the coordinate and compute the derivative of the rotation matrix wrt x
      for(uint p = 0; p < 3; p++)
        {
        // Compute the gradient of the rotation with respect to q[p]
        d_R[p] = d_a1 * q[p] * Qmat +
                 a1 * d_Qmat[p] +
                 d_a2 * q[p] * vnl_matrix_fixed_mat_mat_mult(Qmat, Qmat)
                 + a2 * (vnl_matrix_fixed_mat_mat_mult(d_Qmat[p], Qmat) +
                         vnl_matrix_fixed_mat_mat_mult(Qmat, d_Qmat[p]));
        }
      }
    else
      {
      for(uint p = 0; p < 3; p++)
        d_R[p] = d_Qmat[p];
      }

    // Create a matrix to hold the jacobian
    vnl_matrix<double> jac(m_AffineFn.get_number_of_unknowns(), 6);
    jac.fill(0.0);

    // Zero vector
    Vec3 zero_vec; zero_vec.fill(0.0);
    Mat3 zero_mat; zero_mat.fill(0.0);

    // Fill out the jacobian
    for(uint p = 0; p < 3; p++)
      {
      // Fill the corresponding column
      vnl_vector<double> jac_col_q(m_AffineFn.get_number_of_unknowns());
      flatten_affine_transform(d_R[p], zero_vec, jac_col_q.data_block());
      jac.set_column(p, jac_col_q);

      // Also set column on the right (wrt translation)
      vnl_vector<double> jac_col_b(m_AffineFn.get_number_of_unknowns());
      Vec3 ep; ep.fill(0.0); ep[p] = 1;
      flatten_affine_transform(zero_mat, ep, jac_col_b.data_block());
      jac.set_column(p+3, jac_col_b);
      }

    // Multiply the gradient by the jacobian
    *g = jac.transpose() * g_affine;
    }
  else
    {
    m_AffineFn.compute(x_affine, f, NULL);
    }
}

template <unsigned int VDim, typename TReal>
vnl_vector<double>
GreedyApproach<VDim, TReal>::RigidCostFunction
::GetOptimalParameterScaling(const itk::Size<VDim> &image_dim)
{
  // Initialize the scaling vector
  vnl_vector<double> scaling(this->get_number_of_unknowns());

  // Scaling is harder for rotations. The rotation parameters are in units of
  // radians. We must figure out how many radians are equivalent to a point in
  // the image moving by a single voxel. That actually works out to be 1/dim.

  // So we take the average of the image dimensions and use that as scaling
  double mean_dim = 0;
  for(uint i = 0; i < VDim; i++)
    mean_dim += image_dim[i] / VDim;
  scaling[0] = scaling[1] = scaling[2] = mean_dim;
  scaling[3] = scaling[4] = scaling[5] = 1.0;

  return scaling;
}

template <unsigned int VDim, typename TReal>
vnl_vector<double>
GreedyApproach<VDim, TReal>::RigidCostFunction
::GetCoefficients(LinearTransformType *tran)
{
  // This affine transform is in voxel space. We must first map it into physical
  vnl_vector<double> x_aff_phys = m_AffineFn.GetCoefficients(tran);
  Mat3 A; Vec3 b;
  unflatten_affine_transform(x_aff_phys.data_block(), A, b);

  // Compute polar decomposition of the affine matrix
  vnl_svd<double> svd(A);
  Mat3 R = svd.U() * svd.V().transpose();
  Vec3 q = this->GetAxisAngle(R);

  // Make result
  vnl_vector<double> x(6);
  x[0] = q[0]; x[1] = q[1]; x[2] = q[2];
  x[3] = b[0]; x[4] = b[1]; x[5] = b[2];

  return x;
}

template <unsigned int VDim, typename TReal>
typename GreedyApproach<VDim, TReal>::RigidCostFunction::Vec3
GreedyApproach<VDim, TReal>::RigidCostFunction
::GetAxisAngle(const Mat3 &R)
{
  double eps = 1e-4;
  double f_thresh = cos(eps);

  // Compute the matrix logarithm of R
  double f = (vnl_trace(R) - 1) / 2;
  Vec3 q;
  if(f >= f_thresh)
    {
    q[0] = R(2,1) - R(1,2);
    q[1] = R(0,2) - R(2,0);
    q[2] = R(1,0) - R(0,1);
    q *= 0.5;
    }
  else
    {
    double theta = acos(f);
    double sin_theta = sqrt(1 - f * f);
    q[0] = R(2,1) - R(1,2);
    q[1] = R(0,2) - R(2,0);
    q[2] = R(1,0) - R(0,1);
    q *= theta / (2 * sin_theta);
    }

  return q;
}



template <unsigned int VDim, typename TReal>
vnl_vector<double>
GreedyApproach<VDim, TReal>::RigidCostFunction
::GetRandomCoeff(const vnl_vector<double> &xInit, vnl_random &randy, double sigma_angle, double sigma_xyz,
                 const Vec3 &C_fixed, const Vec3 &C_moving)
{
  // Generate a random axis of rotation. A triple of Gaussian numbers, normalized to 
  // unit length gives a uniform distribution over the sphere
  Vec3 q_axis;
  for(uint d = 0; d < 3; d++)
    q_axis[d] = randy.normal();
  q_axis.normalize();

  // Generate an angle of rotation from the normal distribution (degrees->radians)
  Vec3 q = q_axis * (randy.normal() * sigma_angle * 0.01745329252);

  // Generate a random rotation using given angles
  Mat3 R = this->GetRotationMatrix(q);

  // Generate a rotation matrix for the initial parameters
  Vec3 qInit;
  for(uint d = 0; d < 3; d++)
    qInit[d] = xInit[d];
  Mat3 R_init = this->GetRotationMatrix(qInit);

  // Combined rotation
  Mat3 R_comb = R * R_init;

  // Take the log map
  Vec3 q_comb = this->GetAxisAngle(R_comb);

  // Generate the offset
  Vec3 b = C_moving - R_comb * C_fixed;

  // Apply random offset
  for(uint d = 0; d < 3; d++)
    b[d] += randy.normal() * sigma_xyz;

  // Generate output vector
  vnl_vector<double> x(6);
  x[0] = q_comb[0];
  x[1] = q_comb[1];
  x[2] = q_comb[2];
  x[3] = b[0];
  x[4] = b[1];
  x[5] = b[2];

  return x;
}

template <unsigned int VDim, typename TReal>
typename GreedyApproach<VDim, TReal>::RigidCostFunction::Mat3
GreedyApproach<VDim, TReal>::RigidCostFunction
::GetRotationMatrix(const Vec3 &q)
{
  // Compute theta
  double theta = q.magnitude();

  // Predefine the rotation matrix
  Mat3 R; R.set_identity();

  // Create the Q matrix
  Mat3 Qmat; Qmat.fill(0.0);
  Qmat(0,1) = -q[2]; Qmat(1,0) =  q[2];
  Qmat(0,2) =  q[1]; Qmat(2,0) = -q[1];
  Qmat(1,2) = -q[0]; Qmat(2,1) =  q[0];

  // Compute the square of the matrix
  Mat3 QQ = vnl_matrix_fixed_mat_mat_mult(Qmat, Qmat);

  // When theta = 0, rotation is identity
  double eps = 1e-4;

  if(theta > eps)
    {
    // Compute the constant terms in the Rodriguez formula
    double a1 = sin(theta) / theta;
    double a2 = (1 - cos(theta)) / (theta * theta);

    // Compute the rotation matrix
    R += a1 * Qmat + a2 * QQ;
    }
  else
    {
    R += Qmat;
    }

  return R;
}

template <unsigned int VDim, typename TReal>
void
GreedyApproach<VDim, TReal>::RigidCostFunction
::GetTransform(const vnl_vector<double> &x, LinearTransformType *tran)
{
  // Place parameters into q and b
  Vec3 q, b;
  q[0] = x[0]; q[1] = x[1]; q[2] = x[2];
  b[0] = x[3]; b[1] = x[4]; b[2] = x[5];

  // Get the rotation matrix
  Mat3 R = this->GetRotationMatrix(q);

  // This gives us the physical space affine matrices. Flatten and map to voxel space
  vnl_vector<double> x_aff_phys(m_AffineFn.get_number_of_unknowns());
  flatten_affine_transform(R, b, x_aff_phys.data_block());
  m_AffineFn.GetTransform(x_aff_phys, tran);
}


/*
template <unsigned int VDim, typename TReal>
void
GreedyApproach<VDim, double>::AffineCostFunction
::compute(const vnl_vector<double> &x, double *f, vnl_vector<double> *g)
{
  // Form a matrix/vector from x
  typename TransformType::Pointer tran = TransformType::New();

  // Set the components of the transform
  unflatten_affine_transform(x.data_block(), tran.GetPointer());

  // Compute the gradient
  double val = 0.0;
  if(g)
    {
    typename TransformType::Pointer grad = TransformType::New();
    val = m_OFHelper->ComputeAffineMatchAndGradient(m_Level, tran, grad);
    flatten_affine_transform(grad.GetPointer(), g->data_block());
    }
  else
    {
    val = m_OFHelper->ComputeAffineMatchAndGradient(m_Level, tran, NULL);
    }

  if(f)
    *f = val;
}
*
*
*/

#include "itkTransformFileReader.h"

template <unsigned int VDim, typename TReal>
vnl_matrix<double>
GreedyApproach<VDim, TReal>
::ReadAffineMatrixViaCache(const TransformSpec &ts)
{
  // Physical (RAS) space transform matrix
  vnl_matrix<double> Qp(VDim+1, VDim+1);  Qp.set_identity();

  // An ITK-style transform - forced to floating point here
  typedef itk::MatrixOffsetTransformBase<double, VDim, VDim> TransformType;
  typename TransformType::Pointer itk_tran;

  // See if a transform is already stored in the cache
  typename ImageCache::const_iterator itCache = m_ImageCache.find(ts.filename);
  if(itCache != m_ImageCache.end())
    {
    TransformType *cached = dynamic_cast<TransformType *>(itCache->second);
    if(!cached)
      throw GreedyException("Cached transform %s cannot be cast to type %s",
                            ts.filename.c_str(), typeid(TransformType).name());

    itk_tran = cached;
    }
  else
    {
    // Open the file and read the first line
    std::ifstream fin(ts.filename.c_str());
    std::string header_line, itk_header = "#Insight Transform File";
    std::getline(fin, header_line);

    if(header_line.substr(0, itk_header.size()) == itk_header)
      {
      fin.close();
      try
        {
        // First we try to load the transform using ITK code
        // This code is from c3d_affine_tool
        typedef itk::AffineTransform<double, VDim> AffTran;
        itk::TransformFactory<TransformType>::RegisterTransform();
        itk::TransformFactory<AffTran>::RegisterTransform();

        itk::TransformFileReader::Pointer fltReader = itk::TransformFileReader::New();
        fltReader->SetFileName(ts.filename.c_str());
        fltReader->Update();

        itk::TransformBase *base = fltReader->GetTransformList()->front();
        itk_tran = dynamic_cast<TransformType *>(base);
        }
      catch(...)
        {
        throw GreedyException("Unable to read ITK transform file %s", ts.filename.c_str());
        }
      }
    else
      {
      // Try reading C3D matrix format
      fin.seekg(0);
      for(size_t i = 0; i < VDim+1; i++)
        for(size_t j = 0; j < VDim+1; j++)
          if(fin.good())
            {
            fin >> Qp[i][j];
            }
      fin.close();
      }
    }

  // At this point we might have read the RAS matrix directly, or an ITK transform
  // if the latter, extract the RAS matrix
  if(itk_tran.IsNotNull())
    {
    for(size_t r = 0; r < VDim; r++)
      {
      for(size_t c = 0; c < VDim; c++)
        {
        Qp(r,c) = itk_tran->GetMatrix()(r,c);
        }
      Qp(r,VDim) = itk_tran->GetOffset()[r];
      }

    // RAS - LPI nonsense
    if(VDim == 3)
      {
      Qp(2,0) *= -1; Qp(2,1) *= -1;
      Qp(0,2) *= -1; Qp(1,2) *= -1;
      Qp(0,3) *= -1; Qp(1,3) *= -1;
      }
    }

  // Compute the exponent
  if(ts.exponent == 1.0)
    {
    return Qp;
    }
  else if(ts.exponent == -1.0)
    {
    return vnl_matrix_inverse<double>(Qp);
    }
  else
    {
    throw GreedyException("Transform exponent values of +1 and -1 are the only ones currently supported");
    }

  return Qp;
}




template <unsigned int VDim, typename TReal>
void
GreedyApproach<VDim, TReal>
::WriteAffineMatrixViaCache(
    const std::string &filename, const vnl_matrix<double> &Qp)
{
  // An ITK-style transform - forced to double point here
  typedef itk::MatrixOffsetTransformBase<double, VDim, VDim> TransformType;

  // See if a transform is already stored in the cache
  typename ImageCache::const_iterator itCache = m_ImageCache.find(filename);
  if(itCache != m_ImageCache.end())
    {
    TransformType *cached = dynamic_cast<TransformType *>(itCache->second);
    if(!cached)
      throw GreedyException("Cached transform %s cannot be cast to type %s",
                            filename.c_str(), typeid(TransformType).name());

    // RAS - LPI nonsense
    vnl_matrix<double> Q = Qp;
    if(VDim == 3)
      {
      Q(2,0) *= -1; Q(2,1) *= -1;
      Q(0,2) *= -1; Q(1,2) *= -1;
      Q(0,3) *= -1; Q(1,3) *= -1;
      }

    typename TransformType::MatrixType matrix;
    typename TransformType::OffsetType offset;

    // We have found the output transform and can use it for assignment
    for(size_t r = 0; r < VDim; r++)
      {
      for(size_t c = 0; c < VDim; c++)
        {
        matrix(r, c) = Q(r, c);
        }
      offset[r] = Q(r, VDim);
      }

    cached->SetMatrix(matrix);
    cached->SetOffset(offset);
    }
  else
    {
    std::ofstream matrixFile;
    matrixFile.open(filename.c_str());
    matrixFile << Qp;
    matrixFile.close();
    }
}



template <unsigned int VDim, typename TReal>
template <class TImage>
itk::SmartPointer<TImage>
GreedyApproach<VDim, TReal>
::ReadImageViaCache(const std::string &filename)
{
  // Check the cache for the presence of the image
  typename ImageCache::const_iterator it = m_ImageCache.find(filename);
  if(it != m_ImageCache.end())
    {
    TImage *image = dynamic_cast<TImage *>(it->second);
    if(!image)
      throw GreedyException("Cached image %s cannot be cast to type %s",
                            filename.c_str(), typeid(TImage).name());
    itk::SmartPointer<TImage> pointer = image;

    return pointer;
    }

  // Read the image using ITK reader
  typedef itk::ImageFileReader<TImage> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  itk::SmartPointer<TImage> pointer = reader->GetOutput();
  return pointer;
}


template <unsigned int VDim, typename TReal>
void GreedyApproach<VDim, TReal>
::ReadImages(GreedyParameters &param, OFHelperType &ofhelper)
{
  // If the parameters include a sequence of transforms, apply it first
  VectorImagePointer moving_pre_warp;

  // Read the input images and stick them into an image array
  for(uint i = 0; i < param.inputs.size(); i++)
    {
    // Read fixed and moving images
    CompositeImagePointer imgFix = ReadImageViaCache<CompositeImageType>(param.inputs[i].fixed);
    CompositeImagePointer imgMov = ReadImageViaCache<CompositeImageType>(param.inputs[i].moving);

    // Read the pre-warps (only once)
    if(param.moving_pre_transforms.size() && moving_pre_warp.IsNull())
      {
      ReadTransformChain(param.moving_pre_transforms, imgFix, moving_pre_warp);
      }

    if(moving_pre_warp.IsNotNull())
      {
      // Create an image to store the warp
      CompositeImagePointer warped_moving;
      LDDMMType::alloc_cimg(warped_moving, imgFix,
                            imgMov->GetNumberOfComponentsPerPixel());

      // Interpolate the moving image using the transform chain
      LDDMMType::interp_cimg(imgMov, moving_pre_warp, warped_moving, false, true);

      // Add the image pair to the helper
      ofhelper.AddImagePair(imgFix, warped_moving, param.inputs[i].weight);
      }
    else
      {
      // Add to the helper object
      ofhelper.AddImagePair(imgFix, imgMov, param.inputs[i].weight);
      }
    }

  // Read the fixed-space mask
  if(param.gradient_mask.size())
    {
    typedef typename OFHelperType::FloatImageType MaskType;
    typename MaskType::Pointer imgMask =
        ReadImageViaCache<MaskType>(param.gradient_mask);

    ofhelper.SetGradientMask(imgMask);
    }

  // Read the moving-space mask
  if(param.moving_mask.size())
    {
    typedef typename OFHelperType::FloatImageType MaskType;
    typename MaskType::Pointer imgMovMask =
        ReadImageViaCache<MaskType>(param.moving_mask);

    if(moving_pre_warp.IsNotNull())
      {
      // Create an image to store the warp
      typename MaskType::Pointer warped_moving_mask;
      LDDMMType::alloc_img(warped_moving_mask, moving_pre_warp);

      // Interpolate the moving image using the transform chain
      LDDMMType::interp_img(imgMovMask, moving_pre_warp, warped_moving_mask, false, true);

      // Add the warped mask to the helper
      ofhelper.SetMovingMask(warped_moving_mask);
      }
    else
      {
      // Add the mask to the helper object
      ofhelper.SetMovingMask(imgMovMask);
      }
    }

  // Generate the optimized composite images. For the NCC metric, we add random noise to
  // the composite images, specified in units of the interquartile intensity range.
  double noise = (param.metric == GreedyParameters::NCC) ? param.ncc_noise_factor : 0.0;

  // Build the composite images
  ofhelper.BuildCompositeImages(noise);

  // If the metric is NCC, then also apply special processing to the gradient masks
  if(param.metric == GreedyParameters::NCC)
    ofhelper.DilateCompositeGradientMasksForNCC(array_caster<VDim>::to_itkSize(param.metric_radius));
}

#include <vnl/algo/vnl_lbfgs.h>

template <unsigned int VDim, typename TReal>
vnl_matrix<double>
GreedyApproach<VDim, TReal>
::MapAffineToPhysicalRASSpace(
    OFHelperType &of_helper, int level,
    LinearTransformType *tran)
{
  // Map the transform to NIFTI units
  vnl_matrix<double> T_fix, T_mov, Q, A;
  vnl_vector<double> s_fix, s_mov, p, b;

  GetVoxelSpaceToNiftiSpaceTransform(of_helper.GetReferenceSpace(level), T_fix, s_fix);
  GetVoxelSpaceToNiftiSpaceTransform(of_helper.GetMovingReferenceSpace(level), T_mov, s_mov);

  itk_matrix_to_vnl_matrix(tran->GetMatrix(), A);
  itk_vector_to_vnl_vector(tran->GetOffset(), b);

  Q = T_mov * A * vnl_matrix_inverse<double>(T_fix);
  p = T_mov * b + s_mov - Q * s_fix;

  vnl_matrix<double> Qp(VDim+1, VDim+1);
  Qp.set_identity();
  for(uint i = 0; i < VDim; i++)
    {
    Qp(i, VDim) = p(i);
    for(uint j = 0; j < VDim; j++)
      Qp(i,j) = Q(i,j);
    }

  return Qp;
}

template <unsigned int VDim, typename TReal>
void
GreedyApproach<VDim, TReal>
::MapPhysicalRASSpaceToAffine(
    OFHelperType &of_helper, int level,
    vnl_matrix<double> &Qp,
    LinearTransformType *tran)
{
  // Map the transform to NIFTI units
  vnl_matrix<double> T_fix, T_mov, Q(VDim, VDim), A;
  vnl_vector<double> s_fix, s_mov, p(VDim), b;

  GetVoxelSpaceToNiftiSpaceTransform(of_helper.GetReferenceSpace(level), T_fix, s_fix);
  GetVoxelSpaceToNiftiSpaceTransform(of_helper.GetMovingReferenceSpace(level), T_mov, s_mov);

  for(uint i = 0; i < VDim; i++)
    {
    p(i) = Qp(i, VDim);
    for(uint j = 0; j < VDim; j++)
      Q(i,j) = Qp(i,j);
    }

  // A = vnl_matrix_inverse<double>(T_mov) * (Q * T_fix);
  // b = vnl_matrix_inverse<double>(T_mov) * (p - s_mov + Q * s_fix);
  A=vnl_svd<double>(T_mov).solve(Q * T_fix);
  b=vnl_svd<double>(T_mov).solve(p - s_mov + Q * s_fix);

  typename LinearTransformType::MatrixType tran_A;
  typename LinearTransformType::OffsetType tran_b;

  vnl_matrix_to_itk_matrix(A, tran_A);
  vnl_vector_to_itk_vector(b, tran_b);

  tran->SetMatrix(tran_A);
  tran->SetOffset(tran_b);
}

template <unsigned int VDim, typename TReal>
void
GreedyApproach<VDim, TReal>
::RecordMetricValue(double val)
{
  if(m_MetricLog.size())
    m_MetricLog.back().push_back(val);
}

/**
 * Find a plane of symmetry in an image
 */
/*
template <unsigned int VDim, typename TReal>
vnl_vector<double>
GreedyApproach<VDim, TReal>
::FindSymmetryPlane(ImageType *image, int N, int n_search_pts)
{
  typedef vnl_vector_fixed<double, 3> Vec3;
  typedef vnl_matrix_fixed<double, 3, 3> Mat3;

  // Loop over direction on a sphere, using the Saff & Kuijlaars algorithm
  // https://perswww.kuleuven.be/~u0017946/publications/Papers97/art97a-Saff-Kuijlaars-MI/Saff-Kuijlaars-MathIntel97.pdf
  double phi = 0.0;
  double spiral_const = 3.6 / sqrt(N);
  for(uint k = 0; k < n_sphere_pts; k++)
    {
    // Height of the k-th point
    double cos_theta = -1 * (2 * k) / (N - 1);
    double sin_theta = sqrt(1 - cos_theta * cos_theta);

    // Phase of the k-th point
    if(k > 0 && k < N-1)
      phi = fmod(phi_last + spiral_const / sin_theta, vnl_math::pi * 2);
    else
      phi = 0.0;

    // We now have the polar coordinates of the points, get cartesian coordinates
    Vec3 q;
    q[0] = sin_theta * cos(phi);
    q[1] = sin_theta * sin(phi);
    q[2] = cos_theta;

    // Now q * (x,y,z) = 0 defines a plane through the origin. We will test whether the image
    // is symmetric across this plane. We first construct the reflection matrix
    Mat3 R;
    R(0,0) =  1 - q[0] * q[0]; R(0,1) = -2 * q[1] * q[0]; R(0,2) = -2 * q[2] * q[0];
    R(1,0) = -2 * q[0] * q[1]; R(1,1) =  1 - q[1] * q[1]; R(1,2) = -2 * q[2] * q[1];
    R(2,0) = -2 * q[0] * q[2]; R(2,1) = -2 * q[1] * q[2]; R(2,2) =  1 - q[2] * q[2];

    // We must find the reasonable range of intercepts to search for. An intercept is reasonable
    // if the plane cuts the image volume in at least a 80/20 ratio (let's say)


    // This is a test axis of rotation. We will now measure the symmetry of the image across this axis
    // To do so, we will construct a series of flips across this direction

    }
}
*/

/**
 * This method performs initial alignment by first searching for a plane of symmetry
 * in each image, and then finding the transform between the planes of symmetry.
 *
 * The goal is to have an almost sure-fire method for brain alignment, yet generic
 * enough to work for other data as well.
 */
/*
template <unsigned int VDim, typename TReal>
int GreedyApproach<VDim, TReal>
::SymmetrySearch(GreedyParameters &param, int level, OFHelperType *of_helper)
{

}
*/


template <unsigned int VDim, typename TReal>
int GreedyApproach<VDim, TReal>
::RunAffine(GreedyParameters &param)
{
  // Create an optical flow helper object
  OFHelperType of_helper;

  // Set the scaling factors for multi-resolution
  of_helper.SetDefaultPyramidFactors(param.iter_per_level.size());

  // Add random sampling jitter for affine stability at voxel edges
  of_helper.SetJitterSigma(param.affine_jitter);

  // Read the image pairs to register - this will also build the composite pyramids
  ReadImages(param, of_helper);

  // Matrix describing current transform in physical space
  vnl_matrix<double> Q_physical;

  // The number of resolution levels
  unsigned nlevels = param.iter_per_level.size();

  // Clear the metric log
  m_MetricLog.clear();

  // Iterate over the resolution levels
  for(unsigned int level = 0; level < nlevels; ++level)
    {
    // Add stage to metric log
    m_MetricLog.push_back(std::vector<double>());

    // Define the affine cost function
    AbstractAffineCostFunction *pure_acf, *acf;
    if(param.affine_dof == GreedyParameters::DOF_RIGID)
      {
      RigidCostFunction *rigid_acf = new RigidCostFunction(&param, this, level, &of_helper);
      acf = new ScalingCostFunction(
              rigid_acf,
              rigid_acf->GetOptimalParameterScaling(
                of_helper.GetReferenceSpace(level)->GetBufferedRegion().GetSize()));
      pure_acf = rigid_acf;
      }
    else
      {
      //  PureAffineCostFunction *affine_acf = new PureAffineCostFunction(&param, level, &of_helper);
      PhysicalSpaceAffineCostFunction *affine_acf = new PhysicalSpaceAffineCostFunction(&param, this, level, &of_helper);
      acf = new ScalingCostFunction(
              affine_acf,
              affine_acf->GetOptimalParameterScaling(
                of_helper.GetReferenceSpace(level)->GetBufferedRegion().GetSize()));
      pure_acf = affine_acf;
      }

    // Current transform
    typename LinearTransformType::Pointer tLevel = LinearTransformType::New();

    // Set up the initial transform
    if(level == 0)
      {
      // Get the coefficients corresponding to the identity transform in voxel space
      tLevel->SetIdentity();
      vnl_vector<double> xIdent = acf->GetCoefficients(tLevel);

      // Use the provided initial affine as the starting point
      if(param.affine_init_mode == RAS_FILENAME)
        {
        // Read the initial affine transform from a file
        vnl_matrix<double> Qp = this->ReadAffineMatrixViaCache(param.affine_init_transform);

        // Map this to voxel space
        MapPhysicalRASSpaceToAffine(of_helper, level, Qp, tLevel);
        }
      else if(param.affine_init_mode == RAS_IDENTITY)
        {
        // Physical space transform
        vnl_matrix<double> Qp(VDim+1, VDim+1); Qp.set_identity();

        // Map this to voxel space
        MapPhysicalRASSpaceToAffine(of_helper, level, Qp, tLevel);
        }
      else if(param.affine_init_mode == IMG_CENTERS)
        {
        // Find a translation that maps center voxel of fixed image to the center 
        // voxel of the moving image
        vnl_matrix<double> Qp(VDim+1, VDim+1); Qp.set_identity();
        vnl_vector<double> cfix = GetImageCenterinNiftiSpace(of_helper.GetReferenceSpace(level));
        vnl_vector<double> cmov = GetImageCenterinNiftiSpace(of_helper.GetMovingReferenceSpace(level));

        // TODO: I think that setting the matrix above to affine will break the registration
        // if fixed and moving are in different orientations? Or am I crazy?

        // Compute the transform that takes fixed into moving
        for(uint d = 0; d < VDim; d++)
          Qp(d, VDim) = cmov[d] - cfix[d];

        // Map this to voxel space
        MapPhysicalRASSpaceToAffine(of_helper, level, Qp, tLevel);
        }

      // Get the new coefficients
      vnl_vector<double> xInit = acf->GetCoefficients(tLevel);

      // If the voxel-space transform is identity, apply a little bit of jitter
      if((xIdent - xInit).inf_norm() < 1e-4)
        {
        // Apply jitter
        vnl_random rndy(12345);
        for (unsigned i = 0; i < xInit.size(); i++)
          xInit[i] += rndy.drand32(-0.4, 0.4);

        // Map back into transform format
        acf->GetTransform(xInit, tLevel);
        }

      // If the uses asks for rigid search, do it!
      if(param.rigid_search.iterations > 0)
        {
        // Create a pure rigid acf
        RigidCostFunction search_fun(&param, this, level, &of_helper);

        // Get the parameters corresponding to the current transform
        vnl_vector<double> xRigidInit = search_fun.GetCoefficients(tLevel);

        // Get center of fixed and moving images in physical space
        vnl_vector<double> cfix = GetImageCenterinNiftiSpace(of_helper.GetReferenceSpace(level));
        vnl_vector<double> cmov = GetImageCenterinNiftiSpace(of_helper.GetMovingReferenceSpace(level));

        // At random, try a whole bunch of transforms, around 5 degrees
        vnl_random randy(12345);

        // TODO: make a heap of k best tries
        double fBest;
        vnl_vector<double> xBest = xRigidInit;
        search_fun.compute(xBest, &fBest, NULL);

        // Report the initial best
        std::cout << "Rigid search -> Initial best: " << fBest << " " << xBest << std::endl;

        for(int i = 0; i < param.rigid_search.iterations; i++)
          {
          // Get random coefficient
          // Compute a random rotation
          vnl_vector<double> xTry = search_fun.GetRandomCoeff(xRigidInit, randy,
                                                              param.rigid_search.sigma_angle,
                                                              param.rigid_search.sigma_xyz,
                                                              cfix, cmov);

          // Evaluate this transform
          double f;
          search_fun.compute(xTry, &f, NULL);

          if(f < fBest)
            {
            fBest = f;
            xBest = xTry;
            std::cout << "New best: " << fBest << " " << xBest << std::endl;
            }
          }

        xInit = xBest;
        search_fun.GetTransform(xInit, tLevel);
        }
      }
    else
      {
      // Update the transform from the last level
      MapPhysicalRASSpaceToAffine(of_helper, level, Q_physical, tLevel);
      }

    // Test derivatives
    // Convert to a parameter vector
    vnl_vector<double> xLevel = acf->GetCoefficients(tLevel.GetPointer());

    if(param.flag_debug_deriv)
      {
      // Test the gradient computation
      vnl_vector<double> xGrad(acf->get_number_of_unknowns(), 0.0);
      double f0;
      acf->compute(xLevel, &f0, &xGrad);

      // Propagate the jitter to the transform
      Q_physical = MapAffineToPhysicalRASSpace(of_helper, level, tLevel);
      std::cout << "Initial RAS Transform: " << std::endl << Q_physical  << std::endl;

      printf("ANL gradient: ");
      for (unsigned i = 0; i < xGrad.size(); i++)
        printf("%11.4f ", xGrad[i]);
      printf("\n");

      vnl_vector<double> xGradN(acf->get_number_of_unknowns(), 0.0);
      for(int i = 0; i < acf->get_number_of_unknowns(); i++)
        {
        // double eps = (i % VDim == 0) ? 1.0e-2 : 1.0e-5;
        double eps = param.deriv_epsilon;
        double f1, f2, f3, f4;
        vnl_vector<double> x1 = xLevel, x2 = xLevel, x3 = xLevel, x4 = xLevel;
        x1[i] -= 2 * eps; x2[i] -= eps; x3[i] += eps; x4[i] += 2 * eps;

        // Four-point derivative computation
        acf->compute(x1, &f1, NULL);
        acf->compute(x2, &f2, NULL);
        acf->compute(x3, &f3, NULL);
        acf->compute(x4, &f4, NULL);

        xGradN[i] = (f1 - 8 * f2 + 8 * f3 - f4) / (12 * eps);
        }

      printf("NUM gradient: ");
      for (unsigned i = 0; i < xGradN.size(); i++)
        printf("%11.4f ", xGradN[i]);
      printf("\n");

      std::cout << "f = " << f0 << std::endl;

      acf->GetTransform(xGrad, tLevel.GetPointer());
      std::cout << "A: " << std::endl
                << tLevel->GetMatrix() << std::endl
                << tLevel->GetOffset() << std::endl;

      acf->GetTransform(xGradN, tLevel.GetPointer());
      std::cout << "N: " << std::endl
                << tLevel->GetMatrix() << std::endl
                << tLevel->GetOffset() << std::endl;
      }

    if(param.flag_debug_aff_obj)
      {
      for(int k = -50; k < 50; k++)
        {
        printf("Obj\t%d\t", k);
        for(int i = 0; i < acf->get_number_of_unknowns(); i++)
          {
          vnl_vector<double> xTest = xLevel;
          xTest[i] = xLevel[i] + k * param.deriv_epsilon;
          double f; acf->compute(xTest, &f, NULL);
          printf("%12.8f\t", f);
          }
        printf("\n");
        }
        {
        vnl_vector<double> xTest = xLevel;
          {
          }
        printf("\n");
        }
      }

    // Run the minimization
    if(param.iter_per_level[level] > 0)
      {
      if(param.flag_powell)
        {
        // Set up the optimizer
        vnl_powell *optimizer = new vnl_powell(acf);
        optimizer->set_f_tolerance(1e-9);
        optimizer->set_x_tolerance(1e-4);
        optimizer->set_g_tolerance(1e-6);
        optimizer->set_trace(true);
        optimizer->set_verbose(true);
        optimizer->set_max_function_evals(param.iter_per_level[level]);

        optimizer->minimize(xLevel);
        delete optimizer;

        }
      else
        {
        // Set up the optimizer
        vnl_lbfgs *optimizer = new vnl_lbfgs(*acf);
        optimizer->set_f_tolerance(1e-9);
        optimizer->set_x_tolerance(1e-4);
        optimizer->set_g_tolerance(1e-6);
        optimizer->set_trace(true);
        optimizer->set_max_function_evals(param.iter_per_level[level]);

        optimizer->minimize(xLevel);
        delete optimizer;
        }

      // Did the registration succeed?
      if(xLevel.size() > 0)
        {
        // Get the final transform
        typename LinearTransformType::Pointer tFinal = LinearTransformType::New();
        acf->GetTransform(xLevel, tFinal.GetPointer());
        Q_physical = MapAffineToPhysicalRASSpace(of_helper, level, tFinal);
        }
      else
        {
        // Use the pre-initialization transform parameters
        Q_physical = MapAffineToPhysicalRASSpace(of_helper, level, tLevel);
        }
      }

    std::cout << "Final RAS Transform: " << std::endl << Q_physical << std::endl;

    delete acf;
    delete pure_acf;
    }

  // Write the final affine transform
  this->WriteAffineMatrixViaCache(param.output, Q_physical);
  return 0;
}

#include "itkStatisticsImageFilter.h"


/**
 * This is the main function of the GreedyApproach algorithm
 */
template <unsigned int VDim, typename TReal>
int GreedyApproach<VDim, TReal>
::RunDeformable(GreedyParameters &param)
{
  // Create an optical flow helper object
  OFHelperType of_helper;

  // Set the scaling factors for multi-resolution
  of_helper.SetDefaultPyramidFactors(param.iter_per_level.size());

  // Set the scaling mode depending on the metric
  if(param.metric == GreedyParameters::MAHALANOBIS)
    of_helper.SetScaleFixedImageWithVoxelSize(true);

  // Read the image pairs to register
  ReadImages(param, of_helper);

  // An image pointer desribing the current estimate of the deformation
  VectorImagePointer uLevel = NULL;

  // The number of resolution levels
  unsigned nlevels = param.iter_per_level.size();

  // Iterate over the resolution levels
  for(unsigned int level = 0; level < nlevels; ++level)
    {
    // Reference space
    ImageBaseType *refspace = of_helper.GetReferenceSpace(level);

    // Smoothing factors for this level, in physical units
    typename LDDMMType::Vec sigma_pre_phys =
        of_helper.GetSmoothingSigmasInPhysicalUnits(level, param.sigma_pre.sigma,
                                                    param.sigma_pre.physical_units);

    typename LDDMMType::Vec sigma_post_phys =
        of_helper.GetSmoothingSigmasInPhysicalUnits(level, param.sigma_post.sigma,
                                                    param.sigma_post.physical_units);

    // Report the smoothing factors used
    std::cout << "LEVEL " << level+1 << " of " << nlevels << std::endl;
    std::cout << "  Smoothing sigmas: " << sigma_pre_phys << ", " << sigma_post_phys << std::endl;

    // Set up timers for different critical components of the optimization
    itk::TimeProbe tm_Gradient, tm_Gaussian1, tm_Gaussian2, tm_Iteration, 
      tm_Integration, tm_Update;

    // Intermediate images
    ImagePointer iTemp = ImageType::New();
    VectorImagePointer viTemp = VectorImageType::New();
    VectorImagePointer uk = VectorImageType::New();
    VectorImagePointer uk1 = VectorImageType::New();

    // This is the exponentiated uk, in stationary velocity mode it is uk^(2^N)
    VectorImagePointer uk_exp = VectorImageType::New();

    // A pointer to the full warp image - either uk in greedy mode, or uk_exp in diff demons mdoe
    VectorImageType *uFull;

    // Matrix work image (for Lie Bracket) 
    typedef typename LDDMMType::MatrixImageType MatrixImageType;
    typename MatrixImageType::Pointer work_mat = MatrixImageType::New();

    // Allocate the intermediate data
    LDDMMType::alloc_vimg(uk, refspace);
    if(param.iter_per_level[level] > 0)
      {
      LDDMMType::alloc_img(iTemp, refspace);
      LDDMMType::alloc_vimg(viTemp, refspace);
      LDDMMType::alloc_vimg(uk1, refspace);

      // These are only allocated in diffeomorphic demons mode
      if(param.flag_stationary_velocity_mode)
        {
        LDDMMType::alloc_vimg(uk_exp, refspace);
        LDDMMType::alloc_mimg(work_mat, refspace);
        }
      }

    // Initialize the deformation field from last iteration
    if(uLevel.IsNotNull())
      {
      LDDMMType::vimg_resample_identity(uLevel, refspace, uk);
      LDDMMType::vimg_scale_in_place(uk, 2.0);
      uLevel = uk;
      }
    else if(param.initial_warp.size())
      {
      // The user supplied an initial warp or initial root warp. In this case, we
      // do not start iteration from zero, but use the initial warp to start from
      VectorImagePointer uInit = VectorImageType::New();

      // Read the warp file
      LDDMMType::vimg_read(param.initial_warp.c_str(), uInit );

      // Convert the warp file into voxel units from physical units
      OFHelperType::PhysicalWarpToVoxelWarp(uInit, uInit, uInit);

      // Scale the initial warp by the pyramid level
      LDDMMType::vimg_resample_identity(uInit, refspace, uk);
      LDDMMType::vimg_scale_in_place(uk, 1.0 / (1 << level));
      uLevel = uk;
      itk::Index<VDim> test; test.Fill(24);
      std::cout << "Index 24x24x24 maps to " << uInit->GetPixel(test) << std::endl;
      std::cout << "Index 24x24x24 maps to " << uk->GetPixel(test) << std::endl;
      }
    else if(param.affine_init_mode != VOX_IDENTITY)
      {
      typename LinearTransformType::Pointer tran = LinearTransformType::New();

      if(param.affine_init_mode == RAS_FILENAME)
        {
        // Read the initial affine transform from a file
        vnl_matrix<double> Qp = ReadAffineMatrixViaCache(param.affine_init_transform);

        // Map this to voxel space
        MapPhysicalRASSpaceToAffine(of_helper, level, Qp, tran);
        }
      else if(param.affine_init_mode == RAS_IDENTITY)
        {
        // Physical space transform
        vnl_matrix<double> Qp(VDim+1, VDim+1); Qp.set_identity();

        // Map this to voxel space
        MapPhysicalRASSpaceToAffine(of_helper, level, Qp, tran);
        }

      // Create an initial warp
      OFHelperType::AffineToField(tran, uk);
      uLevel = uk;

      itk::Index<VDim> test; test.Fill(24);
      std::cout << "Index 24x24x24 maps to " << uk->GetPixel(test) << std::endl;
      }

    // Iterate for this level
    for(int iter = 0; iter < param.iter_per_level[level]; iter++)
      {
      // Start the iteration timer
      tm_Iteration.Start();

      // The epsilon for this level
      double eps= param.epsilon_per_level[level];

      // Compute the gradient of objective
      double total_energy;

      // Integrate the total deformation field for this iteration
      if(param.flag_stationary_velocity_mode)
        {
        tm_Integration.Start();

        // This is the exponentiation of the stationary velocity field
        // Take current warp to 'exponent' power - this is the actual warp
        LDDMMType::vimg_exp(uk, uk_exp, viTemp, param.warp_exponent, 1.0);
        uFull = uk_exp;

        tm_Integration.Stop();
        }
      else
        {
        uFull = uk;
        }

      if(param.metric == GreedyParameters::SSD)
        {
        // Begin gradient computation
        tm_Gradient.Start();

        vnl_vector<double> all_metrics =
            of_helper.ComputeOpticalFlowField(level, uFull, iTemp, uk1, eps)  / eps;

        // If there is a mask, multiply the gradient by the mask
        if(param.gradient_mask.size())
          LDDMMType::vimg_multiply_in_place(uk1, of_helper.GetGradientMask(level));

        // End gradient computation
        tm_Gradient.Stop();

        printf("Lev:%2d  Itr:%5d  Met:[", level, iter);
        total_energy = 0.0;
        for (unsigned i = 0; i < all_metrics.size(); i++)
          {
          printf("  %8.6f", all_metrics[i]);
          total_energy += all_metrics[i];
          }
        printf("]  Tot: %8.6f\n", total_energy);
        }

      else if(param.metric == GreedyParameters::MI || param.metric == GreedyParameters::NMI)
        {
        // Begin gradient computation
        tm_Gradient.Start();

        vnl_vector<double> all_metrics =
            of_helper.ComputeMIFlowField(level, param.metric == GreedyParameters::NMI, uFull, iTemp, uk1, eps);

        // If there is a mask, multiply the gradient by the mask
        if(param.gradient_mask.size())
          LDDMMType::vimg_multiply_in_place(uk1, of_helper.GetGradientMask(level));

        // End gradient computation
        tm_Gradient.Stop();

        printf("Lev:%2d  Itr:%5d  Met:[", level, iter);
        total_energy = 0.0;
        for (unsigned i = 0; i < all_metrics.size(); i++)
          {
          printf("  %8.6f", all_metrics[i]);
          total_energy += all_metrics[i];
          }
        printf("]  Tot: %8.6f\n", total_energy);
        }

      else if(param.metric == GreedyParameters::NCC)
        {
        itk::Size<VDim> radius = array_caster<VDim>::to_itkSize(param.metric_radius);

        // Test derivative
        // total_energy = of_helper.ComputeNCCMetricAndGradient(level, uk, uk1, radius, param.epsilon);

        /*
        if(iter == 0)
          {

          // Perform a derivative check!

          itk::Index<VDim> test; test.Fill(24);
          typename VectorImageType::PixelType vtest = uk->GetPixel(test), vv;

          itk::ImageRegion<VDim> region = uk1->GetBufferedRegion();
          // region.ShrinkByRadius(1);

          double eps = param.epsilon;
          for(uint d = 0; d < VDim; d++)
            {
            vv.Fill(0.5); vv[d] -= eps; uk->FillBuffer(vv);
            of_helper.ComputeNCCMetricImage(level, uk, radius, iTemp, uk1, 1.0);

            double a1 = 0.0;
            typedef itk::ImageRegionConstIterator<ImageType> Iter;
            for(Iter it(iTemp, region); !it.IsAtEnd(); ++it)
              {
              a1 += it.Get();
              }


            vv.Fill(0.5); vv[d] += eps; uk->FillBuffer(vv);
            of_helper.ComputeNCCMetricImage(level, uk, radius, iTemp, uk1, 1.0);

            double a2 = 0.0;
            typedef itk::ImageRegionConstIterator<ImageType> Iter;
            for(Iter it(iTemp, region); !it.IsAtEnd(); ++it)
              {
              a2 += it.Get();
              }

            std::cout << "NUM:" << (a2 - a1) / (2*eps) << std::endl;

            }

          vv.Fill(0.5); uk->FillBuffer(vv);
          total_energy = of_helper.ComputeNCCMetricImage(level, uk, radius, iTemp, uk1, 1.0);
          for(uint d = 0; d < VDim; d++)
            {

            double ader = 0.0;
            typedef itk::ImageRegionConstIterator<VectorImageType> Iter;
            for(Iter it(uk1, region); !it.IsAtEnd(); ++it)
              {
              ader += it.Get()[d];
              }

            // itk::Index<VDim> test; test.Fill(24);
            // std::cout << "ANA:" << uk1->GetPixel(test) << std::endl;

            std::cout << "ANA:" << ader << std::endl;
            }
          }
          */

        // Begin gradient computation
        tm_Gradient.Start();

        // Compute the metric - no need to multiply by the mask, this happens already in the NCC metric code
        total_energy = of_helper.ComputeNCCMetricImage(level, uFull, radius, iTemp, uk1, eps) / eps;

        // End gradient computation
        tm_Gradient.Stop();

        printf("Level %5d,  Iter %5d:    Energy = %8.4f\n", level, iter, total_energy);
        fflush(stdout);
        }
      else if(param.metric == GreedyParameters::MAHALANOBIS)
        {
        tm_Gradient.Start();
        double total_energy = of_helper.ComputeMahalanobisMetricImage(level, uFull, iTemp, uk1);
        tm_Gradient.Stop();
        printf("Level %5d,  Iter %5d:    Energy = %8.4f\n", level, iter, total_energy);
        fflush(stdout);
        }

      // Dump the gradient image if requested
      if(param.flag_dump_moving && 0 == iter % param.dump_frequency)
        {
        char fname[256];
        //sprintf(fname, "dump_gradient_lev%02d_iter%04d.nii.gz", level, iter);
        LDDMMType::vimg_write(uk1, fname);
        }

      // We have now computed the gradient vector field. Next, we smooth it
      tm_Gaussian1.Start();
      LDDMMType::vimg_smooth_withborder(uk1, viTemp, sigma_pre_phys, 1);
      tm_Gaussian1.Stop();

      // After smoothing, compute the maximum vector norm and use it as a normalizing
      // factor for the displacement field
      if(param.time_step_mode == GreedyParameters::SCALE)
        LDDMMType::vimg_normalize_to_fixed_max_length(viTemp, iTemp, eps, false);
      else if (param.time_step_mode == GreedyParameters::SCALEDOWN)
        LDDMMType::vimg_normalize_to_fixed_max_length(viTemp, iTemp, eps, true);

      // Dump the smoothed gradient image if requested
      if(param.flag_dump_moving && 0 == iter % param.dump_frequency)
        {
        char fname[256];
        //sprintf(fname, "dump_optflow_lev%02d_iter%04d.nii.gz", level, iter);
        LDDMMType::vimg_write(viTemp, fname);
        }

      // Compute the updated deformation field - in uk1
      tm_Update.Start();
      if(param.flag_stationary_velocity_mode)
        {
        // this is diffeomorphic demons - Vercauteren 2008
        // We now hold the update field in viTemp. This update u should be integrated
        // with the current stationary velocity field such that exp[v'] = exp[v] o exp[u]
        // Vercauteren (2008) suggests using the following expressions
        // v' = v + u (so-so)
        // v' = v + u + [v, u]/2 (this is the Lie bracket)
        
        // Scale the update by 1 / 2^exponent (tiny update, first order approximation)
        LDDMMType::vimg_scale_in_place(viTemp, 1.0 / (2 << param.warp_exponent));

        // Use appropriate update
        if(param.flag_stationary_velocity_mode_use_lie_bracket)
          {
          // Use the Lie Bracket approximation (v + u + [v,u])
          LDDMMType::lie_bracket(uk, viTemp, work_mat, uk1);
          LDDMMType::vimg_scale_in_place(uk1, 0.5); 
          LDDMMType::vimg_add_in_place(uk1, uk);
          LDDMMType::vimg_add_in_place(uk1, viTemp);
          }
        else
          {
          LDDMMType::vimg_copy(uk, uk1);
          LDDMMType::vimg_add_in_place(uk1, viTemp);
          }
        }
      else
        {
        // This is compositive (uk1 = viTemp + uk o viTemp), which is what is done with
        // compositive demons and ANTS
        LDDMMType::interp_vimg(uk, viTemp, 1.0, uk1);
        LDDMMType::vimg_add_in_place(uk1, viTemp);
        }
      tm_Update.Stop();

      // Dump if requested
      if(param.flag_dump_moving && 0 == iter % param.dump_frequency)
        {
        char fname[256];
        //sprintf(fname, "dump_uk1_lev%02d_iter%04d.nii.gz", level, iter);
        LDDMMType::vimg_write(uk1, fname);
        }

      // Another layer of smoothing
      tm_Gaussian2.Start();
      LDDMMType::vimg_smooth_withborder(uk1, uk, sigma_post_phys, 1);
      tm_Gaussian2.Stop();

      tm_Iteration.Stop();
      }

    // Store the end result
    uLevel = uk;

    // Compute the jacobian of the deformation field - but only if we iterated at this level
    if(param.iter_per_level[level] > 0)
      {
      LDDMMType::field_jacobian_det(uk, iTemp);
      TReal jac_min, jac_max;
      LDDMMType::img_min_max(iTemp, jac_min, jac_max);
      printf("END OF LEVEL %5d    DetJac Range: %8.4f  to %8.4f \n", level, jac_min, jac_max);

      // Print timing information
      printf("  Avg. Gradient Time  : %6.4fs  %5.2f%% \n", tm_Gradient.GetMean(), 
             tm_Gradient.GetMean() * 100.0 / tm_Iteration.GetMean());
      printf("  Avg. Gaussian Time  : %6.4fs  %5.2f%% \n", tm_Gaussian1.GetMean() + tm_Gaussian2.GetMean(),
             (tm_Gaussian1.GetMean() + tm_Gaussian2.GetMean()) * 100.0 / tm_Iteration.GetMean());
      printf("  Avg. Integration Time  : %6.4fs  %5.2f%% \n", tm_Integration.GetMean() + tm_Update.GetMean(),
             (tm_Integration.GetMean() + tm_Update.GetMean()) * 100.0 / tm_Iteration.GetMean());
      printf("  Avg. Total Iteration Time : %6.4fs \n", tm_Iteration.GetMean());
      }
    }

  // The transformation field is in voxel units. To work with ANTS, it must be mapped
  // into physical offset units - just scaled by the spacing?
  
  if(param.flag_stationary_velocity_mode)
    {
    // Take current warp to 'exponent' power - this is the actual warp
    VectorImagePointer uLevelExp = LDDMMType::alloc_vimg(uLevel);
    VectorImagePointer uLevelWork = LDDMMType::alloc_vimg(uLevel);
    LDDMMType::vimg_exp(uLevel, uLevelExp, uLevelWork, param.warp_exponent, 1.0);

    // Write the resulting transformation field
    of_helper.WriteCompressedWarpInPhysicalSpace(nlevels - 1, uLevelExp, param.output.c_str(), param.warp_precision);

    if(param.root_warp.size())
      {
      // If asked to write root warp, do so
      of_helper.WriteCompressedWarpInPhysicalSpace(nlevels - 1, uLevel, param.root_warp.c_str(), 0);
      }
    if(param.inverse_warp.size())
      {
      // Compute the inverse (this is probably unnecessary for small warps)
      of_helper.ComputeDeformationFieldInverse(uLevel, uLevelWork, 0);
      of_helper.WriteCompressedWarpInPhysicalSpace(nlevels - 1, uLevelWork, param.inverse_warp.c_str(), param.warp_precision);
      }
    }
  else
    {
    // Write the resulting transformation field
    of_helper.WriteCompressedWarpInPhysicalSpace(nlevels - 1, uLevel, param.output.c_str(), param.warp_precision);

    // If an inverse is requested, compute the inverse using the Chen 2008 fixed method.
    // A modification of this method is that if convergence is slow, we take the square
    // root of the forward transform.
    //
    // TODO: it would be more efficient to check the Lipschitz condition rather than
    // the brute force approach below
    //
    // TODO: the maximum checks should only be done over the region where the warp is
    // not going outside of the image. Right now, they are meaningless and we are doing
    // extra work when computing the inverse.
    if(param.inverse_warp.size())
      {
      // Compute the inverse
      VectorImagePointer uInverse = VectorImageType::New();
      LDDMMType::alloc_vimg(uInverse, uLevel);
      of_helper.ComputeDeformationFieldInverse(uLevel, uInverse, param.warp_exponent);

      // Write the warp using compressed format
      of_helper.WriteCompressedWarpInPhysicalSpace(nlevels - 1, uInverse, param.inverse_warp.c_str(), param.warp_precision);
      }
    }
  return 0;
}

/**
 * This function performs brute force search for similar patches. It generates a discrete displacement
 * field where every pixel in the fixed image is matched to the most similar pixel in the moving image
 * within a certain radius
 */
template <unsigned int VDim, typename TReal>
int GreedyApproach<VDim, TReal>
::RunBrute(GreedyParameters &param)
{
  // Check for valid parameters
  if(param.metric != GreedyParameters::NCC)
    {
    std::cerr << "Brute force search requires NCC metric only" << std::endl;
    return -1;
    }

  if(param.brute_search_radius.size() != VDim)
    {
    std::cerr << "Brute force search radius must be same dimension as the images" << std::endl;
    return -1;
    }

  // Create an optical flow helper object
  OFHelperType of_helper;

  // No multi-resolution
  of_helper.SetDefaultPyramidFactors(1);

  // Read the image pairs to register
  ReadImages(param, of_helper);

  // Reference space
  ImageBaseType *refspace = of_helper.GetReferenceSpace(0);

  // Intermediate images
  VectorImagePointer u_best = VectorImageType::New();
  VectorImagePointer u_curr = VectorImageType::New();
  ImagePointer m_curr = ImageType::New();
  ImagePointer m_best = ImageType::New();

  // Allocate the intermediate data
  LDDMMType::alloc_vimg(u_best, refspace);
  LDDMMType::alloc_vimg(u_curr, refspace);
  LDDMMType::alloc_img(m_best, refspace);
  LDDMMType::alloc_img(m_curr, refspace);

  // Allocate m_best to a negative value
  m_best->FillBuffer(-100.0);

  // Create a neighborhood for computing offsets
  itk::Neighborhood<float, VDim> dummy_nbr;
  itk::Size<VDim> search_rad = array_caster<VDim>::to_itkSize(param.brute_search_radius);
  itk::Size<VDim> metric_rad = array_caster<VDim>::to_itkSize(param.metric_radius);
  dummy_nbr.SetRadius(search_rad);

  // Iterate over all offsets
  for(uint k = 0; k < dummy_nbr.Size(); k++)
    {
    // Get the offset corresponding to this iteration
    itk::Offset<VDim> offset = dummy_nbr.GetOffset(k);

    // Fill the deformation field with this offset
    typename LDDMMType::Vec vec_offset;
    for(uint i = 0; i < VDim; i++)
      vec_offset[i] = offset[i];
    u_curr->FillBuffer(vec_offset);

    // Perform interpolation and metric computation
    of_helper.ComputeNCCMetricImage(0, u_curr, metric_rad, m_curr);

    // Temp: keep track of number of updates
    unsigned long n_updates = 0;

    // Out of laziness, just take a quick pass over the images
    typename VectorImageType::RegionType rgn = refspace->GetBufferedRegion();
    itk::ImageRegionIterator<VectorImageType> it_u(u_best, rgn);
    itk::ImageRegionConstIterator<ImageType> it_m_curr(m_curr, rgn);
    itk::ImageRegionIterator<ImageType> it_m_best(m_best, rgn);
    for(; !it_m_best.IsAtEnd(); ++it_m_best, ++it_m_curr, ++it_u)
      {
      float v_curr = it_m_curr.Value();
      if(v_curr > it_m_best.Value())
        {
        it_m_best.Set(v_curr);
        it_u.Set(vec_offset);
        ++n_updates;
        }
      }

    std::cout << "offset: " << offset << "     updates: " << n_updates << std::endl;
    }

  LDDMMType::vimg_write(u_best, param.output.c_str());
  LDDMMType::img_write(m_best, "mbest.nii.gz");

  return 0;
}


#include "itkWarpVectorImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"


template <unsigned int VDim, typename TReal>
void GreedyApproach<VDim, TReal>
::ReadTransformChain(const std::vector<TransformSpec> &tran_chain,
                     ImageBaseType *ref_space,
                     VectorImagePointer &out_warp)
{
  // Create the initial transform and set it to zero
  out_warp = VectorImageType::New();
  LDDMMType::alloc_vimg(out_warp, ref_space);

  // Read the sequence of transforms
  for(uint i = 0; i < tran_chain.size(); i++)
    {
    // Read the next parameter
    std::string tran = tran_chain[i].filename;

    // Determine if it's an affine transform
    if(itk::ImageIOFactory::CreateImageIO(tran.c_str(), itk::ImageIOFactory::ReadMode))
      {
      // Create a temporary warp
      VectorImagePointer warp_tmp = VectorImageType::New();
      LDDMMType::alloc_vimg(warp_tmp, ref_space);

      // Read the next warp
      VectorImagePointer warp_i = VectorImageType::New();
      LDDMMType::vimg_read(tran.c_str(), warp_i);

      // If there is an exponent on the transform spec, handle it
      if(tran_chain[i].exponent != 1)
        {
        // The exponent may be specified as a negative number, in which case we take the negative
        // input and exponentiate it
        double absexp = fabs(tran_chain[i].exponent);
        double n_real = log(absexp) / log(2.0);
        int n = (int) (n_real + 0.5);
        if(fabs(n - n_real) > 1.0e-4) 
          throw GreedyException("Currently only power of two exponents are supported for warps");

        // Bring the transform into voxel space
        VectorImagePointer warp_exp = LDDMMType::alloc_vimg(warp_i);
        OFHelperType::PhysicalWarpToVoxelWarp(warp_i, warp_i, warp_i);

        // Square the transform N times (in its own space)
        LDDMMType::vimg_exp(warp_i, warp_exp, warp_tmp, n, tran_chain[i].exponent / absexp);

        // Bring the transform back into physical space
        OFHelperType::VoxelWarpToPhysicalWarp(warp_exp, warp_i, warp_i);
        }

      // Now we need to compose the current transform and the overall warp.
      LDDMMType::interp_vimg(warp_i, out_warp, 1.0, warp_tmp, false, true);
      LDDMMType::vimg_add_in_place(out_warp, warp_tmp);
      }
    else
      {
      // Read the transform as a matrix
      vnl_matrix<double> mat = ReadAffineMatrixViaCache(tran_chain[i]);
      vnl_matrix<double>  A = mat.extract(VDim, VDim);
      vnl_vector<double> b = mat.get_column(VDim).extract(VDim), q;

      // TODO: stick this in a filter to take advantage of threading!
      typedef itk::ImageRegionIteratorWithIndex<VectorImageType> IterType;
      for(IterType it(out_warp, out_warp->GetBufferedRegion()); !it.IsAtEnd(); ++it)
        {
        itk::Point<double, VDim> pt, pt2;
        typename VectorImageType::IndexType idx = it.GetIndex();

        // Get the physical position
        // TODO: this calls IsInside() internally, which limits efficiency
        out_warp->TransformIndexToPhysicalPoint(idx, pt);

        // Add the displacement (in DICOM coordinates) and
        for(uint i = 0; i < VDim; i++)
          pt2[i] = pt[i] + it.Value()[i];

        // Switch to NIFTI coordinates
        pt2[0] = -pt2[0]; pt2[1] = -pt2[1];

        // Apply the matrix - get the transformed coordinate in DICOM space
        q = A * pt2.GetVnlVector() + b;
        q[0] = -q[0]; q[1] = -q[1];

        // Compute the difference in DICOM space
        for(uint i = 0; i < VDim; i++)
          it.Value()[i] = q[i] - pt[i];
        }
      }
    }
}

#include "itkBinaryThresholdImageFilter.h"
//#include "itkRecursiveGaussianImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkNaryFunctorImageFilter.h"

template <class TInputImage, class TOutputImage>
class NaryLabelVotingFunctor
{
public:
  typedef NaryLabelVotingFunctor<TInputImage,TOutputImage> Self;
  typedef typename TInputImage::PixelType InputPixelType;
  typedef typename TOutputImage::PixelType OutputPixelType;
  typedef std::vector<OutputPixelType> LabelArray;

  NaryLabelVotingFunctor(const LabelArray &labels)
    : m_LabelArray(labels), m_Size(labels.size()) {}

  NaryLabelVotingFunctor() : m_Size(0) {}


  OutputPixelType operator() (const std::vector<InputPixelType> &pix)
  {
    InputPixelType best_val = pix[0];
    int best_index = 0;
    for(int i = 1; i < m_Size; i++)
      if(pix[i] > best_val)
        {
        best_val = pix[i];
        best_index = i;
        }

    return m_LabelArray[best_index];
  }

  bool operator != (const Self &other)
    { return other.m_LabelArray != m_LabelArray; }

protected:
  LabelArray m_LabelArray;
  int m_Size;
};

#include "itkMeshFileReader.h"
#include "itkMeshFileWriter.h"
#include "itkMesh.h"
#include "itkTransformMeshFilter.h"

template <unsigned int VDim, typename TArray>
class PhysicalCoordinateTransform
{
  static void ras_to_lps(const TArray &src, TArray &trg) {}
  static void lps_to_ras(const TArray &src, TArray &trg) {}
};

template <typename TArray>
class PhysicalCoordinateTransform<2, TArray>
{
public:
  static void ras_to_lps(const TArray &src, TArray &trg)
  {
    trg[0] = -src[0];
    trg[1] = -src[1];
  }

  static void lps_to_ras(const TArray &src, TArray &trg)
  {
    trg[0] = -src[0];
    trg[1] = -src[1];
  }
};

template <typename TArray>
class PhysicalCoordinateTransform<3, TArray>
{
public:
  static void ras_to_lps(const TArray &src, TArray &trg)
  {
    trg[0] = -src[0];
    trg[1] = -src[1];
    trg[2] = src[2];
  }

  static void lps_to_ras(const TArray &src, TArray &trg)
  {
    trg[0] = -src[0];
    trg[1] = -src[1];
    trg[2] = src[2];
  }
};

template <typename TArray>
class PhysicalCoordinateTransform<4, TArray>
{
public:
  static void ras_to_lps(const TArray &src, TArray &trg)
  {
    trg[0] = -src[0];
    trg[1] = -src[1];
    trg[2] = src[2];
    trg[3] = src[3];
  }

  static void lps_to_ras(const TArray &src, TArray &trg)
  {
    trg[0] = -src[0];
    trg[1] = -src[1];
    trg[2] = src[2];
    trg[3] = src[3];
  }
};



template <unsigned int VDim, typename TReal>
class WarpMeshTransformFunctor : public itk::DataObject
{
public:
  typedef WarpMeshTransformFunctor<VDim, TReal>       Self;
  typedef itk::DataObject                             Superclass;
  typedef itk::SmartPointer<Self>                     Pointer;
  typedef itk::SmartPointer<const Self>               ConstPointer;

  itkTypeMacro(WarpMeshTransformFunctor, itk::DataObject)
  itkNewMacro(Self)

  typedef GreedyApproach<VDim, TReal> GreedyAPI;
  typedef typename GreedyAPI::VectorImageType VectorImageType;
  typedef typename GreedyAPI::ImageBaseType ImageBaseType;
  typedef FastLinearInterpolator<VectorImageType, TReal, VDim> FastInterpolator;
  typedef itk::ContinuousIndex<TReal, VDim> CIndexType;
  typedef itk::Point<TReal, VDim> PointType;

  void SetWarp(VectorImageType *warp)
  {
    if(m_Interpolator) delete m_Interpolator;
    m_Interpolator = new FastInterpolator(warp);
  }

  void SetReferenceSpace(ImageBaseType *ref)
  {
    m_ReferenceSpace = ref;
  }

  PointType TransformPoint(const PointType &x)
  {
    // Our convention is to use NIFTI/RAS coordinates for meshes, whereas ITK
    // uses the DICOM/LPS convention. We transform point to LPS first
    PointType x_lps, phi_x;

    PhysicalCoordinateTransform<VDim, PointType>::ras_to_lps(x, x_lps);

    CIndexType cix;
    typename VectorImageType::PixelType vec;
    vec.Fill(0.0);
    m_ReferenceSpace->TransformPhysicalPointToContinuousIndex(x_lps, cix);
    m_Interpolator->Interpolate(cix.GetDataPointer(), &vec);

    for(uint d = 0; d < VDim; d++)
      {
      phi_x[d] = vec[d] + x_lps[d];
      }


    PhysicalCoordinateTransform<VDim, PointType>::lps_to_ras(phi_x, phi_x);

    return phi_x;
  }

protected:

  WarpMeshTransformFunctor() { m_Interpolator = NULL; }
  ~WarpMeshTransformFunctor()
  {
    if(m_Interpolator)
      delete m_Interpolator;
  }

private:

  typename ImageBaseType::Pointer m_ReferenceSpace;
  FastInterpolator *m_Interpolator;

};

/**
 * This code computes the jacobian determinant field for a deformation. The
 * recommended mode for this computation is to take the k-th root of the input
 * transformation and then compose the Jacobians
 */
template <unsigned int VDim, typename TReal>
int GreedyApproach<VDim, TReal>
::RunJacobian(GreedyParameters &param)
{
  // Read the warp as a transform chain
  VectorImagePointer warp;

  // Read the warp file
  LDDMMType::vimg_read(param.jacobian_param.in_warp.c_str(), warp);

  // Convert the warp file into voxel units from physical units
  OFHelperType::PhysicalWarpToVoxelWarp(warp, warp, warp);

  // Compute the root of the warp
  VectorImagePointer root_warp = VectorImageType::New();
  LDDMMType::alloc_vimg(root_warp, warp);

  // Allocate a working warp
  VectorImagePointer work_warp = VectorImageType::New();
  LDDMMType::alloc_vimg(work_warp, warp);

  // Compute the root warp, which is not stored in the variable warp
  OFHelperType::ComputeWarpRoot(warp, root_warp, param.warp_exponent);

  // Initialize empty array of Jacobians
  typedef typename LDDMMType::MatrixImageType JacobianImageType;
  typename JacobianImageType::Pointer jac = JacobianImageType::New();
  LDDMMType::alloc_mimg(jac, warp);

  typename JacobianImageType::Pointer jac_work = JacobianImageType::New();
  LDDMMType::alloc_mimg(jac_work, warp);

  // Compute the Jacobian of the root warp
  LDDMMType::field_jacobian(root_warp, jac);

  // Compute the Jacobian matrix of the root warp; jac[a] = D_a (warp)
  for(int k = 0; k < param.warp_exponent; k++)
    {
    // Compute the composition of the Jacobian with itself
    LDDMMType::jacobian_of_composition(jac, jac, root_warp, jac_work);

    // Swap the pointers, so jac points to the actual composed jacobian
    typename JacobianImageType::Pointer temp = jac_work.GetPointer();
    jac_work = jac.GetPointer();
    jac = temp.GetPointer();

    // Compute the composition of the warp with itself, place into root_warp
    LDDMMType::interp_vimg(root_warp, root_warp, 1.0, work_warp);
    LDDMMType::vimg_add_in_place(root_warp, work_warp);
    }

  // At this point, root_warp should hold the original warp, and jac+I will hold
  // the Jacobian of the original warp. We need to compute the determinant
  ImagePointer jac_det = ImageType::New();
  LDDMMType::alloc_img(jac_det, warp);
  LDDMMType::mimg_det(jac, 1.0, jac_det);

  // Write the computed Jacobian
  LDDMMType::img_write(jac_det, param.jacobian_param.out_det_jac.c_str(), itk::ImageIOBase::FLOAT);
  return 0;
}

/**
 * Run the reslice code - simply apply a warp or set of warps to images
 */
template <unsigned int VDim, typename TReal>
int GreedyApproach<VDim, TReal>
::RunReslice(GreedyParameters &param)
{
  GreedyResliceParameters r_param = param.reslice_param;

  // Check the parameters
  if(!r_param.ref_image.size())
    throw GreedyException("A reference image (-rf) option is required for reslice commands");

  if(r_param.images.size() + r_param.meshes.size() == 0
     && !r_param.out_composed_warp.size()
     && !r_param.out_jacobian_image.size())
    throw GreedyException("No operation specified for reslice mode. "
                          "Use one of -rm, -rs or -rc commands.");

  // Read the fixed as a plain image (we don't care if it's composite)
  ImagePointer ref = ImageType::New();
  LDDMMType::img_read(r_param.ref_image.c_str(), ref);
  itk::ImageBase<VDim> *ref_space = ref;

  // Read the transform chain
  VectorImagePointer warp;
  ReadTransformChain(param.reslice_param.transforms, ref_space, warp);

  // Write the composite warp if requested
  if(r_param.out_composed_warp.size())
    {
    LDDMMType::vimg_write(warp.GetPointer(), r_param.out_composed_warp.c_str(),
                          itk::ImageIOBase::FLOAT);
    }

  // Compute the Jacobian of the warp if requested
  if(r_param.out_jacobian_image.size())
    {
    ImagePointer iTemp = ImageType::New();
    LDDMMType::alloc_img(iTemp, warp);
    LDDMMType::field_jacobian_det(warp, iTemp);

    LDDMMType::img_write(iTemp, r_param.out_jacobian_image.c_str(),
                         itk::ImageIOBase::FLOAT);
    }


  // Process image pairs
  for(uint i = 0; i < r_param.images.size(); i++)
    {
    const char *filename = r_param.images[i].moving.c_str();

    // Handle the special case of multi-label images
    if(r_param.images[i].interp.mode == InterpSpec::LABELWISE)
      {
      // The label image assumed to be an image of shorts
      typedef itk::Image<short, VDim> LabelImageType;
      typedef itk::ImageFileReader<LabelImageType> LabelReaderType;

      // Create a reader
      typename LabelReaderType::Pointer reader = LabelReaderType::New();
      reader->SetFileName(filename);
      reader->Update();
      typename LabelImageType::Pointer moving = reader->GetOutput();

      // Scan the unique labels in the image
      std::set<short> label_set;
      short *labels = moving->GetBufferPointer();
      int n_pixels = moving->GetPixelContainer()->Size();

      // Get the list of unique pixels
      short last_pixel = 0;
      for(int j = 0; j < n_pixels; j++)
        {
        short pixel = labels[j];
        if(last_pixel != pixel || i == 0)
          {
          label_set.insert(pixel);
          last_pixel = pixel;
          if(label_set.size() > 1000)
            throw GreedyException("Label wise interpolation not supported for image %s "
                                  "which has over 1000 distinct labels", filename);
          }
        }

      // Turn this set into an array
      std::vector<short> label_array(label_set.begin(), label_set.end());

      // Create a N-way voting filter
      typedef NaryLabelVotingFunctor<ImageType, LabelImageType> VotingFunctor;
      VotingFunctor vf(label_array);

      typedef itk::NaryFunctorImageFilter<ImageType, LabelImageType, VotingFunctor> VotingFilter;
      typename VotingFilter::Pointer fltVoting = VotingFilter::New();
      fltVoting->SetFunctor(vf);

      // Create a mini-pipeline of streaming filters
      for(uint j = 0; j < label_array.size(); j++)
        {
        // Set up a threshold filter for this label
        typedef itk::BinaryThresholdImageFilter<LabelImageType, ImageType> ThresholdFilterType;
        typename ThresholdFilterType::Pointer fltThreshold = ThresholdFilterType::New();
        fltThreshold->SetInput(moving);
        fltThreshold->SetLowerThreshold(label_array[j]);
        fltThreshold->SetUpperThreshold(label_array[j]);
        fltThreshold->SetInsideValue(1.0);
        fltThreshold->SetOutsideValue(0.0);

        // Set up a smoothing filter for this label
        typedef itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType> SmootherType;
        typename SmootherType::Pointer fltSmooth = SmootherType::New();
        fltSmooth->SetInput(fltThreshold->GetOutput());

        // Work out the sigmas for the filter
        if(r_param.images[i].interp.sigma.physical_units)
          {
          fltSmooth->SetSigma(r_param.images[i].interp.sigma.sigma);
          }
        else
          {
          typename SmootherType::SigmaArrayType sigma_array;
          for(uint d = 0; d < VDim; d++)
            sigma_array[d] = r_param.images[i].interp.sigma.sigma * moving->GetSpacing()[d];
          fltSmooth->SetSigmaArray(sigma_array);
          }

        // TODO: we should really be coercing the output into a vector image to speed up interpolation!
        typedef FastWarpCompositeImageFilter<ImageType, ImageType, VectorImageType> InterpFilter;
        typename InterpFilter::Pointer fltInterp = InterpFilter::New();
        fltInterp->SetMovingImage(fltSmooth->GetOutput());
        fltInterp->SetDeformationField(warp);
        fltInterp->SetUsePhysicalSpace(true);

        fltInterp->Update();

        // Add to the voting filter
        fltVoting->SetInput(j, fltInterp->GetOutput());
        }

      // TODO: test out streaming!
      // Run this big pipeline
      fltVoting->Update();

      // Save
      typedef itk::ImageFileWriter<LabelImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetFileName(r_param.images[i].output.c_str());
      writer->SetUseCompression(true);
      writer->SetInput(fltVoting->GetOutput());
      writer->Update();
      }
    else
      {
      // Read the input image
      CompositeImagePointer moving, warped;
      itk::ImageIOBase::IOComponentType comp = LDDMMType::cimg_read(filename, moving);

      // Allocate the warped image
      LDDMMType::alloc_cimg(warped, ref_space, moving->GetNumberOfComponentsPerPixel());

      // Perform the warp
      LDDMMType::interp_cimg(moving, warp, warped,
                             r_param.images[i].interp.mode == InterpSpec::NEAREST,
                             true, r_param.images[i].interp.outside_value);

      // Write, casting to the input component type
      LDDMMType::cimg_write(warped, r_param.images[i].output.c_str(), comp);
      }
    }

  // Process meshes
  for(uint i = 0; i < r_param.meshes.size(); i++)
    {
    typedef itk::Mesh<TReal, VDim> MeshType;
    typedef itk::MeshFileReader<MeshType> MeshReader;
    typename MeshReader::Pointer reader = MeshReader::New();
    reader->SetFileName(r_param.meshes[i].fixed.c_str());

    typedef WarpMeshTransformFunctor<VDim, TReal> TransformType;
    typename TransformType::Pointer transform = TransformType::New();
    transform->SetWarp(warp);
    transform->SetReferenceSpace(ref);

    typedef itk::TransformMeshFilter<MeshType, MeshType, TransformType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetTransform(transform);
    filter->SetInput(reader->GetOutput());

    typedef itk::MeshFileWriter<MeshType> MeshWriter;
    typename MeshWriter::Pointer writer = MeshWriter::New();
    writer->SetInput(filter->GetOutput());
    writer->SetFileName(r_param.meshes[i].output.c_str());
    writer->Update();
    }



  return 0;
}

template <unsigned int VDim, typename TReal>
void GreedyApproach<VDim, TReal>
::ComputeImageMoments(CompositeImageType *image,
                      const std::vector<double> &weights,
                      VecFx &m1, MatFx &m2)
{
  int n = image->GetNumberOfComponentsPerPixel();
  TReal sum_I = 0.0;
  m1.fill(0.0); m2.fill(0.0);

  typedef itk::ImageRegionConstIteratorWithIndex<CompositeImageType> Iterator;
  for(Iterator it(image, image->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    typedef itk::Point<TReal, VDim> PointType;
    PointType p_lps, p_ras;
    image->TransformIndexToPhysicalPoint(it.GetIndex(), p_lps);
    PhysicalCoordinateTransform<VDim, PointType>::lps_to_ras(p_lps, p_ras);
    VecFx X(p_ras.GetDataPointer());
    MatFx X2 = outer_product(X, X);

    typename CompositeImageType::PixelType pix = it.Get();

    // Just weight the components of intensity by weight vector - this sort of makes sense?
    TReal val = 0.0;
    for(int k = 0; k < n; k++)
      val += weights[k] * pix[k];

    sum_I += val;
    m1 += X * val;
    m2 += X2 * val;
    }

  // Compute the mean and covariance from the sum of squares
  m1 = m1 / sum_I;
  m2 = (m2 - sum_I *  outer_product(m1, m1)) / sum_I;
}

template <unsigned int VDim, typename TReal>
int GreedyApproach<VDim, TReal>
::RunAlignMoments(GreedyParameters &param)
{
  // Create an optical flow helper object
  OFHelperType of_helper;

  // No multi-resolution
  of_helper.SetDefaultPyramidFactors(1);

  // Read the image pairs to register
  ReadImages(param, of_helper);

  // Compute the moments of intertia for the fixed and moving images. For now
  // this is done in an iterator loop, out of laziness. Should be converted to
  // a filter if this whole moments business proves useful
  VecFx m1f, m1m;
  MatFx m2f, m2m;


  std::cout << "--- MATCHING BY MOMENTS OF ORDER " << param.moments_order << " ---" << std::endl;

  ComputeImageMoments(of_helper.GetFixedComposite(0), of_helper.GetWeights(), m1f, m2f);

  std::cout << "Fixed Mean        : " << m1f << std::endl;
  std::cout << "Fixed Covariance  : " << std::endl << m2f << std::endl;

  ComputeImageMoments(of_helper.GetMovingComposite(0), of_helper.GetWeights(), m1m, m2m);

  std::cout << "Moving Mean       : " << m1m << std::endl;
  std::cout << "Moving Covariance : " << std::endl << m2m << std::endl;

  // This flag forces no rotation, only flip
  if(param.moments_order == 1 || param.flag_moments_id_covariance)
    {
    m2f.set_identity();
    m2m.set_identity();
    }

  // Decompose covariance matrices into eigenvectors and eigenvalues
  vnl_vector<TReal> Df, Dm;
  vnl_matrix<TReal> Vf, Vm;
  vnl_symmetric_eigensystem_compute<TReal>(m2f, Vf, Df);
  vnl_symmetric_eigensystem_compute<TReal>(m2m, Vm, Dm);

  // Create a rigid registration problem
  PhysicalSpaceAffineCostFunction cost_fn(&param, this, 0, &of_helper);

  // The best set of coefficients and the associated match value
  vnl_vector<double> xBest;
  TReal xBestMatch = vnl_numeric_traits<TReal>::maxval;

  // Generate all possible flip matrices
  int n_flip = 1 << VDim;
  for(int k_flip = 0; k_flip < n_flip; k_flip++)
    {
    // If using first moments, ignore all flips, only allow identity
    if(param.moments_order == 1 && k_flip != n_flip - 1)
      continue;

    // Generate the flip matrix
    MatFx F(0.0);
    for(uint d = 0; d < VDim; d++)
      F(d,d) = (k_flip & (1 << d)) ? 1 : -1;;

    // Compute the rotation matrix - takes fixed coordinates into moving space
    MatFx R = Vm * F * Vf.transpose();
    VecFx b = m1m - R * m1f;

    vnl_matrix<TReal> A(VDim+1, VDim+1, 0.0);
    A.set_identity();
    A.update(R, 0, 0);
    for(uint d= 0 ;d< VDim;d++)
      A(d,VDim) = b[d];

    // Ignore flips with the wrong determinant
    double det_R = vnl_determinant(R);
    if((param.moments_order == 2 && param.moments_flip_determinant == 1 && det_R < 0) ||
       (param.moments_order == 2 && param.moments_flip_determinant == -1 && det_R > 0))
      {
      continue;
      }

    // Generate affine coefficients from the rotation and shift
    vnl_vector<double> x(cost_fn.get_number_of_unknowns());
    flatten_affine_transform(R, b, x.data_block());

    // Compute similarity
    double f = 0.0;
    cost_fn.compute(x, &f, NULL);

    std::cout << "Metric for flip " << F.get_diagonal() << " : " << f << std::endl;

    // Compare
    if(xBestMatch > f || xBest.size() == 0)
      {
      xBestMatch = f;
      xBest = x;
      }
    }

  // Save the best transform
  typename LinearTransformType::Pointer tran = LinearTransformType::New();
  cost_fn.GetTransform(xBest, tran);
  vnl_matrix<double> Q_physical = MapAffineToPhysicalRASSpace(of_helper, 0, tran);
  this->WriteAffineMatrixViaCache(param.output, Q_physical);

  return 0;
}

/**
 * Post-hoc warp inversion - the Achilles heel of non-symmetric registration :(
 */
template <unsigned int VDim, typename TReal>
int GreedyApproach<VDim, TReal>
::RunInvertWarp(GreedyParameters &param)
{
  // Read the warp as a transform chain
  VectorImagePointer warp;

  // Read the warp file
  LDDMMType::vimg_read(param.invwarp_param.in_warp.c_str(), warp);

  // Convert the warp file into voxel units from physical units
  OFHelperType::PhysicalWarpToVoxelWarp(warp, warp, warp);


  // Compute the inverse of the warp
  VectorImagePointer uInverse = VectorImageType::New();
  LDDMMType::alloc_vimg(uInverse, warp);
  OFHelperType::ComputeDeformationFieldInverse(warp, uInverse, param.warp_exponent, true);

  // Write the warp using compressed format
  OFHelperType::WriteCompressedWarpInPhysicalSpace(uInverse, warp, param.invwarp_param.out_warp.c_str(), param.warp_precision);

  return 0;
}

/**
 * Post-hoc warp root
 */
template <unsigned int VDim, typename TReal>
int GreedyApproach<VDim, TReal>
::RunRootWarp(GreedyParameters &param)
{
  // Read the warp as a transform chain
  VectorImagePointer warp;

  // Read the warp file
  LDDMMType::vimg_read(param.warproot_param.in_warp.c_str(), warp);

  // Convert the warp file into voxel units from physical units
  OFHelperType::PhysicalWarpToVoxelWarp(warp, warp, warp);

  // Allocate the root
  VectorImagePointer warp_root = VectorImageType::New();
  LDDMMType::alloc_vimg(warp_root, warp);

  // Take the n-th root
  OFHelperType::ComputeWarpRoot(warp, warp_root, param.warp_exponent, 1e-6);

  // Write the warp using compressed format
  OFHelperType::WriteCompressedWarpInPhysicalSpace(warp_root, warp, param.warproot_param.out_warp.c_str(), param.warp_precision);

  return 0;
}


template <unsigned int VDim, typename TReal>
void GreedyApproach<VDim, TReal>
::AddCachedInputObject(std::string &string, itk::Object *object)
{
  m_ImageCache[string] = object;
}

template <unsigned int VDim, typename TReal>
const typename GreedyApproach<VDim,TReal>::MetricLogType &
GreedyApproach<VDim,TReal>
::GetMetricLog() const
{
  return m_MetricLog;
}

template <unsigned int VDim, typename TReal>
int GreedyApproach<VDim, TReal>
::Run(GreedyParameters &param)
{
  switch(param.mode)
    {
    case GreedyParameters::GREEDY:
      return Self::RunDeformable(param);
    case GreedyParameters::AFFINE:
      return Self::RunAffine(param);
    case GreedyParameters::BRUTE:
      return Self::RunBrute(param);
    case GreedyParameters::MOMENTS:
      return Self::RunAlignMoments(param);
    case GreedyParameters::RESLICE:
      return Self::RunReslice(param);
    case GreedyParameters::INVERT_WARP:
      return Self::RunInvertWarp(param);
    case GreedyParameters::JACOBIAN_WARP:
      return Self::RunJacobian(param);
    case GreedyParameters::ROOT_WARP:
      return Self::RunRootWarp(param);
    }

  return -1;
}




template class GreedyApproach<2, float>;
template class GreedyApproach<3, float>;
template class GreedyApproach<4, float>;
template class GreedyApproach<2, double>;
template class GreedyApproach<3, double>;
template class GreedyApproach<4, double>;
