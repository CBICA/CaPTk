/*=========================================================================

  Program:   ALFABIS fast medical image registration programs
  Language:  C++
  Website:   github.com/pyushkevich/greedy
  Copyright (c) Paul Yushkevich, University of Pennsylvania. All rights reserved.

  This program is part of ALFABIS: Adaptive Large-Scale Framework for
  Automatic Biomedical Image Segmentation.

  ALFABIS development is funded by the NIH grant R01 EB017255.

  ALFABIS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as publishGed by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ALFABIS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ALFABIS.  If not, see <http://www.gnu.org/licenses/>.

=========================================================================*/
#include "lddmm_data.h"
#include "itkImageRegionIterator.h"
#include "SimpleWarpImageFilter.h"
#include "itkNumericTraitsCovariantVectorPixel.h"
#include "itkOptVectorLinearInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkAddImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkGradientImageFilter.h"
#include "itkUnaryFunctorImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVectorImage.h"
#include "itkDisplacementFieldJacobianDeterminantFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkVectorResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkVectorNearestNeighborInterpolateImageFunction.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkComposeImageFilter.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkTernaryFunctorImageFilter.h"
#include "itkShiftScaleImageFilter.h"

#include "FastWarpCompositeImageFilter.h"

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::alloc_vf(VelocityField &vf, uint nt, ImageBaseType *ref)
{
  vf.resize(nt);
  for(uint i = 0; i < nt; i++)
    alloc_vimg(vf[i], ref);
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::alloc_vimg(VectorImagePointer &img, ImageBaseType *ref)
{
  img = VectorImageType::New();
  img->SetRegions(ref->GetBufferedRegion());
  img->CopyInformation(ref);
  img->Allocate();
  img->FillBuffer(Vec(0.0));
}

template <class TFloat, uint VDim>
typename LDDMMData<TFloat, VDim>::VectorImagePointer
LDDMMData<TFloat, VDim>
::alloc_vimg(ImageBaseType *ref)
{
  VectorImagePointer p = VectorImageType::New();
  alloc_vimg(p, ref);
  return p;
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::alloc_mimg(MatrixImagePointer &img, ImageBaseType *ref)
{
  img = MatrixImageType::New();
  img->SetRegions(ref->GetBufferedRegion());
  img->CopyInformation(ref);
  img->Allocate();
  img->FillBuffer(Mat());
}

template <class TFloat, uint VDim>
typename LDDMMData<TFloat, VDim>::MatrixImagePointer
LDDMMData<TFloat, VDim>
::alloc_mimg(ImageBaseType *ref)
{
  MatrixImagePointer p = MatrixImageType::New();
  alloc_mimg(p, ref);
  return p;
}

template <class TFloat, uint VDim>
void
LDDMMData<TFloat, VDim>
::alloc_cimg(CompositeImagePointer &img, ImageBaseType *ref, int n_comp)
{
  img = CompositeImageType::New();
  img->SetRegions(ref->GetBufferedRegion());
  img->CopyInformation(ref);
  img->SetNumberOfComponentsPerPixel(n_comp);
  img->Allocate();

  typename CompositeImageType::PixelType cpix;
  cpix.SetSize(n_comp);
  cpix.Fill(0.0);
  img->FillBuffer(cpix);
}

template <class TFloat, uint VDim>
typename LDDMMData<TFloat, VDim>::CompositeImagePointer
LDDMMData<TFloat, VDim>
::alloc_cimg(ImageBaseType *ref, int n_comp)
{
  CompositeImagePointer p;
  alloc_cimg(p, ref, n_comp);
  return p;
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::alloc_img(ImagePointer &img, ImageBaseType *ref)
{
  img = ImageType::New();
  img->SetRegions(ref->GetBufferedRegion());
  img->CopyInformation(ref);
  img->Allocate();
  img->FillBuffer(0.0);
}

template <class TFloat, uint VDim>
typename LDDMMData<TFloat, VDim>::ImagePointer
LDDMMData<TFloat, VDim>
::alloc_img(ImageBaseType *ref)
{
  ImagePointer p;
  alloc_img(p, ref);
  return p;
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::init(LDDMMData<TFloat, VDim> &p, 
  ImageType *fix, ImageType *mov, 
  uint nt, double alpha, double gamma, double sigma)
{
  p.fix = fix;
  p.mov = mov;
  p.alpha = alpha;
  p.sigma = sigma;
  p.gamma = gamma;
  p.nt = nt;
  p.dt = 1.0 / (nt - 1.0);
  p.sigma_sq = sigma * sigma;

  // Initialize N and R
  p.r = fix->GetBufferedRegion();
  p.nv = fix->GetBufferedRegion().GetNumberOfPixels();
  for(uint i = 0; i < VDim; i++)
    p.n[i] = p.r.GetSize()[i];

  // Initialize the velocity fields
  alloc_vf(p.v, nt, fix);
  alloc_vf(p.a, nt, fix);
  alloc_vf(p.f, nt, fix);

  // Initialize kernel terms
  alloc_img(p.f_kernel, fix);
  alloc_img(p.f_kernel_sq, fix);

  // Compute these images
  ImageIterator it(p.f_kernel, p.r), itsq(p.f_kernel_sq, p.r);
  for(; !it.IsAtEnd(); ++it, ++itsq)
    {
    TFloat val = 0.0;
    for(uint j = 0; j < VDim; j++)
      val += 1.0 - cos(it.GetIndex()[j] * 2.0 * vnl_math::pi / p.n[j]);
    it.Set(p.gamma + 2.0 * p.alpha * p.nv * val);
    itsq.Set(it.Get() * it.Get());
    }

  // Initialize temporary vector field
  alloc_vimg(p.vtmp, fix);
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::compute_navier_stokes_kernel(ImageType *kernel, double alpha, double gamma)
{
  ImageIterator it(kernel, kernel->GetBufferedRegion());
  itk::Size<VDim> sz = kernel->GetBufferedRegion().GetSize();
  double alpha_scale = 2.0 * alpha * kernel->GetBufferedRegion().GetNumberOfPixels();

  for(; !it.IsAtEnd(); ++it)
    {
    TFloat val = 0.0;
    for(uint j = 0; j < VDim; j++)
      val += 1.0 - cos(it.GetIndex()[j] * 2.0 * vnl_math::pi / sz[j]);
    double k = gamma + alpha_scale * val;
    it.Set(k * k);
    }
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::interp_vimg(VectorImageType *data, VectorImageType *field,
  TFloat def_scale, VectorImageType *out, bool use_nn, bool phys_space)
{
  typedef FastWarpCompositeImageFilter<VectorImageType, VectorImageType, VectorImageType> WF;
  typename WF::Pointer wf = WF::New();
  wf->SetDeformationField(field);
  wf->SetMovingImage(data);
  wf->GraftOutput(out);
  wf->SetDeformationScaling(def_scale);
  wf->SetUseNearestNeighbor(use_nn);
  wf->SetUsePhysicalSpace(phys_space);
  wf->Update();
}

template <class TFloat, uint VDim>
void
LDDMMData<TFloat, VDim>
::vimg_exp(
  const VectorImageType *src, VectorImageType *trg, VectorImageType *work,
  int exponent, TFloat scale)
{
  // Scale the image if needed
  if(scale != 1.0)
    vimg_scale(src, scale, trg);
  else
    vimg_copy(src, trg);

  for(int q = 0; q < exponent; q++)
    {
    interp_vimg(trg, trg, 1.0, work);
    vimg_add_in_place(trg, work);
    }
}

template <class TFloat, uint VDim>
void
LDDMMData<TFloat, VDim>
::vimg_exp_with_jacobian(
  const VectorImageType *src, VectorImageType *trg, VectorImageType *work,
  MatrixImageType *trg_jac, MatrixImageType *work_mat,
  int exponent, TFloat scale)
{
  // Scale the image if needed
  if(scale != 1.0)
    vimg_scale(src, scale, trg);
  else
    vimg_copy(src, trg);

  // Compute the initial Jacobian
  field_jacobian(trg, trg_jac);

  // Perform the exponentiation
  for(int q = 0; q < exponent; q++)
    {
    // Compute the composition of the Jacobian with itself, place in jac_work
    jacobian_of_composition(trg_jac, trg_jac, trg, work_mat);

    // Copy the data (TODO: this is a little wasteful)
    mimg_copy(work_mat, trg_jac);

    // Update the velocity field
    interp_vimg(trg, trg, 1.0, work);
    vimg_add_in_place(trg, work);
    }
}


template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::interp_mimg(MatrixImageType *data, VectorImageType *field,
  MatrixImageType *out, bool use_nn, bool phys_space)
{
  // Decorate the matrix images as multi-component images
  CompositeImagePointer wrap_data = CompositeImageType::New();
  wrap_data->SetRegions(data->GetBufferedRegion());
  wrap_data->CopyInformation(data);
  wrap_data->SetNumberOfComponentsPerPixel(VDim * VDim);
  wrap_data->GetPixelContainer()->SetImportPointer(
    (TFloat *)(data->GetPixelContainer()->GetImportPointer()),
    VDim * VDim * data->GetPixelContainer()->Size(), false);

  // Decorate the output image in the same way
  CompositeImagePointer wrap_out = CompositeImageType::New();
  wrap_out->SetRegions(out->GetBufferedRegion());
  wrap_out->CopyInformation(out);
  wrap_out->SetNumberOfComponentsPerPixel(VDim * VDim);
  wrap_out->GetPixelContainer()->SetImportPointer(
    (TFloat *)(out->GetPixelContainer()->GetImportPointer()),
    VDim * VDim * out->GetPixelContainer()->Size(), false);

  // Perform the interpolation
  LDDMMData<TFloat, VDim>::interp_cimg(wrap_data, field, wrap_out, use_nn, phys_space);
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::interp_img(ImageType *data, VectorImageType *field, ImageType *out, 
  bool use_nn, bool phys_space, TFloat outside_value)
{
  typedef FastWarpCompositeImageFilter<ImageType, ImageType, VectorImageType> WF;
  typename WF::Pointer wf = WF::New();
  wf->SetDeformationField(field);
  wf->SetMovingImage(data);
  wf->GraftOutput(out);
  wf->SetUseNearestNeighbor(use_nn);
  wf->SetUsePhysicalSpace(phys_space);
  wf->SetOutsideValue(outside_value);
  wf->Update();

/*
  // Create a warp filter
  typedef itk::SimpleWarpImageFilter<
    ImageType, ImageType, VectorImageType, TFloat> WarpFilterType;
  typename WarpFilterType::Pointer flt = WarpFilterType::New();

  // Create an interpolation function
  typedef itk::LinearInterpolateImageFunction<ImageType, TFloat> InterpType;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, TFloat> NNInterpType;

  typename InterpType::Pointer func = InterpType::New();
  typename NNInterpType::Pointer funcNN = NNInterpType::New();

  // Graft output of the warp filter
  flt->GraftOutput(out);

  // Set inputs
  flt->SetInput(data);
  if(use_nn)
    flt->SetInterpolator(funcNN);
  else
    flt->SetInterpolator(func);
  flt->SetDeformationField(field);
  flt->SetDeformationScaling(1.0);
  flt->Update();
  */
}

template <class TFloat, uint VDim>
void
LDDMMData<TFloat, VDim>
::interp_cimg(CompositeImageType *data, VectorImageType *field, CompositeImageType *out,
              bool use_nn, bool phys_space, TFloat outside_value)
{
  typedef FastWarpCompositeImageFilter<CompositeImageType, CompositeImageType, VectorImageType> WF;
  typename WF::Pointer wf = WF::New();
  wf->SetDeformationField(field);
  wf->SetMovingImage(data);
  wf->GraftOutput(out);
  wf->SetUseNearestNeighbor(use_nn);
  wf->SetUsePhysicalSpace(phys_space);
  wf->SetOutsideValue(outside_value);
  wf->Update();
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::vimg_add_in_place(VectorImageType *trg, VectorImageType *a)
{
  typedef itk::AddImageFilter<VectorImageType> AddFilter;
  typename AddFilter::Pointer flt = AddFilter::New();
  flt->SetInput(0,trg);
  flt->SetInput(1,a);
  flt->GraftOutput(trg);
  flt->Update();
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::vimg_subtract_in_place(VectorImageType *trg, VectorImageType *a)
{
  typedef itk::SubtractImageFilter<VectorImageType> SubtractFilter;
  typename SubtractFilter::Pointer flt = SubtractFilter::New();
  flt->SetInput(0,trg);
  flt->SetInput(1,a);
  flt->GraftOutput(trg);
  flt->Update();
}

// Scalar math
template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::vimg_multiply_in_place(VectorImageType *trg, ImageType *s)
{
  typedef itk::MultiplyImageFilter<
    VectorImageType, ImageType, VectorImageType> MultiplyFilter;
  typename MultiplyFilter::Pointer flt = MultiplyFilter::New();
  flt->SetInput1(trg);
  flt->SetInput2(s);
  flt->GraftOutput(trg);
  flt->Update();
}

template <class TFloat, uint VDim>
void
LDDMMData<TFloat, VDim>
::img_scale_in_place(ImageType *img, TFloat scale)
{
  typedef itk::ShiftScaleImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetScale(scale);
  filter->SetInput(img);
  filter->GraftOutput(img);
  filter->Update();
}


template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::mimg_multiply_in_place(MatrixImageType *trg, MatrixImageType *s)
{
  typedef itk::MultiplyImageFilter<
    MatrixImageType, MatrixImageType, MatrixImageType> MultiplyFilter;
  typename MultiplyFilter::Pointer flt = MultiplyFilter::New();
  flt->SetInput1(trg);
  flt->SetInput2(s);
  flt->GraftOutput(trg);
  flt->Update();
}

// Scalar math

template <class TFloat, uint VDim>
void
LDDMMData<TFloat, VDim>
::img_add_in_place(ImagePointer &trg, ImageType *a)
{
  typedef itk::AddImageFilter<ImageType> AddFilter;
  typename AddFilter::Pointer flt = AddFilter::New();
  flt->SetInput(0,trg);
  flt->SetInput(1,a);
  flt->GraftOutput(trg);
  flt->Update();
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::img_subtract_in_place(ImagePointer &trg, ImageType *a)
{
  typedef itk::SubtractImageFilter<ImageType> SubtractFilter;
  typename SubtractFilter::Pointer flt = SubtractFilter::New();
  flt->SetInput(0,trg);
  flt->SetInput(1,a);
  flt->GraftOutput(trg);
  flt->Update();
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::img_multiply_in_place(ImagePointer &trg, ImageType *a)
{
  typedef itk::MultiplyImageFilter<ImageType> MultiplyFilter;
  typename MultiplyFilter::Pointer flt = MultiplyFilter::New();
  flt->SetInput(0,trg);
  flt->SetInput(1,a);
  flt->GraftOutput(trg);
  flt->Update();
}

template <class TFloat, uint VDim>
TFloat 
LDDMMData<TFloat, VDim>
::vimg_euclidean_norm_sq(VectorImageType *trg)
{
  // Add all voxels in the image
  double accum = 0.0;
  typedef itk::ImageRegionIterator<VectorImageType> Iter;
  for(Iter it(trg, trg->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    for(uint d = 0; d < VDim; d++)
      accum += it.Value()[d] * it.Value()[d];
    }
  return (TFloat) accum;
}

template <class TFloat, uint VDim>
TFloat 
LDDMMData<TFloat, VDim>
::img_euclidean_norm_sq(ImageType *trg)
{
  // Add all voxels in the image
  double accum = 0.0;
  typedef itk::ImageRegionIterator<ImageType> Iter;
  for(Iter it(trg, trg->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    { accum += it.Value() * it.Value(); }
  return (TFloat) accum;
}


template <class TFloat, uint VDim>
class LinearToConstRectifierFunctor
{
public:
  typedef LinearToConstRectifierFunctor<TFloat, VDim> Self;

  LinearToConstRectifierFunctor() : m_Threshold(0.0), m_Offset(0.0) {}
  LinearToConstRectifierFunctor(TFloat thresh) : m_Threshold(thresh)
    { m_Offset = log(1 + exp(thresh)); }

  TFloat operator() (const TFloat &x)
    { return m_Offset - log(1 + exp(m_Threshold - x)); }

  bool operator== (const Self &other)
    { return m_Threshold == other.m_Threshold; }

  bool operator!= (const Self &other)
    { return m_Threshold != other.m_Threshold; }

protected:
  TFloat m_Threshold, m_Offset;
};

template <class TFloat, uint VDim>
void
LDDMMData<TFloat, VDim>
::img_linear_to_const_rectifier_fn(ImageType *src, ImageType *trg, TFloat thresh)
{
  typedef LinearToConstRectifierFunctor<TFloat, VDim> Functor;
  typedef itk::UnaryFunctorImageFilter<ImageType, ImageType, Functor> Filter;
  typename Filter::Pointer flt = Filter::New();

  Functor func(thresh);
  flt->SetFunctor(func);
  flt->SetInput(trg);
  flt->GraftOutput(trg);
  flt->Update();
}

template <class TFloat, uint VDim>
class LinearToConstRectifierDerivFunctor
{
public:
  typedef LinearToConstRectifierDerivFunctor<TFloat, VDim> Self;

  LinearToConstRectifierDerivFunctor() : m_Threshold(0.0) {}
  LinearToConstRectifierDerivFunctor(TFloat thresh) : m_Threshold(thresh) {}

  TFloat operator() (const TFloat &x)
    { return 1.0 / (1.0 + exp(x - m_Threshold)); }

  bool operator== (const Self &other)
    { return m_Threshold == other.m_Threshold; }

  bool operator!= (const Self &other)
    { return m_Threshold != other.m_Threshold; }

protected:
  TFloat m_Threshold;
};

template <class TFloat, uint VDim>
void
LDDMMData<TFloat, VDim>
::img_linear_to_const_rectifier_deriv(ImageType *src, ImageType *trg, TFloat thresh)
{
  typedef LinearToConstRectifierDerivFunctor<TFloat, VDim> Functor;
  typedef itk::UnaryFunctorImageFilter<ImageType, ImageType, Functor> Filter;
  typename Filter::Pointer flt = Filter::New();

  Functor func(thresh);
  flt->SetFunctor(func);
  flt->SetInput(trg);
  flt->GraftOutput(trg);
  flt->Update();
}

template <class TFloat, uint VDim>
TFloat 
LDDMMData<TFloat, VDim>
::img_voxel_sum(ImageType *trg)
{
  // Add all voxels in the image
  double accum = 0.0;
  typedef itk::ImageRegionIterator<ImageType> Iter;
  for(Iter it(trg, trg->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    accum += it.Value();
  return (TFloat) accum;
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::img_min_max(ImageType *src, TFloat &out_min, TFloat &out_max)
{
  // Add all voxels in the image
  typedef itk::MinimumMaximumImageFilter<ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(src);
  filter->Update();
  out_min = filter->GetMinimum();
  out_max = filter->GetMaximum();
}


template <class TFloat, uint VDim>
class VectorScaleFunctor
{
public:
  VectorScaleFunctor() { this->Scale = 1.0; }
  typedef itk::CovariantVector<TFloat,VDim> Vec;

  Vec operator() (const Vec &x)
    { return x * Scale; }

  bool operator== (const VectorScaleFunctor<TFloat, VDim> &other)
    { return Scale == other.Scale; }

  bool operator!= (const VectorScaleFunctor<TFloat, VDim> &other)
    { return Scale != other.Scale; }

  TFloat Scale;
};

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::vimg_scale_in_place(VectorImageType *trg, TFloat s)
{
  typedef VectorScaleFunctor<TFloat, VDim> Functor;
  typedef itk::UnaryFunctorImageFilter<
    VectorImageType, VectorImageType, Functor> Filter;
  typename Filter::Pointer flt = Filter::New();

  Functor func;
  func.Scale = s;
  flt->SetFunctor(func);
  flt->SetInput(trg);
  flt->GraftOutput(trg);
  flt->Update();
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::vimg_scale(const VectorImageType*src, TFloat s, VectorImageType *trg)
{
  typedef VectorScaleFunctor<TFloat, VDim> Functor;
  typedef itk::UnaryFunctorImageFilter<
    VectorImageType, VectorImageType, Functor> Filter;
  typename Filter::Pointer flt = Filter::New();

  Functor func;
  func.Scale = s;
  flt->SetFunctor(func);
  flt->SetInput(src);
  flt->GraftOutput(trg);
  flt->Update();
}

template <class TFloat, uint VDim>
class VectorScaleAddFunctor
{
public:
  VectorScaleAddFunctor() { this->Scale = 1.0; }
  typedef itk::CovariantVector<TFloat,VDim> Vec;

  Vec operator() (const Vec &x, const Vec &y)
    { return x + y * Scale; }

  bool operator== (const VectorScaleAddFunctor<TFloat, VDim> &other)
    { return Scale == other.Scale; }

  bool operator!= (const VectorScaleAddFunctor<TFloat, VDim> &other)
    { return Scale != other.Scale; }

  TFloat Scale;
};

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::vimg_add_scaled_in_place(VectorImageType *trg, VectorImageType *a, TFloat s)
{
  typedef VectorScaleAddFunctor<TFloat, VDim> Functor;
  typedef itk::BinaryFunctorImageFilter<
    VectorImageType, VectorImageType, VectorImageType, Functor> Filter;
  typename Filter::Pointer flt = Filter::New();

  Functor func;
  func.Scale = s;
  flt->SetFunctor(func);
  flt->SetInput1(trg);
  flt->SetInput2(a);
  flt->GraftOutput(trg);
  flt->Update();
}

template <class TFloat, uint VDim>
class VectorDotProduct
{
public:
  VectorDotProduct() {}
  typedef itk::CovariantVector<TFloat,VDim> Vec;

  TFloat operator() (const Vec &x, const Vec &y)
    {
    TFloat dp = 0.0;
    for(uint d = 0; d < VDim; d++)
      dp += x[d] * y[d];
    return dp;
    }

  bool operator== (const VectorDotProduct<TFloat, VDim> &other)
    { return true; }

  bool operator!= (const VectorDotProduct<TFloat, VDim> &other)
    { return false; }
};

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::vimg_euclidean_inner_product(ImagePointer &trg, VectorImageType *a, VectorImageType *b)
{
  typedef VectorDotProduct<TFloat, VDim> Functor;
  typedef itk::BinaryFunctorImageFilter<
    VectorImageType, VectorImageType, ImageType, Functor> Filter;
  typename Filter::Pointer flt = Filter::New();

  Functor func;
  flt->SetFunctor(func);
  flt->SetInput1(a);
  flt->SetInput2(b);
  flt->GraftOutput(trg);
  flt->Update();
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::compute_semi_lagrangean_a()
{
  for(uint i = 0; i < nt; i++)
    {
    a[i]->FillBuffer(Vec(0.0));
    for (uint j = 0; j < 5; j++)
      {
      interp_vimg(v[i], a[i], -0.5, a[i]);
      vimg_scale_in_place(a[i], dt);
      itk::Index<VDim> x;
      x[0] = 63; x[1] = 63;
      }
    }

}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::integrate_phi_t0()
{
  // Go through and compute phi_t0 for each time point
  for(uint m = 0; m < nt; m++)
    if(m==0)
      {
      f[m]->FillBuffer(Vec(0.0));
      }
    else
      {
      interp_vimg(f[m-1], a[m], -1.0, f[m]);
      vimg_subtract_in_place(f[m], a[m]);
      }
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::integrate_phi_t1()
{
  for(uint m = nt-1; m >= 0; m--)
    {
    if(m == nt-1)
      {
      f[m]->FillBuffer(Vec(0.0));
      }
    else
      {
      interp_vimg(f[m+1], a[m], 1.0, f[m]);
      vimg_add_in_place(f[m], a[m]);
      }
    } 
}

template <class TFloat, uint VDim>
class SetMatrixRowBinaryOperator
{
public:
  typedef SetMatrixRowBinaryOperator<TFloat, VDim> Self;
  typedef LDDMMData<TFloat, VDim> LDDMM;
  typedef typename LDDMM::Vec Vec;
  typedef typename LDDMM::Mat Mat;

  SetMatrixRowBinaryOperator() { m_Row = 0; }
  void SetRow(unsigned int row) { m_Row = row; }
  bool operator != (const Self &other) const { return m_Row != other.m_Row; }

  Mat operator() (const Mat &M, const Vec &V)
    {
    Mat Q;
    for(uint r = 0; r < VDim; r++)
      for(uint c = 0; c < VDim; c++)
        Q(r, c) = (r == m_Row) ? V[c] : M(r,c);
    return Q;
    }

protected:
  unsigned int m_Row;
};

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::field_jacobian(VectorImageType *vec, MatrixImageType *out)
{
  for(uint a = 0; a < VDim; a++)
    {
    // Extract the a'th component of the displacement field
    typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType, ImageType> CompFilterType;
    typename CompFilterType::Pointer comp = CompFilterType::New();
    comp->SetIndex(a);
    comp->SetInput(vec);

    // Compute the gradient of this component
    typedef itk::GradientImageFilter<ImageType, TFloat, TFloat> GradientFilter;
    typename GradientFilter::Pointer grad = GradientFilter::New();
    grad->SetInput(comp->GetOutput());
    grad->SetUseImageSpacingOff();
    grad->SetUseImageDirection(false);

    // Apply to the Jacobian matrix
    typedef SetMatrixRowBinaryOperator<TFloat, VDim> RowOperatorType;
    RowOperatorType rop;
    rop.SetRow(a);
    typedef itk::BinaryFunctorImageFilter<
      MatrixImageType, VectorImageType, MatrixImageType, RowOperatorType> RowFilterType;
    typename RowFilterType::Pointer rof = RowFilterType::New();
    rof->SetInput1(out);
    rof->SetInput2(grad->GetOutput());
    rof->SetFunctor(rop);
    rof->GraftOutput(out);
    rof->Update();
    }
}

template <class TFloat, uint VDim>
class JacobianCompisitionFunctor
{
public:
  typedef typename LDDMMData<TFloat, VDim>::Mat Mat;

  Mat operator() (const Mat &Du_wrp, const Mat &Dv)
    {
    Mat Dw = Dv + Du_wrp * Dv + Du_wrp;
    return Dw;
    }

  bool operator != (const JacobianCompisitionFunctor<TFloat, VDim> &other) {return false; }
};


template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::jacobian_of_composition(
    MatrixImageType *Du, MatrixImageType *Dv, VectorImageType *v, MatrixImageType *out_Dw)
{
  // Interpolate jac_phi by psi and place it into out
  interp_mimg(Du, v, out_Dw);

  // Perform the matrix multiplication and addition
  typedef JacobianCompisitionFunctor<TFloat, VDim> Functor;
  typedef itk::BinaryFunctorImageFilter<MatrixImageType,MatrixImageType,MatrixImageType,Functor> BinaryFilter;
  typename BinaryFilter::Pointer flt = BinaryFilter::New();
  flt->SetInput1(out_Dw);
  flt->SetInput2(Dv);
  flt->GraftOutput(out_Dw);
  flt->Update();
}



template <class TFloat, uint VDim>
class MatrixPlusConstDeterminantFunctor
{
public:
  typedef typename LDDMMData<TFloat, VDim>::Mat Mat;

  TFloat operator() (const Mat &M)
    {
    Mat X = m_LambdaEye;
    X += M;
    return vnl_determinant(X.GetVnlMatrix());
    }

  void SetLambda(TFloat lambda)
    {
    m_LambdaEye.SetIdentity();
    m_LambdaEye *= lambda;
    }

  bool operator != (const MatrixPlusConstDeterminantFunctor<TFloat, VDim> &other) 
    { return m_LambdaEye(0,0) != other.m_LambdaEye(0,0); }

protected:
  Mat m_LambdaEye;

};


template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::mimg_det(MatrixImageType *M, double lambda, ImageType *out_det)
{
  typedef MatrixPlusConstDeterminantFunctor<TFloat, VDim> FunctorType;
  FunctorType functor;
  functor.SetLambda(lambda);
  typedef itk::UnaryFunctorImageFilter<MatrixImageType, ImageType, FunctorType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(M);
  filter->SetFunctor(functor);
  filter->GraftOutput(out_det);
  filter->Update();
}

/**
 * Functor to compute Ax+b
 */
template <class TFloat, uint VDim>
class MatrixVectorMultiplyAndAddVectorFunctor
{
public:
  typedef MatrixVectorMultiplyAndAddVectorFunctor<TFloat, VDim> Self;
  typedef LDDMMData<TFloat, VDim> LDDMMType;
  typedef typename LDDMMType::Mat Mat;
  typedef typename LDDMMType::Vec Vec;

  Vec operator() (const Mat &A, const Vec &x, const Vec &b)
    {
    Vec y = m_Lambda * (A * x) + m_Mu * b;
    return y;
    }

  void SetLambda(TFloat lambda) { m_Lambda = lambda; }
  void SetMu(TFloat mu) { m_Mu = mu; }

  bool operator != (const Self &other) const
    { return m_Lambda != other.m_Lambda || m_Mu != other.m_Mu; }

  bool operator == (const Self &other) const
    { return ! (*this != other); }

protected:
  TFloat m_Lambda, m_Mu;
};

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::mimg_vimg_product_plus_vimg(
    MatrixImageType *A, VectorImageType *x, VectorImageType *b, 
    TFloat lambda, TFloat mu, VectorImageType *out)
{
  typedef MatrixVectorMultiplyAndAddVectorFunctor<TFloat, VDim> Functor;
  Functor functor;
  functor.SetLambda(lambda);
  functor.SetMu(mu);

  typedef itk::TernaryFunctorImageFilter<
    MatrixImageType, VectorImageType, VectorImageType, VectorImageType, 
    Functor> FilterType;

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1(A);
  filter->SetInput2(x);
  filter->SetInput3(b);
  filter->SetFunctor(functor);
  filter->GraftOutput(out);
  filter->Update();
}

#include "LieBracketFilter.h"

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::lie_bracket(VectorImageType *v, VectorImageType *u, MatrixImageType *work, VectorImageType *out)
{
  // Compute Du, place it in work
  field_jacobian(v, work);

  // Multiply by v
  mimg_vimg_product_plus_vimg(work, u, out, 1.0, 0.0, out);

  // Compute Dv, place it in work
  field_jacobian(u, work);

  // Multiply by u and subtract from existing
  mimg_vimg_product_plus_vimg(work, v, out, -1.0, 1.0, out);

  // Alternative approach
  VectorImagePointer alt = alloc_vimg(out);
  
  typedef LieBracketFilter<VectorImageType, VectorImageType> LieBracketFilterType;
  typename LieBracketFilterType::Pointer fltLieBracket = LieBracketFilterType::New();
  fltLieBracket->SetFieldU(v);
  fltLieBracket->SetFieldV(u);
  fltLieBracket->GraftOutput(alt);
  fltLieBracket->Update();

  itk::Index<VDim> idx_probe; for(uint a = 0; a < VDim; a++) idx_probe[a] = out->GetBufferedRegion().GetSize()[a] / 2;
  Vec test1 = out->GetPixel(idx_probe);
  Vec test2 = alt->GetPixel(idx_probe);
  std::cout << "test1 = " << test1 << "   and   test2 = " << test2 << std::endl;
  return;
}


template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::field_jacobian_det(VectorImageType *vec, ImageType *out)
{
  typedef itk::DisplacementFieldJacobianDeterminantFilter<
    VectorImageType, TFloat, ImageType> Filter;
  typename Filter::Pointer filter = Filter::New();
  filter->SetInput(vec);
  filter->SetUseImageSpacingOff();
  filter->GraftOutput(out);
  filter->Update();
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::image_gradient(ImageType *src, VectorImageType *grad)
{
  // Create a gradient image filter
  typedef itk::GradientImageFilter<ImageType, TFloat, TFloat> Filter;
  typename Filter::Pointer flt = Filter::New();
  flt->SetInput(src);
  flt->GraftOutput(grad);
  flt->SetUseImageSpacingOff();
  flt->SetUseImageDirection(false);
  flt->Update();
}

template <class TFloat, uint VDim>
void
LDDMMData<TFloat, VDim>
::img_smooth(ImageType *src, ImageType *trg, double sigma)
{
  // typedef itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType> Filter;
  typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> Filter;
  typename Filter::Pointer flt = Filter::New();
  flt->SetInput(src);
  // flt->SetSigma(sigma);
  flt->SetVariance(sigma * sigma);
  flt->GraftOutput(trg);
  flt->Update();
}

template <class TFloat, uint VDim>
void
LDDMMData<TFloat, VDim>
::vimg_smooth(VectorImageType *src, VectorImageType *trg, double sigma)
{
  Vec sa; sa.Fill(sigma);
  vimg_smooth(src, trg, sa);
}

template <class TFloat, uint VDim>
void
LDDMMData<TFloat, VDim>
::vimg_smooth(VectorImageType *src, VectorImageType *trg, Vec sigma)
{
  typedef itk::SmoothingRecursiveGaussianImageFilter<VectorImageType, VectorImageType> Filter;
  typename Filter::Pointer fltSmooth = Filter::New();
  fltSmooth->SetInput(src);
  fltSmooth->SetSigmaArray(sigma);
  // fltSmooth->SetSigma(sigma);
  // fltSmooth->GraftOutput(trg);
  fltSmooth->Update();

  // TODO: this is a work-around for a stupid bug with this recursive filter. When the data
  // type is float, the filter does not allow me to graft an output
  LDDMMData<TFloat, VDim>::vimg_copy(fltSmooth->GetOutput(), trg);
}

template <class TFloat, uint VDim>
void
LDDMMData<TFloat, VDim>
::vimg_smooth_withborder(VectorImageType *src, VectorImageType *trg, Vec sigma, int border_size)
{

  // Define a region of interest
  RegionType region = src->GetBufferedRegion();
  region.ShrinkByRadius(border_size);

  // Perform smoothing
  vimg_smooth(src, trg, sigma);

  // Clear the border
  Vec zerovec; zerovec.Fill(0);
  typedef itk::ImageRegionIteratorWithIndex<VectorImageType> VIterator;
  for(VIterator it(trg, trg->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    if(!region.IsInside(it.GetIndex()))
      {
      it.Set(zerovec);
      }
    }
}


template <class TFloat, unsigned int VDim>
struct VectorSquareNormFunctor
{
  template <class TVector> TFloat operator() (const TVector &v)
  {
    TFloat norm_sqr = 0.0;
    for(uint i = 0; i < VDim; i++)
      norm_sqr += v[i] * v[i];
    return norm_sqr;
  }
};

template <class TFloat, uint VDim>
void
LDDMMData<TFloat, VDim>
::vimg_norm_min_max(VectorImageType *image, ImageType *normsqr,
                    TFloat &min_norm, TFloat &max_norm)
{
  // Compute the squared norm of the displacement
  typedef VectorSquareNormFunctor<TFloat, VDim> NormFunctor;
  typedef itk::UnaryFunctorImageFilter<VectorImageType, ImageType, NormFunctor> NormFilter;
  typename NormFilter::Pointer fltNorm = NormFilter::New();
  fltNorm->SetInput(image);
  fltNorm->GraftOutput(normsqr);
  fltNorm->Update();

  // Compute the maximum squared norm of the displacement
  img_min_max(normsqr, min_norm, max_norm);

  min_norm = sqrt(min_norm);
  max_norm = sqrt(max_norm);
}

template <class TFloat, uint VDim>
void
LDDMMData<TFloat, VDim>
::vimg_normalize_to_fixed_max_length(VectorImageType *trg, ImageType *normsqr,
                                     double max_displacement, bool scale_down_only)
{
  // Compute the squared norm of the displacement
  typedef VectorSquareNormFunctor<TFloat, VDim> NormFunctor;
  typedef itk::UnaryFunctorImageFilter<VectorImageType, ImageType, NormFunctor> NormFilter;
  typename NormFilter::Pointer fltNorm = NormFilter::New();
  fltNorm->SetInput(trg);
  fltNorm->GraftOutput(normsqr);
  fltNorm->Update();

  // Compute the maximum squared norm of the displacement
  TFloat nsq_min, nsq_max;
  img_min_max(normsqr, nsq_min, nsq_max);

  // Compute the scale functor
  TFloat scale = max_displacement / sqrt(nsq_max);

  // Apply the scale
  if(scale_down_only && scale >= 1.0)
    return;
  else
    vimg_scale_in_place(trg, scale);
}



namespace lddmm_data_io {

template <class TInputImage, class TOutputImage>
void
write_cast(TInputImage *image, const char *filename)
{
  typedef itk::CastImageFilter<TInputImage, TOutputImage> CastType;
  typename CastType::Pointer cast = CastType::New();
  cast->SetInput(image);

  typedef itk::ImageFileWriter<TOutputImage> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput(cast->GetOutput());
  writer->SetFileName(filename);
  writer->SetUseCompression(true);
  writer->Update();
}

template <class TInputImage, class TOutputComponent>
struct image_type_cast
{
};

template <class TPixel, unsigned int VDim, class TOutputComponent>
struct image_type_cast< itk::Image<TPixel, VDim>, TOutputComponent>
{
  typedef itk::Image<TOutputComponent, VDim> OutputImageType;
};

template <class TPixel, unsigned int VDim, class TOutputComponent>
struct image_type_cast< itk::Image<itk::CovariantVector<TPixel>, VDim>, TOutputComponent>
{
  typedef itk::Image<itk::CovariantVector<TOutputComponent>, VDim> OutputImageType;
};

template <class TPixel, unsigned int VDim, class TOutputComponent>
struct image_type_cast< itk::VectorImage<TPixel, VDim>, TOutputComponent>
{
  typedef itk::VectorImage<TOutputComponent, VDim> OutputImageType;
};

template <class TInputImage>
void write_cast_to_iocomp(TInputImage *image, const char *filename,
                          itk::ImageIOBase::IOComponentType comp)
{
  switch(comp)
    {
    case itk::ImageIOBase::UCHAR :
      write_cast<TInputImage, typename image_type_cast<TInputImage, unsigned char>::OutputImageType>(image, filename);
      break;
    case itk::ImageIOBase::CHAR :
      write_cast<TInputImage, typename image_type_cast<TInputImage, char>::OutputImageType>(image, filename);
      break;
    case itk::ImageIOBase::USHORT :
      write_cast<TInputImage, typename image_type_cast<TInputImage, unsigned short>::OutputImageType>(image, filename);
      break;
    case itk::ImageIOBase::SHORT :
      write_cast<TInputImage, typename image_type_cast<TInputImage, short>::OutputImageType>(image, filename);
      break;
    case itk::ImageIOBase::UINT :
      write_cast<TInputImage, typename image_type_cast<TInputImage, unsigned int>::OutputImageType>(image, filename);
      break;
    case itk::ImageIOBase::INT :
      write_cast<TInputImage, typename image_type_cast<TInputImage, int>::OutputImageType>(image, filename);
      break;
    case itk::ImageIOBase::ULONG :
      write_cast<TInputImage, typename image_type_cast<TInputImage, unsigned long>::OutputImageType>(image, filename);
      break;
    case itk::ImageIOBase::LONG :
      write_cast<TInputImage, typename image_type_cast<TInputImage, long>::OutputImageType>(image, filename);
      break;
    case itk::ImageIOBase::FLOAT :
      write_cast<TInputImage, typename image_type_cast<TInputImage, float>::OutputImageType>(image, filename);
      break;
    case itk::ImageIOBase::DOUBLE :
      write_cast<TInputImage, typename image_type_cast<TInputImage, double>::OutputImageType>(image, filename);
      break;
    default:
      typedef itk::ImageFileWriter<TInputImage> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput(image);
      writer->SetFileName(filename);
      writer->SetUseCompression(true);
      writer->Update();
    }
}

} // namespace


template <class TFloat, uint VDim>
typename LDDMMData<TFloat, VDim>::IOComponentType
LDDMMData<TFloat, VDim>
::img_read(const char *fn, ImagePointer &trg)
{
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(fn);
  reader->Update();
  trg = reader->GetOutput();

  return reader->GetImageIO()->GetComponentType();
}

template <class TFloat, uint VDim>
typename LDDMMData<TFloat, VDim>::ImagePointer
LDDMMData<TFloat, VDim>
::img_read(const char *fn)
{
  ImagePointer p;
  img_read(fn, p);
  return p;
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::img_write(ImageType *src, const char *fn, IOComponentType comp)
{
  lddmm_data_io::write_cast_to_iocomp(src, fn, comp);
}

template <class TFloat, uint VDim>
void
LDDMMData<TFloat, VDim>
::cimg_write(CompositeImageType *src, const char *fn, IOComponentType comp)
{
  lddmm_data_io::write_cast_to_iocomp(src, fn, comp);
}


template <class TFloat, uint VDim>
typename LDDMMData<TFloat, VDim>::IOComponentType
LDDMMData<TFloat, VDim>
::vimg_read(const char *fn, VectorImagePointer &trg)
{
  typedef itk::ImageFileReader<VectorImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(fn);
  reader->Update();
  trg = reader->GetOutput();

  return reader->GetImageIO()->GetComponentType();
}

template <class TFloat, uint VDim>
typename LDDMMData<TFloat, VDim>::VectorImagePointer
LDDMMData<TFloat, VDim>
::vimg_read(const char *fn)
{
  VectorImagePointer p;
  vimg_read(fn, p);
  return p;
}



template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::vimg_write(VectorImageType *src, const char *fn, IOComponentType comp)
{
  // Cast to vector image type
  typedef itk::VectorImage<TFloat, VDim> OutputImageType;
  typename OutputImageType::Pointer output = OutputImageType::New();
  output->CopyInformation(src);
  output->SetRegions(src->GetBufferedRegion());
  output->SetNumberOfComponentsPerPixel(VDim);

  // Override the data pointer
  output->GetPixelContainer()->SetImportPointer(
    (TFloat *) src->GetBufferPointer(), 
    VDim * src->GetPixelContainer()->Size(), false);

  // Write
  lddmm_data_io::write_cast_to_iocomp(output.GetPointer(), fn, comp);
}

template <class TFloat, uint VDim>
typename LDDMMData<TFloat, VDim>::IOComponentType
LDDMMData<TFloat, VDim>
::cimg_read(const char *fn, CompositeImagePointer &trg)
{
  typedef itk::ImageFileReader<CompositeImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(fn);
  reader->Update();
  trg = reader->GetOutput();

  return reader->GetImageIO()->GetComponentType();
}


template <class TFloat, uint VDim>
typename LDDMMData<TFloat, VDim>::CompositeImagePointer
LDDMMData<TFloat, VDim>
::cimg_read(const char *fn)
{
  CompositeImagePointer p;
  cimg_read(fn, p);
  return p;
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::vfield_read(uint nt, const char *fnpat, VelocityField &v)
{
  v.clear();
  for(uint i = 0; i < nt; i++)
    {
    char fname[1024];
    //sprintf(fname, fnpat, i);
    VectorImagePointer vt;
    vimg_read(fname, vt);
    v.push_back(vt);
    }
}

template <class TFloat, uint VDim>
void
LDDMMData<TFloat, VDim>
::vimg_copy(const VectorImageType *src, VectorImageType *trg)
{
  typedef itk::CastImageFilter<VectorImageType, VectorImageType> CastFilter;
  typename CastFilter::Pointer fltCast = CastFilter::New();
  fltCast->SetInput(src);
  fltCast->GraftOutput(trg);
  fltCast->Update();
}

template <class TFloat, uint VDim>
void
LDDMMData<TFloat, VDim>
::img_copy(const ImageType *src, ImageType *trg)
{
  typedef itk::CastImageFilter<ImageType, ImageType> CastFilter;
  typename CastFilter::Pointer fltCast = CastFilter::New();
  fltCast->SetInput(src);
  fltCast->GraftOutput(trg);
  fltCast->Update();
}

template <class TFloat, uint VDim>
void
LDDMMData<TFloat, VDim>
::mimg_copy(const MatrixImageType *src, MatrixImageType *trg)
{
  typedef itk::CastImageFilter<MatrixImageType, MatrixImageType> CastFilter;
  typename CastFilter::Pointer fltCast = CastFilter::New();
  fltCast->SetInput(src);
  fltCast->GraftOutput(trg);
  fltCast->Update();
}
template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::img_shrink(ImageType *src, ImageType *trg, int factor)
{
  typedef itk::ShrinkImageFilter<ImageType, ImageType> Filter;
  typename Filter::Pointer filter = Filter::New();
  filter->SetInput(src);
  filter->SetShrinkFactors(factor);
  filter->GraftOutput(trg);
  filter->Update();
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::img_resample_identity(ImageType *src, ImageBaseType *ref, ImageType *trg)
{
  typedef itk::ResampleImageFilter<ImageType, ImageType, TFloat> ResampleFilter;
  typedef itk::IdentityTransform<TFloat, VDim> TranType;
  typedef itk::LinearInterpolateImageFunction<ImageType, TFloat> InterpType;

  typename ResampleFilter::Pointer filter = ResampleFilter::New();
  typename TranType::Pointer tran = TranType::New();
  typename InterpType::Pointer func = InterpType::New();

  filter->SetInput(src);
  filter->SetTransform(tran);
  filter->SetInterpolator(func);
  filter->UseReferenceImageOn();
  filter->SetReferenceImage(ref);
  filter->GraftOutput(trg);
  filter->Update();
}

template <class TFloat, uint VDim>
typename LDDMMData<TFloat, VDim>::ImagePointer
LDDMMData<TFloat, VDim>
::img_downsample(ImageType *src, double factor)
{
  ImagePointer p = ImageType::New();
  img_downsample(src, p, factor);
  return p;
}

template <class TFloat, uint VDim>
void
LDDMMData<TFloat, VDim>
::img_downsample(ImageType *src, ImageType *trg, double factor)
{
  // Begin by smoothing the image
  typedef itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType> SmoothType;
  typename SmoothType::Pointer fltSmooth = SmoothType::New();
  fltSmooth->SetInput(src);
  fltSmooth->SetSigmaArray(0.5 * factor * src->GetSpacing());

  // Now resample the image to occupy the same physical space
  typedef itk::ResampleImageFilter<ImageType, ImageType, TFloat> ResampleFilter;
  typedef itk::IdentityTransform<TFloat, VDim> TranType;
  typedef itk::LinearInterpolateImageFunction<ImageType, TFloat> InterpType;

  typename ResampleFilter::Pointer filter = ResampleFilter::New();
  typename TranType::Pointer tran = TranType::New();
  typename InterpType::Pointer func = InterpType::New();

  // Compute the size of the new image
  typename ImageType::SizeType sz;
  for(uint i = 0; i < VDim; i++)
    sz[i] = (unsigned long) vcl_ceil(src->GetBufferedRegion().GetSize()[i] / factor);

  // Compute the spacing of the new image
  typename ImageType::SpacingType spc_pre = src->GetSpacing();
  typename ImageType::SpacingType spc_post = spc_pre;
  for(size_t i = 0; i < VDim; i++)
    spc_post[i] *= src->GetBufferedRegion().GetSize()[i] * 1.0 / sz[i];

  // Get the bounding box of the input image
  typename ImageType::PointType origin_pre = src->GetOrigin();

  // Recalculate the origin. The origin describes the center of voxel 0,0,0
  // so that as the voxel size changes, the origin will change as well.
  typename ImageType::SpacingType off_pre = (src->GetDirection() * spc_pre) * 0.5;
  typename ImageType::SpacingType off_post = (src->GetDirection() * spc_post) * 0.5;
  typename ImageType::PointType origin_post = origin_pre - off_pre + off_post;

  // Weird - have to allocate the output image?
  trg->SetRegions(sz);
  trg->SetOrigin(origin_post);
  trg->SetSpacing(spc_post);
  trg->SetDirection(src->GetDirection());
  trg->Allocate();

  // Set the image sizes and spacing.
  filter->SetSize(sz);
  filter->SetOutputSpacing(spc_post);
  filter->SetOutputOrigin(origin_post);
  filter->SetOutputDirection(src->GetDirection());
  filter->SetInput(fltSmooth->GetOutput());
  filter->SetTransform(tran);
  filter->SetInterpolator(func);

  filter->GraftOutput(trg);
  filter->Update();
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::vimg_resample_identity(VectorImageType *src, ImageBaseType *ref, VectorImageType *trg)
{
  typedef itk::VectorResampleImageFilter<VectorImageType, VectorImageType, TFloat> ResampleFilter;
  typedef itk::IdentityTransform<TFloat, VDim> TranType;
  typedef itk::OptVectorLinearInterpolateImageFunction<VectorImageType, TFloat> InterpType;

  typename ResampleFilter::Pointer filter = ResampleFilter::New();
  typename TranType::Pointer tran = TranType::New();
  typename InterpType::Pointer func = InterpType::New();

  filter->SetInput(src);
  filter->SetTransform(tran);
  filter->SetInterpolator(func.GetPointer());
  filter->SetSize(ref->GetBufferedRegion().GetSize());
  filter->SetOutputSpacing(ref->GetSpacing());
  filter->SetOutputOrigin(ref->GetOrigin());
  filter->SetOutputDirection(ref->GetDirection());
  filter->GraftOutput(trg);
  filter->Update();
}

template <class TFloat, uint VDim>
typename LDDMMData<TFloat, VDim>::VectorImagePointer
LDDMMData<TFloat, VDim>
::vimg_resample_identity(VectorImageType *src, ImageBaseType *ref)
{
  VectorImagePointer p = VectorImageType::New();
  vimg_resample_identity(src, ref, p);
  return p;
}

template <class TFloat, uint VDim>
void
LDDMMData<TFloat, VDim>
::img_threshold_in_place(ImageType *src, double lt, double ut, double fore, double back)
{
  typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(src);
  filter->GraftOutput(src);
  filter->SetLowerThreshold(lt);
  filter->SetUpperThreshold(ut);
  filter->SetInsideValue(fore);
  filter->SetOutsideValue(back);
  filter->Update();
}

template <class TImage>
struct VoxelToPhysicalFunctor
{

  typedef typename TImage::PixelType PixelType;
  typedef itk::ImageBase<TImage::ImageDimension> ImageBaseType;

  PixelType operator() (const PixelType &x)
  {
    typedef itk::ContinuousIndex<double, TImage::ImageDimension> CIType;
    typedef typename TImage::PointType PtType;
    CIType ci0, ci;
    PtType p0, p;
    for(uint i = 0; i < TImage::ImageDimension; i++)
      {
      ci[i] = x[i];
      ci0[i] = 0.0;
      }

    m_Image->TransformContinuousIndexToPhysicalPoint(ci, p);
    m_Image->TransformContinuousIndexToPhysicalPoint(ci0, p0);
    PixelType y;

    for(uint i = 0; i < TImage::ImageDimension; i++)
      {
      y[i] = p[i] - p0[i];
      }

    return y;
  }

  bool operator != (const VoxelToPhysicalFunctor<TImage> &other)
    { return other.m_Image != m_Image; }

  ImageBaseType *m_Image;
};

template <class TFloat, uint VDim>
void
LDDMMData<TFloat, VDim>
::warp_voxel_to_physical(VectorImageType *src, ImageBaseType *ref_space, VectorImageType *trg)
{
  // Set up functor
  typedef VoxelToPhysicalFunctor<VectorImageType> FunctorType;
  FunctorType fnk;
  fnk.m_Image = ref_space;

  // Set up filter
  typedef itk::UnaryFunctorImageFilter<VectorImageType, VectorImageType, FunctorType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(src);
  filter->GraftOutput(trg);
  filter->SetFunctor(fnk);
  filter->Update();
}



/* =============================== */

#ifdef _LDDMM_FFT_


template <class TFloat, uint VDim>
LDDMMFFTInterface<TFloat, VDim>
::LDDMMFFTInterface(ImageType *ref)
{
  // Work out the data dimensions (large enough, and multiple of four)
  m_Size = ref->GetBufferedRegion().GetSize();
  m_Alloc = m_Size;
  m_Alloc[VDim-1] = 2 * (m_Size[VDim-1] / 2 + 1);

  // Size for calling the plan routines
  int n[VDim];

  // Get the data dimensions
  m_AllocSize = 1; m_DataSize = 1;
  for(uint i = 0; i < VDim; i++)
    {
    m_AllocSize *= m_Alloc[i];
    m_DataSize *= m_Size[i];
    n[i] = m_Size[i];
    }

  // Allocate the complex data (real data is packed in the complex data)
  m_Data = (double *) fftw_malloc(m_AllocSize * sizeof(double));

  // Create plans for forward and inverse transforms
  m_Plan = fftw_plan_dft_r2c(VDim, n, m_Data, (fftw_complex *) m_Data, FFTW_MEASURE);
  m_InvPlan = fftw_plan_dft_c2r(VDim, n, (fftw_complex *) m_Data, m_Data, FFTW_MEASURE);
}

template <class TFloat, uint VDim>
void
LDDMMFFTInterface<TFloat, VDim>
::convolution_fft(
  VectorImageType *img, ImageType *kernel_ft, bool inv_kernel,
  VectorImageType *out)
{
  // Pack the data into m_Data. This requires us to skip a few bytes at
  // the end of each row of data
  uint nskip = m_Alloc[VDim-1] - m_Size[VDim-1];
  uint ncopy = m_Size[VDim-1];
  uint nout = m_Alloc[VDim-1] / 2;
  uint noutskip = kernel_ft->GetBufferedRegion().GetSize()[VDim-1] - nout;
  uint nstrides = m_AllocSize / m_Alloc[VDim-1];

  for(uint d = 0; d < VDim; d++)
    {
    const Vec *src = img->GetBufferPointer();
    double *dst = m_Data;

    // Funky loop
    for(double *end = dst + m_AllocSize; dst < end; dst+=nskip)
      for(double *rowend = dst + ncopy; dst < rowend; dst++, src++)
        *dst = (double) (*src)[d];

    // Execute the plan
    fftw_execute(m_Plan);

    // Multiply or divide the complex values in m_Data by the kernel array
    fftw_complex *c = (fftw_complex *) m_Data;

    // Scaling factor for producing final result
    double scale = 1.0 / m_DataSize;

    /*
    // Precision weirdness (results differ from MATLAB fft in 6th, 7th decimal digit)
    uint tp = (m_Alloc[VDim-1] / 2) * 6 + 8;
    printf("Before scaling, value at %d is %12.12lf, %12.12lf\n",
      tp, c[tp][0], c[tp][1]);
    printf("Kernel value at %d is %12.12lf\n", 128*6+8, kp[128*6+8]);
    */

    // Another such loop
    TFloat *kp = kernel_ft->GetBufferPointer();
    if(inv_kernel)
      {
      for(uint i = 0; i < nstrides; i++)
        {
        for(uint j = 0; j < nout; j++)
          {
          (*c)[0] /= (*kp);
          (*c)[1] /= (*kp);
          c++; kp++;
          }
        kp += noutskip;
        }
      }
    else
      {
      for(uint i = 0; i < nstrides; i++)
        {
        for(uint j = 0; j < nout; j++)
          {
          (*c)[0] *= (*kp);
          (*c)[1] *= (*kp);
          c++; kp++;
          }
        kp += noutskip;
        }
      }

    /*
    fftw_complex *ctest = (fftw_complex *) m_Data;
    printf("After scaling, value at %d is %12.12lf, %12.12lf\n",
      tp, ctest[tp][0], ctest[tp][1]);
    */

    // Inverse transform
    fftw_execute(m_InvPlan);

    // Copy the results to the output image
    const double *osrc = m_Data;
    Vec *odst = out->GetBufferPointer();
    for(uint i = 0; i < nstrides; i++, osrc+=nskip)
      for(uint j = 0; j < ncopy; j++, odst++, osrc++)
        (*odst)[d] = (TFloat) ((*osrc) * scale);

    /*
    odst = out->GetBufferPointer();
    printf("Result %12.12lf\n",  odst[128*6+8][0]);
    */
    }

}

template <class TFloat, uint VDim>
LDDMMFFTInterface<TFloat, VDim>
::~LDDMMFFTInterface()
{
  fftw_destroy_plan(m_Plan);
  fftw_destroy_plan(m_InvPlan);
  fftw_free(m_Data);
}

#endif // _LDDMM_FFT_


template <class TFloat, uint VDim>
LDDMMImageMatchingObjective<TFloat, VDim>
::LDDMMImageMatchingObjective(LDDMM &p)
  : fft(p.fix)
{
  // Allocate intermediate datasets
  LDDMM::alloc_img(Jt0, p.fix);
  LDDMM::alloc_img(Jt1, p.fix);
  LDDMM::alloc_img(DetPhit1, p.fix);
  LDDMM::alloc_vimg(GradJt0, p.fix);
}

template <class TFloat, uint VDim>
TFloat
LDDMMImageMatchingObjective<TFloat, VDim>
::compute_objective_and_gradient(LDDMM &p)
{
  // Compute the regularization energy of v. We can use a[0] for temporary storage
  // e_field = lddmm_vector_field_dot_product(vx, vy, vx, vy, p);
  TFloat e_field = 0.0;
  for(uint m = 0; m < p.nt; m++)
    {
    // a[0] = Lv[m] .* v[m]
    fft.convolution_fft(p.v[m], p.f_kernel_sq, false, p.a[0]);

    // We're sticking the inner product in Jt0
    LDDMM::vimg_euclidean_inner_product(Jt0, p.a[0], p.v[m]); 
    e_field += LDDMM::img_voxel_sum(Jt0) / p.nt;
    }

  // Compute the 'a' array (for semilagrangean scheme)
  p.compute_semi_lagrangean_a();

  // Go through and compute phi_t1 for each time point
  p.integrate_phi_t1();

  // Compute the update for v at each time step
  for(uint m = 0; m < p.nt; m++)
    {
    // Currently, f[m] holds phi_t1[m]. Use it for whatever we need
    // and then replace with phi_t0[m]

    // TODO: for ft00 and ft11, don't waste time on interpolation

    // Jt1 = lddmm_warp_scalar_field(p.I1, ft1x(:,:,it), ft1y(:,:,it), p);
    LDDMM::interp_img(p.mov, p.f[m], Jt1); 

    // detjac_phi_t1 = lddmm_jacobian_determinant(ft1x(:,:,it), ft1y(:,:,it), p);
    LDDMM::field_jacobian_det(p.f[m], DetPhit1);

    // Place phi_t0 into the f array
    if(m==0)
      {
      p.f[m]->FillBuffer(typename LDDMM::Vec(0.0));
      }
    else
      {
      LDDMM::interp_vimg(p.f[m-1], p.a[m], -1.0, p.f[m]);
      LDDMM::vimg_subtract_in_place(p.f[m], p.a[m]);
      }

    // Jt0 = lddmm_warp_scalar_field(p.I0, ft0x(:,:,it), ft0y(:,:,it), p); 
    LDDMM::interp_img(p.fix, p.f[m], Jt0); 

    // [grad_Jt0_x grad_Jt0_y] = gradient(Jt0);
    LDDMM::image_gradient(Jt0, GradJt0);

    // pde_rhs_x = detjac_phi_t1 .* (Jt0 - Jt1) .* grad_Jt0_x; 
    // pde_rhs_y = detjac_phi_t1 .* (Jt0 - Jt1) .* grad_Jt0_y; 

    // Here we do some small tricks. We want to retain Jt0 because it's the warped
    // template image, and we want to retain the difference Jt0-Jt1 = (J0-I1) for
    // calculating the objective at the end. 
    LDDMM::img_subtract_in_place(Jt1, Jt0);           // 'Jt1' stores Jt1 - Jt0 
    LDDMM::img_multiply_in_place(DetPhit1, Jt1);      // 'DetPhit1' stores (det Phi_t1)(Jt1-Jt0)
    LDDMM::vimg_multiply_in_place(GradJt0, DetPhit1); // 'GradJt0' stores  GradJt0 * (det Phi_t1)(Jt1-Jt0)

    // Solve PDE via FFT convolution
    // pde_soln_x = ifft2(fft2(pde_rhs_x) ./ p.f_kernel_sq,'symmetric');
    // pde_soln_y = ifft2(fft2(pde_rhs_y) ./ p.f_kernel_sq,'symmetric');
    fft.convolution_fft(GradJt0, p.f_kernel_sq, true, GradJt0); // 'GradJt0' stores K[ GradJt0 * (det Phi_t1)(Jt1-Jt0) ]

    // dedvx(:,:,it) = dedvx(:,:,it) - 2 * pde_soln_x / p.sigma^2;
    // dedvy(:,:,it) = dedvy(:,:,it) - 2 * pde_soln_y / p.sigma^2;        

    // Store the update in a[m]
    LDDMM::vimg_scale_in_place(GradJt0, 1.0 / p.sigma_sq); // 'GradJt0' stores 1 / sigma^2 K[ GradJt0 * (det Phi_t1)(Jt1-Jt0) ]
    LDDMM::vimg_add_in_place(GradJt0, p.v[m]); // 'GradJt0' stores v + 1 / sigma^2 K[ GradJt0 * (det Phi_t1)(Jt1-Jt0) ]
    LDDMM::vimg_scale(GradJt0, 2.0, p.a[m]); // p.a[m] holds 2 v + 2 / sigma^2 K[ GradJt0 * (det Phi_t1)(Jt1-Jt0) ]
    }

  // Ok, Jt1 currently contains (Jt1-Jt0), we just need to square it.
  TFloat e_image = LDDMM::img_euclidean_norm_sq(Jt1) / p.sigma_sq;

  // Return the energy
  printf("  Energy components: %lf, %lf\n", e_field, e_image);
  return e_field + e_image;
}


template class LDDMMData<float, 2>;
template class LDDMMData<float, 3>;
template class LDDMMData<float, 4>;

template class LDDMMData<double, 2>;
template class LDDMMData<double, 3>;
template class LDDMMData<double, 4>;

#ifdef _LDDMM_FFT_
template class LDDMMFFTInterface<double, 2>;
template class LDDMMFFTInterface<double, 3>;
template class LDDMMFFTInterface<double, 4>;
#endif

template class LDDMMImageMatchingObjective<myreal, 2>;
template class LDDMMImageMatchingObjective<myreal, 3>;

