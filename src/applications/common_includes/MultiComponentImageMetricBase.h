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
#ifndef MULTICOMPONENTIMAGEMETRICBASE_H
#define MULTICOMPONENTIMAGEMETRICBASE_H

#include "itkImageBase.h"
#include "itkImageToImageFilter.h"
#include "itkPoint.h"
#include "itkFixedArray.h"
#include "itkVectorImage.h"
#include "itkMatrixOffsetTransformBase.h"
#include "lddmm_common.h"

/**
 * Default traits for parameterizing the metric filters below
 */
template <class TReal, unsigned int VDim> struct DefaultMultiComponentImageMetricTraits
{
  typedef itk::VectorImage<TReal, VDim> InputImageType;
  typedef itk::Image<TReal, VDim> ScalarImageType;
  typedef itk::Image<itk::CovariantVector<TReal, VDim>, VDim> VectorImageType;

  typedef ScalarImageType MaskImageType;
  typedef VectorImageType DeformationFieldType;
  typedef VectorImageType GradientImageType;
  typedef ScalarImageType MetricImageType;
  typedef itk::MatrixOffsetTransformBase<TReal, VDim, VDim> TransformType;

  typedef TReal RealType;
};


/**
 * \class MultiComponentImageMetricBase
 *
 * Base class for metrics that compute similarity between two multi-component
 * images based on a deformation field. This filter is extended to support
 * normalized cross-correlation and least squares metrics.
 *
 * The metric takes the following inputs:
 *    Fixed multicomponent image (type TMetricTraits::InputImageType)
 *    Moving multicomponent image (type TMetricTraits::InputImageType)
 *    A mask for the fixed image (type TMetricTraits::MaskImageType) [optional]
 *    A deformation field (type TMetricTraits::DeformationFieldType)
 *
 * It produces the following outputs
 *    Metric image (type TMetricTraits::MetricImageType)
 *    Metric image gradient (type TMetricTraits::GradientImageType)
 *    Moving image domain mask (type TMetricTraits::MetricImageType)
 *    Moving image domain mask gradient (type TMetricTraits::GradientImageType)
 */
template <class TMetricTraits>
class ITK_EXPORT MultiComponentImageMetricBase :
    public itk::ImageToImageFilter<typename TMetricTraits::InputImageType,
                                   typename TMetricTraits::MetricImageType>
{
public:
  /** Type definitions from the traits class */
  typedef typename TMetricTraits::InputImageType        InputImageType;
  typedef typename TMetricTraits::MaskImageType         MaskImageType;
  typedef typename TMetricTraits::DeformationFieldType  DeformationFieldType;
  typedef typename TMetricTraits::MetricImageType       MetricImageType;
  typedef typename TMetricTraits::GradientImageType     GradientImageType;
  typedef typename TMetricTraits::TransformType         TransformType;
  typedef typename TMetricTraits::RealType              RealType;

  /** Standard class typedefs. */
  typedef MultiComponentImageMetricBase                            Self;
  typedef itk::ImageToImageFilter<InputImageType,MetricImageType>  Superclass;
  typedef itk::SmartPointer<Self>                                  Pointer;
  typedef itk::SmartPointer<const Self>                            ConstPointer;

  /** Run-time type information (and related methods) */
  itkTypeMacro( MultiComponentImageMetricBase, ImageToImageFilter )

  /** Determine the image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      MetricImageType::ImageDimension );

  /** Typedef to describe the output image region type. */
  typedef typename InputImageType::RegionType OutputImageRegionType;

  /** Inherit some types from the superclass. */
  typedef typename InputImageType::PixelType          InputPixelType;
  typedef typename InputImageType::InternalPixelType  InputComponentType;
  typedef typename MetricImageType::PixelType         MetricPixelType;
  typedef typename MetricImageType::IndexType         IndexType;
  typedef typename MetricImageType::IndexValueType    IndexValueType;
  typedef typename MetricImageType::SizeType          SizeType;
  typedef typename MetricImageType::SpacingType       SpacingType;
  typedef typename MetricImageType::DirectionType     DirectionType;
  typedef typename GradientImageType::PixelType       GradientPixelType;
  typedef itk::ImageBase<ImageDimension>              ImageBaseType;

  typedef typename Superclass::DataObjectIdentifierType DataObjectIdentifierType;

  /** Information from the deformation field class */
  typedef typename DeformationFieldType::Pointer      DeformationFieldPointer;
  typedef typename DeformationFieldType::PixelType    DeformationVectorType;

  /** Weight vector */
  typedef vnl_vector<float>                           WeightVectorType;

  /** Set the fixed image(s) */
  itkNamedInputMacro(FixedImage, InputImageType, "Primary")

  /** Set the moving image(s) and their gradients */
  itkNamedInputMacro(MovingImage, InputImageType, "moving")

  /** Set the optional mask input */
  itkNamedInputMacro(FixedMaskImage, MaskImageType, "fixed_mask")

  /** Set the optional jitter input - for affine images*/
  itkNamedInputMacro(JitterImage, DeformationFieldType, "jitter")

  /**
   * Set the deformation field. If the deformation field is set, the affine
   * transform gradients will not be computed.
   */
  void SetDeformationField(DeformationFieldType *phi)
  {
    this->itk::ProcessObject::SetInput("phi", phi);
    this->UpdateOutputs();
  }

  itkNamedInputGetMacro(DeformationField, DeformationFieldType, "phi")

  /**
   * Set the affine transform. If the affine transform is set, affine transform
   * gradients will be computed, but not the deformation field gradients
   */
  void SetAffineTransform(TransformType *transform)
  {
    this->m_AffineTransform = transform;
    this->m_ComputeAffine = true;
    this->UpdateOutputs();
  }

  itkGetMacro(AffineTransform, TransformType *)

  /** Set the weight vector - for different components in the input image */
  itkSetMacro(Weights, WeightVectorType)
  itkGetConstMacro(Weights, WeightVectorType)

  /** Whether or not the gradient is required */
  itkGetMacro(ComputeGradient, bool)

  /** Whether the transformation is affine */
  itkGetMacro(ComputeAffine, bool)

  /** Specify whether the filter should compute gradients (whether affine or deformable) */
  void SetComputeGradient(bool flag)
  {
    this->m_ComputeGradient = flag;
    this->UpdateOutputs();
  }

  itkBooleanMacro(ComputeMovingDomainMask)

  /** Specify whether the metric should be normalized by the moving image domain */
  void SetComputeMovingDomainMask(bool flag)
  {
    this->m_ComputeMovingDomainMask = flag;
    this->UpdateOutputs();
  }

  /** Get the metric image output - this is the main output */
  itkNamedOutputMacro(MetricOutput, MetricImageType, "Primary")

  /** Get the metric dense gradient output. The gradient may be arbitrarily scaled. */
  itkNamedOutputMacro(DeformationGradientOutput, GradientImageType, "phi_gradient")

  /** Get the gradient of the affine transform */
  itkGetMacro(AffineTransformGradient, TransformType *)

  /**
   * Get the gradient scaling factor. To get the actual gradient of the metric, multiply the
   * gradient output of this filter by the scaling factor. Explanation: for efficiency, the
   * metrics return an arbitrarily scaled vector, such that adding the gradient to the
   * deformation field would INCREASE SIMILARITY. For metrics that are meant to be minimized,
   * this is the opposite of the gradient direction. For metrics that are meant to be maximized,
   * it is the gradient direction.
   */
  virtual double GetGradientScalingFactor() const = 0;

  /** Summary results after running the filter */
  itkGetConstMacro(MetricValue, double)

  /** Get the metric values per component (each component weighted) */
  vnl_vector<double> GetAllMetricValues() const;

protected:
  MultiComponentImageMetricBase();
  ~MultiComponentImageMetricBase() {}

  virtual void VerifyInputInformation() ITK_OVERRIDE {}

  /** It is difficult to compute in advance the input image region
   * required to compute the requested output region. Thus the safest
   * thing to do is to request for the whole input image.
   *
   * For the deformation field, the input requested region
   * set to be the same as that of the output requested region. */
  virtual void GenerateInputRequestedRegion() ITK_OVERRIDE;

  virtual typename itk::DataObject::Pointer MakeOutput(const DataObjectIdentifierType &) ITK_OVERRIDE;

  void UpdateOutputs();
  void ToggleOutput(bool flag, const DataObjectIdentifierType &key);

  virtual void BeforeThreadedGenerateData() ITK_OVERRIDE;
  virtual void AfterThreadedGenerateData() ITK_OVERRIDE;



  // Weight vector
  WeightVectorType                m_Weights;

  bool m_ComputeMovingDomainMask;
  bool m_ComputeGradient;
  bool m_ComputeAffine;

  // Data accumulated for each thread
  struct ThreadData {
    double metric, mask;

    // Component-wise metric values - for reporting
    vnl_vector<double> comp_metric;

    vnl_vector<double> gradient, grad_mask;
    ThreadData() : metric(0.0), mask(0.0),
      gradient(ImageDimension * (ImageDimension+1), 0.0),
      grad_mask(ImageDimension * (ImageDimension+1), 0.0) {}
  };

  // Per-thread data
  std::vector<ThreadData> m_ThreadData;

  // Total accumulated data
  ThreadData m_AccumulatedData;

  // Accumulated metric value
  double m_MetricValue;

  // Affine transform
  typename TransformType::Pointer m_AffineTransform, m_AffineTransformGradient;

private:
  MultiComponentImageMetricBase(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};


/**
 * Flatten an affine transform to a flat array
 */
template<class TFloat, class TFloatArr, unsigned int VDim>
static void flatten_affine_transform(
    const itk::MatrixOffsetTransformBase<TFloat, VDim, VDim> *transform,
    TFloatArr *flat_array)
{
  int pos = 0;
  for(int i = 0; i < VDim; i++)
    {
    flat_array[pos++] = transform->GetOffset()[i];
    for(int j = 0; j < VDim; j++)
      flat_array[pos++] = transform->GetMatrix()(i,j);
    }
}

template<class TFloat, class TFloatArr, unsigned int VDim>
static void flatten_affine_transform(
    const vnl_matrix_fixed<TFloat, VDim, VDim> &matrix,
    const vnl_vector_fixed<TFloat, VDim> &offset,
    TFloatArr *flat_array)
{
  int pos = 0;
  for(int i = 0; i < VDim; i++)
    {
    flat_array[pos++] = offset[i];
    for(int j = 0; j < VDim; j++)
      flat_array[pos++] = matrix(i,j);
    }
}

/**
 * Unflatten a flat array to an affine transform
 */
template<class TFloat, class TFloatArr, unsigned int VDim>
static void unflatten_affine_transform(
   const TFloatArr *flat_array,
   itk::MatrixOffsetTransformBase<TFloat, VDim, VDim> *transform,
   double scaling = 1.0)
{
  typename itk::MatrixOffsetTransformBase<TFloat, VDim, VDim>::MatrixType matrix;
  typename itk::MatrixOffsetTransformBase<TFloat, VDim, VDim>::OffsetType offset;

  int pos = 0;
  for(int i = 0; i < VDim; i++)
    {
    offset[i] = flat_array[pos++] * scaling;
    for(int j = 0; j < VDim; j++)
      matrix(i, j) = flat_array[pos++] * scaling;
    }

  transform->SetMatrix(matrix);
  transform->SetOffset(offset);
}

template<class TFloat, class TFloatArr, unsigned int VDim>
static void unflatten_affine_transform(
    const TFloatArr *flat_array,
    vnl_matrix_fixed<TFloat, VDim, VDim> &matrix,
    vnl_vector_fixed<TFloat, VDim> &offset,
    double scaling = 1.0)
{
  int pos = 0;
  for(int i = 0; i < VDim; i++)
    {
    offset[i] = flat_array[pos++] * scaling;
    for(int j = 0; j < VDim; j++)
      matrix(i, j) = flat_array[pos++] * scaling;
    }
}



template <class TFloat, unsigned int VDim>
static void set_affine_transform(
    const vnl_matrix_fixed<TFloat, VDim, VDim> &matrix,
    const vnl_vector_fixed<TFloat, VDim> &offset,
    itk::MatrixOffsetTransformBase<TFloat, VDim, VDim> *transform)
{
  typename itk::MatrixOffsetTransformBase<TFloat, VDim, VDim>::MatrixType tmatrix;
  typename itk::MatrixOffsetTransformBase<TFloat, VDim, VDim>::OffsetType toffset;

  for(int i = 0; i < VDim; i++)
    {
    toffset[i] = offset[i];
    for(int j = 0; j < VDim; j++)
      tmatrix(i, j) = matrix(i,j);
    }

  transform->SetMatrix(tmatrix);
  transform->SetOffset(toffset);
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "MultiComponentImageMetricBase.txx"
#endif


#endif // MULTICOMPONENTIMAGEMETRICBASE_H
