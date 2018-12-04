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
#ifndef MULTICOMPONENTNCCIMAGEMETRIC_H
#define MULTICOMPONENTNCCIMAGEMETRIC_H

#include "MultiComponentImageMetricBase.h"
#include "itkBarrier.h"

/**
 * Normalized cross-correlation metric. This filter sets up a mini-pipeline with
 * a pre-compute filter that interpolates the moving image, N one-dimensional
 * mean filters, and a post-compute filter that generates the metric and the
 * gradient.
 */
template <class TMetricTraits>
class ITK_EXPORT MultiComponentNCCImageMetric :
    public MultiComponentImageMetricBase<TMetricTraits>
{
public:
  /** Standard class typedefs. */
  typedef MultiComponentNCCImageMetric<TMetricTraits>       Self;
  typedef MultiComponentImageMetricBase<TMetricTraits>      Superclass;
  typedef itk::SmartPointer<Self>                           Pointer;
  typedef itk::SmartPointer<const Self>                     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self)

  /** Run-time type information (and related methods) */
  itkTypeMacro( MultiComponentNCCImageMetric, MultiComponentImageMetricBase )

  /** Typedef to describe the output image region type. */
  typedef typename Superclass::OutputImageRegionType         OutputImageRegionType;

  /** Inherit some types from the superclass. */
  typedef typename Superclass::InputImageType                InputImageType;
  typedef typename Superclass::InputPixelType                InputPixelType;
  typedef typename Superclass::InputComponentType            InputComponentType;
  typedef typename Superclass::MetricImageType               MetricImageType;
  typedef typename Superclass::MetricPixelType               MetricPixelType;
  typedef typename Superclass::GradientImageType             GradientImageType;
  typedef typename Superclass::GradientPixelType             GradientPixelType;
  typedef typename Superclass::MaskImageType                 MaskImageType;


  typedef typename Superclass::IndexType                     IndexType;
  typedef typename Superclass::IndexValueType                IndexValueType;
  typedef typename Superclass::SizeType                      SizeType;
  typedef typename Superclass::SpacingType                   SpacingType;
  typedef typename Superclass::DirectionType                 DirectionType;
  typedef typename Superclass::ImageBaseType                 ImageBaseType;

  /** Information from the deformation field class */
  typedef typename Superclass::DeformationFieldType          DeformationFieldType;

  /** Determine the image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int, InputImageType::ImageDimension );

  /** Set the radius of the cross-correlation */
  itkSetMacro(Radius, SizeType)

  /** Get the radius of the cross-correlation */
  itkGetMacro(Radius, SizeType)

  /**
   * Set this flag to true to compute the approximate gradient of the metric, in the
   * direction of the (fixed - why?) image gradient.
   */
  itkSetMacro(ApproximateGradient, bool)

  /** Get the approximate gradient flag */
  itkGetMacro(ApproximateGradient, bool)

  /**
   * Set the working memory image for this filter. This function should be used to prevent
   * repeated allocation of memory when the metric is created/destructed in a loop. The
   * user can just pass in a pointer to a blank image, the filter will take care of allocating
   * the image as necessary
   */
  itkSetObjectMacro(WorkingImage, InputImageType)

  /**
   * Whether we should reuse the fixed components in the working image. This is true
   * if the filter is being run repeatedly on the same image
   */
  itkSetMacro(ReuseWorkingImageFixedComponents, bool)

  /**
   * Get the gradient scaling factor. To get the actual gradient of the metric, multiply the
   * gradient output of this filter by the scaling factor. Explanation: for efficiency, the
   * metrics return an arbitrarily scaled vector, such that adding the gradient to the
   * deformation field would INCREASE SIMILARITY. For metrics that are meant to be minimized,
   * this is the opposite of the gradient direction. For metrics that are meant to be maximized,
   * it is the gradient direction.
   */
  virtual double GetGradientScalingFactor() const ITK_OVERRIDE { return 1.0; }


protected:
  MultiComponentNCCImageMetric()
    : m_ApproximateGradient(false), m_ReuseWorkingImageFixedComponents(false)
    { m_Radius.Fill(1); }

  ~MultiComponentNCCImageMetric() {}

  // TODO: set up for proper streaming
  // virtual void GenerateInputRequestedRegion();

  virtual void BeforeThreadedGenerateData() ITK_OVERRIDE;
  virtual void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread,
                                    itk::ThreadIdType threadId) ITK_OVERRIDE;

private:
  MultiComponentNCCImageMetric(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  // A pointer to the working image. The user should supply this image in order to prevent
  // unnecessary memory allocation
  typename InputImageType::Pointer m_WorkingImage;

  // Whether we should reuse the fixed components in the working image. This is true
  // if the filter is being run repeatedly on the same image
  bool m_ReuseWorkingImageFixedComponents;

  // Radius of the cross-correlation
  SizeType m_Radius;

  // Whether to compute the approximate gradient of the metric, as in Avants et al., 2008.
  // When not set, the exact gradient is computed
  bool m_ApproximateGradient;
};






/**
 * \class MultiImageNCCPrecomputeFilter
 * \brief Helps compute the NCC metric
 *
 * This filter takes a pair of images plus a warp and computes the components that
 * are used to calculate the cross-correlation metric between them and
 * the gradient. These components are in the form I, I*J, I * gradJ, and
 * so on. These components must then be mean-filtered and combined to get the
 * metric and the gradient.
 *
 * The output of this filter must be a vector image. The input may be a vector image.
 *
 */
template <class TMetricTraits, class TOutputImage>
class ITK_EXPORT MultiImageNCCPrecomputeFilter :
    public itk::ImageToImageFilter<typename TMetricTraits::InputImageType, TOutputImage>
{
public:

  /** The parent/owner class */
  typedef MultiComponentNCCImageMetric<TMetricTraits> ParentType;

  /** Inherit some types from the superclass. */
  typedef typename ParentType::InputImageType         InputImageType;
  typedef typename ParentType::InputPixelType         InputPixelType;
  typedef typename ParentType::InputComponentType     InputComponentType;
  typedef typename ParentType::DeformationVectorType  DeformationVectorType;
  typedef typename ParentType::MetricPixelType        MetricPixelType;
  typedef typename ParentType::GradientPixelType      GradientPixelType;
  typedef typename ParentType::RealType               RealType;
  typedef TOutputImage                                OutputImageType;
  typedef typename OutputImageType::InternalPixelType OutputComponentType;

  /** Standard class typedefs. */
  typedef MultiImageNCCPrecomputeFilter                         Self;
  typedef itk::ImageToImageFilter<InputImageType, TOutputImage> Superclass;
  typedef itk::SmartPointer<Self>                               Pointer;
  typedef itk::SmartPointer<const Self>                         ConstPointer;

  /** Typedef to describe the output image region type. */
  typedef typename Superclass::OutputImageRegionType         OutputImageRegionType;


  typedef typename Superclass::DataObjectIdentifierType DataObjectIdentifierType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self)

  /** Run-time type information (and related methods) */
  itkTypeMacro( MultiImageNCCPrecomputeFilter, ImageToImageFilter )

  /** Determine the image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int, TOutputImage::ImageDimension );

  /** Set the parent class */
  void SetParent(ParentType *parent)
  {
    m_Parent = parent;
  }

  /** Get the number of components in the output */
  int GetNumberOfOutputComponents();

  /** Get the number of components in the output */
  int GetNumberOfFixedOnlyOutputComponents();

  /**
   * Set whether the filter is being run for the first time. When run for the first
   * time, fixed components and products will be stored in the output. Otherwise,
   * these are skipped over in memory
   */
  itkSetMacro(FlagGenerateFixedComponents, bool)


protected:
  MultiImageNCCPrecomputeFilter();
  ~MultiImageNCCPrecomputeFilter() {}

  /** SimpleWarpImageFilter is implemented as a multi-threaded filter.
   * As such, it needs to provide and implementation for
   * ThreadedGenerateData(). */
  virtual void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                            itk::ThreadIdType threadId ) ITK_OVERRIDE;

  /** Set up the output information */
  virtual void GenerateOutputInformation() ITK_OVERRIDE;

  /** Override input checks to allow fixed and moving to be in different space */
  virtual void VerifyInputInformation() ITK_OVERRIDE {}

private:
  MultiImageNCCPrecomputeFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  ParentType *m_Parent;

  bool m_FlagGenerateFixedComponents;
};






#ifndef ITK_MANUAL_INSTANTIATION
#include "MultiComponentNCCImageMetric.txx"
#endif


#endif // MULTICOMPONENTNCCIMAGEMETRIC_H
