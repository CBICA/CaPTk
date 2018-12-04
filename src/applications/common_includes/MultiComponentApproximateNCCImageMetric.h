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
#ifndef MULTICOMPONENTAPPROXIMATENCCIMAGEMETRIC_H
#define MULTICOMPONENTAPPROXIMATENCCIMAGEMETRIC_H

#include "MultiComponentImageMetricBase.h"
#include "itkBarrier.h"

/**
 * Normalized cross-correlation metric similar to the one used in ANTS. The gradient
 * is an approximation, but it seems to work very well.
 */
template <class TMetricTraits>
class ITK_EXPORT MultiComponentApproximateNCCImageMetric :
    public MultiComponentImageMetricBase<TMetricTraits>
{
public:
  /** Standard class typedefs. */
  typedef MultiComponentApproximateNCCImageMetric<TMetricTraits> Self;
  typedef MultiComponentImageMetricBase<TMetricTraits>      Superclass;
  typedef itk::SmartPointer<Self>                           Pointer;
  typedef itk::SmartPointer<const Self>                     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self)

  /** Run-time type information (and related methods) */
  itkTypeMacro( MultiComponentApproximateNCCImageMetric, MultiComponentImageMetricBase )

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
  typedef typename Superclass::RealType                      RealType;

  /** Information from the deformation field class */
  typedef typename Superclass::DeformationFieldType          DeformationFieldType;

  /** Determine the image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int, InputImageType::ImageDimension );

  /** Set the radius of the cross-correlation */
  itkSetMacro(Radius, SizeType)

  /** Get the radius of the cross-correlation */
  itkGetMacro(Radius, SizeType)

  /**
   * Set the working memory image for this filter. This function should be used to prevent
   * repeated allocation of memory when the metric is created/destructed in a loop. The
   * user can just pass in a pointer to a blank image, the filter will take care of allocating
   * the image as necessary
   */
  itkSetObjectMacro(WorkingImage, InputImageType)

  /**
   * Get the gradient scaling factor. To get the actual gradient of the metric, multiply the
   * gradient output of this filter by the scaling factor. Explanation: for efficiency, the
   * metrics return an arbitrarily scaled vector, such that adding the gradient to the
   * deformation field would INCREASE SIMILARITY. For metrics that are meant to be minimized,
   * this is the opposite of the gradient direction. For metrics that are meant to be maximized,
   * it is the gradient direction.
   */
  virtual double GetGradientScalingFactor() const { return 1.0; }


protected:
  MultiComponentApproximateNCCImageMetric()
    { m_Radius.Fill(1); }

  ~MultiComponentApproximateNCCImageMetric() {}

  // TODO: set up for proper streaming
  // virtual void GenerateInputRequestedRegion();

  virtual void BeforeThreadedGenerateData();
  virtual void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread,
                                    itk::ThreadIdType threadId);


private:
  MultiComponentApproximateNCCImageMetric(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  // A pointer to the working image. The user should supply this image in order to prevent
  // unnecessary memory allocation
  typename InputImageType::Pointer m_WorkingImage;

  // Radius of the cross-correlation
  SizeType m_Radius;

  // Whether to compute the approximate gradient of the metric, as in Avants et al., 2008.
  // When not set, the exact gradient is computed
  bool m_ApproximateGradient;
};






/**
 * \class MultiImageApproximateNCCPrecomputeFilter
 * \brief Helps compute the NCC metric
 *
 * This filter takes a pair of images plus a warp and computes the components that
 * are used to calculate the cross-correlation metric between them and
 * the gradient. These components are in the form I, I*J, and
 * so on. These components must then be mean-filtered and combined to get the
 * metric and the gradient.
 *
 * The output of this filter must be a vector image. The input may be a vector image.
 *
 */
template <class TMetricTraits, class TOutputImage>
class ITK_EXPORT MultiImageApproximateNCCPrecomputeFilter :
    public itk::ImageToImageFilter<typename TMetricTraits::InputImageType, TOutputImage>
{
public:

  /** The parent/owner class */
  typedef MultiComponentApproximateNCCImageMetric<TMetricTraits> ParentType;

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
  typedef MultiImageApproximateNCCPrecomputeFilter              Self;
  typedef itk::ImageToImageFilter<InputImageType, TOutputImage> Superclass;
  typedef itk::SmartPointer<Self>                               Pointer;
  typedef itk::SmartPointer<const Self>                         ConstPointer;

  /** Typedef to describe the output image region type. */
  typedef typename Superclass::OutputImageRegionType         OutputImageRegionType;

  typedef typename Superclass::DataObjectIdentifierType DataObjectIdentifierType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self)

  /** Run-time type information (and related methods) */
  itkTypeMacro( MultiImageApproximateNCCPrecomputeFilter, ImageToImageFilter )

  /** Determine the image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int, TOutputImage::ImageDimension );

  /** Set the parent class */
  void SetParent(ParentType *parent)
  {
    m_Parent = parent;
  }

  enum Stage {FIRST, SECOND};

  /** Set the stage of the computation */
  itkSetMacro(Stage, Stage)

  /** Get the number of components in the output */
  int GetNumberOfOutputComponents();


protected:
  MultiImageApproximateNCCPrecomputeFilter();
  ~MultiImageApproximateNCCPrecomputeFilter() {}

  /** SimpleWarpImageFilter is implemented as a multi-threaded filter.
   * As such, it needs to provide and implementation for
   * ThreadedGenerateData(). */
  virtual void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                                    itk::ThreadIdType threadId );

  /** Set up the output information */
  virtual void GenerateOutputInformation();

  /** Override input checks to allow fixed and moving to be in different space */
  virtual void VerifyInputInformation() {}

private:
  MultiImageApproximateNCCPrecomputeFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  ParentType *m_Parent;

  Stage m_Stage;
};






#ifndef ITK_MANUAL_INSTANTIATION
#include "MultiComponentApproximateNCCImageMetric.txx"
#endif


#endif // MULTICOMPONENTAPPROXIMATENCCIMAGEMETRIC_H
