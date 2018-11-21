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
#ifndef FASTWARPCOMPOSITEIMAGEFILTER_H
#define FASTWARPCOMPOSITEIMAGEFILTER_H

#include "lddmm_common.h"
#include "itkImageToImageFilter.h"

/**
 * This is a warp filter that is fast, supports multi-component (vector) images
 * and can be performed using either physical space or voxel-space calculations
 *
 * This filter supports linear interpolation currently
 */
template <class TInputImage, class TOutputImage, class TDeformationField>
class FastWarpCompositeImageFilter
        : public itk::ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  typedef FastWarpCompositeImageFilter<TInputImage,TOutputImage,TDeformationField> Self;
  typedef itk::ImageToImageFilter<TInputImage, TOutputImage>       Superclass;
  typedef itk::SmartPointer<Self>                                  Pointer;
  typedef itk::SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self)

  /** Run-time type information (and related methods) */
  itkTypeMacro( FastWarpCompositeImageFilter, ImageToImageFilter )

  /** Determine the image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension );

  // Lots of typedefs
  typedef TInputImage                                 InputImageType;
  typedef TOutputImage                                OutputImageType;
  typedef TDeformationField                           DeformationFieldType;
  typedef typename Superclass::OutputImageRegionType  OutputImageRegionType;
  typedef typename InputImageType::PixelType          InputPixelType;
  typedef typename InputImageType::InternalPixelType  InputComponentType;
  typedef typename InputImageType::IndexType          IndexType;
  typedef typename InputImageType::IndexValueType     IndexValueType;
  typedef typename InputImageType::SizeType           SizeType;
  typedef typename InputImageType::SpacingType        SpacingType;
  typedef typename InputImageType::DirectionType      DirectionType;
  typedef typename DeformationFieldType::PixelType    DeformationVectorType;
  typedef typename DeformationVectorType::ValueType   DeforamtionScalarType;
  typedef typename OutputImageType::PixelType         OutputPixelType;
  typedef typename OutputImageType::InternalPixelType OutputComponentType;

  typedef itk::ImageBase<ImageDimension>              ImageBaseType;

  /** Set the fixed image */
  itkNamedInputMacro(DeformationField, DeformationFieldType, "Primary")

  /** Set the moving image */
  itkNamedInputMacro(MovingImage, InputImageType, "moving")

  /**
   * Set whether the filter should use physical space calculations, i.e., the
   * displacement field is between physical coordinates of the voxels in the
   * fixed image domain and the physical coordinates of the voxels in the
   * moving image domain. By default, voxel units are used for greater speed.
   */
  itkSetMacro(UsePhysicalSpace, bool)
  itkGetMacro(UsePhysicalSpace, bool)

  /**
   * Set optional scaling for the deformation field
   */
  itkSetMacro(DeformationScaling, DeforamtionScalarType)
  itkGetMacro(DeformationScaling, DeforamtionScalarType)

  /**
   * Should nearest neighbor interpolation be used?
   */
  itkSetMacro(UseNearestNeighbor, bool)
  itkGetMacro(UseNearestNeighbor, bool)

  /**
   * When the warp is slightly outside of the moving image, should we
   * perform linear interpolation between border pixels and outside value (0)
   * or just report the outside value
   */
  itkSetMacro(ExtrapolateBorders, bool)
  itkGetMacro(ExtrapolateBorders, bool)

  /**
   * The outside value is used to fill in parts of the image that lie outside of
   * the domain.
   */
  itkSetMacro(OutsideValue, OutputComponentType);
  itkGetMacro(OutsideValue, OutputComponentType);

protected:

  FastWarpCompositeImageFilter()
  : m_UsePhysicalSpace(false), m_DeformationScaling(1.0),
    m_UseNearestNeighbor(false), m_ExtrapolateBorders(true), m_OutsideValue(0.0) { }

  ~FastWarpCompositeImageFilter() {}

  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                            itk::ThreadIdType threadId ) ITK_OVERRIDE;

  virtual void VerifyInputInformation() ITK_OVERRIDE {}

  virtual void GenerateInputRequestedRegion() ITK_OVERRIDE;

  virtual void GenerateOutputInformation() ITK_OVERRIDE;

  bool m_UsePhysicalSpace, m_UseNearestNeighbor, m_ExtrapolateBorders;
  DeforamtionScalarType m_DeformationScaling;

  OutputComponentType m_OutsideValue;


private:
  FastWarpCompositeImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "FastWarpCompositeImageFilter.txx"
#endif

#endif // FASTWARPCOMPOSITEIMAGEFILTER_H
