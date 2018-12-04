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

#ifndef LINEARTRANSFORMTOWARPFILTER_H
#define LINEARTRANSFORMTOWARPFILTER_H

#include <itkImageToImageFilter.h>
#include "lddmm_common.h"

/**
 * This class transforms a linear transform into a deformation field
 */
template <class TInputImage, class TDeformationField, class TTransform>
class LinearTransformToWarpFilter
    : public itk::ImageToImageFilter<TInputImage,TDeformationField>
{
public:

  /** Standard class typedefs. */
  typedef LinearTransformToWarpFilter<TInputImage,TDeformationField,TTransform> Self;
  typedef itk::ImageToImageFilter<TInputImage,TDeformationField>   Superclass;
  typedef itk::SmartPointer<Self>                                  Pointer;
  typedef itk::SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self)

  /** Run-time type information (and related methods) */
  itkTypeMacro( LinearTransformToWarpFilter, ImageToImageFilter )

  /** Determine the image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension );

  // Lots of typedefs
  typedef TInputImage                                 InputImageType;
  typedef TDeformationField                           DeformationFieldType;
  typedef typename InputImageType::RegionType         OutputImageRegionType;
  typedef typename InputImageType::PixelType          InputPixelType;
  typedef typename InputImageType::InternalPixelType  InputComponentType;
  typedef typename InputImageType::IndexType          IndexType;
  typedef typename InputImageType::IndexValueType     IndexValueType;
  typedef typename InputImageType::SizeType           SizeType;
  typedef typename InputImageType::SpacingType        SpacingType;
  typedef typename InputImageType::DirectionType      DirectionType;
  typedef typename DeformationFieldType::PixelType    DeformationVectorType;
  typedef itk::ImageBase<ImageDimension>              ImageBaseType;

  typedef TTransform                                  TransformType;

  /** Set the fixed image */
  itkNamedInputMacro(FixedImage, InputImageType, "Primary")

  /** Set the moving image */
  itkNamedInputMacro(MovingImage, InputImageType, "moving")

  /** Set the transform */
  itkSetObjectMacro(Transform, TransformType)

  /** Get the transform */
  itkGetObjectMacro(Transform, TransformType)

protected:

  LinearTransformToWarpFilter() {}
  virtual ~LinearTransformToWarpFilter() {}

  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                            itk::ThreadIdType threadId );

  typename TransformType::Pointer m_Transform;

private:
  LinearTransformToWarpFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "LinearTransformToWarpFilter.txx"
#endif

#endif // LINEARTRANSFORMTOWARPFILTER_H

