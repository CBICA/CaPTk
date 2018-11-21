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
#ifndef JACOBIANDETERMINANTIMAGEFILTER_H
#define JACOBIANDETERMINANTIMAGEFILTER_H

#include <itkImageToImageFilter.h>

/** 
 * A class that computes a Jacobian determinant field from a warp field
 */
template <class TInputImage, class TOutputImage>
class JacobianDeterminantImageFilter
        : public itk::ImageToImageFilter<TInputImage, TOutputImage>
{
public:

  typedef JacobianDeterminantImageFilter<TInputImage,TOutputImage> Self;
  typedef itk::ImageToImageFilter<TInputImage, TOutputImage>       Superclass;
  typedef itk::SmartPointer<Self>                                  Pointer;
  typedef itk::SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self)

  /** Run-time type information (and related methods) */
  itkTypeMacro( JacobianDeterminantImageFilter, ImageToImageFilter )

  /** Determine the image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension );

  // Lots of typedefs
  typedef TInputImage                                 InputImageType;
  typedef TOutputImage                                OutputImageType;
  typedef typename InputImageType::RegionType         OutputImageRegionType;
  typedef typename InputImageType::PixelType          InputPixelType;
  typedef typename InputImageType::InternalPixelType  InputComponentType;
  typedef typename InputImageType::IndexType          IndexType;
  typedef typename InputImageType::IndexValueType     IndexValueType;
  typedef typename InputImageType::SizeType           SizeType;
  typedef typename InputImageType::SpacingType        SpacingType;
  typedef typename InputImageType::DirectionType      DirectionType;
  typedef typename OutputImageType::PixelType         OutputPixelType;
  typedef typename OutputImageType::InternalPixelType OutputComponentType;

  typedef itk::ImageBase<ImageDimension>              ImageBaseType;


protected:

  JacobianDeterminantImageFilter() {}
  ~JacobianDeterminantImageFilter() {}

  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                            itk::ThreadIdType threadId );
private:

  JacobianDeterminantImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "JacobianDeterminantImageFilter.txx"
#endif

#endif // JACOBIANDETERMINANTIMAGEFILTER_H
