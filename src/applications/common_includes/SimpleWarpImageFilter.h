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
#ifndef __SimpleWarpImageFilter_h
#define __SimpleWarpImageFilter_h
#include "itkImageBase.h"
#include "itkImageFunction.h"
#include "itkImageToImageFilter.h"
#include "itkPoint.h"
#include "itkFixedArray.h"

namespace itk
{

/** \class SimpleWarpImageFilter
 * \brief Warps an image using an input deformation field (for LDDMM)
 *
 * SimpleWarpImageFilter warps an existing image with respect to
 * a given deformation field.
 *
 * A deformation field is represented as a image whose pixel type is some
 * vector type with at least N elements, where N is the dimension of
 * the input image. The vector type must support element access via operator
 * [].
 *
 * The output image is produced by inverse mapping: the output pixels
 * are mapped back onto the input image. This scheme avoids the creation of
 * any holes and overlaps in the output image.
 *
 * Each vector in the deformation field represent the distance between
 * a geometric point in the input space and a point in the output space 
 * in VOXEL COORDINATES (why? because it's faster!)
 *
 * \f[ p_{in} = p_{out} + d \f]
 *
 * Linear interpolation is used
 *
 * Position mapped to outside of the input image buffer are assigned
 * a edge padding value.
 *
 * This class is templated over the type of the input image, the
 * type of the output image and the type of the deformation field.
 *
 * The input image is set via SetInput. The input deformation field
 * is set via SetDeformationField.
 *
 * This filter is implemented as a multithreaded filter.
 *
 * \warning This filter assumes that the input type, output type
 * and deformation field type all have the same number of dimensions.
 *
 */
template <
  class TInputImage,
  class TOutputImage,
  class TDeformationField,
  class TFloat
  >
class ITK_EXPORT SimpleWarpImageFilter :
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef SimpleWarpImageFilter                              Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( SimpleWarpImageFilter, ImageToImageFilter );

  /** Typedef to describe the output image region type. */
  typedef typename TOutputImage::RegionType OutputImageRegionType;

  /** Inherit some types from the superclass. */
  typedef typename Superclass::InputImageType         InputImageType;
  typedef typename Superclass::InputImagePointer      InputImagePointer;
  typedef typename Superclass::OutputImageType        OutputImageType;
  typedef typename Superclass::OutputImagePointer     OutputImagePointer;
  typedef typename Superclass::InputImageConstPointer InputImageConstPointer;
  typedef typename OutputImageType::IndexType         IndexType;
  typedef typename OutputImageType::IndexValueType    IndexValueType;
  typedef typename OutputImageType::SizeType          SizeType;
  typedef typename OutputImageType::PixelType         PixelType;
  typedef typename OutputImageType::SpacingType       SpacingType;

  /** Determine the image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension );
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension );
  itkStaticConstMacro(DeformationFieldDimension, unsigned int,
                      TDeformationField::ImageDimension );
  /** typedef for base image type at the current ImageDimension */
  typedef ImageBase<itkGetStaticConstMacro(ImageDimension)> ImageBaseType;

  /** Deformation field typedef support. */
  typedef TDeformationField                        DeformationFieldType;
  typedef typename DeformationFieldType::Pointer   DeformationFieldPointer;
  typedef typename DeformationFieldType::PixelType DisplacementType;

  /** Type for representing the direction of the output image */
  typedef typename TOutputImage::DirectionType     DirectionType;

  /** Set the deformation field. */
  void SetDeformationField( const DeformationFieldType * field );

  /** Get a pointer the deformation field. */
  DeformationFieldType * GetDeformationField(void);

  /** Set the edge padding value */
  itkSetMacro( EdgePaddingValue, PixelType );

  /** Get the edge padding value */
  itkGetConstMacro( EdgePaddingValue, PixelType );

  /** Set scaling factor for the deformation field */
  itkSetMacro(DeformationScaling, TFloat);
  itkGetMacro(DeformationScaling, TFloat);

  /** SimpleWarpImageFilter produces an image which is a different
   * size than its input image. As such, it needs to provide an
   * implemenation for GenerateOutputInformation() which set
   * the output information according the OutputSpacing, OutputOrigin
   * and the deformation field's LargestPossibleRegion. */
  virtual void GenerateOutputInformation();

  /** It is difficult to compute in advance the input image region
   * required to compute the requested output region. Thus the safest
   * thing to do is to request for the whole input image.
   *
   * For the deformation field, the input requested region
   * set to be the same as that of the output requested region. */
  virtual void GenerateInputRequestedRegion();

  /** This method is used to set the state of the filter before 
   * multi-threading. */
  virtual void BeforeThreadedGenerateData();

  /** This method is used to set the state of the filter after 
   * multi-threading. */
  virtual void AfterThreadedGenerateData();

  typedef itk::ImageFunction<
    TInputImage,
    typename itk::NumericTraits< typename TInputImage::PixelType >::RealType,
    TFloat> InterpolatorType;

  itkSetObjectMacro(Interpolator, InterpolatorType);
  itkGetObjectMacro(Interpolator, InterpolatorType);

protected:
  SimpleWarpImageFilter();
  ~SimpleWarpImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** SimpleWarpImageFilter is implemented as a multi-threaded filter.
   * As such, it needs to provide and implementation for 
   * ThreadedGenerateData(). */
  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                            ThreadIdType threadId );

private:
  SimpleWarpImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  PixelType                  m_EdgePaddingValue;

  // Scaling for the deformation field
  TFloat m_DeformationScaling;

  // Interpolator
  typename InterpolatorType::Pointer m_Interpolator;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "SimpleWarpImageFilter.txx"
#endif

#endif
