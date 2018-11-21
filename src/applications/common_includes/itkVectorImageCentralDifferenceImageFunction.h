/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkVectorImageCentralDifferenceImageFunction_h
#define itkVectorImageCentralDifferenceImageFunction_h

#include "itkImageFunction.h"
#include "itkMatrix.h"

namespace itk
{

/**
 * \class VectorImageCentralDifferenceImageFunction
 * \brief Calculate the derivative by central differencing.
 *
 * This class is templated over the input image type and
 * the coordinate representation type (e.g. float or double).
 *
 * Possible improvements:
 * - the use of Neighborhood operators may improve efficiency.
 *
 * \author Tom Vercauteren, INRIA & Mauna Kea Technologies
 *
 * This implementation was taken from the Insight Journal paper:
 * http://hdl.handle.net/1926/510
 *
 * \ingroup ImageFunctions
 * \ingroup ITKReview
 *
 * ADAPTED by P Yushkevich to VectorImage
 */
template<
  typename TInputImage,
  typename TCoordRep = float >
class VectorImageCentralDifferenceImageFunction:
  public ImageFunction< TInputImage,
                        typename TInputImage::PixelType,
                        TCoordRep >
{
public:
  typedef typename TInputImage::PixelType InputPixelType;

  /** Dimension underlying input image. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Standard class typedefs. */
  typedef VectorImageCentralDifferenceImageFunction Self;
  typedef ImageFunction< TInputImage,
                         typename TInputImage::PixelType,
                         TCoordRep >       Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(VectorImageCentralDifferenceImageFunction, ImageFunction);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** InputImageType typedef support. */
  typedef TInputImage InputImageType;

  /** OutputType typdef support. */
  typedef typename Superclass::OutputType OutputType;

  /** Index typedef support. */
  typedef typename Superclass::IndexType IndexType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

  /** Point typedef support. */
  typedef typename Superclass::PointType PointType;

  /** Evalulate the image derivative by central differencing at specified index.
   *
   *  No bounds checking is done.
   *  The point is assume to lie within the image buffer.
   *
   *  ImageFunction::IsInsideBuffer() can be used to check bounds before
   *  calling the method. */
  virtual void EvaluateAtIndex(const IndexType & index, OutputType &) const;

  virtual OutputType EvaluateAtIndex(const IndexType & index) const
  {
    OutputType res(ImageDimension * this->GetInputImage()->GetNumberOfComponentsPerPixel());
    this->EvaluateAtIndex(index, res);
    return res;
  }

  /** Evalulate the image derivative by central differencing at non-integer
   *  positions.
   *
   *  No bounds checking is done.
   *  The point is assumed to lie within the image buffer.
   *
   *  ImageFunction::IsInsideBuffer() can be used to check bounds before
   *  calling the method. */
  virtual OutputType Evaluate(const PointType & point) const
  {
    IndexType index;

    this->ConvertPointToNearestIndex(point, index);
    return this->EvaluateAtIndex(index);
  }

  virtual OutputType EvaluateAtContinuousIndex(
    const ContinuousIndexType & cindex) const
  {
    IndexType index;

    this->ConvertContinuousIndexToNearestIndex(cindex, index);
    return this->EvaluateAtIndex(index);
  }

protected:
  VectorImageCentralDifferenceImageFunction();
  ~VectorImageCentralDifferenceImageFunction(){}
  void PrintSelf(std::ostream & os, Indent indent) const;

private:
  VectorImageCentralDifferenceImageFunction(const Self &); //purposely not
                                                      // implemented
  void operator=(const Self &);                       //purposely not
                                                      // implemented

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVectorImageCentralDifferenceImageFunction.hxx"
#endif

#endif
