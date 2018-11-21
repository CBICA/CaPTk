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
#ifndef LieBracketFilter_h
#define LieBracketFilter_h

#include "itkImageToImageFilter.h"
#include "itkCovariantVector.h"
#include "itkImageRegionIterator.h"


/** 
 * \class LieBracketFilter
 * \brief Computes the Lie Bracket of two vector fields.
 *
 * Given vector field U and vector field V, this computes J(u) v - J(v) u
 *
 * Optionally, the result can be added to a third input, x with scaling alpha
 */
template< typename TInputImage, typename TOutputImage>
class ITK_TEMPLATE_EXPORT LieBracketFilter:
  public itk::ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Extract dimension from input image. */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage                       InputImageType;
  typedef typename InputImageType::Pointer  InputImagePointer;
  typedef TOutputImage                      OutputImageType;
  typedef typename OutputImageType::Pointer OutputImagePointer;

  /** Standard class typedefs. */
  typedef LieBracketFilter                                           Self;
  typedef itk::ImageToImageFilter< InputImageType, OutputImageType > Superclass;
  typedef itk::SmartPointer< Self >                                  Pointer;
  typedef itk::SmartPointer< const Self >                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(LieBracketFilter, ImageToImageFilter);

  /** Image typedef support. */
  typedef typename InputImageType::PixelType  InputPixelType;
  typedef typename OutputImageType::PixelType OutputPixelType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  /** First input, image U */
  itkNamedInputMacro(FieldU, InputImageType, "Primary");

  /** Second input, image V */
  itkNamedInputMacro(FieldV, InputImageType, "FieldV");

  /** Third optional input, additive image */
  itkNamedInputMacro(FieldX, InputImageType, "FieldX");

  /** LieBracketFilter needs a larger input requested region than
   * the output requested region.  As such, LieBracketFilter needs
   * to provide an implementation for GenerateInputRequestedRegion()
   * in order to inform the pipeline execution model.
   *
   * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
  virtual void GenerateInputRequestedRegion() ITK_OVERRIDE;

protected:
  LieBracketFilter();
  virtual ~LieBracketFilter();

  void ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread,
                            itk::ThreadIdType threadId) ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(LieBracketFilter);
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "LieBracketFilter.hxx"
#endif

#endif
