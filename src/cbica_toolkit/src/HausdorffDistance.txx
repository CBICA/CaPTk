/**

\brief Implementation taken from https://github.com/InsightSoftwareConsortium/covalic

*/

#ifndef _HausdorffDistanceImageToImageMetric_txx
#define _HausdorffDistanceImageToImageMetric_txx

#include "itkBSplineInterpolateImageFunction.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkMorphologicalGradientImageFilter.h"

#include "vnl/vnl_math.h"

// #include "HausdorffDistanceImageToImageMetric.h"

#include <algorithm>
#include <cmath>
#include <vector>


template <class TFixedImage, class TMovingImage>
HausdorffDistanceImageToImageMetric<TFixedImage, TMovingImage>
::HausdorffDistanceImageToImageMetric()
{
  m_DoBlurring = false;

  m_Percentile = 0.95;
}

template <class TFixedImage, class TMovingImage>
HausdorffDistanceImageToImageMetric<TFixedImage, TMovingImage>
::~HausdorffDistanceImageToImageMetric()
{

}

template <class TFixedImage, class TMovingImage>
void
HausdorffDistanceImageToImageMetric<TFixedImage, TMovingImage>
::SetPercentile(double p)
{
  if (p < 0.0 || p > 1.0)
    itkExceptionMacro(<< "Percentile needs to be in [0, 1]");

  m_Percentile = p;
}

template <class TFixedImage, class TMovingImage>
double
HausdorffDistanceImageToImageMetric<TFixedImage, TMovingImage>
::ComputeMaxDistance(
  const TFixedImage* img1, const TMovingImage* img2) const
{
  typedef itk::Image<float, FixedImageType::ImageDimension> FloatImageType;

  typedef itk::SignedMaurerDistanceMapImageFilter<
    FixedImageType, FloatImageType> DistanceMapFilterType;

  typedef itk::DiscreteGaussianImageFilter<
    FloatImageType, FloatImageType> BlurFilterType;

  // Compute distance transform
  typename FloatImageType::Pointer distMap2;
  {
    typename DistanceMapFilterType::Pointer distanceMapFilter =
      DistanceMapFilterType::New();

    distanceMapFilter->InsideIsPositiveOff();
    distanceMapFilter->SetInput(img2);
    distanceMapFilter->SquaredDistanceOff();
    distanceMapFilter->UseImageSpacingOn();

    distanceMapFilter->Update();

    if (!m_DoBlurring)
    {
      distMap2 = distanceMapFilter->GetOutput();
    }
    else
    {
      typename FloatImageType::SpacingType spacing = img1->GetSpacing();

      double minSpacing = spacing[0];
      for (unsigned int dim = 0; dim < FixedImageType::ImageDimension; dim++)
        if (spacing[dim] < minSpacing)
          minSpacing = spacing[dim];

      typename BlurFilterType::Pointer blurf = BlurFilterType::New();
      blurf->SetInput(distanceMapFilter->GetOutput());
      blurf->SetVariance(1.5 * minSpacing);
      blurf->Update();

      distMap2 = blurf->GetOutput();
    }
  }

  //typedef itk::LinearInterpolateImageFunction<FloatImageType, double>
  typedef itk::BSplineInterpolateImageFunction<FloatImageType, double>
    InterpolatorType;
  typename InterpolatorType::Pointer distInterp2 = InterpolatorType::New();
  distInterp2->SetInputImage(distMap2);
  distInterp2->SetSplineOrder(3);

  // Detect boundary via morphological gradient
  typedef itk::BinaryBallStructuringElement<FixedImagePixelType, TFixedImage::ImageDimension>
    StructElementType;
  typedef
    itk::MorphologicalGradientImageFilter<FixedImageType, FixedImageType,
      StructElementType> EdgeFilterType;

  StructElementType structel;
  structel.SetRadius(1);
  structel.CreateStructuringElement();

  typename EdgeFilterType::Pointer edgef = EdgeFilterType::New();
  edgef->SetInput(img1);
  edgef->SetKernel(structel);
  edgef->Update();

  FixedImagePointer edgeImg1 = edgef->GetOutput();

  std::vector<double> distances;

  typedef itk::ImageRegionConstIteratorWithIndex<FixedImageType>
    FixedIteratorType;
  FixedIteratorType it(edgeImg1, edgeImg1->GetLargestPossibleRegion());

  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    FixedImageIndexType ind = it.GetIndex();

    if (it.Get() == 0 || img1->GetPixel(ind) == 0)
      continue;

     FixedImagePointType p;
     img1->TransformIndexToPhysicalPoint(ind, p);

    if (!distInterp2->IsInsideBuffer(p))
      continue;

    distances.push_back(vnl_math_abs(distInterp2->Evaluate(p)));
  }

  if (distances.size() == 0)
    return vnl_huge_val(1.0);

  std::sort(distances.begin(), distances.end());

  return distances[(int)(m_Percentile*(distances.size() - 1))];
}

template <class TFixedImage, class TMovingImage>
typename HausdorffDistanceImageToImageMetric<TFixedImage, TMovingImage>::MeasureType
HausdorffDistanceImageToImageMetric<TFixedImage, TMovingImage>
::GetValue() const
{
  if (Superclass::m_FixedImage.IsNull() || Superclass::m_MovingImage.IsNull())
    itkExceptionMacro(<< "Need two input classification images");

  // Handle special case where inputs are zeros
  typedef itk::ImageRegionConstIteratorWithIndex<FixedImageType>
    FixedIteratorType;
  FixedIteratorType it(Superclass::m_FixedImage, Superclass::m_FixedImage->GetLargestPossibleRegion());

  double sumFixed = 0;
  double sumMoving = 0;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    sumFixed += it.Get();
    sumMoving += Superclass::m_MovingImage->GetPixel(it.GetIndex());
  }

  if (sumFixed == 0 || sumMoving == 0)
  {
    if (sumFixed == sumMoving)
      return 0.0;
    else
      return vnl_huge_val(1.0);
  }

  // Compute max distances at specified percentile
  double d12 = this->ComputeMaxDistance(
    Superclass::m_FixedImage, Superclass::m_MovingImage);
  double d21 = this->ComputeMaxDistance(
    Superclass::m_MovingImage, Superclass::m_FixedImage);

  if (d12 > d21)
    return d12;
  else
    return d21;

}

#endif