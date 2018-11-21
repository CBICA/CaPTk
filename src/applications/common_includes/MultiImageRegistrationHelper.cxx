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
#include "MultiImageRegistrationHelper.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNumericTraits.h"
#include "itkContinuousIndex.h"
#include "vnl/vnl_math.h"
#include "lddmm_data.h"
#include "MultiImageOpticalFlowImageFilter.h"
#include "MultiComponentNCCImageMetric.h"
#include "MultiComponentApproximateNCCImageMetric.h"
#include "MultiComponentMutualInfoImageMetric.h"
#include "MahalanobisDistanceToTargetWarpMetric.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "OneDimensionalInPlaceAccumulateFilter.h"
#include "itkUnaryFunctorImageFilter.h"
#include "itkImageFileWriter.h"
#include "GreedyException.h"

template <class TFloat, unsigned int VDim>
void
MultiImageOpticalFlowHelper<TFloat, VDim>
::SetDefaultPyramidFactors(int n_levels)
{
  for(int i = n_levels-1; i>=0; --i)
    m_PyramidFactors.push_back(1 << i);
}

template <class TFloat, unsigned int VDim>
void
MultiImageOpticalFlowHelper<TFloat, VDim>
::SetPyramidFactors(const PyramidFactorsType &factors)
{
  m_PyramidFactors = factors;
}

template <class TFloat, unsigned int VDim>
void
MultiImageOpticalFlowHelper<TFloat, VDim>
::AddImagePair(MultiComponentImageType *fixed, MultiComponentImageType *moving, double weight)
{
  // Collect the weights
  for(unsigned i = 0; i < fixed->GetNumberOfComponentsPerPixel(); i++)
    m_Weights.push_back(weight);

  // Store the images
  m_Fixed.push_back(fixed);
  m_Moving.push_back(moving);
}

template <class TFloat, unsigned int VDim>
void
MultiImageOpticalFlowHelper<TFloat, VDim>
::SetJitterSigma(double sigma)
{
  m_JitterSigma = sigma;
}

template <class TFloat, unsigned int VDim>
void
MultiImageOpticalFlowHelper<TFloat, VDim>
::PlaceIntoComposite(FloatImageType *source, MultiComponentImageType *target, int offset)
{
  // We do this using a loop - no threading
  TFloat *src_ptr = source->GetPixelContainer()->GetBufferPointer();
  TFloat *trg_ptr = target->GetPixelContainer()->GetBufferPointer() + offset;

  int trg_comp = target->GetNumberOfComponentsPerPixel();

  int n_voxels = source->GetPixelContainer()->Size();
  TFloat *trg_end = trg_ptr + n_voxels * target->GetNumberOfComponentsPerPixel();

  while(trg_ptr < trg_end)
    {
    *trg_ptr = *src_ptr++;
    trg_ptr += trg_comp;
    }
}

template <class TFloat, unsigned int VDim>
void
MultiImageOpticalFlowHelper<TFloat, VDim>
::PlaceIntoComposite(VectorImageType *source, MultiComponentImageType *target, int offset)
{
  // We do this using a loop - no threading
  VectorType *src_ptr = source->GetPixelContainer()->GetBufferPointer();
  TFloat *trg_ptr = target->GetPixelContainer()->GetBufferPointer() + offset;

  int trg_skip = target->GetNumberOfComponentsPerPixel() - VDim;

  int n_voxels = source->GetPixelContainer()->Size();
  TFloat *trg_end = trg_ptr + n_voxels * target->GetNumberOfComponentsPerPixel();

  while(trg_ptr < trg_end)
    {
    const VectorType &vsrc = *src_ptr++;
    for(int k = 0; k < VDim; k++)
      *trg_ptr++ = vsrc[k];
    trg_ptr += trg_skip;
    }
}

#include <vnl/vnl_random.h>

template <class TFloat, unsigned int VDim>
void
MultiImageOpticalFlowHelper<TFloat, VDim>
::DilateCompositeGradientMasksForNCC(SizeType radius)
{
  typedef LDDMMData<TFloat, VDim> LDDMMType;

  for(int level = 0; level < m_PyramidFactors.size(); level++)
    {
    if(m_GradientMaskComposite[level])
      {
      // Threshold the mask itself
      LDDMMType::img_threshold_in_place(m_GradientMaskComposite[level], 0.5, 1e100, 0.5, 0);

      // Make a copy of the mask
      typename FloatImageType::Pointer mask_copy;
      LDDMMType::alloc_img(mask_copy, m_GradientMaskComposite[level]);
      LDDMMType::img_copy(m_GradientMaskComposite[level], mask_copy);

      // Run the accumulation filter on the mask
      typename FloatImageType::Pointer mask_accum =
          AccumulateNeighborhoodSumsInPlace(mask_copy.GetPointer(), radius);

      // Threshold the mask copy
      LDDMMType::img_threshold_in_place(mask_accum, 0.25, 1e100, 0.5, 0);

      // Add the two images - the result has 1 for the initial mask, 0.5 for the 'outer' mask
      LDDMMType::img_add_in_place(m_GradientMaskComposite[level], mask_accum);
      }
    }
}


template <class TFloat, unsigned int VDim>
void
MultiImageOpticalFlowHelper<TFloat, VDim>
::BuildCompositeImages(double noise_sigma_relative)
{
  typedef LDDMMData<TFloat, VDim> LDDMMType;

  // Offsets into the composite images
  int off_fixed = 0, off_moving = 0;

  // Set up the composite images
  m_FixedComposite.resize(m_PyramidFactors.size());
  m_MovingComposite.resize(m_PyramidFactors.size());

  // Repeat for each of the input images
  for(size_t j = 0; j < m_Fixed.size(); j++)
    {
    // Repeat for each component
    for(unsigned k = 0; k < m_Fixed[j]->GetNumberOfComponentsPerPixel(); k++)
      {
      // Extract the k-th image component from fixed and moving images
      typedef itk::VectorIndexSelectionCastImageFilter<MultiComponentImageType, FloatImageType> ExtractType;
      typename ExtractType::Pointer fltExtractFixed, fltExtractMoving;

      fltExtractFixed = ExtractType::New();
      fltExtractFixed->SetInput(m_Fixed[j]);
      fltExtractFixed->SetIndex(k);
      fltExtractFixed->Update();

      fltExtractMoving = ExtractType::New();
      fltExtractMoving->SetInput(m_Moving[j]);
      fltExtractMoving->SetIndex(k);
      fltExtractMoving->Update();

      // Deal with additive noise
      double noise_sigma_fixed = 0.0, noise_sigma_moving = 0.0;

      if(noise_sigma_relative > 0.0)
        {
        // Figure out the quartiles of the fixed image
        typedef MutualInformationPreprocessingFilter<FloatImageType, FloatImageType> QuantileFilter;
        typename QuantileFilter::Pointer fltQuantileFixed = QuantileFilter::New();
        fltQuantileFixed->SetLowerQuantile(0.01);
        fltQuantileFixed->SetUpperQuantile(0.99);
        fltQuantileFixed->SetInput(fltExtractFixed->GetOutput());
        fltQuantileFixed->Update();
        double range_fixed = fltQuantileFixed->GetUpperQuantileValue(0) - fltQuantileFixed->GetLowerQuantileValue(0);
        noise_sigma_fixed = noise_sigma_relative * range_fixed;

        // Figure out the quartiles of the moving image
        typename QuantileFilter::Pointer fltQuantileMoving = QuantileFilter::New();
        fltQuantileMoving->SetLowerQuantile(0.01);
        fltQuantileMoving->SetUpperQuantile(0.99);
        fltQuantileMoving->SetInput(fltExtractMoving->GetOutput());
        fltQuantileMoving->Update();
        double range_moving = fltQuantileMoving->GetUpperQuantileValue(0) - fltQuantileMoving->GetLowerQuantileValue(0);
        noise_sigma_moving = noise_sigma_relative * range_moving;

        // Report noise levels
        //printf("Noise on image %d component %d: fixed = %g, moving = %g\n", static_cast<int>(j), k, noise_sigma_fixed, noise_sigma_moving);
        }

      // Compute the pyramid for this component
      for(size_t i = 0; i < m_PyramidFactors.size(); i++)
        {
        // Downsample the image to the right pyramid level
        typename FloatImageType::Pointer lFixed, lMoving;
        if (m_PyramidFactors[i] == 1)
          {
          lFixed = fltExtractFixed->GetOutput();
          lMoving = fltExtractMoving->GetOutput();
          }
        else
          {
          lFixed = FloatImageType::New();
          lMoving = FloatImageType::New();
          LDDMMType::img_downsample(fltExtractFixed->GetOutput(), lFixed, m_PyramidFactors[i]);
          LDDMMType::img_downsample(fltExtractMoving->GetOutput(), lMoving, m_PyramidFactors[i]);

          // For the Mahalanobis metric, the fixed image needs to be scaled by the factor of the 
          // pyramid level because it describes voxel coordinates
          if(m_ScaleFixedImageWithVoxelSize)
            LDDMMType::img_scale_in_place(lFixed, 1.0 / m_PyramidFactors[i]);
          }

        // Add some noise to the images
        if(noise_sigma_relative > 0.0)
          {
          vnl_random randy(12345);
          for(size_t i = 0; i < lFixed->GetPixelContainer()->Size(); i++)
            lFixed->GetBufferPointer()[i] += randy.normal() * noise_sigma_fixed;
          for(size_t i = 0; i < lMoving->GetPixelContainer()->Size(); i++)
            lMoving->GetBufferPointer()[i] += randy.normal() * noise_sigma_moving;
          }

        // Compute the gradient of the moving image
        //typename VectorImageType::Pointer gradMoving = VectorImageType::New();
        //LDDMMType::alloc_vimg(gradMoving, lMoving);
        //LDDMMType::image_gradient(lMoving, gradMoving);

        // Allocate the composite images if they have not been allocated
        if(j == 0 && k == 0)
          {
          m_FixedComposite[i] = MultiComponentImageType::New();
          m_FixedComposite[i]->CopyInformation(lFixed);
          m_FixedComposite[i]->SetNumberOfComponentsPerPixel(m_Weights.size());
          m_FixedComposite[i]->SetRegions(lFixed->GetBufferedRegion());
          m_FixedComposite[i]->Allocate();

          m_MovingComposite[i] = MultiComponentImageType::New();
          m_MovingComposite[i]->CopyInformation(lMoving);
          m_MovingComposite[i]->SetNumberOfComponentsPerPixel(m_Weights.size());
          m_MovingComposite[i]->SetRegions(lMoving->GetBufferedRegion());
          m_MovingComposite[i]->Allocate();
          }

        // Pack the data into the fixed and moving composite images
        this->PlaceIntoComposite(lFixed, m_FixedComposite[i], off_fixed);
        this->PlaceIntoComposite(lMoving, m_MovingComposite[i], off_moving);
        }

      // Update the offsets
      off_fixed++;
      off_moving++;
      }
    }

  // Set up the mask pyramid
  m_GradientMaskComposite.resize(m_PyramidFactors.size(), NULL);
  if(m_GradientMaskImage)
    {
    for(size_t i = 0; i < m_PyramidFactors.size(); i++)
      {
      // Downsample the image to the right pyramid level
      if (m_PyramidFactors[i] == 1)
        {
        m_GradientMaskComposite[i] = m_GradientMaskImage;
        }
      else
        {
        m_GradientMaskComposite[i] = FloatImageType::New();

        // Downsampling the mask involves smoothing, so the mask will no longer be binary
        LDDMMType::img_downsample(m_GradientMaskImage, m_GradientMaskComposite[i], m_PyramidFactors[i]);
        LDDMMType::img_threshold_in_place(m_GradientMaskComposite[i], 0.5, 1e100, 1.0, 0.0);
        }      
      }
    }

  // Set up the moving mask pyramid
  m_MovingMaskComposite.resize(m_PyramidFactors.size(), NULL);
  if(m_MovingMaskImage)
    {
    for(size_t i = 0; i < m_PyramidFactors.size(); i++)
      {
      // Downsample the image to the right pyramid level
      if (m_PyramidFactors[i] == 1)
        {
        m_MovingMaskComposite[i] = m_MovingMaskImage;
        }
      else
        {
        m_MovingMaskComposite[i] = FloatImageType::New();

        // Downsampling the mask involves smoothing, so the mask will no longer be binary
        LDDMMType::img_downsample(m_MovingMaskImage, m_MovingMaskComposite[i], m_PyramidFactors[i]);

        // We don't need the moving mask to be binary, we can leave it be floating point...
        // LDDMMType::img_threshold_in_place(m_MovingMaskComposite[i], 0.5, 1e100, 1.0, 0.0);
        }
      }
    }

  // Set up the jitter images
  m_JitterComposite.resize(m_PyramidFactors.size(), NULL);
  if(m_JitterSigma > 0)
    {
    for(size_t i = 0; i < m_PyramidFactors.size(); i++)
      {
      // Get the reference space
      ImageBaseType *base = this->GetReferenceSpace(i);
      VectorImagePointer iJitter = VectorImageType::New();
      iJitter->CopyInformation(base);
      iJitter->SetRegions(base->GetBufferedRegion());
      iJitter->Allocate();

      vnl_random randy(12345);
      typedef itk::ImageRegionIterator<VectorImageType> IterType;
      for(IterType iter(iJitter, iJitter->GetBufferedRegion()); !iter.IsAtEnd(); ++iter)
        {
        for(int k = 0; k < VDim; k++)
          {
          iter.Value()[k] = randy.normal() * m_JitterSigma;
          }
        }

      m_JitterComposite[i] = iJitter;
      }
    }
}

template <class TFloat, unsigned int VDim>
typename MultiImageOpticalFlowHelper<TFloat, VDim>::ImageBaseType *
MultiImageOpticalFlowHelper<TFloat, VDim>
::GetMovingReferenceSpace(int level)
{
  return m_MovingComposite[level];
}

template <class TFloat, unsigned int VDim>
typename MultiImageOpticalFlowHelper<TFloat, VDim>::ImageBaseType *
MultiImageOpticalFlowHelper<TFloat, VDim>
::GetReferenceSpace(int level)
{
  return m_FixedComposite[level];
}

template <class TFloat, unsigned int VDim>
typename MultiImageOpticalFlowHelper<TFloat, VDim>::Vec
MultiImageOpticalFlowHelper<TFloat, VDim>
::GetSmoothingSigmasInPhysicalUnits(int level, double sigma, bool in_physical_units)
{
  Vec sigmas;
  if(in_physical_units)
    {
    sigmas.Fill(sigma * m_PyramidFactors[level]);
    }
  else
    {
    for(int k = 0; k < VDim; k++)
      sigmas[k] = this->GetReferenceSpace(level)->GetSpacing()[k] * sigma;
    }
  return sigmas;
}

template <class TFloat, unsigned int VDim>
vnl_vector<double>
MultiImageOpticalFlowHelper<TFloat, VDim>
::ComputeOpticalFlowField(int level,
                          VectorImageType *def,
                          FloatImageType *out_metric,
                          VectorImageType *out_gradient,
                          double result_scaling)
{
  typedef DefaultMultiComponentImageMetricTraits<TFloat, VDim> TraitsType;
  typedef MultiImageOpticalFlowImageFilter<TraitsType> FilterType;

  typename FilterType::Pointer filter = FilterType::New();

  // Scale the weights by epsilon
  vnl_vector<float> wscaled(m_Weights.size());
  for (unsigned i = 0; i < wscaled.size(); i++)
    wscaled[i] = m_Weights[i] * result_scaling;

  // Run the filter
  filter->SetFixedImage(m_FixedComposite[level]);
  filter->SetMovingImage(m_MovingComposite[level]);
  filter->SetDeformationField(def);
  filter->SetWeights(wscaled);
  filter->SetComputeGradient(true);
  filter->GetMetricOutput()->Graft(out_metric);
  filter->GetDeformationGradientOutput()->Graft(out_gradient);
  filter->Update();

  // Get the vector of the normalized metrics
  return filter->GetAllMetricValues();
}

template <class TFloat, unsigned int VDim>
vnl_vector<double>
MultiImageOpticalFlowHelper<TFloat, VDim>
::ComputeMIFlowField(int level,
                     bool normalized_mutual_information,
                     VectorImageType *def,
                     FloatImageType *out_metric,
                     VectorImageType *out_gradient,
                     double result_scaling)
{
  // Scale the weights by epsilon
  vnl_vector<float> wscaled(m_Weights.size());
  for (unsigned i = 0; i < wscaled.size(); i++)
    wscaled[i] = m_Weights[i] * result_scaling;

  // Set up the mutual information metric
  typedef DefaultMultiComponentMutualInfoImageMetricTraits<TFloat, unsigned char, VDim> TraitsType;
  typedef MultiComponentMutualInfoImageMetric<TraitsType> MetricType;

  typedef itk::VectorImage<unsigned char, VDim> BinnedImageType;
  typedef MutualInformationPreprocessingFilter<MultiComponentImageType, BinnedImageType> BinnerType;

  // TODO: this is utter laziness, get rid of this garbage!
  static typename BinnerType::Pointer binner_fixed;
  static typename BinnerType::Pointer binner_moving;

  if(binner_fixed.IsNull()
     || binner_fixed->GetOutput()->GetBufferedRegion()
     != m_FixedComposite[level]->GetBufferedRegion())
    {
    binner_fixed = BinnerType::New();
    binner_fixed->SetInput(m_FixedComposite[level]);
    binner_fixed->SetBins(128);
    binner_fixed->SetLowerQuantile(0.01);
    binner_fixed->SetUpperQuantile(0.99);
    binner_fixed->SetStartAtBinOne(true);
    binner_fixed->Update();

    binner_moving = BinnerType::New();
    binner_moving->SetInput(m_MovingComposite[level]);
    binner_moving->SetBins(128);
    binner_moving->SetLowerQuantile(0.01);
    binner_moving->SetUpperQuantile(0.99);
    binner_moving->SetStartAtBinOne(true);
    binner_moving->Update();
    }


  typename MetricType::Pointer metric = MetricType::New();

  metric->SetComputeNormalizedMutualInformation(normalized_mutual_information);
  metric->SetFixedImage(binner_fixed->GetOutput());
  metric->SetMovingImage(binner_moving->GetOutput());
  metric->SetDeformationField(def);
  metric->SetWeights(wscaled);
  metric->SetComputeGradient(true);
  metric->GetMetricOutput()->Graft(out_metric);
  metric->GetDeformationGradientOutput()->Graft(out_gradient);
  metric->GetMetricOutput()->Graft(out_metric);
  metric->SetBins(128);
  metric->Update();

  // Process the results
  return metric->GetAllMetricValues();
}

// #undef DUMP_NCC
#define DUMP_NCC 1


template <class TFloat, unsigned int VDim>
typename MultiImageOpticalFlowHelper<TFloat, VDim>::SizeType
MultiImageOpticalFlowHelper<TFloat, VDim>
::AdjustNCCRadius(int level, const SizeType &radius, bool report_on_adjust)
{
  SizeType radius_fix = radius;
  for(int d = 0; d < VDim; d++)
    {
    int sz_d = (int) m_FixedComposite[level]->GetBufferedRegion().GetSize()[d];
    if(radius_fix[d] * 2 + 1 >= sz_d)
      radius_fix[d] = (sz_d - 1) / 2;
    }

  if(report_on_adjust && radius != radius_fix)
    {
    std::cout << "  *** NCC radius adjusted to " << radius_fix
              << " because image too small at level " << level
              << " (" << m_FixedComposite[level]->GetBufferedRegion().GetSize() << ")" << std::endl;
    }

  return radius_fix;
}

template <class TFloat, unsigned int VDim>
double
MultiImageOpticalFlowHelper<TFloat, VDim>
::ComputeNCCMetricImage(int level,
                        VectorImageType *def,
                        const SizeType &radius,
                        FloatImageType *out_metric,
                        VectorImageType *out_gradient,
                        double result_scaling)
{
  typedef DefaultMultiComponentImageMetricTraits<TFloat, VDim> TraitsType;
  typedef MultiComponentNCCImageMetric<TraitsType> FilterType;
  // typedef MultiComponentApproximateNCCImageMetric<TraitsType> FilterType;

  typename FilterType::Pointer filter = FilterType::New();

  // Scale the weights by epsilon
  vnl_vector<float> wscaled(m_Weights.size());
  for (unsigned i = 0; i < wscaled.size(); i++)
    wscaled[i] = m_Weights[i] * result_scaling;

  // Allocate a working image
  if(m_NCCWorkingImage.IsNull())
    m_NCCWorkingImage = MultiComponentImageType::New();

  // Is this the first time that this function is being called with this image?
  bool first_run =
      m_NCCWorkingImage->GetBufferedRegion() != m_FixedComposite[level]->GetBufferedRegion();

  // Check the radius against the size of the image
  SizeType radius_fix = AdjustNCCRadius(level, radius, first_run);

  // Run the filter
  filter->SetFixedImage(m_FixedComposite[level]);
  filter->SetMovingImage(m_MovingComposite[level]);
  filter->SetDeformationField(def);
  filter->SetWeights(wscaled);
  filter->SetComputeGradient(true);
  filter->GetMetricOutput()->Graft(out_metric);
  filter->GetDeformationGradientOutput()->Graft(out_gradient);
  filter->SetRadius(radius_fix);
  filter->SetWorkingImage(m_NCCWorkingImage);
  filter->SetReuseWorkingImageFixedComponents(!first_run);
  filter->SetFixedMaskImage(m_GradientMaskComposite[level]);

  // TODO: support moving masks...
  // filter->SetMovingMaskImage(m_MovingMaskComposite[level]);
  filter->Update();

  // Get the vector of the normalized metrics
  return filter->GetMetricValue();
}

template <class TFloat, unsigned int VDim>
double
MultiImageOpticalFlowHelper<TFloat, VDim>
::ComputeMahalanobisMetricImage(int level, VectorImageType *def, 
  FloatImageType *out_metric, 
  VectorImageType *out_gradient)
{
  typedef DefaultMahalanobisDistanceToTargetMetricTraits<TFloat, VDim> TraitsType;
  typedef MahalanobisDistanceToTargetWarpMetric<TraitsType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  filter->SetFixedImage(m_FixedComposite[level]);
  filter->SetMovingImage(m_MovingComposite[level]);
  filter->SetDeformationField(def);
  filter->SetComputeGradient(true);
  filter->GetMetricOutput()->Graft(out_metric);
  filter->GetDeformationGradientOutput()->Graft(out_gradient);
  filter->SetFixedMaskImage(m_GradientMaskComposite[level]);

  filter->Update();

  return filter->GetMetricValue();
}


template <class TFloat, unsigned int VDim>
double
MultiImageOpticalFlowHelper<TFloat, VDim>
::ComputeAffineMSDMatchAndGradient(int level,
    LinearTransformType *tran,
    FloatImageType *wrkMetric,
    FloatImageType *wrkMask,
    VectorImageType *wrkGradMetric,
    VectorImageType *wrkGradMask,
    VectorImageType *wrkPhi,
    LinearTransformType *grad)
{
  // Scale the weights by epsilon
  vnl_vector<float> wscaled(m_Weights.size());
  for (unsigned i = 0; i < wscaled.size(); i++)
    wscaled[i] = m_Weights[i];

  // Set up the optical flow computation
  typedef DefaultMultiComponentImageMetricTraits<TFloat, VDim> TraitsType;
  typedef MultiImageOpticalFlowImageFilter<TraitsType> MetricType;
  typename MetricType::Pointer metric = MetricType::New();

  metric->SetFixedImage(m_FixedComposite[level]);
  metric->SetMovingImage(m_MovingComposite[level]);
  metric->SetWeights(wscaled);
  metric->SetAffineTransform(tran);
  metric->SetComputeMovingDomainMask(true);
  metric->GetMetricOutput()->Graft(wrkMetric);
  metric->SetComputeGradient(grad != NULL);
  metric->SetFixedMaskImage(m_GradientMaskComposite[level]);
  metric->SetJitterImage(m_JitterComposite[level]);
  metric->Update();

  // TODO: erase this
  /*
  std::cout << "SAVING METRIC, TRAN = " << tran->GetMatrix() << std::endl;
  static int iter = 0;
  std::ostringstream oss; oss << "metric_" << iter << ".nii.gz";
  LDDMMData<TFloat, VDim>::img_write(wrkMetric, oss.str().c_str());
  std::stringstream oss2; oss2 << "metric_mask_" << iter << ".nii.gz";
  LDDMMData<TFloat, VDim>::img_write(wrkMask, oss2.str().c_str());
  ++iter;
  */

  // Process the results
  if(grad)
    {
    grad->SetMatrix(metric->GetAffineTransformGradient()->GetMatrix());
    grad->SetOffset(metric->GetAffineTransformGradient()->GetOffset());
    }

  return metric->GetMetricValue();
}

#include "itkRescaleIntensityImageFilter.h"

template <class TFloat, unsigned int VDim>
double
MultiImageOpticalFlowHelper<TFloat, VDim>
::ComputeAffineMIMatchAndGradient(int level,
                                  bool normalized_mutual_info,
                                  LinearTransformType *tran,
                                  FloatImageType *wrkMetric,
                                  FloatImageType *wrkMask,
                                  VectorImageType *wrkGradMetric,
                                  VectorImageType *wrkGradMask,
                                  VectorImageType *wrkPhi,
                                  LinearTransformType *grad)
{
  // Scale the weights by epsilon
  vnl_vector<float> wscaled(m_Weights.size());
  for (unsigned i = 0; i < wscaled.size(); i++)
    wscaled[i] = m_Weights[i];

  // Set up the mutual information metric
  typedef DefaultMultiComponentMutualInfoImageMetricTraits<TFloat, unsigned char, VDim> TraitsType;
  typedef MultiComponentMutualInfoImageMetric<TraitsType> MetricType;

  typedef itk::VectorImage<unsigned char, VDim> BinnedImageType;
  typedef MutualInformationPreprocessingFilter<MultiComponentImageType, BinnedImageType> BinnerType;

  // TODO: this is utter laziness, get rid of this garbage!
  static typename BinnerType::Pointer binner_fixed;
  static typename BinnerType::Pointer binner_moving;

  if(binner_fixed.IsNull()
     || binner_fixed->GetOutput()->GetBufferedRegion()
     != m_FixedComposite[level]->GetBufferedRegion())
    {
    binner_fixed = BinnerType::New();
    binner_fixed->SetInput(m_FixedComposite[level]);
    binner_fixed->SetBins(128);
    binner_fixed->SetLowerQuantile(0.01);
    binner_fixed->SetUpperQuantile(0.99);
    binner_fixed->SetStartAtBinOne(true);
    binner_fixed->Update();

    binner_moving = BinnerType::New();
    binner_moving->SetInput(m_MovingComposite[level]);
    binner_moving->SetBins(128);
    binner_moving->SetLowerQuantile(0.01);
    binner_moving->SetUpperQuantile(0.99);
    binner_moving->SetStartAtBinOne(true);
    binner_moving->Update();
    }

  typename MetricType::Pointer metric = MetricType::New();

  metric->SetComputeNormalizedMutualInformation(normalized_mutual_info);
  metric->SetFixedImage(binner_fixed->GetOutput());
  metric->SetMovingImage(binner_moving->GetOutput());
  metric->SetWeights(wscaled);
  metric->SetAffineTransform(tran);
  metric->SetComputeMovingDomainMask(true);
  metric->GetMetricOutput()->Graft(wrkMetric);
  metric->SetComputeGradient(grad != NULL);
  metric->SetFixedMaskImage(m_GradientMaskComposite[level]);
  metric->SetBins(128);
  metric->SetJitterImage(m_JitterComposite[level]);
  metric->Update();

  // Process the results
  if(grad)
    {
    grad->SetMatrix(metric->GetAffineTransformGradient()->GetMatrix());
    grad->SetOffset(metric->GetAffineTransformGradient()->GetOffset());
    }

  return metric->GetMetricValue();
}



// TODO: there is a lot of code duplication here!
template <class TFloat, unsigned int VDim>
double
MultiImageOpticalFlowHelper<TFloat, VDim>
::ComputeAffineNCCMatchAndGradient(int level,
                                   LinearTransformType *tran,
                                   const SizeType &radius,
                                   FloatImageType *wrkMetric,
                                   FloatImageType *wrkMask,
                                   VectorImageType *wrkGradMetric,
                                   VectorImageType *wrkGradMask,
                                   VectorImageType *wrkPhi,
                                   LinearTransformType *grad)
{
  // Scale the weights by epsilon
  vnl_vector<float> wscaled(m_Weights.size());
  for (unsigned i = 0; i < wscaled.size(); i++)
    wscaled[i] = m_Weights[i];

  // Allocate a working image
  if(m_NCCWorkingImage.IsNull())
    m_NCCWorkingImage = MultiComponentImageType::New();

  // Set up the optical flow computation
  typedef DefaultMultiComponentImageMetricTraits<TFloat, VDim> TraitsType;
  typedef MultiComponentNCCImageMetric<TraitsType> MetricType;
  typename MetricType::Pointer metric = MetricType::New();

  // Is this the first time that this function is being called with this image?
  bool first_run =
      m_NCCWorkingImage->GetBufferedRegion() != m_FixedComposite[level]->GetBufferedRegion();

  // Check the radius against the size of the image
  SizeType radius_fix = AdjustNCCRadius(level, radius, first_run);

  metric->SetFixedImage(m_FixedComposite[level]);
  metric->SetMovingImage(m_MovingComposite[level]);
  metric->SetWeights(wscaled);
  metric->SetAffineTransform(tran);
  metric->SetComputeMovingDomainMask(false);
  metric->GetMetricOutput()->Graft(wrkMetric);
  metric->SetComputeGradient(grad != NULL);
  metric->SetRadius(radius_fix);
  metric->SetWorkingImage(m_NCCWorkingImage);
  metric->SetReuseWorkingImageFixedComponents(!first_run);
  metric->SetFixedMaskImage(m_GradientMaskComposite[level]);
  metric->SetJitterImage(m_JitterComposite[level]);
  metric->Update();

  // Process the results
  if(grad)
    {
    grad->SetMatrix(metric->GetAffineTransformGradient()->GetMatrix());
    grad->SetOffset(metric->GetAffineTransformGradient()->GetOffset());
    }

  return metric->GetMetricValue();


  /*
  // Scale the weights by epsilon
  vnl_vector<float> wscaled(m_Weights.size());
  for(int i = 0; i < wscaled.size(); i++)
    wscaled[i] = m_Weights[i];

  // Create a deformation field from the affine transform
  typedef LinearTransformToWarpFilter<
      MultiComponentImageType, VectorImageType, LinearTransformType> WarpFilter;
  typename WarpFilter::Pointer warp_source = WarpFilter::New();
  warp_source->SetInput(m_FixedComposite[level]);
  warp_source->SetTransform(tran);
  warp_source->GraftOutput(wrkPhi);

  // Allocate a working image
  if(m_NCCWorkingImage.IsNull())
    m_NCCWorkingImage = MultiComponentImageType::New();

  // Set up the optical flow computation
  typedef DefaultMultiComponentImageMetricTraits<TFloat, VDim> TraitsType;
  typedef MultiComponentNCCImageMetric<TraitsType> MetricType;
  typename MetricType::Pointer metric = MetricType::New();

  metric->SetWorkingImage(m_NCCWorkingImage);
  metric->SetRadius(radius);

  metric->SetFixedImage(m_FixedComposite[level]);
  metric->SetMovingImage(m_MovingComposite[level]);
  metric->SetWeights(wscaled);
  metric->SetDeformationField(warp_source->GetOutput());

  metric->GetMetricOutput()->Graft(wrkMetric);

  metric->SetComputeMovingDomainMask(true);
  metric->GetMovingDomainMaskOutput()->Graft(wrkMask);

  if(grad)
    {
    metric->SetComputeGradient(true);
    metric->GetGradientOutput()->Graft(wrkGradMetric);
    metric->GetMovingDomainMaskGradientOutput()->Graft(wrkGradMask);
    }

  // Use finite differences
  typedef MultiImageAffineMetricFilter<TraitsType> AffineMetricType;
  typename AffineMetricType::Pointer affine_metric = AffineMetricType::New();

  // Run the filter
  affine_metric->SetMetricImage(metric->GetMetricOutput());
  affine_metric->SetMovingDomainMaskImage(metric->GetMovingDomainMaskOutput());

  // TODO: only if gradient!
  if(grad)
    {
    affine_metric->SetComputeGradient(true);
    affine_metric->SetGradientImage(metric->GetGradientOutput());
    affine_metric->SetMovingDomainMaskGradientImage(metric->GetMovingDomainMaskGradientOutput());
    affine_metric->SetGradientScalingFactor(metric->GetGradientScalingFactor());
    }

  affine_metric->Update();

  // Process the results
  if(grad)
    {
    grad->SetMatrix(affine_metric->GetMetricGradient()->GetMatrix());
    grad->SetOffset(affine_metric->GetMetricGradient()->GetOffset());
    }

  // / *
  // LDDMMData<TFloat, VDim>::img_write(wrkMetric, "dump_metric.nii.gz");
  // LDDMMData<TFloat, VDim>::img_write(wrkMask, "dump_mask.nii.gz");
  // LDDMMData<TFloat, VDim>::vimg_write(wrkGradMetric, "dump_grad_metric.nii.gz");
  // LDDMMData<TFloat, VDim>::vimg_write(wrkGradMask, "dump_grad_mask.nii.gz");
  // exit(-1);
  // * /

  return affine_metric->GetMetricValue();
*/
}

template <class TFloat, unsigned int VDim>
void
MultiImageOpticalFlowHelper<TFloat, VDim>
::AffineToField(LinearTransformType *tran, VectorImageType *def)
{
  // TODO: convert this to a filter
  typedef itk::ImageLinearIteratorWithIndex<VectorImageType> IterBase;
  typedef IteratorExtender<IterBase> Iter;
  Iter it(def, def->GetBufferedRegion());
  it.SetDirection(0);

  for(; !it.IsAtEnd(); it.NextLine())
    {
    // Get the pointer to the begin of line
    VectorType *ptr = const_cast<VectorType *>(it.GetPosition());
    VectorType *ptr_end = ptr + def->GetBufferedRegion().GetSize(0);

    // Get the initial index
    typename LinearTransformType::InputPointType pt;
    for(int k = 0; k < VDim; k++)
      pt[k] = it.GetIndex()[k];

    for(; ptr < ptr_end; ++ptr, ++pt[0])
      {
      // Apply transform to the index. TODO: this is stupid, just use an offset
      typename LinearTransformType::OutputPointType pp = tran->TransformPoint(pt);
      for(int k = 0; k < VDim; k++)
        (*ptr)[k] = pp[k] - pt[k];
      }
    }
}


template <class TInputImage, class TOutputImage, class TFunctor>
class UnaryPositionBasedFunctorImageFilter : public itk::ImageToImageFilter<TInputImage, TOutputImage>
{
public:

  typedef UnaryPositionBasedFunctorImageFilter<TInputImage,TOutputImage,TFunctor> Self;
  typedef itk::ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef itk::SmartPointer<Self>                           Pointer;
  typedef itk::SmartPointer<const Self>                     ConstPointer;
  typedef typename Superclass::OutputImageRegionType         OutputImageRegionType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self)

  /** Run-time type information (and related methods) */
  itkTypeMacro( UnaryPositionBasedFunctorImageFilter, itk::ImageToImageFilter )

  /** Determine the image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int, TOutputImage::ImageDimension );

  void SetFunctor(const TFunctor &f) { this->m_Functor = f; }

protected:
  UnaryPositionBasedFunctorImageFilter() {}
  ~UnaryPositionBasedFunctorImageFilter() {}

  virtual void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                                    itk::ThreadIdType threadId) ITK_OVERRIDE 
  {
    typedef itk::ImageRegionConstIteratorWithIndex<TInputImage> InputIter;
    InputIter it_in(this->GetInput(), outputRegionForThread);

    typedef itk::ImageRegionIterator<TOutputImage> OutputIter;
    OutputIter it_out(this->GetOutput(), outputRegionForThread);

    for(; !it_out.IsAtEnd(); ++it_out, ++it_in)
      {
      it_out.Set(m_Functor(it_in.Get(), it_in.GetIndex()));
      }
  }

  TFunctor m_Functor;
};

template <class TWarpImage>
class VoxelToPhysicalWarpFunctor
{
public:
  typedef itk::ImageBase<TWarpImage::ImageDimension> ImageBaseType;
  typedef typename TWarpImage::PixelType VectorType;
  typedef itk::Index<TWarpImage::ImageDimension> IndexType;

  VectorType operator()(const VectorType &v, const IndexType &pos)
  {
    // Get the physical point for the tail of the arrow
    typedef itk::ContinuousIndex<double, TWarpImage::ImageDimension> CIType;
    typedef typename TWarpImage::PointType PtType;

    CIType ia, ib;
    PtType pa, pb;
    for(uint i = 0; i < TWarpImage::ImageDimension; i++)
      {
      ia[i] = pos[i];
      ib[i] = pos[i] + v[i];
      }

    m_Warp->TransformContinuousIndexToPhysicalPoint(ia, pa);
    m_MovingSpace->TransformContinuousIndexToPhysicalPoint(ib, pb);

    VectorType y;
    for(uint i = 0; i < TWarpImage::ImageDimension; i++)
      y[i] = pb[i] - pa[i];

    return y;
  }

  VoxelToPhysicalWarpFunctor(TWarpImage *warp, ImageBaseType *moving)
    : m_Warp(warp), m_MovingSpace(moving) {}

  VoxelToPhysicalWarpFunctor() {}

protected:

  TWarpImage *m_Warp;
  ImageBaseType *m_MovingSpace;
};


template <class TWarpImage>
class PhysicalToVoxelWarpFunctor
{
public:
  typedef itk::ImageBase<TWarpImage::ImageDimension> ImageBaseType;
  typedef typename TWarpImage::PixelType VectorType;
  typedef itk::Index<TWarpImage::ImageDimension> IndexType;

  VectorType operator()(const VectorType &v, const IndexType &pos)
  {
    // Get the voxel offset between the tip of the arrow and the input position
    // Get the physical point for the tail of the arrow
    typedef itk::ContinuousIndex<double, TWarpImage::ImageDimension> CIType;
    typedef typename TWarpImage::PointType PtType;

    CIType ia, ib;
    PtType pa, pb;

    // Get the base physical position
    for(uint i = 0; i < TWarpImage::ImageDimension; i++)
      ia[i] = pos[i];
    m_Warp->TransformContinuousIndexToPhysicalPoint(ia, pa);

    // Compute the tip physical position
    for(uint i = 0; i < TWarpImage::ImageDimension; i++)
      pb[i] = pa[i] + v[i];

    // Map the tip into continuous index
    m_MovingSpace->TransformPhysicalPointToContinuousIndex(pb, ib);

    VectorType y;
    for(uint i = 0; i < TWarpImage::ImageDimension; i++)
      y[i] = ib[i] - ia[i];

    return y;
  }

  PhysicalToVoxelWarpFunctor(TWarpImage *warp, ImageBaseType *moving)
    : m_Warp(warp), m_MovingSpace(moving) {}

  PhysicalToVoxelWarpFunctor() {}

protected:

  TWarpImage *m_Warp;
  ImageBaseType *m_MovingSpace;
};



/**
 * This functor is used to compress a warp before saving it. The input
 * to this functor is a voxel-space warp, and the output is a physical
 * space warp, with the precision of the voxel-space warp reduced to a
 * prescribed value. The functor will also cast the warp to desired
 * output type
 */
template <class TInputWarp, class TOutputWarp>
class CompressWarpFunctor
{
public:
  typedef VoxelToPhysicalWarpFunctor<TInputWarp> PhysFunctor;
  typedef typename PhysFunctor::ImageBaseType ImageBaseType;

  typedef typename TInputWarp::IndexType IndexType;
  typedef typename TInputWarp::PixelType InputVectorType;
  typedef typename TOutputWarp::PixelType OutputVectorType;

  CompressWarpFunctor() {}

  CompressWarpFunctor(TInputWarp *input, ImageBaseType *mov_space, double precision)
    : m_InputWarp(input), m_Precision(precision), m_ScaleFactor(1.0 / m_Precision),
      m_PhysFunctor(input, mov_space) {}

  OutputVectorType operator()(const InputVectorType &v, const IndexType &pos)
  {
    InputVectorType w;

    // Round to precision
    if(m_Precision > 0)
      {
      for(uint i = 0; i < TInputWarp::ImageDimension; i++)
        w[i] = std::floor(v[i] * m_ScaleFactor + 0.5) * m_Precision;
      }
    else
      {
      w = v;
      }

    // Map to physical space
    w = m_PhysFunctor(w, pos);

    // Cast to output type
    InputVectorType y;
    for(uint i = 0; i < TInputWarp::ImageDimension; i++)
      y[i] = static_cast<typename OutputVectorType::ValueType>(w[i]);

    return y;
  }

protected:
  TInputWarp *m_InputWarp;
  double m_Precision, m_ScaleFactor;
  PhysFunctor m_PhysFunctor;
};


template <class TFloat, unsigned int VDim>
void
MultiImageOpticalFlowHelper<TFloat, VDim>
::VoxelWarpToPhysicalWarp(VectorImageType *warp, ImageBaseType *moving_space, VectorImageType *result)
{
  typedef VoxelToPhysicalWarpFunctor<VectorImageType> Functor;
  typedef UnaryPositionBasedFunctorImageFilter<VectorImageType,VectorImageType,Functor> Filter;
  Functor functor(warp, moving_space);

  typename Filter::Pointer filter = Filter::New();
  filter->SetFunctor(functor);
  filter->SetInput(warp);
  filter->GraftOutput(result);
  filter->Update();
}

template <class TFloat, unsigned int VDim>
void
MultiImageOpticalFlowHelper<TFloat, VDim>
::PhysicalWarpToVoxelWarp(VectorImageType *warp, ImageBaseType *moving_space, VectorImageType *result)
{
  typedef PhysicalToVoxelWarpFunctor<VectorImageType> Functor;
  typedef UnaryPositionBasedFunctorImageFilter<VectorImageType,VectorImageType,Functor> Filter;
  Functor functor(warp, moving_space);

  typename Filter::Pointer filter = Filter::New();
  filter->SetFunctor(functor);
  filter->SetInput(warp);
  filter->GraftOutput(result);
  filter->Update();
}

template <class TFloat, unsigned int VDim>
void
MultiImageOpticalFlowHelper<TFloat, VDim>
::DownsampleWarp(VectorImageType *srcWarp, VectorImageType *trgWarp, int srcLevel, int trgLevel)
{
  typedef LDDMMData<TFloat, VDim> LDDMMType;

  // Get the factor by which to downsample
  int src_factor = m_PyramidFactors[srcLevel];
  int trg_factor = m_PyramidFactors[trgLevel];
  if(src_factor < trg_factor)
    {
    // Resample the warp - no smoothing
    LDDMMType::vimg_resample_identity(srcWarp, this->GetReferenceSpace(trgLevel), trgWarp);

    // Scale by the factor
    LDDMMType::vimg_scale_in_place(trgWarp, src_factor / trg_factor);
    }
  else if(src_factor == trg_factor)
    {
    LDDMMType::vimg_copy(srcWarp, trgWarp);
    }
  else
    {
    throw GreedyException("DownsampleWarp called for upsampling");
    }
}

template <class TFloat, unsigned int VDim>
void
MultiImageOpticalFlowHelper<TFloat, VDim>
::WriteCompressedWarpInPhysicalSpace(int level, VectorImageType *warp, const char *filename, double precision)
{
  WriteCompressedWarpInPhysicalSpace(warp, this->GetMovingReferenceSpace(level), filename, precision);
}

template <class TFloat, unsigned int VDim>
void
MultiImageOpticalFlowHelper<TFloat, VDim>
::WriteCompressedWarpInPhysicalSpace(VectorImageType *warp, ImageBaseType *moving_ref_space, const char *filename, double precision)
{
  // Define a _float_ output type, even if working with double precision (less space on disk)
  typedef itk::CovariantVector<float, VDim> OutputVectorType;
  typedef itk::Image<OutputVectorType, VDim> OutputWarpType;
  typedef CompressWarpFunctor<VectorImageType, OutputWarpType> Functor;

  typedef UnaryPositionBasedFunctorImageFilter<VectorImageType,OutputWarpType,Functor> Filter;
  Functor functor(warp, moving_ref_space, precision);

  typename Filter::Pointer filter = Filter::New();
  filter->SetFunctor(functor);
  filter->SetInput(warp);
  filter->Update();

  LDDMMData<float, VDim>::vimg_write(filter->GetOutput(), filename);
}

template <class TFloat, unsigned int VDim>
void 
MultiImageOpticalFlowHelper<TFloat, VDim>
::ComputeWarpSquareRoot(
  VectorImageType *warp, VectorImageType *out, VectorImageType *work, 
  FloatImageType *error_norm, double tol, int max_iter)
{
  typedef LDDMMData<TFloat, VDim> LDDMMType;

  // Use more convenient variables
  VectorImageType *u = warp, *v = out;

  // Initialize the iterate to zero
  v->FillBuffer(typename LDDMMType::Vec(0.0));

  // Perform iteration
  for(int i = 0; i < max_iter; i++)
    {
    // Min and max norm of the error at this iteration
    TFloat norm_max = tol, norm_min = 0.0; 

    // Perform a single iteration
    LDDMMType::interp_vimg(v, v, 1.0, work);   // work = v(v(x))
    LDDMMType::vimg_scale_in_place(work, -1.0); // work = -v(v(x))
    LDDMMType::vimg_add_scaled_in_place(work, v, -1.0); // work = -v(v) - v(v(x))
    LDDMMType::vimg_add_in_place(work, u); // work = u - v - v(v(x)) = f - g o g

    // At this point, 'work' stores the difference between actual warp and square of the
    // square root estimate, i.e., the estimation error
    if(error_norm)
      {
      LDDMMType::vimg_norm_min_max(work, error_norm, norm_min, norm_max);
      std::cout << " " << norm_max << " " << std::endl;
      }

    // Update v - which is being put into the output anyway
    LDDMMType::vimg_add_scaled_in_place(v, work, 0.5);

    // Check the maximum delta
    std::cout << "." << std::flush;

    // Break if the tolerance bound reached
    if(norm_max < tol)
      break;
    }
}


template <class TFloat, unsigned int VDim>
void 
MultiImageOpticalFlowHelper<TFloat, VDim>
::ComputeWarpRoot(VectorImageType *warp, VectorImageType *root, int exponent, TFloat tol, int max_iter)
{
  typedef LDDMMData<TFloat, VDim> LDDMMType;

  // If the exponent is zero, return the image itself
  if(exponent == 0)
    {
    LDDMMType::vimg_copy(warp, root);
    return;
    }

  // Create the current root and the next root
  VectorImagePointer u = VectorImageType::New();
  LDDMMType::alloc_vimg(u, warp);
  LDDMMType::vimg_copy(warp, u);

  // Create a working image
  VectorImagePointer work = VectorImageType::New();
  LDDMMType::alloc_vimg(work, warp);

  // If there is tolerance, create an error norm image
  FloatImagePointer error_norm;
  if(tol > 0.0)
    {
    error_norm = FloatImageType::New();
    LDDMMType::alloc_img(error_norm, warp);
    }
    
  // Compute the square root 'exponent' times
  for(int k = 0; k < exponent; k++)
    {
    // Square root of u goes into root
    ComputeWarpSquareRoot(u, root, work, error_norm, tol, max_iter);
    std::cout << std::endl;

    // Copy root into u
    LDDMMType::vimg_copy(root, u);
    }
}



template <class TFloat, unsigned int VDim>
void
MultiImageOpticalFlowHelper<TFloat, VDim>
::ComputeDeformationFieldInverse(
    VectorImageType *warp, VectorImageType *uInverse, int n_sqrt, bool verbose)
{
  typedef LDDMMData<TFloat, VDim> LDDMMType;

  // Create a copy of the forward warp
  VectorImagePointer uForward = VectorImageType::New();
  LDDMMType::alloc_vimg(uForward, warp);
  LDDMMType::vimg_copy(warp, uForward);

  // Create a working image 
  VectorImagePointer uWork = VectorImageType::New();
  LDDMMType::alloc_vimg(uWork, warp);

  // Take the desired square root of the input warp and place into uForward
  ComputeWarpRoot(warp, uForward, n_sqrt);

  // Clear uInverse
  uInverse->FillBuffer(itk::NumericTraits<typename LDDMMType::Vec>::ZeroValue());

  // At this point, uForward holds the small deformation
  // Try to compute the inverse of the current forward transformation
  for(int i = 0; i < 20; i++)
    {
    // We are using uPhys as temporary storage
    LDDMMType::interp_vimg(uForward, uInverse, 1.0, uWork);
    LDDMMType::vimg_scale_in_place(uWork, -1.0);

    // Compute the maximum change from last iteration
    LDDMMType::vimg_subtract_in_place(uInverse, uWork);

    // std::cout << "inverse iter " << i << " change " << norm_max << std::endl;
    LDDMMType::vimg_copy(uWork, uInverse);
    }

  // Compose the inverses
  for(int i = 0; i < n_sqrt; i++)
    {
    LDDMMType::interp_vimg(uInverse, uInverse, 1.0, uWork);
    LDDMMType::vimg_add_in_place(uInverse, uWork);
    }

  // If verbose, compute the maximum error
  if(verbose)
    {
    FloatImagePointer iNorm = FloatImageType::New();
    LDDMMType::alloc_img(iNorm, uWork);
    LDDMMType::interp_vimg(uInverse, uForward, 1.0, uWork);
    LDDMMType::vimg_add_in_place(uWork, uForward);
    TFloat norm_min, norm_max;
    LDDMMType::vimg_norm_min_max(uWork, iNorm, norm_min, norm_max);
    std::cout << "Warp inverse max residual: " << norm_max << std::endl;
    }
}

// Explicitly instantiate this class
template class MultiImageOpticalFlowHelper<float, 2>;
template class MultiImageOpticalFlowHelper<float, 3>;
template class MultiImageOpticalFlowHelper<float, 4>;
template class MultiImageOpticalFlowHelper<double, 2>;
template class MultiImageOpticalFlowHelper<double, 3>;
template class MultiImageOpticalFlowHelper<double, 4>;
