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
#ifndef __MultiImageRegistrationHelper_h
#define __MultiImageRegistrationHelper_h

#include "itkImageBase.h"
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkMatrixOffsetTransformBase.h"



/**
 * This class is used to perform mean square intensity difference type
 * registration with multiple images. The filter is designed for speed
 * of interpolation.
 */
template <class TFloat, unsigned int VDim>
class MultiImageOpticalFlowHelper
{
public:

  typedef itk::VectorImage<TFloat, VDim> MultiComponentImageType;
  typedef itk::Image<TFloat, VDim> FloatImageType;
  typedef itk::CovariantVector<TFloat, VDim> VectorType;
  typedef itk::Image<VectorType, VDim> VectorImageType;
  typedef itk::ImageBase<VDim> ImageBaseType;

  typedef typename MultiComponentImageType::Pointer MultiComponentImagePointer;
  typedef typename FloatImageType::Pointer FloatImagePointer;
  typedef typename VectorImageType::Pointer VectorImagePointer;

  typedef std::vector<int> PyramidFactorsType;
  typedef itk::Size<VDim> SizeType;
  typedef itk::CovariantVector<TFloat, VDim> Vec;

  typedef itk::MatrixOffsetTransformBase<TFloat, VDim, VDim> LinearTransformType;

  /** Set default (power of two) pyramid factors */
  void SetDefaultPyramidFactors(int n_levels);

  /** Set the pyramid factors - for multi-resolution (e.g., 8,4,2) */
  void SetPyramidFactors(const PyramidFactorsType &factors);

  /** 
   * Set whether the fixed images should be scaled down by the pyramid factors
   * when subsampling. This is needed for the Mahalanobis distance metric, but not for
   * any of the metrics that use image intensities 
   */
  void SetScaleFixedImageWithVoxelSize(bool onoff) { m_ScaleFixedImageWithVoxelSize = onoff; }

  /** Add a pair of multi-component images to the class - same weight for each component */
  void AddImagePair(MultiComponentImageType *fixed, MultiComponentImageType *moving, double weight);

  /** Set the gradient image mask */
  void SetGradientMask(FloatImageType *maskImage) { m_GradientMaskImage = maskImage; }

  /** Set the moving image mask */
  void SetMovingMask(FloatImageType *maskImage) { m_MovingMaskImage = maskImage; }

  /** Set jitter sigma - for jittering image samples in affine mode */
  void SetJitterSigma(double sigma);

  /** Compute the composite image - must be run before any sampling is done */
  void BuildCompositeImages(double noise_sigma_relative = 0.0);

  /**
   * Apply a dilation to the fixed gradient masks - this is used with the NCC metric. The initial
   * user-specified mask is transformed into a dilated mask, with values as follows:
   *   1.0 : voxel is inside the user-specified mask
   *   0.5 : voxel is within radius of the user-specified mask
   *   0.0 : voxel is outside of the user-specified mask
   *
   * NCC metrics exploit this mask format for faster processing - region where the mask is zero are
   * excluded from NCC computation and accumulation
   */
  void DilateCompositeGradientMasksForNCC(SizeType radius);

  /** Get the reference image for level k */
  ImageBaseType *GetReferenceSpace(int level);

  /** Get the reference image for level k */
  ImageBaseType *GetMovingReferenceSpace(int level);

  /** Get the gradient mask at a pyramid level */
  FloatImageType *GetGradientMask(int level) { return m_GradientMaskComposite[level]; }

  /** Get the moving mask at a pyramid level */
  FloatImageType *GetMovingMask(int level) { return m_MovingMaskComposite[level]; }

  /** Get the fixed image at a pyramid level */
  MultiComponentImageType *GetFixedComposite(int level) { return m_FixedComposite[level]; }

  /** Get the moving image at a pyramid level */
  MultiComponentImageType *GetMovingComposite(int level) { return m_MovingComposite[level]; }

  /** Get the smoothing factor for given level based on parameters */
  Vec GetSmoothingSigmasInPhysicalUnits(int level, double sigma, bool in_physical_units);

  /** Get the component weights in the composite */
  const std::vector<double> &GetWeights() const { return m_Weights; }

  /** Perform interpolation - compute [(I - J(Tx)) GradJ(Tx)] */
  vnl_vector<double> ComputeOpticalFlowField(
      int level, VectorImageType *def, FloatImageType *out_metric,
      VectorImageType *out_gradient, double result_scaling = 1.0);

  /** Perform interpolation - compute mutual information metric */
  vnl_vector<double> ComputeMIFlowField(
      int level, bool normalized_mutual_information,
      VectorImageType *def, FloatImageType *out_metric,
      VectorImageType *out_gradient, double result_scaling = 1.0);

  /** Compute the NCC metric without gradient */
  double ComputeNCCMetricImage(int level, VectorImageType *def, const SizeType &radius,
                              FloatImageType *out_metric, VectorImageType *out_gradient = NULL,
                               double result_scaling = 1.0);

  /** Compute the Mahalanobis metric with gradient */
  double ComputeMahalanobisMetricImage(int level, VectorImageType *def, 
                                       FloatImageType *out_metric, 
                                       VectorImageType *out_gradient = NULL);

  /** Compute affine similarity and gradient */
  double ComputeAffineMSDMatchAndGradient(int level, LinearTransformType *tran,
                                          FloatImageType *wrkMetric,
                                          FloatImageType *wrkMask,
                                          VectorImageType *wrkGradMetric,
                                          VectorImageType *wrkGradMask,
                                          VectorImageType *wrkPhi,
                                          LinearTransformType *grad = NULL);


  double ComputeAffineMIMatchAndGradient(int level, bool normalized_mutual_info,
                                         LinearTransformType *tran,
                                         FloatImageType *wrkMetric,
                                         FloatImageType *wrkMask,
                                         VectorImageType *wrkGradMetric,
                                         VectorImageType *wrkGradMask,
                                         VectorImageType *wrkPhi,
                                         LinearTransformType *grad = NULL);

  double ComputeAffineNCCMatchAndGradient(int level, LinearTransformType *tran,
                                          const SizeType &radius,
                                          FloatImageType *wrkMetric,
                                          FloatImageType *wrkMask,
                                          VectorImageType *wrkGradMetric,
                                          VectorImageType *wrkGradMask,
                                          VectorImageType *wrkPhi,
                                          LinearTransformType *grad = NULL);

  static void AffineToField(LinearTransformType *tran, VectorImageType *def);

  void DownsampleWarp(VectorImageType *srcWarp, VectorImageType *trgWarp, int srcLevel, int trgLevel);

  /** Convert a warp to physical space */
  static void VoxelWarpToPhysicalWarp(VectorImageType *warp, ImageBaseType *moving_space, VectorImageType *result);
  static void PhysicalWarpToVoxelWarp(VectorImageType *warp, ImageBaseType *moving_space, VectorImageType *result);

  /* 
   * Write a warp to a file. The warp must be in voxel space, not physical space 
   * this is the static version of this method
   */
  static void WriteCompressedWarpInPhysicalSpace(
    VectorImageType *warp, ImageBaseType *moving_ref_space, const char *filename, double precision);

  /** Write a warp to a file. The warp must be in voxel space, not physical space */
  void WriteCompressedWarpInPhysicalSpace(int level, VectorImageType *warp, const char *filename, double precision);

  /**
   * Invert a deformation field by first dividing it into small transformations using the
   * square root command, and then inverting the small transformations
   */
  static void ComputeDeformationFieldInverse(
    VectorImageType *warp, VectorImageType *result, int n_sqrt, bool verbose = false);

  /**
   * Compute the (2^k)-th root of a warp using an iterative scheme. For each
   * square root computation, the following iteration is used, where f = x + u
   * is the input warp, and g is the square root.
   */
  static void ComputeWarpRoot(
    VectorImageType *warp, VectorImageType *root, int exponent, TFloat tol = 0, int max_iter = 20);

  /**
   * Compute the square root of an input warp f = x + u(x) using an iterative scheme
   *
   *    g[0] = Id
   *    g[t+1] = g[t] + (f - g[t] o g[t]) / 2
   *
   * A working image of the same size as the input and output must be provided
   */
  static void ComputeWarpSquareRoot(
    VectorImageType *warp, VectorImageType *out, VectorImageType *work, 
    FloatImageType *error_norm = NULL, double tol = 0.0, int max_iter = 20);

  MultiImageOpticalFlowHelper() : 
    m_JitterSigma(0.0), m_ScaleFixedImageWithVoxelSize(false) {}

protected:

  // Pyramid factors
  PyramidFactorsType m_PyramidFactors;

  // Weights
  std::vector<double> m_Weights;

  // Vector of images
  typedef std::vector<typename MultiComponentImageType::Pointer> MultiCompImageSet;
  typedef std::vector<typename FloatImageType::Pointer> FloatImageSet;
  typedef std::vector<typename VectorImageType::Pointer> VectorImageSet;

  // Fixed and moving images
  MultiCompImageSet m_Fixed, m_Moving;

  // Composite image at each resolution level
  MultiCompImageSet m_FixedComposite, m_MovingComposite;

  // Working memory image for NCC computation
  typename MultiComponentImageType::Pointer m_NCCWorkingImage;

  // Gradient mask image - used to multiply the gradient
  typename FloatImageType::Pointer m_GradientMaskImage;

  // Moving mask image - used to reduce region where metric is computed
  typename FloatImageType::Pointer m_MovingMaskImage;

  // Mask composites
  FloatImageSet m_GradientMaskComposite, m_MovingMaskComposite;

  // Amount of jitter - for affine only
  double m_JitterSigma;

  // Jitter composite
  VectorImageSet m_JitterComposite;

  void PlaceIntoComposite(FloatImageType *src, MultiComponentImageType *target, int offset);
  void PlaceIntoComposite(VectorImageType *src, MultiComponentImageType *target, int offset);

  // Adjust NCC radius to be smaller than half image size
  SizeType AdjustNCCRadius(int level, const SizeType &radius, bool report_on_adjust);

  // Whether the fixed images should be scaled down by the pyramid factors
  // when subsampling. This is needed for the Mahalanobis distance metric, but not for
  // any of the metrics that use image intensities
  bool m_ScaleFixedImageWithVoxelSize;
};

#endif
