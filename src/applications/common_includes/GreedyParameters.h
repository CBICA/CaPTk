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
#ifndef GREEDYPARAMETERS_H
#define GREEDYPARAMETERS_H

#include <string>
#include <vector>
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"

struct ImagePairSpec
{
  std::string fixed;
  std::string moving;
  double weight;
};

struct SmoothingParameters
{
  double sigma;
  bool physical_units;
  SmoothingParameters(double s, bool pu) : sigma(s), physical_units(pu) {}
  SmoothingParameters() : sigma(0.0), physical_units(true) {}
};

struct RigidSearchSpec
{
  int iterations;
  double sigma_xyz;
  double sigma_angle;

  RigidSearchSpec() : iterations(0), sigma_xyz(0.0), sigma_angle(0.0) {}
};

struct InterpSpec
{
  enum InterpMode { LINEAR, NEAREST, LABELWISE };

  InterpMode mode;
  SmoothingParameters sigma;
  double outside_value;

  InterpSpec() : mode(LINEAR), sigma(0.5, false), outside_value(0.0) {}
};

struct ResliceSpec
{
  std::string moving;
  std::string output;
  InterpSpec interp;
};

struct ResliceMeshSpec
{
  std::string fixed;
  std::string output;
};

struct TransformSpec
{
  // Transform file
  std::string filename;

  // Optional exponent (-1 for inverse, 0.5 for square root)
  double exponent;
};

enum AffineInitMode
{
  VOX_IDENTITY = 0, // Identity mapping in voxel space
  RAS_IDENTITY,     // Identity mapping in physical space (i.e., use headers)
  RAS_FILENAME,     // User-specified matrix in physical space
  IMG_CENTERS,      // Match image centers, identity rotation in voxel space
  IMG_SIDE,         // Match image sides,
  MOMENTS_1,        // Match centers of mass,
  MOMENTS_2         // Match inertia tensors
};

struct GreedyResliceParameters
{
  // For reslice mode
  std::vector<ResliceSpec> images;
  std::vector<ResliceMeshSpec> meshes;

  // Reference image
  std::string ref_image;

  // Chain of transforms
  std::vector<TransformSpec> transforms;

  // Output warp
  std::string out_composed_warp;

  // Output jacobian
  std::string out_jacobian_image;
};

// Parameters for inverse warp command
struct GreedyInvertWarpParameters
{
  std::string in_warp, out_warp;
};

struct GreedyJacobianParameters
{
  std::string in_warp, out_det_jac;
};


// Parameters for inverse warp command
struct GreedyWarpRootParameters
{
  std::string in_warp, out_warp;
};


struct GreedyParameters
{
  enum MetricType { SSD = 0, NCC, MI, NMI, MAHALANOBIS };
  enum TimeStepMode { CONSTANT=0, SCALE, SCALEDOWN };
  enum Mode { GREEDY=0, AFFINE, BRUTE, RESLICE, INVERT_WARP, ROOT_WARP, JACOBIAN_WARP, MOMENTS };
  enum AffineDOF { DOF_RIGID=6, DOF_SIMILARITY=7, DOF_AFFINE=12 };

  std::vector<ImagePairSpec> inputs;
  std::string output;
  unsigned int dim;

  // Output for each iteration. This can be in the format "blah_%04d_%04d.mat" for
  // saving intermediate results into separate files. Or it can point to an object
  // in the GreedyAPI cache
  std::string output_intermediate;

  // Reslice parameters
  GreedyResliceParameters reslice_param;

  // Inversion parameters
  GreedyInvertWarpParameters invwarp_param;

  // Jacobian parameters
  GreedyJacobianParameters jacobian_param;

  // Root warp parameters
  GreedyWarpRootParameters warproot_param;

  // Registration mode
  Mode mode;

  bool flag_dump_moving, flag_debug_deriv, flag_powell;
  int dump_frequency, threads;
  double deriv_epsilon;

  double affine_jitter;

  double background;

  // Smoothing parameters
  SmoothingParameters sigma_pre, sigma_post;

  MetricType metric;
  TimeStepMode time_step_mode;

  // Iterations per level (i.e., 40x40x100)
  std::vector<double> epsilon_per_level;
  
  std::vector<int> iter_per_level;

  std::vector<int> metric_radius;

  std::vector<int> brute_search_radius;

  // List of transforms to apply to the moving image before registration
  std::vector<TransformSpec> moving_pre_transforms;

  // Initial affine transform
  AffineInitMode affine_init_mode;
  AffineDOF affine_dof;
  TransformSpec affine_init_transform;

  // Filename of initial warp
  std::string initial_warp;

  // Mask for gradient computation (fixed mask)
  std::string gradient_mask;

  // Mask for the moving image
  std::string moving_mask;

  // Inverse warp and root warp, for writing in deformable mode
  std::string inverse_warp, root_warp;
  int warp_exponent;

  // Precision for output warps
  double warp_precision;

  // Noise for NCC
  double ncc_noise_factor;

  // Debugging matrices
  bool flag_debug_aff_obj;

  // Rigid search
  RigidSearchSpec rigid_search;

  // Moments of inertia specification
  int moments_flip_determinant;
  int moments_order;
  bool flag_moments_id_covariance;

  // Stationary velocity (diffeomophic demons) mode
  bool flag_stationary_velocity_mode;

  // Whether the lie bracket is used in the y velocity update
  bool flag_stationary_velocity_mode_use_lie_bracket;

  // Floating point precision?
  bool flag_float_math;

  static void SetToDefaults(GreedyParameters &param);
};


#endif // GREEDYPARAMETERS_H
