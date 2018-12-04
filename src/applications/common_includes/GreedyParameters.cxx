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
#include "GreedyParameters.h"



void
GreedyParameters
::SetToDefaults(GreedyParameters &param)
{
	param.dim = 2;
	param.mode = GreedyParameters::GREEDY;

	param.threads = 0;
	param.metric = GreedyParameters::NMI;
	param.affine_init_mode = IMG_CENTERS;
	param.affine_dof = GreedyParameters::DOF_AFFINE;
	param.time_step_mode = GreedyParameters::SCALE;
	param.deriv_epsilon = 1e-4;
	param.warp_exponent = 6;
	param.warp_precision = 0.1;
	param.ncc_noise_factor = 0.001;
	
	param.dump_frequency = 1;
	param.epsilon_per_level = std::vector<double>(1, 1.0);
	param.sigma_pre.sigma = sqrt(3.0);
	param.sigma_pre.physical_units = false;
	param.sigma_post.sigma = sqrt(0.5);
	param.affine_jitter = 0.5;
	param.flag_powell = false;
	param.flag_dump_moving = false;
	param.flag_debug_deriv = false;
	param.flag_debug_aff_obj = false;
	param.flag_float_math = false;
	param.flag_stationary_velocity_mode = false;
	param.flag_stationary_velocity_mode_use_lie_bracket = false;
	param.sigma_post.physical_units = false;
	param.background = 0.0;

	// reslice mode parameters
	InterpSpec interp_current;

	param.iter_per_level.push_back(100);
	param.iter_per_level.push_back(100);

	/*param.iter_per_level.push_back(100);
	param.iter_per_level.push_back(50);
	param.iter_per_level.push_back(10);*/

	// Moments of inertia parameters
	param.moments_flip_determinant = 0;
	param.flag_moments_id_covariance = false;
	param.moments_order = 1;
}
