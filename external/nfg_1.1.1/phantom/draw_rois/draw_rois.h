/*
 *  draw_rois.h
 *  Numerical Fibre Generator
 *
 *  Created by Tom Close on 25/06/08.
 *  Copyright 2008 Tom Close.
 *  Distributed under the GNU General Public Licence.
 *
 *
 *
 *  This file is part of 'Numerical Fibre Generator'.
 *
 *  'Numerical Fibre Generator' is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  'Numerical Fibre Generator' is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with 'Numerical Fibre Generator'.  If not, see <http://www.gnu.org/licenses/>
 *
 */

#include "phantom/shared/strand_collection.h"


#define NUM_VOXELS_DEFAULT 40
#define VOXEL_SIZE_DEFAULT 0.05
#define POINT_DEPTH_L_BOUND_DEFAULT 3
#define POINT_DEPTH_U_BOUND_DEFAULT 6
#define EXTRA_POINT_RADIUS_DEFAULT 0.05
#define VOXEL_RADIAL_L_BOUND_DEFAULT 0.0
#define VOXEL_RADIAL_U_BOUND_DEFAULT 1.0
#define SUBBUNDLE_MASK_DEFAULT 0
#define SAVE_COMBINED_MASK_DEFAULT 1
#define INCLUSION_RADIUS_DEFAULT 0.7
#define NULL_MASK_L_BOUND_DEFAULT 0.94
#define NULL_MASK_U_BOUND_DEFAULT 1.0

int* draw_rois(Strand_collection *c, int num_voxels, double voxel_size, double extra_point_radius, int point_depth_l_bound, int point_depth_u_bound, double voxel_radial_l_bound, double voxel_radial_u_bound, double inclusion_radius, double null_mask_l_bound, double null_mask_u_bound, int subbundle_mask, int *cumulative_num_strands, int *bundle_is_excluded);

void add_to_roi(int *masks, double *masks_closest, int num_voxels, double voxel_size, double voxel_radial_l_bound, double voxel_radial_u_bound,  Segment *segment, double extra_point_radius, int mask_value);

int get_start_roi_id(int bundle_id, int *cumulative_num_strands, int num_strands_in_bundle, int strand_i_in_bundle, int subbundle_mask);

int get_end_roi_id(int bundle_id, int *cumulative_num_strands, int num_strands_in_bundle, int strand_i_in_bundle, int subbundle_mask);