/*
 *  mri_sim.h
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
#include "phantom/mri_sim/strand_collection_stats.h"


#define FA_DEFAULT 0.8

#define DIFFUSIVITY_DEFAULT 0.0009
#define NUM_VOXELS_DEFAULT 20
#define VOXEL_SIZE_DEFAULT 0.1
#define SUBVOXEL_SIZE_DEFAULT 0.0025
#define SPHERE_R_DEFAULT 1.0
#define SAVE_SUBVOXELS_DEFAULT 0
#define SAVE_EXT_SEGMENT_STATS_DEFAULT 0
#define SAVE_EXT_FILL_RADIAL_STATS_DEFAULT 0
#define OUTPUT_FORMAT_DEFAULT "analyze"

#define GRAD_DIR_BLOCK_SIZE 100
#define SUBVOXEL_ROUND_DOWN_WARNING_THRESHOLD 0.00001


float* mri_sim(Strand_collection *c, double fa, double diffusivity, int num_voxels, double voxel_size, int num_subvoxels, int num_grad_directions, double *grad_directions, double *b_values, Strand_collection_stats *stats, int save_subvoxels, double **subvoxel_orientations) ;
