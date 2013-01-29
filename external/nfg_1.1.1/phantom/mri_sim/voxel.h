/*
 *  voxel.h
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
 
 
#ifndef VOXEL_H
#define VOXEL_H 
 
#include "phantom/shared/segment.h" 
#include "phantom/mri_sim/segment_register.h" 
#include "phantom/mri_sim/isotropic_region_register.h"
#include "phantom/mri_sim/subvoxel.h"


typedef struct _voxel {

	double size[3];
	double origin[3];
	int num_subvoxels[3];
	
	Subvoxel *subvoxels;
	double *orientations;

	Segment_register *segment_register;
	Isotropic_region_register *isotropic_region_register;

} Voxel;


void voxel_init(Voxel *voxel, double voxel_size[3], double voxel_origin[3], int num_subvoxels[3]);

void voxel_free(Voxel *voxel);

void register_segment(Voxel *voxel, Segment *segment);

void register_isotropic_region(Voxel *voxel, Isotropic_region *isotropic_region);

void plot_strand_orientations_in_voxel(Voxel *voxel);


#endif
