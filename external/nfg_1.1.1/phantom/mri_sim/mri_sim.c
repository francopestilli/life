/*
 *  mri_sim.c
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



#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <string.h>

#include "phantom/shared/strand.h"
#include "phantom/shared/strand_collection.h"
#include "phantom/shared/segment.h"
#include "phantom/shared/shared.h"
#include "phantom/mri_sim/voxel.h"

#include "phantom/mri_sim/mri_sim.h"
#include "phantom/mri_sim/sim_voxel_intensities.h"


float* mri_sim(Strand_collection *c, double fa, double diffusivity, int num_voxels, double voxel_size, int num_subvoxels, int num_grad_directions, double *grad_directions, double *b_values, Strand_collection_stats *stats, int save_subvoxels, double **subvoxel_orientations) {

	int x, y, z, vox_offset, grad_i, strand_i, ubound_x, ubound_y, ubound_z, lbound_x, lbound_y, lbound_z, isotropic_region_i;
	int image_size;
	float *images;
	double fov;						/*fov is defined from the origin.  i.e. it defines a cube with (2 * fov) length sides.*/
	Voxel *voxels, *voxel;
	Segment *segment;
	Strand *strand;
	Isotropic_region *isotropic_region;
	int total_num_subvoxels;
	double *intensities;
	
	double *subvoxels, dummy;
	
	
	total_num_subvoxels = num_voxels * num_subvoxels;
	
	add_length_curv_stats(stats, c);
	
	if (save_subvoxels) {
	
		*subvoxel_orientations = (double*)calloc(sizeof(double),  total_num_subvoxels * total_num_subvoxels * total_num_subvoxels * 3);
		subvoxels = *subvoxel_orientations;
	} else {
		subvoxels = &dummy;
	}

	
	
	double voxel_origin[3], voxel_dim_sizes[3];
	int num_dim_subvoxels[3];

	fov = voxel_size * ((double)num_voxels) / 2.0;

	image_size = num_voxels * num_voxels * num_voxels;

	voxels = (Voxel*)malloc(sizeof(Voxel) * image_size);

	images = (float*)malloc(sizeof(float) * image_size * num_grad_directions);


	voxel_dim_sizes[X] = voxel_dim_sizes[Y] = voxel_dim_sizes[Z] = voxel_size;

	num_dim_subvoxels[X] = num_dim_subvoxels[Y] = num_dim_subvoxels[Z] = num_subvoxels; 


	for (z = 0; z < num_voxels; z++) {
		for (y = 0; y < num_voxels; y++) {
			for (x = 0; x < num_voxels; x++) {


				voxel_origin[X] = ((double)x) * voxel_size - fov;
				voxel_origin[Y] = ((double)y) * voxel_size - fov;
				voxel_origin[Z] = ((double)z) * voxel_size - fov;
				
				vox_offset = z * num_voxels * num_voxels + y * num_voxels + x;
				voxel_init(&(voxels[vox_offset]), voxel_dim_sizes, voxel_origin, num_dim_subvoxels);
		
			}
		}
	}


	for (isotropic_region_i = 0; isotropic_region_i < c->num_isotropic_regions; isotropic_region_i++) {
		isotropic_region = &(c->isotropic_regions[isotropic_region_i]);
		
		ubound_x = min_int( (int)ceil(  ( isotropic_region->pos[X] + isotropic_region->radius + fov) / voxel_size ), num_voxels);
		ubound_y = min_int( (int)ceil(  ( isotropic_region->pos[Y] + isotropic_region->radius + fov)  / voxel_size ), num_voxels);
		ubound_z = min_int( (int)ceil(  ( isotropic_region->pos[Z] + isotropic_region->radius + fov)  / voxel_size ), num_voxels);		

		lbound_x = max_int( (int)floor(  ( isotropic_region->pos[X] - isotropic_region->radius + fov)  / voxel_size ), 0);
		lbound_y = max_int( (int)floor(  ( isotropic_region->pos[Y] - isotropic_region->radius + fov)  / voxel_size ), 0);
		lbound_z = max_int( (int)floor(  ( isotropic_region->pos[Z] - isotropic_region->radius + fov)  / voxel_size ), 0);
		
		for (z = lbound_z; z < ubound_z; z++) {
			for (y = lbound_y; y < ubound_y; y++) {
				for (x = lbound_x; x < ubound_x; x++) {
		
					vox_offset = z * num_voxels * num_voxels + y * num_voxels + x;
					register_isotropic_region(&(voxels[vox_offset]), isotropic_region);

				}
			}
		}

	}



	for (strand_i = 0; strand_i < c->num_strands; strand_i++) {
	
		strand = &(c->strands[strand_i]);
		segment = strand->start_segment;
		
		while (segment != strand->post_segment) {
		
	
			ubound_x = min_int( (int)ceil(  ( max_dble(segment->start_point.pos[X], segment->end_point.pos[X]) + strand->radius + fov ) / voxel_size ), num_voxels);
			ubound_y = min_int( (int)ceil(  ( max_dble(segment->start_point.pos[Y], segment->end_point.pos[Y]) + strand->radius + fov ) / voxel_size ), num_voxels);
			ubound_z = min_int( (int)ceil(  ( max_dble(segment->start_point.pos[Z], segment->end_point.pos[Z]) + strand->radius + fov ) / voxel_size ), num_voxels);		

			lbound_x = max_int( (int)floor(  ( min_dble(segment->start_point.pos[X], segment->end_point.pos[X]) - strand->radius + fov ) / voxel_size ), 0);
			lbound_y = max_int( (int)floor(  ( min_dble(segment->start_point.pos[Y], segment->end_point.pos[Y]) - strand->radius + fov ) / voxel_size ), 0);
			lbound_z = max_int( (int)floor(  ( min_dble(segment->start_point.pos[Z], segment->end_point.pos[Z]) - strand->radius + fov ) / voxel_size ), 0);

			for (z = lbound_z; z < ubound_z; z++) {
				for (y = lbound_y; y < ubound_y; y++) {
					for (x = lbound_x; x < ubound_x; x++) {
			
						vox_offset = z * num_voxels * num_voxels + y * num_voxels + x;
						register_segment(&(voxels[vox_offset]), segment);

					}
				}
			}

		
			segment = segment->next_segment;
		}
	
	
	}

	int sub_x, sub_y, sub_z;



	for (z = 0; z < num_voxels; z++) {
		for (y = 0; y < num_voxels; y++) {
			for (x = 0; x < num_voxels; x++) {

				vox_offset = z * num_voxels * num_voxels + y * num_voxels + x;
				voxel = &(voxels[vox_offset]);
	
				plot_strand_orientations_in_voxel(voxel);
				

				intensities = sim_voxel_intensities(voxel, num_grad_directions, grad_directions, b_values, fa, diffusivity); 
				
				for (grad_i = 0; grad_i < num_grad_directions; grad_i++) {
					images[grad_i * image_size + vox_offset] = (float)(intensities[grad_i]);				
				}

	
	
				free(intensities);
	
				if (save_subvoxels) {
					
					for (sub_z = 0; sub_z < num_subvoxels; sub_z++) {				
						for (sub_y = 0; sub_y < num_subvoxels; sub_y++) {				
							for (sub_x = 0; sub_x < num_subvoxels; sub_x++) {
								subvoxels[( (z * num_subvoxels + sub_z) * total_num_subvoxels * total_num_subvoxels + (y * num_subvoxels + sub_y) * total_num_subvoxels + (x * num_subvoxels + sub_x)) * 3 + X] = voxel->orientations[(sub_z * num_subvoxels * num_subvoxels + sub_y * num_subvoxels + sub_x) * 3 + X];
								subvoxels[( (z * num_subvoxels + sub_z) * total_num_subvoxels * total_num_subvoxels + (y * num_subvoxels + sub_y) * total_num_subvoxels + (x * num_subvoxels + sub_x)) * 3 + Y] = voxel->orientations[(sub_z * num_subvoxels * num_subvoxels + sub_y * num_subvoxels + sub_x) * 3 + Y];
								subvoxels[( (z * num_subvoxels + sub_z) * total_num_subvoxels * total_num_subvoxels + (y * num_subvoxels + sub_y) * total_num_subvoxels + (x * num_subvoxels + sub_x)) * 3 + Z] = voxel->orientations[(sub_z * num_subvoxels * num_subvoxels + sub_y * num_subvoxels + sub_x) * 3 + Z];
							}
						}
					}
				}

				
				
				add_overlap_stats(stats, voxel, c->sphere_r);
				
		
				voxel_free(voxel);
		
			}
		}
		
		printf("Simulated slice %d intensities\n", z);
		fflush(stdout);

	}
	
	printf("\n");
	fflush(stdout);

	
	free(voxels);
	
	return images;
}


