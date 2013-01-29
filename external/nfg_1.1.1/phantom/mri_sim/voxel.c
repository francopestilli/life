/*
 *  voxel.c
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
#include <math.h>

#include "phantom/shared/segment.h"
#include "phantom/mri_sim/voxel.h"
#include "phantom/mri_sim/segment_register.h" 
#include "phantom/mri_sim/isotropic_region_register.h"
#include "phantom/shared/shared.h"



void voxel_init(Voxel *voxel, double voxel_size[3], double voxel_origin[3], int num_subvoxels[3]) {

	voxel->size[X] = voxel_size[X];
	voxel->size[Y] = voxel_size[Y];
	voxel->size[Z] = voxel_size[Z];
	
	voxel->num_subvoxels[X] = num_subvoxels[X];
	voxel->num_subvoxels[Y] = num_subvoxels[Y];
	voxel->num_subvoxels[Z] = num_subvoxels[Z];	
			
	voxel->origin[X] = voxel_origin[X];
	voxel->origin[Y] = voxel_origin[Y];
	voxel->origin[Z] = voxel_origin[Z];			
			
	voxel->segment_register = NULL;

}


void voxel_free(Voxel *voxel) {

	Segment_register *segment_register;
	int x, y, z, offset;
	
	if (voxel->segment_register != NULL) {
		
		segment_register = voxel->segment_register;
		
		while (segment_register->next != NULL) {
			segment_register = segment_register->next;
			free(segment_register->prev);	
		}
		
		free(segment_register);

	}
	
	for (z = 0; z < voxel->num_subvoxels[Z]; z++) {
		for (y = 0; y < voxel->num_subvoxels[Y]; y++) {
			for (x = 0; x < voxel->num_subvoxels[X]; x++) {
			
				offset = z * voxel->num_subvoxels[Y] * voxel->num_subvoxels[X] + y * voxel->num_subvoxels[X] + x;
				subvoxel_free(&(voxel->subvoxels[offset]));
			}
			
		}
	}
	
	free(voxel->orientations);
	free(voxel->subvoxels);
	
}


void register_segment(Voxel *voxel, Segment *segment) {


	Segment_register *segment_register;
	
	if (voxel->segment_register == NULL) {
		voxel->segment_register = segment_register_alloc(segment);
		voxel->segment_register->prev = NULL;
		
	} else {
	
		segment_register = voxel->segment_register;
		while (segment_register->next != NULL) {
	
			segment_register = segment_register->next;
		}

		segment_register->next = segment_register_alloc(segment);
		segment_register->next->prev = segment_register;

	}

}


void register_isotropic_region(Voxel *voxel, Isotropic_region *isotropic_region) {


	Isotropic_region_register *isotropic_region_register;
	
	if (voxel->isotropic_region_register == NULL) {
		voxel->isotropic_region_register = isotropic_region_register_alloc(isotropic_region);
		voxel->isotropic_region_register->prev = NULL;
		
	} else {
	
		isotropic_region_register = voxel->isotropic_region_register;
		while (isotropic_region_register->next != NULL) {
	
			isotropic_region_register = isotropic_region_register->next;
		}

		isotropic_region_register->next = isotropic_region_register_alloc(isotropic_region);
		isotropic_region_register->next->prev = isotropic_region_register;

	}

}


void plot_strand_orientations_in_voxel(Voxel *voxel) {

	Segment_register *segment_register;
	Isotropic_region_register *isotropic_region_register;
	
	double subvoxel_size[3];
	int ubound_x, ubound_y, ubound_z, lbound_x, lbound_y, lbound_z;
	double subvoxel_centre[3], segment_orientation[3], disp_to_start_point[3], disp_to_end_point[3], in_plane[3], in_plane_length, dist_to_strand, dist_to_region, normal[3], fraction;
	int x, y, z, offset;
	
	Segment *segment;
	Isotropic_region *isotropic_region;

	voxel->orientations = (double*)malloc(sizeof(double) * voxel->num_subvoxels[X] * voxel->num_subvoxels[Y] * voxel->num_subvoxels[Z] * 3);
	voxel->subvoxels = (Subvoxel*)malloc(sizeof(Subvoxel) * voxel->num_subvoxels[X] * voxel->num_subvoxels[Y] * voxel->num_subvoxels[Z]); 


	for (z = 0; z < voxel->num_subvoxels[Z]; z++) {
		for (y = 0; y < voxel->num_subvoxels[Y]; y++) {
			for (x = 0; x < voxel->num_subvoxels[X]; x++) {

				offset = z * voxel->num_subvoxels[X] * voxel->num_subvoxels[Y] + y * voxel->num_subvoxels[X] + x;
				init_subvoxel(&(voxel->subvoxels[offset]), &(voxel->orientations[offset * 3]));

			}
		}
	}
	
	subvoxel_size[X] = voxel->size[X] / ((double)voxel->num_subvoxels[X]);
	subvoxel_size[Y] = voxel->size[Y] / ((double)voxel->num_subvoxels[Y]);
	subvoxel_size[Z] = voxel->size[Z] / ((double)voxel->num_subvoxels[Z]);
	
	isotropic_region_register = voxel->isotropic_region_register;
	
	while (isotropic_region_register != NULL) {
	
		isotropic_region = isotropic_region_register->isotropic_region;

		ubound_x = min_int( (int)ceil(  ( isotropic_region->pos[X] - voxel->origin[X] + isotropic_region->radius ) / subvoxel_size[X]), voxel->num_subvoxels[X]);
		ubound_y = min_int( (int)ceil(  ( isotropic_region->pos[Y] - voxel->origin[Y] + isotropic_region->radius ) / subvoxel_size[Y]), voxel->num_subvoxels[Y]);
		ubound_z = min_int( (int)ceil(  ( isotropic_region->pos[Z] - voxel->origin[Z] + isotropic_region->radius ) / subvoxel_size[Z]), voxel->num_subvoxels[Z]);		

		lbound_x = max_int( (int)floor(  ( isotropic_region->pos[X] - voxel->origin[X] - isotropic_region->radius ) / subvoxel_size[X]  ), 0 );
		lbound_y = max_int( (int)floor(  ( isotropic_region->pos[Y] - voxel->origin[Y] - isotropic_region->radius ) / subvoxel_size[Y]  ), 0 );
		lbound_z = max_int( (int)floor(  ( isotropic_region->pos[Z] - voxel->origin[Z] - isotropic_region->radius ) / subvoxel_size[Z]  ), 0 );

		for (z = lbound_z; z < ubound_z; z++) {
			for (y = lbound_y; y < ubound_y; y++) {
				for (x = lbound_x; x < ubound_x; x++) {

					subvoxel_centre[X] = ((double)(x) + 0.5) * subvoxel_size[X] + voxel->origin[X];
					subvoxel_centre[Y] = ((double)(y) + 0.5) * subvoxel_size[Y] + voxel->origin[Y];
					subvoxel_centre[Z] = ((double)(z) + 0.5) * subvoxel_size[Z] + voxel->origin[Z];
					
					dist_to_region = dist_between_points(isotropic_region->pos, subvoxel_centre);
					
					fraction = dist_to_region / isotropic_region->radius;
					
					if (fraction <= 1.0) { 
						add_isotropic_region(&(voxel->subvoxels[z * voxel->num_subvoxels[X] * voxel->num_subvoxels[Y] + y * voxel->num_subvoxels[X] + x]), fraction, isotropic_region);
					}
				}
			}
		}
		
		isotropic_region_register = isotropic_region_register->next;

					
	}

	segment_register = voxel->segment_register;	

	while (segment_register != NULL) {
	
	
		segment = segment_register->segment;
	
		segment_orientation[X] = segment->disp[X]/segment->length;
		segment_orientation[Y] = segment->disp[Y]/segment->length;
		segment_orientation[Z] = segment->disp[Z]/segment->length;
		
		
		ubound_x = min_int( (int)ceil(  ( max_dble(segment->start_point.pos[X], segment->end_point.pos[X]) - voxel->origin[X] + segment->strand->radius ) / subvoxel_size[X]), voxel->num_subvoxels[X]);
		ubound_y = min_int( (int)ceil(  ( max_dble(segment->start_point.pos[Y], segment->end_point.pos[Y]) - voxel->origin[Y] + segment->strand->radius ) / subvoxel_size[Y]), voxel->num_subvoxels[Y]);
		ubound_z = min_int( (int)ceil(  ( max_dble(segment->start_point.pos[Z], segment->end_point.pos[Z]) - voxel->origin[Z] + segment->strand->radius ) / subvoxel_size[Z]), voxel->num_subvoxels[Z]);		

		lbound_x = max_int( (int)floor(  ( min_dble(segment->start_point.pos[X], segment->end_point.pos[X]) - voxel->origin[X] - segment->strand->radius ) / subvoxel_size[X]  ), 0 );
		lbound_y = max_int( (int)floor(  ( min_dble(segment->start_point.pos[Y], segment->end_point.pos[Y]) - voxel->origin[Y] - segment->strand->radius ) / subvoxel_size[Y]  ), 0 );
		lbound_z = max_int( (int)floor(  ( min_dble(segment->start_point.pos[Z], segment->end_point.pos[Z]) - voxel->origin[Z] - segment->strand->radius ) / subvoxel_size[Z]  ), 0 );

		for (z = lbound_z; z < ubound_z; z++) {
			for (y = lbound_y; y < ubound_y; y++) {
				for (x = lbound_x; x < ubound_x; x++) {
			
					subvoxel_centre[X] = ((double)(x) + 0.5) * subvoxel_size[X] + voxel->origin[X];
					subvoxel_centre[Y] = ((double)(y) + 0.5) * subvoxel_size[Y] + voxel->origin[Y];
					subvoxel_centre[Z] = ((double)(z) + 0.5) * subvoxel_size[Z] + voxel->origin[Z];
					
					disp_to_start_point[X] = subvoxel_centre[X] - segment->start_point.pos[X];
					disp_to_start_point[Y] = subvoxel_centre[Y] - segment->start_point.pos[Y];					
					disp_to_start_point[Z] = subvoxel_centre[Z] - segment->start_point.pos[Z];	
					
					disp_to_end_point[X] = subvoxel_centre[X] - segment->end_point.pos[X];
					disp_to_end_point[Y] = subvoxel_centre[Y] - segment->end_point.pos[Y];					
					disp_to_end_point[Z] = subvoxel_centre[Z] - segment->end_point.pos[Z];
					
					
					if (dot_product(disp_to_start_point, segment->disp) < 0.0) {
						dist_to_strand = sqrt(dot_product(disp_to_start_point, disp_to_start_point));
						
					} else if (dot_product(disp_to_end_point, segment->disp) > 0.0) {
						dist_to_strand = sqrt(dot_product(disp_to_end_point, disp_to_end_point));
					} else {
					
						/* Cross the displacement between the subvoxel centre and the start of the strand segment to get an out of plane vector */ 
						normal[X] = disp_to_start_point[Y] * segment->disp[Z] - disp_to_start_point[Z] * segment->disp[Y];
						normal[Y] = disp_to_start_point[Z] * segment->disp[X] - disp_to_start_point[X] * segment->disp[Z];
						normal[Z] = disp_to_start_point[X] * segment->disp[Y] - disp_to_start_point[Y] * segment->disp[X];
						
						/* Get the vector perpendicular to the segment and in plane with it and the displacement to the subvoxel centre */ 
						in_plane[X] = segment->disp[Y] * normal[Z] - segment->disp[Z] * normal[Y];
						in_plane[Y] = segment->disp[Z] * normal[X] - segment->disp[X] * normal[Z];
						in_plane[Z] = segment->disp[X] * normal[Y] - segment->disp[Y] * normal[X];					
										
						in_plane_length = sqrt(dot_product(in_plane, in_plane));
						
						in_plane[X] = in_plane[X]/in_plane_length;
						in_plane[Y] = in_plane[Y]/in_plane_length;					
						in_plane[Z] = in_plane[Z]/in_plane_length;
						
						dist_to_strand = dot_product(disp_to_start_point, in_plane);
					}
					
					fraction = dist_to_strand/segment->strand->radius;
					
					if (fraction <= 1.0) { 
						add_orientation(&(voxel->subvoxels[z * voxel->num_subvoxels[X] * voxel->num_subvoxels[Y] + y * voxel->num_subvoxels[X] + x]), segment_orientation, fraction, segment);
					}
				}
			}
		}
		
		segment_register = segment_register->next;
	}

}



