/*
 *  draw_roi.c
 *  Numerical Fibre Generator
 *
 *  Created by Tom Close on 11/12/08.
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

#include "phantom/draw_rois/draw_rois.h"


int* draw_rois(Strand_collection *c, int num_voxels, double voxel_size, double extra_point_radius, int point_depth_l_bound, int point_depth_u_bound, double voxel_radial_l_bound, double voxel_radial_u_bound, double inclusion_radius, double null_mask_l_bound, double null_mask_u_bound, int subbundle_mask, int *cumulative_num_strands, int *bundle_is_excluded) {

	Strand *strand;
	Segment *segment, *forward_segment, *backward_segment;
	int *combined_masks, num_elems, strand_i, depth_i, max_bundle_i, x, y, z, offset;
	double *masks_closest, voxel_centre[3], fov, radial_dist;
	int *strand_i_in_bundle;
	int start_roi_id, end_roi_id;
	int in_inclusion_sphere;
	
	
	max_bundle_i = c->bundles[c->num_bundles-1].bundle_i;
	
	strand_i_in_bundle = (int*)calloc(max_bundle_i+1, sizeof(int)); //Used to index each strand within their bundle.
	
	
	//For each bundle, determine whether any of them pass through the inclusion sphere. 
	for (strand_i = 0; strand_i < c->num_strands; strand_i++) {
	
		strand = &(c->strands[strand_i]);
		
		segment = strand->start_segment;
	
		in_inclusion_sphere = 0;
				
		while (segment != NULL) {
		
			if (vector_norm(segment->start_point.pos) <= inclusion_radius) {
				in_inclusion_sphere = 1;
			}
			segment = segment->next_segment;
		}
		
		if (!in_inclusion_sphere) {
			bundle_is_excluded[strand->bundle_i] = 1;
		}
		
	}
	
	
	num_elems = num_voxels * num_voxels * num_voxels;
	
	combined_masks = (int*)malloc(sizeof(int) * num_elems);
	masks_closest = (double*)malloc(sizeof(double) * num_elems);

	fov = ((double)num_voxels) * voxel_size / 2.0;

	for (z = 0; z < num_voxels; z++) {
		for (y = 0; y < num_voxels; y++) {	
			for (x = 0; x < num_voxels; x++) {
				
				offset = z * num_voxels * num_voxels + y * num_voxels + x;
				
				voxel_centre[X] = ((double)x + 0.5) * voxel_size - fov;
				voxel_centre[Y] = ((double)y + 0.5) * voxel_size - fov;
				voxel_centre[Z] = ((double)z + 0.5) * voxel_size - fov;
				
				radial_dist = vector_norm(voxel_centre);
				
				if ( (radial_dist >= null_mask_l_bound) & (radial_dist <= null_mask_u_bound) ) {
					combined_masks[offset] = 1;
				} else {
					combined_masks[offset] = 0;
				}
				
				masks_closest[offset] = 1.0;
			}
		}
		
	}


	for (strand_i = 0; strand_i < c->num_strands; strand_i++) {
	
		strand = &(c->strands[strand_i]);
		
		if (!bundle_is_excluded[strand->bundle_i]) {
			
			
//			if (subbundle_mask) {
//				start_roi_id = cumulative_num_strands[strand->bundle_i] * 2 + strand_i_in_bundle[strand->bundle_i] + 2;
//				end_roi_id = cumulative_num_strands[strand->bundle_i] * 2 + strand->bundle->num_strands + strand_i_in_bundle[strand->bundle_i] + 2;
//			} else {
//				start_roi_id = strand->bundle_i * 2 + 2;
//				end_roi_id = strand->bundle_i * 2 + 3;
//			}

			start_roi_id = get_start_roi_id(strand->bundle_i, cumulative_num_strands, strand->bundle->num_strands, strand_i_in_bundle[strand->bundle_i], subbundle_mask);
			end_roi_id = get_end_roi_id(strand->bundle_i, cumulative_num_strands, strand->bundle->num_strands, strand_i_in_bundle[strand->bundle_i], subbundle_mask);


			strand_i_in_bundle[strand->bundle_i]++;
			
		
			if ( (2 * point_depth_l_bound) > strand->num_control_points) {
				printf("Warning! Point depth %d exceeds half the number of control points in strand %d. Excluding strand from ROI...\n", point_depth_l_bound, strand_i);
			} else {
				
				
				forward_segment = strand->start_segment;
				backward_segment = strand->end_segment;
				
				for (depth_i = 0; depth_i <= point_depth_u_bound; depth_i++) {
					forward_segment = forward_segment->next_segment;
					backward_segment = backward_segment->prev_segment;
					
					if (forward_segment->segment_i >= backward_segment->segment_i)
						break;
					
					if (depth_i >= point_depth_l_bound) {

						/*Add point to start ROI*/
						add_to_roi(combined_masks, masks_closest, num_voxels, voxel_size, voxel_radial_l_bound, voxel_radial_u_bound, forward_segment, extra_point_radius, start_roi_id);

						/*Add point to end ROI*/				
						add_to_roi(combined_masks, masks_closest, num_voxels, voxel_size, voxel_radial_l_bound, voxel_radial_u_bound, backward_segment, extra_point_radius, end_roi_id);
					}	
				}
			}
		}
	}
		
	free(masks_closest);
	free(strand_i_in_bundle);
	return combined_masks;

}

//Very similar to plot_orientations in voxel
void add_to_roi(int *masks, double *masks_closest, int num_voxels, double voxel_size, double voxel_radial_l_bound, double voxel_radial_u_bound,  Segment *segment, double extra_point_radius, int mask_value) {
	
	double l_bound[3], u_bound[3], voxel_centre[3], fraction, normal[3], disp_to_start_point[3], disp_to_end_point[3], dist_to_strand, in_plane[3], in_plane_length, fov, voxel_radial_dist;
	int x, y, z, offset;
	double inflated_radius;
	
	fov = ((double)num_voxels) * voxel_size / 2.0;
	
	inflated_radius = segment->strand->radius * extra_point_radius;
	
	l_bound[X] = max_int((int)floor((min_dble(segment->start_point.pos[X], segment->end_point.pos[X]) - inflated_radius + fov)/voxel_size), 0);
	l_bound[Y] = max_int((int)floor((min_dble(segment->start_point.pos[Y], segment->end_point.pos[Y]) - inflated_radius + fov)/voxel_size), 0);
	l_bound[Z] = max_int((int)floor((min_dble(segment->start_point.pos[Z], segment->end_point.pos[Z]) - inflated_radius + fov)/voxel_size), 0);
	
	u_bound[X] = min_int((int)ceil((max_dble(segment->start_point.pos[X], segment->end_point.pos[X]) + inflated_radius + fov)/voxel_size), num_voxels);
	u_bound[Y] = min_int((int)ceil((max_dble(segment->start_point.pos[Y], segment->end_point.pos[Y]) + inflated_radius + fov)/voxel_size), num_voxels);
	u_bound[Z] = min_int((int)ceil((max_dble(segment->start_point.pos[Z], segment->end_point.pos[Z]) + inflated_radius + fov)/voxel_size), num_voxels);
	
	
	for (z = l_bound[Z]; z < u_bound[Z]; z++) {
		for (y = l_bound[Y]; y < u_bound[Y]; y++) {
			for (x = l_bound[X]; x < u_bound[X]; x++) {
		
				offset = z * num_voxels * num_voxels + y * num_voxels + x;
		
				voxel_centre[X] = ((double)x + 0.5) * voxel_size - fov;
				voxel_centre[Y] = ((double)y + 0.5) * voxel_size - fov;
				voxel_centre[Z] = ((double)z + 0.5) * voxel_size - fov;
				
				
				voxel_radial_dist = vector_norm(voxel_centre);
				
				if (voxel_radial_dist >= voxel_radial_l_bound && voxel_radial_dist <= voxel_radial_u_bound) {
				
					disp_to_start_point[X] = voxel_centre[X] - segment->start_point.pos[X];
					disp_to_start_point[Y] = voxel_centre[Y] - segment->start_point.pos[Y];					
					disp_to_start_point[Z] = voxel_centre[Z] - segment->start_point.pos[Z];	
					
					disp_to_end_point[X] = voxel_centre[X] - segment->end_point.pos[X];
					disp_to_end_point[Y] = voxel_centre[Y] - segment->end_point.pos[Y];					
					disp_to_end_point[Z] = voxel_centre[Z] - segment->end_point.pos[Z];
					
					
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
					
					fraction = dist_to_strand/inflated_radius;
																
														
					if (fraction <= masks_closest[offset]) {
						masks[offset] = mask_value;
						masks_closest[offset] = fraction;
					}
					
				}
				
			}
		}
	}


}



int get_start_roi_id(int bundle_i, int *cumulative_num_strands, int num_strands_in_bundle, int strand_i_in_bundle, int subbundle_mask) {

	int start_roi_id;
	
	if (subbundle_mask) {
		start_roi_id = cumulative_num_strands[bundle_i] * 2 + strand_i_in_bundle + 2;

	} else {
		start_roi_id = bundle_i * 2 + 2;

	}
	
	return start_roi_id;

}

int get_end_roi_id(int bundle_i, int *cumulative_num_strands, int num_strands_in_bundle, int strand_i_in_bundle, int subbundle_mask) {

	int end_roi_id;

	if (subbundle_mask) {
		end_roi_id = cumulative_num_strands[bundle_i] * 2 + num_strands_in_bundle + strand_i_in_bundle + 2;
	} else {
		end_roi_id = bundle_i * 2 + 3;
	}

	return end_roi_id;

}


///* Maps the given bundle index ('bundle_i') to a ROI index, which is the index of the given bundle_i in the 'c->bundles_in_collection' array, multiplied by two (even numbers are start ROIs and odd are end ROIs).*/
//int get_roi_i(Strand_collection *c, int bundle_i) {
//	
//	int roi_i, *bundle_i_address, bundle_i_offset;
//	
//	/* Finds the address of the given bundle index in the 'c->bundles_in_collection' array */
//	bundle_i_address = bsearch(&(bundle_i), c->bundles_in_collection, c->num_bundles, sizeof(int), compare_int);
//	
//	if (bundle_i_address == 0) {
//		return -1; 
//	} 
//	
//	/* The offset of the bundle_i_address from the start of the 'c->bundles_in_collection' array */
//	bundle_i_offset = ((int)bundle_i_address) - ((int)c->bundles_in_collection);
//	
//	roi_i = bundle_i_offset * 2 / sizeof(int);
//	
//	return roi_i;
//}
