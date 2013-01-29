/*
 *  subdiv.c
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
#include "phantom/shared/shared.h"

#include "phantom/resample/resample.h"



void resample_collection(Strand_collection *old_c, Strand_collection *resample_c, double resample_length, double double_back_angle_threshold, double forward_angle_threshold) {

	int strand_i;
	double combined_segment_length;
	int ubound_num_control_points;
	Segment *segment;


	combined_segment_length = 0.0;		//To calculate an upper bound on the total number of points that could be generated for these resampled collection.
	
	for (strand_i=0; strand_i < old_c->num_strands; strand_i++) {

		segment = old_c->strands[strand_i].start_segment;
		
		while (segment != NULL) {
			combined_segment_length += segment->length;
			segment = segment->next_segment;
		}

	}

	
	ubound_num_control_points = (int)ceil(combined_segment_length / resample_length) + 10 * old_c->num_strands; /* Add an extra 10 per strand just to be safe */
	
	collection_alloc(resample_c, old_c->num_strands, ubound_num_control_points, old_c->num_isotropic_regions);
	copy_isotropic_regions(resample_c, old_c->isotropic_regions, old_c->num_isotropic_regions);

	resample_c->num_control_points = 0;

		
	for (resample_c->num_strands = 0; resample_c->num_strands < old_c->num_strands; resample_c->num_strands++) {
	
		resample_strand(&(old_c->strands[resample_c->num_strands]), resample_c, resample_length, double_back_angle_threshold, forward_angle_threshold);
	}	

	resample_c->control_points = (double*)realloc(resample_c->control_points, sizeof(double) * 3 * resample_c->num_control_points);
	resample_c->control_points_grad = (double*)realloc(resample_c->control_points_grad, sizeof(double) * resample_c->num_control_points);
	resample_c->segments = (Segment*)realloc(resample_c->segments, sizeof(Segment) * (resample_c->num_control_points + 3 * resample_c->num_strands)); 

}



void resample_strand(Strand *old_strand, Strand_collection *resample_c, double resample_length, double double_back_angle_threshold, double forward_angle_threshold) {


	Segment *old_segment, *prev_old_segment;
	double a, b, c, k1, k2, k_prev, k;
	double *resample_point, *prev_resample_point, *start_point, *end_point, post_segment_scalar;
	int k1_lies_on_seg, k2_lies_on_seg;
	
	int start_of_strand_control_point_i;
	
	double double_back_angle;
	
	double debug[3];
	double debug_angle, debug_length;
	int debug_count;
	
	
	resample_c->num_strand_control_points[resample_c->num_strands] = 0;
	start_of_strand_control_point_i = resample_c->num_control_points;	

	old_segment = old_strand->pre_segment;

	resample_point = &(resample_c->pre_points[resample_c->num_strands * 3]);
	
	resample_point[X] = old_segment->start_point.pos[X];
	resample_point[Y] = old_segment->start_point.pos[Y];
	resample_point[Z] = old_segment->start_point.pos[Z];

	prev_old_segment = old_segment;
	old_segment = old_strand->start_segment;	
	
	resample_point = &(resample_c->start_points[resample_c->num_strands * 3]);
	
	resample_point[X] = old_segment->start_point.pos[X];
	resample_point[Y] = old_segment->start_point.pos[Y];
	resample_point[Z] = old_segment->start_point.pos[Z];
	
	prev_resample_point = resample_point;
	resample_point = &(resample_c->control_points[resample_c->num_control_points * 3]);
	
	k_prev = 0.0; 

	debug_count = 0;

	while (old_segment != old_strand->post_segment) {
	
		start_point = old_segment->start_point.pos;
		end_point = old_segment->end_point.pos;
	
		while (1) {
	
			/* Solve for the points that lie a distance 'resample_length' away from the previous resampled control point and on the line passing through the current segment in the old strand */
			a = dot_product(end_point, end_point) - 2.0 * dot_product(start_point, end_point) + dot_product(start_point, start_point);
			
			b = 2.0 * dot_product(start_point, end_point) - 2.0 * dot_product(end_point, prev_resample_point) - 2.0 * dot_product(start_point, start_point) + 2.0 * dot_product(start_point, prev_resample_point);
			
			c = dot_product(start_point, start_point) - 2.0 * dot_product(start_point, prev_resample_point) + dot_product(prev_resample_point, prev_resample_point) - resample_length * resample_length;
			
			
			/* k1 and k2 are the scalar displacements along the line-extension of the current old segment from its start point*/
			k1 = (-b + sqrt( b * b - 4.0 * a * c)) / (2.0 * a);
			
			k2 = (-b - sqrt( b * b - 4.0 * a * c)) / (2.0 * a);


			/* Check to see whether the points corresponding to k1 and k2 acutally lie on the segment and are further along than the previous resampled point*/		
			k1_lies_on_seg = (k1 > 0.0) & (k1 <= 1.0);
			
			k2_lies_on_seg = (k2 > 0.0) & (k2 <= 1.0);
		
 
			if ( k1_lies_on_seg & k2_lies_on_seg )
				k = max_dble(k1, k2);
						
			else if ( k1_lies_on_seg & (k1 > k_prev) )
				k = k1;
				
			else if (k2_lies_on_seg & (k2 > k_prev))
				k = k2;
				
			/* Else move on to the next segment in the old strand */	
			else {
				k_prev = 0.0;
				break;
			}
			
			
			resample_point[X] = start_point[X] + k * old_segment->disp[X];
			resample_point[Y] = start_point[Y] + k * old_segment->disp[Y];
			resample_point[Z] = start_point[Z] + k * old_segment->disp[Z];		
												
			resample_c->num_control_points++;
			resample_c->num_strand_control_points[resample_c->num_strands]++;

			k_prev = k;			
			prev_resample_point = resample_point;
			resample_point = &(resample_c->control_points[resample_c->num_control_points * 3]);

		
		}
		
		prev_old_segment = old_segment;
		old_segment = old_segment->next_segment;

		
		/* First check that the next segment doesn't double back on its self, which can be caused by the subdivision method */
		double_back_angle = angle_between_points(prev_resample_point, prev_old_segment->end_point.pos, old_segment->end_point.pos);
		
		if (double_back_angle > double_back_angle_threshold) {

			 while (double_back_angle > forward_angle_threshold) {

				if (old_segment == old_strand->post_segment) {
					break;
				}

				old_segment = old_segment->next_segment;

				double_back_angle = angle_between_points(prev_resample_point, prev_old_segment->end_point.pos, old_segment->start_point.pos);
				
				
			 }
			
			
		}
			


	}
	

	
	/* NB: old_segment is now the post_segment of the old strand */
	
	/* Part of a Slight hack - The post segment of the old strand is rescaled to be the resample length so the last resample is guaranteed to fall on it */
		
	post_segment_scalar = 2.0 * resample_length / old_segment->length;
	
	old_segment->disp[X] *= post_segment_scalar; 
	old_segment->disp[Y] *= post_segment_scalar;
	old_segment->disp[Z] *= post_segment_scalar;
	
	old_segment->end_point.pos[X] = old_segment->start_point.pos[X] + old_segment->disp[X];
	old_segment->end_point.pos[X] = old_segment->start_point.pos[X] + old_segment->disp[X];					
	old_segment->end_point.pos[X] = old_segment->start_point.pos[X] + old_segment->disp[X];
	
	old_segment->length = resample_length;
	
	start_point = old_segment->start_point.pos;
	end_point = old_segment->end_point.pos;
	
	resample_point = &(resample_c->end_points[resample_c->num_strands * 3]);
	
	a = dot_product(end_point, end_point) - 2.0 * dot_product(start_point, end_point) + dot_product(start_point, start_point);
	
	b = 2.0 * dot_product(start_point, end_point) - 2.0 * dot_product(end_point, prev_resample_point) - 2.0 * dot_product(start_point, start_point) + 2.0 * dot_product(start_point, prev_resample_point);
	
	c = dot_product(start_point, start_point) - 2.0 * dot_product(start_point, prev_resample_point) + dot_product(prev_resample_point, prev_resample_point) - resample_length * resample_length;
	
	
	k1 = (-b + sqrt( b * b - 4 * a * c)) / (2 * a);
	
	k2 = (-b - sqrt( b * b - 4 * a * c)) / (2 * a);
	

	
	debug[X] = start_point[X] - prev_resample_point[X];
	debug[Y] = start_point[Y] - prev_resample_point[Y];
	debug[Z] = start_point[Z] - prev_resample_point[Z];	
		

	debug_length = vector_norm(debug);
	debug_angle = angle_between_points(prev_old_segment->start_point.pos, old_segment->start_point.pos, old_segment->end_point.pos);
	
	/* Where the hack comes in.  For some reason neither k1 or k2 fall on the post segment set it to 1 */
	if ( ((k1 < 0) | (k1 > 1))  & ((k2 < 0) | (k2 > 1))) {
		k = 1.0;
	} else {
		k = max_dble(k1, k2);
	}
	
	
	resample_point[X] = start_point[X] + k * old_segment->disp[X];
	resample_point[Y] = start_point[Y] + k * old_segment->disp[Y];
	resample_point[Z] = start_point[Z] + k * old_segment->disp[Z];	
	
	prev_resample_point = resample_point;
	resample_point = &(resample_c->post_points[resample_c->num_strands * 3]);
	
	resample_point[X] = prev_resample_point[X] + (STRAND_EXT_FRACT / old_segment->length) * old_segment->disp[X];
	resample_point[Y] = prev_resample_point[Y] + (STRAND_EXT_FRACT / old_segment->length) * old_segment->disp[Y];
	resample_point[Z] = prev_resample_point[Z] + (STRAND_EXT_FRACT / old_segment->length) * old_segment->disp[Z];
	
	
	resample_c->strand_r[resample_c->num_strands] = old_strand->radius;
	resample_c->bundle_i_of_strand[resample_c->num_strands] = old_strand->bundle_i;
	
	
	
	construct_strand(&(resample_c->strands[resample_c->num_strands]), resample_c->num_strands, resample_c->bundle_i_of_strand[resample_c->num_strands], &(resample_c->control_points[start_of_strand_control_point_i]), &(resample_c->start_points[resample_c->num_strands]), &(resample_c->end_points[resample_c->num_strands]), &(resample_c->pre_points[resample_c->num_strands]), &(resample_c->post_points[resample_c->num_strands]), &(resample_c->segments[start_of_strand_control_point_i + 3 * resample_c->num_strands]), resample_c->num_strand_control_points[resample_c->num_strands], 0.0, resample_c->strand_r[resample_c->num_strands]);

}
