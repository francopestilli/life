/*
 *  strand.c
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "phantom/shared/strand.h"
#include "phantom/shared/segment.h"
#include "phantom/shared/shared.h"

void init_strand(Strand *strand, int strand_i, int bundle_i, int num_control_points, double radius, double *start_point, double *end_point) {
	
	double straight_line_length;

	strand->strand_i = strand_i;
	strand->bundle_i = bundle_i;
	strand->num_control_points = num_control_points;
	strand->radius =  radius;
	
	
	straight_line_length = dist_between_points(start_point, end_point);
	
	strand->natural_seg_length = straight_line_length/((double)(num_control_points+1));
	
	strand->bundle = NULL;
	
}


void construct_strand(Strand *strand, int strand_i, int bundle_i, double *control_points, double *start_point, double *end_point, double *pre_point, double *post_point, Segment segments[], int num_strand_control_points, double sample_density, double strand_r) {

	Segment *segment, *prev_segment;
	Control_point point, next_point;
	
	int segment_i, point_i;
	
	double zeros[3] = {0.0, 0.0, 0.0};
	
	strand->radius = strand_r;

	init_strand(strand, strand_i, bundle_i, num_strand_control_points, strand_r, start_point, end_point);
	
	set_control_point(&point, pre_point, zeros); 
	set_control_point(&next_point, start_point, zeros);
	
	segment = &(segments[0]);
	init_segment(segment, -1, NULL, point, next_point, strand, sample_density);
	
	strand->pre_segment = segment;	
	
	point_i = 0;
	
	for (segment_i = 0; segment_i < strand->num_control_points+1; segment_i++) {
	
		point = next_point;
		prev_segment = segment;
	
		if (segment_i == strand->num_control_points)	{							/* If the segment is the last in the strand -> set the next point to be an endpoint. */
			set_control_point(&next_point, end_point, zeros);
		} else {
			set_control_point(&next_point, (control_points + 3 * point_i++), zeros);
		} 		
	
		segment = &(segments[segment_i+1]);								
		init_segment(segment, segment_i, prev_segment, point, next_point, strand, sample_density);
		

	}
	
	strand->start_segment = strand->pre_segment->next_segment;
	strand->end_segment = segment;
	
	prev_segment = segment;
	point = next_point;
	
	set_control_point(&next_point, post_point, zeros);
	
	segment = &(segments[segment_i+1]);
	
	init_segment(segment, segment_i, prev_segment, point, next_point, strand, sample_density);
	
	strand->post_segment = segment;
	

}

