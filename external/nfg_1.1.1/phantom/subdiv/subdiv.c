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

#include "phantom/subdiv/subdiv.h"


void subdivide_collection(Strand_collection *children, Strand_collection *parents, double strand_r_final) {

	int parent_i;
	int ubound_num_children1, ubound_num_control_points, ubound_num_cols, ubound_num_rows, ubound_new_children;

	ubound_num_children1 = 0;		//An upper bound on the number of children that could be generated
	ubound_num_control_points = 0;		//An upper bound on the total number of points that could be generated for these children.
	
	for (parent_i=0; parent_i < parents->num_strands; parent_i++) {

		ubound_num_cols = 2 * floor(parents->strands[parent_i].radius/strand_r_final) + 1;
		ubound_num_rows = floor(2.0 * parents->strands[parent_i].radius/(sqrt(3) * strand_r_final)) + 1;
		
		ubound_new_children = ubound_num_rows * ubound_num_cols;
		
		ubound_num_children1 += ubound_new_children;
		
		/* add 10 (should only need to add 1) just to be sure and allow leeway for future modifications of the subdivide_strand function */
		ubound_num_control_points += ubound_new_children * (parents->strands[parent_i].num_control_points + 10);
	
	}

	collection_alloc(children, ubound_num_children1, ubound_num_control_points, parents->num_isotropic_regions);
	copy_isotropic_regions(children, parents->isotropic_regions, parents->num_isotropic_regions);

	children->num_strands = 0;
	children->num_control_points = 0;
		
		
	for (parent_i=0; parent_i < parents->num_strands; parent_i++) {
	
		subdivide_strand(&(parents->strands[parent_i]), children, strand_r_final);
	
	}	


	children->strand_r = (double*)realloc(children->strand_r, sizeof(double) * children->num_strands);
	children->num_strand_control_points = (int*)realloc(children->num_strand_control_points, sizeof(int) * children->num_strands);
	children->control_points = (double*)realloc(children->control_points, sizeof(double) * 3 * children->num_control_points);
	children->control_points_grad = (double*)realloc(children->control_points_grad, sizeof(double) * children->num_control_points);
	children->start_points = (double*)realloc(children->start_points, sizeof(double) * 3 * children->num_strands);
	children->end_points = (double*)realloc(children->end_points, sizeof(double) * 3 * children->num_strands);
	children->pre_points = (double*)realloc(children->pre_points, sizeof(double) * 3 * children->num_strands);
	children->post_points = (double*)realloc(children->post_points, sizeof(double) * 3 * children->num_strands);	
	children->strands = (Strand*)realloc(children->strands, sizeof(Strand) * children->num_strands);
	children->segments = (Segment*)realloc(children->segments, sizeof(Segment) * (children->num_control_points + 3 * children->num_strands)); 

}


void subdivide_strand(Strand *parent, Strand_collection *children, double strand_r_final) {


	double *parent_start_point_pos, possible_child_start_point[3];
	double row_basis[3], col_basis[3], row_basis_length, col_basis_length;
	
	Segment *parent_segment; 

	int half_num_rows, half_num_cols, row_i, col_i, start_of_strand_control_point_i;
	
	int num_new_children, total_num_new_children, debug;

	double col_indent;
	double child_disp[3], child_dist, *prev_child_point, *child_point, s[3], tang[3], v[3], mu, child_segment_disp[3], child_segment_length, child_segment_direction[3], child_angle, prev_child_segment_direction[3];
	double tang_norm;
	
	double proj_onto_parent;
	
	/* Get the vector tangent (angle halfway between pre_segment and start_segment) to the parent start point */
	tang[X] = parent->pre_segment->disp[X] / parent->pre_segment->length + parent->start_segment->disp[X] / parent->start_segment->length;
	tang[Y] = parent->pre_segment->disp[Y] / parent->pre_segment->length + parent->start_segment->disp[Y] / parent->start_segment->length;	
	tang[Z] = parent->pre_segment->disp[Z] / parent->pre_segment->length + parent->start_segment->disp[Z] / parent->start_segment->length;		
	
	
	tang_norm = vector_norm(tang);
	
	if (tang_norm < STARTING_ANGLE_THRESHOLD) {
		printf("Error the pre and start segments of strand %d, double back on each other. Not able to proceed with subdivision.  Please check data and try again\n", parent->strand_i);
		exit(EXIT_FAILURE); 
		
	} else if (abs(tang[X]/tang_norm) < LINEAR_DEPEND_THRESHOLD) {
	
		/*Cross the first tangent with [1;0;0] to get a vector perpendicular to the parent->pre_segment to be the row basis vector child start_points;*/
		row_basis[X] = 0;
		row_basis[Y] = tang[Z];
		row_basis[Z] = -tang[Y];
	
	} else {
	
		/*If tangent is colinear with [1;0;0] cross the first tangent with [0;1;0] instead to get a vector perpendicular to the parent->pre_segment to be the row basis vector child start_points;*/
		row_basis[X] = -tang[Z];
		row_basis[Y] = 0;
		row_basis[Z] = tang[X];
	}


	/*Cross the row basis with the tangent vector to get the col basis.*/
	col_basis[X] = row_basis[Y] * tang[Z] - row_basis[Z] * tang[Y];
	col_basis[Y] = row_basis[Z] * tang[X] - row_basis[X] * tang[Z];
	col_basis[Z] = row_basis[X] * tang[Y] - row_basis[Y] * tang[X];
	
	
	
	

	/*Normalize then rescale the basis vectors to the correct row and column increment (2 * strand_r_final and sqrt(3) * strand_r_final respectively)*/
	row_basis_length = sqrt(3) * strand_r_final / sqrt(row_basis[X] * row_basis[X] + row_basis[Y] * row_basis[Y] + row_basis[Z] * row_basis[Z]);
	col_basis_length = 2 * strand_r_final / sqrt(col_basis[X] * col_basis[X] + col_basis[Y] * col_basis[Y] + col_basis[Z] * col_basis[Z]);


			
	row_basis[X] = row_basis[X] * row_basis_length;
	row_basis[Y] = row_basis[Y] * row_basis_length;
	row_basis[Z] = row_basis[Z] * row_basis_length;
	
	col_basis[X] = col_basis[X] * col_basis_length;
	col_basis[Y] = col_basis[Y] * col_basis_length;
	col_basis[Z] = col_basis[Z] * col_basis_length;

	/* The number of children rows and cols tried (to see if they fall in the original strand radius) either side of the parent start_point */
	half_num_rows = floor(parent->radius/(sqrt(3) * strand_r_final)) + 1;
	half_num_cols = floor(parent->radius/(2 * strand_r_final) );
			
	parent_start_point_pos = parent->start_segment->start_point.pos;	
		
	num_new_children = 0;

	for (row_i = -half_num_rows; row_i <= half_num_rows; row_i++) {
	
		col_indent = ((double)(row_i % 2)) / 2;		//indent the colums by 1/2 if it is an odd row. 
		
		for (col_i = -half_num_cols ; col_i <= half_num_cols; col_i++) {
		
			possible_child_start_point[X] = parent_start_point_pos[X] + ((double)row_i) * row_basis[X] + (((double)col_i) + col_indent) * col_basis[X];
			possible_child_start_point[Y] = parent_start_point_pos[Y] + ((double)row_i) * row_basis[Y] + (((double)col_i) + col_indent) * col_basis[Y];
			possible_child_start_point[Z] = parent_start_point_pos[Z] + ((double)row_i) * row_basis[Z] + (((double)col_i) + col_indent) * col_basis[Z];
			
			child_disp[X] = possible_child_start_point[X] - parent_start_point_pos[X];
			child_disp[Y] = possible_child_start_point[Y] - parent_start_point_pos[Y];
			child_disp[Z] = possible_child_start_point[Z] - parent_start_point_pos[Z];
			
			child_dist = sqrt(child_disp[X] * child_disp[X] + child_disp[Y] * child_disp[Y] + child_disp[Z] * child_disp[Z]);
			
			/* If possible child start point lies within parent strand extent then add it to the collection */
			if (child_dist <= parent->radius) {
				
				children->start_points[ (children->num_strands + num_new_children) * 3 + X] = possible_child_start_point[X];
				children->start_points[ (children->num_strands + num_new_children) * 3 + Y] = possible_child_start_point[Y];
				children->start_points[ (children->num_strands + num_new_children) * 3 + Z] = possible_child_start_point[Z];
				
				children->pre_points[ (children->num_strands + num_new_children) * 3 + X] = possible_child_start_point[X] + parent->pre_segment->disp[X];
				children->pre_points[ (children->num_strands + num_new_children) * 3 + Y] = possible_child_start_point[Y] + parent->pre_segment->disp[Y];
				children->pre_points[ (children->num_strands + num_new_children) * 3 + Z] = possible_child_start_point[Z] + parent->pre_segment->disp[Z];
				
				/*children->num_strand_control_points[children->num_strands + num_new_children] = parent->num_control_points;*/
				children->strand_r[children->num_strands + num_new_children] = strand_r_final;
				children->bundle_i_of_strand[children->num_strands + num_new_children] = parent->bundle_i;
			
				num_new_children++;
			}
		}
	}

	total_num_new_children = children->num_strands + num_new_children;
 
	for (; children->num_strands < total_num_new_children; children->num_strands++) {




		parent_segment = parent->start_segment;
		prev_child_point = &(children->start_points[(children->num_strands) * 3]);
		
		prev_child_segment_direction[X] = parent->pre_segment->disp[X]/(parent->pre_segment->length);				
		prev_child_segment_direction[Y] = parent->pre_segment->disp[Y]/(parent->pre_segment->length);				
		prev_child_segment_direction[Z] = parent->pre_segment->disp[Z]/(parent->pre_segment->length);
		
		start_of_strand_control_point_i = children->num_control_points;
		children->num_strand_control_points[children->num_strands] = 0; 
		
		child_point = &(children->control_points[(children->num_control_points) * 3]);
		
		/* While the current parent segment is not the post segment */
		while (parent_segment != parent->post_segment) {
		
		
		
			child_point = &(children->control_points[(children->num_control_points) * 3]);
			
			/* Normalised vector of the segment displacement */
			s[X] = parent_segment->disp[X] / parent_segment->length;
			s[Y] = parent_segment->disp[Y] / parent_segment->length;
			s[Z] = parent_segment->disp[Z] / parent_segment->length;						
						
						
			/* Get the vector tangent (angle halfway between current segment and the next_segment) to the parent strand */
			tang[X] = s[X] + parent_segment->next_segment->disp[X] / parent_segment->next_segment->length;
			tang[Y] = s[Y] + parent_segment->next_segment->disp[Y] / parent_segment->next_segment->length;	
			tang[Z] = s[Z] + parent_segment->next_segment->disp[Z] / parent_segment->next_segment->length;		
      

      
			
      if ((tang[X] == 0) && (tang[Y] == 0.0) && (tang[Z] == 0.0)) 
        printf("Warning! Consecutive segments (%d and %d of strand %d) completely double back on each other (will result in a NaN).\n", parent_segment->segment_i , parent_segment->next_segment->segment_i, parent_segment->strand->strand_i); 
      
			/* The displacement vector between the previous child_point and the next parent control_point */
			v[X] = parent_segment->end_point.pos[X] - prev_child_point[X];
			v[Y] = parent_segment->end_point.pos[Y] - prev_child_point[Y];			
			v[Z] = parent_segment->end_point.pos[Z] - prev_child_point[Z];
			
			/* The distance along s, it intersects with the normal plane at the control point */ 
			mu = dot_product(v,tang)/dot_product(s,tang);
			
			/* Calculate the control point value */
			child_point[X] = prev_child_point[X] + mu * s[X];
			child_point[Y] = prev_child_point[Y] + mu * s[Y];
			child_point[Z] = prev_child_point[Z] + mu * s[Z];
			
			child_segment_disp[X] = child_point[X] - prev_child_point[X];						
			child_segment_disp[Y] = child_point[Y] - prev_child_point[Y];			
			child_segment_disp[Z] = child_point[Z] - prev_child_point[Z];
			
			child_segment_length = vector_norm(child_segment_disp);
			
			child_segment_direction[X] = child_segment_disp[X]/child_segment_length;
			child_segment_direction[Y] = child_segment_disp[Y]/child_segment_length;			
			child_segment_direction[Z] = child_segment_disp[Z]/child_segment_length; 

      
			child_angle = dot_product(child_segment_direction, prev_child_segment_direction);
			
			if ((children->num_strands == 573) && (children->num_strand_control_points[children->num_strands] == 34)) {
				debug = 0;
				debug = 2;
			
			
			}
		
			proj_onto_parent = dot_product(parent_segment->disp, child_segment_disp);
			
			proj_onto_parent /= (parent_segment->length * parent_segment->length);			
					
			
			/* Check to see if the child segment doesn't double back on to itself, which is caused by "wobbling" in the parent strand coupled with a large bundle radius*/
			if (child_angle > 0) {
			

				
		
				/*To ensure an even spacing between the child segments (due to bunching/stretching on the inside/outside of bends) if the new segment is more than 1.5 times the parent segment length then a new control point is inserted in the middle. If the new segment is less than half the parent segment length then the control point is ignored*/

				if (proj_onto_parent > 0.5) {
					children->num_control_points++;
					children->num_strand_control_points[children->num_strands]++;
					prev_child_point = child_point;

					/* Save to check if next segment 'doubles back' */				
					prev_child_segment_direction[X] = child_segment_direction[X];				
					prev_child_segment_direction[Y] = child_segment_direction[Y];				
					prev_child_segment_direction[Z] = child_segment_direction[Z];
				}
			}
			
			/* Iterate */			
			parent_segment = parent_segment->next_segment;
			
		}

		/* Reassign the last generated control to the end_point */
		children->end_points[(children->num_strands * 3) + X] = child_point[X];		
		children->end_points[(children->num_strands * 3) + Y] = child_point[Y];		
		children->end_points[(children->num_strands * 3) + Z] = child_point[Z];
		
		children->num_control_points--;
		children->num_strand_control_points[children->num_strands]--;
		
		/* Add the parent post segment to get the child strand post_point */
		children->post_points[(children->num_strands * 3) + X] = child_point[X] + parent->post_segment->disp[X];		
		children->post_points[(children->num_strands * 3) + Y] = child_point[Y] + parent->post_segment->disp[Y];		
		children->post_points[(children->num_strands * 3) + Z] = child_point[Z] + parent->post_segment->disp[Z];
		
		construct_strand(&(children->strands[children->num_strands]), children->num_strands, parent->bundle_i, &(children->control_points[start_of_strand_control_point_i]), &(children->start_points[children->num_strands]), &(children->end_points[children->num_strands]), &(children->pre_points[children->num_strands]), &(children->post_points[children->num_strands]), &(children->segments[start_of_strand_control_point_i + 3 * children->num_strands]), children->num_strand_control_points[children->num_strands], 0.0, strand_r_final);		
	}


}

