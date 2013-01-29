/*
 *  cost_function.c
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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_multimin.h>


#include "phantom/optimise/cost_function.h"
#include "phantom/shared/segment.h"
#include "phantom/shared/strand.h"
#include "phantom/shared/control_point.h"
#include "phantom/optimise/sample_block.h"
#include "phantom/optimise/reference_block.h"

#include "phantom/shared/shared.h"

/*#define BIG_NUMBER 1e20*/

/* Calculates the cost (the returned variable), and the cost function gradient (control_points_grad) */
double cost_function(
					double *control_points_grad, 
					double *control_points, 
					double *start_points, 
					double *end_points,
					double *pre_points,
					double *post_points,
					int num_strands, 
					int num_control_points, 
					int num_strand_control_points[], 
					double fov, 
					double strand_r[], 
					double sample_density, 
					double cost_term_weights[],
					double cost_term_powers[], 
					Sample_block *samples_store, 
					Reference_block *grid_reference,
					int grid_reference_size,
					Isotropic_region *isotropic_regions,
					int num_isotropic_regions
				) {
				
	int x, y, z, strand_i, segment_i, point_i, cumulative_segment_i, isotropic_region_i;
	
	Control_point point, next_point;
	
	/* These gradients are not actually used but given storage space to simplifiy the looping structure in which the end_points are treated the same as any other control point.  Note that these values can be used if desired */
	double start_points_grad[num_strands][3], end_points_grad[num_strands][3],  pre_points_grad[num_strands][3], post_points_grad[num_strands][3];
	
	
	Strand strands[num_strands], *strand;
	
	/* There is one segment per point plus one at the end of each strand. */
	Segment *segments, start_segments[num_strands], end_segments[num_strands], *segment, *prev_segment;		
	
	
	double repulsion_cost, bend_cost, length_cost, cost;
	
	int print_error;
	int nan_detect;
	
	segments = (Segment*)malloc(sizeof(Segment) * (num_control_points + num_strands));
	
	nan_detect = 0;
	print_error = 1;
	
/*
	if (num_cost_function_calls > max_num_cost_function_calls) {
		// If in calculating the iteration, the minimizer exceeds maximum number of cost function calls, short circuit the cost function calculation until it exceeds its maximum number attempts (set to 100 in 'gsl/multimin/linear_minimize.c').  This is an ad-hoc method of lowering this maximum number of attempts to avoid lengthy attempts to lower the cost function before giving up.
		
		for (point_i = 0; point_i < num_control_points; point_i) {
			control_points_grad[point_i] = BIG_NUMBER;
		}
	
		return BIG_NUMBER;
	}
	
	num_cost_function_calls++;
*/


	/*---------------------------------------------------------------------------------------------
	//
	//	Initialization of various variables
	//
	//---------------------------------------------------------------------------------------------	*/
	
		
	cost =0.0;
			
	reset_sample_block(samples_store);
	
	for (z=0; z < grid_reference_size; z++) {	
		for (y=0; y < grid_reference_size; y++) {
			for (x=0; x < grid_reference_size; x++) {
				reset_reference_block(&(grid_reference[z * grid_reference_size * grid_reference_size + y * grid_reference_size + x])); 
			}
		}
	}

	for (isotropic_region_i = 0; isotropic_region_i < num_isotropic_regions; isotropic_region_i++) {
	
		add_isotropic_region_samples(&(isotropic_regions[isotropic_region_i]), &samples_store, grid_reference, fov, grid_reference_size);
		
	}


	for (point_i=0; point_i<num_control_points; point_i++) {
		control_points_grad[point_i * 3 + X] = 0;
		control_points_grad[point_i * 3 + Y] = 0;	
		control_points_grad[point_i * 3 + Z] = 0;
	}
	
	
	for (strand_i=0; strand_i<num_strands; strand_i++) {
		start_points_grad[strand_i][X] = 0;
		start_points_grad[strand_i][Y] = 0;
		start_points_grad[strand_i][Z] = 0;
		
		end_points_grad[strand_i][X] = 0;
		end_points_grad[strand_i][Y] = 0;
		end_points_grad[strand_i][Z] = 0;

		pre_points_grad[strand_i][X] = 0;
		pre_points_grad[strand_i][Y] = 0;
		pre_points_grad[strand_i][Z] = 0;
		
		post_points_grad[strand_i][X] = 0;
		post_points_grad[strand_i][Y] = 0;
		post_points_grad[strand_i][Z] = 0;
		
	}
	
	
	
	point_i = 0;
	cumulative_segment_i = 0;
	
	/*---------------------------------------------------------------------------------------------
	 *
	 *	Loop through strands and generate segments between the points.
	 *
	 *---------------------------------------------------------------------------------------------*/


	for (strand_i=0; strand_i < num_strands; strand_i++) {
		
		strand = &(strands[strand_i]);
		init_strand(strand, strand_i, strand_i, num_strand_control_points[strand_i], strand_r[strand_i], &(start_points[strand_i * 3]), &(end_points[strand_i * 3]));
		
		set_control_point(&point, &(pre_points[strand_i * 3]), pre_points_grad[strand_i] ); 
		set_control_point(&next_point, &(start_points[strand_i * 3]), start_points_grad[strand_i]);
		
		segment = &(start_segments[strand_i]);
		init_segment(segment, -1, NULL, point, next_point, strand, sample_density);
		
		for (segment_i = 0; segment_i < strand->num_control_points+1; segment_i++) {
		
			point = next_point;
			prev_segment = segment;
		
			/* If the segment is the last in the strand -> set the next point to be an endpoint. */
			if (segment_i == strand->num_control_points)	{							
				set_control_point(&next_point, &(end_points[strand_i * 3]), end_points_grad[strand_i]);
			} else {
			
				set_control_point(&next_point, &(control_points[point_i * 3]), &(control_points_grad[ 3 * (point_i)]) );
				point_i++;
			} 		
		
			segment = &(segments[cumulative_segment_i++]);								
			init_segment(segment, segment_i, prev_segment, point, next_point, strand, sample_density);
		
			length_cost = calc_length_cost(segment, cost_term_weights[2], cost_term_powers[2]);
					
			repulsion_cost = calc_repulsion_cost(segment, &samples_store, grid_reference, cost_term_weights[0], cost_term_powers[0], fov, grid_reference_size);
			
			bend_cost = calc_bend_cost(segment, cost_term_weights[1], cost_term_powers[1]);
			
			cost += repulsion_cost + bend_cost + length_cost;
				
			
		}
		
		prev_segment = segment;
		point = next_point;
		
		set_control_point(&next_point, &(post_points[strand_i* 3]), post_points_grad[strand_i] );
		
		segment = &(end_segments[strand_i]);
		
		init_segment(segment, segment_i, prev_segment, point, next_point, strand, sample_density);

		
		bend_cost = calc_bend_cost(segment, cost_term_weights[1], cost_term_powers[1]);
		
		cost += bend_cost;
		
		
	}		 

	
	free(segments);
	
	printf(".");
	fflush(stdout);
	return cost;

}



double cost_function_f(const gsl_vector *state_vector, void *cost_function_params) {
	
	double cost, *grad;
	Cost_function_params *params;	
	params = (Cost_function_params*)cost_function_params;
	
	grad = (double*)malloc(sizeof(double) * params->num_control_points * 3);

	cost = cost_function(grad, state_vector->data, params->start_points, params->end_points, params->pre_points, params->post_points, params->num_strands, params->num_control_points, params->num_strand_control_points, params->fov, params->strand_r, params->sample_density, params->cost_term_weights, params->cost_term_powers, params->samples_store, params->grid_reference, params->grid_reference_size, params->isotropic_regions, params->num_isotropic_regions);
	
	free(grad);
	
	return cost;

}

void cost_function_df(const gsl_vector *state_vector, void *cost_function_params, gsl_vector *grad_vector) {

	Cost_function_params *params;	
	params = (Cost_function_params*)cost_function_params;
	
	cost_function( grad_vector->data,				
				state_vector->data,				
				params->start_points,			
				params->end_points,
				params->pre_points, 
				params->post_points,				 
				params->num_strands,			
				params->num_control_points,				 
				params->num_strand_control_points,		 
				params->fov,				 
				params->strand_r,				 
				params->sample_density,			 
				params->cost_term_weights,
				params->cost_term_powers, 
				params->samples_store,			 
				params->grid_reference,
				params->grid_reference_size,
				params->isotropic_regions,
				params->num_isotropic_regions						
				);
	
}


void cost_function_fdf(const gsl_vector *state_vector, void *cost_function_params, double *cost, gsl_vector *grad_vector) {

	Cost_function_params *params;	
	params = (Cost_function_params*)cost_function_params;
	
	*cost = cost_function(grad_vector->data, state_vector->data, params->start_points, params->end_points, params->pre_points, params->post_points, params->num_strands, params->num_control_points, params->num_strand_control_points, params->fov, params->strand_r, params->sample_density, params->cost_term_weights, params->cost_term_powers, params->samples_store, params->grid_reference, params->grid_reference_size, params->isotropic_regions, params->num_isotropic_regions);
	
}

void set_cost_function_params(Cost_function_params *params, double *start_points, double *end_points, double *pre_points, double *post_points, int num_strands, int num_control_points, int *num_strand_control_points, double fov, double *strand_r, double sample_density,  double *cost_term_weights, double *cost_term_powers, Sample_block *samples_store, Reference_block *grid_reference, int grid_reference_size, Isotropic_region *isotropic_regions, int num_isotropic_regions){

	params->start_points = start_points;
	params->end_points = end_points;
	params->pre_points = pre_points;
	params->post_points = post_points;
	params->num_strands = num_strands;
	params->num_control_points = num_control_points;
	params->num_strand_control_points = num_strand_control_points;
	params->fov = fov;
	params->strand_r = strand_r;
	params->sample_density = sample_density;
	params->cost_term_weights = cost_term_weights;
	params->cost_term_powers = cost_term_powers;
	params->samples_store = samples_store;
	params->grid_reference = grid_reference;
	params->grid_reference_size = grid_reference_size;
	params->isotropic_regions = isotropic_regions;
	params->num_isotropic_regions = num_isotropic_regions;
}






/************************************************************************************************
********************************* Bending Cost **************************************************
************************************************************************************************/


double calc_bend_cost(Segment *segment, double k, double p) {

	double cost, cost_base, cost_U, cost_V;
	
	double grad_scalar;
	
	double *start, *mid, *end, *start_grad, *mid_grad, *end_grad;
	
	double start_grad_part[3],  mid_grad_part[3], end_grad_part[3];
	
	double start_segment_length, end_segment_length, end_segment_scalar, start_segment_scalar;
	
	double cross_section_area;
	
	cross_section_area = segment->strand->radius * segment->strand->radius;

	start = segment->prev_segment->start_point.pos;
	mid = segment->start_point.pos; 
	end = segment->end_point.pos;		
	
	start_segment_length = segment->prev_segment->length;
	end_segment_length = segment->length;
	
	/*Calculate the derivative in two parts (U & V) then combine using the chain rule. */
	
	/* Calculate the dot product between the first segment and the second segment */
	cost_U = (end[X] - mid[X]) * (mid[X] - start[X]) + (end[Y] - mid[Y]) * (mid[Y] - start[Y]) + (end[Z] - mid[Z]) * (mid[Z] - start[Z]);

	/* Calculate the normalizing factor. */
	cost_V = 1/(end_segment_length * start_segment_length);

	cost_base = cost_U * cost_V;
	
//	cost = k * pow(1-cost_base,p);
	cost = cross_section_area * k * pow(1-cost_base,p);
	
	start_grad = segment->prev_segment->start_point.grad;
	mid_grad = segment->start_point.grad;
	end_grad = segment->end_point.grad;	

	end_segment_scalar = -1 * end_segment_length * end_segment_length * end_segment_length * start_segment_length;
	start_segment_scalar = -1 * start_segment_length * start_segment_length * start_segment_length * end_segment_length;
//	grad_scalar = -1 * k * p * pow(1-cost_base, p-1);
	grad_scalar = -1 * cross_section_area * k * p * pow(1-cost_base, p-1);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 

	start_grad_part[X] = grad_scalar * ( (mid[X] - end[X]) * cost_V   +   cost_U * (start[X] - mid[X]) / start_segment_scalar );
	start_grad_part[Y] = grad_scalar * ( (mid[Y] - end[Y]) * cost_V   +   cost_U * (start[Y] - mid[Y]) / start_segment_scalar );
	start_grad_part[Z] = grad_scalar * ( (mid[Z] - end[Z]) * cost_V   +   cost_U * (start[Z] - mid[Z]) / start_segment_scalar );

	mid_grad_part[X] = grad_scalar * ( (-2*mid[X] + start[X] + end[X]) * cost_V   +   cost_U * ((mid[X] - start[X]) / start_segment_scalar + (mid[X] - end[X]) / end_segment_scalar) );
	mid_grad_part[Y] = grad_scalar * ( (-2*mid[Y] + start[Y] + end[Y]) * cost_V   +   cost_U * ((mid[Y] - start[Y]) / start_segment_scalar + (mid[Y] - end[Y]) / end_segment_scalar) );
	mid_grad_part[Z] = grad_scalar * ( (-2*mid[Z] + start[Z] + end[Z]) * cost_V   +   cost_U * ((mid[Z] - start[Z]) / start_segment_scalar + (mid[Z] - end[Z]) / end_segment_scalar) );
		
	end_grad_part[X] = grad_scalar * ( (mid[X] - start[X]) * cost_V   +   cost_U * (end[X] - mid[X]) / end_segment_scalar );
	end_grad_part[Y] = grad_scalar * ( (mid[Y] - start[Y]) * cost_V   +   cost_U * (end[Y] - mid[Y]) / end_segment_scalar );
	end_grad_part[Z] = grad_scalar * ( (mid[Z] - start[Z]) * cost_V   +   cost_U * (end[Z] - mid[Z]) / end_segment_scalar );
	

	start_grad[X] += start_grad_part[X];
	start_grad[Y] += start_grad_part[Y];
	start_grad[Z] += start_grad_part[Z];
	
	mid_grad[X] += mid_grad_part[X];
	mid_grad[Y] += mid_grad_part[Y];
	mid_grad[Z] += mid_grad_part[Z];
	
	end_grad[X] += end_grad_part[X];
	end_grad[Y] += end_grad_part[Y];
	end_grad[Z] += end_grad_part[Z];
	
	return cost;
	
}
























/************************************************************************************************
********************************** Length Cost **************************************************
************************************************************************************************


*******************************************************************
* Calculate the cost due to the segment length. 
********************************************************************


--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
Atomic Length Cost 
------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------

 
									     p
length_cost :=  k ( (len  - L)/L )


--------------------------------------------------------------------------------------------------------
 dCost/dstart_x 
------------------------------------------------------------------------------------------------------

 
        /          \p
        |len  -   L|
      k |----------|  p (start_x - end_x)
        \    L    /
   ---------------------------------------
                 len    (len    - L)


			  ________________________________________________________________________________________________________________
			 /		2                            2                            2
  len :=	/  end_x  - 2 start_x end_x + start_x  - 2 start_y end_y + start_y + end_y  - 2 start_z end_z + end_z  + start_z
		  \/	


segment->length => len
segment->strand->natural_seg_length => L
length_diff => len - L

*/



double calc_length_cost(Segment *segment, double k, double p) {

	double grad[3];
	double length_diff;
	double length_cost;
	
	double length_grad_scalar;

	double cross_section_area;
	
	cross_section_area = segment->strand->radius * segment->strand->radius;


	length_diff = segment->length - segment->strand->natural_seg_length;	
	
//	length_cost = k * pow(length_diff/segment->strand->natural_seg_length, p);
	length_cost = cross_section_area * k * pow(length_diff/segment->strand->natural_seg_length, p);

	
	/* To avoid a divide by zero error when both the length_cost and the length_diff go to zero */
	if (length_cost != 0) {
//		length_grad_scalar =  k * pow(length_diff,p-1) * p / (segment->length * pow(segment->strand->natural_seg_length, p)) ;
		length_grad_scalar =  cross_section_area * k * pow(length_diff,p-1) * p / (segment->length * pow(segment->strand->natural_seg_length, p)) ;
	} else {
		length_grad_scalar = 0;
	}
	
	grad[X] = -1.0 * segment->disp[X] * length_grad_scalar;
	grad[Y] = -1.0 * segment->disp[Y] * length_grad_scalar;
	grad[Z] = -1.0 * segment->disp[Z] * length_grad_scalar;
	
	segment->start_point.grad[X] += grad[X];
	segment->start_point.grad[Y] += grad[Y];
	segment->start_point.grad[Z] += grad[Z];
	
	segment->end_point.grad[X] -= grad[X];
	segment->end_point.grad[Y] -= grad[Y];
	segment->end_point.grad[Z] -= grad[Z];	
	

	return length_cost;

}




/************************************************************************************************
******************************* Repulsion Cost **************************************************
************************************************************************************************/



double calc_repulsion_cost(Segment *segment, Sample_block **samples_store, Reference_block *grid_reference, double k, double p, double fov, int grid_reference_size) {

	int x, y, z, sample_i;
	
	int cell_ubound[3], cell_lbound[3];
	Sample *new_sample, *exist_sample;	
	Reference_block *cell_reference;
	
	double repulsion_cost, repulsion_cost_tmp;

	
	/*******************************************************************
	* Calculate the cost due to the repulsion between different strands 
	********************************************************************/

	repulsion_cost = 0;
	
	for (sample_i=0; sample_i < segment->num_samples; sample_i++) {
		
		new_sample = new_Sample(samples_store, segment, sample_i);
		
		cell_ubound[X] = ceil( (new_sample->pos[X] + fov + (segment->strand)->radius) * 0.5 * grid_reference_size/fov );
		cell_ubound[Y] = ceil( (new_sample->pos[Y] + fov + (segment->strand)->radius) * 0.5 * grid_reference_size/fov );
		cell_ubound[Z] = ceil( (new_sample->pos[Z] + fov + (segment->strand)->radius) * 0.5 * grid_reference_size/fov );
		
		cell_lbound[X] = floor( (new_sample->pos[X] + fov - (segment->strand)->radius) * 0.5 * grid_reference_size/fov );
		cell_lbound[Y] = floor( (new_sample->pos[Y] + fov - (segment->strand)->radius) * 0.5 * grid_reference_size/fov );
		cell_lbound[Z] = floor( (new_sample->pos[Z] + fov - (segment->strand)->radius) * 0.5 * grid_reference_size/fov );
		
		
		for (z = max_int(cell_lbound[Z], 0); z < min_int(cell_ubound[Z], grid_reference_size); z++) {		
			for (y = max_int(cell_lbound[Y], 0); y < min_int(cell_ubound[Y], grid_reference_size); y++) {
				for (x = max_int(cell_lbound[X], 0); x < min_int(cell_ubound[X], grid_reference_size); x++) {
										
					cell_reference = (grid_reference + z * grid_reference_size * grid_reference_size + y * grid_reference_size + x); 
					
					
					exist_sample = first_reference(&cell_reference);				/* Note that next_sample will loop through the linked-list of Reference_blocks and update the cell_reference pointers to the current Reference_block until it is the last Reference_block in the list. */
					
					while (exist_sample != NULL) {
						
						repulsion_cost_tmp = calc_sample_cost(new_sample, exist_sample, k, p);
						
						repulsion_cost += repulsion_cost_tmp;
						
						
						exist_sample = next_reference(&cell_reference);
					}
					
					add_reference(&cell_reference, new_sample);
				}
			}
			
		}
		
	}	

	
	return repulsion_cost;

}






/************************************************************************************************
********************************** Repulsion cost between 2 samples *****************************
************************************************************************************************/




double calc_sample_cost(Sample *new_sample, Sample *exist_sample, double k, double p) {
	
	double cost,  dist, grad_u_scalar, overlap_dist, overlap;
	double disp[3];
	double cost_u, cost_v, cost_w, grad_u[3], grad_v[3], grad_w[3];
	double new_start_grad[3], new_end_grad[3], exist_start_grad[3], exist_end_grad[3];
	
	cost = 0;

	
	if ((new_sample->segment)->strand != (exist_sample->segment)->strand && exist_sample->last_accessed_sample != new_sample) {
		
		disp[X] = new_sample->pos[X] - exist_sample->pos[X];
		disp[Y] = new_sample->pos[Y] - exist_sample->pos[Y];
		disp[Z] = new_sample->pos[Z] - exist_sample->pos[Z];
		
				
		dist = sqrt(disp[X] * disp[X] + disp[Y] * disp[Y] + disp[Z] * disp[Z]); 
				
		overlap_dist = ((new_sample->segment)->strand)->radius + ((exist_sample->segment)->strand)->radius;
		
		if (dist < overlap_dist) {
						
			overlap = overlap_dist - dist;

			cost_u = k * pow(overlap/overlap_dist, p);
			
			cost_v = new_sample->segment->sample_length;
			
			cost_w = exist_sample->segment->sample_length;
			
			cost = cost_u * cost_v * cost_w;
			
			grad_u_scalar = -p * k * pow (overlap, p-1) / ( pow(overlap_dist,p) * dist); 	                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
			
			grad_u[X] = disp[X] * grad_u_scalar;
			grad_u[Y] = disp[Y] * grad_u_scalar;
			grad_u[Z] = disp[Z] * grad_u_scalar;
					
			/* The partial derivatives of the sample_length have been precalculated for efficiency */		
			grad_v[X] = new_sample->segment->sample_length_grad[X];
			grad_v[Y] = new_sample->segment->sample_length_grad[Y];
			grad_v[Z] = new_sample->segment->sample_length_grad[Z];
			
			grad_w[X] = exist_sample->segment->sample_length_grad[X];
			grad_w[Y] = exist_sample->segment->sample_length_grad[Y];
			grad_w[Z] = exist_sample->segment->sample_length_grad[Z];
			

			/* After the gradient between to samples is calculated it is then assigned back to their generating points. */
			new_start_grad[X] = grad_u[X] * (1.0-new_sample->sample_fract) * cost_v * cost_w    +    cost_u * grad_v[X] * cost_w;
			new_start_grad[Y] = grad_u[Y] * (1.0-new_sample->sample_fract) * cost_v * cost_w    +    cost_u * grad_v[Y] * cost_w;
			new_start_grad[Z] = grad_u[Z] * (1.0-new_sample->sample_fract) * cost_v * cost_w    +    cost_u * grad_v[Z] * cost_w;
			
			new_end_grad[X] = grad_u[X] * new_sample->sample_fract * cost_v * cost_w    -    cost_u * grad_v[X] * cost_w;
			new_end_grad[Y] = grad_u[Y] * new_sample->sample_fract * cost_v * cost_w    -    cost_u * grad_v[Y] * cost_w;
			new_end_grad[Z] = grad_u[Z] * new_sample->sample_fract * cost_v * cost_w    -    cost_u * grad_v[Z] * cost_w;
			
			exist_start_grad[X] = -1*grad_u[X] * (1.0-exist_sample->sample_fract) * cost_v * cost_w    +    cost_u * cost_v * grad_w[X];
			exist_start_grad[Y] = -1*grad_u[Y] * (1.0-exist_sample->sample_fract) * cost_v * cost_w    +    cost_u * cost_v * grad_w[Y];
			exist_start_grad[Z] = -1*grad_u[Z] * (1.0-exist_sample->sample_fract) * cost_v * cost_w    +    cost_u * cost_v * grad_w[Z];
			
			exist_end_grad[X] = -1*grad_u[X] * exist_sample->sample_fract * cost_v * cost_w    -    cost_u * cost_v * grad_w[X];
			exist_end_grad[Y] = -1*grad_u[Y] * exist_sample->sample_fract * cost_v * cost_w    -    cost_u * cost_v * grad_w[Y];
			exist_end_grad[Z] = -1*grad_u[Z] * exist_sample->sample_fract * cost_v * cost_w    -    cost_u * cost_v * grad_w[Z];



			new_sample->segment->start_point.grad[X] += new_start_grad[X];
			new_sample->segment->start_point.grad[Y] += new_start_grad[Y];
			new_sample->segment->start_point.grad[Z] += new_start_grad[Z];
			
			new_sample->segment->end_point.grad[X] += new_end_grad[X];
			new_sample->segment->end_point.grad[Y] += new_end_grad[Y];
			new_sample->segment->end_point.grad[Z] += new_end_grad[Z];
			
			exist_sample->segment->start_point.grad[X] += exist_start_grad[X];
			exist_sample->segment->start_point.grad[Y] += exist_start_grad[Y];
			exist_sample->segment->start_point.grad[Z] += exist_start_grad[Z];
			
			exist_sample->segment->end_point.grad[X] += exist_end_grad[X];
			exist_sample->segment->end_point.grad[Y] += exist_end_grad[Y];
			exist_sample->segment->end_point.grad[Z] += exist_end_grad[Z];

			
		}
		
		exist_sample->last_accessed_sample = new_sample;
		
	}
	
	return cost;
	
}


void add_isotropic_region_samples(Isotropic_region *isotropic_region, Sample_block **samples_store, Reference_block *grid_reference, double fov, int grid_reference_size) {

	int x, y, z;
	int cell_ubound[3], cell_lbound[3];
	Reference_block *cell_reference;
	Sample *new_sample;

	new_sample = new_Sample(samples_store, &(isotropic_region->segment), 0);

	new_sample->pos[X] = isotropic_region->pos[X];
	new_sample->pos[Y] = isotropic_region->pos[Y];
	new_sample->pos[Z] = isotropic_region->pos[Z];	
	
	new_sample->sample_fract = 0;
	new_sample->sample_i = -1;	//For debugging purposes to determine whether the current sample is a isotropic region.
	
	cell_ubound[X] = ceil( (new_sample->pos[X] + fov + isotropic_region->radius) * 0.5 * grid_reference_size/fov );
	cell_ubound[Y] = ceil( (new_sample->pos[Y] + fov + isotropic_region->radius) * 0.5 * grid_reference_size/fov );
	cell_ubound[Z] = ceil( (new_sample->pos[Z] + fov + isotropic_region->radius) * 0.5 * grid_reference_size/fov );
	
	cell_lbound[X] = floor( (new_sample->pos[X] + fov - isotropic_region->radius) * 0.5 * grid_reference_size/fov );
	cell_lbound[Y] = floor( (new_sample->pos[Y] + fov - isotropic_region->radius) * 0.5 * grid_reference_size/fov );
	cell_lbound[Z] = floor( (new_sample->pos[Z] + fov - isotropic_region->radius) * 0.5 * grid_reference_size/fov );
		
		
	for (z = max_int(cell_lbound[Z], 0); z < min_int(cell_ubound[Z], grid_reference_size); z++) {		
		for (y = max_int(cell_lbound[Y], 0); y < min_int(cell_ubound[Y], grid_reference_size); y++) {
			for (x = max_int(cell_lbound[X], 0); x < min_int(cell_ubound[X], grid_reference_size); x++) {
									
				cell_reference = (grid_reference + z * grid_reference_size * grid_reference_size + y * grid_reference_size + x); 
				
				add_reference(&cell_reference, new_sample);
			}
		}
	}

	

}
