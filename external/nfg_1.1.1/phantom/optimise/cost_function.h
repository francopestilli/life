/*
 *  cost_function.h
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

#ifndef COST_FUNCTION_H

#define COST_FUNCTION_H 0

#include <gsl/gsl_multimin.h>

#include "phantom/optimise/sample_block.h"
#include "phantom/optimise/reference_block.h"
#include "phantom/shared/isotropic_region.h"

/*int num_cost_function_calls, max_num_cost_function_calls;*/

typedef struct _cost_function_params {
	double *start_points, *end_points, *pre_points, *post_points, *strand_r;
	
	int num_strands, num_control_points, *num_strand_control_points, grid_reference_size;
	
	double fov, sample_density, *cost_term_weights, *cost_term_powers;
	
	Sample_block *samples_store;
	Reference_block *grid_reference;

	int num_isotropic_regions;
	Isotropic_region *isotropic_regions;

} Cost_function_params;

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
				);



double cost_function_f(const gsl_vector *state_vector, void *params);

void cost_function_df(const gsl_vector *state_vector, void *params, gsl_vector *grad);

void cost_function_fdf(const gsl_vector *state_vector, void *params, double *cost, gsl_vector *grad);

void set_cost_function_params(Cost_function_params *params, double *start_points, double *end_points, double *pre_points, double *post_points, int num_strands, int num_control_points, int *num_strand_control_points, double fov, double *strand_r, double sample_density,  double *cost_term_weights, double *cost_term_powers, Sample_block *samples_store, Reference_block *grid_reference, int grid_reference_size, Isotropic_region *isotropic_regions, int num_isotropic_regions);


double calc_bend_cost(Segment *segment, double k, double p);

double calc_repulsion_cost(Segment *segment, Sample_block **samples_store, Reference_block *grid_reference, double k, double p, double fov, int grid_reference_size);

double calc_length_cost(Segment *segment, double k, double p);

double calc_sample_cost(Sample *new_sample, Sample *exist_sample, double k, double p);

void add_isotropic_region_samples(Isotropic_region *isotropic_region, Sample_block **samples_store, Reference_block *grid_reference, double fov, int grid_reference_size);

#endif
