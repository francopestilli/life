/*
 *  optimise.c
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
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "phantom/shared/shared.h"
#include "phantom/optimise/cost_function.h"

#include "phantom/optimise/optimise.h"

/* Initialisation and iteration of optimisation algorithm */
void optimise(Strand_collection *c, double cost_term_weights[3], double cost_term_powers[3], double sample_density, int grid_reference_size, double initial_step_size, double step_tol, double optimization_tol, double jitter_dev, unsigned int jitter_seed, int max_iterations, int iterations_save_freq, char *iterations_save_path) {



	gsl_multimin_fdfminimizer *optimizer;
	const gsl_multimin_fdfminimizer_type *optimizer_type;
	gsl_multimin_function_fdf cost_function_struct;
	const gsl_vector *state_vector;
	Cost_function_params cost_function_params;
	
	clock_t start, stop;
	double time_for_iteration, total_time, total_time_min;
	
	size_t num_states;
	int iteration, status;
	double max_dx, dx, dx_sum, dx_sum_sq, dx_av, dx_2norm;
	
	unsigned int state_i;
	
	char *iterations_next_save_path, *iterations_next_dirname;
	
	int max_iter_label_width;

	
	gsl_rng *rand_gen;
	
	Sample_block *samples_store; 
	Reference_block *grid_reference;
	
	max_iter_label_width = num_decimal_places(max_iterations);
	
	iterations_next_save_path = (char*)malloc(sizeof(char) * (strlen(iterations_save_path) + 20 + max_iter_label_width));
	
	/*num_cost_function_calls = 0;	*/
				
	/* Allocate memory structures, which will be reused between successive iterations.  Prevents the algorithm from repetitively asking the OS for large blocks of memory, which can become a problem for very large numbers of strands (hence large memory blocks) */
	    					
	allocate_storage_blocks(&samples_store, &grid_reference, grid_reference_size);	
		
	/* Initialise input parameters of the GSL optimiser (implementing the L-BFGS algorithm) */			 
						 		 
	num_states = c->num_control_points * 3;
	optimizer_type = gsl_multimin_fdfminimizer_vector_bfgs2;
	optimizer = gsl_multimin_fdfminimizer_alloc(optimizer_type, num_states);
	
	set_cost_function_params(&cost_function_params, c->start_points, c->end_points, c->pre_points, c->post_points, c->num_strands, c->num_control_points, c->num_strand_control_points, c->fov, c->strand_r, sample_density,  cost_term_weights, cost_term_powers, samples_store, grid_reference, grid_reference_size, c->isotropic_regions, c->num_isotropic_regions); 
	
	cost_function_struct.n = num_states;
	
	/* This is where our custom cost function ('cost_function' in 'cost_function.c') is passed to the optimisation algorithm */ 
	cost_function_struct.f = &cost_function_f;
	cost_function_struct.df = &cost_function_df;
	cost_function_struct.fdf = &cost_function_fdf;
	cost_function_struct.params = (void *)&cost_function_params;
	
	state_vector = gsl_vector_alloc(num_states);
	
	//printf("State vector stride: %d\n", (int)(state_vector->stride));
	
	/* In order to move the state away from any shallow parts of the cost function (which are not ideal points to start the algorithm from), the state is slightly perturbed (by a gaussian distributed random variable) away from any highly regular configuration it may have started in (i.e. all strands are straight (from 'rand_init') or neighbouring strands are hexagonally packed (from 'subdiv') )*/ 
	
	rand_gen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(rand_gen, jitter_seed);
		
	for (state_i=0; state_i < num_states; state_i++) {
		state_vector->data[state_i * state_vector->stride] = c->control_points[state_i] + gsl_rng_uniform(rand_gen) * jitter_dev;
	}
	
	gsl_multimin_fdfminimizer_set(optimizer, &cost_function_struct, state_vector, initial_step_size, step_tol);
	
	
	total_time = 0.0;
	
	iteration = 0;
	printf("\n Iteration\t                Func Val.\t               |Gradient|\t             dx Avg.\t              dx Max\t                |dx|\t     Time taken (s)\tNum. cost function calls\n");
	fflush(stdout);
	printf("------------\t-------------------------\t-------------------------\t--------------------\t--------------------\t--------------------\t-------------------\t------------------------\n");
	fflush(stdout);
	printf("            \t                         \t                         \t                    \t                    \t                    \t                   \t");
	fflush(stdout);


	
	
	iterations_next_dirname = (char*)malloc(sizeof(char) * (17 + max_iter_label_width));

	
	do {
		
		start = clock();
	
		/* Perform an iteration */
		iteration++;
	/*	num_cost_function_calls = 0;*/
		
		gsl_multimin_fdfminimizer_iterate(optimizer);
		
		status = gsl_multimin_test_gradient(optimizer->gradient, optimization_tol);	
		
		max_dx = 0.0;
		

		/* Determine sum statistics about the iteration */
		dx_sum_sq = 0;
		dx_sum = 0;
		for (state_i=0; state_i < num_states; state_i++) {
			dx = optimizer->x->data[state_i] - c->control_points[state_i];
			c->control_points[state_i] = optimizer->x->data[state_i];
			
			if (dx < 0) {
				dx = -dx;
			}
			
			dx_sum += dx;
			dx_sum_sq += dx * dx;
			
			if (dx > max_dx) {
				max_dx = dx;
			}
		}
		
		dx_av = dx_sum/num_states;
		dx_2norm = sqrt(dx_sum_sq);
		
		stop = clock();
	
		time_for_iteration = (double)(stop - start) / CLOCKS_PER_SEC;
	
		total_time += time_for_iteration;
				
		printf("\n%10d\t%25f\t%25f\t%20g\t%20g\t%20g\t%19g\t", iteration, optimizer->f, gsl_blas_dnrm2 (optimizer->gradient), dx_av, max_dx, dx_2norm, time_for_iteration);
		fflush(stdout);
		
		/* Save the configuration after the iteration */

		if ( iterations_save_freq && (iteration % iterations_save_freq == 0)) {
			sprintf(iterations_next_save_path, "%s%citeration_%0*d", iterations_save_path, DIR_SEP, max_iter_label_width, iteration);
						
									
															
			save_collection(c, iterations_next_save_path,0);
		}

		/* Stop the algorithm if there has been no improvement */
		if (max_dx == 0.0) {
			break;
		}

		
	} while (iteration < max_iterations && status == GSL_CONTINUE);
	
	total_time_min = total_time / 60.0;
	
//	if (status == GSL_SUCCESS) {
//		printf("\n\nMinimum was found successfully after %d iterations in %g seconds (%g min):\n", iteration, total_time, total_time_min);
//	} else { 
//		printf("\n\nDid not find minimum after %d iterations in %g seconds (%g min).\n", iteration, total_time, total_time_min);
//	}
//	
	
	printf("\n\nOptimisation terminated after %d iterations in %g seconds (%g minutes).\n",iteration, total_time, total_time_min); 
	
	fflush(stdout);

	/* copy the output of the optimisation back to the original strand collection pointers to give access to the calling function */

	for (state_i=0; state_i < num_states; state_i++) {
		c->control_points[state_i] = gsl_vector_get(optimizer->x, state_i);
	}



	free(iterations_next_dirname);
	
	gsl_multimin_fdfminimizer_free(optimizer);
	
	free_storage_blocks(samples_store, grid_reference, grid_reference_size);

}


void allocate_storage_blocks(Sample_block **start_block, Reference_block **grid_reference, int grid_reference_size) {
	
	int x, y, z;
	
	*start_block = (Sample_block*)malloc(sizeof(Sample_block));
			
	*grid_reference = (Reference_block*)malloc(sizeof(Reference_block)* grid_reference_size * grid_reference_size * grid_reference_size );

	
	init_sample_block(*start_block);
	
	for (z = 0; z < grid_reference_size; z++) {
		for (y = 0; y < grid_reference_size; y++) {
			for (x = 0; x < grid_reference_size; x++) {
				init_reference_block(&(*grid_reference)[grid_reference_size * grid_reference_size * z + grid_reference_size * y + x]);
			}
		}
	}
	
	return;

}


void free_storage_blocks(Sample_block *start_block, Reference_block *grid_reference, int grid_reference_size) {
	
	int x, y, z;
	
	free_sample_blocks(start_block);
	
	for (x = 0; x < grid_reference_size; x++) {
		for (y = 0; y < grid_reference_size; y++) {
			for (z = 0; z < grid_reference_size; z++) {
				free_reference_blocks(&grid_reference[grid_reference_size*grid_reference_size*x + grid_reference_size*y + z]);
			}
		}
	}
	
	free(start_block);
	free(grid_reference);

}


