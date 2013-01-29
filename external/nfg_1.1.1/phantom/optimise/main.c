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

#include "phantom/optimise/optimise.h"





/* Handles the input/output and variable initialisation.  All the 'work' is performed by 'rand_init_collection'*/
int main(int argc, char **argv) {

	int error;
	double sample_density, initial_step_size, step_tol, optimization_tol, fov, jitter_dev, relative_isotropic_repulsion_weight;
	int grid_reference_size, max_iterations, iterations_save_freq;
	unsigned int jitter_seed;


	Strand_collection c;
	DIR *dir;
	FILE *param_file, *file;
	char *param_path, *input_dir_path,  *output_dir_path, *iterations_save_path, *path;
	int stage_num;

	int num_params, num_lines;
	char line[200], key[200];
	double value;

	double cost_term_weights[3] = COST_TERM_WEIGHTS_DEFAULT;
	double cost_term_powers[3] = COST_TERM_POWERS_DEFAULT;
	relative_isotropic_repulsion_weight = RELATIVE_ISOTROPIC_REPULSION_WEIGHT_DEFAULT;
	sample_density = SAMPLE_DENSITY_DEFAULT;
	grid_reference_size = GRID_REFERENCE_SIZE_DEFAULT;
	initial_step_size = INITIAL_STEP_SIZE_DEFAULT;
	step_tol = STEP_TOL_DEFAULT;
	optimization_tol = OPTIMIZATION_TOL_DEFAULT;
	max_iterations = MAX_ITERATIONS_DEFAULT;
/*	max_num_cost_function_calls = MAX_NUM_COST_FUNCTION_CALLS_DEFAULT;*/
	jitter_dev = JITTER_STDEV_DEFAULT;
	jitter_seed = JITTER_SEED_DEFAULT;
	iterations_save_freq = ITERATIONS_SAVE_FREQ_DEFAULT;
	fov = FOV_DEFAULT;

	param_path = "";


	if (argc <= 2 || argc >= 5) {
		printf("\nNumerical Fibre Generator (%s): 'optimise'\n\nNumber of suplied arguments %d.\n\n Usage: \n\t arg[1]: input directory \n\t arg[2]: output directory \n\t arg[3]: parameters file\n\n", VERSION_NUM, (argc-1));
		exit(EXIT_FAILURE);


	} else if (argc == 4) {

		param_path = argv[3];


		if ((param_file = fopen(param_path, "r")) == NULL) {
			printf("Error! Could not open file %s\n", param_path);
			exit(EXIT_FAILURE);
		}


		num_params = 0;
		num_lines = 0;

		printf("\n\nReading parameters file: %s\n\n", param_path);
		fflush(stdout);

		while (fgets(line, 200, param_file) != NULL) {

			if (sscanf(line, "%s %lf", key, &value) != EOF) {

				if (!strcmp(key, "fov")) {
					fov = value;
					printf("Read value for fov, %g\n", value);

				} else if (!strcmp(key, "cost_repulsion_weight")) {
					cost_term_weights[0] = value;
					printf("Read value for cost_repulsion_weight, %g\n", value);

				}  else if (!strcmp(key, "cost_bending_weight")) {
					cost_term_weights[1] = value;
					printf("Read value for cost_bending_weight, %g\n", value);

				}  else if (!strcmp(key, "cost_length_weight")) {
					cost_term_weights[2] = value;
					printf("Read value for cost_length_weight, %g\n", value);

				} else if (!strcmp(key, "cost_repulsion_power")) {
					cost_term_powers[0] = value;
					printf("Read value for cost_repulsion_power, %g\n", value);

				} else if (!strcmp(key, "cost_bending_power")) {
					cost_term_powers[1] = value;
					printf("Read value for cost_bending_power, %g\n", value);

				} else if (!strcmp(key, "cost_length_power")) {
					cost_term_powers[2] = value;
					printf("Read value for cost_length_power, %g\n", value);

				} else if (!strcmp(key, "relative_isotropic_repulsion_weight")) {
					relative_isotropic_repulsion_weight = value;
					printf("Read value for relative_isotropic_repulsion_weight, %g\n", value);

				}  else if (!strcmp(key, "sample_density")) {
					sample_density = value;
					printf("Read value for sample_density, %g\n", value);

				} else if (!strcmp(key, "grid_reference_size")) {
					grid_reference_size = value;
					printf("Read value for grid_reference_size, %g\n", value);

				} else if (!strcmp(key, "initial_step_size")) {
					initial_step_size = value;
					printf("Read value for initial_step_size, %g\n", value);

				} else if (!strcmp(key, "step_tol")) {
					step_tol = value;
					printf("Read value for step_tol, %g\n", value);

				} else if (!strcmp(key, "optimization_tol")) {
					optimization_tol = value;
					printf("Read value for optimization_tol, %g\n", value);

				} else if (!strcmp(key, "jitter_dev")) {
					jitter_dev = value;
					printf("Read value for jitter_dev, %g\n", value);

				} else if (!strcmp(key, "jitter_seed")) {
					jitter_seed = (unsigned int)value;
					printf("Read value for jitter_seed, %u\n", jitter_seed);

				} else if (!strcmp(key, "max_iterations")) {
					max_iterations = (int)value;
					printf("Read value for max_iterations, %d\n", max_iterations);

				} /*else if (!strcmp(key, "max_num_cost_function_calls")) {
					max_num_cost_function_calls = (int)value;
					printf("Read value for max_num_cost_function_calls, %d\n", max_num_cost_function_calls);

				} */ else if (!strcmp(key, "iterations_save_freq")) {

					iterations_save_freq = (int)value;
					printf("Read value for iterations_save_freq, %d\n", iterations_save_freq);

				} else {
					num_params--;
				}
				fflush(stdout);
				num_params++;
			}

			num_lines++;
		}

		printf("\nRead in %d (of possible 16) key/value pairs\n", num_params);

	}

	fflush(stdout);


	input_dir_path = argv[1];
	output_dir_path = argv[2];

	if (load_collection(&c, input_dir_path)) {
		printf("Could not load strands from directory %s, exiting ...\n", input_dir_path);
		exit(EXIT_FAILURE);
	}

	/* If the fov supplied is greater than the furthest extent of the loaded strand then use it otherwise use the furthest extent */
	if (fov > c.fov) {
		c.fov = fov;
	}

	iterations_save_path = (char*)malloc(sizeof(char) * (strlen(output_dir_path) + 20));

	if (iterations_save_freq != 0) {

		if ((dir = opendir(output_dir_path)) == NULL) {
			if ((error = mkdir(output_dir_path, MY_PERMS)) > 0) {
				printf("Error! Could not create directory %s, (Error code: %d)!\n", output_dir_path, error);
				return 1;
			}
		} else {
			closedir(dir);
		}



		sprintf(iterations_save_path, "%s%citerations", output_dir_path, DIR_SEP);

		if ((dir = opendir(iterations_save_path)) == NULL) {
			if ((error = mkdir(iterations_save_path, MY_PERMS)) > 0) {
				printf("Error! Could not create directory %s, (Error code: %d)!\n", argv[0], error);
				exit(EXIT_FAILURE);
			}
		} else {
			closedir(dir);
		}
	}




	printf("\nOptimising strand control points...\n");
	fflush(stdout);

	set_isotropic_regions_weightings(&c, relative_isotropic_repulsion_weight);

	/* The optimisatioin algorithm */
	if (c.num_strands == 0) {
		printf("Warning! No strands to optimise, skipping optimisation\n");
	} else {
		optimise(&c, cost_term_weights, cost_term_powers, sample_density, grid_reference_size, initial_step_size, step_tol, optimization_tol, jitter_dev, jitter_seed, max_iterations, iterations_save_freq, iterations_save_path);
	}

	printf("\nWriting %d strand files to: %s...\n", c.num_strands, output_dir_path);
	fflush(stdout);


	if (save_collection(&c, output_dir_path, 0)) {
		printf("\nCould not save strands to directory %s, exiting ...\n", output_dir_path);
		exit(EXIT_FAILURE);

	}


	stage_num = copy_parameters(input_dir_path, output_dir_path);

	if (stage_num != -1) {
		path = (char*)malloc(sizeof(char) * (strlen(output_dir_path) + 80));

		sprintf(path, "%s%cparameters%c%d_optimise_param.txt", output_dir_path, DIR_SEP, DIR_SEP, stage_num);

		if ( (file =fopen(path,"w")) == NULL) {
			printf("\nCould not open parameters save file %s.!\n", path);
		}

		fprintf(file, "fov %g\n", fov);

		fprintf(file, "cost_repulsion_weight %g\n", cost_term_weights[0]);
		fprintf(file, "cost_bending_weight %g\n", cost_term_weights[1]);
		fprintf(file, "cost_length_weight %g\n", cost_term_weights[2]);

		fprintf(file, "cost_repulsion_power %g\n", cost_term_powers[0]);
		fprintf(file, "cost_bending_power %g\n", cost_term_powers[1]);
		fprintf(file, "cost_length_power %g\n", cost_term_powers[2]);

		fprintf(file, "relative_isotropic_repulsion_weight %g\n", relative_isotropic_repulsion_weight);

		fprintf(file, "sample_density %g\n", sample_density);
		fprintf(file, "grid_reference_size %d\n", grid_reference_size);
		fprintf(file, "initial_step_size %g\n", initial_step_size);
		fprintf(file, "step_tol %g\n", step_tol);
		fprintf(file, "optimization_tol %g\n", optimization_tol);
		fprintf(file, "jitter_dev %g\n", jitter_dev);
		fprintf(file, "jitter_seed %u\n", jitter_seed);
		fprintf(file, "max_iterations %d\n", max_iterations);
		fprintf(file, "iterations_save_freq %d\n", iterations_save_freq);


		if (fclose(file)) {
			printf("Error closing parameters file copy\n");
			exit(EXIT_FAILURE);
		}

		free(path);
	}

	free(iterations_save_path);


	exit(EXIT_SUCCESS);

}

