/*
 *  rand_init.c
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
 */


#include <stdlib.h>
#include <stdio.h>
#include <dirent.h>
#include <sys/stat.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "phantom/rand_init/rand_init.h"




/* Handles the input/output and variable initialisation.  All the 'work' is performed by 'rand_init_collection'*/
int main(int argc, char *argv[]) {


	double sphere_r,control_point_freq, strand_r_lbound, strand_r_ubound, strand_r_buffer_ratio;
	unsigned int seed;
	int max_attempts, error;


	Isotropic_region *isotropic_regions;
	
	double grey_matter_r_lbound, grey_matter_r_ubound, grey_matter_location_lbound, grey_matter_location_ubound;
	double grey_matter_diffusivity, grey_matter_baseline_signal, grey_matter_weighting;
	
	int num_grey_matters, num_loaded_isotropic_regions, num_isotropic_regions;

	Strand_collection c;
	DIR *dir;
	FILE *file;
	char *path, *param_path, *output_dir_path, *isotropic_regions_path, trailing_slash[2];
	

	int num_params, num_lines;
	char line[200], key[200];
	double value;

	gsl_rng *rand_gen;

	sphere_r = SPHERE_R_DEFAULT;
	control_point_freq = CONTROL_POINT_FREQ_DEFAULT;
	strand_r_lbound = STRAND_R_LBOUND_DEFAULT;
	strand_r_ubound = STRAND_R_UBOUND_DEFAULT;
	strand_r_buffer_ratio = STRAND_R_BUFFER_RATIO_DEFAULT;

	num_grey_matters = NUM_GREY_MATTERS_DEFAULT;
	grey_matter_r_lbound = GREY_MATTER_R_LBOUND_DEFAULT;
	grey_matter_r_ubound = GREY_MATTER_R_UBOUND_DEFAULT;
	grey_matter_location_lbound = GREY_MATTER_LOCATION_LBOUND_DEFAULT;
	grey_matter_location_ubound = GREY_MATTER_LOCATION_UBOUND_DEFAULT;
	
	grey_matter_diffusivity = GREY_MATTER_DIFFUSIVITY_DEFAULT;
	grey_matter_baseline_signal = GREY_MATTER_BASELINE_SIGNAL_DEFAULT;
	grey_matter_weighting = GREY_MATTER_WEIGHTING_DEFAULT; 

	seed = time(NULL);
	max_attempts = MAX_ATTEMPTS_DEFAULT;
	

	if (argc <= 1 || argc >= 5) {
		printf("\nNumerical Fibre Generator (%s): 'rand_init'\n\nNumber of suplied arguments %d.\n\nUsage: \n\t arg[1]: output directory \n\t arg[2]: parameters file \n\n", VERSION_NUM, (argc-1));

		exit(EXIT_FAILURE);

	} else if (argc == 2) {
		printf("\n\nNo parameters file was supplied... using defaults \n\n");
		
	} else if (argc >= 3) {
	
		param_path = argv[2];
		
		if ((file = fopen(param_path, "r")) == NULL) {
			printf("Error! Could not open parameters file %s\n", param_path);
			exit(EXIT_FAILURE);
		}

		num_params = 0;
		num_lines = 0;
		
		printf("\n\nReading parameters file: %s\n\n", param_path);
		fflush(stdout);
		
		while (fgets(line, 200, file) != NULL) {
			
			if (sscanf(line, "%s %lf", key, &value) != EOF) {
				
				if (!strcmp(key, "sphere_r")) {
					sphere_r = value;
					printf("Read value for sphere_r, %g\n", value);
					
				} else if (!strcmp(key, "control_point_freq")) {
					control_point_freq = value;
					printf("Read value for control_point_freq, %g\n", value);
				
				} else if (!strcmp(key, "strand_r_lbound")) {
					strand_r_lbound = value;
					printf("Read value for strand_r_lbound, %g\n", value);
				
				} else if (!strcmp(key, "strand_r_ubound")) {
					strand_r_ubound = value;
					printf("Read value for strand_r_ubound, %g\n", value);
				
				} else if (!strcmp(key, "strand_r_buffer_ratio")) {
					strand_r_buffer_ratio = value;
					printf("Read value for strand_r_buffer_ratio, %g\n", value);
				
				} else if (!strcmp(key, "num_grey_matters")) {
					num_grey_matters = (int)value;
					printf("Read value for num_grey_matters, %d\n", num_grey_matters);
				
				} else if (!strcmp(key, "grey_matter_r_lbound")) {
					grey_matter_r_lbound = value;
					printf("Read value for grey_matter_r_lbound, %g\n", value);
				
				} else if (!strcmp(key, "grey_matter_r_ubound")) {
					grey_matter_r_ubound = value;
					printf("Read value for grey_matter_r_ubound, %g\n", value);
				
				} else if (!strcmp(key, "grey_matter_location_lbound")) {
					grey_matter_location_lbound = value;
					printf("Read value for grey_matter_location_lbound, %g\n", value);
				
				} else if (!strcmp(key, "grey_matter_location_ubound")) {
					grey_matter_location_ubound = value;
					printf("Read value for grey_matter_location_ubound, %g\n", value);
				
				} else if (!strcmp(key, "grey_matter_diffusivity")) {
					grey_matter_diffusivity = value;
					printf("Read value for grey_matter_diffusivity, %g\n", value);
				
				} else if (!strcmp(key, "grey_matter_baseline_signal")) {
					grey_matter_baseline_signal = value;
					printf("Read value for grey_matter_baseline_signal, %g\n", value);
				
				} else if (!strcmp(key, "grey_matter_weighting")) {
					grey_matter_weighting = value;
					printf("Read value for grey_matter_weighting, %g\n", value);
				
				} else if (!strcmp(key, "seed")) {
				
					if (sscanf(line, "%*s %u", &seed) != 1) {
						printf("Could not read seed value.  Randomly generating...");
					} 
				
					printf("Read value for seed, %d\n", seed);
				
				} else if (!strcmp(key, "max_attempts")) {
					max_attempts = (int)value;
					printf("Read value for max_attempts, %d\n", max_attempts);
				
				} else {
					num_params--;
				}
				fflush(stdout);
				num_params++;
			}
			
			num_lines++;
		}
		
		printf("\nRead in %d (of 14) key/value pairs\n", num_params);
		
	}
	
	 
	if (argc == 4) {
	
		isotropic_regions_path = argv[3];	
		
		num_loaded_isotropic_regions = count_isotropic_regions(isotropic_regions_path);
		
		num_isotropic_regions = num_loaded_isotropic_regions + num_grey_matters;
		
		isotropic_regions = (Isotropic_region*)malloc(sizeof(Isotropic_region) * num_isotropic_regions);
	
		load_isotropic_regions(isotropic_regions, isotropic_regions_path , num_loaded_isotropic_regions);
	
	} else {
	
		num_loaded_isotropic_regions = 0;
		
		num_isotropic_regions = num_grey_matters;
		
		isotropic_regions = (Isotropic_region*)malloc(sizeof(Isotropic_region) * num_isotropic_regions);
	
	}
	
	fflush(stdout);

	output_dir_path = argv[1];


     rand_gen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(rand_gen, seed);

	if (num_grey_matters > 0) {
		printf("\nGenerating strands...\n");
		fflush(stdout);	
	
		rand_init_grey_matters(&(isotropic_regions[num_loaded_isotropic_regions]), num_grey_matters, grey_matter_r_lbound, grey_matter_r_ubound, grey_matter_location_lbound, grey_matter_location_ubound, grey_matter_diffusivity, grey_matter_baseline_signal, grey_matter_weighting, rand_gen);
	}

	printf("\nGenerating strands...\n");
	fflush(stdout);	

		
	rand_init_collection(&c, sphere_r, control_point_freq, strand_r_lbound, strand_r_ubound, strand_r_buffer_ratio, rand_gen, max_attempts, num_isotropic_regions, isotropic_regions);
	
	free(isotropic_regions);
		
	printf("\nWriting %d strand files to output directory (arg 1): %s\n", c.num_strands, output_dir_path);
	fflush(stdout);
	
	if (save_collection(&c, output_dir_path, 0)) {
		printf("Could not save collection %s, exiting ...", output_dir_path);
		exit(EXIT_FAILURE);
	}
	
	
	int strand_i;

	for (strand_i = 0; strand_i < c.num_strands; strand_i++) {
	
		printf("Created strand %5d - radius %f - length %f\n", strand_i, c.strand_r[strand_i], dist_between_points(&(c.start_points[strand_i * 3]), &(c.end_points[strand_i * 3])) );
		fflush(stdout);
	}
	
	
	if (output_dir_path[strlen(output_dir_path) - 1] == DIR_SEP) {
		strcpy(trailing_slash, "");
	} else {
		strcpy(trailing_slash, DIR_SEP_S);
	}
	
	path = (char*)malloc(sizeof(char) * (strlen(output_dir_path) + 100));	
	sprintf(path, "%s%sparameters%s", output_dir_path, trailing_slash, DIR_SEP_S);
	
	if ((dir = opendir(path)) == NULL) {
		
		if ((error = mkdir(path, MY_PERMS)) > 0) {
			printf("Error! Could not create directory %s, (Error code: %d)!\n", path, error);
			return -1;
		}
	} else {
		closedir(dir);
	}
		
	
	strcat(path, "1_rand_init_param.txt");
	
	if ( (file =fopen(path,"w")) == NULL) {
		printf("\nCould not open file %s.!\n", path);
		exit(EXIT_FAILURE);
	} 

	fprintf(file, "seed %u\n", seed);	
	fprintf(file, "sphere_r %g\n", sphere_r);
	fprintf(file, "control_point_freq %g\n", control_point_freq);
	fprintf(file, "strand_r_lbound %g\n", strand_r_lbound);
	fprintf(file, "strand_r_ubound %g\n", strand_r_ubound);
	fprintf(file, "strand_r_buffer_ratio %g\n", strand_r_buffer_ratio);
	fprintf(file, "max_attempts %d\n", max_attempts);
	fprintf(file, "num_grey_matters %d\n", num_grey_matters);
	fprintf(file, "grey_matter_r_lbound %g\n", grey_matter_r_lbound);
	fprintf(file, "grey_matter_r_ubound %g\n", grey_matter_r_ubound);
	fprintf(file, "grey_matter_location_lbound %g\n", grey_matter_location_lbound);
	fprintf(file, "grey_matter_location_ubound %g\n", grey_matter_location_ubound);
	
	fprintf(file, "grey_matter_diffusivity %g\n", grey_matter_diffusivity);
	fprintf(file, "grey_matter_baseline_signal %g\n", grey_matter_baseline_signal);
	fprintf(file, "grey_matter_weighting %g\n", grey_matter_weighting);
	
	fclose(file);
		
	exit(EXIT_SUCCESS);
}



