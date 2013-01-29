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

#include "phantom/resample/resample.h"


int main(int argc, char *argv[]) {

	int error;
	double resample_length, double_back_angle_threshold, forward_angle_threshold;

	Strand_collection old_c, resample_c;
	DIR *dir;	
	FILE *param_file, *file;
	char *param_path, *input_dir_path,  *output_dir_path, *path;
	

	int num_params, num_lines, stage_num;
	char line[200], key[200];
	double value;


	if (argc != 4) {	
		printf("\nNumerical Fibre Generator (%s): 'resample'\n\nNumber of suplied arguments %d.\n\n Usage: \n\t arg[1]: input directory \n\t arg[2]: output directory \n\t arg[3]: parameters file\n\n", VERSION_NUM, (argc-1));
		exit(EXIT_FAILURE);
	}

	param_path = argv[3];
	
	input_dir_path = argv[1];
	output_dir_path = argv[2];
	
	
	if ((param_file = fopen(param_path, "r")) == NULL) {
		printf("Error! Could not open file %s\n", param_path);
		exit(EXIT_FAILURE);
	}

	if ((dir = opendir(output_dir_path)) == NULL) {
		if ((error = mkdir(output_dir_path, MY_PERMS)) > 0) {
			printf("Error! Could not create directory %s, (Error code: %d)!\n", argv[0], error);
			exit(EXIT_FAILURE);
		}
	} else {
		closedir(dir);
	}	
	
	
	resample_length = RESAMPLE_LENGTH_DEFAULT;
	double_back_angle_threshold = DOUBLE_BACK_ANGLE_THRESHOLD_DEFAULT;
	forward_angle_threshold = FORWARD_ANGLE_THRESHOLD_DEFAULT;
		
	num_params = 0;
	num_lines = 0;
	
	printf("\n\nReading parameters file: %s\n\n", param_path);
	fflush(stdout);
	
	while (fgets(line, 200, param_file) != NULL) {
		
		if (sscanf(line, "%s %lf", key, &value) != EOF) {
			
			if (!strcmp(key, "resample_length")) {
				resample_length = value;
				printf("Read value for resample_length, %g\n", value);
			
			} else if (!strcmp(key, "double_back_angle_threshold")) {
				double_back_angle_threshold = value;
				printf("Read value for double_back_angle_threshold, %g\n", value);
			
			} else if (!strcmp(key, "forward_angle_threshold")) {
				forward_angle_threshold = value;
				printf("Read value for forward_angle_threshold, %g\n", value);
			
			} else {
				num_params--;
				//printf("Did not recognise key, '%s' on line %d (%f) \n", key, num_lines+1, value);
			}
			fflush(stdout);
			num_params++;
		}
		
		num_lines++;
	}
	
	
	
	printf("\nRead in %d (of 3) key/value pairs\n", num_params); 
	fflush(stdout);

	if (load_collection(&old_c, input_dir_path)) {
		printf("Could not load strands from directory %s, exiting ...\n", input_dir_path);
		exit(EXIT_FAILURE);
	
	}
						
	printf("Resampling strands...\n");	
	fflush(stdout);
		
	resample_collection(&old_c, &resample_c, resample_length, double_back_angle_threshold, forward_angle_threshold);
	
	printf("\nWriting %d strand files to: %s\n", resample_c.num_strands, output_dir_path);
	fflush(stdout);
	
	if (save_collection(&resample_c, output_dir_path, 0)) {
		printf("\nCould not save strands to directory %s, exiting ...\n", output_dir_path);
		exit(EXIT_FAILURE);
	
	}
		
		
	/* Copy parameter files to new directory for record of what has been performed */
	
	stage_num = copy_parameters(input_dir_path, output_dir_path);
	
	if (stage_num != -1) {
		path = (char*)malloc(sizeof(char) * (strlen(output_dir_path) + 80));	

		sprintf(path, "%s%cparameters%c%d_resample_param.txt", output_dir_path, DIR_SEP, DIR_SEP, stage_num);
		
		if ( (file =fopen(path,"w")) == NULL) {
			printf("\nCould not open file %s.!\n", path);
			exit(EXIT_SUCCESS);
		} 
										
																			
		fprintf(file, "resample_length %g\n", resample_length);
		fprintf(file, "double_back_angle_threshold %g\n", double_back_angle_threshold);	
		fprintf(file, "forward_angle_threshold %g\n", forward_angle_threshold);			
		
		
		if (fclose(file)) {
			printf("Error closing parameters file copy\n");
		}
		
		free(path);
	}

	exit(EXIT_SUCCESS);

}
