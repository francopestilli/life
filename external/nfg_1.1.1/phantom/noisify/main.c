/*
 *  noisify.c
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
#include <dirent.h>
#include <sys/stat.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "phantom/noisify/noisify.h"
#include "phantom/shared/shared.h"


int main(int argc, char **argv) {

	double noise_level;

	DIR *dir;
	FILE *file;
	struct stat file_status;
	int image_size, count;
	struct dirent *directory_entry;

	char *param_path,  *input_dir_path, *output_dir_path, *path, *path2, *file_ext;

	int num_params, num_lines;
	char line[200], key[200], filename[200], output_format[200];
	double value;

	float *image;
	gsl_rng *rand_gen;


	int stage_num, error;

	rand_gen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(rand_gen, time(NULL));

	if (argc != 4) {

		printf("\nNumerical Fibre Generator (%s): 'noisify'\n\nNumber of suplied arguments %d.\n\n Usage: \n\t arg[1]: input directory \n\t arg[2]: output directory \n\t arg[3]: parameters file\n\n", VERSION_NUM, (argc-1));

		exit(EXIT_FAILURE);
	}

	param_path = argv[3];

	input_dir_path = argv[1];
	output_dir_path = argv[2];

	if ((file = fopen(param_path, "r")) == NULL) {
		printf("Error! Could not open parameters file %s\n", param_path);
		exit(EXIT_FAILURE);
	}


	noise_level = NOISE_LEVEL_DEFAULT;
	strcpy(output_format, OUTPUT_FORMAT_DEFAULT);

	if ((dir = opendir(output_dir_path)) == NULL) {
		if ((error = mkdir(output_dir_path, MY_PERMS)) > 0) {
			printf("Error! Could not create directory %s, (Error code: %d)!\n", argv[0], error);
			exit(EXIT_FAILURE);
		}
	} else {
		closedir(dir);
	}


	num_params = 0;
	while (fgets(line, 200, file) != NULL) {

		if (sscanf(line, "%s %lf", key, &value) != EOF) {

			if (!strcmp(key, "noise_level")) {
				noise_level = value;
				printf("Read value for noise_level, %g\n", noise_level);

			}  else if (!strcmp(key, "output_format")) {

				if (sscanf(line, "%*s %s", output_format) != 1) {
					printf("Could not read output format.  Using default");
					strcpy(output_format, OUTPUT_FORMAT_DEFAULT);
				}

				printf("Read value for output_format, %s\n", output_format);

			} else {
				num_params--;
			}
			fflush(stdout);
			num_params++;
		}

		num_lines++;
	}

	printf("\nRead in %d (of 3) key/value pairs\n", num_params);



	printf("\n\nAdding noise to images...\n");
	fflush(stdout);

	if (strcmp(output_format, "analyze") == 0) {

		count = 0;

		path = (char*)malloc(sizeof(char) * (max_int(strlen(input_dir_path), strlen(output_dir_path)) + 200));


		if ((dir = opendir(input_dir_path)) == NULL) {
			printf("Error! Could not open directory %s!\n", input_dir_path);
			exit(EXIT_FAILURE);
		}

		while ( (directory_entry = readdir(dir)) != NULL) {

			strcpy(filename, directory_entry->d_name);

			file_ext = strrchr(filename, '.');

			if ( (file_ext != NULL) && (strlen(file_ext) == 4)  && ( strcmp(file_ext , ".hdr") == 0) ) {

				sprintf(path, "%s%c%s", input_dir_path, DIR_SEP, filename);

				file_copy(path, output_dir_path);

				strcpy(strrchr(path, '.'), ".img");


				stat(path, &file_status);
				image = (float*)malloc(file_status.st_size);

				image_size = file_status.st_size / sizeof(float);

				if ((file = fopen(path, "r")) == NULL) {
					printf("Error! Could not open file %s\n", path);
					exit(EXIT_FAILURE);
				}


				fread(image, sizeof(char), file_status.st_size, file);

				fclose(file);


				sprintf(path, "%s%c%s", output_dir_path, DIR_SEP, filename);

				strcpy(strrchr(path, '.'), ".img");

				noisify(image, image_size, noise_level, rand_gen);



				if ( (file =fopen(path,"w")) == NULL) {
					printf("\nCould not open file %s.!\n", path);
					exit(EXIT_FAILURE);
				}


				fwrite(image, sizeof(float), image_size , file);

				fclose(file);

				printf("Added Rician noise to %s.\n", path);
				fflush(stdout);
				count++;
			}

		}

		free(path);

	} else {
		printf("Unrecognied file format specified in input parameters file, %s, exiting ...\n\n", output_format);
		exit(EXIT_FAILURE);
	}

	if (count == 0) {
		printf("\n!!! Warning no images were found in directory %s!\n\n", input_dir_path);
		fflush(stdout);
	}


	printf("\n\nAdded noise to %d images.\n", count);
	fflush(stdout);

	/* Copy parameter files to new directory for record of what has been performed */

	stage_num = copy_parameters(input_dir_path, output_dir_path);

	if (stage_num != -1) {

		path = (char*)malloc(sizeof(char) * (max_int(strlen(output_dir_path), strlen(input_dir_path)) + 200));
		path2 = (char*)malloc(sizeof(char) * (max_int(strlen(output_dir_path), strlen(input_dir_path)) + 200));

		sprintf(path, "%s%cparameters%c%d_noisify_param.txt", output_dir_path, DIR_SEP, DIR_SEP, stage_num);

		if ( (file =fopen(path,"w")) == NULL) {
			printf("\nCould not open file %s.!\n", path);
			free(path);
			free(path2);
			exit(EXIT_SUCCESS);
		}


		fprintf(file, "noise_level %g\n", noise_level);
		fprintf(file, "output_format %s\n", output_format);

		if (fclose(file)) {
			printf("Error closing parameters file copy\n");
		}



		free(path);
		free(path2);
	}

	exit(EXIT_SUCCESS);
}

