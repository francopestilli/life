/*
 *  mri_sim.c
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


#include "phantom/mri_sim/mri_sim.h"

int main(int argc, char *argv[]) {

	int error;

	double fa, diffusivity;  
	int num_voxels, num_voxels_pdim[3]; 
	double voxel_size, voxel_size_pdim[3], subvoxel_size;
	char output_format[50];
	double num_subvoxels_dble, dummy;
	int num_subvoxels;

	double sphere_r;
	int max_width, grad_i;

	float *images;
	int num_grad_directions, offset;
	double *grad_directions, *b_values, *grad_directions_tmp, *b_values_tmp, grad_norm;

	Strand_collection c;
	DIR *dir;	
	FILE *param_file, *file;
	char *param_path, *input_dir_path,  *output_dir_path, *path, *grad_directions_path;

	int num_params, num_lines, stage_num;
	char line[200], key[200];
	double value;
	int save_subvoxels, save_ext_segment_stats, save_ext_fill_radial_stats;
	double *subvoxels;
	
	int *debug;

	Strand_collection_stats *stats;

	if (argc != 5) {
		printf("\nNumerical Fibre Generator (%s): 'mri_sim'\n\nNumber of suplied arguments %d.\n\n Usage: \n\t arg[1]: input directory \n\t arg[2]: output directory \n\t arg[3]: gradient directions/b-values file  \n\t arg[4]: parameters file\n\n", VERSION_NUM, (argc-1));
		exit(EXIT_FAILURE);
	}

	param_path = argv[4];
	
	input_dir_path = argv[1];
	output_dir_path = argv[2];
	grad_directions_path = argv[3];
	
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
	
	

	
	fa = FA_DEFAULT;
	diffusivity = DIFFUSIVITY_DEFAULT;
	num_voxels = NUM_VOXELS_DEFAULT;
	voxel_size = VOXEL_SIZE_DEFAULT;
	subvoxel_size = SUBVOXEL_SIZE_DEFAULT;
	sphere_r = SPHERE_R_DEFAULT;
	save_ext_segment_stats = SAVE_EXT_SEGMENT_STATS_DEFAULT;
	save_ext_fill_radial_stats = SAVE_EXT_FILL_RADIAL_STATS_DEFAULT;
	save_subvoxels = SAVE_SUBVOXELS_DEFAULT;
	strcpy(output_format, OUTPUT_FORMAT_DEFAULT); 
	
	
	
	num_params = 0;
	num_lines = 0;
	
	printf("\n\nReading parameters file: %s\n\n", param_path);
	fflush(stdout);
	
	
	while (fgets(line, 200, param_file) != NULL) {
		
		if (sscanf(line, "%s %lf", key, &value) != EOF) {
			
			if (!strcmp(key, "fa")) {
				fa = value;
				printf("Read value for fa, %g\n", fa);
				
			}  else if (!strcmp(key, "diffusivity")) {
				diffusivity = value;
				printf("Read value for diffusivity, %g\n", diffusivity);
			
			} else if (!strcmp(key, "num_voxels")) {
				num_voxels = (int)value;
				printf("Read value for num_voxels, %d\n", num_voxels);
			
			} else if (!strcmp(key, "voxel_size")) {
				voxel_size = value;
				printf("Read value for voxel_size, %g\n", voxel_size);
			
			} else if (!strcmp(key, "subvoxel_size")) {
				subvoxel_size = value;
				printf("Read value for subvoxel_size, %g\n", subvoxel_size);
			
			} else if (!strcmp(key, "sphere_r")) {
				sphere_r = value;
				printf("Read value for sphere_r, %g\n", sphere_r);
			
			} else if (!strcmp(key, "save_subvoxels")) {
				save_subvoxels = (int)value;
				printf("Read value for save_subvoxels, %d\n", save_subvoxels);
			
			} else if (!strcmp(key, "save_ext_segment_stats")) {
				save_ext_segment_stats = (int)value;
				printf("Read value for save_ext_segment_stats, %d\n", save_ext_segment_stats);
			
			} else if (!strcmp(key, "save_ext_fill_radial_stats")) {
				save_ext_fill_radial_stats = (int)value;
				printf("Read value for save_ext_fill_radial_stats, %d\n", save_ext_fill_radial_stats);
			
			} else if (!strcmp(key, "output_format")) {
				
				if (sscanf(line, "%*s %s", output_format) != 1) {
					printf("Could not read output format.  Using default");
					strcpy(output_format, OUTPUT_FORMAT_DEFAULT);
				} 
				
				printf("Read value for output_format, %s\n", output_format);
			
			} else {
				num_params--;
			}
			num_params++;
		}
		
		num_lines++;
	}
	
	printf("\nRead in %d (of 10) key/value pairs\n", num_params); 
	
	
	printf("\nLoading gradient directions...\n");
	fflush(stdout);
	
	
	grad_directions = (double*)malloc(sizeof(double) * GRAD_DIR_BLOCK_SIZE * 3);
	b_values = (double*)malloc(sizeof(double) * GRAD_DIR_BLOCK_SIZE);

	if ((file = fopen(grad_directions_path, "r")) == NULL) {
		printf("Error! Could not open gradient directions file %s\n", grad_directions_path);
		exit(EXIT_FAILURE);
	}

	num_grad_directions = 0;

	while (fgets(line, 200, file) != NULL) {
	
	
		offset = num_grad_directions * 3;
		
		if (sscanf(line, "%lf %lf %lf %lf\n", &(grad_directions[offset + X]), &(grad_directions[offset + Y]), &(grad_directions[offset + Z]), &(b_values[num_grad_directions]) ) != 4) {
			if (line[0] == '\n') {
				continue;
			} else {
				printf("Error!! Could not read the %d th direction/b_value in the gradient directions file %s\n", num_grad_directions+1, grad_directions_path);
				exit(EXIT_FAILURE);
			}
		}
		
		grad_norm = vector_norm(&(grad_directions[offset]));
		
		if (grad_norm != 0) {
			grad_directions[offset + X] = grad_directions[offset + X] / grad_norm; 
			grad_directions[offset + Y] = grad_directions[offset + Y] / grad_norm;
			grad_directions[offset + Z] = grad_directions[offset + Z] / grad_norm;				
		}
		
		num_grad_directions++;
		
		if (num_grad_directions % GRAD_DIR_BLOCK_SIZE * 3 == 0) {
			
			grad_directions_tmp = (double*)malloc(sizeof(double) * GRAD_DIR_BLOCK_SIZE * 3 * (num_grad_directions/GRAD_DIR_BLOCK_SIZE + 1));
			memcpy(grad_directions_tmp, grad_directions, (num_grad_directions-1) * 3 * sizeof(double));
			free(grad_directions);
			grad_directions = grad_directions_tmp;
			
			b_values_tmp = (double*)malloc(sizeof(double) * GRAD_DIR_BLOCK_SIZE * (num_grad_directions/GRAD_DIR_BLOCK_SIZE + 1));
			memcpy(b_values_tmp, b_values, (num_grad_directions-1) * sizeof(double));
			free(b_values);
			b_values = b_values_tmp;
		}
	}
		
	/* Load collection */	
		
	if (load_collection(&c, input_dir_path)) {
		printf("Could not load strands from directory %s, exiting ...\n", input_dir_path);
		exit(EXIT_FAILURE);
	
	}
	
	
	debug = (int*)malloc(sizeof(int));
	
	
	c.sphere_r = sphere_r;
																
	printf("\nSimulating DW-MR images...\n");	
	fflush(stdout);	
	
		
	num_subvoxels_dble = voxel_size / subvoxel_size;
	num_subvoxels = (int)ceil(num_subvoxels_dble);
	
	if (modf(num_subvoxels_dble, &dummy) > SUBVOXEL_ROUND_DOWN_WARNING_THRESHOLD) {
		printf("\nWarning subvoxel size was rounded down to %g in order to fit voxel size %g!\n", voxel_size/((double)num_subvoxels), voxel_size);
		fflush(stdout);
	}
		
	stats = strand_collection_stats_alloc(&c, save_ext_segment_stats, save_ext_fill_radial_stats);	
				
	images = mri_sim(&c, fa, diffusivity, num_voxels, voxel_size, num_subvoxels, num_grad_directions, grad_directions, b_values, stats, save_subvoxels, &subvoxels);
	
	
	path = (char*)malloc(sizeof(char) * (strlen(output_dir_path) + 200));


	if (!strcmp(output_format, "analyze")) {

		num_voxels_pdim[X] = num_voxels_pdim[Y] = num_voxels_pdim[Z] = num_voxels;
		voxel_size_pdim[X] = voxel_size_pdim[Y] = voxel_size_pdim[Z] = voxel_size;


		max_width = num_decimal_places(num_grad_directions);

		for (grad_i = 0; grad_i < num_grad_directions; grad_i++) {

			sprintf(path, "%s%cdwi-%0*d.hdr", output_dir_path, DIR_SEP, max_width, grad_i); 
			
			write_analyze_header(num_voxels_pdim, voxel_size_pdim, ANALYZE_FLOAT, sizeof(float), path);
			
			sprintf(path, "%s%cdwi-%0*d.img", output_dir_path, DIR_SEP, max_width, grad_i); 
			
			
			if ( (file =fopen(path,"wb")) == NULL) {
				printf("\nCould not open file %s.!\n", path);
				exit(EXIT_FAILURE);
			} 
			
			fwrite( &(images[grad_i * num_voxels * num_voxels * num_voxels]), sizeof(float),  num_voxels * num_voxels * num_voxels, file);
			
			fclose(file);
			
			printf("Wrote data file %s.\n", path);
			fflush(stdout);
		}
		

		
		
	} else {
		printf("Unrecognied file format specified in input parameters file, %s, exiting ...\n\n", output_format);
		exit(EXIT_FAILURE);
	}
	
	free(images);	
		
		
	if (save_subvoxels) {
	
		sprintf(path, "%s%csubvoxels", output_dir_path, DIR_SEP);
		
		if ((dir = opendir(path)) == NULL) {
			if ((error = mkdir(path, MY_PERMS)) > 0) {
				printf("Error! Could not create directory %s, (Error code: %d)!\n", path, error);
				exit(EXIT_FAILURE);
			}
		} else {
			closedir(dir);
		}
		
		sprintf(path, "%s%csubvoxels%csubvoxels.dat", output_dir_path, DIR_SEP, DIR_SEP);
		
		
		printf("\nWriting orientations to file...\n");	
		fflush(stdout);
		
		if ( (file =fopen(path,"wb")) == NULL) {
			printf("\nCould not open file %s.!\n", path);
			exit(EXIT_FAILURE);
		} 
		
		
		fwrite(subvoxels, sizeof(double), num_voxels * num_subvoxels * num_voxels * num_subvoxels * num_voxels * num_subvoxels * 3 , file);
		
		fclose(file);
		free(subvoxels);
	}	
		
		
	printf("\n\nSimulated DW-MR images.\n\n\n");		
	fflush(stdout);
			
	sprintf(path, "%s%cstats", output_dir_path, DIR_SEP);
		
	if ((dir = opendir(path)) == NULL) {
		if ((error = mkdir(path, MY_PERMS)) > 0) {
			printf("Error! Could not create directory %s, (Error code: %d)!\n", path, error);
			exit(EXIT_FAILURE);
		}
	} else {
		closedir(dir);
	}
	
		
	save_strand_collection_stats(stats, path);
	
	free(path);
		
		
		
	/* Copy parameter files to new directory for record of what has been performed */
		
	stage_num = copy_parameters(input_dir_path, output_dir_path);
	
	if (stage_num != -1) {
		path = (char*)malloc(sizeof(char) * (strlen(output_dir_path) + 200));	
			
		sprintf(path, "%s%cparameters%c%d_mri_sim_param.txt", output_dir_path, DIR_SEP, DIR_SEP, stage_num);
		
		
		if ( (file =fopen(path,"w")) == NULL) {
			printf("\nCould not open file %s.!\n", path);
			exit(EXIT_SUCCESS);
		} 
										
																			
		fprintf(file, "fa %g\n", fa);
		fprintf(file, "diffusivity %g\n", diffusivity);
		fprintf(file, "num_voxels %d\n", num_voxels);	
		fprintf(file, "voxel_size %g\n", voxel_size);	
		fprintf(file, "subvoxel_size %g\n", subvoxel_size);
		fprintf(file, "save_ext_segment_stats %d\n", save_ext_segment_stats);
		fprintf(file, "save_ext_fill_radial_stats %d\n", save_ext_fill_radial_stats);
		fprintf(file, "save_subvoxels %d\n", save_subvoxels);
		fprintf(file, "sphere_r %g\n", sphere_r);	
		fprintf(file, "output_format %s\n", output_format);	

		
		if (fclose(file)) {
			printf("Error closing parameters file copy\n");
		}
		
		
		sprintf(path, "%s%cparameters%cgrad_directions.txt", output_dir_path, DIR_SEP, DIR_SEP);
		
		
		if ( (file =fopen(path,"w")) == NULL) {
			printf("\nCould not open gradient save file %s.!\n", path);
			exit(EXIT_SUCCESS);
		} 
		free(path);
										
																			
		for (grad_i = 0; grad_i < num_grad_directions; grad_i++) {
			fprintf(file, "%g %g %g %g\n", grad_directions[grad_i * 3 + X], grad_directions[grad_i * 3 + Y], grad_directions[grad_i * 3 + Z], b_values[grad_i]);
		}

		
		if (fclose(file)) {
			printf("Error closing gradient directions file copy\n");

		}
	}
	
	
	exit(EXIT_SUCCESS);

}


