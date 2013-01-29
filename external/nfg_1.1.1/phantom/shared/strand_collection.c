/*
 *  strand_collection.c
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
#include <string.h>
#include <sys/stat.h>
#include <dirent.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "phantom/shared/strand_collection.h"
#include "phantom/shared/isotropic_region.h"
#include "phantom/shared/shared.h"

#define NUM_CONTROL_POINTS_ON_STRAND_BLOCK_SIZE 100

void collection_alloc(Strand_collection *c, int num_strands, int num_control_points, int num_isotropic_regions) {
	
	c->num_strands = num_strands;
	c->num_control_points = num_control_points;
	c->num_isotropic_regions = num_isotropic_regions;
	c->sphere_r = -1.0;
	c->strand_r = (double*)malloc(sizeof(double) * num_strands);
	c->num_strand_control_points = (int*)malloc(sizeof(int) * num_strands);
	c->control_points = (double*)malloc(sizeof(double) * 3 * num_control_points);
	c->control_points_grad = (double*)calloc(num_control_points, sizeof(double));
	c->start_points = (double*)malloc(sizeof(double) * 3 * num_strands);
	c->end_points = (double*)malloc(sizeof(double) * 3 * num_strands);
	c->pre_points = (double*)malloc(sizeof(double) * 3 * num_strands);
	c->post_points = (double*)malloc(sizeof(double) * 3 * num_strands);	
	c->strands = (Strand*)malloc(sizeof(Strand) * num_strands);
	c->segments = (Segment*)malloc(sizeof(Segment) * (num_control_points + 3 * num_strands));
	c->isotropic_regions = (Isotropic_region*)malloc(sizeof(Isotropic_region) * num_isotropic_regions);
	c->bundle_i_of_strand = (int*)malloc(sizeof(int) * num_strands);
	
	//The maximum number of bundles is the number of strands (i.e. each strand is a seperate bundle).
	//c->bundles_in_collection = malloc(sizeof(int) * num_strands);
	c->bundles = (Bundle*)malloc(sizeof(Bundle) * num_strands);
	
}

void collection_free(Strand_collection *c) {

	int bundle_i;

	free(c->control_points);
	free(c->control_points_grad);	
	free(c->start_points);
	free(c->end_points);
	free(c->strand_r);
	free(c->num_strand_control_points);
	free(c->pre_points);	
	free(c->post_points);
	free(c->strands);
	free(c->segments);
	free(c->isotropic_regions);
	free(c->bundle_i_of_strand);
	
	for (bundle_i = 0; bundle_i < c->num_bundles; bundle_i++) {
		free_bundle(&(c->bundles[bundle_i]));
	}
	
	free(c->bundles);
	
}


int save_collection( 
				Strand_collection *c,
				char *output_dir_path,
				int save_start_index
			) {

	FILE *file;
	DIR *dir;
	int strand_i, point_i, error, end_of_strand_triple_i;
	char filename[100];
	char *dir_path, *path;
	int isotropic_region_i;
	int strand_i_max_width, max_bundle_i, bundle_i_max_width;

	point_i = 0;


	dir_path = (char*)malloc((sizeof(char) * strlen(output_dir_path) + 2));	
	strcpy(dir_path, output_dir_path);
	
	if (dir_path[strlen(output_dir_path) - 1] != DIR_SEP) {
		strcat(dir_path, DIR_SEP_S);			
	} 
	

	if ((dir = opendir(dir_path)) == NULL) {
		if ((error = mkdir(dir_path, MY_PERMS)) > 0) {
			printf("Error! Could not create directory %s, (Error code: %d)!\n", dir_path, error);
			exit(EXIT_FAILURE);
		}
	} else {
		closedir(dir);
	}	


	strand_i_max_width = num_decimal_places(c->num_strands);
	
	max_bundle_i = 0;	
	for (strand_i = 0; strand_i < c->num_strands; strand_i++) {
		if (c->bundle_i_of_strand[strand_i] > max_bundle_i) {
			max_bundle_i = c->bundle_i_of_strand[strand_i];
		}
	}
	
	bundle_i_max_width = num_decimal_places(max_bundle_i);
		
	/* Save isotropic_region regions */
	path = (char*)malloc(sizeof(char) * (strlen(dir_path) + strlen("isotropic_regions.txt") + 1));	
	strcpy(path,dir_path);
	strcat(path,"isotropic_regions.txt");
	
	if ( (file =fopen(path,"w")) == NULL) {
		printf("\nCould not open file %s.!\n", path);
		exit(EXIT_FAILURE);
	} 
	
	for (isotropic_region_i = 0; isotropic_region_i < c->num_isotropic_regions; isotropic_region_i++) {
		if (fprintf(file, "%10f\t%10f\t%10f\t%10f\t%10f\t%10f\t%10f\n", c->isotropic_regions[isotropic_region_i].pos[X], c->isotropic_regions[isotropic_region_i].pos[Y], c->isotropic_regions[isotropic_region_i].pos[Z], c->isotropic_regions[isotropic_region_i].radius, c->isotropic_regions[isotropic_region_i].diffusivity, c->isotropic_regions[isotropic_region_i].baseline_signal, c->isotropic_regions[isotropic_region_i].weighting) < 0) {
			printf("Error! Could not print isotropic_region region %d\n", isotropic_region_i);
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(file);
	free(path);

	/* Save strands */
	for (strand_i = 0; strand_i < c->num_strands; strand_i++) {
		
		sprintf(filename, "strand_%0*d-%0*d-r%f.txt", strand_i_max_width, (strand_i + save_start_index), bundle_i_max_width, c->bundle_i_of_strand[strand_i], c->strand_r[strand_i]);
			
		path = (char*)malloc(sizeof(char) * (strlen(dir_path) + strlen(filename) + 1));	
		strcpy(path, dir_path);
		strcat(path, filename);
		
		if ( (file =fopen(path,"w")) == NULL) {
			printf("\nCould not open file %s.!\n", path);
			exit(EXIT_FAILURE);
		} 

		if ( fprintf(file, "%10f\t%10f\t%10f\n", c->pre_points[strand_i * 3 + X], c->pre_points[strand_i * 3 + Y], c->pre_points[strand_i * 3 + Z]) < 0) {
			printf("Error! Could not print pre_point of strand %d\n", strand_i);
			exit(EXIT_FAILURE);
		}


		if ( fprintf(file, "%10f\t%10f\t%10f\n", c->start_points[strand_i * 3 + X], c->start_points[strand_i * 3 + Y], c->start_points[strand_i * 3 + Z]) < 0) {
			printf("Error! Could not print start_point of strand %d\n", strand_i);
			exit(EXIT_FAILURE);
		}

		end_of_strand_triple_i = point_i + c->num_strand_control_points[strand_i];

		for (; point_i < end_of_strand_triple_i; point_i++) {
		
			if ( fprintf(file, "%10f\t%10f\t%10f\n", c->control_points[point_i * 3 + X], c->control_points[point_i * 3 + Y], c->control_points[point_i * 3 + Z]) < 0 ) {
			printf("Error! Could not print control_point %d of strand %d\n", point_i, strand_i);			
				exit(EXIT_FAILURE);
			}
			
		}
		
		
		if ( fprintf(file, "%10f\t%10f\t%10f\n", c->end_points[strand_i * 3 + X], c->end_points[strand_i * 3 + Y], c->end_points[strand_i * 3 + Z]) < 0) { 
			printf("Error! Could not print end_point of strand %d\n", strand_i);		
			exit(EXIT_FAILURE);
		}

		if ( fprintf(file, "%10f\t%10f\t%10f\n", c->post_points[strand_i * 3 + X], c->post_points[strand_i * 3 + Y], c->post_points[strand_i * 3 + Z]) < 0) { 
			printf("Error! Could not print post_point of strand %d\n", strand_i);		
			exit(EXIT_FAILURE);
		}
		
		
		fclose(file);
		free(path);
	}
	
	free(dir_path);
	
			
	return EXIT_SUCCESS;
			
}


int load_collection(
				Strand_collection *c, 
				char *input_dir_path) {
	int error;
	DIR *directory;
	FILE *file;
	struct dirent *directory_entry;
	char *dir_path, *path;
	char filename[100];
	
	char line[100];
	
	int strand_i, saved_strand_i, bundle_i, start_of_strand_triple_i, end_of_strand_triple_i, point_i, dim_i;
	int num_strands, num_control_points, *num_strand_control_points, num_isotropic_regions;
	double strand_r, disp_from_origin; /*sphere_r, new_sphere_r, ;*/
	int  line_num; /*set_sphere_r,*/
	double control_point[3];

			
	int scans;	

		
	error = 0;
	
	dir_path = (char*)malloc((sizeof(char) * strlen(input_dir_path) + 2));	
	strcpy(dir_path, input_dir_path);
	
	if (dir_path[strlen(input_dir_path) - 1] != DIR_SEP) {
		strcat(dir_path, DIR_SEP_S);			
	} 
	
	printf("\nLoading collection, %s ...\n", dir_path);	
	fflush(stdout);
				
	if ((directory = opendir(dir_path)) == NULL) {
		printf("Error! Could not open directory %s!\n", dir_path);
		exit(EXIT_FAILURE);
	}
	
	num_strands = 0;
	num_control_points = 0;
	line_num = 0;
	
	c->fov = 0;
	
	num_isotropic_regions = 0;
	
	num_strand_control_points = (int*)malloc(sizeof(int) * NUM_CONTROL_POINTS_ON_STRAND_BLOCK_SIZE);
	
	while ( (directory_entry = readdir(directory)) != NULL) {
		
		strcpy(filename, directory_entry->d_name);
		
		if (strcmp(filename, "isotropic_regions.txt") == 0) {
		
			path = (char*)malloc(sizeof(char) * (strlen(dir_path) + strlen(filename) + 1));
			strcpy(path, dir_path);
			strcat(path,directory_entry->d_name);
			
			num_isotropic_regions = count_isotropic_regions(path);
			
		
		}
		
		
		if ( sscanf(filename, "strand_%d-%*d-r_%lf.txt", &saved_strand_i, &strand_r) > 0) {
			
			num_strand_control_points[num_strands] = 0;			
			
			path = (char*)malloc(sizeof(char) * (strlen(dir_path) + strlen(filename) + 1));
			strcpy(path, dir_path);
			strcat(path,directory_entry->d_name);
			
			if ((file = fopen(path, "r")) == NULL) {
				printf("Error! Could not open file %s\n", path);
				exit(EXIT_FAILURE);
			}
			
			line_num = 0;
			
			while (fgets(line, 100, file) != NULL) {
			
				if (sscanf(line, "%lf\t%lf\t%lf", &(control_point[X]), &(control_point[Y]), &(control_point[Z])) == 0) {
					printf("Error! Could not read line %d of strand %d\n", line_num, saved_strand_i);
					exit(EXIT_FAILURE);
				}
				
			
				for (dim_i = 0; dim_i < 3; dim_i++) {
					disp_from_origin = fabs(control_point[dim_i]) + strand_r;
					if (disp_from_origin > c->fov) {
						c->fov = disp_from_origin;
					}
				}
			
			
				num_strand_control_points[num_strands]++;
				num_control_points++;
				
				line_num++;				
			}
			
			fclose(file);
			
			/* Start, end, pre and post points are not counted. */
			
			if (num_strand_control_points[num_strands] >= 4) {
				
				num_control_points -= 4;
				num_strand_control_points[num_strands] -= 4;  
				
				num_strands++;
				if (num_strands % NUM_CONTROL_POINTS_ON_STRAND_BLOCK_SIZE == 0)
					num_strand_control_points = (int*)realloc(num_strand_control_points, sizeof(int) * (num_strands + NUM_CONTROL_POINTS_ON_STRAND_BLOCK_SIZE));
			
			} else {
				printf("Error! File %s did not have enough points (min 4) please rectify or remove the file!\n", filename);
				exit(EXIT_FAILURE);
			}
						
			
			
		}
	}
	
	
	if (num_strands == 0) {
		printf("\n!!! Warning no strands loaded from directory %s!\n\n", dir_path);	
		fflush(stdout);	
	}
	
	collection_alloc(c, num_strands, num_control_points, num_isotropic_regions);
				
	strand_i = 0;
	start_of_strand_triple_i = 0;
	point_i = 0;
	
	rewinddir(directory);
	
	while ( (directory_entry = readdir(directory)) != NULL) {
	
		strcpy(filename, directory_entry->d_name);
	
		if (strcmp(filename, "isotropic_regions.txt") == 0) {
			
			path = (char*)malloc(sizeof(char) * (strlen(dir_path) + strlen(filename) + 1));
			strcpy(path, dir_path);
			strcat(path,directory_entry->d_name);
			
			load_isotropic_regions(c->isotropic_regions, path, num_isotropic_regions);
		
		}
	
	
	
		if ((scans = sscanf(filename, "strand_%d-%d-r%lf.txt", &saved_strand_i, &bundle_i, &(c->strand_r[strand_i]))) > 0) {
		
			c->num_strand_control_points[strand_i] = num_strand_control_points[strand_i];
			
			path = (char*)malloc(sizeof(char) * (strlen(dir_path) + strlen(filename) + 1));
			strcpy(path, dir_path);
			strcat(path,directory_entry->d_name);
			
			if ((file = fopen(path, "r")) == NULL) {
				printf("Error! Could not open file %s\n", path);
				exit(EXIT_FAILURE);
			}
			
	
			/* Read pre_points */
			if (fgets(line, 100, file) == NULL) {
				printf("Error! Mismatch between number of points (first time %d, second time 0) in subsequent reads of file %s!\n", (c->num_strand_control_points[strand_i]), filename);
			
				exit(EXIT_FAILURE);
			} else {
				sscanf(line, "%lf\t%lf\t%lf", &(c->pre_points[strand_i * 3 + X]), &(c->pre_points[strand_i * 3 + Y]), &(c->pre_points[strand_i * 3 + Z]));
			}			
							
			/* Read start_points */
			if (fgets(line, 100, file) == NULL) {
				printf("Error! Mismatch between number of points (first time %d, second time 0) in subsequent reads of file %s!\n", (c->num_strand_control_points[strand_i]), filename);
			
				exit(EXIT_FAILURE);
			} else {
				sscanf(line, "%lf\t%lf\t%lf", &(c->start_points[strand_i * 3 + X]), &(c->start_points[strand_i * 3 + Y]), &(c->start_points[strand_i * 3 + Z]));
			}

			
			/* Read control_points */
			start_of_strand_triple_i = point_i;
			end_of_strand_triple_i = start_of_strand_triple_i + c->num_strand_control_points[strand_i];
			
			for (; point_i < end_of_strand_triple_i; point_i++) {
				
				if (fgets(line, 100, file) == NULL) {
					printf("Error! Mismatch between number of points (first time %d, second time %d) in subsequent reads of file %s!\n", (c->num_strand_control_points[strand_i]), (point_i - start_of_strand_triple_i), filename);
					exit(EXIT_FAILURE);
				} else {
					sscanf(line, "%lf\t%lf\t%lf", &(c->control_points[point_i * 3 + X]), &(c->control_points[point_i * 3 + Y]), &(c->control_points[point_i * 3 + Z]));
				}
			}
			
			/* Read end_points */
			if (fgets(line, 100, file) == NULL) {
				printf("Error! Mismatch between number of points (first time %d, second time %d) in subsequent reads of file %s!\n", (c->num_strand_control_points[strand_i]), (point_i - start_of_strand_triple_i), filename);
				exit(EXIT_FAILURE);
			} else {
				sscanf(line, "%lf\t%lf\t%lf", &(c->end_points[strand_i * 3 + X]), &(c->end_points[strand_i * 3 + Y]), &(c->end_points[strand_i * 3 + Z]));
			}
			 
			/* Read post_points */ 
			if (fgets(line, 100, file) == NULL) {
				printf("Error! Mismatch between number of points (first time %d, second time %d) in subsequent reads of file %s!\n", (c->num_strand_control_points[strand_i]), (point_i - start_of_strand_triple_i), filename);
				exit(EXIT_FAILURE);
			} else {
				sscanf(line, "%lf\t%lf\t%lf", &(c->post_points[strand_i * 3 + X]), &(c->post_points[strand_i * 3 + Y]), &(c->post_points[strand_i * 3 + Z]));
			}			 
						 
			construct_strand(&(c->strands[strand_i]), strand_i, bundle_i, &(c->control_points[start_of_strand_triple_i * 3]), &(c->start_points[strand_i * 3]), &(c->end_points[strand_i * 3]), &(c->pre_points[strand_i * 3]), &(c->post_points[strand_i * 3]), &(c->segments[start_of_strand_triple_i + 3 * strand_i]), c->num_strand_control_points[strand_i], 0.0, c->strand_r[strand_i]); 
			
			
			c->bundle_i_of_strand[strand_i] = bundle_i;
			

			fclose(file);
			
			
			free(path);
			
			strand_i++;
		}
		
	}
	

	// Determine the bundle indices that are present and collate them into a order list ('bundles_in_collection').
	


	c->num_bundles = construct_bundles(c->bundles, c->bundle_i_of_strand, c->num_strands, c->strands);

	
	printf("\n\nSuccessfully read %d strands and %d isotropic regions from directory %s\n\n", c->num_strands, c->num_isotropic_regions, input_dir_path);
	fflush(stdout);
	
	closedir(directory);
	
	free(dir_path);
	free(num_strand_control_points);
	
	return EXIT_SUCCESS;	

}


int count_isotropic_regions(char *path) {

	FILE *file;
	int line_num, num_isotropic_regions;
	char line[512];
	double isotropic_region_pos[3], isotropic_region_radius, isotropic_region_diffusivity, isotropic_region_baseline_signal, isotropic_region_weighting;

	
	if ((file = fopen(path, "r")) == NULL) {
		printf("Error! Could not open file %s\n", path);
		exit(EXIT_FAILURE);
	}
	
	line_num = 0;
	num_isotropic_regions = 0;
	
	while (fgets(line, 100, file) != NULL) {
	
		if (sscanf(line, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", &(isotropic_region_pos[X]), &(isotropic_region_pos[Y]), &(isotropic_region_pos[Z]), &isotropic_region_radius, &isotropic_region_diffusivity, &isotropic_region_baseline_signal, &isotropic_region_weighting) == 0) {
			printf("Error! Could not read line %d of isotropic_regions file\n", line_num);
			exit(EXIT_FAILURE);
		}
		
		num_isotropic_regions++;				
		line_num++;				
	}
	
	fclose(file);
	
	return num_isotropic_regions;

}

void load_isotropic_regions(Isotropic_region *destination, char *path, int num_isotropic_regions) {

	FILE *file;
	double isotropic_region_pos[3], isotropic_region_radius, isotropic_region_diffusivity, isotropic_region_baseline_signal, isotropic_region_weighting;
	char line[512];
	int isotropic_region_i;

	if ((file = fopen(path, "r")) == NULL) {
		printf("Error! Could not open file %s\n", path);
		exit(EXIT_FAILURE);
	}
	
	for (isotropic_region_i=0; isotropic_region_i < num_isotropic_regions; isotropic_region_i++) {
	
		if (fgets(line, 100, file) == NULL) {
			printf("Error! Mismatch between number of points (first time %d, second time 0) in subsequent reads of file %s!\n", num_isotropic_regions, path);
		
			exit(EXIT_FAILURE);
		} else {
			sscanf(line, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", &(isotropic_region_pos[X]), &(isotropic_region_pos[Y]), &(isotropic_region_pos[Z]), &isotropic_region_radius, &isotropic_region_diffusivity, &isotropic_region_baseline_signal, &isotropic_region_weighting);
			
			init_isotropic_region(&(destination[isotropic_region_i]), isotropic_region_pos, isotropic_region_radius, isotropic_region_diffusivity, isotropic_region_baseline_signal, isotropic_region_weighting);
		}	
		
	}
				
	fclose(file);


}

void copy_isotropic_regions(Strand_collection *destination, Isotropic_region *isotropic_regions, int num_isotropic_regions) {
	
	destination->isotropic_regions = (Isotropic_region*)realloc(destination->isotropic_regions, sizeof(Isotropic_region) * num_isotropic_regions);
	
	destination->num_isotropic_regions = num_isotropic_regions;

	memcpy(destination->isotropic_regions, isotropic_regions, sizeof(Isotropic_region) * num_isotropic_regions);


}




///* Gives the relative order of the given bundle index with respect to the other bundle indices in the collection.  
//Note that some bundle indices may not be present as they may have been trimmed.  In this case the denoted bundle 
//index does not correspond to the order of the bundle index in the collection.  This function finds the 'relative'
// bundle index, which is the rank of the given bundle index amongst the present bundle indices.  For example, if bundle
// indices 2 and 5 were missing then bundle index 8 would have a relative bundle index 6.*/
//int get_relative_bundle_i(Strand_collection *c, int bundle_i) {
//	
//	int  *bundle_i_address, relative_bundle_i;
//	
//	/* Finds the address of the given bundle index in the 'c->bundles_in_collection' array */
//	bundle_i_address = bsearch(&(bundle_i), c->bundles_in_collection, c->num_bundles, sizeof(int), compare_int);
//	
//	if (bundle_i_address == 0) {
//		return -1;
//	} 
//	
//	/* The offset of the bundle_i_address from the start of the 'c->bundles_in_collection' array */
//	relative_bundle_i = ((int)bundle_i_address) - ((int)c->bundles_in_collection);
//	
//	return relative_bundle_i;
//}


