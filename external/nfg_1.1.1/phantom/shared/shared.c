/*
 *  shared.c
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
#include <sys/stat.h>
#include <dirent.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "phantom/shared/shared.h"






/* start_or_end_intersect is either END_INTERSECT (1) or START_INTERSECT (-1) and determines which of the two roots of mu is desired. */
void sphere_intersect(double *intersect, double *p, double *v, double sphere_r) {	
	
	/* NB p - point, v - vector. Finds the intersection along the line extension of 'v' about 'p' that intersects with a sphere centred on the origin with radius 'sphere_r' */
	
	double mu, a, b, c;
	int sign;

	
	a = dot_product(v, v);
	b = 2 * dot_product(p, v);
	c = dot_product(p, p) - sphere_r * sphere_r;
	
	sign = (b > 0) - (b < 0);
	
	mu = (-b + (double)sign * sqrt(b * b - 4 * a * c))/ (2 * a);
	
	intersect[X] = p[X] + mu * v[X];
	intersect[Y] = p[Y] + mu * v[Y];
	intersect[Z] = p[Z] + mu * v[Z];

}  

double dot_product(double *v1, double *v2) {
	double product;
	product = v1[X] * v2[X] + v1[Y] * v2[Y] + v1[Z] * v2[Z];
	return product;
}

double* cross_product(double *v1, double *v2) {

	double *output;

	output = (double*)malloc(sizeof(double) * 3 );

	output[X] = v1[Y] * v2[Z] - v1[Z] * v2[Y];
	output[Y] = v1[Z] * v2[X] - v1[X] * v2[Z];
	output[Z] = v1[X] * v2[Y] - v1[Y] * v2[X];
	
	return output;
}

double vector_norm(double *v) {

	return sqrt(v[X] * v[X] + v[Y] * v[Y] + v[Z] * v[Z]);

}

void normalize(double *v) {
	
	double norm;
	
	norm = vector_norm(v);
	
	v[X] = v[X] / norm;
	v[Y] = v[Y] / norm;
	v[Z] = v[Z] / norm;
	
}

double max_dble(double x, double y) {
	if (x > y) { 
		return x;
	} else {
		return y;
	}
}

double min_dble(double x, double y) {
	if (x < y) { 
		return x;
	} else {
		return y;
	}
}

int max_int(int x, int y) {
	if (x > y) { 
		return x;
	} else {
		return y;
	}
}

int min_int(int x, int y) {
	if (x < y) { 
		return x;
	} else {
		return y;
	}
}

void cart_to_spher(double spher[3], double cart[3]) {
	
	spher[R] = sqrt(cart[X] * cart[X] + cart[Y] * cart[Y] + cart[Z] * cart[Z]);
	
	if (spher[R] != 0.0) {
		spher[EL] = acos(cart[Z] / spher[R]);
		spher[AZ] = atan2(cart[Y], cart[X]);  
	} else {
		spher[EL] = 0.0;
		spher[AZ] = 0.0;
	}
}
 
void spher_to_cart(double cart[3], double spher[3]) {
	
	cart[X] = spher[R] * cos(spher[AZ]) * sin(spher[EL]);
	cart[Y] = spher[R] * sin(spher[AZ]) * sin(spher[EL]);
	cart[Z] = spher[R] * cos(spher[EL]);

}

int num_decimal_places(int number) {
	int num_places, ten_power;
	
	num_places = 1;
	ten_power = 10;
	while ((number / ten_power) != 0) {
		num_places++;
		ten_power *= 10;
	} 
	return num_places;
}

double radius_of_curvature(double *start, double *middle, double *end) {

/*	double radius, secant_a[3], secant_b[3], *out_of_plane, *in_plane_perpen_a, centre_disp[3], midpoint_a[3],midpoint_b[3], midpoint_disp[3], alpha;*/
	
	double radius, seg_disp[3];
	double chord;
	double chord_disp[3];
	double height;
	double height_disp[3];    
	double mu;
	
	
	/* The radius of curvature is approximated by calculating the radius of an arc whose base is the 'chord' between the start and end points and its height
	* is the shortest distance from this chord to the corner of the path*/ 
	
	/*The first segment displacement*/
	seg_disp[X] = middle[X] - start[X];
	seg_disp[Y] = middle[Y] - start[Y];
	seg_disp[Z] = middle[Z] - start[Z];		
	
	/* The chord which forms the base of the final arc */ 
	chord_disp[X] = end[X] - start[X];
	chord_disp[Y] = end[Y] - start[Y];
	chord_disp[Z] = end[Z] - start[Z];
	
	/*its length*/
	chord = vector_norm(chord_disp);
	
	/* Normalise the chord displacemnt */
	chord_disp[X] = chord_disp[X]/ chord;
	chord_disp[Y] = chord_disp[Y]/ chord;
	chord_disp[Z] = chord_disp[Z]/ chord;
	
	
	/*The projection distance of the first segment onto the direction of chord*/
	mu = dot_product(seg_disp, chord_disp); 

	/* The shortest displacement from the corner to the chord.*/ 
	height_disp[X] = seg_disp[X] - chord_disp[X] * mu;
	height_disp[Y] = seg_disp[Y] - chord_disp[Y] * mu;	
	height_disp[Z] = seg_disp[Z] - chord_disp[Z] * mu;	
	
	/*The 'height' of the arc.*/
	height = vector_norm(height_disp);
	
	/*From the intersecting chord theorem, the radius of the approximated arc*/
	radius = height / 2 + chord * chord / ( 8 * height );
	
	
	return radius;
	

}

double dist_between_points(double *point_1, double *point_2) {

	return sqrt((point_1[X] - point_2[X]) * (point_1[X] - point_2[X]) + (point_1[Y] - point_2[Y]) * (point_1[Y] - point_2[Y]) + (point_1[Z] - point_2[Z]) * (point_1[Z] - point_2[Z])); 
		

}

double angle_between_points(double *start, double *middle, double *end) {

	double angle, dp, v1[3], v2[3];
	
	v1[X] = middle[X] - start[X];
	v1[Y] = middle[Y] - start[Y];
	v1[Z] = middle[Z] - start[Z];
	
	v2[X] = end[X] - middle[X];
	v2[Y] = end[Y] - middle[Y];
	v2[Z] = end[Z] - middle[Z];
	
	normalize(v1);
	normalize(v2);
	
	dp = dot_product(v1,v2);
	
	
	/* Sometimes due to rounding errors the dot product can exceed 1.0.  Guard against this case. */
	if (dp >= 1.0) {	
		angle = 0;
	} else {
		angle = acos(dp);
	}
	
	angle = angle * 180/PI;

	
	return angle;
	

}

int copy_parameters(char *src_dirname, char *dest_dirname) {

	DIR *src_dir, *dest_dir;
	struct dirent *directory_entry;
	
	char *src_param_dirname, *dest_param_dirname, *path, filename[100], *appendage;
	int stage_num, error;
	
	stage_num = 0;
	
	
	dest_param_dirname = (char*)malloc(sizeof(char) * strlen(dest_dirname) + 100);

	
	sprintf(dest_param_dirname, "%s%cparameters", dest_dirname, DIR_SEP);
	
	if ((dest_dir = opendir(dest_param_dirname)) == NULL) {
		if ((error = mkdir(dest_param_dirname, MY_PERMS)) > 0) {
			printf("Could not create parameters directory %s, (Error code: %d)!\n", dest_param_dirname, error);
			return -1;
		}
	} else {
		closedir(dest_dir);
	}
	
	
	
	src_param_dirname = (char*)malloc(strlen(src_dirname) + 100);

	sprintf(src_param_dirname, "%s%cparameters", src_dirname, DIR_SEP);
	path = (char*)malloc(strlen(src_dirname) + 200);
	
	if ((src_dir = opendir(src_param_dirname)) != NULL) {		

		while ( (directory_entry = readdir(src_dir)) != NULL) {
		
			strcpy(filename, directory_entry->d_name);
			
			appendage = strrchr(filename, '_');
			
			if ((appendage != NULL) && (strcmp(appendage, "_param.txt") == 0)) {

				sprintf(path, "%s%c%s", src_param_dirname, DIR_SEP, filename);
				
				file_copy(path, dest_param_dirname);

				stage_num++;
				
			} else if (strcmp(filename, "grad_directions.txt") == 0) {
			
				sprintf(path, "%s%c%s", src_param_dirname, DIR_SEP, filename);
				
				file_copy(path, dest_param_dirname);
			
			}
			
			
		}

	
		free(dest_param_dirname);
		free(path);
	

	} else {
		printf("\nNo parameters directory found in input directory %s!\n", src_param_dirname);
		free(src_param_dirname);
	
		return 0;
	}

	free(src_param_dirname);
	
	return stage_num+1;
}

int file_copy(char *src_path, char *dest_dir)  {
	
	FILE *src, *dest;
	char *dest_path, ch, *src_filename;

	if ( (src_filename = strrchr(src_path, DIR_SEP) ) == NULL) {
		printf("Source path is not complete (i.e. Could not find '/' ('\' on windows) in the path)!, source path = %s\n", src_path);
		return 1;
	}


	dest_path = (char*)malloc(sizeof(char) * (strlen(dest_dir) + strlen(src_filename) + 10 ));	
	
	strcpy(dest_path, dest_dir);
	strcat(dest_path, src_filename);

	/* open source file */
	if((src = fopen(src_path, "rb"))==NULL) {
		printf("Cannot open source file, %s.\n", src_path);
		return 1;
	}

	/* open destination file */
	if((dest = fopen(dest_path, "wb"))==NULL) {
		printf("Cannot open destination file, %s.\n", dest_path);
		return 1;
	}

	/* copy the file */
	while(!feof(src)) {
		ch = fgetc(src);
		if(ferror(src)) {
			 printf("Error reading source file.\n");
			 return 1;
		}
		
		if(!feof(src)) fputc(ch, dest); {
			if(ferror(dest)) {
			 printf("Error writing destination file.\n");
			 return 1;
			}
		}
	}

	if(fclose(src)==EOF) {
		printf("Error closing source file.\n");
		return 1;
	}

	if(fclose(dest)==EOF) {
		printf("Error closing destination file.\n");
		return 1;
	}

	free(dest_path);

	return 0;
}

void read_analyze_header(int num_voxels[3], double voxel_size[3], char *input_path) {
	FILE *file;
	int dim_i;
	char header[348];


	if ( (file =fopen(input_path,"rb")) == NULL) {
		printf("\nCould not open file %s.!\n", input_path);
		exit(EXIT_FAILURE);
	} 


        if (!fread(header, sizeof(char), 348, file)) {
           printf("\nNo elements were read from file %s.!\n", input_path);
           exit(EXIT_FAILURE);
        }
	
	fclose(file);
	
	for (dim_i = 0; dim_i < 3; dim_i++) {
		num_voxels[dim_i] = *((int*) (header + 42 + 2 * dim_i) );
		voxel_size[dim_i] = *((double*) (header + 80 + 4 * dim_i) );
	}

	return;
}



void write_analyze_header(int num_voxels[3], double voxel_size[3], int analyze_data_type, int data_length, char *output_path) {
/* data_type can be either 'f' for float or 'd' for int */

	FILE *file;
	int dim_i;
	char header[348];
	
	/* initialise all characters to 0 */
	memset(header, 0, 348);
	
	/* Required opening of an analyze file. */
	*((int*) header) = 348;
	
	/* Datatype of file. */
	strncpy(header + 4, "dsr      \0", 10);

	/* Database name */
	strncpy(header + 14, "randomly generated\0", 18);
	
	/* Unknown yet required cryptic number */
	*((int*) (header + 32) ) = 16384;
	
	/* Unknown and maybe not even required. */
	strncpy(header + 38, "r0", 2);
	
	/* Number of dimensions */
	*((short*) (header + 40) ) = 3;
	
	/*First 3 Dimension lengths */
	for (dim_i = 0; dim_i < 3; dim_i++) 
		*((short*) (header + 42 + 2 * dim_i) ) = num_voxels[dim_i];
	
	/* 4th dimension length */
	*((short*) (header + 48) ) = 1;
	
	/*Data type and length */
	*((short*) (header + 70) ) = analyze_data_type;
	*((short*) (header + 72) ) = data_length * 8;

	/* Voxel dimensions */
	for (dim_i = 0; dim_i < 3; dim_i++) 
		*((float*) (header + 80 + 4 * dim_i) ) = voxel_size[dim_i];
	
	/* Length of voxels in 4th dimension (not applicable in this case) */
	*((float*) (header + 92) ) = 1;
	
	/* Scale */
	*((float*) (header + 112) ) = 1.0;
	
	/* Offset */
	*((float*) (header + 116) ) = 0.0;
	
	/* Description */
	strncpy(header + 148, "Generated from NFG by Tom Close, Brain Research Institute, Melbourne, Australia", 80);
	
	strncpy(header + 228, "none                   \0", 24);

	/* Open the file and write the header array.*/
	
	if ( (file =fopen(output_path,"wb")) == NULL) {
		printf("\nCould not open file %s.!\n", output_path);
		exit(EXIT_FAILURE);
	} 


        if (!fwrite(header, sizeof(char), 348, file)) {
           printf("\nNo elements were written from file %s.!\n", output_path);
           exit(EXIT_FAILURE);
        }
	
	
	fclose(file);

	return;
}


int compare_int (const void *a, const void *b)
{
  int temp = *((int*)a) - *((int*)b);
  if (temp > 0)
    return 1;
  else if (temp < 0)
    return -1;
  else
    return 0;
}

