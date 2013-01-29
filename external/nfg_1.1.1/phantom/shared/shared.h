/*
 *  shared.h
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

#ifndef SHARED_H
#define SHARED_H

#define VERSION_NUM "v1.11"


#ifdef __MINGW32__
	#define mkdir(a,b) mkdir(a)
	#define MY_PERMS 0
#else
     #define MY_PERMS (S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)
#endif


#define DIR_SEP '/'
#define DIR_SEP_S "/"

#define PI 3.14159265359

/* Define the order of cartesian coordinates in an array of 3 */
#define X 0
#define Y 1
#define Z 2

/* Define the order of spherical coordinates in an array of 3 */
#define EL 0		/* Elevation */
#define AZ 1		/*Azimuth */
#define R 2			/*Radius */

/* Used when loading strands from file and determining when the deviation from the specified sphere radius is not due to rounding errors. */
#define SPHERE_R_TOL 0.005
#define STRAND_EXT_FRACT 0.01

/* Used in writing Analyze headers */
#define ANALYZE_FLOAT 16
#define ANALYZE_INT 8
#define ANALYZE_UCHAR 2

double dist_between_points(double *p1, double *p2);

double angle_between_points(double *start, double *middle, double *end);

void sphere_intersect(double *intersect, double *p, double *v, double sphere_r);

double dot_product(double *v1, double *v2);

double* cross_product(double *v1, double *v2);

double min_dble(double x, double y);

double max_dble(double x, double y);

int min_int(int x, int y);

int max_int(int x, int y);

double vector_norm(double *v);

void normalize(double *v);

void cart_to_spher(double spher[3], double cart[3]);
 
void spher_to_cart(double cart[3], double spher[3]);

int num_decimal_places(int number);

double radius_of_curvature(double *start, double *middle, double *end);

int copy_parameters(char *src_dirname, char *dest_dirname);

int file_copy(char *src_path, char *dest_dir);

void write_analyze_header(int num_voxels[3], double voxel_size[3], int analyze_data_type, int data_length, char *output_path);

void read_analyze_header(int num_voxels[3], double voxel_size[3], char *input_path);

int compare_int (const void *a, const void *b);

#endif
