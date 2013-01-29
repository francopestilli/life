/*
 *  strand_collection.h
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

#ifndef _ISOTROPIC_REGION_H
#include "phantom/shared/isotropic_region.h"
#endif

#ifndef _STRAND_COLLECTION_H
#define _STRAND_COLLECTION_H 0

#include "phantom/shared/strand.h"
#include "phantom/shared/bundle.h"

typedef struct _strand_collection {
	
	int num_strands;
	int num_control_points;
	int *num_strand_control_points;
	double *strand_r;
	
	double *control_points;
	double *control_points_grad;
	
	double *start_points;
	double *end_points;
	
	double *pre_points;	
	double *post_points;

	double sphere_r;
	double fov;
	
	Strand *strands;
	Segment *segments;
	Isotropic_region *isotropic_regions;
	
	int num_isotropic_regions;


	
	int num_bundles;
	Bundle *bundles;

	int *bundle_i_of_strand;

	/* Note that not all bundle indices will make it to the final collection as they may have been trimmed */
	//int *bundles_in_collection;
	
} Strand_collection;



void collection_alloc(Strand_collection *c, int num_strands, int num_control_points, int num_isotropic_regions);

void collection_free(Strand_collection *c);

int save_collection(Strand_collection *c, char *dir_path, int save_start_index);
			
int load_collection(Strand_collection *c, char *dir_path);

int count_isotropic_regions(char *path);

void load_isotropic_regions(Isotropic_region *destination, char *file_path, int num_isotropic_regions);

void copy_isotropic_regions(Strand_collection *dest_collection, Isotropic_region *isotropic_regions, int num_isotropic_regions);

#endif
