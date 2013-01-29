/*
 *  subvoxel.c
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

#include "phantom/mri_sim/subvoxel.h"
#include "phantom/shared/shared.h"

void init_subvoxel(Subvoxel *subvoxel, double *orientation) {
		
	subvoxel->orientation = orientation;
	
	subvoxel->orientation[X] = 0.0;
	subvoxel->orientation[Y] = 0.0;
	subvoxel->orientation[Z] = 0.0;
	
	subvoxel->closest_fraction = 1.0;
	subvoxel->closest_segment = NULL;
	subvoxel->closest_isotropic_region = NULL;
	
	subvoxel->overlap_strands = NULL;
	
}



void subvoxel_free(Subvoxel *subvoxel) {
		
	Overlap_strand *current, *prev;

	current = subvoxel->overlap_strands;
	
	while (current != NULL) {
		prev = current;
		current = current->next;
		free(prev);
	}
	
}


void add_orientation(Subvoxel *subvoxel, double orientation[3], double fraction, Segment *segment) {
	
	
	if (fraction <= subvoxel->closest_fraction) {
	
		subvoxel->orientation[X] = orientation[X];
		subvoxel->orientation[Y] = orientation[Y];
		subvoxel->orientation[Z] = orientation[Z];
		
		subvoxel->closest_fraction = fraction;
		subvoxel->closest_segment = segment;
		
	}
	
	
	add_fill_segment(&(subvoxel->overlap_strands), fraction, segment);
	
}

void add_isotropic_region(Subvoxel *subvoxel, double fraction, Isotropic_region *isotropic_region) {


	if (fraction <= subvoxel->closest_fraction) {
	
		subvoxel->orientation[X] = -99;
		subvoxel->orientation[Y] = -99;
		subvoxel->orientation[Z] = -99;
	
		subvoxel->closest_fraction = fraction;
		subvoxel->closest_isotropic_region = isotropic_region;
	
	}

}
