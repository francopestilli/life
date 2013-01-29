/*
 *  isotropic_region.c
 *  Numerical Fibre Generator
 *
 *  Created by Tom Close on 11/01/09.
 *  Copyright 2009 Tom Close.
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

#include <math.h>
#include <stdio.h>


#include "phantom/shared/segment.h"
#include "phantom/shared/strand.h"
#include "phantom/shared/shared.h"

#include "phantom/shared/isotropic_region.h"


void init_isotropic_region(Isotropic_region *isotropic_region, double pos[3], double radius, double diffusivity, double baseline_signal, double weighting){ 
	
	/* This function may be a bit confusing to follow as it is utilizing the structures used to store the samples of the strand segments to hold the isotropic_region regions.  This is so the isotropic_region regions can be preloaded as samples into the cost function calculation.  It only has set up the fields that are required in calculating the sample cost but they are stored within 'strand' and 'segment' structs to fit with the legacy code.  */
	
	/* Set the radius of the isotropic_region region */
	isotropic_region->strand.radius = radius;
	
	/* Set the isotropic_region struct's strand as the isotropic_region struct's segment's strand */ 
	isotropic_region->segment.strand = &(isotropic_region->strand); 
	
	/* Provide a more logical and direct way of accessing this the isotropic_region's radius */
	isotropic_region->radius = radius;
	
	
	/* The overlap penalty for a strand passing through a isotropic_region region can be arbitrarily weighted within the legacy framework by utilizing the weighting on the segment length each sample is supposed to represent.  By setting the sample length to be longer the weighting for the isotropic_region region, which has no length, can be increased. Should be set to the inverse of the sample_density parameter for equal weighting between strands and isotropic_region regions. */

	isotropic_region->segment.sample_length_grad[X] = 0.0;
	isotropic_region->segment.sample_length_grad[Y] = 0.0;
	isotropic_region->segment.sample_length_grad[Z] = 0.0;
	isotropic_region->weighting = weighting;
	
	
	/* Note that the gradients of the isotropic_region regions are not used. However since the sample fraction is 0 the end_point gradient will never be set to the control_point_grad will contain the gradient value of the isotropic_region region if it could move.*/
	isotropic_region->segment.start_point.grad = isotropic_region->control_point_grad;
	isotropic_region->segment.end_point.grad = isotropic_region->control_point_grad;
	isotropic_region->segment.num_samples = 1;

	/* Initialise the segment start and end point positions to any old array. As these positions will be accessed in the 'new_Sample' function, when the sample is generated from the isotropic region, to generate the sample's position before this position is subsequently overwritten*/
	isotropic_region->segment.start_point.pos = isotropic_region->control_point_grad;
	isotropic_region->segment.end_point.pos = isotropic_region->control_point_grad;
	
	isotropic_region->diffusivity = diffusivity;
	isotropic_region->baseline_signal = baseline_signal;
	
	/* Set the isotropic_region's position */
	
	isotropic_region->pos[X] = pos[X];
	isotropic_region->pos[Y] = pos[Y];
	isotropic_region->pos[Z] = pos[Z];
	
	
}

void set_isotropic_regions_weightings(Strand_collection *c, double relative_isotropic_repulsion_weight) {

	int isotropic_region_i;
	
	for (isotropic_region_i = 0; isotropic_region_i < c->num_isotropic_regions; isotropic_region_i++) {
		c->isotropic_regions[isotropic_region_i].segment.sample_length = c->isotropic_regions[isotropic_region_i].weighting * relative_isotropic_repulsion_weight;
	}

}

