/*
 *  isotropic_region.h
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

#ifndef _ISOTROPIC_REGION_H
#define _ISOTROPIC_REGION_H

#include "phantom/shared/segment.h"
#include "phantom/shared/strand.h"

typedef struct _isotropic_region {
	
	Strand strand;
	Segment segment;
	double control_point_grad[3];
	
	double diffusivity;
	double baseline_signal;  //Change in signal baseline due to differring T2 decay. Note that diffusion signal baseline for strands is 1.
	double radius;
	double weighting;
		
	double pos[3];
		
} Isotropic_region;

#include "phantom/shared/strand_collection.h"

#endif

#ifndef ISOTROPIC_REGION_H
#define ISOTROPIC_REGION_H

void init_isotropic_region(Isotropic_region *isotropic_region, double pos[3], double radius, double diffusivity, double baseline_signal, double weighting);

void set_isotropic_regions_weightings(Strand_collection *c, double relative_isotropic_repulsion_weight);




#endif
