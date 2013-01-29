/*
 *  subvoxel.h
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
 *  it under the term\s of the GNU General Public License as published by
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

#ifndef SUBVOXEL_H
#define SUBVOXEL_H

#include "phantom/shared/segment.h"
#include "phantom/mri_sim/overlap_strands.h"
#include "phantom/shared/isotropic_region.h"


typedef struct _subvoxel {

	double *orientation;
	double closest_fraction;
	Segment *closest_segment;
	Isotropic_region *closest_isotropic_region;
	
	Overlap_strand *overlap_strands;

} Subvoxel;


void init_subvoxel(Subvoxel *subvoxel, double *orientation);

void subvoxel_free(Subvoxel *subvoxel);

void add_orientation(Subvoxel *subvoxel, double orientation[3], double fraction, Segment *segment);

void add_isotropic_region(Subvoxel *subvoxel, double fraction, Isotropic_region *isotropic_region);

#endif
