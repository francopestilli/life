/*
 *  strand.h
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







#ifndef _STRAND_H
#define _STRAND_H


typedef struct _strand {
	
	int strand_i;
	int bundle_i;

	int num_control_points;
	double radius;	
	
	double natural_seg_length;
	
	struct _segment *start_segment;
	struct _segment *end_segment;
	
	struct _segment *pre_segment;
	struct _segment *post_segment;
	
	struct _bundle *bundle;
	
} Strand;



#endif


#ifndef STRAND_H
#define STRAND_H

#include "phantom/shared/segment.h"
#include "phantom/shared/control_point.h"
#include "phantom/shared/shared.h"


void init_strand(Strand *strand, int strand_i, int bundle_i, int num_control_points, double radius, double *start_point, double *end_point);

void construct_strand(Strand *strand, int strand_i, int bundle_i, double *control_point, double *start_point, double *end_point, double *pre_point, double *post_point, Segment segments[], int num_strand_control_points, double sample_density, double strand_r);


#endif

