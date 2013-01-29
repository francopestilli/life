/*
 *  segment.h
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



#ifndef _SEGMENT_H
#define _SEGMENT_H


#include "phantom/shared/control_point.h"

typedef struct _segment {
	
	int segment_i;
	struct _strand *strand;
	Control_point start_point, end_point;
		
	int num_samples;
	double length;
	double sample_length;
	double disp[3];
	double sample_length_grad[3];
	
	struct _segment *prev_segment, *next_segment;
	
} Segment;


#ifndef _STRAND_H
#include "phantom/shared/strand.h"
#endif

#endif


#ifndef SEGMENT_H
#define SEGMENT_H

#include "phantom/shared/strand.h"

void init_segment(Segment *segment, int segment_i, Segment *prev_segment, Control_point point, Control_point next_point, Strand *strand, double sample_density); 

void print_segment(Segment *segment, char indent[]);

#endif
