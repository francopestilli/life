/*
 *  segment.c
 *  Numerical Fibre Generator
 *
 *  Created by Tom Close, on 25/06/08.
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


#include <math.h>
#include <stdio.h>
#include <time.h>

#include "phantom/shared/segment.h"
#include "phantom/shared/shared.h"

#define COINCIDE_DISP 1e-4

void init_segment(Segment *segment, int segment_i, Segment *prev_segment, Control_point point, Control_point next_point, Strand *strand, double sample_density) { 	
	
	double grad_scalar;
	
	
	segment->disp[X] =  next_point.pos[X] - point.pos[X];
	segment->disp[Y] =  next_point.pos[Y] - point.pos[Y];
	segment->disp[Z] =  next_point.pos[Z] - point.pos[Z];
	
	if ( (segment->disp[X] == 0.0) & (segment->disp[X] == 0.0) & (segment->disp[X] == 0.0)) {
	
		segment->disp[X] = COINCIDE_DISP;
		next_point.pos[X] += COINCIDE_DISP;
		
	}
	
	
	segment->segment_i = segment_i;

	segment->length = sqrt(segment->disp[X] * segment->disp[X] + segment->disp[Y] * segment->disp[Y] + segment->disp[Z] * segment->disp[Z]); 
	
	
	segment->num_samples = ceil(segment->length * sample_density);

	/* differs from the inverse of 'sample_density' due to the ceiling num_samples to an int. */	
	segment->sample_length = segment->length/((double)segment->num_samples);				

	if (segment->length !=0) {
		grad_scalar = 1/(segment->length * (double)segment->num_samples);
	} else {
		grad_scalar = 0;
	}
	
	segment->sample_length_grad[X] = -1.0 * segment->disp[X] * grad_scalar;
	segment->sample_length_grad[Y] = -1.0 * segment->disp[Y] * grad_scalar;
	segment->sample_length_grad[Z] = -1.0 * segment->disp[Z] * grad_scalar;

	
	segment->strand = strand;
	
	segment->start_point = point;
	segment->end_point = next_point;
	
	segment->prev_segment = prev_segment;
	segment->next_segment = NULL;
	
	if (segment->prev_segment != NULL) {
		segment->prev_segment->next_segment = segment;
	}
	
	if (segment_i==0) {
		strand->start_segment = segment;
	}
	
	if (segment_i == strand->num_control_points) {
		strand->end_segment = segment;
	}
	
	
}




void print_segment(Segment *segment, char indent[]) {
	
	printf("\n%s-------- Start Segment --------\n\n", indent);
	printf("%sSegment %d:\n", indent, segment->segment_i);
	
	printf("\n%sStart Control_point:\n", indent);
	print_control_point(&segment->start_point, indent);
	printf("\n%sEnd Control_point:\n   ", indent);
	print_control_point(&segment->end_point, indent);	
	printf("\n%sNum samples %d, length %g, sample_length %g\n", indent, segment->num_samples, segment->length, segment->sample_length); 
	printf("\n%sDisp: [%g,%g,%g]", indent, segment->disp[X], segment->disp[Y], segment->disp[Z]);
	printf("\n%s-------- End Segment --------\n\n", indent);
}


