/*
 *  sample.c
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



#include <stdio.h>
#include <math.h>

#include "phantom/optimise/sample.h"
#include "phantom/shared/segment.h"

#define QUAD_COST 1

void init_sample(Sample *sample, int sample_i, Segment *segment) {
	
	double sample_fract;
	
	sample_fract = (double)sample_i/(double)(segment->num_samples);
	sample->sample_fract = sample_fract;
	
	sample->pos[X] = (segment->start_point).pos[X] * ( 1 - sample_fract) + (segment->end_point).pos[X] * sample_fract;
	sample->pos[Y] = (segment->start_point).pos[Y] * ( 1 - sample_fract) + (segment->end_point).pos[Y] * sample_fract;
	sample->pos[Z] = (segment->start_point).pos[Z] * ( 1 - sample_fract) + (segment->end_point).pos[Z] * sample_fract;
	
	sample->sample_i = sample_i; 
	sample->segment = segment;


}


void print_sample(Sample *sample) {
	
	printf("\n\n-------- Start Sample --------\n\n");
	printf("Strand: %d, Segment %d, Sample %d: \n", sample->segment->strand->strand_i, sample->segment->segment_i, sample->sample_i);
	
	print_segment(sample->segment, "\t");
	
	printf("Coord [%g, %g, %g]\n\n", sample->pos[X],sample->pos[Y], sample->pos[Z]);
	
	if (sample->last_accessed_sample != NULL ) {
		printf("Last Accessed sample [Strand %d, Seg %d, Sample %d]", sample->last_accessed_sample->segment->strand->strand_i, sample->last_accessed_sample->segment->segment_i, sample->last_accessed_sample->sample_i);
	}
	printf("\n\n-------- End Sample --------\n\n");
}

