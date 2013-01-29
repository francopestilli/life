/*
 *  sample.h
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
#include "phantom/shared/segment.h"
#endif

#ifndef SAMPLE_H
#define SAMPLE_H 0

typedef struct _sample {
	
	int sample_i;
	
	double sample_fract;
	Segment *segment;
	
	double pos[3];

	
	struct _sample *last_accessed_sample;		
	
} Sample;





void init_sample(Sample *sample, int sample_i, Segment *segment);

void print_sample(Sample *sample);

#endif
