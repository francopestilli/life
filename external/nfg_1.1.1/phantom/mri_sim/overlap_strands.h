/*
 *  overlap_strands.h
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


#ifndef OVERLAP_STRANDS_H
#define OVERLAP_STRANDS_H

#include "phantom/shared/segment.h"
#include "phantom/shared/strand.h"

typedef struct _fill_strand {
	
	Strand *strand;
	Segment *closest_segment;
	double closest_fraction;
	
	struct _fill_strand *next;
	
} Overlap_strand;
void fill_strand_alloc(Overlap_strand **fill_strand, Strand *strand);

void add_fill_segment(Overlap_strand **overlap_strands, double fraction, Segment *segment);
#endif

