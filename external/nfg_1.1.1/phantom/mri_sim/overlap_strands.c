/*
 *  overlap_strands.c
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

#include "phantom/mri_sim/overlap_strands.h"

void fill_strand_alloc(Overlap_strand **f, Strand *strand) {

	Overlap_strand *fill_strand;

	*f = (Overlap_strand*)malloc(sizeof(Overlap_strand));
	
	fill_strand = *f;

	fill_strand->strand = strand;
	fill_strand->closest_fraction = 1.0;
	fill_strand->next = NULL;

}

void add_fill_segment(Overlap_strand **overlap_strands, double fraction, Segment *segment) {

		Overlap_strand *current;

		if (*overlap_strands == NULL) {
		
			fill_strand_alloc(overlap_strands, segment->strand);
			current = *overlap_strands;
		
		} else {
						
			current = *overlap_strands;
										
			while (current->strand != segment->strand) {
			
				if (current->next == NULL) {
					fill_strand_alloc(&(current->next), segment->strand);				
				}
			
				current = current->next;
				
			}
		}
		
		if (fraction <= current->closest_fraction) {
		
			current->closest_segment = segment;
			current->closest_fraction = fraction;
			
		}
		
		
}

