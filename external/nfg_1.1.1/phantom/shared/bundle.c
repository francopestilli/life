/*
 *  bundle.c
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "phantom/shared/strand.h"
#include "phantom/shared/segment.h"
#include "phantom/shared/shared.h"
#include "phantom/shared/bundle.h"

#define BUNDLE_STRAND_ALLOC_BLOCK 10


int construct_bundles(Bundle *bundles, int *bundle_i_of_strand, int num_strands, Strand *strands) {
	
	Strand *strand;
	int strand_i, prev_strand_i, i, num_bundles, is_new;
	Bundle *bundle;
	
	
	num_bundles = 0;
	
	for (strand_i = 0; strand_i < num_strands; strand_i++) {

		strand = &(strands[strand_i]);

		is_new = 1;

		for (prev_strand_i = 0; prev_strand_i < strand_i; prev_strand_i++) {
			if (strand->bundle_i == bundle_i_of_strand[prev_strand_i])
				is_new = 0;
		}
		
		if (is_new) {
//			c->bundles_in_collection[num_bundles] = strand->bundle_i;

			bundle = &(bundles[num_bundles]);
			
			init_bundle(bundle, strand);	
				
			num_bundles++;
			
		} else {
			
			for (i = 0; i < num_bundles; i++) {
			
				bundle = &(bundles[i]);
				
				if (bundle->bundle_i == strand->bundle_i) {
					add_strand_to_bundle(bundle, strand);
					break;
				}
			
			}
			
		}
	}


	qsort(bundles, num_bundles, sizeof(Bundle), compare_bundle_indices);

	return num_bundles;

}

void init_bundle(Bundle *bundle, Strand *strand) {

	bundle->bundle_i = strand->bundle_i;
	
	bundle->strands = (Strand**)malloc(sizeof(Strand*) * BUNDLE_STRAND_ALLOC_BLOCK);
	
	bundle->strands[0] = strand;
	strand->bundle = bundle;
	
	bundle->num_strands = 1;
}


void free_bundle(Bundle *bundle) {

	free(bundle->strands);

}

void add_strand_to_bundle(Bundle *bundle, Strand *strand) {

	Strand **temp_strands;

	bundle->strands[bundle->num_strands] = strand;
	strand->bundle = bundle;
	
	bundle->num_strands++;
	
	if ((bundle->num_strands % BUNDLE_STRAND_ALLOC_BLOCK) == 0) {
		temp_strands = (Strand**)malloc(sizeof(Strand*) * bundle->num_strands);
		memcpy(temp_strands, bundle->strands, (sizeof(Strand*) * bundle->num_strands));

		bundle->strands = (Strand**)realloc(bundle->strands, (sizeof(Strand*) * (bundle->num_strands + BUNDLE_STRAND_ALLOC_BLOCK)));
		memcpy(bundle->strands, temp_strands, (sizeof(Strand*) * bundle->num_strands));
		free(temp_strands);
	}

}





int compare_bundle_indices(const void *a, const void *b) {

	return ( ((Bundle*)a)->bundle_i - ((Bundle*)b)->bundle_i);


}


double calculate_cross_sectional_area(Bundle *bundle) {

	Strand *strand;
	int strand_i;
	double cross_sectional_area;

	cross_sectional_area = 0;

	for (strand_i = 0; strand_i < bundle->num_strands; strand_i++) {
	
		strand = bundle->strands[strand_i];
		
		cross_sectional_area += strand->radius * strand->radius * PI;
	
	}

	return cross_sectional_area;

}

double calculate_average_strand_length(Bundle *bundle) {

	Strand *strand;
	Segment *segment;
	int strand_i;
	double average_strand_length, strand_length;

	average_strand_length = 0;

	for (strand_i = 0; strand_i < bundle->num_strands; strand_i++) {
	
		strand = bundle->strands[strand_i];
		
		strand_length = 0;
		
		segment = strand->start_segment;
		
		while (segment != strand->post_segment) {
			
			strand_length += vector_norm(segment->disp);
			segment = segment->next_segment;		
		}
		

		average_strand_length += strand_length;

	}

	average_strand_length /= bundle->num_strands;
	
	return average_strand_length;


}


int compare_bundle_cross_sectional_area(const void *a, const void *b) {

	double a_cross_section, b_cross_section;
	
	a_cross_section = calculate_cross_sectional_area((Bundle*)a);
	b_cross_section = calculate_cross_sectional_area((Bundle*)b);
	
	return a_cross_section - b_cross_section;


}

