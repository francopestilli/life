/*
 *  bundle.h
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





#ifndef BUNDLE_H
#define BUNDLE_H

typedef struct _bundle {
	
	int bundle_i;
	int num_strands;
	struct _strand **strands;	
	
} Bundle;




int construct_bundles(Bundle *bundles, int *bundle_i_of_strand, int num_strands, Strand *strands);

void init_bundle(Bundle *bundle, Strand *strand);

void free_bundle(Bundle *bundle);

void add_strand_to_bundle(Bundle *bundle, Strand *strand);

int compare_bundle_indices(const void *a, const void *b);

double calculate_cross_sectional_area(Bundle *bundle);

double calculate_average_strand_length(Bundle *bundle);

int compare_bundle_cross_sectional_area(const void *a, const void *b);




#endif

