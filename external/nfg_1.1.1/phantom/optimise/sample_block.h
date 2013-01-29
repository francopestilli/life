/*
 *  sample_block.h
 *  Numerical Fibre Generator
 *
 *  The sample block structure is one item in a linked list of sample blocks used to 
 *  dynamically allocate memory to newly generated samples as required.
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

#ifndef SAMPLE_BLOCK_H
#define SAMPLE_BLOCK_H 0

#include "phantom/optimise/sample.h"
#include "phantom/shared/segment.h"

#define SAMPLE_BLOCK_SIZE 68

typedef struct _sample_block {
	
	//The index of the next free sample in the 'sample' array.
	int sample_count;
	
	//Used when iterating through the existing values.
	int next_i;
	
	//Create an array to hold the sample coordinates, the cost_function gradient at those coordinates and its sample fraction (the fraction between its two generating points it lies) 
	Sample samples[SAMPLE_BLOCK_SIZE];
	
		
	int init_check;
	int next_block_open;
	
	//Reference to the next sample block in the linked list (will be NULL if the sample block is the last in the list).
	struct _sample_block *next_block;

	
} Sample_block;

void reset_sample_block(Sample_block *sample_block);

void init_sample_block(Sample_block *sample_block);

Sample* new_Sample(Sample_block **sample_block, Segment *segment, int sample_i);

Sample* next_sample(Sample_block **sample_block);

Sample* first_sample(Sample_block **sample_block);

void free_sample_blocks(Sample_block *start_block);

#endif
