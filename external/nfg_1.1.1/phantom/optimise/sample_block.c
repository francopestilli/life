/*
 *  sample_block.c
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


#include <stdlib.h>
#include <stdio.h>

#include "phantom/optimise/sample.h"
#include "phantom/shared/segment.h"
#include "phantom/optimise/sample_block.h"



//------------------------------------------------------------------------//
// Initialize a newly allocated sample block.
//------------------------------------------------------------------------//


void reset_sample_block(Sample_block *sample_block) {
	if (sample_block->init_check != -7) {
	
		printf("Sample block was never initiated\n");
		exit(0);
	}
	
	sample_block->next_i = 0;
	sample_block->sample_count = 0;
	sample_block->next_block_open = 0;
}

void init_sample_block(Sample_block *sample_block) {



	sample_block->next_i = 0;
	sample_block->sample_count = 0;
	sample_block->next_block = NULL;
	sample_block->init_check = -7;
	sample_block->next_block_open = 0;	
}

//------------------------------------------------------------------------//
// Append a new sample to the end of the sample block
//------------------------------------------------------------------------//

Sample* new_Sample(Sample_block **sample_block, Segment *segment, int sample_i) {
	
	if ((*sample_block)->sample_count >= SAMPLE_BLOCK_SIZE) {
	
		(*sample_block)->next_block_open = 1;
	
		if ( (*sample_block)->next_block == NULL ) {
	
			(*sample_block)->next_block = (Sample_block*)malloc(sizeof(Sample_block));
			(*sample_block) = (*sample_block)->next_block; 
			init_sample_block(*sample_block);
		} else {
			(*sample_block) = (*sample_block)->next_block;
			reset_sample_block(*sample_block);
		}
			
	}
	
	init_sample(&((*sample_block)->samples[(*sample_block)->sample_count]), sample_i, segment);
	
	return &((*sample_block)->samples[(*sample_block)->sample_count++]);
	
}


//------------------------------------------------------------------------//
// Iterate through all the samples in the chain of sample blocks
//------------------------------------------------------------------------//


Sample* next_sample(Sample_block **sample_block) { 
	
	Sample *sample;
	
	if ((*sample_block)->next_i >= SAMPLE_BLOCK_SIZE) {
		if ((*sample_block)->next_block_open) {
			(*sample_block) = (*sample_block)->next_block;
			(*sample_block)->next_i = 0;
		} else {
			return NULL;
		}
	} 
	
	if ((*sample_block)->next_i >= (*sample_block)->sample_count ) {
		sample = NULL;
	} else {
		sample = &((*sample_block)->samples[(*sample_block)->next_i++]);
	}		
	
	return sample;
}

//------------------------------------------------------------------------//
// Get the first sample in the block.
//------------------------------------------------------------------------//


Sample* first_sample(Sample_block **sample_block) {
	(*sample_block)->next_i = 0;
	return next_sample(sample_block);
}

//------------------------------------------------------------------------//
// Free all the memory allocated in the sample block chain
//------------------------------------------------------------------------//


void free_sample_blocks(Sample_block *start_block) {
	
	Sample_block *curr_block;
	Sample_block *prev_block;
	
	if (start_block->next_block != NULL) {
		curr_block = start_block->next_block;
		
		while (curr_block->next_block != NULL) {
			prev_block = curr_block;
			curr_block = curr_block->next_block;
			
			free(prev_block);
		}
		
		free(curr_block);
		
	}
}



