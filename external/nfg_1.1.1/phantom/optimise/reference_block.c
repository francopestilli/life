/*
 *  reference_block.c
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
#include <stdio.h>

#include "phantom/optimise/reference_block.h"
#include "phantom/optimise/sample.h"


//------------------------------------------------------------------------//
// Initialize a newly allocated reference block.
//------------------------------------------------------------------------//

void reset_reference_block(Reference_block *reference_block) {  //For reuse of the memory allocated during previous iteration resets all the values to their initial values except the 'next_block' pointers to the next block in the chain.
	if (reference_block->init_check != -7) {
	
		printf("Reference block was never initiated\n");
		exit(0);
	}


	
	reference_block->next_i = 0;
	reference_block->reference_count = 0;
	reference_block->next_block_open = 0;

	
}

void init_reference_block(Reference_block *reference_block) {

	reference_block->next_i = 0;
	reference_block->reference_count = 0;
	reference_block->next_block_open = 0;
	reference_block->next_block = NULL;
	reference_block->init_check = -7;

}





//------------------------------------------------------------------------//
// Append a new reference to the end of the sample block
//------------------------------------------------------------------------//

void add_reference(Reference_block **reference_block, Sample *sample) {
	
	if ((*reference_block)->reference_count >= REFERENCE_BLOCK_SIZE) {
	
		(*reference_block)->next_block_open = 1;
		
		
		if ((*reference_block)->next_block == NULL) {
			(*reference_block)->next_block = (Reference_block*)malloc(sizeof(Reference_block));
			*reference_block = (*reference_block)->next_block; 
			init_reference_block(*reference_block);
		} else {
			*reference_block = (*reference_block)->next_block;
			reset_reference_block(*reference_block);
		}
	}

	(*reference_block)->references[(*reference_block)->reference_count++] = sample;
	
}

//------------------------------------------------------------------------//
// Iterate through all the references in the chain of sample blocks
//------------------------------------------------------------------------//

Sample* next_reference(Reference_block **reference_block) { 
	
	Sample *reference;
	
	if ((*reference_block)->next_i >= REFERENCE_BLOCK_SIZE) {
		
		
		if ((*reference_block)->next_block_open) {
			*reference_block = (*reference_block)->next_block;
			(*reference_block)->next_i=0;
		
		} else {
			return NULL;
		}
		
	} 
	
	
	if ((*reference_block)->next_i >= (*reference_block)->reference_count ) {
		reference = NULL;
	} else {
		reference = (*reference_block)->references[(*reference_block)->next_i++];
	}		
	
	return reference;
}


//------------------------------------------------------------------------//
// Get the first sample in the reference block chain.
//------------------------------------------------------------------------//

Sample* first_reference(Reference_block **reference_block) {
	(*reference_block)->next_i=0;
	return next_reference(reference_block);
}

//------------------------------------------------------------------------//
// Free all the memory allocated in the reference block chain
//------------------------------------------------------------------------//

void free_reference_blocks(Reference_block *start_block) {
	
	Reference_block *curr_block;
	Reference_block *prev_block;
	
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

