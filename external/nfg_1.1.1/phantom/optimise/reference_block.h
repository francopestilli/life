/*
 *  reference_block.h
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

#ifndef REFERENCE_BLOCK_H
#define REFERENCE_BLOCK_H 0

#include "phantom/optimise/sample.h"


#define REFERENCE_BLOCK_SIZE 1000

typedef struct _reference_block {
	
	//The index of the next free reference in the 'reference reference' array.
	int reference_count;
	
	int next_i;
	
	//Create an array to hold the reference references
	Sample *references[REFERENCE_BLOCK_SIZE];
	
	int init_check;
	
	int next_block_open;
			
	//Reference to the next reference reference block in the linked list (will be NULL if the reference reference block is the last in the list).
	struct _reference_block *next_block;
	

	
} Reference_block;



void reset_reference_block(Reference_block *reference_block);

void init_reference_block(Reference_block *reference_block);

void add_reference(Reference_block **reference_block, Sample *sample);

Sample* next_reference(Reference_block **reference_block);

Sample* first_reference(Reference_block **reference_block); 


void free_reference_blocks(Reference_block *start_block);

#endif
