/*
 *  segment_register.c
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

#include "phantom/mri_sim/segment_register.h"
#include "phantom/shared/segment.h"


Segment_register* segment_register_alloc(Segment *segment) {

	Segment_register *segment_register;
	
	segment_register = (Segment_register*)malloc(sizeof(Segment_register));
	
	segment_register->segment = segment;
	segment_register->next = NULL;
	
	return segment_register;
	
	
}

void segment_register_free(Segment_register *segment_register) {

	free(segment_register);

}
