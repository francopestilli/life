/*
 *  strand_section.c
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

#include "phantom/trim/strand_section.h"
#include "phantom/shared/shared.h"
#include "phantom/shared/segment.h"



Strand_section* strand_section_alloc(Strand_section *prev_section, Strand *strand) {

	Strand_section *new_section;
	
	new_section = (Strand_section*)malloc(sizeof(Strand_section));
	prev_section->next_section = new_section;
	new_section->prev_section = prev_section;
	
	strand_section_init(new_section, strand, 0);
	
	return new_section;
} 


void strand_section_init(Strand_section *section, Strand *strand, int is_first_section) {
	section->next_section = NULL;
	section->num_control_points = 0;
	section->in_sphere = 0;
	section->length_accepted = 0;
	section->strand = strand;
	section->is_first_section = is_first_section;


}

void strand_section_free(Strand_section *section) {
	
	while (section->next_section != NULL) {
	
		section = section->next_section;
		free(section->prev_section);
	}
	
	free(section);

}

void set_section_entry(Strand_section *section, Segment *segment, double sphere_r) {

	sphere_intersect(section->entry_point, segment->start_point.pos, segment->disp, sphere_r);
	section->pre_point = segment->start_point.pos;
	section->in_sphere = 1;
	section->start_segment = segment->next_segment;


}


void set_section_exit(Strand_section *section, Segment *segment, double sphere_r) {

	sphere_intersect(section->exit_point, segment->start_point.pos, segment->disp, sphere_r);
	section->post_point = segment->end_point.pos;
	section->in_sphere = 0;

}

double section_length(Strand_section *section) {
	
	return dist_between_points(section->entry_point, section->exit_point);

}
