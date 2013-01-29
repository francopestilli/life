/*
 *  trim.c
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
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <string.h>

#include "phantom/shared/strand.h"
#include "phantom/shared/strand_collection.h"
#include "phantom/shared/segment.h"
#include "phantom/shared/shared.h"
#include "phantom/trim/strand_section.h"

#include "phantom/trim/trim.h"


#define NEW_STRAND_SECTION 1


void trim(Strand_collection *trim_c, Strand_collection *c, double new_sphere_r, double length_reject_threshold, int save_seperate_bundles) {


	int strand_i, num_trimmed_strands, num_trim_control_points, start_of_strand_control_point_i, control_point_i;
	Strand_section **strand_sections, *strand_section;
	Segment *segment;
	Strand *strand;
	int debug;
	
	double *extra_points;
	int num_extra_points;
	
	strand_sections = (Strand_section**)malloc(sizeof(Strand_section*) * c->num_strands);

	extra_points = (double*)malloc(sizeof(double) * (3 * c->num_strands * 2));

	num_trimmed_strands = 0;
	num_trim_control_points = 0;
	num_extra_points = 0;
	
	for (strand_i = 0; strand_i < c->num_strands; strand_i++) {

	
		strand = &(c->strands[strand_i]);
	
		
		strand_section = strand_sections[strand_i] = NULL;

		segment = strand->pre_segment;
	
		
		if (vector_norm(segment->start_point.pos) < new_sphere_r) {

			strand_section = strand_sections[strand_i] = (Strand_section*)malloc(sizeof(Strand_section));
			strand_section->prev_section = NULL;
			strand_section_init(strand_section, strand, NEW_STRAND_SECTION);
			
			strand_section->pre_point = &(extra_points[num_extra_points * 3]);
			num_extra_points++;
			
			strand_section->pre_point[X] = segment->start_point.pos[X] - segment->disp[X];  
			strand_section->pre_point[Y] = segment->start_point.pos[Y] - segment->disp[Y];
			strand_section->pre_point[Z] = segment->start_point.pos[Z] - segment->disp[Z];						
			
			strand_section->entry_point[X] = segment->start_point.pos[X];
			strand_section->entry_point[Y] = segment->start_point.pos[Y];
			strand_section->entry_point[Z] = segment->start_point.pos[Z];
			strand_section->in_sphere = 1;
			strand_section->start_segment = segment->next_segment;
		} 
		
		
		while (segment->next_segment != NULL) {
		
			segment = segment->next_segment;
			if ( ( (strand_section == NULL) || !(strand_section->in_sphere) ) && ( vector_norm(segment->start_point.pos) < new_sphere_r ) ) {
				
				
				if (strand_section == NULL) {
					strand_section = strand_sections[strand_i] = (Strand_section*)malloc(sizeof(Strand_section));
					strand_section->prev_section = NULL;
					strand_section_init(strand_section, strand, NEW_STRAND_SECTION);
				
				} else {
					strand_section = strand_section_alloc(strand_section, strand);
				}
				
				set_section_entry(strand_section, segment->prev_segment, new_sphere_r);
				
			} else if ( (strand_section != NULL) && strand_section->in_sphere) {
			 
				if (vector_norm(segment->start_point.pos) > new_sphere_r ) {
				
					set_section_exit(strand_section, segment->prev_segment, new_sphere_r);
				
					if ( section_length(strand_section) >= length_reject_threshold)  {
						strand_section->length_accepted = 1;
						num_trimmed_strands++;
					}
					
				} else {
					strand_section->num_control_points++;
					num_trim_control_points++;
				}
					
			} 
			
		
		}
		
		if ((strand_section != NULL) && strand_section->in_sphere) {
		
			strand_section->exit_point[X] = segment->end_point.pos[X];
			strand_section->exit_point[Y] = segment->end_point.pos[Y];
			strand_section->exit_point[Z] = segment->end_point.pos[Z];
			
			strand_section->post_point = &(extra_points[num_extra_points * 3]);
			num_extra_points++;
			
			strand_section->post_point[X] = segment->end_point.pos[X] + segment->disp[X];  
			strand_section->post_point[Y] = segment->end_point.pos[Y] + segment->disp[Y];
			strand_section->post_point[Z] = segment->end_point.pos[Z] + segment->disp[Z];	
				
		
			if ( section_length(strand_section) >= length_reject_threshold)  {
				strand_section->length_accepted = 1;
				num_trimmed_strands++;
			}
		
		}
	
		
		
	}
	

	collection_alloc(trim_c, num_trimmed_strands, num_trim_control_points, c->num_isotropic_regions);
	
	copy_isotropic_regions(trim_c, c->isotropic_regions, c->num_isotropic_regions);
	
	trim_c->num_strands = 0;
	trim_c->num_control_points = 0;
	trim_c->num_bundles = c->num_bundles;
	
	trim_c->sphere_r = new_sphere_r;
	trim_c->fov = new_sphere_r;

	for (strand_i = 0; strand_i < c->num_strands; strand_i++) {
	
	
		if (strand_i == 24) {
			debug = 0;
		}
	
		strand_section = strand_sections[strand_i];
	
		while (strand_section != NULL) {
		
			if (strand_section->length_accepted) {
		
				trim_c->start_points[(trim_c->num_strands * 3) + X] = strand_section->entry_point[X]; 
				trim_c->start_points[(trim_c->num_strands * 3) + Y] = strand_section->entry_point[Y];
				trim_c->start_points[(trim_c->num_strands * 3) + Z] = strand_section->entry_point[Z];						
			
				trim_c->end_points[(trim_c->num_strands * 3) + X] = strand_section->exit_point[X]; 
				trim_c->end_points[(trim_c->num_strands * 3) + Y] = strand_section->exit_point[Y];
				trim_c->end_points[(trim_c->num_strands * 3) + Z] = strand_section->exit_point[Z];		
			
			
				if (strand_section->pre_point != NULL) {
					trim_c->pre_points[(trim_c->num_strands * 3) + X] = strand_section->pre_point[X];
					trim_c->pre_points[(trim_c->num_strands * 3) + Y] = strand_section->pre_point[Y];
					trim_c->pre_points[(trim_c->num_strands * 3) + Z] = strand_section->pre_point[Z];
					
				} else {
				
					trim_c->pre_points[(trim_c->num_strands * 3) + X] = strand_section->entry_point[X] * (1.0 + STRAND_EXT_FRACT);
					trim_c->pre_points[(trim_c->num_strands * 3) + Y] = strand_section->entry_point[Y] * (1.0 + STRAND_EXT_FRACT);
					trim_c->pre_points[(trim_c->num_strands * 3) + Z] = strand_section->entry_point[Z] * (1.0 + STRAND_EXT_FRACT);
					
				}
			
			
				if (strand_section->post_point != NULL) {
					trim_c->post_points[(trim_c->num_strands * 3) + X] = strand_section->post_point[X];
					trim_c->post_points[(trim_c->num_strands * 3) + Y] = strand_section->post_point[Y];
					trim_c->post_points[(trim_c->num_strands * 3) + Z] = strand_section->post_point[Z];
					
				} else {
				
					trim_c->post_points[(trim_c->num_strands * 3) + X] = strand_section->exit_point[X] * (1.0 + STRAND_EXT_FRACT);
					trim_c->post_points[(trim_c->num_strands * 3) + Y] = strand_section->exit_point[Y] * (1.0 + STRAND_EXT_FRACT);
					trim_c->post_points[(trim_c->num_strands * 3) + Z] = strand_section->exit_point[Z] * (1.0 + STRAND_EXT_FRACT);
					
				}
			
			
				trim_c->num_strand_control_points[trim_c->num_strands] = strand_section->num_control_points;
				trim_c->strand_r[trim_c->num_strands] = strand_section->strand->radius;
				
				if (save_seperate_bundles) {
				
					if (strand_section->is_first_section) {
						trim_c->bundle_i_of_strand[trim_c->num_strands] = strand_section->strand->bundle_i;
					} else {
						trim_c->bundle_i_of_strand[trim_c->num_strands] = trim_c->num_bundles;
						trim_c->num_bundles++;
					}
				
				} else {
					trim_c->bundle_i_of_strand[trim_c->num_strands] = strand_section->strand->bundle_i;
				}
			
				segment = strand_section->start_segment;
			
				start_of_strand_control_point_i = trim_c->num_control_points;
			
				for (control_point_i = 0; control_point_i < strand_section->num_control_points; control_point_i++) {
					
					trim_c->control_points[trim_c->num_control_points * 3 + X] = segment->start_point.pos[X]; 
					trim_c->control_points[trim_c->num_control_points * 3 + Y] = segment->start_point.pos[Y];
					trim_c->control_points[trim_c->num_control_points * 3 + Z] = segment->start_point.pos[Z];
					
					trim_c->num_control_points++;
					
					segment = segment->next_segment;
				}
				
				debug = start_of_strand_control_point_i + 3 * trim_c->num_strands;
				
				construct_strand(&(trim_c->strands[trim_c->num_strands]), trim_c->num_strands, trim_c->bundle_i_of_strand[trim_c->num_strands], &(trim_c->control_points[start_of_strand_control_point_i]), &(trim_c->start_points[trim_c->num_strands]), &(trim_c->end_points[trim_c->num_strands]), &(trim_c->pre_points[trim_c->num_strands]), &(trim_c->post_points[trim_c->num_strands]), &(trim_c->segments[start_of_strand_control_point_i + 3 * trim_c->num_strands]), trim_c->num_strand_control_points[trim_c->num_strands], 0.0, trim_c->strand_r[trim_c->num_strands]);	
				
				
				trim_c->num_strands++;
					
			}
		
			strand_section = strand_section->next_section;
		}
	
		
			
		free(strand_sections[strand_i]);
	}
	
	free(strand_sections);
	free(extra_points);
}


