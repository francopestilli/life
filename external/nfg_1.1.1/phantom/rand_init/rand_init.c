/*
 *  rand_init.c
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
 */


#include <stdlib.h>
#include <stdio.h>
#include <dirent.h>
#include <sys/stat.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "phantom/shared/strand_collection.h"
#include "phantom/shared/shared.h"
#include "phantom/rand_init/rand_init.h"

 
/* Randomly initialises a collection of strands by generating pairs of points on a sphere then joining them with evenly spaced control points */
int rand_init_collection(Strand_collection *c, double sphere_r, double control_point_freq, double strand_r_lbound, double strand_r_ubound, double strand_r_buffer_ratio, gsl_rng *rand_gen, int max_attempts, int num_isotropic_regions, Isotropic_region *isotropic_regions) {

	int ubound_num_strands, ubound_num_control_points, num_attempts;
	int reject;
	
	double overlap_dist;

	int point_i, strand_i, start_of_strand_triple_i;
	double segment_disp[3], num_segment_disps;
	double seperation_threshold;

	

	ubound_num_strands = (int)( sphere_r * sphere_r / (strand_r_lbound * strand_r_lbound) + 0.5) * 2; 
	
	ubound_num_control_points = (int)(ubound_num_strands * 2.0 * sphere_r *control_point_freq);
	
	collection_alloc(c, ubound_num_strands, ubound_num_control_points, num_isotropic_regions);
	
	copy_isotropic_regions(c, isotropic_regions, num_isotropic_regions);
	
	/* Randomly generate the 'end-points' (including both start_points and end_points) and radii for the strands to be generated */
	
	c->num_strands = 0;
	c->num_control_points = 0;
	
	num_attempts = 0;
	
	while (num_attempts < max_attempts && c->num_strands < ubound_num_strands) {
		
		c->strand_r[c->num_strands] = gsl_ran_flat(rand_gen, strand_r_lbound, strand_r_ubound);
					
		rand_gen_point_on_sphere(&(c->start_points[c->num_strands * 3]), sphere_r, rand_gen);
		rand_gen_point_on_sphere(&(c->end_points[c->num_strands * 3]), sphere_r, rand_gen);		
		
		
		/* Rejection sampling procedure to favour longer strands and thus better fill the sphere controlior */
		seperation_threshold = gsl_ran_flat(rand_gen, -1.0, 1.0) * sphere_r * sphere_r;
		reject = seperation_threshold < dot_product(&(c->start_points[c->num_strands * 3]), &(c->end_points[c->num_strands * 3]));
		
		if ( !reject ) {
		
			for (strand_i=0; strand_i < c->num_strands; strand_i++) {
				
				/* Ensure that no two 'end points' (both start_points and end_points) are closer than the sum of their respective strands' radii scaled by the 'strand_r_buffer_ratio' */ 
				overlap_dist = strand_r_buffer_ratio * (c->strand_r[c->num_strands] + c->strand_r[strand_i]);
				
				if (	(dist_between_points(&(c->start_points[c->num_strands * 3]), &(c->start_points[strand_i * 3])) < overlap_dist) |
					(dist_between_points(&(c->end_points[c->num_strands * 3]), &(c->end_points[strand_i * 3])) < overlap_dist) |
					(dist_between_points(&(c->start_points[c->num_strands * 3]), &(c->end_points[strand_i * 3])) < overlap_dist) |
					(dist_between_points(&(c->end_points[c->num_strands * 3]), &(c->start_points[strand_i * 3])) < overlap_dist) )
				 {
					reject = 1;
					break;
				}
			}
		}
		
		
		if (reject) {
			num_attempts++;
		} else {
		
			c->num_strand_control_points[c->num_strands] = (int)floor(dist_between_points(&(c->start_points[c->num_strands * 3]), &(c->end_points[c->num_strands * 3])) *control_point_freq); 
			c->num_control_points += c->num_strand_control_points[c->num_strands];
			c->bundle_i_of_strand[strand_i] = strand_i;
			c->num_strands++;
			num_attempts = 0;
		}

	}
	
	
	/* Construct the strand collection including the joining control points from the randomly generated 'end-points' and strand radii. */
	
	point_i = 0;
	
	for (strand_i = 0; strand_i < c->num_strands; strand_i++) {
	
		start_of_strand_triple_i = point_i;
		
		segment_disp[X] =  (c->end_points[strand_i * 3 + X] - c->start_points[strand_i * 3 + X]) / (c->num_strand_control_points[strand_i] + 1);
		segment_disp[Y] =  (c->end_points[strand_i * 3 + Y] - c->start_points[strand_i * 3 + Y]) / (c->num_strand_control_points[strand_i] + 1);
		segment_disp[Z] =  (c->end_points[strand_i * 3 + Z] - c->start_points[strand_i * 3 + Z]) / (c->num_strand_control_points[strand_i] + 1);

	
		for ( ; point_i < (start_of_strand_triple_i + c->num_strand_control_points[strand_i]); point_i ++) {
		
			num_segment_disps = (double)(point_i - start_of_strand_triple_i + 1);
		
			c->control_points[point_i * 3 + X] = c->start_points[strand_i * 3 + X] + num_segment_disps * segment_disp[X];
			c->control_points[point_i * 3 + Y] = c->start_points[strand_i * 3 + Y] + num_segment_disps * segment_disp[Y];
			c->control_points[point_i * 3 + Z] = c->start_points[strand_i * 3 + Z] + num_segment_disps * segment_disp[Z];
		
		}
	
	
		/* Define the reference point for the bending cost of the first segment by extending the start_point/end_point outwards from the sphere.  If you want to use a surface other than a sphere this is how you would define the perpendicular direction.  By defining the pre_point/post_point extended from their start_point/end_point in the perpendicular direction.  */
			
		c->pre_points[strand_i * 3 + X] = c->start_points[strand_i * 3 + X] * (1.0 + STRAND_EXT_FRACT);
		c->pre_points[strand_i * 3 + Y] = c->start_points[strand_i * 3 + Y] * (1.0 + STRAND_EXT_FRACT);
		c->pre_points[strand_i * 3 + Z] = c->start_points[strand_i * 3 + Z] * (1.0 + STRAND_EXT_FRACT);

		c->post_points[strand_i * 3 + X] = c->end_points[strand_i * 3 + X] * (1.0 + STRAND_EXT_FRACT);
		c->post_points[strand_i * 3 + Y] = c->end_points[strand_i * 3 + Y] * (1.0 + STRAND_EXT_FRACT);
		c->post_points[strand_i * 3 + Z] = c->end_points[strand_i * 3 + Z] * (1.0 + STRAND_EXT_FRACT); 
	
		construct_strand(&(c->strands[strand_i]), strand_i, strand_i, &(c->control_points[start_of_strand_triple_i * 3]), &(c->start_points[strand_i * 3]), &(c->end_points[strand_i * 3]), &(c->pre_points[strand_i * 3]), &(c->post_points[strand_i * 3]), &(c->segments[start_of_strand_triple_i + 3 * strand_i]), c->num_strand_control_points[strand_i], 0.0, c->strand_r[strand_i]); 

	}


	return 0;
}

void rand_init_grey_matters(Isotropic_region *isotropic_regions, int num_grey_matters, double grey_matter_r_lbound, double grey_matter_r_ubound, double grey_matter_location_lbound, double grey_matter_location_ubound, double grey_matter_diffusivity, double grey_matter_baseline_signal, double grey_matter_weighting, gsl_rng *rand_gen) {

	int grey_matter_i;
	Isotropic_region *grey_matter;
	double point_on_sphere[3], location_radius, region_radius;

	for (grey_matter_i = 0; grey_matter_i < num_grey_matters; grey_matter_i++) {
		grey_matter = &(isotropic_regions[grey_matter_i]);
		
		location_radius = gsl_ran_flat(rand_gen, grey_matter_location_lbound, grey_matter_location_ubound);
		region_radius = gsl_ran_flat(rand_gen, grey_matter_r_lbound, grey_matter_r_ubound);
		
		rand_gen_point_on_sphere(point_on_sphere, location_radius, rand_gen);
		
		init_isotropic_region(grey_matter, point_on_sphere, region_radius, grey_matter_diffusivity, grey_matter_baseline_signal, grey_matter_weighting);
	 
	}


}

/* As the name suggests, randomly generates a point on a sphere */
void rand_gen_point_on_sphere(double *point, double sphere_r, const gsl_rng *random_gen) {
	
	double length;

	point[X] = gsl_ran_gaussian(random_gen, 1.0);
	point[Y] = gsl_ran_gaussian(random_gen, 1.0);
	point[Z] = gsl_ran_gaussian(random_gen, 1.0);
	
	length = sqrt(point[X] * point[X] + point[Y] * point[Y] + point[Z] * point[Z]);
	
	point[X] = sphere_r * point[X]/length;
	point[Y] = sphere_r * point[Y]/length;
	point[Z] = sphere_r * point[Z]/length;
	
	
}


