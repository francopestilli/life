/*
 *  strand_collection_stats.h
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

#ifndef COLLECTION_STATS_H
#define COLLECTION_STATS_H

#include "phantom/shared/segment.h"
#include "phantom/shared/strand_collection.h"

#include "phantom/mri_sim/segment_stats.h"
#include "phantom/mri_sim/voxel.h"


#define RADIUS_CURVATURE_UBOUND 1e30


/*#define SAVE_RADIAL_FILL*/

typedef struct _strand_collection_stats {

	Strand_collection *collection;
	
	int save_ext_segment_stats;
	int save_ext_fill_radial_stats;
		
	double length_avg;
	double angle_avg;
	double radius_curv_avg;
	
	double length_stdev;
	double angle_stdev;
	double radius_curv_stdev;
	
	/* double length_max;
	double length_min;
	double angle_max;
	double radius_curv_min; */
		
	int fill_gross;
	int fill_net;
	
	double overlap_fraction_stdev;
	double overlap_fraction_max;
				
	int fill_gross_in_sphere;
	int fill_net_in_sphere;	
	int subvoxel_in_sphere;
		
	double *fill_gross_radial_dist;
	double *fill_net_radial_dist;
	double *subvoxel_in_sphere_radial_dist;	
		
	Segment_stats **segment_stats;		/* Statistics on each segment organised into strands.*/

} Strand_collection_stats;


Strand_collection_stats *strand_collection_stats_alloc(Strand_collection *collection, int save_ext_segment_stats, int save_ext_fill_radial_stats);

void strand_collection_stats_free(Strand_collection_stats *stats);

void add_overlap_stats(Strand_collection_stats *stats, Voxel *voxel, double sphere_fov);

void add_length_curv_stats(Strand_collection_stats *stats, Strand_collection *c);

int save_strand_collection_stats(Strand_collection_stats *stats, char *output_dir_path);

int compare_dble(const void *a, const void *b);

#endif
