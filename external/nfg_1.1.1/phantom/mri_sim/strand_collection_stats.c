/*
 *  strand_collection_stats.c
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
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>

#include "phantom/mri_sim/strand_collection_stats.h"
#include "phantom/mri_sim/segment_stats.h"
#include "phantom/mri_sim/overlap_strands.h"
#include "phantom/shared/shared.h"

#define RADIAL_DISTR_BLOCK_SIZE 1028

Strand_collection_stats *strand_collection_stats_alloc(Strand_collection *collection, int save_ext_segment_stats, int save_ext_fill_radial_stats) {

	int strand_i;

	Strand_collection_stats *stats;
	
	stats = (Strand_collection_stats*)malloc(sizeof(Strand_collection_stats));
	
	
	stats->save_ext_segment_stats = save_ext_segment_stats; 
	stats->save_ext_fill_radial_stats = save_ext_fill_radial_stats; 
	
	stats->collection = collection;
	
	stats->fill_net = 0;
	stats->fill_gross = 0;
	
	stats->subvoxel_in_sphere = 0;
	stats->fill_net_in_sphere = 0;
	stats->fill_gross_in_sphere = 0;
	
	stats->angle_avg = 0.0;
	stats->length_avg = 0.0;
	stats->radius_curv_avg = 0.0;
	/*
	stats->length_max = 0.0;
	stats->length_min = collection->sphere_r * 2.0;	
	stats->angle_max = 0.0;
	stats->radius_curv_min = RADIUS_CURVATURE_UBOUND;*/
	
	stats->angle_stdev = 0.0;
	stats->length_stdev = 0.0;
	stats->radius_curv_stdev = 0.0;
	
	
	if (save_ext_fill_radial_stats) {
		stats->fill_gross_radial_dist = (double*)malloc(sizeof(double) * RADIAL_DISTR_BLOCK_SIZE);
		stats->fill_net_radial_dist = (double*)malloc(sizeof(double) * RADIAL_DISTR_BLOCK_SIZE);
		stats->subvoxel_in_sphere_radial_dist = (double*)malloc(sizeof(double) * RADIAL_DISTR_BLOCK_SIZE);		/* Could be calculated analytically but it is easier to store while going. */
	}
	
	
	stats->segment_stats = (Segment_stats**)malloc(sizeof(Segment_stats*) * collection->num_strands);
	
	for (strand_i = 0; strand_i < collection->num_strands; strand_i++) {
		stats->segment_stats[strand_i] = (Segment_stats*)malloc(sizeof(Segment_stats) * (collection->num_strand_control_points[strand_i]+1));	/* It is the size of the number of points on the strand plus one for the extra segment. */
	}
	
	return stats;


}


void strand_collection_stats_free(Strand_collection_stats *stats) {

	int strand_i;

	for (strand_i = 0; strand_i < stats->collection->num_strands; strand_i++) {
		free(stats->segment_stats[strand_i]);
	}
	
	free(stats->segment_stats);
	
	if (stats->save_ext_fill_radial_stats) {
		free(stats->fill_gross_radial_dist);
		free(stats->fill_net_radial_dist);
		free(stats->subvoxel_in_sphere_radial_dist);
	}

	free(stats);

}

void add_overlap_stats(Strand_collection_stats *stats, Voxel *voxel, double sphere_fov) {



/**********************************************************************
*  Calculate strand overlap and volume fill ratios.	
***********************************************************************/


	int x, y, z, offset;
	Subvoxel *subvoxel;
	Overlap_strand *overlap_strand;

	double radial_disp[3], origin[3], radial_dist, subvoxel_size[3];

	origin[X] = origin[Y] = origin[Z] = 0.0;
	
	
	subvoxel_size[X] = voxel->size[X]/((double)(voxel->num_subvoxels[X])); 
	subvoxel_size[Y] = voxel->size[Y]/((double)(voxel->num_subvoxels[Y]));
	subvoxel_size[Z] = voxel->size[Z]/((double)(voxel->num_subvoxels[Z]));	

	
	for (z = 0; z < voxel->num_subvoxels[Z]; z++) {
		for (y = 0; y < voxel->num_subvoxels[Y]; y++) {
			for (x = 0; x < voxel->num_subvoxels[X]; x++) {			
			
				offset = z * voxel->num_subvoxels[X] * voxel->num_subvoxels[Y] + y * voxel->num_subvoxels[X] + x;

				radial_disp[X] = ((double)x + 0.5) * subvoxel_size[X] + voxel->origin[X];
				radial_disp[Y] = ((double)y + 0.5) * subvoxel_size[Y] + voxel->origin[Y];
				radial_disp[Z] = ((double)z + 0.5) * subvoxel_size[Z] + voxel->origin[Z];
												
				radial_dist = dist_between_points(radial_disp, origin);
				
				
				if (radial_dist <= sphere_fov) {
					stats->subvoxel_in_sphere++;
				}
				
			
				if (stats->save_ext_fill_radial_stats) {
					
				
					if (radial_dist <= sphere_fov) {
						if ( (stats->subvoxel_in_sphere != 0) && (stats->subvoxel_in_sphere % RADIAL_DISTR_BLOCK_SIZE == 0) ) {
							stats->subvoxel_in_sphere_radial_dist = (double*)realloc(stats->subvoxel_in_sphere_radial_dist, sizeof(double) * (stats->subvoxel_in_sphere + RADIAL_DISTR_BLOCK_SIZE));
						}
						stats->subvoxel_in_sphere_radial_dist[stats->subvoxel_in_sphere] = radial_dist;
					}
				}				
				
				subvoxel = &(voxel->subvoxels[offset]);
				
				if (subvoxel->closest_segment != NULL) {
					
					
					if (stats->save_ext_fill_radial_stats) {
						if (stats->save_ext_fill_radial_stats) {
							if ( (stats->fill_net != 0) & (stats->fill_net % RADIAL_DISTR_BLOCK_SIZE == 0) ) {
								stats->fill_net_radial_dist = (double*)realloc(stats->fill_net_radial_dist, sizeof(double) * (stats->fill_net + RADIAL_DISTR_BLOCK_SIZE));
							}
							
							stats->fill_net_radial_dist[stats->fill_net] = radial_dist;
						}
					}
					
					
					stats->fill_net++;
					
					if (radial_dist <= sphere_fov) {
						stats->fill_net_in_sphere++;
					}
					
					stats->segment_stats[subvoxel->closest_segment->strand->strand_i][subvoxel->closest_segment->segment_i].spatial_extent_net++; 
				
				} else if (subvoxel->closest_isotropic_region != NULL) {
					
					stats->fill_net++;
					
					if (radial_dist <= sphere_fov) {
						stats->fill_net_in_sphere++;
					}
					
				}
			
				overlap_strand = subvoxel->overlap_strands;
						
				while (overlap_strand != NULL) {
	
					if (stats->save_ext_fill_radial_stats) {
		
						if ( (stats->fill_gross != 0) & (stats->fill_gross % RADIAL_DISTR_BLOCK_SIZE == 0) ) {
							stats->fill_gross_radial_dist = (double*)realloc(stats->fill_gross_radial_dist, sizeof(double) * (stats->fill_gross + RADIAL_DISTR_BLOCK_SIZE));
						}
						
						stats->fill_gross_radial_dist[stats->fill_gross] = radial_dist;
					}
					
					
					stats->fill_gross++;
	
					if (radial_dist <= sphere_fov) {
						stats->fill_gross_in_sphere++;
					}
	
					stats->segment_stats[overlap_strand->strand->strand_i][overlap_strand->closest_segment->segment_i].spatial_extent_gross++;	
	
					overlap_strand = overlap_strand->next; 

				}
				
				

				
							
			}
		}
	}


	overlap_strand = overlap_strand;
}


void add_length_curv_stats(Strand_collection_stats *stats, Strand_collection *c) {

	int strand_i, segment_i;
	Segment *segment;

	Segment_stats *segment_stats;
	double *start, *middle, *end;
	int num_segments;

	double length_diff, angle_diff, radius_curv_diff;

/**********************************************************************
*  Calculate strand segments length, angle and radius of curvature.	
***********************************************************************/

	num_segments = 0;
	
	for (strand_i = 0; strand_i < stats->collection->num_strands; strand_i++) {
		
		segment = stats->collection->strands[strand_i].start_segment;
		
		start = segment->start_point.pos;
		middle = segment->end_point.pos;
		end = segment->next_segment->end_point.pos;
		
		while (segment->next_segment != NULL) {
			segment_stats = &(stats->segment_stats[strand_i][segment->segment_i]);
			
			init_segment_stats(segment_stats, segment);
			
			segment_stats->segment = segment;
			
			start = segment->start_point.pos;
			middle = segment->end_point.pos;
			end = segment->next_segment->end_point.pos;
			
			segment_stats->angle = angle_between_points(start, middle, end);
			segment_stats->length = dist_between_points(start, middle);
			segment_stats->radius_curv = radius_of_curvature(start, middle, end);
						
			stats->angle_avg += segment_stats->angle;
			stats->length_avg += segment_stats->length;
/*
			if (segment_stats->radius_curv <  stats->radius_curv_min) {
				stats->radius_curv_min = segment_stats->radius_curv;
			}			
			
			if (segment_stats->angle >  stats->angle_max) {
				stats->angle_max = segment_stats->angle;
			}
			
			if (segment_stats->length >  stats->length_max) {
				stats->length_max = segment_stats->length;
			} else if (segment_stats->length <  stats->length_min) {
				stats->length_min = segment_stats->length;
			}
*/			
			/* Cap the stats->radius curvature at RADIUS_CURVATURE_UBOUND */
			if ( (stats->radius_curv_avg >= RADIUS_CURVATURE_UBOUND) | (segment_stats->radius_curv >= RADIUS_CURVATURE_UBOUND) ) {
				stats->radius_curv_avg = RADIUS_CURVATURE_UBOUND;
			} else {
				stats->radius_curv_avg += segment_stats->radius_curv;
			}
			
			
			segment = segment->next_segment;
			num_segments++;
		}
		
	//	segment_stats = &(stats->segment_stats[strand_i][segment->segment_i]);
	//	segment_stats->length = dist_between_points(middle, end);

	}



	stats->length_avg /= (double)num_segments;
	stats->angle_avg /= (double)num_segments;
	
	
	if (stats->radius_curv_avg >= RADIUS_CURVATURE_UBOUND) {
		stats->radius_curv_avg = RADIUS_CURVATURE_UBOUND;
	} else {
		stats->radius_curv_avg /= (double)num_segments;
	}


	for (strand_i = 0; strand_i < stats->collection->num_strands; strand_i++) {
	
		for (segment_i = 1; segment_i < stats->collection->num_strand_control_points[strand_i]; segment_i++) {
		
			length_diff = stats->segment_stats[strand_i][segment_i].length - stats->length_avg;
			angle_diff = stats->segment_stats[strand_i][segment_i].angle - stats->angle_avg;
			radius_curv_diff = stats->segment_stats[strand_i][segment_i].radius_curv - stats->radius_curv_avg;
			
			stats->length_stdev += length_diff * length_diff;
			stats->angle_stdev += angle_diff * angle_diff;
			stats->radius_curv_stdev += radius_curv_diff * radius_curv_diff;

		}
		
		length_diff = stats->segment_stats[strand_i][segment_i].length - stats->length_avg;
		stats->length_stdev += length_diff * length_diff;
		
	}

	
	stats->length_stdev = sqrt(stats->length_stdev/(double)num_segments);
	stats->angle_stdev = sqrt(stats->angle_stdev/(double)num_segments);
	stats->radius_curv_stdev = sqrt(stats->radius_curv_stdev/(double)num_segments);	

}


int save_strand_collection_stats(Strand_collection_stats *stats, char *output_dir_path) {

	FILE *file;
	DIR *dir;
	char *path;
	
	int strand_i, segment_i; 
		
	double overlap_fraction, fill_fraction;
	
	int output_i, segment_stats_size;
	int *fill_gross_output,	*fill_net_output, *new_strand_markers_output;
	double *overlap_fraction_output, *angle_output, *length_output, *radius_curv_output;
	
	int segment_fill_gross, segment_fill_net;
	double segment_overlap_fraction;
	
	int error;
	
	
	double overlap_fraction_min;
	double angle_min;
	double length_min;
	double radius_curv_min;

	double overlap_fraction_max;
	double angle_max;
	double length_max;
	double radius_curv_max;
	
	double overlap_fraction_median;
	double overlap_fraction_95conf_low;
	double overlap_fraction_95conf_high;
	double overlap_fraction_1quartile;
	double overlap_fraction_3quartile;
	
	double angle_median;
	double angle_95conf_low;
	double angle_95conf_high;
	double angle_1quartile;
	double angle_3quartile;
	
	double radius_curv_median;
	double radius_curv_95conf_low;
	double radius_curv_95conf_high;
	double radius_curv_1quartile;
	double radius_curv_3quartile;
	
	double length_median;
	double length_95conf_low;
	double length_95conf_high;
	double length_1quartile;
	double length_3quartile;

	int landmark_index;
	double landmark_index_unrounded, remaining_fraction;
	
	error = 0;
	
	
	printf("\nSaving statistics ...\n");
	
	overlap_fraction =  ((double)(stats->fill_gross_in_sphere - stats->fill_net_in_sphere)) / (double)(stats->fill_gross_in_sphere);
	fill_fraction = ((double)(stats->fill_net_in_sphere)) / ((double) stats->subvoxel_in_sphere);
	
	
	segment_stats_size = 0;
	stats->overlap_fraction_max = 0.0;
	
	for (strand_i = 0; strand_i < stats->collection->num_strands; strand_i++) {
	
		/*First and last segments are excluded as they are likely to have just been trimmed */
		for (segment_i = 1; segment_i < (stats->collection->num_strand_control_points[strand_i]-1); segment_i++) {
		
			segment_fill_gross = stats->segment_stats[strand_i][segment_i].spatial_extent_gross;
			segment_fill_net = stats->segment_stats[strand_i][segment_i].spatial_extent_net;
		
			if (segment_fill_gross == 0) {
				segment_overlap_fraction = 0;
			} else {
				segment_overlap_fraction = ((double)(segment_fill_gross - segment_fill_net )) / (double)segment_fill_gross; 
			}
			
			if (segment_overlap_fraction > stats->overlap_fraction_max) {
				stats->overlap_fraction_max = segment_overlap_fraction;
			}
			
			stats->overlap_fraction_stdev += (segment_overlap_fraction - overlap_fraction) * (segment_overlap_fraction - overlap_fraction);
			segment_stats_size++;
		
		}
		
	}

	
	stats->overlap_fraction_stdev = sqrt(stats->overlap_fraction_stdev/(double)segment_stats_size);




	if ((dir = opendir(output_dir_path)) == NULL) {
		if ((error = mkdir(output_dir_path, MY_PERMS)) > 0) {
			printf("Error! Could not create directory %s, (Error code: %d)!\n", output_dir_path, error);
			return 1;
		}
	} else {
		closedir(dir);
	}	


	path = (char*)malloc(sizeof(char) * (strlen(output_dir_path) + 100));



	
	
	if (stats->save_ext_fill_radial_stats) {
		
		printf("\nSaving extended filling statistics as function of radius ...\n");

		
		/*Gross fill radial distances.*/

		sprintf(path, "%s%cfill_gross_radial_dist.double", output_dir_path, DIR_SEP);

		if ( (file =fopen(path,"wb")) == NULL) {
			printf("\nCould not open file %s.!\n", path);
			return 1;
		} 

		fwrite(stats->fill_gross_radial_dist, sizeof(double), stats->fill_gross, file);
		
		fclose(file);


		//Net fill radial distances.

		sprintf(path, "%s%cfill_net_radial_dist.double", output_dir_path, DIR_SEP);

		if ( (file =fopen(path,"wb")) == NULL) {
			printf("\nCould not open file %s.!\n", path);
			return 1;
		} 

		fwrite(stats->fill_net_radial_dist, sizeof(double), stats->fill_net, file);
		
		fclose(file);


		/*Subvoxels in sphere radial distances.*/

		sprintf(path, "%s%csubvoxel_in_sphere_radial_dist.double", output_dir_path, DIR_SEP);

		if ( (file =fopen(path,"wb")) == NULL) {
			printf("\nCould not open file %s.!\n", path);
			return 1;
		} 

		fwrite(stats->subvoxel_in_sphere_radial_dist, sizeof(double), stats->subvoxel_in_sphere, file);

		fclose(file);

	}
	
	
	fill_gross_output = (int*)malloc(sizeof(int) * segment_stats_size);
	fill_net_output = (int*)malloc(sizeof(int) * segment_stats_size);
	length_output = (double*)malloc(sizeof(double) * segment_stats_size);
	angle_output = (double*)malloc(sizeof(double) * segment_stats_size);
	radius_curv_output = (double*)malloc(sizeof(double) * segment_stats_size);
	new_strand_markers_output = (int*)malloc(sizeof(int) * segment_stats_size);
	overlap_fraction_output = (double*)malloc(sizeof(double) * segment_stats_size);

	/*Gross fill per segment.*/
	
	output_i = 0;

	for (strand_i = 0; strand_i < stats->collection->num_strands; strand_i++) {
	
		/*First and last segments are excluded as they are likely to have just been trimmed */
		for (segment_i = 1; segment_i <  (stats->collection->num_strand_control_points[strand_i]-1); segment_i++) {
		
			fill_gross_output[output_i] = stats->segment_stats[strand_i][segment_i].spatial_extent_gross;
			fill_net_output[output_i] = stats->segment_stats[strand_i][segment_i].spatial_extent_net;
			
			if (fill_gross_output[output_i] == 0) {
				overlap_fraction_output[output_i] = 0.0;
			} else {
				overlap_fraction_output[output_i] = ((double)(fill_gross_output[output_i] - fill_net_output[output_i])) / (double)(fill_gross_output[output_i]);
			}
			
			angle_output[output_i] = stats->segment_stats[strand_i][segment_i].angle;
			length_output[output_i] = stats->segment_stats[strand_i][segment_i].length;
			radius_curv_output[output_i] = stats->segment_stats[strand_i][segment_i].radius_curv;		

			if (segment_i == 1) {
				new_strand_markers_output[output_i] = 1;
			} else {
				new_strand_markers_output[output_i] = 0;
			}
			
			output_i++;
		}
		
	}
	
		
	if (stats->save_ext_segment_stats) {	
		
		printf("\nSaving extended segment statistics ...\n");
		
		sprintf(path, "%s%csegment_fill_gross.int", output_dir_path, DIR_SEP);
		
		if ( (file =fopen(path,"wb")) == NULL) {
			printf("\nCould not open file %s.!\n", path);
			return 1;
		} 
		
		fwrite(fill_gross_output, sizeof(int), segment_stats_size, file);
		
		fclose(file);
		
		
		
		sprintf(path, "%s%csegment_fill_net.int", output_dir_path, DIR_SEP);
		
		if ( (file =fopen(path,"wb")) == NULL) {
			printf("\nCould not open file %s.!\n", path);
			return 1;
		} 
		
		fwrite(fill_net_output, sizeof(int), segment_stats_size, file);
		
		fclose(file);
		
		sprintf(path, "%s%csegment_overlap_fraction.double", output_dir_path, DIR_SEP);
		
		if ( (file =fopen(path,"wb")) == NULL) {
			printf("\nCould not open file %s.!\n", path);
			return 1;
		} 
		
		fwrite(overlap_fraction_output, sizeof(double), segment_stats_size, file);
		
		fclose(file);	
		
		sprintf(path, "%s%csegment_length.double", output_dir_path, DIR_SEP);
		
		if ( (file =fopen(path,"wb")) == NULL) {
			printf("\nCould not open file %s.!\n", path);
			return 1;
		} 
		
		fwrite(length_output, sizeof(double), segment_stats_size, file);
		
		fclose(file);
		
		
		sprintf(path, "%s%csegment_angle.double", output_dir_path, DIR_SEP);
		
		if ( (file =fopen(path,"wb")) == NULL) {
			printf("\nCould not open file %s.!\n", path);
			return 1;
		} 
		
		fwrite(angle_output, sizeof(double), segment_stats_size, file);
		
		fclose(file);

		
		sprintf(path, "%s%csegment_radius_curv.double", output_dir_path, DIR_SEP);
		
		if ( (file =fopen(path,"wb")) == NULL) {
			printf("\nCould not open file %s.!\n", path);
			return 1;
		} 
		
		fwrite(radius_curv_output, sizeof(double), segment_stats_size, file);
		
		fclose(file);
		
		sprintf(path, "%s%cnew_strand_markers.int", output_dir_path, DIR_SEP);
	
		if ( (file =fopen(path,"wb")) == NULL) {
			printf("\nCould not open file %s.!\n", path);
			return 1;
		} 
		
		fwrite(new_strand_markers_output, sizeof(int), segment_stats_size, file);
		
		fclose(file);
		
		

			
	}

	
	qsort(overlap_fraction_output, segment_stats_size, sizeof(double), compare_dble);
	qsort(angle_output, segment_stats_size, sizeof(double), compare_dble);
	qsort(length_output, segment_stats_size, sizeof(double), compare_dble);
	qsort(radius_curv_output, segment_stats_size, sizeof(double), compare_dble);
	
	
	overlap_fraction_min =  overlap_fraction_output[0];
	angle_min = angle_output[0];
	length_min = length_output[0];
	radius_curv_min = radius_curv_output[0];
	
	overlap_fraction_max =  overlap_fraction_output[segment_stats_size-1];
	angle_max = angle_output[segment_stats_size-1];
	length_max = length_output[segment_stats_size-1];
	radius_curv_max = radius_curv_output[segment_stats_size-1];
	
	
	landmark_index_unrounded = ((double)segment_stats_size) * 0.5;
	landmark_index = (int)floor(landmark_index_unrounded);
	remaining_fraction = landmark_index_unrounded - (double)landmark_index;
	
	overlap_fraction_median = ( 1 - remaining_fraction)  * overlap_fraction_output[landmark_index] + remaining_fraction * overlap_fraction_output[landmark_index + 1];
	angle_median = ( 1 - remaining_fraction)  * angle_output[landmark_index] + remaining_fraction * angle_output[landmark_index + 1];
	length_median = ( 1 - remaining_fraction)  * length_output[landmark_index] + remaining_fraction * length_output[landmark_index + 1];
	radius_curv_median = ( 1 - remaining_fraction)  * radius_curv_output[landmark_index] + remaining_fraction * radius_curv_output[landmark_index + 1];
	
	
	landmark_index_unrounded = ((double)segment_stats_size) * 0.025;
	landmark_index = (int)floor(landmark_index_unrounded);
	remaining_fraction = landmark_index_unrounded - (double)landmark_index;

	overlap_fraction_95conf_low = ( 1 - remaining_fraction)  * overlap_fraction_output[landmark_index] + remaining_fraction * overlap_fraction_output[landmark_index + 1];
	angle_95conf_low = ( 1 - remaining_fraction)  * angle_output[landmark_index]+ remaining_fraction * angle_output[landmark_index + 1];
	length_95conf_low = ( 1 - remaining_fraction)  * length_output[landmark_index]+ remaining_fraction * length_output[landmark_index + 1];
	radius_curv_95conf_low = ( 1 - remaining_fraction)  * radius_curv_output[landmark_index] + remaining_fraction * radius_curv_output[landmark_index + 1];

	landmark_index_unrounded = ((double)segment_stats_size) * 0.975;
	landmark_index = (int)floor(landmark_index_unrounded);
	remaining_fraction = landmark_index_unrounded - (double)landmark_index;

	overlap_fraction_95conf_high = ( 1 - remaining_fraction)  * overlap_fraction_output[landmark_index] + remaining_fraction * overlap_fraction_output[landmark_index + 1];
	angle_95conf_high = ( 1 - remaining_fraction)  * angle_output[landmark_index]+ remaining_fraction * angle_output[landmark_index + 1];
	length_95conf_high = ( 1 - remaining_fraction)  * length_output[landmark_index]+ remaining_fraction * length_output[landmark_index + 1];
	radius_curv_95conf_high = ( 1 - remaining_fraction)  * radius_curv_output[landmark_index] + remaining_fraction * radius_curv_output[landmark_index + 1];
	
	
	landmark_index_unrounded = ((double)segment_stats_size) * 0.25;
	landmark_index = (int)floor(landmark_index_unrounded);
	remaining_fraction = landmark_index_unrounded - (double)landmark_index;

	overlap_fraction_1quartile = ( 1 - remaining_fraction)  * overlap_fraction_output[landmark_index] + remaining_fraction * overlap_fraction_output[landmark_index + 1];
	angle_1quartile = ( 1 - remaining_fraction)  * angle_output[landmark_index]+ remaining_fraction * angle_output[landmark_index + 1];
	length_1quartile = ( 1 - remaining_fraction)  * length_output[landmark_index]+ remaining_fraction * length_output[landmark_index + 1];
	radius_curv_1quartile = ( 1 - remaining_fraction)  * radius_curv_output[landmark_index] + remaining_fraction * radius_curv_output[landmark_index + 1];	
	

	landmark_index_unrounded = ((double)segment_stats_size) * 0.75;
	landmark_index = (int)floor(landmark_index_unrounded);
	remaining_fraction = landmark_index_unrounded - (double)landmark_index;

	overlap_fraction_3quartile = ( 1 - remaining_fraction)  * overlap_fraction_output[landmark_index] + remaining_fraction * overlap_fraction_output[landmark_index + 1];
	angle_3quartile = ( 1 - remaining_fraction)  * angle_output[landmark_index] + remaining_fraction * angle_output[landmark_index + 1];
	length_3quartile = ( 1 - remaining_fraction)  * length_output[landmark_index] + remaining_fraction * length_output[landmark_index + 1];
	radius_curv_3quartile = ( 1 - remaining_fraction)  * radius_curv_output[landmark_index] + remaining_fraction * radius_curv_output[landmark_index + 1];	


	free(fill_gross_output);
	free(fill_net_output);		
	free(angle_output);
	free(length_output);
	free(radius_curv_output);
	free(new_strand_markers_output);
	
	
	
	sprintf(path, "%s%cstats.txt", output_dir_path, DIR_SEP);
	
	/*Overall statistics*/



	if ( (file =fopen(path,"w")) == NULL) {
		printf("\nCould not open file %s.!\n", path);
		return 1;
	} 

	fprintf(file, "Fill fraction of sphere: %g\n", fill_fraction);
	fprintf(file, "Total overlap fraction (in sphere): %g\n", overlap_fraction);


	
	fprintf(file, "\n");	
	
	fprintf(file, "Segment overlap fraction min: %g\n", overlap_fraction_min);
	fprintf(file, "Segment overlap fraction max: %g\n", overlap_fraction_max);
	fprintf(file, "Segment overlap fraction median: %g\n", overlap_fraction_median);
	fprintf(file, "Segment overlap fraction 2.5th percentile: %g\n", overlap_fraction_95conf_low);
	fprintf(file, "Segment overlap fraction 97.5th percentile: %g\n", overlap_fraction_95conf_high);	
	fprintf(file, "Segment overlap fraction 1st quartile: %g\n", overlap_fraction_1quartile);
	fprintf(file, "Segment overlap fraction 3rd quartile: %g\n", overlap_fraction_3quartile);
	fprintf(file, "Segment overlap fraction standard deviation: %g\n", stats->overlap_fraction_stdev);

	fprintf(file, "\n");
	
	fprintf(file, "Angle min: %g\n", angle_min);
	fprintf(file, "Angle max: %g\n", angle_max);
	fprintf(file, "Angle median: %g\n", angle_median);
	fprintf(file, "Angle 2.5th percentile: %g\n", angle_95conf_low);
	fprintf(file, "Angle 97.5th percentile: %g\n", angle_95conf_high);
	fprintf(file, "Angle 1st quartile: %g\n", angle_1quartile);
	fprintf(file, "Angle 3rd quartile: %g\n", angle_3quartile);	
	fprintf(file, "Angle average: %g\n", stats->angle_avg);
	fprintf(file, "Angle standard deviation: %g\n", stats->angle_stdev);
	
	fprintf(file, "\n");	

	fprintf(file, "Length max: %g\n", length_max); 
	fprintf(file, "Length min: %g\n", length_min);
	fprintf(file, "Length median: %g\n", length_median);
	fprintf(file, "Length 2.5th percentile: %g\n", length_95conf_low);
	fprintf(file, "Length 97.5th percentile: %g\n", length_95conf_high);	
	fprintf(file, "Length 1st quartile: %g\n", length_1quartile);
	fprintf(file, "Length 3rd quartile: %g\n", length_3quartile);
	fprintf(file, "Length average: %g\n", stats->length_avg);
	fprintf(file, "Length standard deviation: %g\n", stats->length_stdev);

	fprintf(file, "\n");	
	
	
	fprintf(file, "Radius of curvature min: %g\n", radius_curv_min);
	fprintf(file, "Radius of curvature max: %g\n", radius_curv_max);
	fprintf(file, "Radius of curvature median: %g\n", radius_curv_median);
	fprintf(file, "Radius of curvature 2.5th percentile: %g\n", radius_curv_95conf_low);
	fprintf(file, "Radius of curvature 97.5th percentile: %g\n", radius_curv_95conf_high);	
	fprintf(file, "Radius of curvature 1st quartile: %g\n", radius_curv_1quartile);
	fprintf(file, "Radius of curvature 3rd quartile: %g\n", radius_curv_3quartile);
	fprintf(file, "Radius of curvature average: %g\n", stats->radius_curv_avg);	
	fprintf(file, "Radius of curvature standard deviation: %g\n", stats->radius_curv_stdev);


	fprintf(file, "\n");

	fprintf(file, "Total grid subvoxels in sphere: %d\n", stats->subvoxel_in_sphere);	
	fprintf(file, "Net grid subvoxels filled (in sphere): %d\n", stats->fill_net_in_sphere);					
	fprintf(file, "Gross grid subvoxels filled (in sphere): %d\n", stats->fill_gross_in_sphere);

	fprintf(file, "\n");

	
	fclose(file);
	
	
	free(path);	



	return error;

}

int compare_dble(const void *a_ptr, const void *b_ptr) {
	
	int comparison;
	double a, b;
	
	a = *(double*)a_ptr;
	b = *(double*)b_ptr;
		
	if (a == b)
		comparison = 0;
	else if (a > b)
		comparison = 1;
	else 
		comparison = -1;
	
	
	return comparison;
	
}

