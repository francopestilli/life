/*
 *  draw_roi.c
 *  Numerical Fibre Generator
 *
 *  Created by Tom Close on 11/12/08.
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


#include "phantom/draw_rois/draw_rois.h"



int main(int argc, char *argv[]) {

	int error;

	Strand_collection c;
	int num_voxels, point_depth_l_bound, point_depth_u_bound;
	double voxel_size, extra_point_radius;
	DIR *dir;
	FILE *param_file, *file;
	char *param_path, *input_dir_path,  *output_dir_path, *path;

//	unsigned char *seperate_masks;
	int *combined_masks, num_elems,  bundle_i, num_voxels_pdim[3];
//	int elem_i, max_width, seperate_masks,bundle_i;

	double voxel_size_pdim[3], voxel_radial_l_bound, voxel_radial_u_bound, inclusion_radius, null_mask_l_bound, null_mask_u_bound;

	int save_combined_mask, subbundle_mask;
	int max_bundle_i, *cumulative_num_strands, *bundle_is_excluded;

	int num_params, num_lines, stage_num;
	char line[200], key[200];
	double value;

	Bundle *bundle;
	double cross_sectional_area, average_strand_length;
	int start_roi_id, end_roi_id;
	int i;


	if (argc != 4) {

		printf("\nNumerical Fibre Generator (%s): 'trim'\n\nNumber of suplied arguments %d.\n\n Usage: \n\t arg[1]: input directory \n\t arg[2]: output directory \n\t arg[3]: parameters file\n\n", VERSION_NUM, (argc-1));
		exit(EXIT_FAILURE);
	}

	param_path = argv[3];

	input_dir_path = argv[1];
	output_dir_path = argv[2];


	if ((param_file = fopen(param_path, "r")) == NULL) {
		printf("Error! Could not open file %s\n", param_path);
		exit(EXIT_FAILURE);
	}

	if ((dir = opendir(output_dir_path)) == NULL) {
		if ((error = mkdir(output_dir_path, MY_PERMS)) > 0) {
			printf("Error! Could not create directory %s, (Error code: %d)!\n", argv[0], error);
			exit(EXIT_FAILURE);
		}
	} else {
		closedir(dir);
	}


	num_voxels = NUM_VOXELS_DEFAULT;
	voxel_size = VOXEL_SIZE_DEFAULT;
	point_depth_l_bound = POINT_DEPTH_L_BOUND_DEFAULT;
	point_depth_u_bound = POINT_DEPTH_U_BOUND_DEFAULT;
	extra_point_radius = EXTRA_POINT_RADIUS_DEFAULT;
	voxel_radial_l_bound = VOXEL_RADIAL_L_BOUND_DEFAULT;
	voxel_radial_u_bound = VOXEL_RADIAL_U_BOUND_DEFAULT;
//	save_combined_mask = SAVE_COMBINED_MASK_DEFAULT;
	subbundle_mask = SUBBUNDLE_MASK_DEFAULT;
	inclusion_radius = INCLUSION_RADIUS_DEFAULT;
	null_mask_l_bound = NULL_MASK_L_BOUND_DEFAULT;
	null_mask_u_bound = NULL_MASK_U_BOUND_DEFAULT;


	num_params = 0;
	num_lines = 0;

	printf("\n\nReading parameters file: %s\n\n", param_path);
	fflush(stdout);


	while (fgets(line, 200, param_file) != NULL) {

		if (sscanf(line, "%s %lf", key, &value) != EOF) {

			if (!strcmp(key, "num_voxels")) {
				num_voxels = (int)value;
				printf("Read value for num_voxels, %d\n", num_voxels);

			} else if (!strcmp(key, "voxel_size")) {
				voxel_size = value;
				printf("Read value for voxel_size, %g\n", voxel_size);

			} else if (!strcmp(key, "extra_point_radius")) {
				extra_point_radius = value;
				printf("Read value for extra_point_radius, %g\n", extra_point_radius);

			} else if (!strcmp(key, "point_depth_l_bound")) {
				point_depth_l_bound = (int)value;
				printf("Read value for point_depth_l_bound, %d\n", point_depth_l_bound);

			} else if (!strcmp(key, "point_depth_u_bound")) {
				point_depth_u_bound = (int)value;
				printf("Read value for point_depth_u_bound, %d\n", point_depth_u_bound);

			} else if (!strcmp(key, "voxel_radial_l_bound")) {
				voxel_radial_l_bound = value;
				printf("Read value for voxel_radial_l_bound, %g\n", voxel_radial_l_bound);

			} else if (!strcmp(key, "voxel_radial_u_bound")) {
				voxel_radial_u_bound = value;
				printf("Read value for voxel_radial_u_bound, %g\n", voxel_radial_u_bound);

			} else if (!strcmp(key, "inclusion_radius")) {
				inclusion_radius = value;
				printf("Read value for inclusion_radius, %g\n", inclusion_radius);

			} else if (!strcmp(key, "null_mask_l_bound")) {
				null_mask_l_bound = value;
				printf("Read value for null_mask_l_bound, %g\n", null_mask_l_bound);

			} else if (!strcmp(key, "null_mask_u_bound")) {
				null_mask_u_bound = value;
				printf("Read value for null_mask_u_bound, %g\n", null_mask_u_bound);

			} else if (!strcmp(key, "subbundle_mask")) {
				subbundle_mask = (int)value;
				printf("Read value for subbundle_mask, %d\n", subbundle_mask);

//			} else if (!strcmp(key, "save_combined_mask")) {
//				save_combined_mask = (int)value;
//				printf("Read value for save_combined_mask, %d\n", save_combined_mask);

			} else {
				num_params--;
			}
			num_params++;
		}

		num_lines++;
	}



	printf("\nRead in %d (of 12) key/value pairs\n", num_params);
	fflush(stdout);


	if (load_collection(&c, input_dir_path)) {
		printf("Could not load strands from directory %s, exiting ...\n", input_dir_path);
		exit(EXIT_FAILURE);

	}


	num_voxels_pdim[X] = num_voxels_pdim[Y] = num_voxels_pdim[Z] = num_voxels;
	voxel_size_pdim[X] = voxel_size_pdim[Y] = voxel_size_pdim[Z] = voxel_size;


	if (c.num_bundles == 0) {
		printf("Error! the number of bundles in the collection is 0!");
		exit(EXIT_FAILURE);
	}


	max_bundle_i = c.bundles[c.num_bundles-1].bundle_i;
	cumulative_num_strands = (int*)malloc(sizeof(int) * (max_bundle_i + 2));

	//Get the cumulative number of strands after each bundle is added.  This is used to index the ROIs.
	cumulative_num_strands[0] = 0;

	for (bundle_i = 0; bundle_i < max_bundle_i; bundle_i++) {

		cumulative_num_strands[bundle_i + 1] = cumulative_num_strands[bundle_i] + c.bundles[bundle_i].num_strands;
	}

	bundle_is_excluded = (int*)calloc(max_bundle_i+1, sizeof(int)); //If the bundle has a strand that does not pass though the inclusion sphere it is excluded.


	printf("\nDrawing ROIs...\n");
	fflush(stdout);

	combined_masks = draw_rois(&c, num_voxels, voxel_size, extra_point_radius, point_depth_l_bound, point_depth_u_bound, voxel_radial_l_bound, voxel_radial_u_bound, inclusion_radius, null_mask_l_bound, null_mask_u_bound, subbundle_mask, cumulative_num_strands, bundle_is_excluded);

	num_elems = num_voxels * num_voxels * num_voxels;

	path = (char*)malloc(sizeof(char) * (strlen(output_dir_path) + 200));


	sprintf(path, "%s%crois.hdr", output_dir_path, DIR_SEP);

	write_analyze_header(num_voxels_pdim, voxel_size_pdim, ANALYZE_INT, sizeof(int), path);

	sprintf(path, "%s%crois.img", output_dir_path, DIR_SEP);

	if ( (file =fopen(path,"wb")) == NULL) {
		printf("\nCould not open file %s.!\n", path);
		exit(EXIT_FAILURE);
	}

	fwrite(combined_masks, sizeof(int), num_elems, file);

	fclose(file);


	sprintf(path, "%s%crois_info.txt", output_dir_path, DIR_SEP);

	if ( (file =fopen(path,"w")) == NULL) {
		printf("\nCould not open file %s.!\n", path);
		exit(EXIT_FAILURE);
	}


	fprintf(file, "Bundle index\tNum. strands\tCross-sect. area\tAvg. strand length\tStart ROI id\tEnd ROI id\n");
  fprintf(file, "----------------------------------------------------------------------------------------------------------------\n");


	for (i = 0; i < c.num_bundles; i++) {

		bundle = &(c.bundles[i]);

		if (!(bundle_is_excluded[bundle->bundle_i])) {

			cross_sectional_area = calculate_cross_sectional_area(bundle);
			average_strand_length = calculate_average_strand_length(bundle);

			start_roi_id = get_start_roi_id(bundle->bundle_i, cumulative_num_strands, bundle->num_strands, 0, subbundle_mask);
			end_roi_id = get_end_roi_id(bundle->bundle_i, cumulative_num_strands, bundle->num_strands, 0, subbundle_mask);

			fprintf(file, "%d\t\t%d\t\t%10g\t\t%10g\t\t%d\t\t%d\n", bundle->bundle_i, bundle->num_strands, cross_sectional_area, average_strand_length, start_roi_id, end_roi_id);
		}

	}

	fclose(file);


	printf("Wrote mask to file %s.\n", path);
	fflush(stdout);


	free(combined_masks);


	/* Copy parameter files to new directory for record of what has been performed */

	stage_num = copy_parameters(input_dir_path, output_dir_path);

	if (stage_num != -1) {


		sprintf(path, "%s%cparameters%cdraw_rois_param.txt", output_dir_path, DIR_SEP, DIR_SEP);

		if ( (file =fopen(path,"w")) == NULL) {
			printf("\nCould not open parameters save file %s.!\n", path);
			exit(EXIT_SUCCESS);
		}


		fprintf(file, "num_voxels %d\n", num_voxels);
		fprintf(file, "voxel_size %g\n", voxel_size);
		fprintf(file, "extra_point_radius %g\n", extra_point_radius);
		fprintf(file, "point_depth_l_bound %d\n", point_depth_l_bound);
		fprintf(file, "point_depth_u_bound %d\n", point_depth_u_bound);
		fprintf(file, "voxel_radial_l_bound %g\n", voxel_radial_l_bound);
		fprintf(file, "voxel_radial_u_bound %g\n", voxel_radial_u_bound);
		fprintf(file, "inclusion_radius %g\n", inclusion_radius);
		fprintf(file, "null_mask_l_bound %g\n", null_mask_l_bound);
		fprintf(file, "null_mask_u_bound %g\n", null_mask_u_bound);
		fprintf(file, "subbundle_mask %d\n", subbundle_mask);
		fprintf(file, "save_combined_mask %d\n", save_combined_mask);

		if (fclose(file)) {
			printf("Error closing parameters file copy\n");
		}
	}

	free(path);
	free(cumulative_num_strands);
	free(bundle_is_excluded);
	exit(EXIT_SUCCESS);

}


