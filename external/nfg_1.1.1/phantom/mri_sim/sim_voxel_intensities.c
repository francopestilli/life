/*
 *  sim_voxel_intensities.c
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


#include <math.h>
#include <stdlib.h>
#include <stdio.h>


#include "phantom/mri_sim/voxel.h"
#include "phantom/mri_sim/subvoxel.h"
#include "phantom/mri_sim/sim_voxel_intensities.h"


/* Override this function to change the way the DW-MR signal is simulated.  'voxel' is a struct defined in 'Voxel.h'.  The imaging gradient directions are defined in cartesian coordinates.  fa, b_value and diffusivity are all self explanatory. */
double* sim_voxel_intensities(Voxel *voxel, int num_grad_directions, double *grad_directions,  double *b_values, double fa, double diffusivity) {
	
	double *intensities, tensor[3][3], isotropic_diffusivity, baseline_signal;
	int x, y, z, grad_i, offset;
	Subvoxel *subvoxel;

	intensities  = (double*)calloc(sizeof(double), num_grad_directions);
	
	for (z = 0; z < voxel->num_subvoxels[Z]; z++) {
		for (y = 0; y < voxel->num_subvoxels[Y]; y++) {
			for (x = 0; x < voxel->num_subvoxels[X]; x++) {
			
				offset = z * voxel->num_subvoxels[Y] * voxel->num_subvoxels[X] + y * voxel->num_subvoxels[X] + x;
				subvoxel = &(voxel->subvoxels[offset]);
				
				if (subvoxel->closest_isotropic_region != NULL) {
					grad_i = 1;
				}
				
			
				if (subvoxel->closest_segment != NULL) {

					create_tensor(tensor, fa, diffusivity, subvoxel->orientation);	
						
					for (grad_i = 0; grad_i < num_grad_directions; grad_i++) {
						intensities[grad_i] += sample_tensor(tensor, b_values[grad_i], &(grad_directions[grad_i * 3]));
					}
				
				} else if (subvoxel->closest_isotropic_region != NULL) {		
						
					/* Defined isotropic regions can have their own b=0 signal (defined relative to the white matter b=0 signal) and diffusivity*/
					isotropic_diffusivity = subvoxel->closest_isotropic_region->diffusivity;
					baseline_signal = subvoxel->closest_isotropic_region->baseline_signal;
							
					for (grad_i = 0; grad_i < num_grad_directions; grad_i++) {
						intensities[grad_i] += baseline_signal * exp(-1 * b_values[grad_i] * isotropic_diffusivity);
					}			
						
				} else { 
							
					/* Voids within the structure take on the same b=0 signal and diffusivity as the white matter strands*/	
					for (grad_i = 0; grad_i < num_grad_directions; grad_i++) {
						intensities[grad_i] += exp(-1 * b_values[grad_i] * diffusivity);
					}			
										
				}
			
			}
		}
	}
	
	for (grad_i = 0; grad_i < num_grad_directions; grad_i++) {
		intensities[grad_i] /= voxel->num_subvoxels[X] * voxel->num_subvoxels[Y] * voxel->num_subvoxels[Z];
	}
	
	return intensities;
	
}

void create_tensor(double tensor[3][3], double fa, double diffusivity, double eig_v[3]) {
	

	double a, dr, dp;
/*	
	if (vector_norm(eig_v) == 0.0) {				//Is an empty subvoxel -> fill with isotropic tensor
	
		tensor[X][X] = diffusivity;
		tensor[X][Y] = 0.0;
		tensor[X][Z] = 0.0;
		tensor[Y][X] = 0.0;
		tensor[Y][Y] = diffusivity;
		tensor[Y][Z] = 0.0;
		tensor[Z][X] = 0.0;
		tensor[Z][Y] = 0.0;
		tensor[Z][Z] = diffusivity;
	
	} else {									//Is an nonempty subvoxel -> fill with anisotropic tensor
*/	

		a = fa / sqrt(3-2 * fa * fa);
		dr = diffusivity * (1.0 - a);
		dp = diffusivity * (1.0 + 2.0 * a);
		
		
		tensor[X][X] = dr + (dp - dr) * eig_v[X] * eig_v[X];
		tensor[Y][Y] = dr + (dp - dr) * eig_v[Y] * eig_v[Y];
		tensor[Z][Z] = dr + (dp - dr) * eig_v[Z] * eig_v[Z];
		tensor[X][Y] = tensor[Y][X] = (dp - dr) * eig_v[X] * eig_v[Y];
		tensor[X][Z] = tensor[Z][X] = (dp - dr) * eig_v[X] * eig_v[Z];
		tensor[Y][Z] = tensor[Z][Y] = (dp - dr) * eig_v[Y] * eig_v[Z];
		
/*	}*/

}


double sample_tensor(double tensor[3][3], double b_value, double *grad_direction) {

	double intensity;

	intensity = exp(-b_value * ( 
						grad_direction[X] * grad_direction[X] * tensor[X][X] +
						grad_direction[Y] * grad_direction[Y] * tensor[Y][Y] +
						grad_direction[Z] * grad_direction[Z] * tensor[Z][Z] +
						2.0 * grad_direction[X] * grad_direction[Y] * tensor[X][Y] +
						2.0 * grad_direction[X] * grad_direction[Z] * tensor[X][Z] +
						2.0 * grad_direction[Y] * grad_direction[Z] * tensor[Y][Z]));   
						
	return intensity;

}

