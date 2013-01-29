/*
 *  sim_voxel_intensities.h
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

#ifndef SIM_VOXEL_INTENSITY_H 
#define SIM_VOXEL_INTENSITY_H

double* sim_voxel_intensities(Voxel *voxel, int num_grad_directions, double *grad_directions, double *b_values, double fa, double diffusivity);

double sample_tensor(double tensor[3][3], double b_value, double *grad_direction);

void create_tensor(double tensor[3][3], double fa, double diffusivity, double eig_v[3]);

#endif
