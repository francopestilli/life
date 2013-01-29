/*
 *  rand_init.h
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

#ifndef _RAND_INIT_H_
#define _RAND_INIT_H_

#include <gsl/gsl_rng.h>
#include "phantom/shared/strand_collection.h"

#define SPHERE_R_DEFAULT 1.0
#define CONTROL_POINT_FREQ_DEFAULT 50.0
#define STRAND_R_LBOUND_DEFAULT 0.075
#define STRAND_R_UBOUND_DEFAULT 0.275
#define STRAND_R_BUFFER_RATIO_DEFAULT 1.1
#define MAX_ATTEMPTS_DEFAULT 1000000
#define GREY_MATTER_R_LBOUND_DEFAULT 0.15
#define GREY_MATTER_R_UBOUND_DEFAULT 0.6
#define GREY_MATTER_LOCATION_LBOUND_DEFAULT 1.5
#define GREY_MATTER_LOCATION_UBOUND_DEFAULT 2.5
#define NUM_GREY_MATTERS_DEFAULT 0

#define GREY_MATTER_DIFFUSIVITY_DEFAULT 0.0009
#define GREY_MATTER_BASELINE_SIGNAL_DEFAULT 1.0
#define GREY_MATTER_WEIGHTING_DEFAULT 0.5

/* Randomly initialises a collection of strands by generating pairs of points on a sphere then joining them with evenly spaced control points */
int rand_init_collection(Strand_collection *c, double sphere_r, double control_point_freq, double strand_r_lbound, double strand_r_ubound, double strand_r_buffer_ratio, gsl_rng *rand_gen, int max_attempts, int num_isotropic_regions, Isotropic_region *isotropic_regions);

/* As the name suggests, randomly generates a point on a sphere */
void rand_gen_point_on_sphere(double *point, double sphere_r, const gsl_rng *random_gen);


void rand_init_grey_matters(Isotropic_region *isotropic_regions, int num_grey_matters, double grey_matter_r_lbound, double grey_matter_r_ubound, double grey_matter_location_lbound, double grey_matter_location_ubound, double grey_matter_diffusivity, double grey_matter_baseline_signal, double grey_matter_weighting, gsl_rng *rand_gen);
	
#endif
