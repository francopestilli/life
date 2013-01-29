/*
 *  noisify.c
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
 
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <sys/stat.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "phantom/shared/shared.h"
#include "phantom/noisify/noisify.h"




void noisify(float *image, int image_size, double noise_level, gsl_rng *rand_gen) {

	double noise[2];
	int image_i;
	
	for (image_i = 0; image_i < image_size; image_i++) {
		noise[X] = (double)image[image_i] + gsl_ran_gaussian(rand_gen, noise_level);
		noise[Y] = gsl_ran_gaussian(rand_gen, noise_level);		

		image[image_i] = (float)sqrt(noise[X] * noise[X] + noise[Y] * noise[Y]);
		
		
	}
	
 
}
