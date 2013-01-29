/*
 *  optimise.h
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

#ifndef OPTIMISE_H
#define OPTIMISE_H

#include "phantom/shared/strand_collection.h"
#include "phantom/optimise/sample_block.h"
#include "phantom/optimise/reference_block.h"

#define COST_TERM_WEIGHTS_DEFAULT {1.0e5, 1.0, 20.0}
#define COST_TERM_POWERS_DEFAULT {2.0, 2.0, 2.0}
#define RELATIVE_ISOTROPIC_REPULSION_WEIGHT_DEFAULT 1.0
#define SAMPLE_DENSITY_DEFAULT 100.0
#define GRID_REFERENCE_SIZE_DEFAULT 10
#define INITIAL_STEP_SIZE_DEFAULT 1.0
#define STEP_TOL_DEFAULT  0.1
#define OPTIMIZATION_TOL_DEFAULT 1.0e-2
#define MAX_ITERATIONS_DEFAULT 100
/*#define MAX_NUM_COST_FUNCTION_CALLS_DEFAULT 30*/
#define ITERATIONS_SAVE_FREQ_DEFAULT 0
#define JITTER_STDEV_DEFAULT 0.01
#define JITTER_SEED_DEFAULT 1193125662
#define FOV_DEFAULT 1.2


void optimise(Strand_collection *c, double cost_term_weights[3], double cost_term_powers[3], double sample_density, int grid_reference_size, double initial_step_size, double step_tol, double optimization_tol, double jitter_dev, unsigned int jitter_seed, int max_iterations, int iterations_save_freq, char *iterations_save_path) ;


void allocate_storage_blocks(Sample_block **start_block, Reference_block **grid_reference, int grid_reference_size);

void free_storage_blocks(Sample_block *start_block, Reference_block *grid_reference, int grid_reference_size);


#endif
