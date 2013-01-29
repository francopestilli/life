/*
 *  trim.h
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

#ifndef TRIM_H
#define TRIM_H

#include "phantom/shared/strand_collection.h"


#define NEW_SPHERE_R_DEFAULT 1.0
#define LENGTH_REJECT_THRESHOLD_DEFAULT 0.2
#define SAVE_SEPERATE_BUNDLES_DEFAULT 0

void trim(Strand_collection *trim_c, Strand_collection *c, double new_sphere_r, double length_reject_threshold, int save_seperate_bundles);


#endif