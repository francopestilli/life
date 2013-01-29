/*
 *  resample.h
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
 
#ifndef RESAMPLE_H
#define RESAMPLE_H

#include "phantom/shared/strand_collection.h"
#include "phantom/shared/strand.h"


#define RESAMPLE_LENGTH_DEFAULT 0.02
#define DOUBLE_BACK_ANGLE_THRESHOLD_DEFAULT 150
#define FORWARD_ANGLE_THRESHOLD_DEFAULT 80 
 
void resample_collection(Strand_collection *old_c, Strand_collection *resample_c, double resample_length, double double_back_angle_threshold, double forward_angle_threshold);

void resample_strand(Strand *old_strand, Strand_collection *resample_c, double resample_length, double double_back_angle_threshold, double forward_angle_threshold); 

#endif

