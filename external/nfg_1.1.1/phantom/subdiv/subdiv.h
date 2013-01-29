/*
 *  subdiv.h
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
 
#ifndef SUBDIVIDE_H
#define SUBDIVIDE_H


#include "phantom/shared/strand_collection.h"
#include "phantom/shared/strand.h"
 
#define STARTING_ANGLE_THRESHOLD 0.001
#define LINEAR_DEPEND_THRESHOLD 0.95

#define STRAND_R_FINAL_DEFAULT 0.02
 
void subdivide_collection(Strand_collection *children, Strand_collection *parents, double strand_r_final);

void subdivide_strand(Strand *parent, Strand_collection *children, double strand_r_final); 

#endif

