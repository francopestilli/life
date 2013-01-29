/*
 *  point.c
 *  Numerical Fibre Generator
 *
 *  Created by Tom Close, on 25/06/08.
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

#include "phantom/shared/control_point.h"


void set_control_point(Control_point *control_point, double *pos, double *grad) {
	
	control_point->pos = pos;
	control_point->grad = grad;
	
}


void print_control_point(Control_point *control_point, char indent[]) {
	 printf("%sCoord [%g,%g,%g], Grad [%g,%g,%g]\n", indent, control_point->pos[0], control_point->pos[1], control_point->pos[2], control_point->grad[0], control_point->grad[1], control_point->grad[2]);
 }
