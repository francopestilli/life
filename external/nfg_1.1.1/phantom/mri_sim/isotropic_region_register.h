/*
 *  isotropic_region_register.h
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

#ifndef ISOTROPIC_REGION_REGISTER_H
#define ISOTROPIC_REGION_REGISTER_H

#include "phantom/shared/isotropic_region.h"

typedef struct _isotropic_region_register {

	Isotropic_region *isotropic_region;
	struct _isotropic_region_register *next;
	struct _isotropic_region_register *prev;

} Isotropic_region_register;


Isotropic_region_register* isotropic_region_register_alloc(Isotropic_region *isotropic_region);

void isotropic_region_register_free(Isotropic_region_register *isotropic_region_register);

#endif
