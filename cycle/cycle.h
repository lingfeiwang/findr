/* Copyright 2016, 2017 Lingfei Wang
 * 
 * This file is part of Findr.
 * 
 * Findr is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Findr is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 * 
 * You should have received a copy of the GNU Affero General Public License
 * along with Findr.  If not, see <http://www.gnu.org/licenses/>.
 */
/* This lib contains the general definitions of cycle detection routines.
 */
 
#ifndef _HEADER_LIB_CYCLE_H_
#define _HEADER_LIB_CYCLE_H_
#include "../base/config.h"
#include <stdlib.h>
#include "vg.h"

#ifdef __cplusplus
extern "C"
{
#endif

#define CYCLEF(X)	cycle_vg_ ## X

#ifdef __cplusplus
}
#endif
#endif
