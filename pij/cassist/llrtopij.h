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
/* This file contains the conversion from log likelihood ratio to probabilities
 *
 */

#ifndef _HEADER_LIB_PIJ_CASSIST_LLRTOPIJ_H_
#define _HEADER_LIB_PIJ_CASSIST_LLRTOPIJ_H_
#include "../../base/config.h"
#include "../../base/types.h"
#include "llrtopij_a.h"
#ifdef __cplusplus
extern "C"
{
#endif


/* Always return probability of step 1 is 1. This is useful when best eQTL are already selected in advance.
 */
int pij_cassist_llrtopij1_1(VECTORF* p1);











#ifdef __cplusplus
}
#endif
#endif
