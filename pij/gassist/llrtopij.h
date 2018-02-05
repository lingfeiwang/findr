/* Copyright 2016-2018 Lingfei Wang
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

#ifndef _HEADER_LIB_PIJ_GASSIST_LLRTOPIJ_H_
#define _HEADER_LIB_PIJ_GASSIST_LLRTOPIJ_H_
#include "../../base/config.h"
#include "../../base/gsl/histogram.h"
#include "../../base/types.h"
#ifdef __cplusplus
extern "C"
{
#endif



/* Converts four LLRs into probabilities together.
 * Uses pij_gassist_llrtopij1 to pij_gassist_llrtopij5.
 * See above functions for parameter definitions.
 * h:		Null histograms. 0 to 3 for tests 2 to 5.
 * Return: 0 if all functions are successful.
 */
int pij_gassist_llrtopijs(const MATRIXG* g,VECTORF* p1,MATRIXF* p2,MATRIXF* p3,MATRIXF* p4,MATRIXF* p5,size_t nv,const gsl_histogram * const * h[4],char nodiag,long nodiagshift);











#ifdef __cplusplus
}
#endif
#endif
