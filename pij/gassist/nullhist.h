/* Copyright 2016-2018, 2020 Lingfei Wang
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
/* This file provides the analytical method to construct histograms
 * for the null pdf of LLRs of genotype assisted pij inference.
 */

#ifndef _HEADER_LIB_PIJ_GASSIST_NULLHIST_H_
#define _HEADER_LIB_PIJ_GASSIST_NULLHIST_H_
#include "../../base/config.h"
#include "../../base/gsl/histogram.h"
#include "../../base/types.h"
#ifdef __cplusplus
extern "C"
{
#endif

/* Produce null histograms for all tests (2 to 5).
 * h:		Output location of null histograms. 0 to 3 for tests 2 to 5.
 * nt:		Number of targets
 * ns:		Number of samples
 * nv:		Number of values, = number of alleles + 1
 * dmax:	Maximum value of all LLRs, for histogram construction.
 * 			It can be larger than the maximum of d, if memlimit is not infinite.
 * 			0 to 4 for tests 1 to 5.
 * Return:	0 on success and 1 otherwise
 */
int pij_gassist_nullhists(gsl_histogram** h[4],size_t nt,size_t ns,size_t nv,const FTYPE dmax[4]);





















#ifdef __cplusplus
}
#endif
#endif
