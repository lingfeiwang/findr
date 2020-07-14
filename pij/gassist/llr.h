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
/* This part contains the log likelihood ratio calculations.
 */

#ifndef _HEADER_LIB_PIJ_GASSIST_LLR_H_
#define _HEADER_LIB_PIJ_GASSIST_LLR_H_
#include "../../base/config.h"
#include "../../base/types.h"
#ifdef __cplusplus
extern "C"
{
#endif

/* Multithread calculation of log likelihood ratios for 5 tests.
 * g:		MATRIXF (ng,ns) Full genotype data matrix
 * t:		MATRIXF (ng,ns) Supernormalized transcript data matrix of A
 * t2:		MATRIXF (nt,ns) Supernormalized transcript data matrix of B
 * llr1:	VECTORF (ng). Log likelihood ratios for test 1. Tests E->A v.s. E  A.
 * llr2:	MATRIXF (ng,nt). Log likelihood ratios for test 2. Tests E->B v.s. E  B.
 * llr3:	MATRIXF (ng,nt). Log likelihood ratios for test 3. Tests E->A->B v.s. E->A->B with E->B.
 * llr4:	MATRIXF (ng,nt). Log likelihood ratios for test 4. Tests E->A->B with E->B v.s. E->A  B.
 * llr5:	MATRIXF (ng,nt). Log likelihood ratios for test 5. Tests E->A->B with E->B v.s. A<-E->B.
 * nv:		Number of possible values for each genotype
 * Return:	0 on success.
 */
int pij_gassist_llr(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,VECTORF* llr1,MATRIXF* llr2,MATRIXF* llr3,MATRIXF* llr4,MATRIXF* llr5,size_t nv);































#ifdef __cplusplus
}
#endif
#endif
