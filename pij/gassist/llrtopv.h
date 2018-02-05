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
/* This file contains the conversion from log likelihood ratio to p-values for discrete anchors, e.g. genotypes
 *
 */
#ifndef _HEADER_LIB_PIJ_GASSIST_LLRTOPV_H_
#define _HEADER_LIB_PIJ_GASSIST_LLRTOPV_H_
#include "../../base/config.h"
#include "../../base/types.h"
#ifdef __cplusplus
extern "C"
{
#endif

/* Convert log likelihood ratios into p-values for matrix in multi thread
 * p1:		(ng)
 * p2:		(ng,nt)
 * p3:		(ng,nt)
 * p4:		(ng,nt)
 * p5:		(ng,nt)	Data for 5 tests from p1 to p5, as input for log likelihood ratios,
 			and also as output for converted p-values of LLRs of each test.
 * g:		(ng,ns)	Original genotype matrix, used for analytical calculation
 * 			of null distribution. Every element=0,1,...,nv-1.
 			Has matching rows with data p1 to p5..
 * nv:		Maximum number of values each g may take.
 * Return:	0 if success.
 */
int pij_gassist_llrtopvs(VECTORF* p1,MATRIXF* p2,MATRIXF* p3,MATRIXF* p4,MATRIXF* p5,const MATRIXG* g,size_t nv);






































#ifdef __cplusplus
}
#endif
#endif
