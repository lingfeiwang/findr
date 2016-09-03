/* Copyright 2016 Lingfei Wang
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

#ifndef _HEADER_LIB_PIJ_GASSIST_LLRTOPIJ_A_H_
#define _HEADER_LIB_PIJ_GASSIST_LLRTOPIJ_A_H_
#include "../../base/config.h"
#include "../../base/types.h"
#include "../llrtopij_a.h"
#ifdef __cplusplus
extern "C"
{
#endif

/* Convert real log likelihood ratios into probability functions.
 * This function converts every A in hypothesis (E->A->B) separately.
 * Suppose there are ng (E,A) pairs and nt Bs, this function converts ng times,
 * each for one (E,A) pair but all Bs.
 * d:		(ng,nt)	Input log likelihood ratios for construction of 
 * 			histograms and calculation of probability of true hypothesis.
 * dconv:	(ng,nx)	Actual log likelihood ratios to be converted to
 * 			probabilities using histogram constructed base on d.
 * ans:		(ng,nx) Output matrix for probabilities.
 * g:		(ng,ns)	Original genotype matrix, used for analytical calculation
 * 			of null distribution. Every element=0,1,...,nv-1.
 * nv:		Maximum number of values each g may take.
 * n1c,
 * n1d,
 * n2c,
 * n2d:		Parameters to specify null distribution. See pij_llrtopij_a_nullhist
 * nodiag:	If diagonal elements of d should be removed in construction of real
 * 			histogram. This should be set to true (!=0) when t is identical with
 * 			the top rows of t2 (in calculation of llr).
 * Return:	0 if success.
 */
int pij_gassist_llrtopij_a_convert(const MATRIXF* d,const MATRIXF* dconv,MATRIXF* ans,const MATRIXG* g,size_t nv,long n1c,size_t n1d,long n2c,size_t n2d,char nodiag);

/* Functions to convert LLR of specific steps into probabilities.
 * Uses pij_llrtopij_a_convert with different settings of n1d and n2d.
 * Function name suffices indicate which LLR to convert.
 */
static inline int pij_gassist_llrtopij1_a(const MATRIXF* d,const MATRIXF* dconv,MATRIXF* ans,const MATRIXG* g,size_t nv,char nodiag);
static inline int pij_gassist_llrtopij2_a(const MATRIXF* d,const MATRIXF* dconv,MATRIXF* ans,const MATRIXG* g,size_t nv,char nodiag);
static inline int pij_gassist_llrtopij3_a(const MATRIXF* d,const MATRIXF* dconv,MATRIXF* ans,const MATRIXG* g,size_t nv,char nodiag);
static inline int pij_gassist_llrtopij4_a(const MATRIXF* d,const MATRIXF* dconv,MATRIXF* ans,const MATRIXG* g,size_t nv,char nodiag);
static inline int pij_gassist_llrtopij5_a(const MATRIXF* d,const MATRIXF* dconv,MATRIXF* ans,const MATRIXG* g,size_t nv,char nodiag);

/* Converts four LLRs into probabilities together.
 * Uses pij_gassist_llrtopij1_a to pij_gassist_llrtopij5_a.
 * See above functions for parameter definitions.
 * Return: 0 if all functions are successful.
 */
int pij_gassist_llrtopijs_a(const MATRIXG* g,const VECTORF* llr1,const MATRIXF* llr2,const MATRIXF* llr3,const MATRIXF* llr4,const MATRIXF* llr5,VECTORF* p1,MATRIXF* p2,MATRIXF* p3,MATRIXF* p4,MATRIXF* p5,size_t nv,char nodiag);


/**********************************************************************
 * Static functions
 **********************************************************************/

static inline int pij_gassist_llrtopij1_a(const MATRIXF* d,const MATRIXF* dconv,MATRIXF* ans,const MATRIXG* g,size_t nv,char nodiag)
{
	LOG(9,"Converting LLR to probabilities for step 1 on per A basis.")
	return pij_gassist_llrtopij_a_convert(d,dconv,ans,g,nv,1,1,1,g->size2-2,nodiag);
}

static inline int pij_gassist_llrtopij2_a(const MATRIXF* d,const MATRIXF* dconv,MATRIXF* ans,const MATRIXG* g,size_t nv,char nodiag)
{
	LOG(9,"Converting LLR to probabilities for step 2 on per A basis.")
	return pij_gassist_llrtopij_a_convert(d,dconv,ans,g,nv,1,1,1,g->size2-2,nodiag);
}

static inline int pij_gassist_llrtopij3_a(const MATRIXF* d,const MATRIXF* dconv,MATRIXF* ans,const MATRIXG* g,size_t nv,char nodiag)
{
	LOG(9,"Converting LLR to probabilities for step 3 on per A basis.")
	if(pij_gassist_llrtopij_a_convert(d,dconv,ans,g,nv,1,1,1,g->size2-3,nodiag))
		return 1;
	MATRIXFF(scale)(ans,-1);
	MATRIXFF(add_constant)(ans,1);
	return 0;
}

static inline int pij_gassist_llrtopij4_a(const MATRIXF* d,const MATRIXF* dconv,MATRIXF* ans,const MATRIXG* g,size_t nv,char nodiag)
{
	LOG(9,"Converting LLR to probabilities for step 4 on per A basis.")
	return pij_gassist_llrtopij_a_convert(d,dconv,ans,g,nv,1,2,1,g->size2-3,nodiag);
}

static inline int pij_gassist_llrtopij5_a(const MATRIXF* d,const MATRIXF* dconv,MATRIXF* ans,const MATRIXG* g,size_t nv,char nodiag)
{
	LOG(9,"Converting LLR to probabilities for step 5 on per A basis.")
	if(pij_gassist_llrtopij_a_convert(d,dconv,ans,g,nv,0,1,1,g->size2-3,nodiag))
		return 1;
	return 0;
}









































#ifdef __cplusplus
}
#endif
#endif
