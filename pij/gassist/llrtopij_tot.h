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
 */

#ifndef _HEADER_LIB_PIJ_GASSIST_LLRTOPIJ_TOT_H_
#define _HEADER_LIB_PIJ_GASSIST_LLRTOPIJ_TOT_H_
#include "../../base/config.h"
#include "../../base/gsl/histogram.h"
#include "../../base/types.h"
#include "../llrtopij_tot.h"

#ifdef __cplusplus
extern "C"
{
#endif

/* The recommended wrapper conversion methods of step 1/2/3 data. Uses analytical generation of pdfs of null distribution.
 * g:		(ng,ns)	Genotype data.
 * llr:		Log likelihood ratios of real data
 * p:		Output matrix of converted probabilities
 * nv:		number of possible genotype values for each SNP
 * nratio:	if not NULL, return the ratio of null distribution in real data.
 * Return:	0 on success
 */
int pij_gassist_llrtopij1_tot(const MATRIXG* g,const VECTORF* llr,VECTORF* p,size_t nv,double* nratio);
int pij_gassist_llrtopij2c_tot(const MATRIXG* g,const VECTORF* llr,VECTORF* p,size_t nv,double* nratio);
int pij_gassist_llrtopij2b_tot(const MATRIXG* g,const VECTORF* llr,VECTORF* p,size_t nv,double* nratio);
int pij_gassist_llrtopij3_tot(const MATRIXG* g,const VECTORF* llr,VECTORF* p,size_t nv,double* nratio);


/* Calculates null log likelihood ratios from permuted data, and then convert log likelihood ratios of real data
 * into probabilities. For step 3, permute A and B.
 * g:		(ng,ns)	Full genotype data
 * t:		(ng,ns) Full transcript data of A. Each rows best eQTL must be the same row of g
 * t2:		(nt,ns) Full transcript data of B.
 * llr1:	(ng) Log likelihood ratios of real data for step 1
 * llr2b:	(ng,nt) Log likelihood ratios of real data for step 2 bold
 * llr2c:	(ng,nt) Log likelihood ratios of real data for step 2 conservative
 * llr3:	(ng,nt) Log likelihood ratios of real data for step 3
 * p1:		(ng) Output for converted probabilities of step 1
 * p2b:		(ng,nt) Output for converted probabilities of step 2 bold
 * p2c:		(ng,nt) Output for converted probabilities of step 2 conservative
 * p3:		(ng,nt) Output for converted probabilities of step 3
 * ans:		(ng,nt) Output matrix of converted probabilities
 * nv:		number of possible genotype values for each SNP
 * nodiag:	When the top ng rows of t2 is exactly t, diagonals of p2 and p3 are meaningless.
 *			In this case, set nodiag to 1 to avoid inclusion of NANs. For nodiag=0, t and t2
 *			should not have any identical genes.
 * Return:	0 on success
 */
int pij_gassist_llrtopijs_tot(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,const VECTORF* llr1,const MATRIXF* llr2b,const MATRIXF* llr2c,const MATRIXF* llr3,VECTORF* p1,MATRIXF* p2b,MATRIXF* p2c,MATRIXF* p3,size_t nv,char nodiag);

/* Calculates null log likelihood ratios from permuted data, and then convert log likelihood ratios of real data
 * into probabilities. Combines results from pij_gassist_llrtopijs. For more information, see pij_gassist_llrtopijs.
 * ansb:	(ng,nt) Output for converted probabilities with step 2 bold.
 * ansc:	(ng,nt) Output for converted probabilities with step 2 conservative.
 * Return:	0 on success.
 */
int pij_gassist_llrtopij_tot(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,const VECTORF* llr1,const MATRIXF* llr2b,const MATRIXF* llr2c,const MATRIXF* llr3,MATRIXF* ansb,MATRIXF* ansc,size_t nv,char nodiag);






























#ifdef __cplusplus
}
#endif
#endif
