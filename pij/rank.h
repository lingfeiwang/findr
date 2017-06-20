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
/* This file contains rank-based pij prediction without genotype information.
 * Input expression data are first supernormalized so only rank information
 * remains. Then different prediction method can be applied for pij.
 *
 * Currently only one method is provided. It first calculates the log likelihood
 * ratio (LLR) between null A   B and alternative A---B hypotheses. The LLR
 * is then converted into probability of alternative hypothesis per A.
 * The probability is regarded as pij. This is in function pij_rank_a.
 */

#ifndef _HEADER_LIB_PIJ_RANK_H_
#define _HEADER_LIB_PIJ_RANK_H_
#include "../base/config.h"
#include "../base/types.h"
#ifdef __cplusplus
extern "C"
{
#endif

/* Calculate p-values of A  B against A--B based on LLR distributions of real data 
 * and null hypothesis.
 * t:		(ng,ns) Expression data for A
 * t2:		(nt,ns) Expression data for B
 * p:		(ng,nt) Output for p-values of A--B is false
 * memlimit:Specifies approximate memory usage. Function can fail if memlimit is too small. For unlimited memory, set memlimit=-1.	
 * Return:	0 if succeed.
 */
int pij_rank_pv(const MATRIXF* t,const MATRIXF* t2,MATRIXF* p,size_t memlimit);

/* Calculate probabilities of A--B based on LLR distributions of real data 
 * and null hypothesis. Conversion from LLR to probability is per A.
 * See pij_rank_pij_any.
 */
int pij_rank_a(const MATRIXF* t,const MATRIXF* t2,MATRIXF* p,char nodiag,size_t memlimit);
//Duplicate with a different name
int pij_rank(const MATRIXF* t,const MATRIXF* t2,MATRIXF* p,char nodiag,size_t memlimit);
















#ifdef __cplusplus
}
#endif
#endif
