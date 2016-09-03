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

/* Convert LLR into probabilities per A. Uses pij_llrtopij_a_convert.
 * d:		(ng,nx) Source real LLRs to compare with null LLRs
 * dconv:	(ng,nt) LLRs to convert to probabilities. Can differ from d
 * ans:		(ng,nt) Output location of converted probabilities.
 * ns:		Number of samples, used for calculation of null distribution
 * nodiag:	Whether diagonal elements of d should be ignored when converting
 * 			to probabilities
 * Return:	0 if succeed.
 */
int pij_rank_llrtopij_a(const MATRIXF* d,const MATRIXF* dconv,MATRIXF* ans,size_t ns,char nodiag);

/* Calculate probabilities of A--B based on LLR distributions of real data 
 * and null hypothesis. Conversion from LLR to probability is per A.
 * See pij_rank_pij_any.
 */
int pij_rank_a(const MATRIXF* t,const MATRIXF* t2,MATRIXF* p,char nodiag);
















#ifdef __cplusplus
}
#endif
#endif
