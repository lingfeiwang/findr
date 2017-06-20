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
/* This file contains the conversion from log likelihood ratio to p-values
 *
 */
#ifndef _HEADER_LIB_PIJ_LLRTOPV_H_
#define _HEADER_LIB_PIJ_LLRTOPV_H_
#include "../base/config.h"
#include "../base/types.h"
#include "nulldist.h"
#ifdef __cplusplus
extern "C"
{
#endif


/* Converts a vector of log likelihood ratios into p-values with the same null distribution.
 * Single thread.
 * For null distribution, see pij_nulldist_cdfQ.
 * p:	data as input for LLR and output for p-values
 * n1,
 * n2:	Null distribution parameters.
 */
static inline void pij_llrtopv_block(VECTORF* p,size_t n1,size_t n2);
// Converts a matrix with the same null distribution in single thread
static inline void pij_llrtopvm_block(MATRIXF* p,size_t n1,size_t n2);
// Converts a matrix with the same null distribution in multi threads
void pij_llrtopvm(MATRIXF* p,size_t n1,size_t n2);



static inline void pij_llrtopv_block(VECTORF* p,size_t n1,size_t n2)
{
	size_t i;
	for(i=0;i<p->size;i++)
		VECTORFF(set)(p,i,(FTYPE)pij_nulldist_cdfQ(VECTORFF(get)(p,i),n1,n2));
}

static inline void pij_llrtopvm_block(MATRIXF* p,size_t n1,size_t n2)
{
	size_t i,j;
	for(i=0;i<p->size1;i++)
		for(j=0;j<p->size2;j++)
			MATRIXFF(set)(p,i,j,(FTYPE)pij_nulldist_cdfQ(MATRIXFF(get)(p,i,j),n1,n2));
}

#ifdef __cplusplus
}
#endif
#endif
