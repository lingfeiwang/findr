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
/* This file contains the conversion from log likelihood ratio to p-values for continuous anchors
 *
 */

#ifndef _HEADER_LIB_PIJ_CASSIST_LLRTOPV_H_
#define _HEADER_LIB_PIJ_CASSIST_LLRTOPV_H_
#include "../../base/config.h"
#include "../../base/types.h"
#include "../llrtopv.h"
#ifdef __cplusplus
extern "C"
{
#endif

/* Converts log likelihood ratios into p-values for continuous assisted causal inference test for each test separately.
 * d:	MATRIXF of any size, as input of LLRs and also output of corresponding p-values
 * ns:	Number of samples, to be used to calculate the null distribution
 */
static inline void pij_cassist_llrtopv1(VECTORF* d,size_t ns)
{
	assert(ns>3);
	pij_llrtopv_block(d,1,ns-2);
}

static inline void pij_cassist_llrtopv2(MATRIXF* d,size_t ns)
{
	assert(ns>3);
	pij_llrtopvm(d,1,ns-2);
}

static inline void pij_cassist_llrtopv3(MATRIXF* d,size_t ns)
{
	assert(ns>3);
	pij_llrtopvm(d,1,ns-3);
}

static inline void pij_cassist_llrtopv4(MATRIXF* d,size_t ns)
{
	assert(ns>3);
	pij_llrtopvm(d,2,ns-3);
}

static inline void pij_cassist_llrtopv5(MATRIXF* d,size_t ns)
{
	assert(ns>3);
	pij_llrtopvm(d,1,ns-3);
}

/* Converts log likelihood ratios into p-values for continuous assisted causal inference test for all tests together
 * p1:		(ng)
 * p2:		(ng,nt)
 * p3:		(ng,nt)
 * p4:		(ng,nt)
 * p5:		(ng,nt)	Data for 5 tests from p1 to p5, as input for log likelihood ratios,
 			and also as output for converted p-values.
 * ns:		Number of samples.
 */
static inline void pij_cassist_llrtopvs(VECTORF* p1,MATRIXF* p2,MATRIXF* p3,MATRIXF* p4,MATRIXF* p5,size_t ns)
{
	pij_cassist_llrtopv1(p1,ns);
	pij_cassist_llrtopv2(p2,ns);
	pij_cassist_llrtopv3(p3,ns);
	pij_cassist_llrtopv4(p4,ns);
	pij_cassist_llrtopv5(p5,ns);
}






































#ifdef __cplusplus
}
#endif
#endif
