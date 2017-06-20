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
#include "../../base/config.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "../../base/logger.h"
#include "../llrtopij.h"
#include "llrtopij_a.h"
#include "llrtopij.h"
#pragma GCC diagnostic ignored "-Wunused-parameter"

/* Functions to convert LLR of specific steps into probabilities.
 * Uses pij_llrtopij_a_convert with different settings of n1d and n2d.
 * Function name suffices indicate which LLR to convert.
 */
static inline int pij_cassist_llrtopij1_a(VECTORF* d)
{
	LOG(9,"Converting LLR to probabilities for step 1 on per A basis.")
	return pij_cassist_llrtopij1_1(d);
}

static inline int pij_cassist_llrtopij2_a(MATRIXF* d,size_t ns,char nodiag)
{
	LOG(9,"Converting LLR to probabilities for step 2 on per A basis.")
	assert(ns>2);
	return pij_llrtopij_a_convert_single_self(d,1,ns-2,nodiag,0);
}

static inline int pij_cassist_llrtopij3_a(MATRIXF* d,size_t ns,char nodiag)
{
	LOG(9,"Converting LLR to probabilities for step 3 on per A basis.")
	assert(ns>3);
	if(pij_llrtopij_a_convert_single_self(d,1,ns-3,nodiag,0))
		return 1;
	MATRIXFF(scale)(d,-1);
	MATRIXFF(add_constant)(d,1);
	return 0;
}

static inline int pij_cassist_llrtopij4_a(MATRIXF* d,size_t ns,char nodiag)
{
	LOG(9,"Converting LLR to probabilities for step 4 on per A basis.")
	assert(ns>3);
	return pij_llrtopij_a_convert_single_self(d,2,ns-3,nodiag,0);
}

static inline int pij_cassist_llrtopij5_a(MATRIXF* d,size_t ns,char nodiag)
{
	LOG(9,"Converting LLR to probabilities for step 5 on per A basis.")
	assert(ns>3);
	return pij_llrtopij_a_convert_single_self(d,1,ns-3,nodiag,0);
}


int pij_cassist_llrtopijs_a(VECTORF* p1,MATRIXF* p2,MATRIXF* p3,MATRIXF* p4,MATRIXF* p5,size_t ns,char nodiag)
{
	int	ret=0,ret2=0;
	
	if(ns<4)
	{
		LOG(0,"Cannot convert log likelihood ratios to probabilities. Needs at least 4 samples.")
		return 1;
	}
	ret=ret||(ret2=pij_cassist_llrtopij2_a(p2,ns,nodiag));
	if(ret2)
		LOG(1,"Failed to convert log likelihood ratios to probabilities in step 2.")
	//For p1, if nodiag, copy p2 data, otherwise set all to 1.
	if(nodiag)
	{
		VECTORFF(view)	vv;
		vv=MATRIXFF(diagonal)(p2);
		ret=ret||(ret2=VECTORFF(memcpy)(p1,&vv.vector));
	}
	else
		ret=ret||(ret2=pij_cassist_llrtopij1_1(p1));
	if(ret2)
		LOG(1,"Failed to convert log likelihood ratios to probabilities in step 1.")
	ret=ret||(ret2=pij_cassist_llrtopij3_a(p3,ns,nodiag));
	if(ret2)
		LOG(1,"Failed to convert log likelihood ratios to probabilities in step 3.")
	ret=ret||(ret2=pij_cassist_llrtopij4_a(p4,ns,nodiag));
	if(ret2)
		LOG(1,"Failed to convert log likelihood ratios to probabilities in step 4.")
	ret=ret||(ret2=pij_cassist_llrtopij5_a(p5,ns,nodiag));
	if(ret2)
		LOG(1,"Failed to convert log likelihood ratios to probabilities in step 5.")
	return ret;
}











