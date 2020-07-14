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
#include "../../base/config.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "../../base/logger.h"
#include "llrtopij.h"
#include "../llrtopij.h"






/* Always return probability of step 1 is 1. This is useful when best eQTL are already selected in advance.
 */
static inline int pij_cassist_llrtopij1_1(VECTORF* p1)
{
	LOG(9,"Converting LLR to probabilities for step 1. Filling with 1.")
	VECTORFF(set_all)(p1,1);
	return 0;
}

/* Functions to convert LLR of specific steps into probabilities.
 * Uses pij_llrtopij_convert with different settings of n1d and n2d.
 * Function name suffices indicate which LLR to convert.
 */
static inline int pij_cassist_llrtopij1(VECTORF* d)
{
	LOG(9,"Converting LLR to probabilities for step 1 on per A basis.")
	return pij_cassist_llrtopij1_1(d);
}

static inline int pij_cassist_llrtopij2(MATRIXF* d,size_t ns,char nodiag)
{
	LOG(9,"Converting LLR to probabilities for step 2 on per A basis.")
	assert(ns>2);
	return pij_llrtopij_convert_single_self(d,1,ns-2,nodiag,0);
}

static inline int pij_cassist_llrtopij3(MATRIXF* d,size_t ns,char nodiag)
{
	LOG(9,"Converting LLR to probabilities for step 3 on per A basis.")
	assert(ns>3);
	if(pij_llrtopij_convert_single_self(d,1,ns-3,nodiag,0))
		return 1;
	MATRIXFF(scale)(d,-1);
	MATRIXFF(add_constant)(d,1);
	return 0;
}

static inline int pij_cassist_llrtopij4(MATRIXF* d,size_t ns,char nodiag)
{
	LOG(9,"Converting LLR to probabilities for step 4 on per A basis.")
	assert(ns>3);
	return pij_llrtopij_convert_single_self(d,2,ns-3,nodiag,0);
}

static inline int pij_cassist_llrtopij5(MATRIXF* d,size_t ns,char nodiag)
{
	LOG(9,"Converting LLR to probabilities for step 5 on per A basis.")
	assert(ns>3);
	return pij_llrtopij_convert_single_self(d,1,ns-3,nodiag,0);
}


int pij_cassist_llrtopijs(VECTORF* p1,MATRIXF* p2,MATRIXF* p3,MATRIXF* p4,MATRIXF* p5,size_t ns,char nodiag)
{
	int	ret=0,ret2=0;
	
	if(ns<4)
	{
		LOG(0,"Cannot convert log likelihood ratios to probabilities. Needs at least 4 samples.")
		return 1;
	}
	ret=ret||(ret2=pij_cassist_llrtopij2(p2,ns,nodiag));
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
	ret=ret||(ret2=pij_cassist_llrtopij3(p3,ns,nodiag));
	if(ret2)
		LOG(1,"Failed to convert log likelihood ratios to probabilities in step 3.")
	ret=ret||(ret2=pij_cassist_llrtopij4(p4,ns,nodiag));
	if(ret2)
		LOG(1,"Failed to convert log likelihood ratios to probabilities in step 4.")
	ret=ret||(ret2=pij_cassist_llrtopij5(p5,ns,nodiag));
	if(ret2)
		LOG(1,"Failed to convert log likelihood ratios to probabilities in step 5.")
	return ret;
}













