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
#include "../base/config.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include "../base/gsl/blas.h"
#include "../base/random.h"
#include "../base/const.h"
#include "../base/logger.h"
#include "../base/macros.h"
#include "../base/data_process.h"
#include "../base/supernormalize.h"
#include "../base/threading.h"
#include "../pij/llrtopij_a.h"
#include "rank.h"

void pij_rank_llr_block(const MATRIXF* t,const MATRIXF* t2,MATRIXF* llr)
{
	size_t	i,j;
	size_t	ng=t->size1;
	size_t	nt=t2->size1;
#ifndef NDEBUG
	size_t	ns=t->size2;
#endif
	assert((t2->size2==ns));
	assert((llr->size1==ng));
	assert((llr->size2==nt));
	MATRIXFF(cov2_bounded)(t,t2,llr);
	MATRIXFF(mul_elements)(llr,llr);
	MATRIXFF(scale)(llr,-1);
	MATRIXFF(add_constant)(llr,1);
	for(i=0;i<ng;i++)
		for(j=0;j<nt;j++)
			MATRIXFF(set)(llr,i,j,(FTYPE)log(MATRIXFF(get)(llr,i,j)));

	MATRIXFF(scale)(llr,-0.5);
	//Bounding from 0
	MATRIXFF(bound_below)(llr,0);
}

void pij_rank_llr(const MATRIXF* t,const MATRIXF* t2,MATRIXF* llr)
{
	assert((t->size2==t2->size2)&&(llr->size1==t->size1)&&(llr->size2==t2->size1));
	#pragma omp parallel
	{
		size_t	n1,n2;
		
		threading_get_startend(t->size1,&n1,&n2);
		if(n2>n1)
		{
			MATRIXFF(const_view) mvt=MATRIXFF(const_submatrix)(t,n1,0,n2-n1,t->size2);
			MATRIXFF(view)	mvllr;
			mvllr=MATRIXFF(submatrix)(llr,n1,0,n2-n1,llr->size2);
			pij_rank_llr_block(&mvt.matrix,t2,&mvllr.matrix);
		}
	}
}

int pij_rank_llrtopij_a(const MATRIXF* d,const MATRIXF* dconv,MATRIXF* ans,size_t ns,char nodiag)
{
	LOG(9,"Converting LLR to probabilities on per A basis.")
	return pij_llrtopij_a_convert_single(d,dconv,ans,1,ns-2,nodiag);
}

int pij_rank_pij_any(const MATRIXF* t,const MATRIXF* t2,MATRIXF* p,char nodiag,int (*pij)(const MATRIXF*,const MATRIXF*,MATRIXF*,size_t,char))
{
#define	CLEANUP			CLEANMATF(tnew)CLEANMATF(tnew2)\
						CLEANMATF(llr)
	MATRIXF		*tnew,*tnew2;			//(nt,ns) Supernormalized transcript matrix
	MATRIXF		*llr;
	VECTORFF(view)	vv;
	int			ret;
	size_t		ng,nt,ns;
	
	ng=t->size1;
	nt=t2->size1;
	ns=t->size2;

	tnew=tnew2=llr=0;
	
	//Validation
	assert((t2->size2==ns)&&(p->size1==ng)&&(p->size2==nt));
	
	llr=MATRIXFF(alloc)(ng,nt);
	tnew=MATRIXFF(alloc)(ng,ns);
	tnew2=MATRIXFF(alloc)(nt,ns);
	if(!(llr&&tnew&&tnew2))
		ERRRET("Not enough memory.")

	//Step 1: Supernormalization
	LOG(9,"Supernormalizing...")
	MATRIXFF(memcpy)(tnew,t);
	ret=supernormalize_byrow(tnew);
	MATRIXFF(memcpy)(tnew2,t2);
	ret=ret||supernormalize_byrow(tnew2);
	if(ret)
		ERRRET("Supernormalization failed.")
	
	//Step 2: Log likelihood ratios from nonpermuted data
	LOG(9,"Calculating real log likelihood ratios...")
	pij_rank_llr(tnew,tnew2,llr);
	if(nodiag)
	{
		vv=MATRIXFF(diagonal)(llr);
		VECTORFF(set_zero)(&vv.vector);
	}
	//Step 3: Convert log likelihood ratios to probabilities
	if(pij(llr,llr,p,ns,nodiag))
		LOG(4,"Failed to convert log likelihood ratios to probabilities.")

	//Cleanup
	CLEANUP
	return ret;
#undef	CLEANUP		
}

int pij_rank_a(const MATRIXF* t,const MATRIXF* t2,MATRIXF* p,char nodiag)
{
	return pij_rank_pij_any(t,t2,p,nodiag,pij_rank_llrtopij_a);
}
















