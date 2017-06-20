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
#include <assert.h>
#include "../../base/logger.h"
#include "../../base/macros.h"
#include "../../base/const.h"
#include "../../base/supernormalize.h"
#include "../../base/threading.h"
#include "llr.h"
#include "llrtopij.h"
#include "llrtopv.h"
#include "cassist.h"


int pijs_cassist_pv(const MATRIXF* g,const MATRIXF* t,const MATRIXF* t2,VECTORF* p1,MATRIXF* p2,MATRIXF* p3,MATRIXF* p4,MATRIXF* p5,size_t memlimit)
{
#define	CLEANUP			CLEANMATF(gnew)CLEANMATF(tnew)CLEANMATF(tnew2)
	MATRIXF		*gnew;			//(ng,ns) Supernormalized transcript matrix
	MATRIXF		*tnew,*tnew2;	//(nt,ns) Supernormalized transcript matrix
	int			ret;
	size_t		ns=g->size2;
#ifndef NDEBUG
	size_t		nt;
	size_t		ng=g->size1;
	
	nt=t2->size1;
	ns=g->size2;
#endif

	gnew=tnew=tnew2=0;

	//Validation
	assert(!((t->size1!=ng)||(t->size2!=ns)||(t2->size2!=ns)
		||(p1&&(p1->size!=ng))
		||(p2&&((p2->size1!=ng)||(p2->size2!=nt)))
		||(p3&&((p3->size1!=ng)||(p3->size2!=nt)))
		||(p4&&((p4->size1!=ng)||(p4->size2!=nt)))
		||(p5&&((p5->size1!=ng)||(p5->size2!=nt)))));
	assert(memlimit);
	
	if(ns<4)
		ERRRET("Cannot compute p-values with fewer than 4 samples.")

	{
		size_t mem1;
		mem1=(4*t->size1*t->size2+2*t2->size1*t2->size2+p1->size+p2->size1*p2->size2*4)*sizeof(FTYPE);
		if(memlimit<=mem1)
			ERRRET("Memory limit lower than minimum memory needed. Try increasing your memory usage limit.")
		LOG(10,"Memory limit: %lu bytes.",memlimit)
	}
	
	gnew=MATRIXFF(alloc)(g->size1,g->size2);
	tnew=MATRIXFF(alloc)(t->size1,t->size2);
	tnew2=MATRIXFF(alloc)(t2->size1,t2->size2);
	if(!(gnew&&tnew&&tnew2))
		ERRRET("Not enough memory.")

	//Step 1: Supernormalization
	LOG(9,"Supernormalizing...")
	MATRIXFF(memcpy)(gnew,g);
	ret=supernormalizea_byrow(gnew);
	MATRIXFF(memcpy)(tnew,t);
	ret=ret||supernormalizea_byrow(tnew);
	MATRIXFF(memcpy)(tnew2,t2);
	ret=ret||supernormalizea_byrow(tnew2);
	if(ret)
		ERRRET("Supernormalization failed.")

	//Step 2: Log likelihood ratios from nonpermuted data
	LOG(9,"Calculating real log likelihood ratios...")
	pij_cassist_llr(gnew,tnew,tnew2,p1,p2,p3,p4,p5);
	//Step 3: Convert log likelihood ratios to p-values
	LOG(9,"Converting log likelihood ratios into p-values...")
	pij_cassist_llrtopvs(p1,p2,p3,p4,p5,ns);
	//Cleanup
	CLEANUP
	return ret;
#undef	CLEANUP		
}

int pijs_cassist_a(const MATRIXF* g,const MATRIXF* t,const MATRIXF* t2,VECTORF* p1,MATRIXF* p2,MATRIXF* p3,MATRIXF* p4,MATRIXF* p5,char nodiag,size_t memlimit)
{
#define	CLEANUP			CLEANMATF(gnew)CLEANMATF(tnew)CLEANMATF(tnew2)
	MATRIXF		*gnew;			//(ng,ns) Supernormalized transcript matrix
	MATRIXF		*tnew,*tnew2;	//(nt,ns) Supernormalized transcript matrix
	VECTORFF(view)	vv;
	int			ret;
	size_t		ns=g->size2;
#ifndef NDEBUG
	size_t		nt;
	size_t		ng=g->size1;
	
	nt=t2->size1;
	ns=g->size2;
#endif

	gnew=tnew=tnew2=0;

	//Validation
	assert(!((t->size1!=ng)||(t->size2!=ns)||(t2->size2!=ns)
		||(p1&&(p1->size!=ng))
		||(p2&&((p2->size1!=ng)||(p2->size2!=nt)))
		||(p3&&((p3->size1!=ng)||(p3->size2!=nt)))
		||(p4&&((p4->size1!=ng)||(p4->size2!=nt)))
		||(p5&&((p5->size1!=ng)||(p5->size2!=nt)))));
	assert(memlimit);
	
	if(ns<4)
		ERRRET("Cannot compute probabilities with fewer than 4 samples.")
	//Defaults to 8GB memory usage
	{
		size_t mem1;
		mem1=(4*t->size1*t->size2+2*t2->size1*t2->size2+p1->size+p2->size1*p2->size2*4)*sizeof(FTYPE);
		if(memlimit<=mem1)
			ERRRET("Memory limit lower than minimum memory needed. Try increasing your memory usage limit.")
		LOG(10,"Memory limit: %lu bytes.",memlimit)
	}
	
	gnew=MATRIXFF(alloc)(g->size1,g->size2);
	tnew=MATRIXFF(alloc)(t->size1,t->size2);
	tnew2=MATRIXFF(alloc)(t2->size1,t2->size2);
	if(!(gnew&&tnew&&tnew2))
		ERRRET("Not enough memory.")

	//Step 1: Supernormalization
	LOG(9,"Supernormalizing...")
	MATRIXFF(memcpy)(gnew,g);
	ret=supernormalizea_byrow(gnew);
	MATRIXFF(memcpy)(tnew,t);
	ret=ret||supernormalizea_byrow(tnew);
	MATRIXFF(memcpy)(tnew2,t2);
	ret=ret||supernormalizea_byrow(tnew2);
	if(ret)
		ERRRET("Supernormalization failed.")

	//Step 2: Log likelihood ratios from nonpermuted data
	LOG(9,"Calculating real log likelihood ratios...")
	pij_cassist_llr(gnew,tnew,tnew2,p1,p2,p3,p4,p5);
	//Step 3: Convert log likelihood ratios to probabilities
	if((ret=pij_cassist_llrtopijs_a(p1,p2,p3,p4,p5,ns,nodiag)))
		LOG(4,"Failed to convert all log likelihood ratios to probabilities.")
	if(nodiag)
	{
		vv=MATRIXFF(diagonal)(p2);
		VECTORFF(set_zero)(&vv.vector);
		vv=MATRIXFF(diagonal)(p3);
		VECTORFF(set_zero)(&vv.vector);
		vv=MATRIXFF(diagonal)(p4);
		VECTORFF(set_zero)(&vv.vector);
		vv=MATRIXFF(diagonal)(p5);
		VECTORFF(set_zero)(&vv.vector);
	}
	
	//Cleanup
	CLEANUP
	return ret;
#undef	CLEANUP		
}

int pijs_cassist(const MATRIXF* g,const MATRIXF* t,const MATRIXF* t2,VECTORF* p1,MATRIXF* p2,MATRIXF* p3,MATRIXF* p4,MATRIXF* p5,char nodiag,size_t memlimit)
{
	return pijs_cassist_a(g,t,t2,p1,p2,p3,p4,p5,nodiag,memlimit);
}

/* Estimates the probability p(E->A->B) from genotype and expression data. Combines results
 * from any pij_cassist_pijs. For more information, see pij_cassist_pijs_ab.
 * ans:	(ng,nt) Output matrix for probabilities. ans[A,B] is p(E->A->B).
 * pijs:	Function to calculate pijs.
 */
static int pij_cassist_any(const MATRIXF* g,const MATRIXF* t,const MATRIXF* t2,MATRIXF* ans,char nodiag,int (*pijs)(const MATRIXF*,const MATRIXF*,const MATRIXF*,VECTORF*,MATRIXF*,MATRIXF*,MATRIXF*,MATRIXF*,char,size_t),size_t memlimit)
{
#define	CLEANUP			CLEANVECF(p1)CLEANMATF(p2)CLEANMATF(p3)CLEANMATF(p4)
	VECTORF	*p1;
	MATRIXF	*p2,*p3,*p4;
	size_t	ng=g->size1;
	size_t	nt=t2->size1;

	assert(g&&t&&t2&&ans&&pijs);
	assert((g->size2==t->size2)&&(g->size2==t2->size2));
	assert((t->size1==ng)&&(ans->size1==ng)&&(ans->size2==nt));
	p1=VECTORFF(alloc)(ng);
	p2=MATRIXFF(alloc)(ng,nt);
	p3=MATRIXFF(alloc)(ng,nt);
	p4=MATRIXFF(alloc)(ng,nt);
	if(!(p1&&p2&&p3&&p4))
		ERRRET("Not enough memory.")
	if(pijs(g,t,t2,p1,p2,p3,p4,ans,nodiag,memlimit))
		ERRRET("pij_cassist_pijs failed.")
		
	//Combine tests
	#pragma omp parallel
	{
		size_t	ng1,ng2;
		MATRIXFF(view)	mva,mv2,mv4;
		threading_get_startend(g->size1,&ng1,&ng2);
		if(ng1<ng2)
		{
			mva=MATRIXFF(submatrix)(ans,ng1,0,ng2-ng1,ans->size2);
			mv2=MATRIXFF(submatrix)(p2,ng1,0,ng2-ng1,p2->size2);
			mv4=MATRIXFF(submatrix)(p4,ng1,0,ng2-ng1,p4->size2);
			MATRIXFF(mul_elements)(&mva.matrix,&mv2.matrix);
			MATRIXFF(add)(&mva.matrix,&mv4.matrix);
			MATRIXFF(scale)(&mva.matrix,0.5);
		}
	}	
	//Cleanup
	CLEANUP
	return 0;
#undef	CLEANUP
}

int pij_cassist_a(const MATRIXF* g,const MATRIXF* t,const MATRIXF* t2,MATRIXF* ans,char nodiag,size_t memlimit)
{
	return pij_cassist_any(g,t,t2,ans,nodiag,pijs_cassist_a,memlimit);
}

int pij_cassist(const MATRIXF* g,const MATRIXF* t,const MATRIXF* t2,MATRIXF* ans,char nodiag,size_t memlimit)
{
	return pij_cassist_a(g,t,t2,ans,nodiag,memlimit);
}

int pij_cassist_trad(const MATRIXF* g,const MATRIXF* t,const MATRIXF* t2,MATRIXF* ans,char nodiag,size_t memlimit)
{
#define	CLEANUP			CLEANVECF(p1)CLEANMATF(p2)CLEANMATF(p4)CLEANMATF(p5)
	VECTORF	*p1;
	MATRIXF	*p2,*p4,*p5;
	size_t	ng=g->size1;
	size_t	nt=t2->size1;

	assert(g&&t&&t2&&ans);
	assert((g->size2==t->size2)&&(g->size2==t2->size2));
	assert((t->size1==ng)&&(ans->size1==ng)&&(ans->size2==nt));
	p1=VECTORFF(alloc)(ng);
	p2=MATRIXFF(alloc)(ng,nt);
	p4=MATRIXFF(alloc)(ng,nt);
	p5=MATRIXFF(alloc)(ng,nt);
	if(!(p1&&p2&&p5&&p4))
		ERRRET("Not enough memory.")
	if(pijs_cassist_a(g,t,t2,p1,p2,ans,p4,p5,nodiag,memlimit))
		ERRRET("pij_cassist_pijs failed.")
		
	//Combine tests
	MATRIXFF(mul_elements)(ans,p2);
	
	//Cleanup
	CLEANUP
	return 0;
#undef	CLEANUP
}















