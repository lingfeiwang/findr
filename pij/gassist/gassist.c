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
#include "../../base/config.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "../../base/logger.h"
#include "../../base/macros.h"
#include "../../base/const.h"
#include "../../base/supernormalize.h"
#include "llr.h"
#include "llrtopij.h"
#include "gassist.h"


int pijs_gassist_tot(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,VECTORF* p1,MATRIXF* p2,MATRIXF* p3,MATRIXF* p4,MATRIXF* p5,size_t nv,char nodiag)
{
#define	CLEANUP			CLEANMATF(tnew)CLEANMATF(tnew2)\
						CLEANVECF(llr1)CLEANMATF(llr2)CLEANMATF(llr3)\
						CLEANMATF(llr4)CLEANMATF(llr5)
	MATRIXF		*tnew,*tnew2;	//(nt,ns) Supernormalized transcript matrix
	VECTORF		*llr1;
	MATRIXF		*llr2,*llr3,*llr4,*llr5;
	VECTORFF(view)	vv;
	int			ret;
#ifndef NDEBUG
	size_t		ng,nt,ns;
	
	ng=g->size1;
	nt=t2->size1;
	ns=g->size2;
#endif

	tnew=tnew2=llr2=llr3=llr4=llr5=0;
	llr1=0;
	
	//Validation
	assert(!((t->size1!=ng)||(t->size2!=ns)||(t2->size2!=ns)
		||(p1&&(p1->size!=ng))
		||(p2&&((p2->size1!=ng)||(p2->size2!=nt)))
		||(p3&&((p3->size1!=ng)||(p3->size2!=nt)))
		||(p4&&((p4->size1!=ng)||(p4->size2!=nt)))
		||(p5&&((p5->size1!=ng)||(p5->size2!=nt)))));
	assert(!(nv>CONST_NV_MAX));
	
	llr1=VECTORFF(alloc)(g->size1);
	tnew=MATRIXFF(alloc)(t->size1,t->size2);
	tnew2=MATRIXFF(alloc)(t2->size1,t2->size2);
	llr2=MATRIXFF(alloc)(t->size1,t2->size1);
	llr3=MATRIXFF(alloc)(t->size1,t2->size1);
	llr4=MATRIXFF(alloc)(t->size1,t2->size1);
	llr5=MATRIXFF(alloc)(t->size1,t2->size1);
	if(!(llr1&&tnew&&tnew2&&llr2&&llr3&&llr4&&llr5))
		ERRRET("Not enough memory.")

	//Step 1: Supernormalization
	LOG(9,"Supernormalizing...")
	MATRIXFF(memcpy)(tnew,t);
	ret=supernormalizea_byrow(tnew);
	MATRIXFF(memcpy)(tnew2,t2);
	ret=ret||supernormalizea_byrow(tnew2);
	if(ret)
		ERRRET("Supernormalization failed.")
	
	//Step 2: Log likelihood ratios from nonpermuted data
	LOG(9,"Calculating real log likelihood ratios...")
	if(pij_gassist_llr(g,tnew,tnew2,llr1,llr2,llr3,llr4,llr5,nv))
		ERRRET("pij_gassist_llr failed.")
	if(nodiag)
	{
		vv=MATRIXFF(diagonal)(llr2);
		VECTORFF(set_zero)(&vv.vector);
		vv=MATRIXFF(diagonal)(llr3);
		VECTORFF(set_zero)(&vv.vector);
		vv=MATRIXFF(diagonal)(llr4);
		VECTORFF(set_zero)(&vv.vector);
		vv=MATRIXFF(diagonal)(llr5);
		VECTORFF(set_zero)(&vv.vector);
	}
	//Step 3: Convert log likelihood ratios to probabilities with permuted data
	if(pij_gassist_llrtopijs_tot(g,tnew,tnew2,llr1,llr2,llr3,llr4,llr5,p1,p2,p3,p4,p5,nv,nodiag))
		LOG(4,"Failed to convert all log likelihood ratios to probabilities.")

	//Cleanup
	CLEANUP
	return ret;
#undef	CLEANUP		
}

int pijs_gassist_a(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,VECTORF* p1,MATRIXF* p2,MATRIXF* p3,MATRIXF* p4,MATRIXF* p5,size_t nv,char nodiag)
{
#define	CLEANUP			CLEANMATF(tnew)CLEANMATF(tnew2)\
						CLEANVECF(llr1)CLEANMATF(llr2)CLEANMATF(llr3)\
						CLEANMATF(llr4)CLEANMATF(llr5)
	MATRIXF		*tnew,*tnew2;	//(nt,ns) Supernormalized transcript matrix
	VECTORF		*llr1;
	MATRIXF		*llr2,*llr3,*llr4,*llr5;
	VECTORFF(view)	vv;
	int			ret;
#ifndef NDEBUG
	size_t		ng,nt,ns;
	
	ng=g->size1;
	nt=t2->size1;
	ns=g->size2;
#endif

	tnew=tnew2=llr2=llr3=llr4=llr5=0;
	llr1=0;
	
	//Validation
	assert(!((t->size1!=ng)||(t->size2!=ns)||(t2->size2!=ns)
		||(p1&&(p1->size!=ng))
		||(p2&&((p2->size1!=ng)||(p2->size2!=nt)))
		||(p3&&((p2->size1!=ng)||(p2->size2!=nt)))
		||(p4&&((p2->size1!=ng)||(p2->size2!=nt)))
		||(p5&&((p3->size1!=ng)||(p3->size2!=nt)))));
	assert(!(nv>CONST_NV_MAX));
	
	llr1=VECTORFF(alloc)(g->size1);
	tnew=MATRIXFF(alloc)(t->size1,t->size2);
	tnew2=MATRIXFF(alloc)(t2->size1,t2->size2);
	llr2=MATRIXFF(alloc)(t->size1,t2->size1);
	llr3=MATRIXFF(alloc)(t->size1,t2->size1);
	llr4=MATRIXFF(alloc)(t->size1,t2->size1);
	llr5=MATRIXFF(alloc)(t->size1,t2->size1);
	if(!(llr1&&tnew&&tnew2&&llr2&&llr3&&llr4&&llr5))
		ERRRET("Not enough memory.")

	//Step 1: Supernormalization
	LOG(9,"Supernormalizing...")
	MATRIXFF(memcpy)(tnew,t);
	ret=supernormalizea_byrow(tnew);
	MATRIXFF(memcpy)(tnew2,t2);
	ret=ret||supernormalizea_byrow(tnew2);
	if(ret)
		ERRRET("Supernormalization failed.")
	
	//Step 2: Log likelihood ratios from nonpermuted data
	LOG(9,"Calculating real log likelihood ratios...")
	if(pij_gassist_llr(g,tnew,tnew2,llr1,llr2,llr3,llr4,llr5,nv))
		ERRRET("pij_gassist_llr failed.")
	if(nodiag)
	{
		vv=MATRIXFF(diagonal)(llr2);
		VECTORFF(set_zero)(&vv.vector);
		vv=MATRIXFF(diagonal)(llr3);
		VECTORFF(set_zero)(&vv.vector);
		vv=MATRIXFF(diagonal)(llr4);
		VECTORFF(set_zero)(&vv.vector);
		vv=MATRIXFF(diagonal)(llr5);
		VECTORFF(set_zero)(&vv.vector);
	}
	//Step 3: Convert log likelihood ratios to probabilities with permuted data
	if(pij_gassist_llrtopijs_a(g,llr1,llr2,llr3,llr4,llr5,p1,p2,p3,p4,p5,nv,nodiag))
		LOG(4,"Failed to convert all log likelihood ratios to probabilities.")

	//Cleanup
	CLEANUP
	return ret;
#undef	CLEANUP		
}

/* Estimates the probability p(E->A->B) from genotype and expression data. Combines results
 * from any pij_gassist_pijs. For more information, see pij_gassist_pijs_ab.
 * ans:	(ng,nt) Output matrix for probabilities. ans[A,B] is p(E->A->B).
 * pijs:	Function to calculate pijs.
 */
static int pij_gassist_any(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,MATRIXF* ans,size_t nv,char nodiag,int (*pijs)(const MATRIXG*,const MATRIXF*,const MATRIXF*,VECTORF*,MATRIXF*,MATRIXF*,MATRIXF*,MATRIXF*,size_t,char))
{
#define	CLEANUP			CLEANVECF(p1)CLEANMATF(p2)CLEANMATF(p3)CLEANMATF(p4)
	VECTORF	*p1;
	MATRIXF	*p2,*p3,*p4;
	size_t	ng=g->size1;
	size_t	nt=t2->size1;

	assert(g&&t&&t2&&ans&&pijs);
	assert((g->size2==t->size2)&&(g->size2==t2->size2));
	assert((t->size1==ng)&&(ans->size1==ng)&&(ans->size2==nt)&&(nv>1));
	p1=VECTORFF(alloc)(ng);
	p2=MATRIXFF(alloc)(ng,nt);
	p3=MATRIXFF(alloc)(ng,nt);
	p4=MATRIXFF(alloc)(ng,nt);
	if(!(p1&&p2&&p3&&p4))
		ERRRET("Not enough memory.")
	if(pijs(g,t,t2,p1,p2,p3,p4,ans,nv,nodiag))
		ERRRET("pij_gassist_pijs failed.")
		
	//Combine tests
	{
		FTYPE	s;
		size_t	i;
		VECTORFF(view)	vv;
		MATRIXFF(mul_elements)(ans,p2);
		s=MATRIXFF(max)(ans)+MATRIXFF(max)(p4);
		MATRIXFF(add)(ans,p4);
		MATRIXFF(scale)(ans,1./s);
		for(i=0;i<g->size1;i++)
		{
			vv=MATRIXFF(row)(ans,i);
			VECTORFF(scale)(&vv.vector,VECTORFF(get)(p1,i));
		}
	}
	
	//Cleanup
	CLEANUP
	return 0;
#undef	CLEANUP
}

int pij_gassist_a(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,MATRIXF* ans,size_t nv,char nodiag)
{
	return pij_gassist_any(g,t,t2,ans,nv,nodiag,pijs_gassist_a);
}

int pij_gassist_tot(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,MATRIXF* ans,size_t nv,char nodiag)
{
	return pij_gassist_any(g,t,t2,ans,nv,nodiag,pijs_gassist_tot);
}


















