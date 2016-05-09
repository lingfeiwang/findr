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


int pijs_gassist_tot(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,VECTORF* p1,MATRIXF* p2b,MATRIXF* p2c,MATRIXF* p3,size_t nv,char nodiag)
{
#define	CLEANUP			CLEANMATF(tnew)CLEANMATF(tnew2)\
						CLEANVECF(llr1)CLEANMATF(llr2b)CLEANMATF(llr2c)\
						CLEANMATF(llr3)
	MATRIXF		*tnew,*tnew2;			//(nt,ns) Supernormalized transcript matrix
	VECTORF		*llr1;
	MATRIXF		*llr2b,*llr2c,*llr3;
	VECTORFF(view)	vv;
	int			ret;
#ifndef NDEBUG
	size_t		ng,nt,ns;
	
	ng=g->size1;
	nt=t2->size1;
	ns=g->size2;
#endif

	tnew=tnew2=llr2b=llr2c=llr3=0;
	llr1=0;
	
	//Validation
	assert(!((t->size1!=ng)||(t->size2!=ns)||(t2->size2!=ns)
		||(p1&&(p1->size!=ng))
		||(p2b&&((p2b->size1!=ng)||(p2b->size2!=nt)))
		||(p2c&&((p2c->size1!=ng)||(p2c->size2!=nt)))
		||(p3&&((p3->size1!=ng)||(p3->size2!=nt)))));
	assert(!(nv>CONST_NV_MAX));
	
	llr1=VECTORFF(alloc)(g->size1);
	tnew=MATRIXFF(alloc)(t->size1,t->size2);
	tnew2=MATRIXFF(alloc)(t2->size1,t2->size2);
	llr2b=MATRIXFF(alloc)(t->size1,t2->size1);
	llr2c=MATRIXFF(alloc)(t->size1,t2->size1);
	llr3=MATRIXFF(alloc)(t->size1,t2->size1);
	if(!(llr1&&tnew&&tnew2&&llr2b&&llr2c&&llr3))
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
	if(pij_gassist_llr(g,tnew,tnew2,llr1,llr2b,llr2c,llr3,nv))
		ERRRET("pij_gassist_llr failed.")
	if(nodiag)
	{
		vv=MATRIXFF(diagonal)(llr2b);
		VECTORFF(set_zero)(&vv.vector);
		vv=MATRIXFF(diagonal)(llr2c);
		VECTORFF(set_zero)(&vv.vector);
		vv=MATRIXFF(diagonal)(llr3);
		VECTORFF(set_zero)(&vv.vector);
	}
	//Step 3: Convert log likelihood ratios to probabilities with permuted data
	if(pij_gassist_llrtopijs_tot(g,tnew,tnew2,llr1,llr2b,llr2c,llr3,p1,p2b,p2c,p3,nv,nodiag))
		LOG(4,"Failed to convert all log likelihood ratios to probabilities.")

	//Cleanup
	CLEANUP
	return ret;
#undef	CLEANUP		
}

int pijs_gassist_a(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,VECTORF* p1,MATRIXF* p2b,MATRIXF* p2c,MATRIXF* p3,size_t nv,char nodiag)
{
#define	CLEANUP			CLEANMATF(tnew)CLEANMATF(tnew2)\
						CLEANVECF(llr1)CLEANMATF(llr2b)CLEANMATF(llr2c)\
						CLEANMATF(llr3)
	MATRIXF		*tnew,*tnew2;			//(nt,ns) Supernormalized transcript matrix
	VECTORF		*llr1;
	MATRIXF		*llr2b,*llr2c,*llr3;
	VECTORFF(view)	vv;
	int			ret;
#ifndef NDEBUG
	size_t		ng,nt,ns;
	
	ng=g->size1;
	nt=t2->size1;
	ns=g->size2;
#endif

	tnew=tnew2=llr2b=llr2c=llr3=0;
	llr1=0;
	
	//Validation
	assert(!((t->size1!=ng)||(t->size2!=ns)||(t2->size2!=ns)
		||(p1&&(p1->size!=ng))
		||(p2b&&((p2b->size1!=ng)||(p2b->size2!=nt)))
		||(p2c&&((p2c->size1!=ng)||(p2c->size2!=nt)))
		||(p3&&((p3->size1!=ng)||(p3->size2!=nt)))));
	assert(!(nv>CONST_NV_MAX));
	
	llr1=VECTORFF(alloc)(g->size1);
	tnew=MATRIXFF(alloc)(t->size1,t->size2);
	tnew2=MATRIXFF(alloc)(t2->size1,t2->size2);
	llr2b=MATRIXFF(alloc)(t->size1,t2->size1);
	llr2c=MATRIXFF(alloc)(t->size1,t2->size1);
	llr3=MATRIXFF(alloc)(t->size1,t2->size1);
	if(!(llr1&&tnew&&tnew2&&llr2b&&llr2c&&llr3))
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
	if(pij_gassist_llr(g,tnew,tnew2,llr1,llr2b,llr2c,llr3,nv))
		ERRRET("pij_gassist_llr failed.")
	if(nodiag)
	{
		vv=MATRIXFF(diagonal)(llr2b);
		VECTORFF(set_zero)(&vv.vector);
		vv=MATRIXFF(diagonal)(llr2c);
		VECTORFF(set_zero)(&vv.vector);
		vv=MATRIXFF(diagonal)(llr3);
		VECTORFF(set_zero)(&vv.vector);
	}
	//Step 3: Convert log likelihood ratios to probabilities with permuted data
	if(pij_gassist_llrtopijs_a(g,llr1,llr2b,llr2c,llr3,p1,p2b,p2c,p3,nv,nodiag))
		LOG(4,"Failed to convert all log likelihood ratios to probabilities.")

	//Cleanup
	CLEANUP
	return ret;
#undef	CLEANUP		
}


int pij_gassist_any(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,MATRIXF* ansb,MATRIXF* ansc,size_t nv,char nodiag,int (*pijs)(const MATRIXG*,const MATRIXF*,const MATRIXF*,VECTORF*,MATRIXF*,MATRIXF*,MATRIXF*,size_t,char))
{
#define	CLEANUP			CLEANVECF(p1)CLEANMATF(p3)
	VECTORF	*p1;
	MATRIXF	*p3;
	VECTORFF(view)	vv;
	size_t	i;
	
	p1=VECTORFF(alloc)(g->size1);
	p3=MATRIXFF(alloc)(g->size1,t2->size1);
	if(!(p1&&p3))
		ERRRET("Not enough memory.")
	if(pijs(g,t,t2,p1,ansb,ansc,p3,nv,nodiag))
		ERRRET("pij_gassist_pijs failed.")

	for(i=0;i<g->size1;i++)
	{
		vv=MATRIXFF(row)(ansb,i);
		VECTORFF(scale)(&vv.vector,VECTORFF(get)(p1,i));
		vv=MATRIXFF(row)(ansc,i);
		VECTORFF(scale)(&vv.vector,VECTORFF(get)(p1,i));
	}
	MATRIXFF(mul_elements)(ansb,p3);
	MATRIXFF(mul_elements)(ansc,p3);
	
	//Cleanup
	CLEANUP
	return 0;
#undef	CLEANUP
}

int pij_gassist_a(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,MATRIXF* ansb,MATRIXF* ansc,size_t nv,char nodiag)
{
	return pij_gassist_any(g,t,t2,ansb,ansc,nv,nodiag,pijs_gassist_a);
}

int pij_gassist_tot(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,MATRIXF* ansb,MATRIXF* ansc,size_t nv,char nodiag)
{
	return pij_gassist_any(g,t,t2,ansb,ansc,nv,nodiag,pijs_gassist_tot);
}


















