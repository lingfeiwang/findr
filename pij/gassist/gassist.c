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
#include "../../base/config.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "../../base/logger.h"
#include "../../base/macros.h"
#include "../../base/const.h"
#include "../../base/supernormalize.h"
#include "../../base/threading.h"
#include "../../base/data_process.h"
#include "../llrtopij.h"
#include "llr.h"
#include "llrtopv.h"
#include "llrtopij.h"
#include "nullhist.h"
#include "gassist.h"

int pijs_gassist_pv(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,VECTORF* p1,MATRIXF* p2,MATRIXF* p3,MATRIXF* p4,MATRIXF* p5,size_t nv,size_t memlimit)
{
#define	CLEANUP			CLEANMATF(tnew)CLEANMATF(tnew2)
	MATRIXF		*tnew,*tnew2;	//(nt,ns) Supernormalized transcript matrix
	MATRIXFF(view)	mvt,mvp2,mvp3,mvp4,mvp5;
	VECTORFF(view)	vvp1;
	int			ret;
	size_t		i,ng,ns,ngnow,nsplit;
#ifndef NDEBUG
	size_t		nt;
	
	nt=t2->size1;
#endif
	ns=g->size2;
	ng=g->size1;

	tnew=tnew2=0;

	//Validation
	assert(!((t->size1!=ng)||(t->size2!=ns)||(t2->size2!=ns)
		||(p1&&(p1->size!=ng))
		||(p2&&((p2->size1!=ng)||(p2->size2!=nt)))
		||(p3&&((p2->size1!=ng)||(p2->size2!=nt)))
		||(p4&&((p2->size1!=ng)||(p2->size2!=nt)))
		||(p5&&((p3->size1!=ng)||(p3->size2!=nt)))));
	assert(!(nv>CONST_NV_MAX));
	assert(memlimit);
	if(ns<4)
		ERRRET("Needs at least 4 samples to compute p-values.")

	{
		size_t mem1,mem2;
		mem1=g->size1*g->size2*sizeof(GTYPE)+(2*t->size1*t->size2+2*t2->size1*t2->size2+p1->size+p2->size1*p2->size2*4)*sizeof(FTYPE);
		mem2=t2->size1*2*nv*sizeof(FTYPE);
		if((memlimit<=mem1)||!(nsplit=(memlimit-mem1)/mem2))
			ERRRET("Memory limit lower than minimum memory needed. Try increasing your memory usage limit.")
		LOG(10,"Memory limit: %lu bytes.",memlimit)
		nsplit=(size_t)ceil((float)ng/ceil((float)ng/(float)nsplit));
		//if(nsplit<ng)
			LOG(9,"Splitting %lu primary targets into groups of about size %lu.",ng,nsplit)
	}
	
	tnew=MATRIXFF(alloc)(t->size1,t->size2);
	tnew2=MATRIXFF(alloc)(t2->size1,t2->size2);
	if(!(tnew&&tnew2))
		ERRRET("Not enough memory.")

	//Step 1: Supernormalization
	LOG(9,"Supernormalizing...")
	MATRIXFF(memcpy)(tnew,t);
	ret=supernormalizea_byrow(tnew);
	MATRIXFF(memcpy)(tnew2,t2);
	ret=ret||supernormalizea_byrow(tnew2);
	if(ret)
		ERRRET("Supernormalization failed.")
	
	for(i=0;i<ng;i+=nsplit)
	{
		ngnow=GSL_MIN(ng-i,nsplit);

		MATRIXGF(const_view) mvg=MATRIXGF(const_submatrix)(g,i,0,ngnow,g->size2);
		mvt=MATRIXFF(submatrix)(tnew,i,0,ngnow,tnew->size2);
		vvp1=VECTORFF(subvector)(p1,i,ngnow);
		mvp2=MATRIXFF(submatrix)(p2,i,0,ngnow,p2->size2);
		mvp3=MATRIXFF(submatrix)(p3,i,0,ngnow,p3->size2);
		mvp4=MATRIXFF(submatrix)(p4,i,0,ngnow,p4->size2);
		mvp5=MATRIXFF(submatrix)(p5,i,0,ngnow,p5->size2);
		//Step 2: Log likelihood ratios from nonpermuted data
		LOG(9,"Calculating real log likelihood ratios...")
		if(pij_gassist_llr(&mvg.matrix,&mvt.matrix,tnew2,&vvp1.vector,&mvp2.matrix,&mvp3.matrix,&mvp4.matrix,&mvp5.matrix,nv))
			ERRRET("pij_gassist_llr failed.")
		//Step 3: Convert log likelihood ratios to p-values
		LOG(9,"Converting log likelihood ratios into p-values...")
		if((ret=pij_gassist_llrtopvs(&vvp1.vector,&mvp2.matrix,&mvp3.matrix,&mvp4.matrix,&mvp5.matrix,&mvg.matrix,nv)))
			LOG(4,"Failed to convert all log likelihood ratios to p-values.")
	}

	//Cleanup
	CLEANUP
	return ret;
#undef	CLEANUP		
}

int pijs_gassist(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,VECTORF* p1,MATRIXF* p2,MATRIXF* p3,MATRIXF* p4,MATRIXF* p5,size_t nv,char nodiag,size_t memlimit)
{
#define	CLEANUP			CLEANMATF(tnew)CLEANMATF(tnew2)for(i=0;i<4;i++){if(hnull[i])for(j=0;j<nv-1;j++)CLEANHIST(hnull[i][j]);CLEANMEM(hnull[i]);}
	MATRIXF			*tnew,*tnew2;	//(nt,ns) Supernormalized transcript matrix
	MATRIXFF(view)	mvt,mvp2,mvp3,mvp4,mvp5;
	VECTORFF(view)	vv,vvp1;
	gsl_histogram**	hnull[4]={0,0,0,0};
	int				ret;
	size_t			i,j,ng,nt,ngnow,nsplit,ns;
	
	nt=t2->size1;
	ng=g->size1;
	ns=g->size2;

	tnew=tnew2=0;

	//Validation
	assert(!((t->size1!=ng)||(t->size2!=ns)||(t2->size2!=ns)
		||(p1&&(p1->size!=ng))
		||(p2&&((p2->size1!=ng)||(p2->size2!=nt)))
		||(p3&&((p2->size1!=ng)||(p2->size2!=nt)))
		||(p4&&((p2->size1!=ng)||(p2->size2!=nt)))
		||(p5&&((p3->size1!=ng)||(p3->size2!=nt)))));
	assert(!(nv>CONST_NV_MAX));
	assert(memlimit);
	if(ns<4)
		ERRRET("Needs at least 4 samples to compute probabilities.")
	//Defaults to 8GB memory usage
	{
		size_t mem1,mem2;
		mem1=g->size1*g->size2*sizeof(GTYPE)+(2*t->size1*t->size2+2*t2->size1*t2->size2+p1->size+p2->size1*p2->size2*4)*sizeof(FTYPE);
		mem2=t2->size1*2*nv*sizeof(FTYPE);
		if((memlimit<=mem1)||!(nsplit=(memlimit-mem1)/mem2))
			ERRRET("Memory limit lower than minimum memory needed. Try increasing your memory usage limit.")
		LOG(10,"Memory limit: %lu bytes.",memlimit)
		nsplit=(size_t)ceil((float)ng/ceil((float)ng/(float)nsplit));
		if(nsplit<ng)
			LOG(9,"Splitting %lu primary targets into groups of about size %lu.",ng,nsplit)
	}
	
	tnew=MATRIXFF(alloc)(t->size1,t->size2);
	tnew2=MATRIXFF(alloc)(t2->size1,t2->size2);
	if(!(tnew&&tnew2))
		ERRRET("Not enough memory.")

	//Check for identical rows in input data
	{
		VECTORFF(view) vbuff1=MATRIXFF(column)(tnew,0);
		VECTORFF(view) vbuff2=MATRIXFF(row)(tnew2,0);
		MATRIXFF(cmprow)(t,t2,&vbuff1.vector,&vbuff2.vector,nodiag,1);
	}

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
	for(i=0;i<ng;i+=nsplit)
	{
		ngnow=GSL_MIN(ng-i,nsplit);

		MATRIXGF(const_view) mvg=MATRIXGF(const_submatrix)(g,i,0,ngnow,g->size2);
		mvt=MATRIXFF(submatrix)(tnew,i,0,ngnow,tnew->size2);
		vvp1=VECTORFF(subvector)(p1,i,ngnow);
		mvp2=MATRIXFF(submatrix)(p2,i,0,ngnow,p2->size2);
		mvp3=MATRIXFF(submatrix)(p3,i,0,ngnow,p3->size2);
		mvp4=MATRIXFF(submatrix)(p4,i,0,ngnow,p4->size2);
		mvp5=MATRIXFF(submatrix)(p5,i,0,ngnow,p5->size2);
		if(pij_gassist_llr(&mvg.matrix,&mvt.matrix,tnew2,&vvp1.vector,&mvp2.matrix,&mvp3.matrix,&mvp4.matrix,&mvp5.matrix,nv))
			ERRRET("pij_gassist_llr failed.")
	}
	
	//Step 3: Obtain null histograms
	{
		FTYPE			dmax[4];
		dmax[0]=pij_llrtopij_llrmatmax(p2,nodiag);
		dmax[1]=pij_llrtopij_llrmatmax(p3,nodiag);
		dmax[2]=pij_llrtopij_llrmatmax(p4,nodiag);
		dmax[3]=pij_llrtopij_llrmatmax(p5,nodiag);
		if(!(dmax[0]&&dmax[1]&&dmax[2]&&dmax[3]))
			ERRRET("Negative or NAN found in LLR.")
		if(pij_gassist_nullhists(hnull,nt,ns,nv,dmax))
			ERRRET("Failed to construct null histograms.")
	}

	//Step 4: Convert log likelihood ratios to probabilities
	for(i=0;i<ng;i+=nsplit)
	{
		ngnow=GSL_MIN(ng-i,nsplit);

		MATRIXGF(const_view) mvg=MATRIXGF(const_submatrix)(g,i,0,ngnow,g->size2);
		mvt=MATRIXFF(submatrix)(tnew,i,0,ngnow,tnew->size2);
		vvp1=VECTORFF(subvector)(p1,i,ngnow);
		mvp2=MATRIXFF(submatrix)(p2,i,0,ngnow,p2->size2);
		mvp3=MATRIXFF(submatrix)(p3,i,0,ngnow,p3->size2);
		mvp4=MATRIXFF(submatrix)(p4,i,0,ngnow,p4->size2);
		mvp5=MATRIXFF(submatrix)(p5,i,0,ngnow,p5->size2);
		if((ret=pij_gassist_llrtopijs(&mvg.matrix,&vvp1.vector,&mvp2.matrix,&mvp3.matrix,&mvp4.matrix,&mvp5.matrix,nv,(const gsl_histogram* const **)hnull,nodiag,(long)i)))
			LOG(4,"Failed to convert all log likelihood ratios to probabilities.")
		if(nodiag)
		{
			vv=MATRIXFF(superdiagonal)(&mvp2.matrix,i);
			VECTORFF(set_zero)(&vv.vector);
			vv=MATRIXFF(superdiagonal)(&mvp3.matrix,i);
			VECTORFF(set_zero)(&vv.vector);
			vv=MATRIXFF(superdiagonal)(&mvp4.matrix,i);
			VECTORFF(set_zero)(&vv.vector);
			vv=MATRIXFF(superdiagonal)(&mvp5.matrix,i);
			VECTORFF(set_zero)(&vv.vector);
		}
	}

	//Cleanup
	CLEANUP
	return ret;
#undef	CLEANUP		
}

int pij_gassist(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,MATRIXF* ans,size_t nv,char nodiag,size_t memlimit)
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
	if(pijs_gassist(g,t,t2,p1,p2,p3,p4,ans,nv,nodiag,memlimit))
		ERRRET("pij_gassist_pijs failed.")
		
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

int pij_gassist_trad(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,MATRIXF* ans,size_t nv,char nodiag,size_t memlimit)
{
#define	CLEANUP			CLEANVECF(p1)CLEANMATF(p2)CLEANMATF(p4)CLEANMATF(p5)
	VECTORF	*p1;
	MATRIXF	*p2,*p4,*p5;
	size_t	ng=g->size1;
	size_t	nt=t2->size1;

	assert(g&&t&&t2&&ans);
	assert((g->size2==t->size2)&&(g->size2==t2->size2));
	assert((t->size1==ng)&&(ans->size1==ng)&&(ans->size2==nt)&&(nv>1));
	p1=VECTORFF(alloc)(ng);
	p2=MATRIXFF(alloc)(ng,nt);
	p4=MATRIXFF(alloc)(ng,nt);
	p5=MATRIXFF(alloc)(ng,nt);
	if(!(p1&&p2&&p5&&p4))
		ERRRET("Not enough memory.")
	if(pijs_gassist(g,t,t2,p1,p2,ans,p4,p5,nv,nodiag,memlimit))
		ERRRET("pij_gassist_pijs failed.")
		
	//Combine tests
	MATRIXFF(mul_elements)(ans,p2);
	
	//Cleanup
	CLEANUP
	return 0;
#undef	CLEANUP
}















