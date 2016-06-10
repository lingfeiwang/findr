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
#include <math.h>
#include <float.h>
#include <string.h>
#include "../../base/gsl/math.h"
#include "../../base/gsl/histogram.h"
#include "../../base/gsl/blas.h"
#include "../../base/logger.h"
#include "../../base/macros.h"
#include "../../base/data_process.h"
#include "../../base/histogram.h"
#include "../nullmodeler.h"
#include "../nullsampler.h"
#include "../nullhist.h"
#include "../llrtopij.h"
#include "llrtopij.h"

#pragma GCC diagnostic ignored "-Wunused-parameter"


int pij_gassist_llrtopij1_1(VECTORF* p1)
{
	LOG(9,"Converting LLR to probabilities for step 1. Filling with 1.")
	VECTORFF(set_all)(p1,1);
	return 0;
}

int pij_gassist_llrtopij_getnvs(const MATRIXG* g,VECTORG* ans,size_t nv)
{
#define	CLEANUP	AUTOFREEVEC(cond)
	size_t		i,j,ng,ns;
	GTYPE	*gptr;
	
	assert(g->size1==ans->size);
	ng=g->size1;
	ns=g->size2;
	AUTOALLOCVECUC(cond,nv,200)
	if(!cond)
		ERRRET("Not enough memory.")
	VECTORGF(set_zero)(ans);
	for(i=0;i<ng;i++)
	{
		VECTORUCF(set_zero)(cond);
		for(j=0;j<ns;j++)
			VECTORUCF(set)(cond,MATRIXGF(get)(g,i,j),1);
		gptr=VECTORGF(ptr)(ans,i);
		for(j=0;j<nv;j++)
			gptr[0]=(GTYPE)(gptr[0]+VECTORUCF(get)(cond,j));
	}
	CLEANUP
	return 0;
#undef	CLEANUP
}

int pij_gassist_llrtopijs_nv(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,const VECTORF* llr1,const MATRIXF* llr2,const MATRIXF* llr3,VECTORF* p1,MATRIXF* p2,MATRIXF* p3,size_t nv,size_t nthread,char nodiag)
{
#define	CLEANUP	CLEANVECF(vp2)CLEANVECF(vp3)CLEANVECF(tp)CLEANVECF(tp2)\
				CLEANVECG(vnv)CLEANVECUC(vselect)CLEANMATG(gs)

	int		ret;
	size_t	ng,nt,ns,n,nc,i,ngnow;
	double	nratio;		//Ratio of null distribution
	VECTORF *vp2,*vp3,*tp,*tp2;
	VECTORG	*vnv;
	VECTORUC	*vselect;
	MATRIXG	*gs;
	MATRIXGF(view)	mgs;
	VECTORFF(view)	vv1,vv2;
	MATRIXFF(view)	mv;

	//Initialize
	vp2=vp3=tp=tp2=0;
	ng=t->size1;
	nt=t2->size1;
	ns=t->size2;
	n=ng*nt;
	if(nodiag)
		n-=GSL_MIN(ng,nt);
	vp2=VECTORFF(alloc)(n);
	vp3=VECTORFF(alloc)(n);
	tp=VECTORFF(alloc)(n);
	tp2=VECTORFF(alloc)(n);
	vnv=VECTORGF(alloc)(ng);
	vselect=VECTORUCF(alloc)(ng);
	gs=MATRIXGF(alloc)(ng,ns);
	
	if(!(vp2&&vp3&&tp&&tp2&&vnv&&vselect))
		ERRRET("Not enough memory.")
	
	if(VECTORFF(first_nan)(llr1)>=0)
		LOG(4,"Infinity found for log likelihood ratio, possibly because data contains fully correlated/anticorrelated columns. This may affect downstream analyses.")
	ret=pij_gassist_llrtopij1_1(p1);
	if(ret)
		ERRRET("Failed to calculate probabilities in step 1.")
	
	
	ret=pij_gassist_llrtopij_getnvs(g,vnv,nv);
	if(ret)
		ERRRET("Failed to get value counts for each genotype.")
	
	MATRIXFF(set_all)(p2,0);
	MATRIXFF(set_all)(p3,0);
	for(i=0;i<ng;i++)
		if(VECTORGF(get)(vnv,i)<=1)
			ERRRET("Genotype has only one value at location "PRINTFSIZET".",i);
	for(i=2;i<=nv;i++)
	{
		VECTORGF(eq)(vnv,vselect,(GTYPE)i);
		ngnow=MATRIXGF(rows_save)(g,gs,vselect);
		if(!ngnow)
			continue;
		mgs=MATRIXGF(submatrix)(gs,0,0,ngnow,ns);
		if(nodiag)
			nc=MATRIXFF(rows_save_nodiag)(llr2,tp,vselect);
		else
		{
			mv=MATRIXFF(view_vector)(tp,ng,nt);
			MATRIXFF(rows_save)(llr2,&mv.matrix,vselect);
			nc=ngnow*nt;
		}
		vv1=VECTORFF(subvector)(tp,0,nc);
		vv2=VECTORFF(subvector)(vp2,0,nc);
		if(VECTORFF(first_nan)(&vv1.vector)>=0)
		{
			LOG(4,"Infinity found for log likelihood ratio of nv="PRINTFSIZET", possibly because data contains fully correlated/anticorrelated columns. This may affect downstream analyses.",i)
			VECTORFF(set_nan)(&vv1.vector,0);
		}
		ret=pij_gassist_llrtopij2c_tot(&mgs.matrix,&vv1.vector,&vv2.vector,nv,&nratio);
		if(ret)
			ERRRET("Failed to calculate probabilities in step 2 for nv="PRINTFSIZET".",i)
		if(nodiag)
			MATRIXFF(rows_load_nodiag)(vp2,p2,vselect);
		else
		{
			mv=MATRIXFF(view_vector)(vp2,ng,nt);
			MATRIXFF(rows_load)(&mv.matrix,p2,vselect);
		}

		if(nodiag)
			MATRIXFF(rows_save_nodiag)(llr3,tp2,vselect);
		else
		{
			mv=MATRIXFF(view_vector)(tp2,ng,nt);
			MATRIXFF(rows_save)(llr3,&mv.matrix,vselect);
		}
		vv1=VECTORFF(subvector)(tp2,0,nc);
		vv2=VECTORFF(subvector)(vp3,0,nc);
		if(VECTORFF(first_nan)(&vv1.vector)>=0)
		{
			LOG(4,"Infinity found for log likelihood ratio of nv="PRINTFSIZET", possibly because data contains fully correlated/anticorrelated columns. This may affect downstream analyses.",i)
			VECTORFF(set_nan)(&vv1.vector,0);
		}
		ret=pij_gassist_llrtopij3_tot(&mgs.matrix,&vv1.vector,&vv2.vector,nv,&nratio);
		if(ret)
			ERRRET("Failed to calculate probabilities in step 3 for nv="PRINTFSIZET".",i)
		if(nodiag)
			MATRIXFF(rows_load_nodiag)(vp3,p3,vselect);
		else
		{
			mv=MATRIXFF(view_vector)(vp3,ng,nt);
			MATRIXFF(rows_load)(&mv.matrix,p3,vselect);
		}
	}
	
	
	//Free memory
	CLEANUP
	return 0;
#undef	CLEANUP
}















