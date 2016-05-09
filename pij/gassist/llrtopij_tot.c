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
#include "../../base/threading.h"
#include "../../base/histogram.h"
#include "../nullhist.h"
#include "nullhist.h"
#include "llrtopij_tot.h"
#include "llrtopij.h"
#pragma GCC diagnostic ignored "-Wunused-parameter"


int pij_gassist_llrtopij1_tot(const MATRIXG* g,const VECTORF* llr,VECTORF* p,size_t nv,double* nratio)
{
	const struct histogram_unequalbins_exp_param			ph={g->size2-nv,1};
	const struct pij_gassist_nullhist_analytical_pdf_param	pnh={g,nv,5};
	return pij_llrtopij_tot_convert(llr,llr,p,histogram_unequalbins_exp,pij_gassist_nullhist_analytical1_pdf,&ph,&pnh,nratio);
}

int pij_gassist_llrtopij2b_tot(const MATRIXG* g,const VECTORF* llr,VECTORF* p,size_t nv,double* nratio)
{
	const struct histogram_unequalbins_exp_param			ph={g->size2-nv,1};
	const struct pij_gassist_nullhist_analytical_pdf_param	pnh={g,nv,5};
	return pij_llrtopij_tot_convert(llr,llr,p,histogram_unequalbins_exp,pij_gassist_nullhist_analytical2b_pdf,&ph,&pnh,nratio);
}

int pij_gassist_llrtopij2c_tot(const MATRIXG* g,const VECTORF* llr,VECTORF* p,size_t nv,double* nratio)
{
	const struct histogram_unequalbins_exp_param			ph={g->size2-nv,1};
	const struct pij_gassist_nullhist_analytical_pdf_param	pnh={g,nv,5};
	return pij_llrtopij_tot_convert(llr,llr,p,histogram_unequalbins_exp,pij_gassist_nullhist_analytical2c_pdf,&ph,&pnh,nratio);
}

int pij_gassist_llrtopij3_tot(const MATRIXG* g,const VECTORF* llr,VECTORF* p,size_t nv,double* nratio)
{
#define	CLEANUP
	const struct histogram_unequalbins_exp_param			ph={g->size2-nv,1};
	const struct pij_gassist_nullhist_analytical_pdf_param	pnh={g,nv,5};
	if(pij_llrtopij_tot_convert(llr,llr,p,histogram_unequalbins_exp,pij_gassist_nullhist_analytical3_pdf,&ph,&pnh,nratio))
		ERRRET("Conversion failed.")
	VECTORFF(scale)(p,-1);
	VECTORFF(add_constant)(p,1);
	return 0;	
#undef	CLEANUP
}

int pij_gassist_llrtopijs_tot(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,const VECTORF* llr1,const MATRIXF* llr2b,const MATRIXF* llr2c,const MATRIXF* llr3,VECTORF* p1,MATRIXF* p2b,MATRIXF* p2c,MATRIXF* p3,size_t nv,char nodiag)
{
#define	CLEANUP	CLEANVECF(vp2)CLEANVECF(vp3)CLEANVECF(tp)CLEANVECF(tp2)
	int		ret;
	size_t	ng,n;
	double	nratio;		//Ratio of null distribution
	VECTORF *vp2,*vp3,*tp,*tp2;

	//Initialize
	vp2=vp3=tp=tp2=0;
	ng=t->size1;
	n=ng*t2->size1;
	if(nodiag)
		n-=GSL_MIN(ng,t2->size1);
	vp2=VECTORFF(alloc)(n);
	vp3=VECTORFF(alloc)(n);
	tp=VECTORFF(alloc)(n);
	tp2=VECTORFF(alloc)(n);
	if(!(vp2&&vp3&&tp&&tp2))
		ERRRET("Not enough memory.")
	
	if(VECTORFF(first_nan)(llr1)>=0)
		LOG(4,"Infinity found for log likelihood ratio, possibly because data contains fully correlated/anticorrelated columns. This may affect downstream analyses.")
	ret=pij_gassist_llrtopij1_1(p1);
	if(ret)
		ERRRET("Failed to calculate probabilities in step 1.")
		
	if(nodiag)
		MATRIXFF(flatten_nodiag)(llr2b,tp);
	else
		MATRIXFF(flatten)(llr2b,tp);
	if(VECTORFF(first_nan)(tp)>=0)
		LOG(4,"Infinity found for log likelihood ratio, possibly because data contains fully correlated/anticorrelated columns. This may affect downstream analyses.")
	ret=pij_gassist_llrtopij2c_tot(g,tp,vp2,nv,&nratio);
	if(ret)
		ERRRET("Failed to calculate probabilities in step 2.")
	if(nodiag)
		VECTORFF(wrap_nodiag)(vp2,p2b);
	else
		VECTORFF(wrap)(vp2,p2b);
	
	if(nodiag)
		MATRIXFF(flatten_nodiag)(llr2c,tp);
	else
		MATRIXFF(flatten)(llr2c,tp);
	if(VECTORFF(first_nan)(tp)>=0)
		LOG(4,"Infinity found for log likelihood ratio, possibly because data contains fully correlated/anticorrelated columns. This may affect downstream analyses.")
	ret=pij_gassist_llrtopij2c_tot(g,tp,vp2,nv,&nratio);
	if(ret)
		ERRRET("Failed to calculate probabilities in step 2.")
	if(nodiag)
		VECTORFF(wrap_nodiag)(vp2,p2c);
	else
		VECTORFF(wrap)(vp2,p2c);
	
	if(nodiag)
		MATRIXFF(flatten_nodiag)(llr3,tp2);
	else
		MATRIXFF(flatten)(llr3,tp2);
	if(VECTORFF(first_nan)(tp2)>=0)
	{
		LOG(4,"Infinity found for log likelihood ratio, possibly because data contains fully correlated/anticorrelated columns. This may affect downstream analyses.")
		VECTORFF(set_nan)(tp2,0);
	}
	ret=pij_gassist_llrtopij3_tot(g,tp2,vp3,nv,&nratio);
	if(ret)
		ERRRET("Failed to calculate probabilities in step 3.")
	if(nodiag)
		VECTORFF(wrap_nodiag)(vp3,p3);
	else
		VECTORFF(wrap)(vp3,p3);
	
	//Free memory
	CLEANUP
	return 0;
#undef	CLEANUP
}

int pij_gassist_llrtopij_tot(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,const VECTORF* llr1,const MATRIXF* llr2b,const MATRIXF* llr2c,const MATRIXF* llr3,MATRIXF* ansb,MATRIXF* ansc,size_t nv,char nodiag)
{
#define	CLEANUP	CLEANVECF(p1)CLEANMATF(p3)
	int		ret;
	size_t	i,ng;
	VECTORF *p1;
	MATRIXF	*p3;
	VECTORFF(view)	vv;

	//Initialize
	ng=g->size1;
	p1=VECTORFF(alloc)(llr1->size);
	p3=MATRIXFF(alloc)(llr2b->size1,llr2b->size2);
	if(!(p1&&p3))
		ERRRET("Not enough memory.")

	ret=pij_gassist_llrtopijs_tot(g,t,t2,llr1,llr2b,llr2c,llr3,p1,ansb,ansc,p3,nv,nodiag);
	if(ret)
		ERRRET("Failed to convert log likelihood ratios to probabilities.")

	//Combine probabilities to full:
	MATRIXFF(mul_elements)(ansb,p3);
	MATRIXFF(mul_elements)(ansc,p3);
	for(i=0;i<ng;i++)
	{
		vv=MATRIXFF(row)(ansb,i);
		VECTORFF(scale)(&vv.vector,VECTORFF(get)(p1,i));
		vv=MATRIXFF(row)(ansc,i);
		VECTORFF(scale)(&vv.vector,VECTORFF(get)(p1,i));
	}

	//Free memory
	CLEANUP
	return 0;
#undef	CLEANUP
}

















