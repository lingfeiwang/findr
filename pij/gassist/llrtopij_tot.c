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

int pij_gassist_llrtopij2_tot(const MATRIXG* g,const VECTORF* llr,VECTORF* p,size_t nv,double* nratio)
{
	const struct histogram_unequalbins_exp_param			ph={g->size2-nv,1};
	const struct pij_gassist_nullhist_analytical_pdf_param	pnh={g,nv,5};
	return pij_llrtopij_tot_convert(llr,llr,p,histogram_unequalbins_exp,pij_gassist_nullhist_analytical2_pdf,&ph,&pnh,nratio);
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

int pij_gassist_llrtopij4_tot(const MATRIXG* g,const VECTORF* llr,VECTORF* p,size_t nv,double* nratio)
{
	const struct histogram_unequalbins_exp_param			ph={g->size2-nv,1};
	const struct pij_gassist_nullhist_analytical_pdf_param	pnh={g,nv,5};
	return pij_llrtopij_tot_convert(llr,llr,p,histogram_unequalbins_exp,pij_gassist_nullhist_analytical4_pdf,&ph,&pnh,nratio);
}

int pij_gassist_llrtopij5_tot(const MATRIXG* g,const VECTORF* llr,VECTORF* p,size_t nv,double* nratio)
{
	const struct histogram_unequalbins_exp_param			ph={g->size2-nv,1};
	const struct pij_gassist_nullhist_analytical_pdf_param	pnh={g,nv,5};
	return pij_llrtopij_tot_convert(llr,llr,p,histogram_unequalbins_exp,pij_gassist_nullhist_analytical5_pdf,&ph,&pnh,nratio);
}

int pij_gassist_llrtopijs_tot(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,const VECTORF* llr1,const MATRIXF* llr2,const MATRIXF* llr3,const MATRIXF* llr4,const MATRIXF* llr5,VECTORF* p1,MATRIXF* p2,MATRIXF* p3,MATRIXF* p4,MATRIXF* p5,size_t nv,char nodiag)
{
#define	CLEANUP	CLEANVECF(vp)CLEANVECF(tp)
	int		ret;
	size_t	ng,n;
	double	nratio;		//Ratio of null distribution
	VECTORF *vp,*tp;

	//Initialize
	vp=tp=0;
	ng=t->size1;
	n=ng*t2->size1;
	if(nodiag)
		n-=GSL_MIN(ng,t2->size1);
	vp=VECTORFF(alloc)(n);
	tp=VECTORFF(alloc)(n);
	if(!(vp&&tp))
		ERRRET("Not enough memory.")
	
	if(VECTORFF(first_nan)(llr1)>=0)
		LOG(4,"Infinity found for log likelihood ratio, possibly because data contains fully correlated/anticorrelated columns. This may affect downstream analyses.")
	ret=pij_gassist_llrtopij1_1(p1);
	if(ret)
		ERRRET("Failed to calculate probabilities in step 1.")
		
	if(nodiag)
		MATRIXFF(flatten_nodiag)(llr2,tp);
	else
		MATRIXFF(flatten)(llr2,tp);
	if(VECTORFF(first_nan)(tp)>=0)
		LOG(4,"Infinity found for log likelihood ratio, possibly because data contains fully correlated/anticorrelated columns. This may affect downstream analyses.")
	ret=pij_gassist_llrtopij2_tot(g,tp,vp,nv,&nratio);
	if(ret)
		ERRRET("Failed to calculate probabilities in step 2.")
	if(nodiag)
		VECTORFF(wrap_nodiag)(vp,p2);
	else
		VECTORFF(wrap)(vp,p2);
	
	if(nodiag)
		MATRIXFF(flatten_nodiag)(llr3,tp);
	else
		MATRIXFF(flatten)(llr3,tp);
	if(VECTORFF(first_nan)(tp)>=0)
	{
		LOG(4,"Infinity found for log likelihood ratio, possibly because data contains fully correlated/anticorrelated columns. This may affect downstream analyses.")
		VECTORFF(set_nan)(tp,0);
	}
	ret=pij_gassist_llrtopij3_tot(g,tp,vp,nv,&nratio);
	if(ret)
		ERRRET("Failed to calculate probabilities in step 3.")
	if(nodiag)
		VECTORFF(wrap_nodiag)(vp,p3);
	else
		VECTORFF(wrap)(vp,p3);
	
	if(nodiag)
		MATRIXFF(flatten_nodiag)(llr4,tp);
	else
		MATRIXFF(flatten)(llr4,tp);
	if(VECTORFF(first_nan)(tp)>=0)
		LOG(4,"Infinity found for log likelihood ratio, possibly because data contains fully correlated/anticorrelated columns. This may affect downstream analyses.")
	ret=pij_gassist_llrtopij4_tot(g,tp,vp,nv,&nratio);
	if(ret)
		ERRRET("Failed to calculate probabilities in step 4.")
	if(nodiag)
		VECTORFF(wrap_nodiag)(vp,p4);
	else
		VECTORFF(wrap)(vp,p4);
	
	if(nodiag)
		MATRIXFF(flatten_nodiag)(llr4,tp);
	else
		MATRIXFF(flatten)(llr4,tp);
	if(VECTORFF(first_nan)(tp)>=0)
		LOG(4,"Infinity found for log likelihood ratio, possibly because data contains fully correlated/anticorrelated columns. This may affect downstream analyses.")
	ret=pij_gassist_llrtopij5_tot(g,tp,vp,nv,&nratio);
	if(ret)
		ERRRET("Failed to calculate probabilities in step 5.")
	if(nodiag)
		VECTORFF(wrap_nodiag)(vp,p5);
	else
		VECTORFF(wrap)(vp,p5);
	
	//Free memory
	CLEANUP
	return 0;
#undef	CLEANUP
}
















