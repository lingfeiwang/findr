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
#include <assert.h>
#include <stdint.h>
#include "../../base/gsl/blas.h"
#include "../../base/gsl/math.h"
#include "../../base/gsl/sf.h"
#include "../../base/logger.h"
#include "../../base/macros.h"
#include "../../base/histogram.h"
#include "../../base/math.h"
#include "../../base/data_process.h"
#include "../nulldist.h"
#include "nulldist.h"
#pragma GCC diagnostic ignored "-Wunused-parameter"



/*************************************************************
 * Generic functions for any step
 *************************************************************/


void pij_gassist_nulldist_mixed_pdf_buffed(const MATRIXG* g,size_t nv,const VECTORD* loc,VECTORD* ans,VECTORD* vb2,VECTORG* vb3,VECTORD* vb4,VECTORUC* vb5,MATRIXD* mb1,void (*func)(size_t,size_t,const VECTORD*,MATRIXD*,VECTORD*))
{
	size_t	ns=g->size2;

	assert(loc->size&&(ans->size==loc->size)&&(vb2->size==loc->size)&&(mb1->size2==loc->size));
	assert((nv>1)&&(mb1->size1==(nv-1)));
	
	//First calculate separate log pdfs
	func(ns,nv,loc,mb1,vb2);

	//Count alleles for each genotype
	MATRIXGF(countv_byrow_buffed)(g,vb3,vb5);
	VECTORGF(count_ratio_d)(vb3,vb4);
	
	//Calculate averaged pdf w.r.t genotype
	if(VECTORDF(get)(vb4,1)!=0)
		LOG(1,"Found genotype that has single value across all samples. Ignored.")
	
	VECTORDF(const_view) vvc=VECTORDF(const_subvector)(vb4,2,(size_t)nv-1);
	gsl_blas_dgemv(CblasTrans,1/(1-VECTORDF(get)(vb4,1)),mb1,&vvc.vector,0,ans);
}

int pij_gassist_nulldist_mixed_pdf(const MATRIXG* g,size_t nv,const VECTORD* loc,VECTORD* ans,void (*func)(size_t,size_t,const VECTORD*,MATRIXD*,VECTORD*))
{
#define CLEANUP	AUTOFREEVEC(vb2)AUTOFREEVEC(vb4)AUTOFREEVEC(vb4)
	size_t		ng=g->size1,nd=loc->size;

	VECTORG *vb3=VECTORGF(alloc)(ng);
	MATRIXD	*mb1=MATRIXDF(alloc)(nv-1,nd);
	AUTOALLOCVECD(vb2,nd,10000)
	AUTOALLOCVECD(vb4,nv+1,1000)
	AUTOALLOCVECUC(vb5,nv,1000)
	
	if(!(vb2&&vb3&&vb4&&vb5&&mb1))
		ERRRET("Not enough memory.")
	
	pij_gassist_nulldist_mixed_pdf_buffed(g,nv,loc,ans,vb2,vb3,vb4,vb5,mb1,func);
	CLEANUP
	return 0;
#undef	CLEANUP
}

int pij_gassist_nulldist_nullhist_mixed_pdf0(const MATRIXG* g,size_t nv,const double* restrict range,size_t nbin,double* restrict hist,void* param,void (*func)(size_t,size_t,const VECTORD*,MATRIXD*,VECTORD*))
{
#define CLEANUP	AUTOFREEVEC(vloc)
	int	ret;
	VECTORDF(view)	vvh=VECTORDF(view_array)(hist,nbin);
	AUTOALLOCVECD(vloc,nbin,1000)
	
	assert(!param);
	if(!vloc)
		ERRRET("Not enough memory.")
	//Obtain bin central values
	{
		VECTORDF(const_view) vv1=VECTORDF(const_view_array)(range+1,nbin);
		VECTORDF(memcpy)(vloc,&vv1.vector);
	}
	{
		VECTORDF(const_view) vv1=VECTORDF(const_view_array)(range,nbin);
		VECTORDF(add)(vloc,&vv1.vector);
	}
	VECTORDF(scale)(vloc,0.5);
	ret=pij_gassist_nulldist_mixed_pdf(g,nv,vloc,&vvh.vector,func);
	CLEANUP
	return ret;
#undef	CLEANUP
}

int pij_gassist_nulldist_nullhist_mixed_pdf(const MATRIXG* g,size_t nv,const double* restrict range,size_t nbin,double* restrict hist,void* param,void (*func)(size_t,size_t,const VECTORD*,MATRIXD*,VECTORD*))
{
#define CLEANUP	CLEANVECD(loc)CLEANVECD(val)
	VECTORD	*loc,*val;
	VECTORDF(view)	vvh=VECTORDF(view_array)(hist,nbin);
	size_t		*n=param;
	size_t		i;
	size_t		nsp;
	
	assert(n&&(*n<10));
	if(!*n)
		return pij_gassist_nulldist_nullhist_mixed_pdf0(g,nv,range,nbin,hist,param,func);
	nsp=(size_t)1<<(*n-1);
	loc=VECTORDF(alloc)(nbin*nsp);
	val=VECTORDF(alloc)(nbin*nsp);
	if(!(loc&&val))
		ERRRET("Not enough memory.")
	
	//Construct bin ranges
	{
		VECTORDF(const_view) vvc=VECTORDF(const_view_array)(range,nbin+1);
		histogram_finer_central(&vvc.vector,loc,nsp);
	}
	
	//Calculate bin values
	if(pij_gassist_nulldist_mixed_pdf(g,nv,loc,val,func))
	{
		CLEANUP
		return 1;
	}
	
	//Shrink to output
	VECTORDF(set_zero)(&vvh.vector);
	for(i=0;i<nsp;i++)
	{
		VECTORDF(const_view)	vvc=VECTORDF(const_subvector_with_stride)(val,i,nsp,nbin);
		VECTORDF(add)(&vvh.vector,&vvc.vector);
	}
	VECTORDF(scale)(&vvh.vector,1/gsl_blas_dasum(&vvh.vector));

	CLEANUP
	return 0;
#undef	CLEANUP
}


/*************************************************************
 * Specific functions for each step
 *************************************************************/

void pij_gassist_nulldist1_calcpdf_buffed(size_t ns,size_t nv,const VECTORD* loc,MATRIXD* ans,VECTORD* vb2)
{
	pij_nulldist_calcpdf_buffed(1,nv,ns-1,loc,ans,vb2);
}

void pij_gassist_nulldist2b_calcpdf_buffed(size_t ns,size_t nv,const VECTORD* loc,MATRIXD* ans,VECTORD* vb2)
{
 	pij_nulldist_calcpdf_buffed(2,nv+1,ns-1,loc,ans,vb2);
}

void pij_gassist_nulldist3_calcpdf_buffed(size_t ns,size_t nv,const VECTORD* loc,MATRIXD* ans,VECTORD* vb2)
{
	pij_nulldist_calcpdf_buffed(1,nv,ns-2,loc,ans,vb2);
}

int pij_gassist_nulldist1_mixed_pdf(const VECTORD* loc,VECTORD* ans,const void* param)
{
	const struct pij_gassist_nulldist_mixed_pdf_data	*p=param;
	return pij_gassist_nulldist_mixed_pdf(p->g,p->nv,loc,ans,pij_gassist_nulldist1_calcpdf_buffed);
}

int pij_gassist_nulldist2b_mixed_pdf(const VECTORD* loc,VECTORD* ans,const void* param)
{
	const struct pij_gassist_nulldist_mixed_pdf_data	*p=param;
	return pij_gassist_nulldist_mixed_pdf(p->g,p->nv,loc,ans,pij_gassist_nulldist2b_calcpdf_buffed);
}

int pij_gassist_nulldist3_mixed_pdf(const VECTORD* loc,VECTORD* ans,const void* param)
{
	const struct pij_gassist_nulldist_mixed_pdf_data	*p=param;
	return pij_gassist_nulldist_mixed_pdf(p->g,p->nv,loc,ans,pij_gassist_nulldist3_calcpdf_buffed);
}

int pij_gassist_nulldist_nullhist1_pdf(const MATRIXG* g,size_t nv,const double* restrict range,size_t nbin,double* restrict hist,void* param)
{
	return pij_gassist_nulldist_nullhist_mixed_pdf(g,nv,range,nbin,hist,param,pij_gassist_nulldist1_calcpdf_buffed);
}

int pij_gassist_nulldist_nullhist2b_mixed_pdf(const MATRIXG* g,size_t nv,const double* restrict range,size_t nbin,double* restrict hist,void* param)
{
	return pij_gassist_nulldist_nullhist_mixed_pdf(g,nv,range,nbin,hist,param,pij_gassist_nulldist2b_calcpdf_buffed);
}

int pij_gassist_nulldist_nullhist3_mixed_pdf(const MATRIXG* g,size_t nv,const double* restrict range,size_t nbin,double* restrict hist,void* param)
{
	return pij_gassist_nulldist_nullhist_mixed_pdf(g,nv,range,nbin,hist,param,pij_gassist_nulldist3_calcpdf_buffed);
}




















