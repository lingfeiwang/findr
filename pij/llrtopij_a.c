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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../base/gsl/math.h"
#include "../base/gsl/histogram.h"
#include "../base/gsl/blas.h"
#include "../base/logger.h"
#include "../base/macros.h"
#include "../base/data_process.h"
#include "../base/histogram.h"
#include "nulldist.h"
#include "llrtopij_a.h"
#include "llrtopij.h"

gsl_histogram* pij_llrtopij_a_nullhist_single(double dmax,size_t nd,size_t n1,size_t n2)
{
#define	CLEANUP	CLEANHIST(h)
	struct pij_nulldist_pdfs_param param={n1,n2};
	size_t	nbin;
	gsl_histogram *h=0;
	
	assert(n1&&n2);
	dmax*=(1+1E-6);
	nbin=histogram_unequalbins_param_count(nd);
	if(nbin<5)
		ERRRETV(0,"Determined %lu bins constructed. Bin count too small.",nbin)
	else if(nbin<10)
		LOG(5,"Determined %lu bins, smaller than recommended minimum bin count (10).",nbin)
	else
		LOG(9,"Determined %lu bins.",nbin)
	h=gsl_histogram_alloc(nbin);
	if(!h)
		ERRRETV(0,"Not enough memory.")
	//Null density histogram
	gsl_histogram_set_ranges_uniform(h,0,dmax);
	//Set null histogram ranges
	if(histogram_unequalbins_fromnullpdfs(nbin,h->range,pij_nulldist_pdfs,&param))
		ERRRETV(0,"histogram_unequalbins_fromnullpdfs failed.")
	//Calculate null density histogram
	if(pij_nulldist_hist_pdf(h->range,nbin,h->bin,param.n1,param.n2,5))
		ERRRETV(0,"pij_nulldist_hist_pdf failed.")
	return h;
#undef	CLEANUP
}

gsl_histogram** pij_llrtopij_a_nullhist(double dmax,size_t nv,size_t ns,size_t nd,long n1d,long n2d)
{
#define	CLEANUP	if(h){for(i=0;i<nv-1;i++)CLEANHIST(h[i])free(h);h=0;}
// 	struct pij_nulldist_cdf_callback_param cdfparam;
	struct pij_nulldist_pdfs_param param;
	size_t	nbin,i;
	int		ret;
	gsl_histogram **h;
	
	assert((2+n1d>0)&&((long)(ns-nv)+n2d>0));
	assert(nv>=2);
	dmax*=(1+1E-6);
	h=calloc(nv-1,sizeof(gsl_histogram*));
	if(!h)
		ERRRETV(0,"Not enough memory.")
	nbin=histogram_unequalbins_param_count(nd);
	if(nbin<5)
		ERRRETV(0,"Determined %lu bins constructed. Bin count too small.",nbin)
	else if(nbin<10)
		LOG(5,"Determined %lu bins, smaller than recommended minimum bin count (10).",nbin)
	else
		LOG(9,"Determined %lu bins.",nbin)
	ret=1;
	for(i=0;i<nv-1;i++)
		ret=ret&&(h[i]=gsl_histogram_alloc(nbin));
	if(!ret)
		ERRRETV(0,"Not enough memory.")
	//Null density histogram
	for(i=0;i<nv-1;i++)
	{
		gsl_histogram_set_ranges_uniform(h[i],0,dmax);
		//Set null histogram ranges
		param.n1=(size_t)((long)(i+2)+n1d);
		param.n2=(size_t)((long)(ns-i-2)+n2d);
		if(histogram_unequalbins_fromnullpdfs(nbin,h[i]->range,pij_nulldist_pdfs,&param))
			ERRRETV(0,"histogram_unequalbins_fromnullpdfs failed.")
		//Calculate null density histogram
		if(pij_nulldist_hist_pdf(h[i]->range,nbin,h[i]->bin,param.n1,param.n2,5))
			ERRRETV(0,"pij_nulldist_hist_pdf failed.")
	}
	return h;
#undef	CLEANUP
}

int pij_llrtopij_a_convert_single(const MATRIXF* d,const MATRIXF* dconv,MATRIXF* ans,size_t n1,size_t n2,char nodiag)
{
#define	CLEANUP	CLEANHIST(hreal)CLEANHIST(hc)\
				CLEANHIST(h)CLEANVECD(vb1)CLEANVECD(vb2)
	size_t		ng=d->size1;
	size_t		j,k,nbin;
	gsl_histogram *h,*hreal,*hc;
	VECTORD		*vb1,*vb2;
	
	h=0;
	hreal=hc=0;
	vb1=0;
	vb2=0;
	//Validity checks
	assert((dconv->size1==ng)&&(ans->size1==ng)&&(ans->size2==dconv->size2));

	//Construct null density histograms
	{
		FTYPE		dmin,dmax;
		MATRIXFF(minmax)(d,&dmin,&dmax);
		if((!(dmin>=0))||gsl_isnan(dmax)||gsl_isinf(dmax))
			ERRRET("Negative or NAN found in input data. It may invalidate follow up analysis. This may be due to incorrect previous steps.")
		h=pij_llrtopij_a_nullhist_single((double)dmax,d->size2,n1,n2);
		if(!h)
			ERRRET("pij_llrtopij_a_nullhist_single failed.")
		nbin=h->n;
	}
	
	//Prepare for real histogram
	hreal=gsl_histogram_clone(h);
	hc=gsl_histogram_alloc(hreal->n+2);
	if(!(hreal&&hc))
		ERRRET("Not enough memory.");
	if(pij_llrtopij_convert_histograms_make_buffs(nbin,&vb1,&vb2))
		ERRRET("pij_llrtopij_convert_histograms_make_buffs failed.")
	
	//Conversion
	{
		VECTORDF(view)	vv1,vvreal;
		VECTORFF(view)	vva;
		VECTORD	*vwidth,*vnull;
		vwidth=VECTORDF(alloc)(nbin);
		vnull=VECTORDF(alloc)(nbin);
		if(!(vwidth&&vnull))
		{
			CLEANVECD(vwidth)CLEANVECD(vnull)
			ERRRET("Not enough memory.")
		}
		vv1=VECTORDF(view_array)(h->range+1,nbin);
		VECTORDF(memcpy)(vwidth,&vv1.vector);
		vv1=VECTORDF(view_array)(h->range,nbin);
		VECTORDF(sub)(vwidth,&vv1.vector);
		vvreal=VECTORDF(view_array)(hreal->bin,nbin);
		vv1=VECTORDF(view_array)(h->bin,nbin);
		for(j=0;j<ng;j++)
		{
			VECTORFF(const_view)	vvd=MATRIXFF(const_row)(dconv,j);
			VECTORDF(memcpy)(vnull,&vv1.vector);
			memcpy(hreal->range,h->range,(nbin+1)*sizeof(double));
			memset(hreal->bin,0,nbin*sizeof(double));
			//Construct real histogram
			if(nodiag&&j<d->size2)
			{
				for(k=0;k<j;k++)
					gsl_histogram_increment(hreal,MATRIXFF(get)(d,j,k));
				for(k=j+1;k<d->size2;k++)
					gsl_histogram_increment(hreal,MATRIXFF(get)(d,j,k));
				VECTORDF(scale)(&vvreal.vector,1./(double)(d->size2-1));
			}
			else
			{
				for(k=0;k<d->size2;k++)
					gsl_histogram_increment(hreal,MATRIXFF(get)(d,j,k));
				VECTORDF(scale)(&vvreal.vector,1./(double)(d->size2));
			}
			//Convert to density histogram
			VECTORDF(div)(&vvreal.vector,vwidth);
			//Convert to probability central histogram
			pij_llrtopij_convert_histograms_buffed(hreal,vnull,hc,vb1,vb2);
			//Convert likelihoods to probabilities
			vva=MATRIXFF(row)(ans,j);
			pij_llrtopij_histogram_interpolate_linear(hc,&vvd.vector,&vva.vector);
		}
		CLEANVECD(vwidth)CLEANVECD(vnull)
	}
	CLEANUP
	return 0;
#undef	CLEANUP
}













