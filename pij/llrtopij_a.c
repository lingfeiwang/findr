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
#include "../base/config.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../base/threading.h"
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

/* Construct one null histogram for a specific genotype value count.
 * The function calculates the null density histogram for random variable:
 * x=-0.5*log(1-z1/(z1+z2)), where z1~chi2(n1),z2~chi2(n2),
 * Histogram bin count and width are automatically determined
 * from real data count (nd).
 * For bin range settings, see histogram_unequalbins_fromnullcdf.
 * For null density histogram from pdf, see pij_nulldist_hist_pdf.
 * dmax:	Specifies the histogram bound as [0,dmax).
 * nd:		Count of real data to form real histograms. This is used to
 * 			automatically decide number of bins and widths.
 * n1,
 * n2:		Parameters of null distribution.
 * Return:	Constructed null distribution histograms with preset
 * 			bin ranges and values as density.
 */
static gsl_histogram* pij_llrtopij_a_nullhist_single(double dmax,size_t nd,size_t n1,size_t n2)
{
#define	CLEANUP	CLEANHIST(h)
	struct pij_nulldist_pdfs_param param={n1,n2};
	size_t	nbin;
	gsl_histogram *h=0;
	
	assert(n1&&n2);
	dmax*=(1+1E-6);
	nbin=histogram_unequalbins_param_count(nd);
	if(nbin<5)
		ERRRETV(0,"Determined "PRINTFSIZET" bins constructed. Bin count too small.",nbin)
	else if(nbin<10)
		LOG(5,"Determined "PRINTFSIZET" bins, smaller than recommended minimum bin count (10).",nbin)
	else
		LOG(10,"Determined "PRINTFSIZET" bins.",nbin)
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

gsl_histogram** pij_llrtopij_a_nullhist(double dmax,size_t nv,size_t nd,long n1c,size_t n1d,long n2c,size_t n2d)
{
#define	CLEANUP	if(h){for(i=0;i<nv-1;i++)CLEANHIST(h[i])free(h);h=0;}
	struct pij_nulldist_pdfs_param param;
	size_t	nbin,i;
	int		ret;
	gsl_histogram **h;
	
	assert(nv>=2);
	dmax*=(1+1E-6);
	CALLOCSIZE(h,nv-1);
	if(!h)
		ERRRETV(0,"Not enough memory.")
	nbin=histogram_unequalbins_param_count(nd);
	if(nbin<5)
		ERRRETV(0,"Determined "PRINTFSIZET" bins constructed. Bin count too small.",nbin)
	else if(nbin<10)
		LOG(5,"Determined "PRINTFSIZET" bins, smaller than recommended minimum bin count (10).",nbin)
	else
		LOG(10,"Determined "PRINTFSIZET" bins.",nbin)
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
		param.n1=(size_t)((long)i*n1c+(long)n1d);
		param.n2=(size_t)(-(long)i*n2c+(long)n2d);
		if(histogram_unequalbins_fromnullpdfs(nbin,h[i]->range,pij_nulldist_pdfs,&param))
			ERRRETV(0,"histogram_unequalbins_fromnullpdfs failed.")
		//Calculate null density histogram
		if(pij_nulldist_hist_pdf(h[i]->range,nbin,h[i]->bin,param.n1,param.n2,5))
			ERRRETV(0,"pij_nulldist_hist_pdf failed.")
	}
	return h;
#undef	CLEANUP
}

int pij_llrtopij_a_convert_single(const MATRIXF* d,const MATRIXF* dconv,MATRIXF* ans,size_t n1,size_t n2,char nodiag,long nodiagshift)
{
#define	CLEANUP	CLEANHIST(hreal)CLEANHIST(hc)\
				CLEANHIST(h)CLEANVECD(vb1)CLEANVECD(vb2)
	size_t		ng=d->size1;
	long		k;
	size_t		j,nbin;
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
			memcpy(hreal->range,h->range,(nbin+1)*sizeof(*hreal->range));
			memset(hreal->bin,0,nbin*sizeof(*hreal->bin));
			//Construct real histogram
			if(nodiag&&((long)j+nodiagshift>=0)&&((long)j+nodiagshift<(long)d->size2))
			{
				for(k=(long)j+nodiagshift-1;k>=0;k--)
					gsl_histogram_increment(hreal,MATRIXFF(get)(d,j,(size_t)k));
				for(k=(long)j+nodiagshift+1;k<(long)d->size2;k++)
					gsl_histogram_increment(hreal,MATRIXFF(get)(d,j,(size_t)k));
				VECTORDF(scale)(&vvreal.vector,1./(double)(d->size2-1));
			}
			else
			{
				for(k=0;k<(long)d->size2;k++)
					gsl_histogram_increment(hreal,MATRIXFF(get)(d,j,(size_t)k));
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

int pij_llrtopij_a_convert_single_self(MATRIXF* d,size_t n1,size_t n2,char nodiag,long nodiagshift)
{
#define	CLEANUP	CLEANAMHIST(hreal,nth)CLEANAMHIST(hc,nth)\
				CLEANHIST(h)CLEANMATD(mb1)CLEANMATD(mb2)CLEANMATD(mnull)CLEANMATF(mb3)CLEANVECD(vwidth)
		
	size_t		ng=d->size1;
	size_t		i,nbin;
	gsl_histogram	*h;
	MATRIXD		*mb1,*mb2,*mnull;
	MATRIXF		*mb3;
	VECTORD		*vwidth;
	VECTORDF(view)	vv1;
	size_t		nth;
	
	h=0;
	mb1=mb2=mnull=0;
	mb3=0;
	vwidth=0;
	//Validity checks
	{
		int	nth0=omp_get_max_threads();
		assert(nth0>0);
		nth=(size_t)nth0;
	}
	
	
	AUTOCALLOC(gsl_histogram*,hreal,nth,64)
	AUTOCALLOC(gsl_histogram*,hc,nth,64)
	if(!(hreal&&hc))
		ERRRET("Not enough memory.");

	//Construct null density histograms
	{
		FTYPE		dmin,dmax;
		if(nodiag)
			MATRIXFF(minmax_nodiag)(d,&dmin,&dmax,nodiagshift);
		else
			MATRIXFF(minmax)(d,&dmin,&dmax);
		if((!(dmin>=0))||gsl_isnan(dmax))
			ERRRET("Negative or NAN found in input data. It may invalidate follow up analysis. This may be due to incorrect previous steps.")
		if(gsl_isinf(dmax))
		{
			LOG(5,"INF found in input data. It may invalidate follow up analysis. This may be due to incorrect previous steps or duplicate rows (by Spearman correlation).")
			MATRIXFF(set_inf)(d,-1);
			if(nodiag)
				MATRIXFF(minmax_nodiag)(d,&dmin,&dmax,nodiagshift);
			else
				MATRIXFF(minmax)(d,&dmin,&dmax);
			MATRIXFF(set_value)(d,-1,dmax);
		}
		h=pij_llrtopij_a_nullhist_single((double)dmax,d->size2,n1,n2);
		if(!h)
			ERRRET("pij_llrtopij_a_nullhist_single failed.")
		nbin=h->n;
	}
	//Memory allocation
	{
		size_t	n1,n2;
		pij_llrtopij_convert_histograms_get_buff_sizes(nbin,&n1,&n2);
		mb1=MATRIXDF(alloc)(nth,n1);
		mb2=MATRIXDF(alloc)(nth,n2);
		mnull=MATRIXDF(alloc)(nth,nbin);
		mb3=MATRIXFF(alloc)(nth,d->size2);
		vwidth=VECTORDF(alloc)(nbin);
		if(!(mb1&&mb2&&mnull&&mb3&&vwidth))
			ERRRET("Not enough memory.")
	}
	
	//Prepare for real histogram
	{
		int	ret;
		for(i=0,ret=1;i<nth;i++)
		{
			hreal[i]=gsl_histogram_clone(h);
			hc[i]=gsl_histogram_alloc(hreal[i]->n+2);
			ret=ret&&hreal[i]&&hc[i];
		}
		if(!ret)
			ERRRET("Not enough memory.");
	}
	
	//Conversion
	vv1=VECTORDF(view_array)(h->range+1,nbin);
	VECTORDF(memcpy)(vwidth,&vv1.vector);
	vv1=VECTORDF(view_array)(h->range,nbin);
	VECTORDF(sub)(vwidth,&vv1.vector);
	vv1=VECTORDF(view_array)(h->bin,nbin);
	#pragma omp parallel
	{
		size_t	ng1,ng2,id;
		size_t	j;
		long	k;
		VECTORDF(view)	vvreal,vvnull,vvb1,vvb2;
		VECTORFF(view)	vvb3,vva;
	
		id=(size_t)omp_get_thread_num();
		vvreal=VECTORDF(view_array)(hreal[id]->bin,nbin);
		vvnull=MATRIXDF(row)(mnull,id);
		vvb1=MATRIXDF(row)(mb1,id);
		vvb2=MATRIXDF(row)(mb2,id);
		vvb3=MATRIXFF(row)(mb3,id);
		
		threading_get_startend(ng,&ng1,&ng2);
		
		for(j=ng1;j<ng2;j++)
		{
			MATRIXFF(get_row)(&vvb3.vector,d,j);
			VECTORDF(memcpy)(&vvnull.vector,&vv1.vector);
			memcpy(hreal[id]->range,h->range,(nbin+1)*sizeof(*hreal[id]->range));
			memset(hreal[id]->bin,0,nbin*sizeof(*hreal[id]->bin));
			//Construct real histogram
			if(nodiag&&((long)j+nodiagshift>=0)&&((long)j+nodiagshift<(long)d->size2))
			{
				for(k=(long)j+nodiagshift-1;k>=0;k--)
					gsl_histogram_increment(hreal[id],MATRIXFF(get)(d,j,(size_t)k));
				for(k=(long)j+nodiagshift+1;k<(long)d->size2;k++)
					gsl_histogram_increment(hreal[id],MATRIXFF(get)(d,j,(size_t)k));
				VECTORDF(scale)(&vvreal.vector,1./(double)(d->size2-1));
			}
			else
			{
				for(k=0;k<(long)d->size2;k++)
					gsl_histogram_increment(hreal[id],MATRIXFF(get)(d,j,(size_t)k));
				VECTORDF(scale)(&vvreal.vector,1./(double)(d->size2));
			}					

			//Convert to density histogram
			VECTORDF(div)(&vvreal.vector,vwidth);
			//Convert to probability central histogram
			pij_llrtopij_convert_histograms_buffed(hreal[id],&vvnull.vector,hc[id],&vvb1.vector,&vvb2.vector);
			//Convert likelihoods to probabilities
			vva=MATRIXFF(row)(d,j);
			pij_llrtopij_histogram_interpolate_linear(hc[id],&vvb3.vector,&vva.vector);
		}
	}
	
	CLEANUP
	return 0;
#undef	CLEANUP
}



























































