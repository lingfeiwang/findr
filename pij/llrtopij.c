/* Copyright 2016-2018, 2020 Lingfei Wang
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
#include <math.h>
#include <float.h>
#include <string.h>
#include "../base/gsl/math.h"
#include "../base/gsl/histogram.h"
#include "../base/gsl/blas.h"
#include "../base/logger.h"
#include "../base/macros.h"
#include "../base/data_process.h"
#include "../base/histogram.h"
#include "../base/threading.h"
#include "nullhist.h"
#include "llrtopij.h"

// #pragma GCC diagnostic ignored "-Wunused-parameter"

/* Smoothen true ratio histogram to make it monotonically increasing
 * Method:	1.	Construct increasing upper bound histogram (hup)
 * 				so every bin is max(self,left neighbor of self)
 * 			2.	Construct increasing lower bound histogram (hlow)
 * 				so every bin is min(self,right neighbor of self)
 * 			3.	Return (hup+hlow)/2
 * h:		[n] double histogram to be smoothened
 * 			Input as histogram of true ratio and,
 * 			output as smoothened histogram of true ratio
 * n:		size of histogram
 * Return:	0 on success.
 */
static void pij_llrtopij_histogram_force_incremental_buffed(gsl_histogram *h,VECTORD* hup)
{
	VECTORDF(view)	hlow=VECTORDF(view_array)(h->bin,h->n);
	size_t	i;
	
	assert(hup->size==h->n);
	VECTORDF(memcpy)(hup,&hlow.vector);
	//Step 1
	for(i=1;i<hup->size;i++)
		VECTORDF(set)(hup,i,GSL_MAX(VECTORDF(get)(hup,i-1),VECTORDF(get)(hup,i)));
	//Step 2
	for(i=hlow.vector.size-1;i;i--)
		VECTORDF(set)(&hlow.vector,i-1,GSL_MIN(VECTORDF(get)(&hlow.vector,i-1),VECTORDF(get)(&hlow.vector,i)));
	//Step 3
	VECTORDF(add)(&hlow.vector,hup);
	VECTORDF(scale)(&hlow.vector,0.5);
}

/* Smoothen 1D histogram with buffer provided as the following:
 * 1. Calculate bin differences.
 * 2. Extrapolates with fixed boundary values to ensure same size of histogram
 * 3. Convolution with Gaussian filter.
 * 4. Calculate cumulative sum
 * 5. Scale and shift to original head-tail positions
 * WARNING:	Gaussian filter applied using bin index as distance measure,
 * 			not the actual histogram range.
 * h:		[n] double histogram to be smoothened
 * n:		size of histogram
 * sigma:	sigma of Gaussian filter
 * ncut:	size of Gaussian convolution vector is 2*ncut+1
 * vlarge:	[n+2*ncut-1] Buffer for data before convolution
 * vconv:	[2*ncut+1] Buffer for convolution mask.
 * Return:	0 on success
 */
static void pij_llrtopij_histogram_smoothen_buffed(gsl_histogram *h,double sigma,size_t ncut,VECTORD *vlarge,VECTORD *vconv)
{
	VECTORDF(view)	vv1,vv2;
	double	tmp,dmin,ddiff,t1,t2;
	size_t	i;

	assert(vlarge->size==h->n+2*ncut-1);
	assert(vconv->size==2*ncut+1);
	//Construct convolution vector
	VECTORDF(set)(vconv,ncut,1);
	for(i=1;i<=ncut;i++)
	{
		tmp=exp(-gsl_pow_2((double)ncut/sigma));
		VECTORDF(set)(vconv,ncut-i,tmp);
		VECTORDF(set)(vconv,ncut+i,tmp);
	}
	
	//Normalize convolution vector
	for(i=0,tmp=0;i<=2*ncut;i++)
		tmp+=VECTORDF(get)(vconv,i);
	VECTORDF(scale)(vconv,1/tmp);
	
	//Construct difference vector
	vv1=VECTORDF(subvector)(vlarge,ncut-1,h->n);
	vv2=VECTORDF(view_array)(h->bin,h->n);
	VECTORDF(minmax)(&vv2.vector,&dmin,&ddiff);
	ddiff-=dmin;
	if(!ddiff)
		return;
	VECTORDF(memcpy)(&vv1.vector,&vv2.vector);
	VECTORDF(diff)(&vv1.vector);

	//Construct extended unconvolved difference vector
	tmp=VECTORDF(get)(vlarge,ncut);
	vv1=VECTORDF(subvector)(vlarge,0,ncut);
	VECTORDF(set_all)(&vv1.vector,tmp);
	tmp=VECTORDF(get)(vlarge,vlarge->size-1-ncut);
	vv1=VECTORDF(subvector)(vlarge,vlarge->size-ncut,ncut);
	VECTORDF(set_all)(&vv1.vector,tmp);
	
	//Convolution
	for(i=1;i<h->n;i++)
	{
		vv1=VECTORDF(subvector)(vlarge,i-1,2*ncut+1);
		gsl_blas_ddot(&vv1.vector,vconv,h->bin+i);
	}
	
	//Cumulative sum and recover original mean and difference
	VECTORDF(set)(&vv2.vector,0,0);
	VECTORDF(cumsum)(&vv2.vector);
	VECTORDF(minmax)(&vv2.vector,&t1,&t2);
	VECTORDF(scale)(&vv2.vector,ddiff/(t2-t1));
	VECTORDF(add_constant)(&vv2.vector,dmin-t1*ddiff/(t2-t1));
}

/* Construct central value histogram from bounded histogram for interpolation.
 * In central value histogram, bin[i] is the value at range[i].
 * h:		(nbin) Input bounded histogram
 * hc:		(nbin+2) Output central value histogram
 */
static void pij_llrtopij_histogram_to_central(const gsl_histogram *h,gsl_histogram* hc)
{
	size_t	n;
	VECTORDF(view)	vv1,vv2;
	
	assert(hc->n==h->n+2);
	n=h->n;
	memcpy(hc->range+1,h->range+1,n*sizeof(*hc->range));
	vv1=VECTORDF(view_array)(hc->range+1,n+1);
	vv2=VECTORDF(view_array)(h->range,n+1);
	VECTORDF(add)(&vv1.vector,&vv2.vector);
	VECTORDF(scale)(&vv1.vector,0.5);
	hc->range[0]=h->range[0];
	hc->range[n+1]=h->range[n];
	hc->range[n+2]=hc->range[n+1]+fabs(hc->range[n+1]);
	memcpy(hc->bin+1,h->bin,h->n*sizeof(*hc->bin));
	//Use fixed boundary condition
	hc->bin[0]=hc->bin[1];
	hc->bin[h->n+1]=hc->bin[h->n];
}

void pij_llrtopij_histogram_interpolate_linear(const gsl_histogram *hc,const VECTORF* d,VECTORF* ans)
{
	size_t	i;
	size_t	loc;
	FTYPE	f;
	
	for(i=0;i<d->size;i++)
	{
		f=VECTORFF(get)(d,i);
		if(f<=hc->range[0])
			VECTORFF(set)(ans,i,(FTYPE)hc->bin[0]);
		else if(f>=hc->range[hc->n])
			VECTORFF(set)(ans,i,(FTYPE)hc->bin[hc->n-1]);
		else
		{
			gsl_histogram_find(hc,f,&loc);
			VECTORFF(set)(ans,i,(FTYPE)(hc->bin[loc]+(f-hc->range[loc])*(hc->bin[loc+1]-hc->bin[loc])/(hc->range[loc+1]-hc->range[loc])));
		}
	}
}

void pij_llrtopij_convert_histograms_get_buff_sizes(size_t n,size_t *n1,size_t *n2)
{
	size_t ncut;

	ncut=GSL_MIN(n/2,50);
	*n1=n+2*ncut-1;
	*n2=2*ncut+1;
}

int pij_llrtopij_convert_histograms_make_buffs(size_t n,VECTORD** vb1,VECTORD** vb2)
{
#define	CLEANUP	CLEANVECD(*vb1)CLEANVECD(*vb2)
	size_t n1,n2;
	
	pij_llrtopij_convert_histograms_get_buff_sizes(n,&n1,&n2);
	*vb1=VECTORDF(alloc)(n1);
	*vb2=VECTORDF(alloc)(n2);
	if(!(*vb1&&*vb2))
		ERRRET("Not enough memory.")
	return 0;
#undef	CLEANUP
}

void pij_llrtopij_convert_histograms_buffed(gsl_histogram* hreal,VECTORD* vnull,gsl_histogram* hc,VECTORD* vb1,VECTORD* vb2)
{
	size_t	minpos;		//Position of minimum true density		
	double	ndens;
	double	psmooth;
	size_t	nbin=hreal->n,ncut,nb;
	size_t	i;
	VECTORDF(view)	vv1,vvreal,vv2,vvnull,vv3;

	//Trim trailing 0s in histogram	
	for(nb=nbin-1;nb&&(!hreal->bin[nb]);nb--);
	nb++;
	//If all in smallest histogram, then output zeros
	if(nb==1)
	{
		vv1=VECTORDF(view_array)(hreal->bin,hreal->n);
		VECTORDF(set_all)(&vv1.vector,0);
		pij_llrtopij_histogram_to_central(hreal,hc);
		return;
	}

	ncut=GSL_MIN(nb/2,50);
	assert(vb1->size>=hreal->n+2*ncut-1);
	assert(vb2->size>=2*ncut+1);
	assert((vnull->size==nbin)&&(hc->n==nbin+2));
	vvreal=VECTORDF(view_array)(hreal->bin,nb);
	vvnull=VECTORDF(subvector)(vnull,0,nb);
	for(i=1;i<nb;i++)
		if(!hreal->bin[i])
			hreal->bin[i]=hreal->bin[i-1];

	//Find minimum location
	VECTORDF(div)(&vvnull.vector,&vvreal.vector);
	//Fill zeros with the one before
	minpos=VECTORDF(max_index)(&vvnull.vector);
	ndens=1/VECTORDF(get)(&vvnull.vector,minpos);
	if(ndens==0)
	{
		vv1=VECTORDF(view_array)(hreal->bin,hreal->n);
		VECTORDF(set_all)(&vv1.vector,1);
		pij_llrtopij_histogram_to_central(hreal,hc);
		return;
	}
	//Calculate true distribution ratio
	VECTORDF(scale)(&vvnull.vector,-ndens);
	VECTORDF(add_constant)(&vvnull.vector,1);
	if(nb<nbin)
	{
		vv1=VECTORDF(subvector)(vnull,nb,nbin-nb);
		VECTORDF(set_all)(&vv1.vector,VECTORDF(get)(&vvnull.vector,nb-1));
	}
	
	//Trim true ratio so bin location<=minpos all become zero
	vv1=VECTORDF(subvector)(&vvnull.vector,0,minpos+1);
	VECTORDF(set_zero)(&vv1.vector);
	vvreal=VECTORDF(view_array)(hreal->bin,nbin);
	VECTORDF(memcpy)(&vvreal.vector,vnull);
	
	//Convert true ratio to monotonic function and smoothening
	vv3=VECTORDF(subvector)(vb1,0,hreal->n);
	pij_llrtopij_histogram_force_incremental_buffed(hreal,&vv3.vector);
	//Slightly smoothen before convert to monotonic function
	psmooth=GSL_MIN((double)ncut/3.,30.);
	vv1=VECTORDF(subvector)(vb1,0,hreal->n+2*ncut-1);
	vv2=VECTORDF(subvector)(vb2,0,2*ncut+1);
	pij_llrtopij_histogram_smoothen_buffed(hreal,psmooth,ncut,&vv1.vector,&vv2.vector);
	//Convert bounded histogram to central histogram
	pij_llrtopij_histogram_to_central(hreal,hc);
	
}

int pij_llrtopij_convert_histograms(gsl_histogram* hreal,VECTORD* vnull,gsl_histogram* hc)
{
#define	CLEANUP
	VECTORD	*vb1,*vb2;
	if(pij_llrtopij_convert_histograms_make_buffs(hreal->n,&vb1,&vb2))
		ERRRET("pij_llrtopij_convert_histograms_make_buffs failed.")
	pij_llrtopij_convert_histograms_buffed(hreal,vnull,hc,vb1,vb2);
	CLEANVECD(vb1)
	CLEANVECD(vb2)
	return 0;
#undef	CLEANUP
}

FTYPE pij_llrtopij_llrmatmax(MATRIXF* d,char nodiag)
{
	FTYPE	dmin,dmax;
	if(nodiag)
		MATRIXFF(minmax_nodiag)(d,&dmin,&dmax,0);
	else
		MATRIXFF(minmax)(d,&dmin,&dmax);
	if((!(dmin>=0))||gsl_isnan(dmax))
	{
		LOG(1,"Negative or NAN found in input data. It may invalidate follow up analysis. This may be due to incorrect previous steps.")
		return 0;
	}
	if(gsl_isinf(dmax))
	{
		LOG(5,"INF found in input data. It may invalidate follow up analysis. This may be due to incorrect previous steps or duplicate rows. Now regard INFs as largest non-INF value.")
		MATRIXFF(set_inf)(d,-1);
		if(nodiag)
			MATRIXFF(minmax_nodiag)(d,&dmin,&dmax,0);
		else
			dmax=MATRIXFF(max)(d);
		MATRIXFF(set_value)(d,-1,dmax);
	}
	return dmax;
}

FTYPE pij_llrtopij_llrvecmax(VECTORF* d)
{
	FTYPE	dmin,dmax;
	VECTORFF(minmax)(d,&dmin,&dmax);
	if((!(dmin>=0))||gsl_isnan(dmax))
	{
		LOG(1,"Negative or NAN found in input data. It may invalidate follow up analysis. This may be due to incorrect previous steps.")
		return 0;
	}
	if(gsl_isinf(dmax))
	{
		LOG(5,"INF found in input data. It may invalidate follow up analysis. This may be due to incorrect previous steps or duplicate rows. Now regard INFs as largest non-INF value.")
		VECTORFF(set_inf)(d,-1);
		dmax=VECTORFF(max)(d);
		VECTORFF(set_value)(d,-1,dmax);
	}
	return dmax;
}


int pij_llrtopij_convert_single(const MATRIXF* d,const MATRIXF* dconv,MATRIXF* ans,size_t n1,size_t n2,char nodiag,long nodiagshift)
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
		h=pij_nullhist_single((double)dmax,d->size2,n1,n2);
		if(!h)
			ERRRET("pij_nullhist_single failed.")
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

int pij_llrtopij_convert_single_self(MATRIXF* d,size_t n1,size_t n2,char nodiag,long nodiagshift)
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
		h=pij_nullhist_single((double)dmax,d->size2,n1,n2);
		if(!h)
			ERRRET("pij_nullhist_single failed.")
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




















