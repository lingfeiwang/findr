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
#include "nullmodeler.h"
#include "nullsampler.h"
#include "nullhist.h"
#include "llrtopij.h"
#pragma GCC diagnostic ignored "-Wunused-parameter"

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














