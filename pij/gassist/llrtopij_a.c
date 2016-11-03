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
#include <string.h>
#include "../../base/gsl/math.h"
#include "../../base/gsl/histogram.h"
#include "../../base/gsl/blas.h"
#include "../../base/logger.h"
#include "../../base/threading.h"
#include "../../base/macros.h"
#include "../../base/data_process.h"
#include "../../base/histogram.h"
#include "../nulldist.h"
#include "../llrtopij.h"
#include "llrtopij_a.h"
#include "llrtopij.h"
#pragma GCC diagnostic ignored "-Wunused-parameter"

int pij_gassist_llrtopij_a_convert(const MATRIXF* d,const MATRIXF* dconv,MATRIXF* ans,const MATRIXG* g,size_t nv,long n1c,size_t n1d,long n2c,size_t n2d,char nodiag,long nodiagshift)
{
#define	CLEANUP	CLEANVECG(vcount)CLEANAMHIST(hreal,nth)CLEANAMHIST(hc,nth)\
				CLEANMHIST(h,(nv-1))CLEANMATD(mb1)CLEANMATD(mb2)CLEANMATD(mnull)CLEANVECD(vwidth)
		
	VECTORG		*vcount;
	size_t		ng=g->size1;
	size_t		i,nbin;
	gsl_histogram **h;
	//gsl_histogram **hreal,**hc;
	MATRIXD		*mb1,*mb2,*mnull;
	VECTORD		*vwidth;
	VECTORDF(view)	vv1;
	size_t		nth;
	
	h=0;
	mb1=mb2=mnull=0;
	vwidth=0;
	vcount=0;
	//Validity checks
	assert((dconv->size1==ng)&&(ans->size1==ng)&&(d->size1==ng)&&(dconv->size2==ans->size2));

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
		MATRIXFF(minmax)(d,&dmin,&dmax);
		if((!(dmin>=0))||gsl_isnan(dmax)||gsl_isinf(dmax))
			ERRRET("Negative or NAN found in input data. It may invalidate follow up analysis. This may be due to incorrect previous steps.")
		h=pij_llrtopij_a_nullhist((double)dmax,nv,d->size2,n1c,n1d,n2c,n2d);
		if(!h)
			ERRRET("pij_llrtopij_a_nullhist failed.")
		nbin=h[0]->n;
	}
	//Memory allocation
	{
		size_t	n1,n2;
		pij_llrtopij_convert_histograms_get_buff_sizes(nbin,&n1,&n2);
		mb1=MATRIXDF(alloc)(nth,n1);
		mb2=MATRIXDF(alloc)(nth,n2);
		mnull=MATRIXDF(alloc)(nth,nbin);
		vwidth=VECTORDF(alloc)(nbin);
		if(!(mb1&&mb2&&mnull&&vwidth))
			ERRRET("Not enough memory.")
	}
	
	//Prepare for real histogram
	{
		int	ret;
		for(i=0,ret=1;i<nth;i++)
		{
			hreal[i]=gsl_histogram_clone(h[0]);
			hc[i]=gsl_histogram_alloc(hreal[i]->n+2);
			ret=ret&&hreal[i]&&hc[i];
		}
		vcount=VECTORGF(alloc)(ng);
		if(!(ret&&vcount))
			ERRRET("Not enough memory.");
	}
	
	{
		VECTORUC	*vb4=VECTORUCF(alloc)(nv);
		if(!vb4)
			ERRRET("Not enough memory.");	
		MATRIXGF(countv_byrow_buffed)(g,vcount,vb4);
		CLEANVECUC(vb4)
	}
	
	
	//Conversion
	vv1=VECTORDF(view_array)(h[0]->range+1,nbin);
	VECTORDF(memcpy)(vwidth,&vv1.vector);
	vv1=VECTORDF(view_array)(h[0]->range,nbin);
	VECTORDF(sub)(vwidth,&vv1.vector);
	
	
	//Conversion
	for(i=2;i<=nv;i++)
	{
		vv1=VECTORDF(view_array)(h[i-2]->range+1,nbin);
		VECTORDF(memcpy)(vwidth,&vv1.vector);
		vv1=VECTORDF(view_array)(h[i-2]->range,nbin);
		VECTORDF(sub)(vwidth,&vv1.vector);
		vv1=VECTORDF(view_array)(h[i-2]->bin,nbin);
		#pragma omp parallel
		{
			size_t	ng1,ng2,id;
			size_t	j;
			long	k;
			VECTORDF(view)	vvreal,vvnull,vvb1,vvb2;
			VECTORFF(view)	vva;
		
			id=(size_t)omp_get_thread_num();
			vvreal=VECTORDF(view_array)(hreal[id]->bin,nbin);
			vvnull=MATRIXDF(row)(mnull,id);
			vvb1=MATRIXDF(row)(mb1,id);
			vvb2=MATRIXDF(row)(mb2,id);
			threading_get_startend(ng,&ng1,&ng2);
			
			for(j=ng1;j<ng2;j++)
				if(VECTORGF(get)(vcount,j)==i)
				{
			//MATRIXFF(get_row)(&vvb3.vector,d,j);
			
					VECTORFF(const_view)	vvd=MATRIXFF(const_row)(dconv,j);
					VECTORDF(memcpy)(&vvnull.vector,&vv1.vector);
					memcpy(hreal[id]->range,h[i-2]->range,(nbin+1)*sizeof(*hreal[id]->range));
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
					vva=MATRIXFF(row)(ans,j);
					pij_llrtopij_histogram_interpolate_linear(hc[id],&vvd.vector,&vva.vector);
				}
		}
	}
	CLEANUP
	return 0;
#undef	CLEANUP
}

int pij_gassist_llrtopij_a_convert_self(MATRIXF* d,const MATRIXG* g,size_t nv,long n1c,size_t n1d,long n2c,size_t n2d,char nodiag,long nodiagshift)
{
#define	CLEANUP	CLEANVECG(vcount)CLEANAMHIST(hreal,nth)CLEANAMHIST(hc,nth)\
				CLEANMHIST(h,(nv-1))CLEANMATD(mb1)CLEANMATD(mb2)CLEANMATD(mnull)CLEANMATF(mb3)CLEANVECD(vwidth)
		
	VECTORG		*vcount;
	size_t		ng=g->size1;
	size_t		i,nbin;
	gsl_histogram **h;
	//gsl_histogram **hreal,**hc;
	MATRIXD		*mb1,*mb2,*mnull;
	MATRIXF		*mb3;
	VECTORD		*vwidth;
	VECTORDF(view)	vv1;
	size_t		nth;
	
	h=0;
	mb1=mb2=mnull=0;
	mb3=0;
	vwidth=0;
	vcount=0;
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
		FTYPE		dmin=FTYPE_MAX,dmax=FTYPE_MIN;
		if(nodiag)
		{
			FTYPE	dmin1,dmax1;
			size_t	j;
			long	tl;
			VECTORFF(view)	vvl;
			for(j=0;j<d->size1;j++)
			{
				tl=(long)j+nodiagshift;
				if((tl>0)&&(tl+1<(long)d->size2))
				{
					vvl=MATRIXFF(subrow)(d,j,0,(size_t)tl);
					VECTORFF(minmax)(&vvl.vector,&dmin1,&dmax1);
					dmin=GSL_MIN(dmin,dmin1);
					dmax=GSL_MAX(dmax,dmax1);
					vvl=MATRIXFF(subrow)(d,j,(size_t)tl+1,d->size2-(size_t)tl-1);
					VECTORFF(minmax)(&vvl.vector,&dmin1,&dmax1);
				}
				else if((tl<0)||((size_t)tl>=d->size2))
				{
					vvl=MATRIXFF(row)(d,j);
					VECTORFF(minmax)(&vvl.vector,&dmin1,&dmax1);
				}
				else
				{
					vvl=MATRIXFF(subrow)(d,j,tl?0:1,d->size2-1);
					VECTORFF(minmax)(&vvl.vector,&dmin1,&dmax1);
				}
				dmin=GSL_MIN(dmin,dmin1);
				dmax=GSL_MAX(dmax,dmax1);
			}
		}
		else
			MATRIXFF(minmax)(d,&dmin,&dmax);
		if((!(dmin>=0))||gsl_isnan(dmax)||gsl_isinf(dmax))
			ERRRET("Negative or NAN found in input data. It may invalidate follow up analysis. This may be due to incorrect previous steps.")
		h=pij_llrtopij_a_nullhist((double)dmax,nv,d->size2,n1c,n1d,n2c,n2d);
		if(!h)
			ERRRET("pij_llrtopij_a_nullhist failed.")
		nbin=h[0]->n;
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
			hreal[i]=gsl_histogram_clone(h[0]);
			hc[i]=gsl_histogram_alloc(hreal[i]->n+2);
			ret=ret&&hreal[i]&&hc[i];
		}
		vcount=VECTORGF(alloc)(ng);
		if(!(ret&&vcount))
			ERRRET("Not enough memory.");
	}
	
	{
		VECTORUC	*vb4=VECTORUCF(alloc)(nv);
		if(!vb4)
			ERRRET("Not enough memory.");	
		MATRIXGF(countv_byrow_buffed)(g,vcount,vb4);
		CLEANVECUC(vb4)
	}
	
	
	//Conversion
	vv1=VECTORDF(view_array)(h[0]->range+1,nbin);
	VECTORDF(memcpy)(vwidth,&vv1.vector);
	vv1=VECTORDF(view_array)(h[0]->range,nbin);
	VECTORDF(sub)(vwidth,&vv1.vector);
	
	
	//Conversion
	for(i=2;i<=nv;i++)
	{
		vv1=VECTORDF(view_array)(h[i-2]->range+1,nbin);
		VECTORDF(memcpy)(vwidth,&vv1.vector);
		vv1=VECTORDF(view_array)(h[i-2]->range,nbin);
		VECTORDF(sub)(vwidth,&vv1.vector);
		vv1=VECTORDF(view_array)(h[i-2]->bin,nbin);
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
				if(VECTORGF(get)(vcount,j)==i)
				{
					MATRIXFF(get_row)(&vvb3.vector,d,j);
					VECTORDF(memcpy)(&vvnull.vector,&vv1.vector);
					memcpy(hreal[id]->range,h[i-2]->range,(nbin+1)*sizeof(*hreal[id]->range));
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
	}
	CLEANUP
	return 0;
#undef	CLEANUP
}

/* Functions to convert LLR of specific steps into probabilities.
 * Uses pij_llrtopij_a_convert with different settings of n1d and n2d.
 * Function name suffices indicate which LLR to convert.
 */
static inline int pij_gassist_llrtopij1_a(MATRIXF* d,const MATRIXG* g,size_t nv,char nodiag,long nodiagshift)
{
	LOG(9,"Converting LLR to probabilities for step 1 on per A basis.")
	return pij_gassist_llrtopij_a_convert_self(d,g,nv,1,1,1,g->size2-2,nodiag,nodiagshift);
}

static inline int pij_gassist_llrtopij2_a(MATRIXF* d,const MATRIXG* g,size_t nv,char nodiag,long nodiagshift)
{
	LOG(9,"Converting LLR to probabilities for step 2 on per A basis.")
	return pij_gassist_llrtopij_a_convert_self(d,g,nv,1,1,1,g->size2-2,nodiag,nodiagshift);
}

static inline int pij_gassist_llrtopij3_a(MATRIXF* d,const MATRIXG* g,size_t nv,char nodiag,long nodiagshift)
{
	LOG(9,"Converting LLR to probabilities for step 3 on per A basis.")
	if(pij_gassist_llrtopij_a_convert_self(d,g,nv,1,1,1,g->size2-3,nodiag,nodiagshift))
		return 1;
	MATRIXFF(scale)(d,-1);
	MATRIXFF(add_constant)(d,1);
	return 0;
}

static inline int pij_gassist_llrtopij4_a(MATRIXF* d,const MATRIXG* g,size_t nv,char nodiag,long nodiagshift)
{
	LOG(9,"Converting LLR to probabilities for step 4 on per A basis.")
	return pij_gassist_llrtopij_a_convert_self(d,g,nv,1,2,1,g->size2-3,nodiag,nodiagshift);
}

static inline int pij_gassist_llrtopij5_a(MATRIXF* d,const MATRIXG* g,size_t nv,char nodiag,long nodiagshift)
{
	LOG(9,"Converting LLR to probabilities for step 5 on per A basis.")
	if(pij_gassist_llrtopij_a_convert_self(d,g,nv,0,1,1,g->size2-3,nodiag,nodiagshift))
		return 1;
	return 0;
}


int pij_gassist_llrtopijs_a(const MATRIXG* g,VECTORF* p1,MATRIXF* p2,MATRIXF* p3,MATRIXF* p4,MATRIXF* p5,size_t nv,char nodiag,long nodiagshift)
{
	int	ret=0,ret2=0;
	ret=ret||(ret2=pij_gassist_llrtopij2_a(p2,g,nv,0,0));
	if(ret2)
		LOG(1,"Failed to log likelihood ratios to probabilities in step 2.")
	//For p1, if nodiag, copy p2 data, otherwise set all to 1.
	if(nodiag)
	{
		VECTORFF(view)	vv;
		vv=MATRIXFF(superdiagonal)(p2,(size_t)nodiagshift);
		ret=(ret2=VECTORFF(memcpy)(p1,&vv.vector));
	}
	else
		ret=(ret2=pij_gassist_llrtopij1_1(p1));
	if(ret2)
		LOG(1,"Failed to log likelihood ratios to probabilities in step 1.")
	ret=ret||(ret2=pij_gassist_llrtopij3_a(p3,g,nv,nodiag,nodiagshift));
	if(ret2)
		LOG(1,"Failed to log likelihood ratios to probabilities in step 3.")
	ret=ret||(ret2=pij_gassist_llrtopij4_a(p4,g,nv,nodiag,nodiagshift));
	if(ret2)
		LOG(1,"Failed to log likelihood ratios to probabilities in step 4.")
	ret=ret||(ret2=pij_gassist_llrtopij5_a(p5,g,nv,nodiag,nodiagshift));
	if(ret2)
		LOG(1,"Failed to log likelihood ratios to probabilities in step 5.")
	return ret;
}











