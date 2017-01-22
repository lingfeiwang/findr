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
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "gsl/math.h"
#include "gsl/blas.h"
#include "gsl/sort.h"
#include "logger.h"
#include "random.h"
#include "macros.h"
#include "general_alg.h"
#include "histogram.h"

/*********************************************************
 * Histogram bin range constructions
 *********************************************************/

// Uniform bins
double* histogram_uniformbins(double dmin,double dmax,size_t n)
{
	double* ans;
	size_t	i;
	double	step;
	
	MALLOCSIZE(ans,n+1);
	if(!ans)
	{
		LOG(1,"Not enough memory.")
		return 0;
	}
	step=(dmax-dmin)/(double)n;
	ans[0]=dmin;
	for(i=1;i<=n;i++)
		ans[i]=ans[i-1]+step;
	histogram_fluc_binrange(ans,n,1E-5);
	return ans;
}

// Bins of (near) equal count
int histogram_equalbins_fromnullpdfs(size_t n,double* binrange,void (*pdfs)(const VECTORD*,VECTORD*,const void*),const void* param)
{
#define	CLEANUP	CLEANVECD(vc)CLEANVECD(vval)AUTOFREEVEC(vstep)

	size_t	i;
	size_t	nb;
	const size_t	nbextend=1000;		//Interpolation count within each bin
	const double	eps=1E-5;			//Relative precision
	VECTORD			*vc,*vval;
	AUTOALLOCVECD(vstep,nbextend,10000)
	double	step,dnext,dlast,dnow,diffmax,t1,t2;
	size_t	inext,iaim;
	VECTORDF(view)	vv1,vv2;
	MATRIXDF(view)	mv1;

	assert(n);
	//Extended bin range
	nb=nbextend*n;
	vc=VECTORDF(alloc)(nb);
	vval=VECTORDF(alloc)(nb);
	
	if(!(vc&&vval&&vstep))
		ERRRETV(0,"Not enough memory.")
	
	//Initial uniform bin range
	step=(binrange[n]-binrange[0])/(double)n;
	for(i=1;i<n;i++)
		binrange[i]=binrange[i-1]+step;
	//Initialize step vector for interpolated locations
	VECTORDF(set)(vstep,0,0.5);
	for(i=1;i<nbextend;i++)
		VECTORDF(set)(vstep,i,VECTORDF(get)(vstep,i-1)+1);
	VECTORDF(scale)(vstep,1./(double)nbextend);
	
	diffmax=eps*2;
	while(diffmax>eps)
	{
		//Calculate extended central locations
		//First fill start locations
		for(i=0;i<n;i++)
		{
			vv1=VECTORDF(subvector)(vc,i*nbextend,nbextend);
			VECTORDF(set_all)(&vv1.vector,binrange[i]);
		}
		//Then calculate bin widths
		vv1=VECTORDF(view_array)(binrange+1,n);
		vv2=VECTORDF(subvector)(vval,0,n);
		VECTORDF(memcpy)(&vv2.vector,&vv1.vector);
		vv1=VECTORDF(view_array)(binrange,n);
		VECTORDF(sub)(&vv2.vector,&vv1.vector);
		//Finally rank 1 update
		mv1=MATRIXDF(view_vector)(vc,n,nbextend);
		gsl_blas_dger(1,&vv2.vector,vstep,&mv1.matrix);

		//Calculate bin pdf values
		pdfs(vc,vval,param);
		//Calculate bin width in vc
		for(i=0;i<n;i++)
		{
			vv1=VECTORDF(subvector)(vc,i*nbextend,nbextend);
			VECTORDF(set_all)(&vv1.vector,binrange[i+1]-binrange[i]);
		}
		VECTORDF(scale)(vc,1./(double)nbextend);
		//Convert to histogram count
		VECTORDF(mul)(vval,vc);
		
		//Convert to cumulative sum and scale to n at maximum
		for(i=1;i<nb;i++)
			VECTORDF(ptr)(vval,i)[0]+=VECTORDF(get)(vval,i-1);
		VECTORDF(scale)(vval,(double)n/VECTORDF(get)(vval,nb-1));
		//Store original bin location
		vv1=VECTORDF(subvector_with_stride)(vc,0,nbextend,n);
		vv2=VECTORDF(view_array)(binrange,n);
		VECTORDF(memcpy)(&vv1.vector,&vv2.vector);
		
		//Reconstruct bins
		for(i=0,dlast=0,dnext=1,inext=1,diffmax=0;inext<n;i++)
		{
			dnow=VECTORDF(get)(vval,i);
			if(dnow>dnext)
			{
				step=VECTORDF(get)(vc,(i/nbextend)*nbextend+1)/(dnow-dlast);
				iaim=(size_t)floor(dnow);
				if(iaim>=n)
					iaim=n-1;
				while(inext<=iaim)
				{
					//Part 1: histogram bin start location
					t1=VECTORDF(get)(vc,(i/nbextend)*nbextend)
						//Part 2: interpolation bin start location
						+VECTORDF(get)(vc,(i/nbextend)*nbextend+1)*(double)(i%nbextend)
						//Part 3: increment inside interpolation bin
						+((double)inext-dlast)*step;
					//Store maximum relative variation in bin width
					t2=fabs(t1-binrange[inext])/(binrange[inext+1]-binrange[inext]);
					if(t2>diffmax)
						diffmax=t2;
					binrange[inext++]=t1;
				}
				dnext=(double)inext;
			}
			dlast=dnow;
		}
	}
	CLEANUP
	return 0;
#undef	CLEANUP
}

/* Bins of with smooth transition from bins of equal count on left side
 * to bins of equal width on right side.
 */
size_t histogram_unequalbins_param_count(size_t n)
{
	size_t	t;
	t=(size_t)floor(exp(log((double)n)/2.5));
	if(t>100)
		t=100;
	return t;
}

int histogram_unequalbins_fromnullpdfs(size_t n,double* binrange,void (*pdfs)(const VECTORD*,VECTORD*,const void*),const void* param)
{
#define	CLEANUP
	//First produce bins of equal density
	if(histogram_equalbins_fromnullpdfs(n,binrange,pdfs,param))
		ERRRET("histogram_equalbins_fromnullpdfs failed.")
	//Then convert to unequal bins
	histogram_unequalbins_fromequalbins(n,binrange);
	return 0;
#undef	CLEANUP
}

void histogram_unequalbins_fromnullcdf(size_t n,double* binrange,double (*cdf)(double,const void*),const void* param)
{
	//Obtain quantiles
	histogram_equalbins_fromnullcdf(n,binrange[0],binrange[n],cdf,param,(binrange[n]-binrange[0])*1E-6,binrange+1);
	//Then convert to unequal bins
	histogram_unequalbins_fromequalbins(n,binrange);
}

/* Bins of with sharp transition from bins of equal count on left side
 * to bins of equal width on right side.
 */
gsl_histogram* histogram_unequalbins_exact(const VECTORF* d,const void* param)
{
#define	CLEANUP	CLEANHIST(h)CLEANMEM(binr)
	int	ret;
	const struct histogram_unequalbins_exact_param *p=param;
	gsl_histogram	*h=0;
	double			*binr=0;
	size_t			nbin,nbin1;
	double			nbinsplit;
	
	ret=histogram_unequalbins_param_sizing(d,p->ebwrate,&nbinsplit,&nbin1);
	if(ret)
		ERRRETV(0,"Failed to auto-discover histogram paramters: not enough data.")
	binr=histogram_unequalbins_exact_vf(d,nbinsplit,&nbin1,&nbin);
	if(!binr)
		ERRRETV(0,"Bin range construction failed for "PRINTFSIZET" bins of equal count for %F quantile.",nbin1,nbinsplit)
	LOG(9,"Histogram range constructed: "PRINTFSIZET" bins of equal count for %F quantile, "PRINTFSIZET" bins of equal width for the rest.",nbin1,nbinsplit,nbin-nbin1)
	h=gsl_histogram_alloc(nbin);
	if(!h)
		ERRRETV(0,"Not enough memory.")
	gsl_histogram_set_ranges(h,binr,h->n+1);
	CLEANMEM(binr)
	return h;
#undef	CLEANUP
}

gsl_histogram* histogram_unequalbins_exp(const VECTORF* d,const void* param)
{
#define	CLEANUP	CLEANHIST(h)CLEANMEM(binr)
	int	ret;
	const struct histogram_unequalbins_exp_param *p=param;
	gsl_histogram	*h=0;
	double			*binr=0;
	size_t			nbin,nbin1;
	double			nbinsplit;
	
	ret=histogram_unequalbins_param_sizing(d,p->ebwrate,&nbinsplit,&nbin1);
	if(ret)
		ERRRETV(0,"Failed to auto-discover histogram paramters: not enough data.")
	binr=histogram_unequalbins_exp_vf(d,p->nsmnv,nbinsplit,nbin1,&nbin);
	if(!binr)
		ERRRETV(0,"Bin range construction failed for "PRINTFSIZET" bins of equal count for %F quantile.",nbin1,nbinsplit)
	LOG(9,"Histogram range constructed: "PRINTFSIZET" bins of equal count for %F quantile, "PRINTFSIZET" bins of equal width for the rest.",nbin1,nbinsplit,nbin-nbin1)
	h=gsl_histogram_alloc(nbin);
	if(!h)
		ERRRETV(0,"Not enough memory.")
	gsl_histogram_set_ranges(h,binr,h->n+1);
	CLEANMEM(binr)
	return h;
#undef	CLEANUP
}

int histogram_unequalbins_param_sizing(const VECTORF* d,double ebwrate,double * restrict nbinsplit, size_t* restrict nbin1)
{
	double rate=ebwrate*5;
	int		nbin;
	*nbinsplit=1-rate/sqrt((double)d->size);
	*nbinsplit=GSL_MIN(*nbinsplit,0.9995);
	nbin=(int)(sqrt(sqrt((double)d->size))-rate);
	nbin=GSL_MIN(nbin,100);
	*nbin1=(size_t)nbin;
	return (nbin<4)||(*nbinsplit<0);
}

double* histogram_unequalbins_exp_vf(const VECTORF* d,size_t nsmnv,double nbinsplit,size_t nbin1,size_t* restrict nbin)
{
#define	CLEANUP	CLEANHIST(h1)

	size_t	i;
	size_t	nb,nbe;
	const size_t	nbr=100;
	const size_t	nbextend=50;
	double*	restrict ans;
	gsl_histogram	*h1;
	size_t	head;
	double	tunit,tremain,tadd,tloc;
	double	tdensity,tlack,twidth;
	VECTORDF(view)	vv;

	assert(d&&nsmnv&&(nbinsplit>0)&&(nbinsplit<1)&&nbin1&&nbin);
	*nbin=0;
	//Extended bin range
	nb=nbr*nbin1+1;
	nbe=nb+nbextend*nbin1;
	h1=gsl_histogram_alloc(nbe);
	MALLOCSIZE(ans,nbin1+1);
	if(!(h1&&ans))
		ERRRETV(0,"Not enough memory.")
	
	//Construct initial bin range
	for(i=0;i<nb;i++)
		h1->range[i]=(double)i;
	vv=VECTORDF(view_array)(h1->range+1,nb-1);
	VECTORDF(scale)(&vv.vector,-nbinsplit/(double)nb);
	VECTORDF(add_constant)(&vv.vector,1);
	for(i=1;i<nb;i++)
		h1->range[i]=log(h1->range[i]);
	VECTORDF(scale)(&vv.vector,-1/(double)nsmnv);
	//Construct initial extended bin range
	twidth=h1->range[nb-1]-h1->range[nb-2];
	for(i=nb;i<h1->n;i++)
		h1->range[i]=h1->range[i-1]+twidth;
	h1->range[h1->n]=INFINITY;
	
	//Sow histogram
	memset(h1->bin,0,h1->n*sizeof(*h1->bin));
	for(i=0;i<d->size;i++)
		gsl_histogram_increment(h1,VECTORFF(get)(d,i));
	
	//Form new bin range, assuming uniform within bin
	tremain=0;
	head=1;
	ans[0]=0;
	tunit=floor((double)d->size*nbinsplit/(double)nbin1);
	for(i=0;(i<h1->n)&&(head<=nbin1);i++)
	{
		tloc=0;
		tadd=gsl_histogram_get(h1,i);
		tremain+=tadd;
		//Continue accumulate
		if(tremain<tunit)
			continue;
		//Form new bin
		tremain-=tadd;
		twidth=h1->range[i+1]-h1->range[i];
		tdensity=tadd/twidth;
		while((tremain+tadd>=tunit)&&(head<=nbin1))
		{
			//Update accumulation
			tlack=(tunit-tremain);
			tloc+=tlack/tdensity;
			tadd-=tlack;
			tremain=0;
			//Construct new bin
			ans[head++]=h1->range[i]+tloc;
		}
		tremain+=tadd;
	}
	CLEANHIST(h1)
	
	//Calculate uniform bin part and arrange for output
	{
		double	vmax;
		size_t	nbin2;
		
		vmax=(double)VECTORFF(max)(d);
		twidth=ans[nbin1]-ans[nbin1-1];
		nbin2=(size_t)ceil((vmax-ans[nbin1])/twidth);
		
		//Arrange for output
		*nbin=nbin1+nbin2;
		ans=realloc(ans,(*nbin+1)*sizeof(*ans));
		if(!ans)
			ERRRETV(0,"Not enough memory.")
		for(i=nbin1+1;i<=*nbin;i++)
			ans[i]=ans[i-1]+twidth;
	}
	
	histogram_fluc_binrange(ans,*nbin,1E-5);
	return ans;
#undef	CLEANUP
}


double* histogram_unequalbins_exact_vf(const VECTORF* d,double nbinsplit,size_t* restrict nbin1,size_t* restrict nbin)
{
	size_t		i,nbin2;
	AUTOALLOCVECF(vs,d->size,10000)
	double*	restrict	ans;	//Return data
	double		step;
	size_t		l1,l2;

	*nbin=0;
	if(!vs)
	{
		LOG(1,"Not enough memory.")
		return 0;
	}

	VECTORFF(memcpy)(vs,d);
	CONCATENATE2(gsl_sort_vector,FTYPE_SUF)(vs);
	//Determine total number of bins
	l1=(size_t)floor((double)d->size*nbinsplit);
	step=(double)l1/(double)*nbin1;
	l2=(size_t)floor((double)l1-step);
	nbin2=(size_t)ceil((VECTORFF(get)(vs,vs->size-1)-VECTORFF(get)(vs,l1))/(VECTORFF(get)(vs,l1)-VECTORFF(get)(vs,l2)));
	MALLOCSIZE(ans,*nbin1+nbin2+1);
	if(!ans)
	{
		LOG(1,"Not enough memory.")
		AUTOFREEVEC(vs)
		return 0;
	}
	//Setting bin values
	for(i=0;i<*nbin1;i++)
		ans[i]=VECTORFF(get)(vs,(size_t)floor(step*(double)i));
	//Remove duplicates
	*nbin1=remove_sorted_duplicates(ans,*nbin1);
	*nbin=*nbin1+nbin2;
		
	step=(VECTORFF(get)(vs,vs->size-1)-VECTORFF(get)(vs,l1))/(double)nbin2;
	for(i=0;i<nbin2;i++)
		ans[i+*nbin1]=VECTORFF(get)(vs,l1)+(double)i*step;
	ans[*nbin1+nbin2]=VECTORFF(get)(vs,vs->size-1);
	
	histogram_fluc_binrange(ans,*nbin,1E-5);
	AUTOFREEVEC(vs)
	return ans;
}



/*********************************************************
 * Histogram bin range operations
 *********************************************************/
 
void histogram_fluc_binrange(double* restrict range,size_t n,double diff)
{
	size_t i;
	//Weakly fluctuate bin edges
	for(i=n-1;i;i--)
		range[i]+=GSL_MIN(range[i]-range[i-1],range[i+1]-range[i])*(random_uniform()-0.5)*diff;
	//Extend boundary bin edges
	range[0]-=(range[1]-range[0])*diff;
	range[n]+=(range[n]-range[n-1])*diff;
}

void histogram_finer_range(const VECTORD* hold,VECTORD* hnew,size_t nsp)
{
	size_t	nbin=hold->size-1;
	size_t	i;
	VECTORDF(view) vv1,vv2,vv3;
	
	assert(nsp>1);
	assert(hnew->size==nbin*nsp+1);

	//For existing ones
	vv1=VECTORDF(subvector_with_stride)(hnew,0,nsp,nbin+1);
	VECTORDF(memcpy)(&vv1.vector,hold);
	
	//For new ones
	//Obtain difference
	vv1=VECTORDF(subvector_with_stride)(hnew,0,nsp,nbin);
	vv2=VECTORDF(subvector_with_stride)(hnew,nsp,nsp,nbin);
	vv3=VECTORDF(subvector_with_stride)(hnew,nsp-1,nsp,nbin);
	VECTORDF(memcpy)(&vv3.vector,&vv2.vector);
	VECTORDF(sub)(&vv3.vector,&vv1.vector);
	VECTORDF(scale)(&vv3.vector,1/(double)nsp);
	//Cumulative add up difference
	for(i=1;i<nsp-1;i++)
	{
		vv1=VECTORDF(subvector_with_stride)(hnew,i-1,nsp,nbin);
		vv2=VECTORDF(subvector_with_stride)(hnew,i,nsp,nbin);
		VECTORDF(memcpy)(&vv2.vector,&vv1.vector);
		VECTORDF(add)(&vv2.vector,&vv3.vector);
	}
	vv2=VECTORDF(subvector_with_stride)(hnew,nsp,nsp,nbin);
	VECTORDF(sub)(&vv3.vector,&vv2.vector);
	VECTORDF(scale)(&vv3.vector,-1);
}

void histogram_finer_central(const VECTORD* hold,VECTORD* hnew,size_t nsp)
{
	size_t	nbin=hold->size-1;
	size_t	i;
	VECTORDF(view) vv1,vv2,vv3;
	VECTORDF(const_view)	vvc1=VECTORDF(const_subvector)(hold,0,nbin);
	VECTORDF(const_view)	vvc2=VECTORDF(const_subvector)(hold,1,nbin);
	
	assert(nsp>1);
	assert(hnew->size==nbin*nsp);
	
	//Obtain difference
	vv3=VECTORDF(subvector_with_stride)(hnew,nsp-1,nsp,nbin);
	VECTORDF(memcpy)(&vv3.vector,&vvc2.vector);
	VECTORDF(sub)(&vv3.vector,&vvc1.vector);
	VECTORDF(scale)(&vv3.vector,1/(double)nsp);
	//Obtain first central bins
	vv1=VECTORDF(subvector_with_stride)(hnew,0,nsp,nbin);
	VECTORDF(memcpy)(&vv1.vector,&vvc1.vector);
	gsl_blas_daxpy(0.5,&vv3.vector,&vv1.vector);
	//Obtain other central bins
	for(i=1;i<nsp-1;i++)
	{
		vv1=VECTORDF(subvector_with_stride)(hnew,i-1,nsp,nbin);
		vv2=VECTORDF(subvector_with_stride)(hnew,i,nsp,nbin);
		VECTORDF(memcpy)(&vv2.vector,&vv1.vector);
		VECTORDF(add)(&vv2.vector,&vv3.vector);
	}
	//Obtain final central bin
	VECTORDF(sub)(&vv3.vector,&vvc2.vector);
	VECTORDF(scale)(&vv3.vector,-1);
}


void histogram_unequalbins_fromequalbins(size_t n,double* binrange)
{
	size_t	i;
	double	twidth,width,t1,t2;
	VECTORDF(view)	vv;
	
	//Convert to bin widths
	twidth=binrange[n]-binrange[0];
	width=log(twidth/(double)n);
	for(i=n;i;i--)
		binrange[i]-=binrange[i-1];
	//bin width[i]=(old width)^((n-2-i)/(n-2))(mean width)^(i/(n-2))
	t2=((double)(n-1));
	for(i=2;i<n;i++)
		binrange[i]=exp((log(binrange[i])*(t2-(double)(i-1))+width*((double)(i-1)))/t2);
	binrange[n]=twidth/(double)n;
	vv=VECTORDF(view_array)(binrange+1,n);
	t1=twidth/gsl_blas_dasum(&vv.vector);
	VECTORDF(scale)(&vv.vector,t1);
	for(i=1;i<n;i++)
		binrange[i]+=binrange[i-1];
	binrange[n]=binrange[0]+twidth;
}


























