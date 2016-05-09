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
#include "../../base/macros.h"
#include "../../base/data_process.h"
#include "../../base/histogram.h"
#include "../nulldist.h"
#include "../llrtopij.h"
#include "llrtopij_a.h"
#include "llrtopij.h"


int pij_gassist_llrtopij_a_convert(const MATRIXF* d,const MATRIXF* dconv,MATRIXF* ans,const MATRIXG* g,size_t nv,long n1d,long n2d,char nodiag)
{
#define	CLEANUP	CLEANVECG(vcount)CLEANHIST(hreal)CLEANHIST(hc)\
				if(h){for(i=0;i<nv-1;i++)CLEANHIST(h[i])free(h);h=0;}\
				CLEANVECD(vb1)CLEANVECD(vb2)
	VECTORG		*vcount;
	size_t		ng=g->size1,ns=g->size2;
	size_t		i,j,k,nbin;
	gsl_histogram **h,*hreal,*hc;
	VECTORD		*vb1,*vb2;
	
	h=0;
	hreal=hc=0;
	vcount=0;
	vb1=0;
	vb2=0;
	//Validity checks
	assert((dconv->size1==ng)&&(ans->size1==ng)&&(d->size1==ng)&&(dconv->size2==ans->size2));

	//Construct null density histograms
	{
		FTYPE		dmin,dmax;
		MATRIXFF(minmax)(d,&dmin,&dmax);
		if(!(dmin>=0))
			ERRRET("Negative or NAN found in input data. It may invalidate follow up analysis. This may be due to incorrect previous steps.")
		h=pij_llrtopij_a_nullhist((double)dmax,nv,ns,d->size2,n1d,n2d);
		if(!h)
			ERRRET("pij_llrtopij_a_nullhist failed.")
		nbin=h[0]->n;
	}
	
	//Prepare for real histogram
	hreal=gsl_histogram_clone(h[0]);
	hc=gsl_histogram_alloc(hreal->n+2);
	vcount=VECTORGF(alloc)(ng);
	if(!(hreal&&hc&&vcount))
		ERRRET("Not enough memory.");
	if(pij_llrtopij_convert_histograms_make_buffs(nbin,&vb1,&vb2))
		ERRRET("pij_llrtopij_convert_histograms_make_buffs failed.")
	{
		VECTORUC	*vb3=VECTORUCF(alloc)(nv);
		if(!vb3)
			ERRRET("Not enough memory.");	
		MATRIXGF(countv_byrow_buffed)(g,vcount,vb3);
		CLEANVECUC(vb3)
	}
	
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
		vvreal=VECTORDF(view_array)(hreal->bin,nbin);
		for(i=2;i<=nv;i++)
		{
			vv1=VECTORDF(view_array)(h[i-2]->range+1,nbin);
			VECTORDF(memcpy)(vwidth,&vv1.vector);
			vv1=VECTORDF(view_array)(h[i-2]->range,nbin);
			VECTORDF(sub)(vwidth,&vv1.vector);
			for(j=0;j<ng;j++)
				if(VECTORGF(get)(vcount,j)==i)
				{
				
					VECTORFF(const_view)	vvd=MATRIXFF(const_row)(dconv,j);
					vv1=VECTORDF(view_array)(h[i-2]->bin,nbin);
					VECTORDF(memcpy)(vnull,&vv1.vector);
					memcpy(hreal->range,h[i-2]->range,(nbin+1)*sizeof(double));
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
						VECTORDF(scale)(&vvreal.vector,1./(double)d->size2);
					}
					//Convert to density histogram
					VECTORDF(div)(&vvreal.vector,vwidth);
					//Convert to probability central histogram
					pij_llrtopij_convert_histograms_buffed(hreal,vnull,hc,vb1,vb2);
					//Convert likelihoods to probabilities
					vva=MATRIXFF(row)(ans,j);
					pij_llrtopij_histogram_interpolate_linear(hc,&vvd.vector,&vva.vector);
				}
		}
		CLEANVECD(vwidth)CLEANVECD(vnull)
	}
	CLEANUP
	return 0;
#undef	CLEANUP
}

int pij_gassist_llrtopijs_a(const MATRIXG* g,const VECTORF* llr1,const MATRIXF* llr2b,const MATRIXF* llr2c,const MATRIXF* llr3,VECTORF* p1,MATRIXF* p2b,MATRIXF* p2c,MATRIXF* p3,size_t nv,char nodiag)
{
	int	ret=0,ret2;
	ret2=pij_gassist_llrtopij1_1(p1);
	ret=ret||ret2;
	if(ret2)
		LOG(1,"Failed to log likelihood ratios to probabilities in step 1.")
	ret2=pij_gassist_llrtopij2b_a(llr2b,llr2b,p2b,g,nv,nodiag);
	ret=ret||ret2;
	if(ret2)
		LOG(1,"Failed to log likelihood ratios to probabilities in step 2b.")
	ret2=pij_gassist_llrtopij2c_a(llr2c,llr2c,p2c,g,nv,nodiag);
	ret=ret||ret2;
	if(ret2)
		LOG(1,"Failed to log likelihood ratios to probabilities in step 2c.")
	ret2=pij_gassist_llrtopij3_a(llr3,llr3,p3,g,nv,nodiag);
	ret=ret||ret2;
	if(ret2)
		LOG(1,"Failed to log likelihood ratios to probabilities in step 3.")
	return ret;
}



















