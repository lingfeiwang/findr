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
#include "../llrtopij.h"
#include "llrtopij.h"



/* Always return probability of step 1 is 1. This is useful when best eQTL are already selected in advance.
 */
static inline int pij_gassist_llrtopij1_1(VECTORF* p1)
{
	LOG(9,"Converting LLR to probabilities for step 1. Filling with 1.")
	VECTORFF(set_all)(p1,1);
	return 0;
}

/* Convert real log likelihood ratios into probability functions.
 * This function converts every A in hypothesis (E->A->B) separately.
 * Suppose there are ng (E,A) pairs and nt Bs, this function converts ng times,
 * each for one (E,A) pair but all Bs.
 * d:		(ng,nt)	Input log likelihood ratios for construction of 
 * 			histograms and calculation of probability of true hypothesis.
 * g:		(ng,ns)	Original genotype matrix, used for analytical calculation
 * 			of null distribution. Every element=0,1,...,nv-1.
 * h:		[nv-1]. Null histogram of the specific test.
 * 			Output of pij_nullhist.
 * nv:		Maximum number of values each g may take.
 * nodiag:	If diagonal elements of d should be removed in construction of real
 * 			histogram. This should be set to true (!=0) when t is identical with
 * 			the top rows of t2 (in calculation of llr).
 * Return:	0 if success.
 */
static int pij_gassist_llrtopij_convert_self(MATRIXF* d,const MATRIXG* g,const gsl_histogram * const * h, size_t nv,char nodiag,long nodiagshift)
{
#define	CLEANUP	CLEANVECG(vcount)CLEANAMHIST(hreal,nth)CLEANAMHIST(hc,nth)\
				CLEANMATD(mb1)CLEANMATD(mb2)CLEANMATD(mnull)CLEANMATF(mb3)CLEANVECD(vwidth)

	VECTORG		*vcount;
	size_t		ng=g->size1;
	size_t		i,nbin;
	//gsl_histogram **hreal,**hc;
	MATRIXD		*mb1,*mb2,*mnull;
	MATRIXF		*mb3;
	VECTORD		*vwidth;
	VECTORDF(view)	vv1;
	size_t		nth;

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
	nbin=h[0]->n;
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

int pij_gassist_llrtopijs(const MATRIXG* g,VECTORF* p1,MATRIXF* p2,MATRIXF* p3,MATRIXF* p4,MATRIXF* p5,size_t nv,const gsl_histogram * const * h[4],char nodiag,long nodiagshift)
{
	int	ret=0,ret2=0;
	if(g->size2<=3)
	{
		LOG(0,"Needs at least 4 samples to compute probabilities.")
		return 1;
	}
	ret=ret||(ret2=pij_gassist_llrtopij_convert_self(p2,g,h[0],nv,nodiag,nodiagshift));
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
	ret=ret||(ret2=pij_gassist_llrtopij_convert_self(p3,g,h[1],nv,nodiag,nodiagshift));
	if(ret2)
		LOG(1,"Failed to log likelihood ratios to probabilities in step 3.")
	MATRIXFF(scale)(p3,-1);
	MATRIXFF(add_constant)(p3,1);
	ret=ret||(ret2=pij_gassist_llrtopij_convert_self(p4,g,h[2],nv,nodiag,nodiagshift));
	if(ret2)
		LOG(1,"Failed to log likelihood ratios to probabilities in step 4.")
	ret=ret||(ret2=pij_gassist_llrtopij_convert_self(p5,g,h[3],nv,nodiag,nodiagshift));
	if(ret2)
		LOG(1,"Failed to log likelihood ratios to probabilities in step 5.")
	return ret;
}












