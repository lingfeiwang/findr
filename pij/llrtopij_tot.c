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
#include <math.h>
#include <float.h>
#include <string.h>
#include "../base/gsl/math.h"
#include "../base/gsl/histogram.h"
#include "../base/gsl/blas.h"
#include "../base/logger.h"
#include "../base/macros.h"
#include "../base/data_process.h"
#include "../base/threading.h"
#include "../base/histogram.h"
#include "nullhist.h"
#include "llrtopij_tot.h"
#include "llrtopij.h"
#pragma GCC diagnostic ignored "-Wunused-parameter"

int pij_llrtopij_tot_convert(const VECTORF* d,const VECTORF* dconv,VECTORF* ans,gsl_histogram* (*fh)(const VECTORF*,const void*),int (*fnh)(const void*,gsl_histogram*),const void* ph,const void* pnh,double* nullratio)
{
#define	CLEANUP	CLEANHIST(hist)CLEANHIST(histc)AUTOFREE(histn)AUTOFREEVEC(vbuff)

	gsl_histogram	*hist,*histc;
// 	double*	histn;	Normal histogram of null density
	double*	tp;
	VECTORDF(view)	vv;
	size_t			i;
	int				ret;

	hist=fh(d,ph);
	if(!hist)
	{
		LOG(1,"Histogram bin range construction failed.")
		return 1;
	}
	AUTOALLOC(double,histn,hist->n,1000)
	histc=gsl_histogram_alloc(hist->n+2);
	AUTOALLOCVECD(vbuff,hist->n,1000)
	if(!(histn&&histc&&vbuff))
		ERRRET("Not enough memory.")
	//Real histogram	
	for(i=0;i<d->size;i++)
			gsl_histogram_increment(hist,(double)VECTORFF(get)(d,i));

	//Null histogram
	tp=histn;
	histn=hist->bin;
	hist->bin=tp;
	ret=fnh(pnh,hist);
	tp=histn;
	histn=hist->bin;
	hist->bin=tp;
	if(ret)
		ERRRET("Custom null histogram function failed.")
	
	//Calculate ratio of true distribution
	{
		VECTORDF(view)	vvnull=VECTORDF(view_array)(histn,hist->n);
		//convert to density histogram
		vv=VECTORDF(view_array)(hist->range+1,hist->n);
		VECTORDF(memcpy)(vbuff,&vv.vector);
		vv=VECTORDF(view_array)(hist->range,hist->n);
		VECTORDF(sub)(vbuff,&vv.vector);
		VECTORDF(div)(&vvnull.vector,vbuff);
		vv=VECTORDF(view_array)(hist->bin,hist->n);
		VECTORDF(scale)(&vv.vector,1./(double)d->size);
		VECTORDF(div)(&vv.vector,vbuff);
		if(pij_llrtopij_convert_histograms(hist,&vvnull.vector,histc))
			ERRRET("Failed in pij_llrtopij_convert_histograms.")
	}
	//Convert likelihoods to probabilities
	pij_llrtopij_histogram_interpolate_linear(histc,dconv,ans);
	
	if(nullratio)
	{
		*nullratio=1-BLASF(asum)(ans)/(FTYPE)ans->size;
		LOG(9,"Null distribution ratio: %G",*nullratio)
	}
	CLEANUP
	return 0;
#undef CLEANUP
}





















