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
#include <string.h>
#include "../base/gsl/histogram.h"
#include "../base/logger.h"
#include "../base/macros.h"
#include "../base/histogram.h"
#include "nulldist.h"
#include "nullhist.h"



gsl_histogram* pij_nullhist_single(double dmax,size_t nd,size_t n1,size_t n2)
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

gsl_histogram** pij_nullhist(double dmax,size_t nv,size_t nd,long n1c,size_t n1d,long n2c,size_t n2d)
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






































