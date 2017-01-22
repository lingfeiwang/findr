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
#include <assert.h>
#include <sys/time.h>
#include "../base/gsl/blas.h"
#include "../base/logger.h"
#include "../base/macros.h"
#include "../base/histogram.h"
#include "../base/threading.h"
#include "nullhist.h"
#include "nulldist.h"

/* Null histogram generation of sample-model method on single thread.
 * g,t,t2,sampler,modeler,ps:	See pij_nullhist_sample_model_param.
 * pmc:		Model container object. See definition in nullmodel.h.
 * Return:	0 if success.
 */
static int pij_nullhist_sample_model_single(const void* data,const struct pij_nullsampler* sampler,const struct pij_nullmodeler* modeler,const void* ps,void* pmc)
{
#define	CLEANUP	if(ptrs){sampler->close(ptrs);ptrs=0;}\
				if(pm){modeler->close(pm);pm=0;}
	
	struct timeval	tv;			
	void	*ptrs,*pm;
	size_t	n1,n2;
	size_t	i;
	
	assert(sampler&&modeler);
	ptrs=pm=0;
	//Initialize sampler with timely random seed.
	gettimeofday(&tv,NULL);
	ptrs=sampler->init(data,ps,(long unsigned int)tv.tv_usec);
	if(!ptrs)
		ERRRET("Sampler construction failed.")
	//Obtain dimensions of random samples.
	sampler->dimensions(ptrs,&n1,&n2);
	assert(n1&&n2);

	AUTOALLOCVECF(vd,n1,30000)	
	if(!vd)
		ERRRET("Not enough memory.")
	//Initialize modeler object given model container object.
	pm=modeler->init(pmc,data,n1,n2);
	if(!pm)
	{
		AUTOFREEVEC(vd)
		ERRRET("Modeler construction failed.")
	}
	//Sample
	for(i=0;i<n2;i++)
	{
		sampler->sample(ptrs,vd);
		modeler->input(pm,vd);
	}
	//Synced updated of model container.
	#pragma omp critical
		modeler->merge(pm,pmc);
	CLEANUP
	AUTOFREEVEC(vd)
	return 0;
#undef	CLEANUP
}

int	pij_nullhist_sample_model(const void* param,gsl_histogram* h)
{
#define	CLEANUP	if(container){\
				p->modeler->close_container(container);\
				container=0;}
	const struct pij_nullhist_sample_model_param *p=param;
	int		ret;
	void*	container;
	

	//Initialize model container
	container=p->modeler->init_container(p->d,p->pm,h);
	if(!container)
		ERRRET("Modeler container construction failed.")

	#pragma omp parallel shared(ret)
	{
		char d[p->dsize];
		//Initializing block selection
		if(p->partition(p->d,d,(size_t)omp_get_thread_num(),(size_t)omp_get_num_threads()))
		{
			int	retth;
		//Calculate null histograms
			retth=pij_nullhist_sample_model_single(d,p->sampler,p->modeler,p->ps,container);
			#pragma omp atomic
				ret+=retth;
		}
	}
	if(ret)
		ERRRET("Failed to construct histogram.")

	//Output from model container
	ret=p->modeler->output(container,h);
	if(ret)
		ERRRET("Failed to output histogram.")
	CLEANUP
	return 0;
#undef	CLEANUP
}

int	pij_nullhist_analytical_pdf(const void* param,gsl_histogram* h)
{
#define	CLEANUP	CLEANVECD(loc)CLEANVECD(val)
	const struct pij_nullhist_analytical_pdf_param	*p=param;

	VECTORD	*loc,*val;
	VECTORDF(view)	vvh=VECTORDF(view_array)(h->bin,h->n);
	size_t		i;
	size_t		nsp;
	double		v;
	
	assert((p->nsplit<10));
	nsp=(size_t)1<<(p->nsplit);
	loc=VECTORDF(alloc)(h->n*nsp);
	val=VECTORDF(alloc)(h->n*nsp);
	if(!(loc&&val))
		ERRRET("Not enough memory.")
	
	//Construct bin ranges
	if(nsp==1)
	{
		VECTORDF(const_view) vv1=VECTORDF(const_view_array)(h->range,h->n);
		VECTORDF(const_view) vv2=VECTORDF(const_view_array)(h->range+1,h->n);
		VECTORDF(memcpy)(loc,&vv1.vector);
		VECTORDF(add)(loc,&vv2.vector);
		VECTORDF(scale)(loc,0.5);
	}
	else
	{
		VECTORDF(const_view) vvc=VECTORDF(const_view_array)(h->range,h->n+1);
		histogram_finer_central(&vvc.vector,loc,nsp);
	}
	//Calculate bin densities
	if(p->func(loc,val,p->param))
	{
		CLEANUP
		return 1;
	}
	
	//Shrink to output
	if(nsp==1)
	{
		VECTORDF(memcpy)(&vvh.vector,val);
	}
	else
	{
		VECTORDF(set_zero)(&vvh.vector);
		for(i=0;i<nsp;i++)
		{
			VECTORDF(const_view)	vvc=VECTORDF(const_subvector_with_stride)(val,i,nsp,h->n);
			VECTORDF(add)(&vvh.vector,&vvc.vector);
		}
	}
	//Convert to bin values from densities
	{
		VECTORDF(const_view)	vv1=VECTORDF(const_view_array)(h->range,h->n);
		VECTORDF(const_view)	vv2=VECTORDF(const_view_array)(h->range+1,h->n);
		VECTORDF(view)			vv3=VECTORDF(subvector)(loc,0,h->n);
		VECTORDF(memcpy)(&vv3.vector,&vv2.vector);
		VECTORDF(sub)(&vv3.vector,&vv1.vector);
		VECTORDF(mul)(&vvh.vector,&vv3.vector);
	}

	//Scale to unity.
	v=gsl_blas_dasum(&vvh.vector);
	VECTORDF(scale)(&vvh.vector,1/v);
	CLEANUP
	return 0;
#undef	CLEANUP
}




























