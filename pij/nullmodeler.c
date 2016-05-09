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
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "../base/gsl/blas.h"
#include "../base/logger.h"
#include "../base/macros.h"
#include "nullmodeler.h"
#pragma GCC diagnostic ignored "-Wunused-parameter"

void* pij_nullmodeler_exp_init_container(const void* d,const void* pm,const gsl_histogram* h)
{
	//Model container same with modeler
	return pij_nullmodeler_exp_init(0,d,0,0);
}

int pij_nullmodeler_exp_output(const void* c,gsl_histogram* h)
{
#define	CLEANUP	AUTOFREEVEC(vwidth)
	const struct pij_nullmodeler_exp_state* p=c;
	double	v;
	size_t	i;
	AUTOALLOCVECD(vwidth,h->n,1000)
	VECTORDF(view)	vv1,vv2;

	if(!vwidth)
		ERRRET("Not enough memory.")

	//Calculate central exponential values as density.
	vv1=VECTORDF(view_array)(h->range,h->n);
	vv2=VECTORDF(view_array)(h->bin,h->n);
	VECTORDF(memcpy)(&vv2.vector,&vv1.vector);
	vv1=VECTORDF(view_array)(h->range+1,h->n);
	VECTORDF(add)(&vv2.vector,&vv1.vector);
	VECTORDF(scale)(&vv2.vector,-(double)(p->n)/(2*p->v));
	for(i=0;i<h->n;i++)
 		h->bin[i]=exp(h->bin[i]);
	//Bin=density*width
	VECTORDF(memcpy)(vwidth,&vv1.vector);
	vv1=VECTORDF(view_array)(h->range,h->n);
	VECTORDF(sub)(vwidth,&vv1.vector);
	VECTORDF(mul)(&vv2.vector,vwidth);
	//Normalize to total 1.
	v=gsl_blas_dasum(&vv2.vector);
	VECTORDF(scale)(&vv2.vector,1/v);
	CLEANUP
	return 0;
#undef	CLEANUP
}

void pij_nullmodeler_exp_close_container(void* c)
{
	pij_nullmodeler_exp_close(c);
}

void* pij_nullmodeler_exp_init(const void* c,const void* d,size_t n1,size_t n2)
{
	struct pij_nullmodeler_exp_state* p;
	p=malloc(sizeof(struct pij_nullmodeler_exp_state));
	if(!p)
		return 0;
	p->n=0;
	p->v=0;
	return p;
}
void pij_nullmodeler_exp_input(void* m,const VECTORF* data)
{
	struct pij_nullmodeler_exp_state* p=m;
 	p->v+=BLASF(asum)(data);
	p->n+=data->size;
}

void pij_nullmodeler_exp_merge(const void* m,void* c)
{
	const struct pij_nullmodeler_exp_state* pm=m;
	struct pij_nullmodeler_exp_state* pc=c;
	pc->n+=pm->n;
	pc->v+=pm->v;
}

void pij_nullmodeler_exp_close(void* m)
{
	if(m)
		free(m);
}








void* pij_nullmodeler_naive_init_container(const void* d,const void* pm,const gsl_histogram* h)
{
#define	CLEANUP	if(p){CLEANHIST(p->h)free(p);p=0;}
	struct pij_nullmodeler_naive_state *p;
	p=calloc(1,sizeof(struct pij_nullmodeler_naive_state));
	if(!p)
		ERRRETV(0,"Not enough memory.")
	p->h=gsl_histogram_alloc(h->n);
	if(!p->h)
		ERRRETV(0,"Not enough memory.")
	//Initialize histogram with same ranges.
	gsl_histogram_set_ranges(p->h,h->range,h->n+1);
	return p;
#undef	CLEANUP
}

int pij_nullmodeler_naive_output(const void* c,gsl_histogram* h)
{
	const struct pij_nullmodeler_naive_state* p=c;
	VECTORDF(view)	vv=VECTORDF(view_array)(h->bin,h->n);
	
	memcpy(h->bin,p->h->bin,h->n*sizeof(double));
	//Scale to unit total probability.
	VECTORDF(scale)(&vv.vector,1/gsl_blas_dasum(&vv.vector));
	return 0;
}

void pij_nullmodeler_naive_close_container(void* c)
{
	struct pij_nullmodeler_naive_state* p=c;
	CLEANHIST(p->h)
	free(p);
}

void* pij_nullmodeler_naive_init(const void* c,const void* d,size_t n1,size_t n2)
{
	const struct pij_nullmodeler_naive_state* p=c;
	return pij_nullmodeler_naive_init_container(d,0,p->h);
}

void pij_nullmodeler_naive_input(void* m,const VECTORF* data)
{
	struct pij_nullmodeler_naive_state* p=m;
	size_t	i;
	for(i=0;i<data->size;i++)
		gsl_histogram_increment(p->h,VECTORFF(get)(data,i));
}

void pij_nullmodeler_naive_merge(const void* m,void* c)
{
	const struct pij_nullmodeler_naive_state* pm=m;
	struct pij_nullmodeler_naive_state* pc=c;
	gsl_histogram_add(pc->h,pm->h);
}

void pij_nullmodeler_naive_close(void* m)
{
	pij_nullmodeler_naive_close_container(m);
}
































