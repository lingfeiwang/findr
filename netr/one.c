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
#include <stdlib.h>
#include <time.h>
#include "../base/gsl/math.h"
#include "../base/gsl/permutation.h"
#include "../base/gsl/sort.h"
#include "../base/logger.h"
#include "../base/macros.h"
#include "../base/data_process.h"
#include "../cycle/cycle.h"
#include "one.h"


size_t netr_one_greedy(const MATRIXF* p,MATRIXUC* net,size_t nam,size_t nimax,size_t nomax)
{
#define CLEANUP	CLEANVECF(v)CYCLEF(free)(&cs);CLEANPERM(perm)
#define	TOID(N,V1,V2)	V1=(N)/(n-1);V2=(N)%(n-1);if((V2)>=(V1))V2++;

	struct CYCLEF(system)	cs;
	VECTORF*	v=0;
	gsl_permutation*	perm=0;
	int	ret;
	size_t	n,na,i,ntot;



	//Initialize
	n=p->size1;
	assert((n==p->size2)&&(n==net->size1)&&(n==net->size2));
	assert(nimax&&nomax&&nam);
	ntot=n*(n-1);
	nam=GSL_MIN(ntot/2,nam);
	{
		size_t	t1;
		t1=GSL_MIN(nimax,nomax);
		if(nam/n>=t1)
			nam=t1*n;
	}
	ret=CYCLEF(init)(&cs,n,nam);
	if(ret)
		ERRRETV(0,"Failed to initialize cycle detection.")
	cs.nim=nimax;
	cs.nom=nomax;
	v=VECTORFF(alloc)(ntot);
	perm=gsl_permutation_alloc(ntot);
	if(!(v&&perm))
		ERRRETV(0,"Not enough memory.")
	
	//Obtain edge order
	MATRIXFF(flatten_nodiag)(p,v);
	ret=CONCATENATE3(gsl_sort_vector,FTYPE_SUF,_index)(perm,v);
	if(ret)
		ERRRETV(0,"Failed to sort vector.")
	CLEANVECF(v)
	
	//Add edges
	for(i=0,na=0;(i<ntot)&&(na<nam);i++)
	{
		size_t	v0,v1,v2;
		v0=gsl_permutation_get(perm,ntot-i-1);
		TOID(v0,v1,v2)
		ret=cycle_vg_add(&cs,v1,v2);
		na+=!ret;
	}
	CYCLEF(extract_graph)(&cs,net);
	
	CLEANUP
	return na;
#undef	TOID
#undef	CLEANUP
}

size_t netr_one_greedy_info(const MATRIXF* p,MATRIXL* net,MATRIXD* time,size_t nam,size_t nimax,size_t nomax)
{
#define CLEANUP	CLEANVECF(v)CYCLEF(free)(&cs);CLEANPERM(perm)
#define	TOID(N,V1,V2)	V1=(N)/(n-1);V2=(N)%(n-1);if((V2)>=(V1))V2++;

	struct CYCLEF(system)	cs;
	VECTORF*	v=0;
	gsl_permutation*	perm=0;
	int	ret;
	size_t	n,na,i,ntot;
	int	sign[2]={1,-1};
	clock_t	cstart,cnow;

	//Initialize
	n=p->size1;
	assert((n==p->size2)&&(n==net->size1)&&(n==net->size2));
	assert(nimax&&nomax);
	ntot=n*(n-1);
	nam=GSL_MIN(ntot/2,nam);
	ret=CYCLEF(init)(&cs,n,nam);
	if(ret)
		ERRRETV(0,"Failed to initialize cycle detection.")	
	cs.nim=nimax;
	cs.nom=nomax;
	v=VECTORFF(alloc)(ntot);
	perm=gsl_permutation_alloc(ntot);
	if(!(v&&perm))
		ERRRETV(0,"Not enough memory.")
	
	//Obtain edge order
	MATRIXFF(flatten_nodiag)(p,v);
	ret=CONCATENATE3(gsl_sort_vector,FTYPE_SUF,_index)(perm,v);
	if(ret)
		ERRRETV(0,"Failed to sort vector.")
	CLEANVECF(v)
	
	//Add edges
	MATRIXLF(set_zero)(net);
	cstart=clock();
	MATRIXDF(set_all)(time,(double)cstart);
	for(i=0,na=0;(i<ntot)&&(na<nam);i++)
	{
		size_t	v0,v1,v2;
		v0=gsl_permutation_get(perm,ntot-i-1);
		TOID(v0,v1,v2)
		ret=cycle_vg_add(&cs,v1,v2);
		cnow=clock();
		MATRIXLF(set)(net,v1,v2,(int)(i+1)*sign[ret]);
		MATRIXDF(set)(time,v1,v2,(double)cnow);
		na+=!ret;
	}
	MATRIXDF(add_constant)(time,(double)-cstart);
	MATRIXDF(scale)(time,1/(double)CLOCKS_PER_SEC);
	
	CLEANUP
	return na;
#undef	TOID
#undef	CLEANUP
}























