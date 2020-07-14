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
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include "gsl/sort.h"
#include "logger.h"
#include "macros.h"
#include "threading.h"
#include "data_process.h"
#include "supernormalize.h"

void supernormalize_byrow_single_buffed(MATRIXF* m,gsl_permutation *p1,const FTYPE* restrict Pinv)
{
	size_t i,j;
	
	for(j=0;j<m->size1;j++)
	{
		VECTORFF(view) vvs=MATRIXFF(row)(m,j);
		
		//Rank
		CONCATENATE3(gsl_sort_vector,FTYPE_SUF,_index)(p1,&(vvs.vector));
		//Distribution
		for(i=0;i<m->size2;i++)
			MATRIXFF(set)(m,j,gsl_permutation_get(p1,i),Pinv[i]);
	}
	//Normalize again for unit variance
	MATRIXFF(normalize_row)(m);
}

int supernormalize_byrow_single(MATRIXF* m,const FTYPE* restrict Pinv)
{
	gsl_permutation *p1;
	
	p1=gsl_permutation_alloc(m->size2);
	if(!p1)
	{
		LOG(1,"Can't allocate permutations.")
		return 1;
	}
	supernormalize_byrow_single_buffed(m,p1,Pinv);
	gsl_permutation_free(p1);
	return 0;
}

void supernormalize_byrow_buffed(MATRIXF* m,gsl_permutation * const *p,FTYPE* Pinv)
{
	size_t	nth=(size_t)omp_get_max_threads();
	LOG(10,"Supernormalization started for matrix size ("PRINTFSIZET"*"PRINTFSIZET") on "PRINTFSIZET" threads.",m->size1,m->size2,nth)
	supernormalize_Pinv(m->size2,Pinv);

	#pragma omp parallel
	{
		size_t	nid=(size_t)omp_get_thread_num();
		size_t	n1,n2;
		MATRIXFF(view)	mv;
		
		threading_get_startend(m->size1,&n1,&n2);
		if(n2>n1)
		{
			mv=MATRIXFF(submatrix)(m,n1,0,n2-n1,m->size2);
			supernormalize_byrow_single_buffed(&mv.matrix,p[nid],Pinv);
		}
	}

	LOG(10,"Supernormalization completed.")
}

int supernormalize_byrow(MATRIXF* m)
{
#define CLEANUP	for(i=0;i<nth;i++)CLEANPERM(p[i])CLEANMEM(Pinv)

	size_t	nth=(size_t)omp_get_max_threads();
	size_t	i;
	int		ret;
	FTYPE	*Pinv;
	gsl_permutation* p[nth];
	
	
	MALLOCSIZE(Pinv,m->size2);
	ret=!!Pinv;
	for(i=0;i<nth;i++)
	{
		p[i]=gsl_permutation_alloc(m->size2);
		ret=ret&&p[i];
	}
	
	if(!ret)
		ERRRET("Not enough memory.")
	supernormalize_byrow_buffed(m,p,Pinv);
	CLEANUP
	return 0;
#undef	CLEANUP
}

int supernormalizef_byrow_single(MATRIXF* m,const FTYPE* restrict Pinv,FTYPE fluc)
{
	int	ret;
	ret=supernormalizea_byrow_single(m,Pinv);
	MATRIXFF(fluc)(m,fluc);
	MATRIXFF(normalize_row)(m);
	return ret;
}

void supernormalizef_byrow_buffed(MATRIXF* m,gsl_permutation * const *p,FTYPE* Pinv,FTYPE fluc)
{
	supernormalize_byrow_buffed(m,p,Pinv);
	MATRIXFF(fluc)(m,fluc);
	MATRIXFF(normalize_row)(m);
}

int supernormalizef_byrow(MATRIXF* m,FTYPE fluc)
{
	int		ret;
	ret=supernormalize_byrow(m);
	MATRIXFF(fluc)(m,fluc);
	MATRIXFF(normalize_row)(m);
	return ret;
}

void supernormalizer_byrow_single_buffed(MATRIXF* m,gsl_permutation *p1,VECTORF* vb,const gsl_rng* r)
{
	size_t i,j;
	VECTORFF(view) vvs;
	
	for(j=0;j<m->size1;j++)
	{
		//Random data
		for(i=0;i<m->size2;i++)
			VECTORFF(set)(vb,i,(FTYPE)random_gaussian_any(r,1));
		CONCATENATE2(gsl_sort_vector,FTYPE_SUF)(vb);
		
		//Rank
		vvs=MATRIXFF(row)(m,j);
		CONCATENATE3(gsl_sort_vector,FTYPE_SUF,_index)(p1,&(vvs.vector));
		//Distribution
		for(i=0;i<m->size2;i++)
			MATRIXFF(set)(m,j,gsl_permutation_get(p1,i),VECTORFF(get)(vb,i));
	}
	//Normalize again for unit variance
	MATRIXFF(normalize_row)(m);
}

void supernormalizer_byrow_buffed(MATRIXF* m,MATRIXF* mb,gsl_permutation * const *p,gsl_rng * const* rng)
{
	size_t	nth=(size_t)omp_get_max_threads();
	LOG(10,"Randomized normalization started for matrix size ("PRINTFSIZET"*"PRINTFSIZET") on "PRINTFSIZET" threads.",m->size1,m->size2,nth)

	#pragma omp parallel
	{
		size_t	nid=(size_t)omp_get_thread_num();
		size_t	n1,n2;
		MATRIXFF(view)	mv;
		VECTORFF(view)	vv;
		
		threading_get_startend(m->size1,&n1,&n2);
		if(n2>n1)
		{
			mv=MATRIXFF(submatrix)(m,n1,0,n2-n1,m->size2);
			vv=MATRIXFF(row)(mb,nid);
			supernormalizer_byrow_single_buffed(&mv.matrix,p[nid],&vv.vector,rng[nid]);
		}
	}

	LOG(10,"Randomized normalization completed.")
}

int supernormalizer_byrow(MATRIXF* m)
{
#define CLEANUP	for(i=0;i<nth;i++){CLEANPERM(p[i])\
				if(r[i]){random_free_any(r[i]);r[i]=0;}}\
				CLEANMATF(mb)

	size_t	nth=(size_t)omp_get_max_threads();
	size_t	i;
	int		ret;
	gsl_permutation* p[nth];
	MATRIXF	*mb;
	gsl_rng	*r[nth];
	
	mb=MATRIXFF(alloc)(nth,m->size2);
	ret=!!mb;
	for(i=0;i<nth;i++)
	{
		p[i]=gsl_permutation_alloc(m->size2);
		r[i]=random_new();
		ret=ret&&p[i]&&r[i];
	}
	if(!ret)
		ERRRET("Not enough memory.")
	random_seed_any(r[0],(size_t)time(NULL));
	for(i=1;i<nth;i++)
		random_seed_any(r[i],gsl_rng_get(r[0])+2314653);
	supernormalizer_byrow_buffed(m,mb,p,r);
	CLEANUP
	return 0;
#undef	CLEANUP
}


















































































































