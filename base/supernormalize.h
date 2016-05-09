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
/* This is the header file for supernormalization, i.e. transforming
 * samples of a variable to normal distribution N(0,1). Two method
 * are provided: deterministic and random.
 */

#ifndef _HEADER_LIB_SUPERNORMALIZE_H_
#define _HEADER_LIB_SUPERNORMALIZE_H_
#include "config.h"
#include "gsl/permutation.h"
#include "gsl/cdf.h"
#include "random.h"
#include "types.h"
#ifdef __cplusplus
extern "C"
{
#endif

/**********************************************************************
 * Deterministic supernormalization
 **********************************************************************/

//Struct of info pack for each supernormalize thread
struct supernormalize_byrow_threadinfo{
	MATRIXF*		m;		//matrix to be supernormalized
	const FTYPE*	Pinv;	//Inverse cdf
};

/* Supernormalize matrix per row with single thread and buff provided.
 * Supernormalization takes place by converting the existing data into a normal distribution
 * with 0 mean and 1 variance. Due to numerical errors, their values may be inexact. This is performed
 * by first converting data into their ranking, and assign new values according to the cummulative
 * distribution function of the respective fraction. After that, a normalization is perform to scale
 * the new data into 0 mean and 1 variance.
 * m:		Matrix to be supernormalized. Overwrites data.
 * p1:		Permutation objects for ranking conversion
 * Pinv:	Inverse transformation from ranking to normal distribution
 * 			(precalculated CDF values of normal distribution of the respective ranking)
 */
void supernormalize_byrow_single_buffed(MATRIXF* m,gsl_permutation *p1,const FTYPE* restrict Pinv);

/* Supernormalize matrix per row with single thread and buff provided.
 * See supernormalize_byrow_single_buffed for detail.
 * m:		Matrix to be supernormalized. Overwrites data.
 * Pinv:	Inverse transformation from ranking to normal distribution
 * 			(precalculated CDF values of normal distribution of the respective ranking)
 * Return:	0 on success.
 */
static inline int supernormalize_byrow_single(MATRIXF* m,const FTYPE* restrict Pinv);

/* Obtain Inverse CDF for normal distribution of n fractiles.
 * n:		n
 * Pinv:	Inverse CDF of normal distribution. Return[i]=CDF^(-1)((i+1)/(n+1)).
 */
static inline void supernormalize_Pinv(size_t n,FTYPE*restrict Pinv);


/* Supernormalizes and overwrites each row of matrix m.
 * Supernormalize into 0 mean and 1 variance, and fulfills normal distribution
 * Therefore numbers are assigned purely according to the rankings.
 * Uses multiple threads
 * Ties are ordered sequentially by GSL (potential increased correlation between rows)
 * With or without buffer included:
 * m:	(n1,n2) Matrix to be supernormalized
 * p:	(nth) permutation buffer
 * Pinv:Buffer to calculate and place inverse CDF
 * nth:	Number of threads.
 * Return:	0 if success.
 */
void supernormalize_byrow_buffed(MATRIXF* m,gsl_permutation * const *p,FTYPE* Pinv);
int supernormalize_byrow(MATRIXF* m);


/**********************************************************************
 * Random supernormalization
 **********************************************************************/

//Struct of info pack for each supernormalize thread
struct supernormalizer_byrow_threadinfo{
	MATRIXF	*m;		//matrix to be supernormalized
	//VECTORF	*vb;	//Random data buff
	//gsl_rng	*r;		//Random Number generator
};

//Check their supernormalize counterparts for definition.
void supernormalizer_byrow_buffed(MATRIXF* m,MATRIXF* mb,gsl_permutation * const *p,gsl_rng * const* rng);

int supernormalizer_byrow(MATRIXF* m);


/**********************************************************************
 * Static functions
 **********************************************************************/

static inline int supernormalize_byrow_single(MATRIXF* m,const FTYPE* restrict Pinv)
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

static inline void supernormalize_Pinv(size_t n,FTYPE* restrict Pinv)
{
	size_t	i;

	for(i=0;i<n;i++)
		Pinv[i]=(FTYPE)gsl_cdf_gaussian_Pinv(((FTYPE)(i+1))/(FTYPE)(n+1),1);
}

#ifdef __cplusplus
}
#endif
#endif


























