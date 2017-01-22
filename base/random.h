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
/* This file contains the low level randomization routines.
 */

#ifndef _HEADER_LIB_RANDOM_H_
#define _HEADER_LIB_RANDOM_H_
#include "config.h"
#include <time.h>
#include "gsl/rng.h"
#include "gsl/randist.h"
#include "logger.h"
#include "types.h"
#ifdef __cplusplus
extern "C"
{
#endif

extern gsl_rng* random_gen;

#define	random_new()		gsl_rng_alloc(gsl_rng_taus2)

static inline void random_init_any(gsl_rng** rng)
{
	*rng=random_new();
	if(!(*rng))
		LOG(1,"Can't allocate random number generator.")
}
#define	random_init()			random_init_any(&random_gen)

#define	random_seed_any(r,s)	gsl_rng_set(r,s)
#define	random_seed(s)			random_seed_any(random_gen,s)

#define	random_free_any(r)		gsl_rng_free(r)
#define	random_free()			random_free_any(&random_gen)

#define	random_seed_now_any(r)	random_seed_any(r,(unsigned long int)time(NULL))
#define	random_seed_now()		random_seed_now_any(random_gen)

// Generate uniformly distributed random number
#define random_uniform_any(r)	gsl_rng_uniform(r)
#define random_uniform()		random_uniform_any(random_gen)
#define random_uniformi_any(r,n)	gsl_rng_uniform_int(r,n)
#define random_uniformi(n)			random_uniformi_any(random_gen,n)
// Generate gaussian distributed random number
#define random_gaussian_any(r,sigma)	gsl_ran_gaussian(r,sigma)
#define random_gaussian(sigma)		random_gaussian_any(random_gen,sigma)

// Randomly shuffle items
#define random_shufflevf_any(r,f)	gsl_ran_shuffle(r,(f)->data,(f)->size,(f)->stride*sizeof(FTYPE))
#define random_shufflevf(f)			random_shufflevf_any(random_gen,f)

//Random shuffle gsl_permutation
#define random_shuffle_any(r,f)		gsl_ran_shuffle(r,(f)->data,(f)->size,sizeof(size_t))
#define random_shuffle(f)			random_shuffle_any(random_gen,f)


#ifdef __cplusplus
}
#endif
#endif
