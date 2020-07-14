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
// This file contains the functions for multi-threading
#ifndef _HEADER_LIB_THREADING_H_
#define _HEADER_LIB_THREADING_H_
#include "config.h"
#include <stdlib.h>
#include <omp.h>

#ifdef __cplusplus
extern "C"
{
#endif


/* Calculate the split position of the big problem into smaller ones.
 * ntotal:	Total size of the problem
 * nthread:	Total number of threads
 * x:		ID of current thread
 * Return:	The start position of problem id by thread x
 */
static inline size_t threading_get_start_bare(size_t ntotal,size_t nthread,size_t x);

/* Calculate the start and end position of the big problem for any thread.
 * ntotal:	Total size of the problem
 * start,
 * end:		Return location of start and end positions for any thread
 * id:		ID of current thread
 * ida:		Total number of threads
 */
static inline void threading_get_startend_from(size_t ntotal,size_t *start,size_t *end,size_t id,size_t ida);

/* Calculate the start and end position of the big problem for current thread (with openMP)
 * ntotal:	Total size of the problem
 * start,
 * end:		Return location of start and end positions for current thread
 */
static inline void threading_get_startend(size_t ntotal,size_t *start,size_t *end);


static inline size_t threading_get_start_bare(size_t ntotal,size_t nthread,size_t x)
{
	size_t i,j;
	i=ntotal/nthread;
	j=ntotal-i*nthread;
	if(j>x)
		j=x;
	j+=i*x;
	if(j>ntotal)
		j=ntotal;
	return j;
}

static inline void threading_get_startend_from(size_t ntotal,size_t *start,size_t *end,size_t id,size_t ida)
{
	*start=threading_get_start_bare(ntotal,ida,id);
	*end=threading_get_start_bare(ntotal,ida,id+1);
}

static inline void threading_get_startend(size_t ntotal,size_t *start,size_t *end)
{
	size_t	id=(size_t)omp_get_thread_num();
	size_t	ida=(size_t)omp_get_num_threads();
	threading_get_startend_from(ntotal,start,end,id,ida);
}

#ifdef __cplusplus
}
#endif

#endif



























