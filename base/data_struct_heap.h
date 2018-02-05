/* Copyright 2016-2018 Lingfei Wang
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
/* This lib contains definition of heap data structure
 */
 
#ifndef _HEADER_LIB_DATA_STRUCT_HEAP_H_
#define _HEADER_LIB_DATA_STRUCT_HEAP_H_
#include "config.h"
#include <stdlib.h>
#include <assert.h>

#ifdef __cplusplus
extern "C"
{
#endif


#define HTYPE	size_t

//Incremental heap
struct data_heap
{
	size_t nmax;
	size_t n;
	HTYPE* restrict	d;
};

int data_heap_init(struct data_heap* h,size_t nmax);
void data_heap_free(struct data_heap* h);
void data_heap_empty(struct data_heap* h);
int data_heap_push(struct data_heap* h, HTYPE d);
HTYPE data_heap_pop(struct data_heap* h);
// int data_heap_popto(struct data_heap* h, HTYPE* d);
static inline HTYPE data_heap_get(const struct data_heap* h,size_t n);
static inline HTYPE data_heap_top(const struct data_heap* h);

//Decremental heap
#define data_heapdec data_heap
#define data_heapdec_init data_heap_init
#define data_heapdec_free data_heap_free
#define data_heapdec_empty data_heap_empty
int data_heapdec_push(struct data_heapdec* h, HTYPE d);
HTYPE data_heapdec_pop(struct data_heapdec* h);
#define data_heapdec_get data_heap_get
#define data_heapdec_top data_heap_top


static inline HTYPE data_heap_get(const struct data_heap* h,size_t n)
{
	assert(h->n>n);
	return h->d[n];
}

static inline HTYPE data_heap_top(const struct data_heap* h)
{
	return data_heap_get(h,0);
}

#ifdef __cplusplus
}
#endif
#endif
