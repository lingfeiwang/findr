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
/* This lib contains definition of linked list data structure
 */
 
#ifndef _HEADER_LIB_DATA_STRUCT_LL_H_
#define _HEADER_LIB_DATA_STRUCT_LL_H_
#include "config.h"
#include <stdlib.h>
#include <assert.h>
#include "logger.h"

#ifdef __cplusplus
extern "C"
{
#endif


//Linked list for size_t
struct data_ll
{
	//Max number of items
	size_t	nmax;
	//Current number of items
	size_t	n;
	/* Data and links
	 * d[2*i] is the child and d[2*i+1] is data.
	 * Item i has child j (at d[2*j] and d[2*j+1]) if d[2*i]=j.
	 * For data[2*i]=-1 is no child.
	 */
	size_t* restrict	d;
};

int data_ll_init(struct data_ll* ll,size_t nmax);
void data_ll_free(struct data_ll* ll);
void data_ll_empty(struct data_ll* ll);
//Insert entry with value val with no parent. Returns id.
static inline size_t data_ll_insert(struct data_ll* ll,size_t val);
//Insert entry with value val with parent id. Returns self id.
static inline size_t data_ll_insert_after(struct data_ll* ll,size_t id,size_t val);
/* Insert entry with value val with child id. Returns self id.
 * NOTE: Does not fix child of father of id.
 */
static inline size_t data_ll_insert_before(struct data_ll* ll,size_t id,size_t val);
//Return child id
static inline size_t data_ll_child(const struct data_ll* ll,size_t id);
//Return value
static inline size_t data_ll_val(const struct data_ll* ll,size_t id);

static inline size_t data_ll_insert(struct data_ll* ll,size_t val)
{
	size_t loc;
	
	if(ll->n==ll->nmax)
	{
		LOG(5,"Linked list insertion failed: linked list full.")
		return (size_t)-1;
	}
	loc=2*ll->n;
	ll->d[loc+1]=val;
	return ll->n++;
}

static inline size_t data_ll_insert_after(struct data_ll* ll,size_t id,size_t val)
{
	size_t loc;

	loc=data_ll_insert(ll,val);
	if(loc==(size_t)-1)
		return loc;
	ll->d[2*loc]=ll->d[2*id];
	ll->d[2*id]=loc;
	return loc;
}

static inline size_t data_ll_insert_before(struct data_ll* ll,size_t id,size_t val)
{
	size_t loc;
	
	loc=data_ll_insert(ll,val);
	if(loc==(size_t)-1)
		return loc;
	ll->d[2*loc]=id;
	return loc;
}

static inline size_t data_ll_child(const struct data_ll* ll,size_t id)
{
	assert(id<ll->n);
	return ll->d[2*id];
}

static inline size_t data_ll_val(const struct data_ll* ll,size_t id)
{
	assert(id<ll->n);
	return ll->d[2*id+1];
}


#ifdef __cplusplus
}
#endif
#endif
