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
#include <string.h>
#include "logger.h"
#include "macros.h"
#include "data_struct_heap.h"

int data_heap_init(struct data_heap* h,size_t nmax)
{
	h->nmax=nmax;
	h->n=0;
	MALLOCSIZE(h->d,nmax);
	return !h->d;
}

void data_heap_free(struct data_heap* h)
{
	h->nmax=0;
	if(h->d)
		free(h->d);
	h->d=0;
}

void data_heap_empty(struct data_heap* h)
{
	h->n=0;
}

int data_heap_push(struct data_heap* h, HTYPE d)
{
	size_t c,p;
	
	if(h->nmax==h->n)
	{
		LOG(10,"Heap push failed: heap full.")
		return 1;
	}
	h->d[h->n]=d;
	c=h->n++;
	while(c)
	{
		p=(c-1)/2;
		if(d>h->d[p])
			return 0;
		h->d[c]=h->d[p];
		h->d[p]=d;
		c=p;
	}
	return 0;
}

HTYPE data_heap_pop(struct data_heap* h)
{
	size_t	p,c2,cm,pn;
	HTYPE	v,ret;
	
	assert(h->n);
	ret=h->d[0];
	h->d[0]=h->d[--(h->n)];
	if(!h->n)
		return ret;
	p=0;
	pn=h->n/2;
	while(p+1<pn)
	{
		cm=2*p+1;
		c2=cm+1;
		if(h->d[c2]<h->d[cm])
			cm=c2;
		if(h->d[p]<h->d[cm])
			return ret;
		v=h->d[p];
		h->d[p]=h->d[cm];
		h->d[cm]=v;
		p=cm;
	}
	if(p+1==pn)
	{
		cm=2*p+1;
		c2=cm+1;
		if((c2<h->n)&&(h->d[c2]<h->d[cm]))
			cm=c2;
		if(h->d[p]<h->d[cm])
			return ret;
		v=h->d[p];
		h->d[p]=h->d[cm];
		h->d[cm]=v;
	}
	return ret;
}

int data_heapdec_push(struct data_heapdec* h, HTYPE d)
{
	size_t c,p;
	
	if(h->nmax==h->n)
	{
		LOG(10,"Heap push failed: heap full.")
		return 1;
	}
	h->d[h->n]=d;
	c=h->n++;
	while(c)
	{
		p=(c-1)/2;
		if(d<h->d[p])
			return 0;
		h->d[c]=h->d[p];
		h->d[p]=d;
		c=p;
	}
	return 0;
}

HTYPE data_heapdec_pop(struct data_heapdec* h)
{
	size_t	p,c2,cm,pn;
	HTYPE	v,ret;
	
	assert(h->n);
	ret=h->d[0];
	h->d[0]=h->d[--(h->n)];
	if(!h->n)
		return ret;
	p=0;
	pn=h->n/2;
	while(p+1<pn)
	{
		cm=2*p+1;
		c2=cm+1;
		if(h->d[c2]>h->d[cm])
			cm=c2;
		if(h->d[p]>h->d[cm])
			return ret;
		v=h->d[p];
		h->d[p]=h->d[cm];
		h->d[cm]=v;
		p=cm;
	}
	if(p+1==pn)
	{
		cm=2*p+1;
		c2=cm+1;
		if((c2<h->n)&&(h->d[c2]>h->d[cm]))
			cm=c2;
		if(h->d[p]>h->d[cm])
			return ret;
		v=h->d[p];
		h->d[p]=h->d[cm];
		h->d[cm]=v;
	}
	return ret;
}



