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
#include "config.h"
#include <string.h>
#include "macros.h"
#include "data_struct_ll.h"

int data_ll_init(struct data_ll* ll,size_t nmax)
{
	assert(ll);
	ll->nmax=nmax;
	ll->n=0;
	MALLOCSIZE(ll->d,2*nmax);
	if(!ll->d)
	{
		LOG(1,"Not enough memory.")
		return 1;
	}
	memset(ll->d,-1,2*nmax*sizeof(*ll->d));
	return 0;
}

void data_ll_free(struct data_ll* ll)
{
	assert(ll);
	if(ll->d)
	{
		free(ll->d);
		ll->d=0;
	}
	ll->nmax=0;
}

void data_ll_empty(struct data_ll* ll)
{
	assert(ll);
	ll->n=0;
	memset(ll->d,-1,2*ll->nmax*sizeof(*ll->d));
}



