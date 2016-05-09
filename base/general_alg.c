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
#include "config.h"
#include <stdlib.h>
#include <string.h>
#include "logger.h"
#include "macros.h"
#include "general_alg.h"

size_t break_string_into_substrings(char** s,const char* restrict sep,size_t nmax)
{
	size_t	i,l,n;
	
	l=strlen(s[0]);
	for(i=0,n=1;i<l;)
	{
		i+=strcspn(s[0]+i,sep);
		if(i>=l)
			break;
		if(n>=nmax)
			return nmax+1;
		s[0][i++]=0;
		i+=strspn(s[0]+i,sep);
		if(i>=l)
			break;
		s[n++]=s[0]+i;
	}
	return n;
}

long count_instance_fromfile(FILE* f,size_t size,size_t num,int (*func)(const void*))
{
#define	CLEANUP	AUTOFREE(dat)
	size_t	i,n;
	AUTOALLOC(char,dat,size,10000)
	size_t	sret;

	if(!dat)
	{
		LOG(3,"Can't allocate memory of size %lu.",size)
		return -1;
	}
	for(i=0,n=0;i<num;i++)
	{
		sret=fread(dat,size,1,f);
		if(sret!=1)
		{
			CLEANUP
			LOG(3,"Can't read file.")
			return -1;
		}
		n+=(func(dat)!=0);
	}
	CLEANUP
	return (long)n;
}































