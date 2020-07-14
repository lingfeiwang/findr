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
/* This lib contains general algorithms helpful in other c programs,
 * such as breaking strings, counting instances, and binary search.
 */

#ifndef _HEADER_LIB_GENERAL_ALG_H_
#define _HEADER_LIB_GENERAL_ALG_H_
#include "config.h"
#include <stdio.h>
#ifdef __cplusplus
extern "C"
{
#endif

/* Categorize data according to categorical information into separate arrays.
 * s:	Source of data.
 * c:	Categorical information. Each element contains a category of the corresponding element of s.
 * d:	Destination of categorization. Element i with c[i]=j is put into d[j].
 *		Modified value after this function indicates size of outcome arrays.
 * n:	Size of s and c.
 */
static inline void general_alg_categorize(const size_t* restrict s,const unsigned char* restrict c,size_t* restrict* restrict d,size_t n);

/* Categorize data according to embedded categorical information into separate arrays.
 * s:	Source of data.
 * c:	Categorical information. Each element contains a category of the corresponding element of s.
 * d:	Destination of categorization. Element i with c[s[i]]=j is put into d[j].
 *		Modified value after this function indicates size of outcome arrays.
 * n:	Size of s and c.
 */
static inline void general_alg_categorize_embed(const size_t* restrict s,const unsigned char* restrict c,size_t* restrict* restrict d,size_t n);

/* Removes duplicates in a sorted array of double, and shifts unique values to
 * the front of the array.
 * a:	array
 * n:	size of array
 * Return:	Size of new array
 */
static inline size_t remove_sorted_duplicates(double* restrict a,size_t n);






static inline void general_alg_categorize(const size_t* restrict s,const unsigned char* restrict c,size_t* restrict* restrict d,size_t n)
{
	size_t	i;
	for(i=0;i<n;i++)
		*(d[c[i]]++)=s[i];
}

static inline void general_alg_categorize_embed(const size_t* restrict s,const unsigned char* restrict c,size_t* restrict* restrict d,size_t n)
{
	size_t	i;
	for(i=0;i<n;i++)
		*(d[c[s[i]]]++)=s[i];
}

static inline size_t remove_sorted_duplicates(double* restrict a,size_t n)
{
	size_t p1,p2;
	
	if(!n)
		return 0;
	p1=p2=0;
	for(p1=p2=0;p2<n-1;p2++)
		if(a[p2]!=a[p2+1])
			a[p1++]=a[p2];
	a[p1++]=a[n-1];
	return p1;
}






#ifdef __cplusplus
}
#endif
#endif
