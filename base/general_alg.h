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

/* Break a long string into substrings.
 * The initial pointer s[0] points to the long string
 * It is broken by the any of the separator chars sep.
 * Multiple continous separators are treated as one.
 * The broken strings are stored in s[0,1,2,...].
 * Maximum broken string count is nmax.
 * Return the number of strings broken into.
 * If return value > nmax, the last string s[nmax-1] contains
 * separators because insufficient breaking.
 */
size_t break_string_into_substrings(char** s,const char* restrict sep,size_t nmax);

/* Counts the number of instances that gives the criterion function a return value of true.
 * Parameters:
 * inst:	the pointer to the first item of the instance array
 * size:	size of one instance
 * num:		the number of instances to be tested
 * func:	the criterion function to test the instances.
 * Return:	the number of instances that gives true for the function func, for a total of num
 * 			instances, starting from inst, each of which has size bytes.
 */
static inline size_t count_instance(const void* inst,size_t size,size_t num,int (*func)(const void*));

/* Similar with above (count_instance), but processes a file handle.
 * Return:	the number of instances satisfying the criterion function func if succeed,
 * 		or -1 if fail.
 */
long count_instance_fromfile(FILE* f,size_t size,size_t num,int (*func)(const void*));

/* Search for left insert position of item in ascending array.
 * x:	item to search for position
 * a:	ascending array 
 * n:	size of a
 * Return:	left insert position of x (=0,...,n)
 */
static inline size_t bsearch_d(double x,const double* restrict a,size_t n);

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






static inline size_t count_instance(const void* inst,size_t size,size_t num,int (*func)(const void*))
{
	size_t i,n;
	for(i=0,n=0;i<num;i++)
		n+=(func(((char*)inst)+i*size)!=0);
	return n;
}

static inline size_t bsearch_d(double x,const double* restrict a,size_t n)
{
	size_t head,end,mid;
	
	head=0;
	end=n;
	while(head<end)
	{
		mid=(head+end)/2;
		if(a[mid]>=x)
			end=mid;
		else if(a[mid]<x)
			head=mid+1;
	}
	return head;
}

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
