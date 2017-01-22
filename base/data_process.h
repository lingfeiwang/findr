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
/* This lib contains the data processing routines.
 * Such as normalization, covariance, off-diagonal matrix flatten, etc.
 */

#ifndef _HEADER_LIB_DATA_PROCESS_H_
#define _HEADER_LIB_DATA_PROCESS_H_
#include "config.h"
#include <string.h>
#include <assert.h>
#include "gsl/math.h"
#include "gsl/blas.h"
#include "gsl/permutation.h"
#include "types.h"
#include "random.h"
#ifdef __cplusplus
extern "C"
{
#endif

/* Constructs data matrix from existing data array (data) with row count (nrow),
 * column count (ncol). The matrix information is put into existing matrix struct
 * (dest).
 */
void MATRIXFF(from_dense)(MATRIXF* dest,const FTYPE* restrict data,size_t nrow,size_t ncol);

/* Constructs data matrxi from existing dense data file handle (f) with row count
 * (nrow), column cout (ncol). The matrix pointer is returned on success, or 0 on fail.
 */
MATRIXF* MATRIXFF(from_densefile)(FILE* f,size_t nrow,size_t ncol);

/* Normalize matrix (m) for each row, with buffer vectors v1, v2, of sizes
 * m->size1 and m->size2 respectively.
 * m:	(n1,n2) Matrix to be normalized row by row.
 * v1:	(n1) Buffer
 * v2:	(n2) Buffer. Must be all 1 initially
 */
void MATRIXFF(normalize_row_buffed)(MATRIXF* m,VECTORF* v1,const VECTORF* v2);

/* Normalize matrix (m) for each row, with temporarily allocated buffer.
 * Invokes MATRIXFF(normalize_row_buffed). Return 0 on success.
 */
int MATRIXFF(normalize_row)(MATRIXF* m);

/* Calculates the covariance of matrix (m) and store the results in the
 * existing proper sized matrix cov. cov=m*m^T.
 */
static inline void MATRIXFF(cov1)(const MATRIXF* m,MATRIXF* cov);

/* Calculates the covariance between matrices m1 and m2,
 * and store the results in the
 * existing proper sized matrix cov. cov=m1*m2^T.
 */
static inline void MATRIXFF(cov2)(const MATRIXF* m1,const MATRIXF* m2,MATRIXF* cov);

/* Calculates the row-by-row covariance between matrices m1 and m2,
 * and store the results in the
 * existing proper sized vector cov.
 * m1:	(n1,n2)
 * m2:	(n1,n2)
 * cov:	(n1). cov[i]=m1[i]*m2[i]^T/n2
 */
static inline void MATRIXFF(cov2_1v1)(const MATRIXF* m1,const MATRIXF* m2,VECTORF* cov);

/* Calculates the covariance matrix of m, and store the result into cov.
 * Results are bounded at -1 to 1.
 */
static inline void MATRIXFF(cov1_bounded)(const MATRIXF* m,MATRIXF* cov);

/* Calculates the covariance matrix between m1 and m2, and store the result into cov.
 * Results are bounded at -1 to 1.
 */
static inline void MATRIXFF(cov2_bounded)(const MATRIXF* m1,const MATRIXF* m2,MATRIXF* cov);

/* Calculates the row-by-row covariance between matrices m1 and m2,
 * and store the results in the
 * existing proper sized vector cov. Results are capped within [-1,1]
 * m1:	(n1,n2)
 * m2:	(n1,n2)
 * cov:	(n1). cov[i]=m1[i]*m2[i]^T/n2
 */
static inline void MATRIXFF(cov2_1v1_bounded)(const MATRIXF* m1,const MATRIXF* m2,VECTORF* cov);

/* Flattens diagonal matrix into vector.
 * m:	diagonal matrix, size (n1,n2)
 * v:	vector for output, size (n1*n2)
 * Return:	0 on success.
 */
static inline void MATRIXFF(flatten)(const MATRIXF* m,VECTORF* v);

/* Flattens diagonal matrix into vector, but skipping diagonal components.
 * m:	diagonal matrix, size (n1,n2)
 * v:	vector for output, size n1*n2-min(n1,n2)
 * Return:	0 on success.
 */
void MATRIXFF(flatten_nodiag)(const MATRIXF* m,VECTORF* v);

/* Wraps vector into diagonal matrix
 * v:	vector, size n1*n2
 * m:	diagonal matrix for output, size (n1,n2)
 * Return:	0 on success
 */
static inline void VECTORFF(wrap)(const VECTORF* v,MATRIXF* m);

/* Wraps vector into diagonal matrix, but skipping diagonal components.
 * v:	vector, size n1*n2-min(n1,n2)
 * m:	diagonal matrix for output, size (n1,n2)
 * Return:	0 on success
 */
void VECTORFF(wrap_nodiag)(const VECTORF* v,MATRIXF* m);

/* Strips diagonal elements (i=j) of (row-first) matrix, by overwriting diagonal elements
 * with end elements of the same matrix. So order of elements is changed.
 * d:	pointer of matrix data, row-first
 * n1:	row number of matrix
 * n2:	column number of matrix
 * s:	size of each element
 * Return:	remaining number of elements of d after stripping. They are stored
 * 			as the frontmost elements in d. Or 0 if failed.
 */
static inline size_t strip_diag(void* restrict data,size_t n1,size_t n2,size_t s);

/* Bounds vector from below so it is never below the value.
 */
static inline void VECTOROF(bound_below)(VECTORO* m,float v);
static inline void VECTORDF(bound_below)(VECTORD* m,double v);

/* Bounds vector from above so it is never above the value.
 */
static inline void VECTOROF(bound_above)(VECTORO* m,float v);
static inline void VECTORDF(bound_above)(VECTORD* m,double v);

/* Bounds vector from above and below so it is always between the values.
 */
static inline void VECTOROF(bound_both)(VECTORO* m,float vlow,float vhigh);
static inline void VECTORDF(bound_both)(VECTORD* m,double vlow,double vhigh);

/* Bounds matrix from below so it is never below the value.
 */
static inline void MATRIXOF(bound_below)(MATRIXO* m,float v);
static inline void MATRIXDF(bound_below)(MATRIXD* m,double v);

/* Bounds matrix from above so it is never above the value.
 */
static inline void MATRIXOF(bound_above)(MATRIXO* m,float v);
static inline void MATRIXDF(bound_above)(MATRIXD* m,double v);

/* Bounds matrix from above and below so it is always between the values.
 */
static inline void MATRIXOF(bound_both)(MATRIXO* m,float vlow,float vhigh);
static inline void MATRIXDF(bound_both)(MATRIXD* m,double vlow,double vhigh);

/* Sets all elements of vector that satisfies the condition to one single value.
 * v:		vector
 * func:	condition to satisfy to change value
 * val:		new value
 */
static inline void VECTOROF(set_cond)(VECTORO* v,int (*func)(float),float val);
static inline void VECTORDF(set_cond)(VECTORD* v,int (*func)(double),double val);
/* Sets all elements of matrix that satisfies the condition to one single value.
 * m:		matrix
 * func:	condition to satisfy to change value
 * val:		new value
 */
static inline void MATRIXOF(set_cond)(MATRIXO* v,int (*func)(float),float val);
static inline void MATRIXDF(set_cond)(MATRIXD* v,int (*func)(double),double val);

// Sets all inf or nan to new value
static inline void VECTOROF(set_inf)(VECTORO* v,float val);
static inline void VECTOROF(set_nan)(VECTORO* v,float val);
static inline void VECTORDF(set_inf)(VECTORD* v,double val);
static inline void VECTORDF(set_nan)(VECTORD* v,double val);
static inline void MATRIXOF(set_inf)(MATRIXO* m,float val);
static inline void MATRIXOF(set_nan)(MATRIXO* m,float val);
static inline void MATRIXDF(set_inf)(MATRIXD* m,double val);
static inline void MATRIXDF(set_nan)(MATRIXD* m,double val);

/* Locates the first not number in a vector/matrix, and return its location.
 * If not found, return -1.
 */
static inline long VECTORFF(first_nan)(const VECTORF *restrict v);
static inline long MATRIXFF(first_nan)(const MATRIXF *restrict m);

/* Count the number of values for each row
 * g:	(ng,ns)	Data
 * ans:	(ng) Output vector for number of values for row.
 * vb:	(nv) Buffer vector.
 * nv:	Max number of values each element can take. g[i,j]=0,...,nv-1
 */
void MATRIXGF(countv_byrow_buffed)(const MATRIXG* g,VECTORG* ans,VECTORUC* vb);

/* Obtain occurance ratio for each possible value in a vector.
 * d:	(n)	Data
 * ans:	(nv+1) Output vector for occurance ratio of each value.
 * nv:	Max number of values each element can take. d[i]=0,...,nv
 */ 
void VECTORGF(count_ratio_d)(const VECTORG* d,VECTORD* ans);

/* Scale each row of matrix by different values from the vector.
 * m:	(n1,n2) Matrix to be scaled.
 * v:	(n1)
 * After scaling, new m[i,j]=m[i,j]*v[i]
 */
static inline void MATRIXFF(scale_row)(MATRIXF* m,const VECTORF* v);

/* Scale each row of matrix by different values from the vector and function.
 * m:	(n1,n2) Matrix to be scaled.
 * v:	(n1)
 * f:	Function to apply on vector before scaling
 * After scaling, new m[i,j]=m[i,j]*f(v[i])
 */
static inline void MATRIXFF(scale_row_func)(MATRIXF* m,const VECTORF* v,FTYPE (*f)(FTYPE));

/* Scale each column of matrix by different values from the vector.
 * m:	(n1,n2) Matrix to be scaled.
 * v:	(n2)
 * After scaling, new m[i,j]=m[i,j]*v[j]
 */
static inline void MATRIXFF(scale_column)(MATRIXF* m,const VECTORF* v);

/* Scale each column of matrix by different values from the vector and function.
 * m:	(n1,n2) Matrix to be scaled.
 * v:	(n2)
 * f:	Function to apply on vector before scaling
 * After scaling, new m[i,j]=m[i,j]*f(v[j])
 */
static inline void MATRIXFF(scale_column_func)(MATRIXF* m,const VECTORF* v,FTYPE (*f)(FTYPE));

/* Select rows of matrix according to index vector and save into new matrix.
 * d:	(n1,n2) Source matrix.
 * dest:(nx,n2) Destination matrix.
 * c:	(n1) Condition vector.
 * Return: n as number of selected rows.
 * The rows of d whose corresponding elements in c are nonzero are selected
 * and put into dest in original order, filling first n rows of dest.
 * dest size must be large enough (nx>=n).
 */
size_t	MATRIXGF(rows_save)(const MATRIXG* d,MATRIXG* dest,const VECTORUC* c);
size_t	MATRIXFF(rows_save)(const MATRIXF* d,MATRIXF* dest,const VECTORUC* c);

/* Select rows of matrix according to index vector, discard diagonal elements,
 * and save into new vector.
 * d:	(n1,n2) Source matrix.
 * dest:(nx) Destination vector.
 * c:	(n1) Condition vector.
 * Return: n as number of selected elements.
 * The rows of d whose corresponding elements in c are nonzero are selected
 * and put into dest in original order after discarding diagonal elements.
 * This fills first n elements of dest.
 * dest size must be large enough (nx>=n). However always n<=n1*n2-min(n1,n2).
 */
size_t	MATRIXFF(rows_save_nodiag)(const MATRIXF* d,VECTORF* dest,const VECTORUC* c);

/* Select rows of matrix according to index vector and load from source matrix.
 * d:	(nx,n2) Source matrix.
 * dest:(n1,n2) Destination matrix.
 * c:	(n1) Condition vector.
 * Return: n as number of selected rows.
 * The rows of dest whose corresponding elements in c are nonzero are selected
 * and loaded from d in original order, using first n rows of d.
 * d size must be large enough (nx>=n).
 */
size_t	MATRIXGF(rows_load)(const MATRIXG* d,MATRIXG* dest,const VECTORUC* c);
size_t	MATRIXFF(rows_load)(const MATRIXF* d,MATRIXF* dest,const VECTORUC* c);

/* Select rows of matrix according to index vector and load from source vector,
 * while skipping diagonal elements.
 * d:	(nx) Source vector.
 * dest:(n1,n2) Destination matrix.
 * c:	(n1) Condition vector.
 * Return: n as number of selected elements.
 * The rows of dest whose corresponding elements in c are nonzero are selected
 * and loaded from d in original order after skipping diagonal elements.
 * This loads first n elements of d.
 * d size must be large enough (nx>=n). However always n<=n1*n2-min(n1,n2).
 */
size_t	MATRIXFF(rows_load_nodiag)(const VECTORF* d,MATRIXF* dest,const VECTORUC* c);

/* Permute matrix columns with given permutation and buffer provided.
 * After permutation new_matrix[i]=old_matrix[perm[i]]
 * m:	(n1,n2)
 * perm:Permutation
 * b1:	(n1) Buffer
 * b2:	(n2) Buffer
 */
void MATRIXFF(permute_column_buffed)(MATRIXF* m,const gsl_permutation* perm,VECTORF* b1,VECTORUC* b2);

/* Permute matrix rows with given permutation and buffer provided.
 * After permutation new_matrix[i]=old_matrix[perm[i]]
 * m:	(n1,n2)
 * perm:Permutation
 * b1:	(n2) Buffer
 * b2:	(n1) Buffer
 */
void MATRIXFF(permute_row_buffed)(MATRIXF* m,const gsl_permutation* perm,VECTORF* b1,VECTORUC* b2);

/* Calculates equality relations for every element and save answer in another
 * vector.
 * d:	(n) Source vector
 * ans:	(n) Answer vector
 * val:	Value to compare against.
 * ans[i]=(d[i]==val)
 */
static inline void VECTORGF(eq)(const VECTORG *d,VECTORUC* ans,GTYPE val);

/* Calculates difference on vector.
 * dnew[n]=dold[n]-dold[n-1], dnew[0]=dold[0].
 * d:	(n) Source and destination vector
 */
static inline void VECTORDF(diff)(VECTORD* d);
static inline void VECTOROF(diff)(VECTORF* d);

/* Calculates cumulative sum on vector.
 * d:	(n) Source and destination vector
 */
static inline void VECTORDF(cumsum)(VECTORD* d);
static inline void VECTOROF(cumsum)(VECTORF* d);

/* Randomly fluctuate every element in matrix relatively.
 * m:		Matrix to fluctuate
 * fluc:	Relative amplitude of fluctuation. Every element x is replaced with
 * 			x*(1+y*fluc), where y is uniformly distributed in [-1,1)
 */
static inline void MATRIXFF(fluc)(MATRIXF* m,FTYPE fluc);

/* Randomly fluctuate every element in vector relatively.
 * v:		Vector to fluctuate
 * fluc:	Relative amplitude of fluctuation. Every element x is replaced with
 * 			x*(1+y*fluc), where y is uniformly distributed in [-1,1)
 */
static inline void VECTORFF(fluc)(VECTORF* v,FTYPE fluc);

static inline void MATRIXFF(flatten)(const MATRIXF* m,VECTORF* v)
{
	size_t	i,n1,n2;
	VECTORFF(view)	vv2;
	
	//Initialize
	n1=m->size1;
	n2=m->size2;
	assert(v->size==n1*n2);
	
	for(i=0;i<n1;i++)
	{
		VECTORFF(const_view) vv1=MATRIXFF(const_row)(m,i);
		vv2=VECTORFF(subvector)(v,i*n2,n2);
		VECTORFF(memcpy)(&vv2.vector,&vv1.vector);
	}
}

static inline void VECTORFF(wrap)(const VECTORF* v,MATRIXF* m)
{
	size_t	n1,n2;

	//Initialize
	n1=m->size1;
	n2=m->size2;
	assert(v->size==n1*n2);
	MATRIXFF(const_view)	mv=MATRIXFF(const_view_vector)(v,n1,n2);
	MATRIXFF(memcpy)(m,&mv.matrix);
}

static inline size_t strip_diag(void* restrict data,size_t n1,size_t n2,size_t s)
{
	char*	d=data;
	size_t	i,end,nc;
	
	nc=GSL_MIN(n1,n2);
	end=n1==n2?n1*n2-1:n1*n2;
	for(i=0;i<nc;i++)
	{
		if(i*(n2+1)>=end)
			break;
		memcpy(d+i*(n2+1)*s,d+(--end)*s,s);
		if((end-1)/n2==(end-1)%n2)
			end--;
	}
	return end;
}

static inline void VECTOROF(bound_below)(VECTORO* m,float v)
{
	size_t	i;
	float	f;
	for(i=0;i<m->size;i++)
	{
		f=VECTOROF(get)(m,i);
		if(f<v)
			VECTOROF(set)(m,i,v);
	}
}

static inline void VECTOROF(bound_above)(VECTORO* m,float v)
{
	size_t	i;
	float	f;
	for(i=0;i<m->size;i++)
	{
		f=VECTOROF(get)(m,i);
		if(f>v)
			VECTOROF(set)(m,i,v);
	}
}

static inline void VECTOROF(bound_both)(VECTORO* m,float vlow,float vhigh)
{
	size_t	i;
	float	f;
	for(i=0;i<m->size;i++)
	{
		f=VECTOROF(get)(m,i);
		if(f<vlow)
			VECTOROF(set)(m,i,vlow);
		else if(f>vhigh)
			VECTOROF(set)(m,i,vhigh);
	}
}

static inline void MATRIXOF(bound_below)(MATRIXO* m,float v)
{
	size_t	i,j;
	float	f;
	for(i=0;i<m->size1;i++)
		for(j=0;j<m->size2;j++)
		{
			f=MATRIXOF(get)(m,i,j);
			if(f<v)
				MATRIXOF(set)(m,i,j,v);
		}
}

static inline void MATRIXOF(bound_above)(MATRIXO* m,float v)
{
	size_t	i,j;
	float	f;
	for(i=0;i<m->size1;i++)
		for(j=0;j<m->size2;j++)
		{
			f=MATRIXOF(get)(m,i,j);
			if(f>v)
				MATRIXOF(set)(m,i,j,v);
		}
}

static inline void MATRIXOF(bound_both)(MATRIXO* m,float vlow,float vhigh)
{
	size_t	i,j;
	float	f;
	for(i=0;i<m->size1;i++)
		for(j=0;j<m->size2;j++)
		{
			f=MATRIXOF(get)(m,i,j);
			if(f<vlow)
				MATRIXOF(set)(m,i,j,vlow);
			else if(f>vhigh)
				MATRIXOF(set)(m,i,j,vhigh);
		}
}

static inline void VECTORDF(bound_below)(VECTORD* m,double v)
{
	size_t	i;
	double	f;
	for(i=0;i<m->size;i++)
	{
		f=VECTORDF(get)(m,i);
		if(f<v)
			VECTORDF(set)(m,i,v);
	}
}

static inline void VECTORDF(bound_above)(VECTORD* m,double v)
{
	size_t	i;
	double	f;
	for(i=0;i<m->size;i++)
	{
		f=VECTORDF(get)(m,i);
		if(f>v)
			VECTORDF(set)(m,i,v);
	}
}

static inline void VECTORDF(bound_both)(VECTORD* m,double vlow,double vhigh)
{
	size_t	i;
	double	f;
	for(i=0;i<m->size;i++)
	{
		f=VECTORDF(get)(m,i);
		if(f<vlow)
			VECTORDF(set)(m,i,vlow);
		else if(f>vhigh)
			VECTORDF(set)(m,i,vhigh);
	}
}

static inline void MATRIXDF(bound_below)(MATRIXD* m,double v)
{
	size_t	i,j;
	double	f;
	for(i=0;i<m->size1;i++)
		for(j=0;j<m->size2;j++)
		{
			f=MATRIXDF(get)(m,i,j);
			if(f<v)
				MATRIXDF(set)(m,i,j,v);
		}
}

static inline void MATRIXDF(bound_above)(MATRIXD* m,double v)
{
	size_t	i,j;
	double	f;
	for(i=0;i<m->size1;i++)
		for(j=0;j<m->size2;j++)
		{
			f=MATRIXDF(get)(m,i,j);
			if(f>v)
				MATRIXDF(set)(m,i,j,v);
		}
}

static inline void MATRIXDF(bound_both)(MATRIXD* m,double vlow,double vhigh)
{
	size_t	i,j;
	double	f;
	for(i=0;i<m->size1;i++)
		for(j=0;j<m->size2;j++)
		{
			f=MATRIXDF(get)(m,i,j);
			if(f<vlow)
				MATRIXDF(set)(m,i,j,vlow);
			else if(f>vhigh)
				MATRIXDF(set)(m,i,j,vhigh);
		}
}

static inline void VECTOROF(set_cond)(VECTORO* v,int (*func)(float),float val)
{
	size_t	i;
	for(i=0;i<v->size;i++)
		if(func(VECTOROF(get)(v,i)))
			VECTOROF(set)(v,i,val);
}

static inline void VECTORDF(set_cond)(VECTORD* v,int (*func)(double),double val)
{
	size_t	i;
	for(i=0;i<v->size;i++)
		if(func(VECTORDF(get)(v,i)))
			VECTORDF(set)(v,i,val);
}

static inline void MATRIXOF(set_cond)(MATRIXO* m,int (*func)(float),float val)
{
	size_t	i,j;
	for(i=0;i<m->size1;i++)
		for(j=0;j<m->size2;j++)
			if(func(MATRIXOF(get)(m,i,j)))
				MATRIXOF(set)(m,i,j,val);
}

static inline void MATRIXDF(set_cond)(MATRIXD* m,int (*func)(double),double val)
{
	size_t	i,j;
	for(i=0;i<m->size1;i++)
		for(j=0;j<m->size2;j++)
			if(func(MATRIXDF(get)(m,i,j)))
				MATRIXDF(set)(m,i,j,val);
}

static inline int gsl_isinf_f(float f)
{
	return gsl_isinf(f);
}

static inline int gsl_isnan_f(float f)
{
	return gsl_isnan(f);
}

static inline void VECTOROF(set_inf)(VECTORO* v,float val)
{
	VECTOROF(set_cond)(v,gsl_isinf_f,val);
}

static inline void VECTOROF(set_nan)(VECTORO* v,float val)
{
	VECTOROF(set_cond)(v,gsl_isnan_f,val);
}

static inline void VECTORDF(set_inf)(VECTORD* v,double val)
{
	VECTORDF(set_cond)(v,gsl_isinf,val);
}

static inline void VECTORDF(set_nan)(VECTORD* v,double val)
{
	VECTORDF(set_cond)(v,gsl_isnan,val);
}

static inline void MATRIXOF(set_inf)(MATRIXO* m,float val)
{
	MATRIXOF(set_cond)(m,gsl_isinf_f,val);
}

static inline void MATRIXOF(set_nan)(MATRIXO* m,float val)
{
	MATRIXOF(set_cond)(m,gsl_isnan_f,val);
}

static inline void MATRIXDF(set_inf)(MATRIXD* m,double val)
{
	MATRIXDF(set_cond)(m,gsl_isinf,val);
}

static inline void MATRIXDF(set_nan)(MATRIXD* m,double val)
{
	MATRIXDF(set_cond)(m,gsl_isnan,val);
}

static inline void MATRIXFF(cov1)(const MATRIXF* m,MATRIXF* cov)
{
	size_t i,j;
	
	BLASF(syrk)(CblasUpper,CblasNoTrans,(FTYPE)1./((FTYPE)(m->size2)),m,0,cov);
	//Set bottom half
	for(i=m->size1-1;i;i--)
		for(j=0;j<i;j++)
			MATRIXFF(set)(cov,i,j,MATRIXFF(get)(cov,j,i));
}

static inline void MATRIXFF(cov2)(const MATRIXF* m1,const MATRIXF* m2,MATRIXF* cov)
{
	BLASF(gemm)(CblasNoTrans,CblasTrans,(FTYPE)1./((FTYPE)(m1->size2)),m1,m2,0,cov);
}

static inline void MATRIXFF(cov2_1v1)(const MATRIXF* m1,const MATRIXF* m2,VECTORF* cov)
{
	size_t	i,n;
	n=m1->size1;
	assert((m2->size1==n)&&(cov->size==n)&&(m1->size2==m2->size2));
	
	for(i=0;i<n;i++)
	{
		VECTORFF(const_view) vv1=MATRIXFF(const_row)(m1,i);
		VECTORFF(const_view) vv2=MATRIXFF(const_row)(m2,i);
		BLASF(dot)(&vv1.vector,&vv2.vector,VECTORFF(ptr)(cov,i));
	}
	VECTORFF(scale)(cov,1/((FTYPE)(m1->size2)));
}

static inline void MATRIXFF(cov1_bounded)(const MATRIXF* m,MATRIXF* cov)
{
	MATRIXFF(cov1)(m,cov);
	MATRIXFF(bound_both)(cov,-1,1);
}

static inline void MATRIXFF(cov2_bounded)(const MATRIXF* m1,const MATRIXF* m2,MATRIXF* cov)
{
	MATRIXFF(cov2)(m1,m2,cov);
	MATRIXFF(bound_both)(cov,-1,1);
}

static inline void MATRIXFF(cov2_1v1_bounded)(const MATRIXF* m1,const MATRIXF* m2,VECTORF* cov)
{
	MATRIXFF(cov2_1v1)(m1,m2,cov);
	VECTORFF(bound_both)(cov,-1,1);
}

static inline long VECTORFF(first_nan)(const VECTORF *restrict v)
{
	size_t i;
	for(i=0;i<v->size;i++)
		if(gsl_isnan(VECTORFF(get)(v,i)))
			return (int)i;
	return -1;
}

static inline long MATRIXFF(first_nan)(const MATRIXF *restrict m)
{
	size_t i,j;
	for(i=0;i<m->size1;i++)
		for(j=0;j<m->size2;j++)
			if(gsl_isnan(MATRIXFF(get)(m,i,j)))
				return (long)(i*m->size2+j);
	return -1;
}

static inline void MATRIXFF(scale_row)(MATRIXF* m,const VECTORF* v)
{
	size_t	i;
	VECTORFF(view) vv;
	
	assert(m->size1==v->size);
	for(i=0;i<m->size1;i++)
	{
		vv=MATRIXFF(row)(m,i);
		VECTORFF(scale)(&vv.vector,VECTORFF(get)(v,i));
	}
}

static inline void MATRIXFF(scale_row_func)(MATRIXF* m,const VECTORF* v,FTYPE (*f)(FTYPE))
{
	size_t	i;
	VECTORFF(view) vv;
	
	assert(m->size1==v->size);
	for(i=0;i<m->size1;i++)
	{
		vv=MATRIXFF(row)(m,i);
		VECTORFF(scale)(&vv.vector,f(VECTORFF(get)(v,i)));
	}
}


static inline void MATRIXFF(scale_column)(MATRIXF* m,const VECTORF* v)
{
	size_t	i;
	VECTORFF(view) vv;
	
	assert(m->size2==v->size);
	for(i=0;i<m->size2;i++)
	{
		vv=MATRIXFF(column)(m,i);
		VECTORFF(scale)(&vv.vector,VECTORFF(get)(v,i));
	}
}

static inline void MATRIXFF(scale_column_func)(MATRIXF* m,const VECTORF* v,FTYPE (*f)(FTYPE))
{
	size_t	i;
	VECTORFF(view) vv;
	
	assert(m->size2==v->size);
	for(i=0;i<m->size2;i++)
	{
		vv=MATRIXFF(column)(m,i);
		VECTORFF(scale)(&vv.vector,f(VECTORFF(get)(v,i)));
	}
}

static inline void VECTORGF(eq)(const VECTORG *d,VECTORUC* ans,GTYPE val)
{
	size_t i;
	
	assert(d->size==ans->size);
	for(i=0;i<d->size;i++)
		VECTORUCF(set)(ans,i,VECTORGF(get)(d,i)==val);
}

static inline void VECTORDF(diff)(VECTORD* d)
{
	size_t	i;
	for(i=d->size-1;i;i--)
		VECTORDF(ptr)(d,i)[0]-=VECTORDF(get)(d,i-1);
}

static inline void VECTOROF(diff)(VECTORF* d)
{
	size_t	i;
	for(i=d->size-1;i;i--)
		VECTORFF(ptr)(d,i)[0]-=VECTORFF(get)(d,i-1);
}

static inline void VECTORDF(cumsum)(VECTORD* d)
{
	size_t	i;
	for(i=1;i<d->size;i++)
		VECTORDF(ptr)(d,i)[0]+=VECTORDF(get)(d,i-1);
}

static inline void VECTOROF(cumsum)(VECTORF* d)
{
	size_t	i;
	for(i=1;i<d->size;i++)
		VECTORFF(ptr)(d,i)[0]+=VECTORFF(get)(d,i-1);
}

static inline void MATRIXFF(fluc)(MATRIXF* m,FTYPE fluc)
{
	size_t	i,j;
	for(i=0;i<m->size1;i++)
		for(j=0;j<m->size2;j++)
			MATRIXFF(ptr)(m,i,j)[0]*=(FTYPE)(1+(2*random_uniform()-1)*fluc);
}

static inline void VECTORFF(fluc)(VECTORF* v,FTYPE fluc)
{
	size_t	i;
	for(i=0;i<v->size;i++)
		VECTORFF(ptr)(v,i)[0]*=(FTYPE)(1+(2*random_uniform()-1)*fluc);
}












































#ifdef __cplusplus
}
#endif
#endif
