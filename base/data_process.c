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
#include <stdio.h>
#include <math.h>
#include <tgmath.h>
#include "types.h"
#include "gsl/sort.h"
#include "logger.h"
#include "macros.h"
#include "data_process.h"
#pragma GCC diagnostic ignored "-Wconversion"

void MATRIXFF(from_dense)(MATRIXF* dest,const FTYPE* restrict data,size_t nrow,size_t ncol)
{
	MATRIXFF(view) mview=MATRIXFF(view_array)((FTYPE*)data,nrow,ncol);
	memcpy(dest,&(mview.matrix),sizeof(*dest));
}

MATRIXF* MATRIXFF(from_densefile)(FILE* f,size_t nrow,size_t ncol)
{
	MATRIXF* m;
	int		ret;
	m=MATRIXFF(alloc)(nrow,ncol);
	ret=MATRIXFF(fread)(f,m);
	if(ret)
	{
		LOG(3,"Failed to read data matrix file.")
		MATRIXFF(free)(m);
		return 0;
	}
	return m;
}

void MATRIXFF(normalize_row_buffed)(MATRIXF* m,VECTORF* v1,const VECTORF* v2)
{
	VECTORFF(view)	vv;
	size_t	i;
	FTYPE	result;
	
	//Calculate mean
	//VECTORFF(set_all)(v2,1);
	BLASF(gemv)(CblasNoTrans,1,m,v2,0,v1);
	VECTORFF(scale)(v1,1./(double)m->size2);
	for(i=0;i<m->size1;i++)
	{
		//Subtract mean
		vv=MATRIXFF(row)(m,i);
		VECTORFF(add_constant)(&(vv.vector),-VECTORFF(get)(v1,i));
		//Calculate variance
		BLASF(dot)(&(vv.vector),&(vv.vector),&result);
		result=1/((float)sqrt(result/(float)m->size2));
		//Scale variance
		VECTORFF(scale)(&(vv.vector),result);
	}
}

int MATRIXFF(normalize_row)(MATRIXF* m)
{
#define	CLEANUP	AUTOFREEVEC(v1)AUTOFREEVEC(v2)
	AUTOALLOCVECF(v1,m->size1,10000)
	AUTOALLOCVECF(v2,m->size2,10000)
	
	if(!(v1&&v2))
	{
		LOG(2,"Failed to allocate memory.")
		CLEANUP
		return -1;
	}
	VECTORFF(set_zero)(v1);
	VECTORFF(set_all)(v2,1);
	MATRIXFF(normalize_row_buffed)(m,v1,v2);
	CLEANUP
	return 0;
}

void MATRIXFF(flatten_nodiag)(const MATRIXF* m,VECTORF* v)
{
	size_t	i,n1,n2,nc;
	VECTORFF(view)	vv2,vv4;
	MATRIXFF(view)	mv1;
	
	//Initialize
	n1=m->size1;
	n2=m->size2;
	nc=GSL_MIN(n1,n2);
	assert(v->size==n1*n2-nc);
	
	{
		VECTORFF(const_view) vv1=MATRIXFF(const_subrow)(m,0,1,n2-1);
		vv2=VECTORFF(subvector)(v,0,n2-1);
		VECTORFF(memcpy)(&vv2.vector,&vv1.vector);
	}
	for(i=1;i<nc-1;i++)
	{
		VECTORFF(const_view) vv1=MATRIXFF(const_subrow)(m,i,0,i);
		vv2=VECTORFF(subvector)(v,i*(n2-1),i);
		VECTORFF(memcpy)(&vv2.vector,&vv1.vector);
		VECTORFF(const_view) vv3=MATRIXFF(const_subrow)(m,i,i+1,n2-i-1);
		vv2=VECTORFF(subvector)(v,i*(n2-1)+i,n2-i-1);
		VECTORFF(memcpy)(&vv2.vector,&vv3.vector);
	}
	{
		VECTORFF(const_view) vv1=MATRIXFF(const_subrow)(m,nc-1,0,nc-1);
		vv2=VECTORFF(subvector)(v,(nc-1)*(n2-1),nc-1);
		VECTORFF(memcpy)(&vv2.vector,&vv1.vector);
	}
	if(nc==n1)
	{
		if(n1!=n2)
		{
			VECTORFF(const_view) vv1=MATRIXFF(const_subrow)(m,nc-1,nc,n2-nc);
			vv2=VECTORFF(subvector)(v,(nc-1)*n2,n2-nc);
			VECTORFF(memcpy)(&vv2.vector,&vv1.vector);
		}
		return;
	}
	vv4=VECTORFF(subvector)(v,(n2-1)*n2,(n1-n2)*n2);
	mv1=MATRIXFF(view_vector)(&vv4.vector,n1-n2,n2);
	MATRIXFF(const_view) mv2=MATRIXFF(const_submatrix)(m,n2,0,n1-n2,n2);
	MATRIXFF(memcpy)(&mv1.matrix,&mv2.matrix);
}

void VECTORFF(wrap_nodiag)(const VECTORF* v,MATRIXF* m)
{
	size_t	i,n1,n2,nc;
	VECTORFF(view)	vv1;
	MATRIXFF(view)	mv1;
	
	//Initialize
	n1=m->size1;
	n2=m->size2;
	nc=GSL_MIN(n1,n2);
	assert(v->size==n1*n2-nc);
	
	{
		vv1=MATRIXFF(subrow)(m,0,1,n2-1);
		VECTORFF(const_view) cvv2=VECTORFF(const_subvector)(v,0,n2-1);
		VECTORFF(memcpy)(&vv1.vector,&cvv2.vector);
	}
	for(i=1;i<nc-1;i++)
	{
		vv1=MATRIXFF(subrow)(m,i,0,i);
		VECTORFF(const_view) cvv2=VECTORFF(const_subvector)(v,i*(n2-1),i);
		VECTORFF(memcpy)(&vv1.vector,&cvv2.vector);
		vv1=MATRIXFF(subrow)(m,i,i+1,n2-i-1);
		VECTORFF(const_view) cvv4=VECTORFF(const_subvector)(v,i*(n2-1)+i,n2-i-1);
		VECTORFF(memcpy)(&vv1.vector,&cvv4.vector);
	}
	{
		vv1=MATRIXFF(subrow)(m,nc-1,0,nc-1);
		VECTORFF(const_view) cvv2=VECTORFF(const_subvector)(v,(nc-1)*(n2-1),nc-1);
		VECTORFF(memcpy)(&vv1.vector,&cvv2.vector);
	}
	if(nc==n1)
	{
		if(n1!=n2)
		{
			vv1=MATRIXFF(subrow)(m,nc-1,nc,n2-nc);
			VECTORFF(const_view) cvv2=VECTORFF(const_subvector)(v,(nc-1)*n2,n2-nc);
			VECTORFF(memcpy)(&vv1.vector,&cvv2.vector);
		}
		return;
	}
	VECTORFF(const_view) cvv4=VECTORFF(const_subvector)(v,(n2-1)*n2,(n1-n2)*n2);
	MATRIXFF(const_view) cmv1=MATRIXFF(const_view_vector)(&cvv4.vector,n1-n2,n2);
	mv1=MATRIXFF(submatrix)(m,n2,0,n1-n2,n2);
	MATRIXFF(memcpy)(&mv1.matrix,&cmv1.matrix);
}

void MATRIXGF(countv_byrow_buffed)(const MATRIXG* g,VECTORG* ans,VECTORUC* vb)
{
	size_t	ng=g->size1;
	size_t	ns=g->size2;
	size_t	i,j;
	unsigned char * restrict p;
	
	VECTORGF(set_zero)(ans);
	for(i=0;i<ng;i++)
	{
		VECTORUCF(set_zero)(vb);
		p=VECTORGF(ptr)(ans,i);
		for(j=0;j<ns;j++)
			VECTORUCF(set)(vb,MATRIXGF(get)(g,i,j),1);
		for(j=0;j<vb->size;j++)
			p[0]+=VECTORUCF(get)(vb,j);
	}
}

void VECTORGF(count_ratio_d)(const VECTORG* d,VECTORD* ans)
{
	size_t	i;
	
	VECTORDF(set_zero)(ans);
	for(i=0;i<d->size;i++)
		VECTORDF(ptr)(ans,(size_t)VECTORGF(get)(d,i))[0]+=(double)1.;
	VECTORDF(scale)(ans,(double)1./(double)(d->size));
}

size_t	MATRIXGF(rows_save)(const MATRIXG* d,MATRIXG* dest,const VECTORUC* c)
{
	size_t	n,i;
	VECTORGF(view)	vv;
	assert((d->size2==dest->size2)&&(d->size1==c->size));
	for(i=0,n=0;i<d->size1;i++)
		if(VECTORUCF(get)(c,i))
		{
			vv=MATRIXGF(row)(dest,n++);
			MATRIXGF(get_row)(&vv.vector,d,i);
		}
	return n;
}

size_t	MATRIXFF(rows_save)(const MATRIXF* d,MATRIXF* dest,const VECTORUC* c)
{
	size_t	n,i;
	VECTORFF(view)	vv;
	assert((d->size2==dest->size2)&&(d->size1==c->size));
	for(i=0,n=0;i<d->size1;i++)
		if(VECTORUCF(get)(c,i))
		{
			vv=MATRIXFF(row)(dest,n++);
			MATRIXFF(get_row)(&vv.vector,d,i);
		}
	return n;
}

size_t	MATRIXFF(rows_save_nodiag)(const MATRIXF* d,VECTORF* dest,const VECTORUC* c)
{
	size_t	n,i,ng,ns;
	VECTORFF(view)	vv;
	assert((d->size1==c->size));
	
	ng=d->size1;
	ns=d->size2;
	for(i=0,n=0;i<ng;i++)
		if(VECTORUCF(get)(c,i))
		{
			if(i<ns)
			{
				if(i)
				{
					vv=VECTORFF(subvector)(dest,n,i);
					VECTORFF(const_view) vvc=MATRIXFF(const_subrow)(d,i,0,i);
					VECTORFF(memcpy)(&vv.vector,&vvc.vector);
					n+=i;
				}
				if(i<ns-1)
				{
					vv=VECTORFF(subvector)(dest,n,ns-i-1);
					VECTORFF(const_view) vvc=MATRIXFF(const_subrow)(d,i,i+1,ns-i-1);
					VECTORFF(memcpy)(&vv.vector,&vvc.vector);
					n+=ns-i-1;
				}
			}
			else
			{
				vv=VECTORFF(subvector)(dest,n,ns);
				MATRIXFF(get_row)(&vv.vector,d,i);
				n+=ns;
			}
		}
	return n;
}


size_t	MATRIXGF(rows_load)(const MATRIXG* d,MATRIXG* dest,const VECTORUC* c)
{
	size_t	n,i;
	VECTORGF(view)	vv;
	assert((d->size2==dest->size2)&&(dest->size1==c->size));
	for(i=0,n=0;i<dest->size1;i++)
		if(VECTORUCF(get)(c,i))
		{
			vv=MATRIXGF(row)(dest,i);
			MATRIXGF(get_row)(&vv.vector,d,n++);
		}
	return n;
}

size_t	MATRIXFF(rows_load)(const MATRIXF* d,MATRIXF* dest,const VECTORUC* c)
{
	size_t	n,i;
	VECTORFF(view)	vv;
	assert((d->size2==dest->size2)&&(dest->size1==c->size));
	for(i=0,n=0;i<dest->size1;i++)
		if(VECTORUCF(get)(c,i))
		{
			vv=MATRIXFF(row)(dest,i);
			MATRIXFF(get_row)(&vv.vector,d,n++);
		}
	return n;
}

size_t	MATRIXFF(rows_load_nodiag)(const VECTORF* d,MATRIXF* dest,const VECTORUC* c)
{
	size_t	n,i,ng,ns;
	VECTORFF(view)	vv;
	assert((dest->size1==c->size));
	
	ng=dest->size1;
	ns=dest->size2;
	for(i=0,n=0;i<ng;i++)
		if(VECTORUCF(get)(c,i))
		{
			if(i<ns)
			{
				if(i)
				{
					VECTORFF(const_view) vvc=VECTORFF(const_subvector)(d,n,i);
					vv=MATRIXFF(subrow)(dest,i,0,i);
					VECTORFF(memcpy)(&vv.vector,&vvc.vector);
					n+=i;
				}
				if(i<ns-1)
				{
					VECTORFF(const_view) vvc=VECTORFF(const_subvector)(d,n,ns-i-1);
					vv=MATRIXFF(subrow)(dest,i,i+1,ns-i-1);
					VECTORFF(memcpy)(&vv.vector,&vvc.vector);
					n+=ns-i-1;
				}
			}
			else
			{
				VECTORFF(const_view) vvc=VECTORFF(const_subvector)(d,n,ns);
				MATRIXFF(set_row)(dest,i,&vvc.vector);
				n+=ns;
			}
		}
	return n;
}

void MATRIXFF(permute_column_buffed)(MATRIXF* m,const gsl_permutation* perm,VECTORF* b1,VECTORUC* b2)
{
	size_t	i,j,k;
	const size_t	*p=perm->data;
	VECTORFF(view)	vv1,vv2;
	
	VECTORUCF(set_zero)(b2);
	for(i=0;i<m->size2;i++)
		if((!VECTORUCF(get)(b2,i))&&((k=p[i])!=i))
		{
			VECTORUCF(set)(b2,i,1);
			j=i;
			vv1=MATRIXFF(column)(m,j);
			VECTORFF(memcpy)(b1,&vv1.vector);
			do
			{
				VECTORUCF(set)(b2,k,1);
				vv2=MATRIXFF(column)(m,k);
				VECTORFF(memcpy)(&vv1.vector,&vv2.vector);
				vv1=vv2;
				j=k;
				k=p[j];
				
			}while(k!=i);
			VECTORFF(memcpy)(&vv1.vector,b1);
		}
}

void MATRIXFF(permute_row_buffed)(MATRIXF* m,const gsl_permutation* perm,VECTORF* b1,VECTORUC* b2)
{
	size_t	i,j,k;
	const size_t	*p=perm->data;
	VECTORFF(view)	vv1,vv2;
	
	VECTORUCF(set_zero)(b2);
	for(i=0;i<m->size1;i++)
		if((!VECTORUCF(get)(b2,i))&&((k=p[i])!=i))
		{
			VECTORUCF(set)(b2,i,1);
			j=i;
			vv1=MATRIXFF(row)(m,j);
			VECTORFF(memcpy)(b1,&vv1.vector);
			do
			{
				VECTORUCF(set)(b2,k,1);
				vv2=MATRIXFF(row)(m,k);
				VECTORFF(memcpy)(&vv1.vector,&vv2.vector);
				vv1=vv2;
				j=k;
				k=p[j];			
			}while(k!=i);
			VECTORFF(memcpy)(&vv1.vector,b1);
		}
}

































