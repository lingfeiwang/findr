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
#include "../../base/config.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../../base/logger.h"
#include "../../base/threading.h"
#include "../../base/macros.h"
#include "../llrtopv.h"

/* Convert log likelihood ratios into p-values for matrix
 * d:		(ng,nt)	Data, as input for log likelihood ratios,
 			and also as output for converted p-values.
 * g:		(ng,ns)	Original genotype matrix, used for analytical calculation
 * 			of null distribution. Every element=0,1,...,nv-1.
 			Has matching rows with data d.
 * nv:		Maximum number of values each g may take.
 * n1c,
 * n1d,
 * n2c,
 * n2d:		Parameters to specify null distribution. See pij_llrtopij_a_nullhist
 * Return:	0 if success.
 */
static int pij_gassist_llrtopv_block(MATRIXF* d,const MATRIXG* g,size_t nv,long n1c,size_t n1d,long n2c,size_t n2d)
{
#define	CLEANUP	AUTOFREE(nexist)

	size_t		i,j,nvr;
	VECTORFF(view)	vv;
	assert(d->size1==g->size1);
	assert(MATRIXGF(max)(g)<nv);
	
	AUTOALLOC(unsigned char,nexist,nv,16)
	if(!nexist)
		ERRRET("Not enough memory.");
	
	for(i=0;i<d->size1;i++)
	{
		//Count for number of genotypes
		memset(nexist,0,nv*sizeof(nexist[0]));
		for(j=0;j<g->size2;j++)
			nexist[MATRIXGF(get)(g,i,j)]=1;
		nvr=0;
		for(j=0;j<nv;j++)
			nvr+=nexist[j];
		assert(nvr>1);
		nvr-=2;
		assert(((long)nvr*n1c+(long)n1d>0)&&((long)n2d>(long)nvr*n2c));
		vv=MATRIXFF(row)(d,i);
		pij_llrtopv_block(&vv.vector,(size_t)((long)nvr*n1c+(long)n1d),(size_t)((long)n2d-(long)nvr*n2c));
	}
	
	CLEANUP
	return 0;
#undef	CLEANUP
}

/* Convert log likelihood ratios into p-values for vector
 * d:		(ng)	Data, as input for log likelihood ratios,
 			and also as output for converted p-values.
 * g:		(ng,ns)	Original genotype matrix, used for analytical calculation
 * 			of null distribution. Every element=0,1,...,nv-1.
 			Has matching rows with data d.
 * nv:		Maximum number of values each g may take.
 * n1c,
 * n1d,
 * n2c,
 * n2d:		Parameters to specify null distribution. See pij_llrtopij_a_nullhist
 * Return:	0 if success.
 */
static inline int pij_gassist_llrtopv_vec_block(VECTORF* d,const MATRIXG* g,size_t nv,long n1c,size_t n1d,long n2c,size_t n2d)
{
	MATRIXFF(view)	mv=MATRIXFF(view_vector)(d,d->size,1);
	assert(d->size==g->size1);
	assert((long)n2d>n2c*(long)nv);
	return pij_gassist_llrtopv_block(&mv.matrix,g,nv,n1c,n1d,n2c,n2d);
}

/* Convert log likelihood ratios into p-values for matrix in single thread
 * p1:		(ng)
 * p2:		(ng,nt)
 * p3:		(ng,nt)
 * p4:		(ng,nt)
 * p5:		(ng,nt)	Data for 5 tests from p1 to p5, as input for log likelihood ratios,
 			and also as output for converted p-values.
 * g:		(ng,ns)	Original genotype matrix, used for analytical calculation
 * 			of null distribution. Every element=0,1,...,nv-1.
 			Has matching rows with data p1 to p5..
 * nv:		Maximum number of values each g may take.
 * Return:	0 if success.
 */
static int pij_gassist_llrtopvs_block(VECTORF* p1,MATRIXF* p2,MATRIXF* p3,MATRIXF* p4,MATRIXF* p5,const MATRIXG* g,size_t nv)
{
	int	ret=0,ret2=0;
	assert((p1->size==g->size1)&&(p2->size1==g->size1)&&(p3->size1==g->size1)&&(p4->size1==g->size1)&&(p5->size1==g->size1));
	assert((p2->size2==p3->size2)&&(p2->size2==p4->size2)&&(p2->size2==p5->size2));

	ret=ret||(ret2=pij_gassist_llrtopv_vec_block(p1,g,nv,1,1,1,g->size2-2));
	if(ret2)
		LOG(1,"Failed to log likelihood ratios to p-values in step 1.")
	ret=ret||(ret2=pij_gassist_llrtopv_block(p2,g,nv,1,1,1,g->size2-2));
	if(ret2)
		LOG(1,"Failed to log likelihood ratios to p-values in step 2.")
	ret=ret||(ret2=pij_gassist_llrtopv_block(p3,g,nv,1,1,1,g->size2-3));
	if(ret2)
		LOG(1,"Failed to log likelihood ratios to p-values in step 3.")
	ret=ret||(ret2=pij_gassist_llrtopv_block(p4,g,nv,1,2,1,g->size2-3));
	if(ret2)
		LOG(1,"Failed to log likelihood ratios to p-values in step 4.")
	ret=ret||(ret2=pij_gassist_llrtopv_block(p5,g,nv,0,1,1,g->size2-3));
	if(ret2)
		LOG(1,"Failed to log likelihood ratios to p-values in step 5.")
	return ret;
}

int pij_gassist_llrtopvs(VECTORF* p1,MATRIXF* p2,MATRIXF* p3,MATRIXF* p4,MATRIXF* p5,const MATRIXG* g,size_t nv)
{
	int	ret=0;
	assert((p1->size==g->size1)&&(p2->size1==g->size1)&&(p3->size1==g->size1)&&(p4->size1==g->size1)&&(p5->size1==g->size1));
	assert((p2->size2==p3->size2)&&(p2->size2==p4->size2)&&(p2->size2==p5->size2));
	
	if(g->size2<nv+1)
	{
		LOG(0,"Needs sample size at least number of alleles + 2 to compute p-values.")
		return 1;
	}
	
	#pragma omp parallel
	{
		size_t	ng1,ng2;
		int		ret2=0;

		threading_get_startend(p1->size,&ng1,&ng2);
		if(ng2>ng1)
		{
			size_t			dn=ng2-ng1;
			VECTORFF(view)	vv1;
			MATRIXFF(view)	mv2,mv3,mv4,mv5;
			MATRIXGF(const_view)	mvg=MATRIXGF(const_submatrix)(g,ng1,0,dn,g->size2);
			vv1=VECTORFF(subvector)(p1,ng1,dn);
			mv2=MATRIXFF(submatrix)(p2,ng1,0,dn,p2->size2);
			mv3=MATRIXFF(submatrix)(p3,ng1,0,dn,p2->size2);
			mv4=MATRIXFF(submatrix)(p4,ng1,0,dn,p2->size2);
			mv5=MATRIXFF(submatrix)(p5,ng1,0,dn,p2->size2);
			ret2=pij_gassist_llrtopvs_block(&vv1.vector,&mv2.matrix,&mv3.matrix,&mv4.matrix,&mv5.matrix,&mvg.matrix,nv);
		}
		#pragma omp critical
			ret=ret||ret2;
	}
	return ret;
}
















