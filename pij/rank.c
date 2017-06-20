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
#include "../base/config.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include "../base/gsl/blas.h"
#include "../base/random.h"
#include "../base/const.h"
#include "../base/logger.h"
#include "../base/macros.h"
#include "../base/data_process.h"
#include "../base/supernormalize.h"
#include "../base/threading.h"
#include "llrtopij_a.h"
#include "llrtopv.h"
#include "rank.h"

/* Calculates the log likelihood ratio correlated v.s. uncorrelated models.
 * Uses GSL BLAS.
 */
static void pij_rank_llr_block(const MATRIXF* t,const MATRIXF* t2,MATRIXF* llr)
{
	size_t	i,j;
	size_t	ng=t->size1;
	size_t	nt=t2->size1;
#ifndef NDEBUG
	size_t	ns=t->size2;
#endif
	assert((t2->size2==ns));
	assert((llr->size1==ng));
	assert((llr->size2==nt));
	MATRIXFF(cov2_bounded)(t,t2,llr);
	MATRIXFF(mul_elements)(llr,llr);
	MATRIXFF(scale)(llr,-1);
	MATRIXFF(add_constant)(llr,1);
	for(i=0;i<ng;i++)
		for(j=0;j<nt;j++)
			MATRIXFF(set)(llr,i,j,(FTYPE)log(MATRIXFF(get)(llr,i,j)));

	MATRIXFF(scale)(llr,-0.5);
	//Bounding from 0
	MATRIXFF(bound_below)(llr,0);
}

/* Multithread calculation of log likelihood ratio
 * t:	(ng,ns) Full transcript data matrix of A
 * t2:	(nt,ns) Full transcript data matrix of B
 * llr:	(ng,nt). Log likelihood ratios for test.
 */
static void pij_rank_llr(const MATRIXF* t,const MATRIXF* t2,MATRIXF* llr)
{
	assert((t->size2==t2->size2)&&(llr->size1==t->size1)&&(llr->size2==t2->size1));
	#pragma omp parallel
	{
		size_t	n1,n2;
		
		threading_get_startend(t->size1,&n1,&n2);
		if(n2>n1)
		{
			MATRIXFF(const_view) mvt=MATRIXFF(const_submatrix)(t,n1,0,n2-n1,t->size2);
			MATRIXFF(view)	mvllr;
			mvllr=MATRIXFF(submatrix)(llr,n1,0,n2-n1,llr->size2);
			pij_rank_llr_block(&mvt.matrix,t2,&mvllr.matrix);
		}
	}
}

/* Converts log likelihood ratios into p-values for ranked correlation test
 * d:	MATRIXF of any size, as input of LLRs and also output of corresponding p-values
 * ns:	Number of samples, to be used to calculate the null distribution
 */
static inline void pij_rank_llrtopv(MATRIXF* d,size_t ns)
{
	assert(ns>2);
	pij_llrtopvm(d,1,ns-2);
}

int pij_rank_pv(const MATRIXF* t,const MATRIXF* t2,MATRIXF* p,size_t memlimit)
{
#define	CLEANUP		CLEANMATF(tnew)CLEANMATF(tnew2)
	MATRIXF		*tnew,*tnew2;			//(nt,ns) Supernormalized transcript matrix
	int			ret;
	size_t		ng,nt,ns;
	
	ng=t->size1;
	nt=t2->size1;
	ns=t->size2;

	tnew=tnew2=0;
	
	//Validation
	assert((t2->size2==ns)&&(p->size1==ng)&&(p->size2==nt)&&memlimit);
	if(ns<3)
		ERRRET("Needs at least 3 samples to compute p-values.")
	
	{
		size_t mem;
		mem=(2*t->size1*t->size2+2*t2->size1*t2->size2+p->size1*p->size2)*FTYPEBITS/8;
		if(memlimit<=mem)
			ERRRET("Memory limit lower than minimum memory needed. Try increasing your memory usage limit.")
		LOG(10,"Memory limit: %lu bytes.",memlimit)
	}

	tnew=MATRIXFF(alloc)(ng,ns);
	tnew2=MATRIXFF(alloc)(nt,ns);
	if(!(tnew&&tnew2))
		ERRRET("Not enough memory.")

	//Step 1: Supernormalization
	LOG(9,"Supernormalizing...")
	MATRIXFF(memcpy)(tnew,t);
	ret=supernormalizea_byrow(tnew);
	MATRIXFF(memcpy)(tnew2,t2);
	ret=ret||supernormalizea_byrow(tnew2);
	if(ret)
		ERRRET("Supernormalization failed.")

	//Step 2: Log likelihood ratios from nonpermuted data
	LOG(9,"Calculating real log likelihood ratios...")
	pij_rank_llr(tnew,tnew2,p);
	//Step 3: Convert log likelihood ratios to probabilities
	LOG(9,"Converting likelihood ratios into p-values...")
	pij_rank_llrtopv(p,ns);

	//Cleanup
	CLEANUP
	return 0;
#undef	CLEANUP		
}

/* Convert LLR into probabilities per A. Uses pij_llrtopij_a_convert.
 * ans:		(ng,nt) Source real LLRs to compare with null LLRs,
 * 			also output location of converted probabilities.
 * ns:		Number of samples, used for calculation of null distribution
 * nodiag:	Whether diagonal elements of d should be ignored when converting
 * 			to probabilities
 * nodiagshift:	Offdiagonal shift.
 * 				For nodiagshift>0/<0, use upper/lower diagonals of corresponding id.
 * Return:	0 if succeed.
 */
static int pij_rank_llrtopij_a(MATRIXF* ans,size_t ns,char nodiag,long nodiagshift)
{
	LOG(9,"Converting LLR to probabilities on per A basis.")
	if(ns<=2)
	{
		LOG(0,"Needs at least 3 samples to compute probabilities.")
		return 1;
	}
	return pij_llrtopij_a_convert_single_self(ans,1,ns-2,nodiag,nodiagshift);
}

/* Calculate probabilities of A--B based on LLR distributions of real data 
 * and null hypothesis.
 * t:		(ng,ns) Expression data for A
 * t2:		(nt,ns) Expression data for B
 * p:		(ng,nt) Output for probabilities A--B is true
 * nodiag:	When the top ng rows of t2 is exactly t, diagonals of pij are meaningless.
 *			In this case, set nodiag to 1 to avoid inclusion of NANs. For nodiag=0, t and t2
 *			should not have any identical genes.
 * pij:		Function to convert LLR to probabilities, such as pij_rank_llrtopij_a
 * memlimit:Specifies approximate memory usage. Function can fail if memlimit is too small. For large dataset, memory usage will be reduced by spliting t into smaller chunks and infer separately. For unlimited memory, set memlimit=-1.	
 * Return:	0 if succeed.
 */
static int pij_rank_any(const MATRIXF* t,const MATRIXF* t2,MATRIXF* p,char nodiag,int (*pij)(MATRIXF*,size_t,char,long),size_t memlimit)
{
#define	CLEANUP		CLEANMATF(tnew)CLEANMATF(tnew2)
	MATRIXF		*tnew,*tnew2;			//(nt,ns) Supernormalized transcript matrix
	VECTORFF(view)	vv;
	int			ret;
	size_t		ng,nt,ns;
	
	ng=t->size1;
	nt=t2->size1;
	ns=t->size2;

	tnew=tnew2=0;
	
	//Validation
	assert((t2->size2==ns)&&(p->size1==ng)&&(p->size2==nt)&&memlimit);
	
	if(ns<=2)
		ERRRET("Needs at least 3 samples to compute probabilities.")
	{
		size_t mem;
		mem=(2*t->size1*t->size2+2*t2->size1*t2->size2+p->size1*p->size2)*FTYPEBITS/8;
		if(memlimit<=mem)
			ERRRET("Memory limit lower than minimum memory needed. Try increasing your memory usage limit.")
		LOG(10,"Memory limit: %lu bytes.",memlimit)
	}

	tnew=MATRIXFF(alloc)(ng,ns);
	tnew2=MATRIXFF(alloc)(nt,ns);
	if(!(tnew&&tnew2))
		ERRRET("Not enough memory.")

	//Step 1: Supernormalization
	LOG(9,"Supernormalizing...")
	MATRIXFF(memcpy)(tnew,t);
	ret=supernormalizea_byrow(tnew);
	MATRIXFF(memcpy)(tnew2,t2);
	ret=ret||supernormalizea_byrow(tnew2);
	if(ret)
		ERRRET("Supernormalization failed.")

	//Step 2: Log likelihood ratios from nonpermuted data
	LOG(9,"Calculating real log likelihood ratios...")
	pij_rank_llr(tnew,tnew2,p);
	if(nodiag)
	{
		vv=MATRIXFF(diagonal)(p);
		VECTORFF(set_zero)(&vv.vector);
	}
	//Step 3: Convert log likelihood ratios to probabilities
	if((ret=pij(p,ns,nodiag,0)))
		LOG(1,"Failed to convert log likelihood ratios to probabilities.")

	//Cleanup
	CLEANUP
	return ret;
#undef	CLEANUP		
}

int pij_rank_a(const MATRIXF* t,const MATRIXF* t2,MATRIXF* p,char nodiag,size_t memlimit)
{
	return pij_rank_any(t,t2,p,nodiag,pij_rank_llrtopij_a,memlimit);
}

int pij_rank(const MATRIXF* t,const MATRIXF* t2,MATRIXF* p,char nodiag,size_t memlimit)
{
	return pij_rank_a(t,t2,p,nodiag,memlimit);
}












