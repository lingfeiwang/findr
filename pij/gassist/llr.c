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
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include "../../base/gsl/blas.h"
#include "../../base/random.h"
#include "../../base/const.h"
#include "../../base/logger.h"
#include "../../base/macros.h"
#include "../../base/data_process.h"
#include "../../base/threading.h"
#include "llr.h"


//Struct for parameters of pij_gassist_llr_block_buffed (below)
struct pij_gassist_llr_block_buffed_params{
	//Number of (E,A) pairs
	size_t			ng;
	//Number of possible genotypes per SNP
	size_t			nv;
	//(nv,ng) Genotype ratio matrix
	const MATRIXF*	mratio;
	//(nv,ng) Mean per genotype matrix
	const MATRIXF*	mmean1;
	//[nv](ng,nt) Mean per genotype matrix
	const MATRIXF**	mmean2;
	//(ng) Output vector for log likelihood ratio 1
	VECTORF* 		llr1;		
	//(ng,nt) Output matrix for log likelihood ratio 2
	MATRIXF*		llr2;
	//(ng,nt) Output matrix for log likelihood ratio 3
	MATRIXF*		llr3;
	//(ng,nt) Output matrix for log likelihood ratio 4
	MATRIXF*		llr4;
	//(ng,nt) Output matrix for log likelihood ratio 5,
	//Also used as correlation matrix for input
	MATRIXF*		llr5;
	//[nv](ng,nt) Buffer matrix
	MATRIXF**		mb1;
};

void pij_gassist_llr_ratioandmean_1v1(const MATRIXG* g,const MATRIXF* t1,const MATRIXF* t2,MATRIXF* mratio,MATRIXF* mmean1,MATRIXF* mmean2,size_t nv)
{
	size_t	ng=g->size1;
	size_t	ns=g->size2;
	size_t	i,j;
	size_t	val;
	FTYPE	f1;
	
	MATRIXFF(set_zero)(mratio);
	MATRIXFF(set_zero)(mmean1);
	MATRIXFF(set_zero)(mmean2);
	
	//Calculating sums
	for(i=0;i<ng;i++)
		for(j=0;j<ns;j++)
		{
			val=(size_t)MATRIXGF(get)(g,i,j);
			MATRIXFF(ptr)(mratio,i,val)[0]++;
			*MATRIXFF(ptr)(mmean1,i,val)+=MATRIXFF(get)(t1,i,j);
			*MATRIXFF(ptr)(mmean2,i,val)+=MATRIXFF(get)(t2,i,j);
		}
	//Convert to means
	for(i=0;i<ng;i++)
		for(j=0;j<nv;j++)
		{
			f1=MATRIXFF(get)(mratio,i,j);
			*MATRIXFF(ptr)(mmean1,i,j)/=f1+FTYPE_MIN;
			*MATRIXFF(ptr)(mmean2,i,j)/=f1+FTYPE_MIN;
		}
	MATRIXFF(scale)(mratio,1/(FTYPE)ns);
}

/* Calculates the ratio and mean of all transcripts (t) for genes (g) with existing buffers.
 * g:		MATRIXG (ng,ns) genotype data, for multiple SNP and samples. Each element takes the value 0 to nv-1
 * t:		MATRIXF (ng,ns) of transcript data for A
 * t2:		MATRIXF (nt,ns) of transcript data for B
 * mratio:	MATRIXF (nv,ng) of the ratio of samples for each SNP type. For return purpose.
 * mmean1:	MATRIXF (nv,ng) of the means of each transcript among the samples of a specific SNP type. For return purpose.
 * mmean2:	MATRIXF[nv] (ng,nt) of the means of each transcript among the samples of a specific SNP type. For return purpose.
 * nv:		number of possible values of g.
 * mb1:		MATRIXF[nv] (ng,ns) Buffer matrix
 * vb:		buffer. const VECTORF (ns). Must be set to 1 for all elements.
 */
static void pij_gassist_llr_ratioandmean_buffed(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,MATRIXF* mratio,MATRIXF* mmean1,MATRIXF** mmean2,size_t nv,MATRIXF** mb1,const VECTORF* vb)
{
	size_t	ng=g->size1;
	size_t	ns=t->size2;
	size_t	i,j;
	VECTORFF(view) vv;
	FTYPE	t1;
		
	//Initialize matrices for ratio and mean.
	for(i=0;i<nv;i++)
		MATRIXFF(set_zero)(mb1[i]);
	for(i=0;i<ng;i++)
		for(j=0;j<ns;j++)
			MATRIXFF(set)(mb1[MATRIXGF(get)(g,i,j)],i,j,1);
	
	for(i=0;i<nv;i++)
	{
		//Count into (mratio)
		vv=MATRIXFF(row)(mratio,i);
		BLASF(gemv)(CblasNoTrans,1,mb1[i],vb,0,&(vv.vector));
		//Sum for t2
		BLASF(gemm)(CblasNoTrans,CblasTrans,1,mb1[i],t2,0,mmean2[i]);
		//Sum for t1
		vv=MATRIXFF(row)(mmean1,i);
		MATRIXFF(mul_elements)(mb1[i],t);
		BLASF(gemv)(CblasNoTrans,1,mb1[i],vb,0,&vv.vector);
		//Mean=Sum/Count
		for(j=0;j<ng;j++)
		{
			t1=MATRIXFF(get)(mratio,i,j);
			MATRIXFF(ptr)(mmean1,i,j)[0]/=t1+FTYPE_MIN;
			vv=MATRIXFF(row)(mmean2[i],j);
			VECTORFF(scale)(&(vv.vector),1/(t1+FTYPE_MIN));
		}
	}

	//Ratio=count/total count
	MATRIXFF(scale)(mratio,((FTYPE)1)/(FTYPE)ns);
}


/* Calculates the ratio and mean of all transcripts (t) for genes (g) with existing buffers.
 * g:		MATRIXG (ng,ns) genotype data, for multiple SNP and samples. Each element takes the value 0 to nv-1
 * t:		MATRIXF (nt,ns) of transcript data.
 * t2:		MATRIXF (nt,ns) of transcript data for B
 * mratio:	MATRIXF (nv,ng) of the ratio of samples for each SNP type. For return purpose.
 * mmean1:	MATRIXF (nv,ng) of the means of each transcript among the samples of a specific SNP type. For return purpose.
 * mmean2:	MATRIXF[nv] (ng,nt) of the means of each transcript among the samples of a specific SNP type. For return purpose.
 * nv:		number of possible values of g.
 * vb:		buffer. const VECTORF (ns). Must be set to 1 for all elements.
 * Return:	0 on success
 */
static int pij_gassist_llr_ratioandmean(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,MATRIXF* mratio,MATRIXF* mmean1,MATRIXF** mmean2,size_t nv,const VECTORF* vb)
{
#define	CLEANUP			CLEANAMMATF(mb1,nv)
	size_t	i;
	int		ret;
	
	//Memory allocation
	AUTOCALLOC(MATRIXF*,mb1,nv,200)
	if(!mb1)
		ERRRET("Not enough memory.")
	ret=1;
	for(i=0;i<nv;i++)
		ret=ret&&(mb1[i]=MATRIXFF(alloc)(g->size1,g->size2));
	if(!ret)
		ERRRET("Not enough memory.")
		
	pij_gassist_llr_ratioandmean_buffed(g,t,t2,mratio,mmean1,mmean2,nv,mb1,vb);
	CLEANUP
	return 0;
#undef	CLEANUP
}


/* Calculates the 5 log likelihood ratios of Trigger with nonpermuted data in block form:
 * 1. E->A v.s. E no relation with A
 * 2. A<-E->B with A--B v.s. E->A<-B
 * 3. E->A->B v.s. A<-E->B with A--B
 * 4. A<-E->B with A->B v.s. E->A
 * 5. A<-E->B with A->B v.s. A<-E->B
 * Uses GSL BLAS.
 * Note: for each row, g must be the best eQTL of t of the same row.
 */
static void pij_gassist_llr_block_buffed(const struct pij_gassist_llr_block_buffed_params* p)
{
	size_t	i,j;
	FTYPE	f1;
	VECTORFF(view)	vv;
	size_t	ng=p->ng;
#ifndef NDEBUG
	size_t	nt=p->llr4->size2;
#endif
	assert(p->nv&&(p->mratio->size1==p->nv)&&(p->mmean1->size1==p->nv));	
	assert(ng&&(p->mratio->size2==ng)&&(p->mmean1->size2==ng)
		&&(p->mmean2[0]->size1==ng)&&(p->llr5->size1==ng)&&(p->llr1->size==ng)
		&&(p->llr4->size1==ng)&&(p->llr2->size1==ng)&&(p->llr3->size1==ng)
		&&(p->mb1[0]->size1==ng));
	assert(nt&&(p->mmean2[0]->size2==nt)&&(p->llr5->size2==nt)
		&&(p->llr2->size2==nt)&&(p->llr3->size2==nt)&&(p->mb1[0]->size2==nt));
	
	for(i=0;i<p->nv;i++)
	{
		MATRIXFF(memcpy)(p->mb1[i],p->mmean2[i]);
		for(j=0;j<ng;j++)
		{
			f1=MATRIXFF(get)(p->mratio,i,j);
			vv=MATRIXFF(row)(p->mb1[i],j);
			VECTORFF(scale)(&vv.vector,f1);
		}
	}
	
	//llr4=llr2=1-sum_alpha f_{alpha i}mu_{alpha ij}^2, using llr3 as buff
	MATRIXFF(set_zero)(p->llr4);
	for(i=0;i<p->nv;i++)
	{
		MATRIXFF(memcpy)(p->llr3,p->mb1[i]);
		MATRIXFF(mul_elements)(p->llr3,p->mmean2[i]);
		MATRIXFF(add)(p->llr4,p->llr3);
	}
	MATRIXFF(scale)(p->llr4,-1);
	MATRIXFF(add_constant)(p->llr4,1);
	MATRIXFF(memcpy)(p->llr2,p->llr4);
	//llr1=1-sum_alpha f_{alpha i}mu_{alpha ii}^2, using llr3[:,0] as buff
	VECTORFF(set_zero)(p->llr1);
	vv=MATRIXFF(column)(p->llr3,0);
	for(i=0;i<p->nv;i++)
	{
		VECTORFF(const_view) vv2=MATRIXFF(const_row)(p->mratio,i);
		VECTORFF(const_view) vv3=MATRIXFF(const_row)(p->mmean1,i);
		VECTORFF(memcpy)(&vv.vector,&vv2.vector);
		VECTORFF(mul)(&vv.vector,&vv3.vector);
		VECTORFF(mul)(&vv.vector,&vv3.vector);
		VECTORFF(add)(p->llr1,&vv.vector);
	}
	VECTORFF(scale)(p->llr1,-1);
	VECTORFF(add_constant)(p->llr1,1);
		
	//mb1[i][j,k]=mb1[i][j,k]*mmean1[i][j]=f_{ij}mu_{ijj}mu_{ijk}
	for(i=0;i<p->nv;i++)
		for(j=0;j<ng;j++)
		{
			vv=MATRIXFF(row)(p->mb1[i],j);
			VECTORFF(scale)(&vv.vector,MATRIXFF(get)(p->mmean1,i,j));
		}
	//mb1[0]=sum_alpha f_{alpha i}mu_{alpha ii}mu_{alpha ij}
	for(i=1;i<p->nv;i++)
		MATRIXFF(add)(p->mb1[0],p->mb1[i]);
	//mb1[0]=(rho_{ij}-sum_alpha f_{alpha i}mu_{alpha ii}mu_{alpha ij})^2
	MATRIXFF(sub)(p->mb1[0],p->llr5);
	MATRIXFF(mul_elements)(p->mb1[0],p->mb1[0]);
	
	//llr4=(1-sum_alpha f_{alpha i}mu_{alpha ii}^2)(1-sum_alpha f_{alpha i}mu_{alpha ij}^2)
	for(i=0;i<ng;i++)
	{
		vv=MATRIXFF(row)(p->llr4,i);
		VECTORFF(scale)(&vv.vector,VECTORFF(get)(p->llr1,i));
	}
	//llr4=ML3=(1-sum_alpha f_{alpha i}mu_{alpha ii}^2)(1-sum_alpha f_{alpha i}mu_{alpha ij}^2)-(rho_{ij}-sum_alpha f_{alpha i}mu_{alpha ii}mu_{alpha ij})^2
	MATRIXFF(sub)(p->llr4,p->mb1[0]);
	
	//mb1[1]=(1-rho_{ij}^2)
	MATRIXFF(memcpy)(p->mb1[1],p->llr5);
	MATRIXFF(mul_elements)(p->mb1[1],p->llr5);
	MATRIXFF(scale)(p->mb1[1],-1);
	MATRIXFF(add_constant)(p->mb1[1],1);
	
	//ALL log
	for(i=0;i<ng;i++)
	{
		VECTORFF(set)(p->llr1,i,(FTYPE)log(VECTORFF(get)(p->llr1,i)));
		for(j=0;j<p->llr4->size2;j++)
		{
			MATRIXFF(set)(p->llr4,i,j,(FTYPE)log(MATRIXFF(get)(p->llr4,i,j)));
			MATRIXFF(set)(p->llr2,i,j,(FTYPE)log(MATRIXFF(get)(p->llr2,i,j)));
			MATRIXFF(set)(p->mb1[0],i,j,(FTYPE)log(MATRIXFF(get)(p->mb1[0],i,j)));
			MATRIXFF(set)(p->mb1[1],i,j,(FTYPE)log(MATRIXFF(get)(p->mb1[1],i,j)));
		}
	}
	
	
	//llr4=log((1-sum_alpha f_{alpha i}mu_{alpha ii}^2)(1-sum_alpha f_{alpha i}mu_{alpha ij}^2)-(rho_{ij}-sum_alpha f_{alpha i}mu_{alpha ii}mu_{alpha ij})^2)-llr1
	for(i=0;i<ng;i++)
	{
		vv=MATRIXFF(row)(p->llr4,i);
		VECTORFF(add_constant)(&vv.vector,-VECTORFF(get)(p->llr1,i));
	}
	
	//llr3 final
	MATRIXFF(memcpy)(p->llr3,p->llr4);
	MATRIXFF(sub)(p->llr3,p->mb1[1]);
	MATRIXFF(scale)(p->llr3,-0.5);
	//llr5 final
	MATRIXFF(memcpy)(p->llr5,p->llr4);
	MATRIXFF(sub)(p->llr5,p->llr2);
	MATRIXFF(scale)(p->llr5,-0.5);

	//llr4 final
	MATRIXFF(scale)(p->llr4,-0.5);
	
	//llr2 final
	MATRIXFF(scale)(p->llr2,-0.5);
	
	//llr1 final
	VECTORFF(scale)(p->llr1,-0.5);
	
	//Bounding from 0
	VECTORFF(bound_below)(p->llr1,0);
	MATRIXFF(bound_below)(p->llr2,0);
	MATRIXFF(bound_below)(p->llr3,0);
	MATRIXFF(bound_below)(p->llr4,0);
	MATRIXFF(bound_below)(p->llr5,0);
}

/* Wrapper of pij_gassist_llr_block_buffed. Performs memory allocation and pre-calculations of
 * categorical mean and ratio, and the covariance matrix before invoking pij_gassist_llr_block_buffed
 * g:		MATRIXF (ng,ns) Full genotype data matrix
 * t:		MATRIXF (ng,ns) Supernormalized transcript data matrix for A
 * t2:		MATRIXF (nt,ns) Supernormalized transcript data matrix for B
 * nv:		Number of possible values for each genotype
 * llr1:	VECTORF (end-start). Log likelihood ratios for test 1.
 * llr2:	MATRIXF (end-start,nt). Log likelihood ratios for test 2.
 * llr3:	MATRIXF (end-start,nt). Log likelihood ratios for test 3.
 * llr4:	MATRIXF (end-start,nt). Log likelihood ratios for test 4.
 * llr5:	MATRIXF (end-start,nt). Log likelihood ratios for test 5.
 * vb1:		VECTORF (ns). Constant buffer vector, must be set to 1 initially.
 * Return:	0 on success.
 * Notes:	1.	block range only applicable to A, all other transcripts as B are always considered.
 * 			2.	for each row, g must be the best eQTL of t of the same row.
 * 			3.	Definitions:	ng: number of SNPs=number of transcripts for A
 * 								nt: number of transcripts for B
 * 								ns: number of samples
 */
static int pij_gassist_llr_block(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,size_t nv,VECTORF* llr1,MATRIXF* llr2,MATRIXF* llr3,MATRIXF* llr4,MATRIXF* llr5,const VECTORF* vb1)
{
#define	CLEANUP	CLEANMATF(mratio)CLEANMATF(mmean1)CLEANAMMATF(mmean2,nv)CLEANAMMATF(mmb1,nv)
	size_t	i,j;
	int		ret;
	size_t	ng=g->size1;
	size_t	nt=t2->size1;
	MATRIXF	*mratio,*mmean1;	//Buffer matrix (nv,ng)
	struct pij_gassist_llr_block_buffed_params llp;
	
	//Memory allocation
	mmean1=mratio=0;
	AUTOCALLOC(MATRIXF*,mmean2,nv,200)
	AUTOCALLOC(MATRIXF*,mmb1,nv,200)
	if(!(mmean2&&mmb1))
		ERRRET("Not enough memory.")
	for(i=0,j=1;i<nv;i++)
		j=j&&(mmean2[i]=MATRIXFF(alloc)(ng,nt))&&(mmb1[i]=MATRIXFF(alloc)(ng,nt));
	mratio=MATRIXFF(alloc)(nv,ng);
	mmean1=MATRIXFF(alloc)(nv,ng);
	if(!(mratio&&mmean1&&j))
		ERRRET("Not enough memory.")
	
	//Calculate ratio and mean
	ret=pij_gassist_llr_ratioandmean(g,t,t2,mratio,mmean1,mmean2,nv,vb1);
	if(ret)
		ERRRET("Not enough memory.")
	//Calculate covariance
	MATRIXFF(cov2_bounded)(t,t2,llr5);
	//Initialize parameter pack
	llp.ng=ng;
	llp.nv=nv;
	llp.mratio=(const MATRIXF*)mratio;
	llp.mmean1=(const MATRIXF*)mmean1;
	llp.mmean2=(const MATRIXF**)mmean2;
	llp.llr1=llr1;
	llp.llr2=llr2;
	llp.llr3=llr3;
	llp.llr4=llr4;
	llp.llr5=llr5;
	llp.mb1=mmb1;
	
	pij_gassist_llr_block_buffed(&llp);
	
	CLEANUP
	return 0;
#undef CLEANUP
}

int pij_gassist_llr(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,VECTORF* llr1,MATRIXF* llr2,MATRIXF* llr3,MATRIXF* llr4,MATRIXF* llr5,size_t nv)
{
#define	CLEANUP			CLEANVECF(vbuff1)
	VECTORF	*vbuff1=0;		//(ns) Const buffer, set to all 1.
	int		ret;
#ifndef NDEBUG
	size_t	ng,nt,ns;

	ng=g->size1;
	nt=t2->size1;
	ns=t->size2;
#endif

	//Validation
	assert(!((g->size2!=ns)||(t2->size2!=ns)||(t->size1!=ng)||(llr2->size1!=ng)||(llr2->size2!=nt)||(llr3->size1!=ng)||(llr3->size2!=nt)||(llr4->size1!=ng)||(llr4->size2!=nt)||(llr5->size1!=ng)||(llr5->size2!=nt)));
	assert(!(llr1->size!=ng));
	assert(!(nv>CONST_NV_MAX));
	assert(nv>MATRIXGF(max)(g));
	
	{
		GTYPE	tg=MATRIXGF(max)(g);
		if(tg>=nv)
			ERRRET("Maximum genotype value "PRINTFSIZET" exceeds the stated maximum possible value "PRINTFSIZET". Please check your input genotype matrix and allele count.",tg,nv-1)
	}
	//Buff unit vector
	vbuff1=VECTORFF(alloc)(g->size2);
	if(!(vbuff1))
		ERRRET("Not enough memory.")
	VECTORFF(set_all)(vbuff1,1);
	
	ret=0;
	#pragma omp parallel
	{
		size_t	n1,n2;
		int		retth;
		
		threading_get_startend(t->size1,&n1,&n2);
		if(n2>n1)
		{
			MATRIXGF(const_view) mvg=MATRIXGF(const_submatrix)(g,n1,0,n2-n1,g->size2);
			MATRIXFF(const_view) mvt=MATRIXFF(const_submatrix)(t,n1,0,n2-n1,t->size2);
			VECTORFF(view)	vvllr1;
			MATRIXFF(view)	mvllr2,mvllr3,mvllr4,mvllr5;
			vvllr1=VECTORFF(subvector)(llr1,n1,n2-n1);
			mvllr2=MATRIXFF(submatrix)(llr2,n1,0,n2-n1,llr2->size2);
			mvllr3=MATRIXFF(submatrix)(llr3,n1,0,n2-n1,llr3->size2);
			mvllr4=MATRIXFF(submatrix)(llr4,n1,0,n2-n1,llr4->size2);
			mvllr5=MATRIXFF(submatrix)(llr5,n1,0,n2-n1,llr5->size2);
			retth=pij_gassist_llr_block(&mvg.matrix,&mvt.matrix,t2,nv,&vvllr1.vector,&mvllr2.matrix,&mvllr3.matrix,&mvllr4.matrix,&mvllr5.matrix,vbuff1);
			#pragma omp atomic
			ret+=retth;
		}
	}

	if(ret)
		ERRRET("Failed to calculate nonpermuted log likelihood ratios.")

	//Cleanup
	CLEANUP
	return 0;
#undef	CLEANUP		
}

