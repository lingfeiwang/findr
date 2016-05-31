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

void pij_gassist_llr_ratioandmean_buffed(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,MATRIXF* mratio,MATRIXF* mmean1,MATRIXF** mmean2,size_t nv,MATRIXF** mb1,const VECTORF* vb)
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

int pij_gassist_llr_ratioandmean(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,MATRIXF* mratio,MATRIXF* mmean1,MATRIXF** mmean2,size_t nv,const VECTORF* vb)
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

void pij_gassist_llr_block_buffed(const struct pij_gassist_llr_block_buffed_params* p)
{
	size_t	i,j;
	FTYPE	f1;
	VECTORFF(view)	vv;
	size_t	ng=p->ng;
#ifndef NDEBUG
	size_t	nt=p->llr2b->size2;
#endif
	assert(p->nv&&(p->mratio->size1==p->nv)&&(p->mmean1->size1==p->nv));	
	assert(ng&&(p->mratio->size2==ng)&&(p->mmean1->size2==ng)
		&&(p->mmean2[0]->size1==ng)&&(p->mcov->size1==ng)&&(p->llr1->size==ng)
		&&(p->llr2b->size1==ng)&&(p->llr2c->size1==ng)&&(p->llr3->size1==ng)
		&&(p->mb1[0]->size1==ng));
	assert(nt&&(p->mmean2[0]->size2==nt)&&(p->mcov->size2==nt)
		&&(p->llr2c->size2==nt)&&(p->llr3->size2==nt)&&(p->mb1[0]->size2==nt));
	//mb1=f_{alpha i}*mu_{alpha ij}
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
	
	//llr2b=llr2c=1-sum_alpha f_{alpha i}mu_{alpha ij}^2, using llr3 as buff
	MATRIXFF(set_zero)(p->llr2b);
	for(i=0;i<p->nv;i++)
	{
		MATRIXFF(memcpy)(p->llr3,p->mb1[i]);
		MATRIXFF(mul_elements)(p->llr3,p->mmean2[i]);
		MATRIXFF(add)(p->llr2b,p->llr3);
	}
	MATRIXFF(scale)(p->llr2b,-1);
	MATRIXFF(add_constant)(p->llr2b,1);
	MATRIXFF(memcpy)(p->llr2c,p->llr2b);
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
	MATRIXFF(sub)(p->mb1[0],p->mcov);
	MATRIXFF(mul_elements)(p->mb1[0],p->mb1[0]);
	
	//llr2b=(1-sum_alpha f_{alpha i}mu_{alpha ii}^2)(1-sum_alpha f_{alpha i}mu_{alpha ij}^2)
	for(i=0;i<ng;i++)
	{
		vv=MATRIXFF(row)(p->llr2b,i);
		VECTORFF(scale)(&vv.vector,VECTORFF(get)(p->llr1,i));
	}
	//llr2b=ML3=(1-sum_alpha f_{alpha i}mu_{alpha ii}^2)(1-sum_alpha f_{alpha i}mu_{alpha ij}^2)-(rho_{ij}-sum_alpha f_{alpha i}mu_{alpha ii}mu_{alpha ij})^2
	MATRIXFF(sub)(p->llr2b,p->mb1[0]);
	
	//mb1[1]=(1-rho_{ij}^2)
	MATRIXFF(memcpy)(p->mb1[1],p->mcov);
	MATRIXFF(mul_elements)(p->mb1[1],p->mcov);
	MATRIXFF(scale)(p->mb1[1],-1);
	MATRIXFF(add_constant)(p->mb1[1],1);
	
	//ALL log
	for(i=0;i<ng;i++)
	{
		VECTORFF(set)(p->llr1,i,(FTYPE)log(VECTORFF(get)(p->llr1,i)));
		for(j=0;j<p->llr2b->size2;j++)
		{
			MATRIXFF(set)(p->llr2b,i,j,(FTYPE)log(MATRIXFF(get)(p->llr2b,i,j)));
			MATRIXFF(set)(p->llr2c,i,j,(FTYPE)log(MATRIXFF(get)(p->llr2c,i,j)));
			MATRIXFF(set)(p->mb1[0],i,j,(FTYPE)log(MATRIXFF(get)(p->mb1[0],i,j)));
			MATRIXFF(set)(p->mb1[1],i,j,(FTYPE)log(MATRIXFF(get)(p->mb1[1],i,j)));
		}
	}
	
	//llr2b=log((1-sum_alpha f_{alpha i}mu_{alpha ii}^2)(1-sum_alpha f_{alpha i}mu_{alpha ij}^2)-(rho_{ij}-sum_alpha f_{alpha i}mu_{alpha ii}mu_{alpha ij})^2)-llr1
	for(i=0;i<ng;i++)
	{
		vv=MATRIXFF(row)(p->llr2b,i);
		VECTORFF(add_constant)(&vv.vector,-VECTORFF(get)(p->llr1,i));
	}
	
	//llr3 final (Michoel's version)
	MATRIXFF(memcpy)(p->llr3,p->llr2b);
	MATRIXFF(sub)(p->llr3,p->mb1[1]);
	MATRIXFF(scale)(p->llr3,-0.5);

	//llr2b final
	MATRIXFF(scale)(p->llr2b,-0.5);
	
	//llr2c final
	MATRIXFF(scale)(p->llr2c,-0.5);
	
	//llr1 final
	VECTORFF(scale)(p->llr1,-0.5);
	
	//Bounding from 0
	VECTORFF(bound_below)(p->llr1,0);
	MATRIXFF(bound_below)(p->llr2b,0);
	MATRIXFF(bound_below)(p->llr2c,0);
	MATRIXFF(bound_below)(p->llr3,0);
}

int pij_gassist_llr_block(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,size_t nv,VECTORF* llr1,MATRIXF* llr2b,MATRIXF* llr2c,MATRIXF* llr3,const VECTORF* vb1)
{
#define	CLEANUP	CLEANMATF(mratio)CLEANMATF(mcov)CLEANMATF(mmean1)CLEANAMMATF(mmean2,nv)CLEANAMMATF(mmb1,nv)
	size_t	i,j;
	int		ret;
	size_t	ng=g->size1;
	size_t	nt=t2->size1;
	MATRIXF	*mratio,*mcov,*mmean1;	//Buffer matrix (nv,ng),(ng,nt),(nv,ng)
	struct pij_gassist_llr_block_buffed_params llp;
	
	//Memory allocation
	mmean1=mratio=mcov=0;
	AUTOCALLOC(MATRIXF*,mmean2,nv,200)
	AUTOCALLOC(MATRIXF*,mmb1,nv,200)
	if(!(mmean2&&mmb1))
		ERRRET("Not enough memory.")
	for(i=0,j=1;i<nv;i++)
		j=j&&(mmean2[i]=MATRIXFF(alloc)(ng,nt))&&(mmb1[i]=MATRIXFF(alloc)(ng,nt));
	mmean1=MATRIXFF(alloc)(nv,ng);
	mratio=MATRIXFF(alloc)(nv,ng);
	mcov=MATRIXFF(alloc)(ng,nt);
	if(!(mratio&&mcov&&j))
		ERRRET("Not enough memory.")
	
	//Calculate ratio and mean
	ret=pij_gassist_llr_ratioandmean(g,t,t2,mratio,mmean1,mmean2,nv,vb1);
	if(ret)
		ERRRET("Not enough memory.")
	//Calculate covariance
	MATRIXFF(cov2_bounded)(t,t2,mcov);
	//Initialize parameter pack
	llp.ng=ng;
	llp.nv=nv;
	llp.mratio=(const MATRIXF*)mratio;
	llp.mmean1=(const MATRIXF*)mmean1;
	llp.mmean2=(const MATRIXF**)mmean2;
	llp.mcov=mcov;
	llp.llr1=llr1;
	llp.llr2b=llr2b;
	llp.llr2c=llr2c;
	llp.llr3=llr3;
	llp.mb1=mmb1;
	
	pij_gassist_llr_block_buffed(&llp);
	
	CLEANUP
	return 0;
#undef CLEANUP
}


void pij_gassist_llr2c_buffed(const MATRIXG* g,const MATRIXF* t,const MATRIXF* tp,size_t nv,MATRIXF* llr2,MATRIXF* mratio,MATRIXF* mmean1,MATRIXF** mmean2,MATRIXF** mmb1,const VECTORF* vb1)
{
	size_t	i,j;
	FTYPE	f1;
	VECTORFF(view)	vv;
	size_t	ng=g->size1;

	//ratio and mean
	pij_gassist_llr_ratioandmean_buffed(g,t,tp,mratio,mmean1,mmean2,nv,mmb1,vb1);
	
	//Calculation
	MATRIXFF(set_zero)(llr2);
	for(i=0;i<nv;i++)
	{
		MATRIXFF(mul_elements)(mmean2[i],mmean2[i]);
		for(j=0;j<ng;j++)
		{
			f1=MATRIXFF(get)(mratio,i,j);
			vv=MATRIXFF(row)(mmean2[i],j);
			VECTORFF(scale)(&vv.vector,f1);
		}
		MATRIXFF(sub)(llr2,mmean2[i]);
	}
	MATRIXFF(add_constant)(llr2,1);

	//log
	for(i=0;i<ng;i++)
		for(j=0;j<llr2->size2;j++)
			MATRIXFF(set)(llr2,i,j,(FTYPE)log(MATRIXFF(get)(llr2,i,j)));
	
	//llr2 final
	MATRIXFF(scale)(llr2,-0.5);
	MATRIXFF(bound_below)(llr2,0);
}

int pij_gassist_llr2c(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,size_t nv,MATRIXF* llr2)
{
#define	CLEANUP CLEANVECF(vb1)\
				CLEANMATF(mratio)CLEANMATF(mmean1)\
				CLEANAMMATF(mmean2,nv)CLEANAMMATF(mmb1,nv)
	size_t	i,j;
	size_t	ng=g->size1;
	size_t	ns=g->size2;	//Number of samples
	size_t	nt=t2->size1;	//Number of B's before permutations
	MATRIXF	*mratio,*mmean1;
	VECTORF *vb1;
	
	vb1=0;
	mratio=mmean1=0;
	
	AUTOCALLOC(MATRIXF*,mmean2,nv,200)
	AUTOCALLOC(MATRIXF*,mmb1,nv,200)
	if(!(mmean2&&mmb1))
		ERRRET("Not enough memory.")
	
	for(i=0,j=1;i<nv;i++)
	{
		mmean2[i]=MATRIXFF(alloc)(ng,nt);
		mmb1[i]=MATRIXFF(alloc)(ng,ns);
		j=j&&mmean2[i]&&mmb1[i];
	}
	vb1=VECTORFF(alloc)(ns);
	mratio=MATRIXFF(alloc)(nv,ng);
	mmean1=MATRIXFF(alloc)(nv,ng);
	if(!(j&&vb1&&mratio))
		ERRRET("Not enough memory.")
	VECTORFF(set_all)(vb1,1);
	
	pij_gassist_llr2c_buffed(g,t,t2,nv,llr2,mratio,mmean1,mmean2,mmb1,vb1);
	
	CLEANUP
	return 0;
#undef CLEANUP
}

void pij_gassist_llr2b_buffed(const MATRIXG* g,const MATRIXF* t,const MATRIXF* tp,size_t nv,MATRIXF* llr2,MATRIXF* mratio,MATRIXF* mmean1,MATRIXF** mmean2,MATRIXF** mmb1,MATRIXF** mmb2,MATRIXF* mb1,MATRIXF* mb2,MATRIXF* mb3,const VECTORF* vb1,VECTORF* vb2)
{
	const struct pij_gassist_llr_block_buffed_params p={g->size1,nv,mratio,mmean1,(const MATRIXF**)mmean2,mb1,vb2,llr2,mb2,mb3,mmb2};
	pij_gassist_llr_ratioandmean_buffed(g,t,tp,mratio,mmean1,mmean2,nv,mmb1,vb1);
	MATRIXFF(cov2_bounded)(t,tp,mb1);
	pij_gassist_llr_block_buffed(&p);
}

int pij_gassist_llr2b(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,size_t nv,MATRIXF* llr2)
{
#define	CLEANUP CLEANVECF(vb1)CLEANVECF(vb2)\
				CLEANMATF(mratio)CLEANMATF(mmean1)CLEANMATF(mb1)CLEANMATF(mb2)\
				CLEANMATF(mb3)\
				CLEANMMATF(mmean2,nv)CLEANMMATF(mmb1,nv)CLEANMMATF(mmb2,nv)
	size_t	i,j;
	size_t	ng=g->size1;
	size_t	ns=g->size2;	//Number of samples
	size_t	nt=t2->size1;	//Number of B's before permutations
	MATRIXF	*mratio,*mmean1,**mmean2,**mmb1,**mmb2,*mb1,*mb2,*mb3;
	VECTORF *vb1,*vb2;
	
	vb1=vb2=0;
	mratio=mmean1=mb1=mb2=mb3=0;
	mmean2=mmb1=mmb2=0;

	CALLOCSIZE(mmean2,nv);
	CALLOCSIZE(mmb1,nv);
	CALLOCSIZE(mmb2,nv);
	if(!(mmean2&&mmb1&&mmb2))
		ERRRET("Not enough memory.")
	
	for(i=0,j=1;i<nv;i++)
	{
		mmean2[i]=MATRIXFF(alloc)(ng,nt);
		mmb1[i]=MATRIXFF(alloc)(ng,ns);
		mmb2[i]=MATRIXFF(alloc)(ng,nt);
		j=j&&mmean2[i]&&mmb1[i]&&mmb2[i];
	}
	vb1=VECTORFF(alloc)(ns);
	vb2=VECTORFF(alloc)(ng);
	mratio=MATRIXFF(alloc)(nv,ng);
	mmean1=MATRIXFF(alloc)(nv,ng);
	mb1=MATRIXFF(alloc)(ng,nt);
	mb2=MATRIXFF(alloc)(ng,nt);
	mb3=MATRIXFF(alloc)(ng,nt);
	if(!(j&&vb1&&vb2&&mratio&&mmean1&&mb1&&mb2&&mb3))
		ERRRET("Not enough memory.")
	VECTORFF(set_all)(vb1,1);
	
	pij_gassist_llr2b_buffed(g,t,t2,nv,llr2,mratio,mmean1,mmean2,mmb1,mmb2,mb1,mb2,mb3,vb1,vb2);
	
	CLEANUP
	return 0;
#undef CLEANUP
}



FTYPE pij_gassist_llr3_one_from_stats_buffed(FTYPE fcov,const VECTORF* vratio,const VECTORF* vmean1,const VECTORF* vmean2,VECTORF* vb1)
{
	FTYPE f1,f2,f3;
	//vb1=f_alpha mu_{alpha 1}
	VECTORFF(memcpy)(vb1,vratio);
	VECTORFF(mul)(vb1,vmean1);
	//f1=1-sum_alpha f_alpha mu_{alpha 1}^2
	BLASF(dot)(vb1,vmean1,&f1);
	f1=1-f1;
	//f2=(rho-sum_alpha f_alpha mu_{alpha 1}mu_{alpha 2})^2
	BLASF(dot)(vb1,vmean2,&f2);
	f2-=fcov;
	f2*=f2;
	//vb1=f_alpha mu_{alpha 2}
	VECTORFF(memcpy)(vb1,vratio);
	VECTORFF(mul)(vb1,vmean2);
	//f3=1-sum_alpha f_alpha mu_{alpha 2}^2
	BLASF(dot)(vb1,vmean2,&f3);
	f3=1-f3;
	//f3=log((1-sum_alpha f_alpha mu_{alpha 1}^2)(1-sum_alpha f_alpha mu_{alpha 2}^2)-(rho-sum_alpha f_alpha mu_{alpha 1}mu_{alpha 2})^2)
	f3=(FTYPE)log(f1*f3-f2);
	//f1=log(1-sum_alpha f_alpha mu_{alpha 1}^2)+log(1-rho^2)
	f1=(FTYPE)(log(f1)+log(1-fcov*fcov));
	//f3=final
	f3=(f3-f1)/2;
	if(f3>0)
		f3=0;
	return f3;
}

void pij_gassist_llr3_E1_buffed(const VECTORG* g,const VECTORF* t1,const MATRIXF* t2,size_t nv,VECTORF* llr3,const VECTORF* vb1,VECTORF *vratio,VECTORF* vmean1,MATRIXF *mmean,VECTORF *vcov,VECTORF *vb2,VECTORF *vb3,MATRIXF *mb1)
{
	size_t	i;
	size_t	ns=g->size;
	size_t	nt=t2->size1;
	FTYPE	fb1;
	
	//Categorical ratio and mean
	MATRIXFF(set_zero)(mb1);
	for(i=0;i<ns;i++)
		MATRIXFF(set)(mb1,(size_t)VECTORGF(get)(g,i),i,1);
	VECTORFF(set_all)(vratio,FTYPE_MIN);
	BLASF(gemv)(CblasNoTrans,1,mb1,vb1,1,vratio);
	BLASF(gemv)(CblasNoTrans,1,mb1,t1,0,vmean1);
	VECTORFF(div)(vmean1,vratio);
	BLASF(gemm)(CblasNoTrans,CblasTrans,1,t2,mb1,0,mmean);
	for(i=0;i<nv;i++)
	{
		VECTORFF(view)	vv=MATRIXFF(column)(mmean,i);
		VECTORFF(scale)(&vv.vector,1/(VECTORFF(get)(vratio,i)));
	}
	VECTORFF(scale)(vratio,1/(FTYPE)ns);
	//Covariance
	BLASF(gemv)(CblasNoTrans,1,t2,t1,0,vcov);
	VECTORFF(scale)(vcov,1/(FTYPE)ns);

	//Buffer vb2=f_alpha mu_{alpha 1}
	VECTORFF(memcpy)(vb2,vratio);
	VECTORFF(mul)(vb2,vmean1);
	//vb3=(rho-sum_alpha f_alpha mu_{alpha 1}mu_{alpha 2})^2
	VECTORFF(memcpy)(vb3,vcov);
	BLASF(gemv)(CblasNoTrans,-1,mmean,vb2,1,vb3);
	VECTORFF(mul)(vb3,vb3);
	//fb1=1-sum_alpha f_alpha mu_{alpha 1}^2
	BLASF(dot)(vmean1,vb2,&fb1);
	fb1=1-fb1;
	//vcov=log(1-rho^2)
	VECTORFF(mul)(vcov,vcov);
	VECTORFF(scale)(vcov,-1);
	VECTORFF(add_constant)(vcov,1);
	for(i=0;i<nt;i++)
		VECTORFF(set)(vcov,i,(FTYPE)log(VECTORFF(get)(vcov,i)));
	//llr3=log(1-sum_alpha f_alpha mu_{alpha 1}^2)+log(1-rho^2)
	VECTORFF(set_all)(llr3,(FTYPE)log(fb1));
	VECTORFF(add)(llr3,vcov);
	
	//vcov=1-sum_alpha f_alpha mu_{alpha 2}^2
	MATRIXFF(mul_elements)(mmean,mmean);
	VECTORFF(set_all)(vcov,1);
	BLASF(gemv)(CblasNoTrans,-1,mmean,vratio,1,vcov);
	
	//vcov=log((1-sum_alpha f_alpha mu_{alpha 1}^2)(1-sum_alpha f_alpha mu_{alpha 2}^2)-(rho-sum_alpha f_alpha mu_{alpha 1}mu_{alpha 2})^2)
	VECTORFF(scale)(vcov,fb1);
	VECTORFF(sub)(vcov,vb3);
	for(i=0;i<nt;i++)
		VECTORFF(set)(vcov,i,(FTYPE)log(VECTORFF(get)(vcov,i)));
	//Final llr3
	VECTORFF(sub)(llr3,vcov);
	VECTORFF(scale)(llr3,-0.5);
	VECTORFF(bound_below)(llr3,0);
}

void pij_gassist_llr3_buffed(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,size_t nv,MATRIXF* llr3,const VECTORF* vb1,MATRIXF *mratio,MATRIXF *mmean1,MATRIXF **mmean2,MATRIXF *mcov,MATRIXF **mmb1,MATRIXF **mmb2)
{
	size_t	i,j;
	size_t	nt=t2->size1;
	size_t	ng=g->size1;
	FTYPE	f1;
	VECTORFF(view)	vv,vv2,vv3;
	
	//Calculate ratio and mean
	pij_gassist_llr_ratioandmean_buffed(g,t,t2,mratio,mmean1,mmean2,nv,mmb2,vb1);

	//mmb1=f_{alpha i}mu_{alpha ij}
	for(i=0;i<nv;i++)
	{
		MATRIXFF(memcpy)(mmb1[i],mmean2[i]);
		for(j=0;j<ng;j++)
		{
			f1=MATRIXFF(get)(mratio,i,j);
			vv=MATRIXFF(row)(mmb1[i],j);
			VECTORFF(scale)(&vv.vector,f1);
		}
	}
	
	//llr3=1-sum_alpha f_{alpha i}mu_{alpha ij}^2, use mcov as buffer
	MATRIXFF(set_all)(llr3,1);
	for(i=0;i<nv;i++)
	{
		MATRIXFF(memcpy)(mcov,mmb1[i]);
		MATRIXFF(mul_elements)(mcov,mmean2[i]);
		MATRIXFF(sub)(llr3,mcov);
	}
	
	//Calculate variance
	MATRIXFF(cov2_bounded)(t,t2,mcov);

	//mmb1[i][j,k]=mmb1[i][j,k]*mmean1[i][j]=f_{ij}mu_{ijj}mu_{ijk}
	for(i=0;i<nv;i++)
		for(j=0;j<ng;j++)
		{
			vv=MATRIXFF(row)(mmb1[i],j);
			VECTORFF(scale)(&vv.vector,MATRIXFF(get)(mmean1,i,j));
		}
	//mmb1[0]=sum_alpha f_{alpha i}mu_{alpha ii}mu_{alpha ij}
	for(i=1;i<nv;i++)
		MATRIXFF(add)(mmb1[0],mmb1[i]);
	//mmb1[0]=(rho_{ij}-sum_alpha f_{alpha i}mu_{alpha ii}mu_{alpha ij})^2
	MATRIXFF(sub)(mmb1[0],mcov);
	MATRIXFF(mul_elements)(mmb1[0],mmb1[0]);
	
	//mmb1[1][:,0]=1-sum_alpha f_{alpha i}mu_{alpha ii}^2
	vv=MATRIXFF(column)(mmb1[1],0);
	VECTORFF(set_zero)(&vv.vector);
	vv2=MATRIXFF(column)(mmb1[1],1);
	for(i=0;i<nv;i++)
	{
		vv3=MATRIXFF(row)(mratio,i);
		VECTORFF(memcpy)(&vv2.vector,&vv3.vector);
		vv3=MATRIXFF(row)(mmean1,i);
		VECTORFF(mul)(&vv2.vector,&vv3.vector);
		VECTORFF(mul)(&vv2.vector,&vv3.vector);
		VECTORFF(add)(&vv.vector,&vv2.vector);
	}
	VECTORFF(scale)(&vv.vector,-1);
	VECTORFF(add_constant)(&vv.vector,1);
	
	//mcov=log(1-rho_{ij}^2)
	MATRIXFF(mul_elements)(mcov,mcov);
	MATRIXFF(scale)(mcov,-1);
	MATRIXFF(add_constant)(mcov,1);
	for(i=0;i<ng;i++)
		for(j=0;j<nt;j++)
			MATRIXFF(set)(mcov,i,j,(FTYPE)log(MATRIXFF(get)(mcov,i,j)));
	
	//mcov=log(1-rho_{ij}^2)+log(1-sum_alpha f_{alpha i}mu_{alpha ii}^2)
	for(i=0;i<ng;i++)
	{
		f1=(FTYPE)log(MATRIXFF(get)(mmb1[1],i,0));
		vv=MATRIXFF(row)(mcov,i);
		VECTORFF(add_constant)(&vv.vector,f1);
	}
	
	//llr3=(1-sum_alpha f_{alpha i}mu_{alpha ii}^2)(1-sum_alpha f_{alpha i}mu_{alpha ij}^2)
	for(i=0;i<ng;i++)
	{
		f1=MATRIXFF(get)(mmb1[1],i,0);
		vv=MATRIXFF(row)(llr3,i);
		VECTORFF(scale)(&vv.vector,f1);
	}
	
	//Final steps
	MATRIXFF(sub)(llr3,mmb1[0]);
	for(i=0;i<ng;i++)
		for(j=0;j<nt;j++)
			MATRIXFF(set)(llr3,i,j,(FTYPE)log(MATRIXFF(get)(llr3,i,j)));
	MATRIXFF(sub)(llr3,mcov);
	MATRIXFF(scale)(llr3,0.5);
	MATRIXFF(bound_above)(llr3,0);
}

int pij_gassist_llr(const MATRIXG* g,const MATRIXF* t,const MATRIXF* t2,VECTORF* llr1,MATRIXF* llr2b,MATRIXF* llr2c,MATRIXF* llr3,size_t nv)
{
#define	CLEANUP			CLEANVECF(vbuff1)
	VECTORF	*vbuff1;		//(ns) Const buffer, set to all 1.
	int		ret;
	//struct pij_gassist_llr_threadinfo tinfo;
#ifndef NDEBUG
	size_t	ng,nt,ns;

	ng=g->size1;
	nt=t2->size1;
	ns=t->size2;
#endif

	vbuff1=0;
	//Validation
	assert(!((g->size2!=ns)||(t2->size2!=ns)||(t->size1!=ng)
		||(llr2b->size1!=ng)||(llr2c->size1!=ng)||(llr2b->size2!=nt)
		||(llr2c->size2!=nt)||(llr3->size1!=ng)||(llr3->size2!=nt)));
	assert(!(llr1->size!=ng));
	assert(!(nv>CONST_NV_MAX));
	assert(nv>MATRIXGF(max)(g));
	
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
			MATRIXFF(view)	mvllr2b,mvllr2c,mvllr3;
			vvllr1=VECTORFF(subvector)(llr1,n1,n2-n1);
			mvllr2b=MATRIXFF(submatrix)(llr2b,n1,0,n2-n1,llr2b->size2);
			mvllr2c=MATRIXFF(submatrix)(llr2c,n1,0,n2-n1,llr2c->size2);
			mvllr3=MATRIXFF(submatrix)(llr3,n1,0,n2-n1,llr3->size2);
			retth=pij_gassist_llr_block(&mvg.matrix,&mvt.matrix,t2,nv,&vvllr1.vector,&mvllr2b.matrix,&mvllr2c.matrix,&mvllr3.matrix,vbuff1);
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

