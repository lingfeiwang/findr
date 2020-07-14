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


/* Calculates the 5 log likelihood ratios of Trigger with nonpermuted data in block form:
 * 1. E->A v.s. E no relation with A
 * 2. A<-E->B with A--B v.s. E->A<-B
 * 3. E->A->B v.s. A<-E->B with A--B
 * 4. A<-E->B with A->B v.s. E->A
 * 5. A<-E->B with A->B v.s. A<-E->B
 * Uses GSL BLAS.
 * Note: for each row, g must be the best eQTL of t of the same row.
 */
static void pij_cassist_llr_block(const MATRIXF* g,const MATRIXF* t,const MATRIXF* t2,VECTORF* llr1,MATRIXF* llr2,MATRIXF* llr3,MATRIXF* llr4,MATRIXF* llr5)
{
	size_t	i,j;
	VECTORFF(view)	vv;
	size_t	ng=g->size1;
	size_t	nt=t2->size1;
	
	assert(ng&&(t->size1==ng)&&(llr1->size==ng)&&(llr2->size1==ng)&&(llr3->size1==ng)
		&&(llr4->size1==ng)&&(llr5->size1==ng));
	assert(nt&&(llr2->size2==nt)&&(llr3->size2==nt)&&(llr4->size2==nt)&&(llr5->size2==nt));
	assert(g->size2&&(t->size2==g->size2)&&(t2->size2==g->size2));
	
	//llr1=rho_EA
	MATRIXFF(cov2_1v1_bounded)(g,t,llr1);
	//llr2=rho_EB
	MATRIXFF(cov2_bounded)(g,t2,llr2);
	//llr5=rho_AB
	MATRIXFF(cov2_bounded)(t,t2,llr5);

	//llr4=rho_EA*rho_EB
	MATRIXFF(memcpy)(llr4,llr2);
	for(i=0;i<ng;i++)
	{
		vv=MATRIXFF(row)(llr4,i);
		VECTORFF(scale)(&vv.vector,VECTORFF(get)(llr1,i));
	}
	//llr4=(rho_AB-rho_EA*rho_EB)^2
	MATRIXFF(sub)(llr4,llr5);
	MATRIXFF(mul_elements)(llr4,llr4);
	//llr2=1-rho_EB^2
	MATRIXFF(mul_elements)(llr2,llr2);
	MATRIXFF(scale)(llr2,-1);
	MATRIXFF(add_constant)(llr2,1);
	//llr5=1-rho_AB^2
	MATRIXFF(mul_elements)(llr5,llr5);
	MATRIXFF(scale)(llr5,-1);
	MATRIXFF(add_constant)(llr5,1);
	//llr1=1-rho_EA^2
	VECTORFF(mul)(llr1,llr1);
	VECTORFF(scale)(llr1,-1);
	VECTORFF(add_constant)(llr1,1);
	//llr4=(1-rho_EA^2)*(1-rho_EB^2)-(rho_AB-rho_EA*rho_EB)^2, using llr3 as buffer
	MATRIXFF(memcpy)(llr3,llr2);
	for(i=0;i<ng;i++)
	{
		vv=MATRIXFF(row)(llr3,i);
		VECTORFF(scale)(&vv.vector,VECTORFF(get)(llr1,i));
	}
	MATRIXFF(sub)(llr4,llr3);
	MATRIXFF(scale)(llr4,-1);
	MATRIXFF(bound_below)(llr4,FTYPE_MIN);
	
	//ALL log
	for(i=0;i<ng;i++)
	{
		VECTORFF(set)(llr1,i,(FTYPE)log(VECTORFF(get)(llr1,i)));
		for(j=0;j<nt;j++)
		{
			MATRIXFF(set)(llr2,i,j,(FTYPE)log(MATRIXFF(get)(llr2,i,j)));
			MATRIXFF(set)(llr4,i,j,(FTYPE)log(MATRIXFF(get)(llr4,i,j)));
			MATRIXFF(set)(llr5,i,j,(FTYPE)log(MATRIXFF(get)(llr5,i,j)));
		}
	}
	//llr4=llr4 before scaling -0.5
	for(i=0;i<ng;i++)
	{
		FTYPE v;
		v=VECTORFF(get)(llr1,i);
		vv=MATRIXFF(row)(llr4,i);
		VECTORFF(add_constant)(&vv.vector,-v);
	}
	//llr3=llr3 before scaling -0.5
	MATRIXFF(memcpy)(llr3,llr4);
	MATRIXFF(sub)(llr3,llr5);
	//llr5=llr5 before scaling -0.5
	MATRIXFF(memcpy)(llr5,llr4);
	MATRIXFF(sub)(llr5,llr2);

	//Finalize all for coefficients
	VECTORFF(scale)(llr1,-0.5);
	MATRIXFF(scale)(llr2,-0.5);
	MATRIXFF(scale)(llr3,-0.5);
	MATRIXFF(scale)(llr4,-0.5);
	MATRIXFF(scale)(llr5,-0.5);
	
	//Bounding from 0
	VECTORFF(bound_below)(llr1,0);
	MATRIXFF(bound_below)(llr2,0);
	MATRIXFF(bound_below)(llr3,0);
	MATRIXFF(bound_below)(llr4,0);
	MATRIXFF(bound_below)(llr5,0);
}

void pij_cassist_llr(const MATRIXF* g,const MATRIXF* t,const MATRIXF* t2,VECTORF* llr1,MATRIXF* llr2,MATRIXF* llr3,MATRIXF* llr4,MATRIXF* llr5)
{
#ifndef NDEBUG
	size_t	ng,nt,ns;

	ng=g->size1;
	nt=t2->size1;
	ns=t->size2;
#endif

	//Validation
	assert(!((g->size2!=ns)||(t2->size2!=ns)||(t->size1!=ng)||(llr2->size1!=ng)||(llr2->size2!=nt)||(llr3->size1!=ng)||(llr3->size2!=nt)||(llr4->size1!=ng)||(llr4->size2!=nt)||(llr5->size1!=ng)||(llr5->size2!=nt)));
	assert(!(llr1->size!=ng));
	
	#pragma omp parallel
	{
		size_t	n1,n2;
		threading_get_startend(t->size1,&n1,&n2);
		if(n2>n1)
		{
			MATRIXFF(const_view) mvg=MATRIXFF(const_submatrix)(g,n1,0,n2-n1,g->size2);
			MATRIXFF(const_view) mvt=MATRIXFF(const_submatrix)(t,n1,0,n2-n1,t->size2);
			VECTORFF(view)	vvllr1;
			MATRIXFF(view)	mvllr2,mvllr3,mvllr4,mvllr5;
			vvllr1=VECTORFF(subvector)(llr1,n1,n2-n1);
			mvllr2=MATRIXFF(submatrix)(llr2,n1,0,n2-n1,llr2->size2);
			mvllr3=MATRIXFF(submatrix)(llr3,n1,0,n2-n1,llr3->size2);
			mvllr4=MATRIXFF(submatrix)(llr4,n1,0,n2-n1,llr4->size2);
			mvllr5=MATRIXFF(submatrix)(llr5,n1,0,n2-n1,llr5->size2);
			pij_cassist_llr_block(&mvg.matrix,&mvt.matrix,t2,&vvllr1.vector,&mvllr2.matrix,&mvllr3.matrix,&mvllr4.matrix,&mvllr5.matrix);
		}
	}
}

