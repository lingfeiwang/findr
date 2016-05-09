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
/* This file contains functions to be called by R
 * Function names are self explanatory.
 */

#include "../base/config.h"
#include <stdlib.h>
#include <stdio.h>
#include "../base/logger.h"
#include "../base/random.h"
#include "../base/macros.h"
#include "../base/lib.h"
#include "../pij/gassist/gassist.h"
#include "../pij/rank.h"

void external_R_lib_init(const int *loglv,const int *rs0,const int *nthread)
{
	lib_init((unsigned char)(*loglv),(unsigned long)(*rs0),(size_t)(*nthread));
}

void external_R_lib_name(const char** ans)
{
	*ans=lib_name();
}

void external_R_lib_version(const char** ans)
{
	*ans=lib_version();
}

void external_R_pijs_gassist_any(const int *ng,const int *nt,const int *ns,const int* g,const double* t,const double* t2,double* p1,double* p2b,double* p2c,double* p3,const int* nv,const int* nodiag,int *ret,int (*func)(const MATRIXG*,const MATRIXF*,const MATRIXF*,VECTORF*,MATRIXF*,MATRIXF*,MATRIXF*,size_t,char))
{
#define	CLEANUP	CLEANMATG(mg)CLEANMATF(mt)CLEANMATF(mt2)CLEANVECF(vp1)\
				CLEANMATF(mp2b)CLEANMATF(mp2c)CLEANMATF(mp3)
	size_t	i,j;
	char	nd=(char)(*nodiag);
	size_t	ngv,ntv,nsv,nvv;
	ngv=(size_t)*ng;
	ntv=(size_t)*nt;
	nsv=(size_t)*ns;
	nvv=(size_t)*nv;
	MATRIXG *mg;
	MATRIXF	*mt,*mt2,*mp2b,*mp2c,*mp3;
	VECTORF	*vp1;
	
	//Construct GTYPE matrix for g
	mg=MATRIXGF(alloc)(ngv,nsv);
	mt=MATRIXFF(alloc)(ngv,nsv);
	mt2=MATRIXFF(alloc)(ntv,nsv);
	vp1=VECTORFF(alloc)(ngv);
	mp2b=MATRIXFF(alloc)(ngv,ntv);
	mp2c=MATRIXFF(alloc)(ngv,ntv);
	mp3=MATRIXFF(alloc)(ngv,ntv);
	if(!(mg&&mt&&mt2&&vp1&&mp2b&&mp2c&&mp3))
	{
		LOG(1,"Not enough memory.")
		CLEANUP
		*ret=1;
		return;
	}
	
	//Copy data, R uses column major
	for(i=0;i<ngv;i++)
		for(j=0;j<nsv;j++)
		{
			MATRIXGF(set)(mg,i,j,(GTYPE)(g[j*ngv+i]));
			MATRIXFF(set)(mt,i,j,(FTYPE)(t[j*ngv+i]));
		}	
	for(i=0;i<ntv;i++)
		for(j=0;j<nsv;j++)
			MATRIXFF(set)(mt2,i,j,(FTYPE)(t2[j*ntv+i]));
	
	//Calculation
	*ret=func(mg,mt,mt2,vp1,mp2b,mp2c,mp3,nvv,nd);
	//Copy data back
	if(!*ret)
	{
		for(i=0;i<ngv;i++)
		{
			p1[i]=(double)VECTORFF(get)(vp1,i);
			for(j=0;j<ntv;j++)
			{
				p2b[i*ngv+j]=(double)MATRIXFF(get)(mp2b,j,i);
				p2c[i*ngv+j]=(double)MATRIXFF(get)(mp2c,j,i);
				p3[i*ngv+j]=(double)MATRIXFF(get)(mp3,j,i);
			}
		}
	}
	CLEANUP
#undef CLEANUP
}


void external_R_pijs_gassist_a(const int *ng,const int *nt,const int *ns,const int* g,const double* t,const double* t2,double* p1,double* p2b,double* p2c,double* p3,const int* nv,const int* nodiag,int *ret)
{
	external_R_pijs_gassist_any(ng,nt,ns,g,t,t2,p1,p2b,p2c,p3,nv,nodiag,ret,pijs_gassist_a);
}

void external_R_pijs_gassist_tot(const int *ng,const int *nt,const int *ns,const int* g,const double* t,const double* t2,double* p1,double* p2b,double* p2c,double* p3,const int* nv,const int* nodiag,int *ret)
{
	external_R_pijs_gassist_any(ng,nt,ns,g,t,t2,p1,p2b,p2c,p3,nv,nodiag,ret,pijs_gassist_tot);
}

void external_R_pij_rank_a(const int *ng,const int *nt,const int *ns,const double* t,const double* t2,double* p,const int* nodiag,int *ret)
{
#define	CLEANUP	CLEANMATF(mt)CLEANMATF(mt2)CLEANMATF(mp)
	size_t	i,j;
	char	nd=(char)(*nodiag);
	size_t	ngv,ntv,nsv;
	ngv=(size_t)*ng;
	ntv=(size_t)*nt;
	nsv=(size_t)*ns;
	MATRIXF	*mt,*mt2,*mp;
	
	mt=MATRIXFF(alloc)(ngv,nsv);
	mt2=MATRIXFF(alloc)(ntv,nsv);
	mp=MATRIXFF(alloc)(ngv,ntv);
	if(!(mt&&mt2&&mp))
	{
		LOG(1,"Not enough memory.")
		CLEANUP
		*ret=1;
		return;
	}

	//Copy data, R uses column major
	for(i=0;i<ngv;i++)
		for(j=0;j<nsv;j++)
			MATRIXFF(set)(mt,i,j,(FTYPE)(t[j*ngv+i]));
	for(i=0;i<ntv;i++)
		for(j=0;j<nsv;j++)
			MATRIXFF(set)(mt2,i,j,(FTYPE)(t2[j*ntv+i]));
	
	//Calculation
	*ret=pij_rank_a(mt,mt2,mp,nd);
	//Copy data back
	if(!*ret)
		for(i=0;i<ngv;i++)
			for(j=0;j<ntv;j++)
				p[i*ngv+j]=(double)MATRIXFF(get)(mp,j,i);
	CLEANUP
#undef CLEANUP
}





























