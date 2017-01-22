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

void external_R_pijs_gassist_any(const int *ng,const int *nt,const int *ns,const int* g,const double* t,const double* t2,double* p1,double* p2,double* p3,double* p4,double* p5,const int* nv,const int* nodiag,int *ret,int (*func)(const MATRIXG*,const MATRIXF*,const MATRIXF*,VECTORF*,MATRIXF*,MATRIXF*,MATRIXF*,MATRIXF*,size_t,char,size_t))
{
#define	CLEANUP	CLEANMATG(mg)CLEANMATF(mt)CLEANMATF(mt2)CLEANVECF(vp1)\
				CLEANMATF(mp2)CLEANMATF(mp3)CLEANMATF(mp4)CLEANMATF(mp5)
	size_t	i,j;
	char	nd=(char)(*nodiag);
	size_t	ngv,ntv,nsv,nvv;
	ngv=(size_t)*ng;
	ntv=(size_t)*nt;
	nsv=(size_t)*ns;
	nvv=(size_t)*nv;
	MATRIXG *mg;
	MATRIXF	*mt,*mt2,*mp2,*mp3,*mp4,*mp5;
	VECTORF	*vp1;
	
	//Construct GTYPE matrix for g
	mg=MATRIXGF(alloc)(ngv,nsv);
	mt=MATRIXFF(alloc)(ngv,nsv);
	mt2=MATRIXFF(alloc)(ntv,nsv);
	vp1=VECTORFF(alloc)(ngv);
	mp2=MATRIXFF(alloc)(ngv,ntv);
	mp3=MATRIXFF(alloc)(ngv,ntv);
	mp4=MATRIXFF(alloc)(ngv,ntv);
	mp5=MATRIXFF(alloc)(ngv,ntv);
	if(!(mg&&mt&&mt2&&vp1&&mp2&&mp3&&mp4&&mp5))
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
	*ret=func(mg,mt,mt2,vp1,mp2,mp3,mp4,mp5,nvv,nd,(size_t)-1);
	//Copy data back
	if(!*ret)
	{
		for(i=0;i<ngv;i++)
		{
			p1[i]=(double)VECTORFF(get)(vp1,i);
			for(j=0;j<ntv;j++)
			{
				p2[j*ngv+i]=(double)MATRIXFF(get)(mp2,i,j);
				p3[j*ngv+i]=(double)MATRIXFF(get)(mp3,i,j);
				p4[j*ngv+i]=(double)MATRIXFF(get)(mp4,i,j);
				p5[j*ngv+i]=(double)MATRIXFF(get)(mp5,i,j);
			}
		}
	}
	CLEANUP
#undef CLEANUP
}


void external_R_pijs_gassist(const int *ng,const int *nt,const int *ns,const int* g,const double* t,const double* t2,double* p1,double* p2,double* p3,double* p4,double* p5,const int* nv,const int* nodiag,int *ret)
{
	external_R_pijs_gassist_any(ng,nt,ns,g,t,t2,p1,p2,p3,p4,p5,nv,nodiag,ret,pijs_gassist);
}

void external_R_pij_gassist_any(const int *ng,const int *nt,const int *ns,const int* g,const double* t,const double* t2,double* p,const int* nv,const int* nodiag,int *ret,int (*func)(const MATRIXG*,const MATRIXF*,const MATRIXF*,MATRIXF*,size_t,char,size_t))
{
#define	CLEANUP	CLEANMATG(mg)CLEANMATF(mt)CLEANMATF(mt2)CLEANMATF(mp)
	size_t	i,j;
	char	nd=(char)(*nodiag);
	size_t	ngv,ntv,nsv,nvv;
	ngv=(size_t)*ng;
	ntv=(size_t)*nt;
	nsv=(size_t)*ns;
	nvv=(size_t)*nv;
	MATRIXG *mg;
	MATRIXF	*mt,*mt2,*mp;
	
	//Construct GTYPE matrix for g
	mg=MATRIXGF(alloc)(ngv,nsv);
	mt=MATRIXFF(alloc)(ngv,nsv);
	mt2=MATRIXFF(alloc)(ntv,nsv);
	mp=MATRIXFF(alloc)(ngv,ntv);
	if(!(mg&&mt&&mt2&&mp))
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
	*ret=func(mg,mt,mt2,mp,nvv,nd,(size_t)-1);
	//Copy data back
	if(!*ret)
	{
		for(i=0;i<ngv;i++)
			for(j=0;j<ntv;j++)
				p[j*ngv+i]=(double)MATRIXFF(get)(mp,i,j);
	}
	CLEANUP
#undef CLEANUP
}

void external_R_pij_gassist(const int *ng,const int *nt,const int *ns,const int* g,const double* t,const double* t2,double* p,const int* nv,const int* nodiag,int *ret)
{
	external_R_pij_gassist_any(ng,nt,ns,g,t,t2,p,nv,nodiag,ret,pij_gassist);
}

void external_R_pij_gassist_trad(const int *ng,const int *nt,const int *ns,const int* g,const double* t,const double* t2,double* p,const int* nv,const int* nodiag,int *ret)
{
	external_R_pij_gassist_any(ng,nt,ns,g,t,t2,p,nv,nodiag,ret,pij_gassist_trad);
}

void external_R_pij_rank(const int *ng,const int *nt,const int *ns,const double* t,const double* t2,double* p,const int* nodiag,int *ret)
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
	*ret=pij_rank(mt,mt2,mp,nd,(size_t)-1);
	//Copy data back
	if(!*ret)
		for(i=0;i<ngv;i++)
			for(j=0;j<ntv;j++)
				p[j*ngv+i]=(double)MATRIXFF(get)(mp,i,j);
	CLEANUP
#undef CLEANUP
}





























