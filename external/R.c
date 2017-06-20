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
#include "../pij/cassist/cassist.h"
#include "../pij/rank.h"
#include "../netr/one.h"

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

void external_R_pijs_gassist_pv(const int *ng,const int *nt,const int *ns,const int* g,const double* t,const double* t2,double* p1,double* p2,double* p3,double* p4,double* p5,const int* nv,int *ret)
{
#define	CLEANUP	CLEANMATG(mg)CLEANMATF(mt)CLEANMATF(mt2)CLEANVECF(vp1)\
				CLEANMATF(mp2)CLEANMATF(mp3)CLEANMATF(mp4)CLEANMATF(mp5)
	size_t	i,j;
	size_t	ngv,ntv,nsv,nvv;
	ngv=(size_t)*ng;
	ntv=(size_t)*nt;
	nsv=(size_t)*ns;
	nvv=(size_t)*nv;
	MATRIXG *mg;
	MATRIXF	*mt,*mt2,*mp2,*mp3,*mp4,*mp5;
	VECTORF	*vp1;
	
	LOG(12,"R interface for external_R_pijs_gassist_pv: nt=%i, nt2=%i, ns=%i, nv=%i",*ng,*nt,*ns,*nv)
	
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
	
	*ret=pijs_gassist_pv(mg,mt,mt2,vp1,mp2,mp3,mp4,mp5,nvv,(size_t)-1);
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
	LOG(12,"R interface for external_R_pijs_gassist: nt=%i, nt2=%i, ns=%i, nv=%i, nodiag=%i",*ng,*nt,*ns,*nv,*nodiag)
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
	LOG(12,"R interface for external_R_pij_gassist: nt=%i, nt2=%i, ns=%i, nv=%i, nodiag=%i",*ng,*nt,*ns,*nv,*nodiag)
	external_R_pij_gassist_any(ng,nt,ns,g,t,t2,p,nv,nodiag,ret,pij_gassist);
}

void external_R_pij_gassist_trad(const int *ng,const int *nt,const int *ns,const int* g,const double* t,const double* t2,double* p,const int* nv,const int* nodiag,int *ret)
{
	LOG(12,"R interface for external_R_pij_gassist_trad: nt=%i, nt2=%i, ns=%i, nv=%i, nodiag=%i",*ng,*nt,*ns,*nv,*nodiag)
	external_R_pij_gassist_any(ng,nt,ns,g,t,t2,p,nv,nodiag,ret,pij_gassist_trad);
}

void external_R_pijs_cassist_pv(const int *ng,const int *nt,const int *ns,const double* g,const double* t,const double* t2,double* p1,double* p2,double* p3,double* p4,double* p5,int *ret)
{
#define	CLEANUP	CLEANMATF(mg)CLEANMATF(mt)CLEANMATF(mt2)CLEANVECF(vp1)\
				CLEANMATF(mp2)CLEANMATF(mp3)CLEANMATF(mp4)CLEANMATF(mp5)
	size_t	i,j;
	size_t	ngv,ntv,nsv;
	ngv=(size_t)*ng;
	ntv=(size_t)*nt;
	nsv=(size_t)*ns;
	MATRIXF *mg;
	MATRIXF	*mt,*mt2,*mp2,*mp3,*mp4,*mp5;
	VECTORF	*vp1;
	
	LOG(12,"R interface for external_R_pijs_cassist_pv: nt=%i, nt2=%i, ns=%i",*ng,*nt,*ns)
	
	mg=MATRIXFF(alloc)(ngv,nsv);
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
			MATRIXFF(set)(mg,i,j,(FTYPE)(g[j*ngv+i]));
			MATRIXFF(set)(mt,i,j,(FTYPE)(t[j*ngv+i]));
		}	
	for(i=0;i<ntv;i++)
		for(j=0;j<nsv;j++)
			MATRIXFF(set)(mt2,i,j,(FTYPE)(t2[j*ntv+i]));
	
	//Calculation
	
	*ret=pijs_cassist_pv(mg,mt,mt2,vp1,mp2,mp3,mp4,mp5,(size_t)-1);
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

void external_R_pijs_cassist_any(const int *ng,const int *nt,const int *ns,const double* g,const double* t,const double* t2,double* p1,double* p2,double* p3,double* p4,double* p5,const int* nodiag,int *ret,int (*func)(const MATRIXF*,const MATRIXF*,const MATRIXF*,VECTORF*,MATRIXF*,MATRIXF*,MATRIXF*,MATRIXF*,char,size_t))
{
#define	CLEANUP	CLEANMATF(mg)CLEANMATF(mt)CLEANMATF(mt2)CLEANVECF(vp1)\
				CLEANMATF(mp2)CLEANMATF(mp3)CLEANMATF(mp4)CLEANMATF(mp5)
	size_t	i,j;
	char	nd=(char)(*nodiag);
	size_t	ngv,ntv,nsv;
	ngv=(size_t)*ng;
	ntv=(size_t)*nt;
	nsv=(size_t)*ns;
	MATRIXF *mg;
	MATRIXF	*mt,*mt2,*mp2,*mp3,*mp4,*mp5;
	VECTORF	*vp1;
	
	mg=MATRIXFF(alloc)(ngv,nsv);
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
			MATRIXFF(set)(mg,i,j,(FTYPE)(g[j*ngv+i]));
			MATRIXFF(set)(mt,i,j,(FTYPE)(t[j*ngv+i]));
		}	
	for(i=0;i<ntv;i++)
		for(j=0;j<nsv;j++)
			MATRIXFF(set)(mt2,i,j,(FTYPE)(t2[j*ntv+i]));
	
	//Calculation
	*ret=func(mg,mt,mt2,vp1,mp2,mp3,mp4,mp5,nd,(size_t)-1);
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

void external_R_pijs_cassist(const int *ng,const int *nt,const int *ns,const double* g,const double* t,const double* t2,double* p1,double* p2,double* p3,double* p4,double* p5,const int* nodiag,int *ret)
{
	LOG(12,"R interface for external_R_pijs_cassist: nt=%i, nt2=%i, ns=%i, nodiag=%i",*ng,*nt,*ns,*nodiag)
	external_R_pijs_cassist_any(ng,nt,ns,g,t,t2,p1,p2,p3,p4,p5,nodiag,ret,pijs_cassist);
}

void external_R_pij_cassist_any(const int *ng,const int *nt,const int *ns,const double* g,const double* t,const double* t2,double* p,const int* nodiag,int *ret,int (*func)(const MATRIXF*,const MATRIXF*,const MATRIXF*,MATRIXF*,char,size_t))
{
#define	CLEANUP	CLEANMATF(mg)CLEANMATF(mt)CLEANMATF(mt2)CLEANMATF(mp)
	size_t	i,j;
	char	nd=(char)(*nodiag);
	size_t	ngv,ntv,nsv;
	ngv=(size_t)*ng;
	ntv=(size_t)*nt;
	nsv=(size_t)*ns;
	MATRIXF *mg;
	MATRIXF	*mt,*mt2,*mp;
	
	mg=MATRIXFF(alloc)(ngv,nsv);
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
			MATRIXFF(set)(mg,i,j,(FTYPE)(g[j*ngv+i]));
			MATRIXFF(set)(mt,i,j,(FTYPE)(t[j*ngv+i]));
		}	
	for(i=0;i<ntv;i++)
		for(j=0;j<nsv;j++)
			MATRIXFF(set)(mt2,i,j,(FTYPE)(t2[j*ntv+i]));
	
	//Calculation
	*ret=func(mg,mt,mt2,mp,nd,(size_t)-1);
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

void external_R_pij_cassist(const int *ng,const int *nt,const int *ns,const double* g,const double* t,const double* t2,double* p,const int* nodiag,int *ret)
{
	LOG(12,"R interface for external_R_pij_cassist: nt=%i, nt2=%i, ns=%i, nodiag=%i",*ng,*nt,*ns,*nodiag)
	external_R_pij_cassist_any(ng,nt,ns,g,t,t2,p,nodiag,ret,pij_cassist);
}

void external_R_pij_cassist_trad(const int *ng,const int *nt,const int *ns,const double* g,const double* t,const double* t2,double* p,const int* nodiag,int *ret)
{
	LOG(12,"R interface for external_R_pij_cassist_trad: nt=%i, nt2=%i, ns=%i, nodiag=%i",*ng,*nt,*ns,*nodiag)
	external_R_pij_cassist_any(ng,nt,ns,g,t,t2,p,nodiag,ret,pij_cassist_trad);
}

void external_R_pij_rank_pv(const int *ng,const int *nt,const int *ns,const double* t,const double* t2,double* p,int *ret)
{
#define	CLEANUP	CLEANMATF(mt)CLEANMATF(mt2)CLEANMATF(mp)
	LOG(12,"R interface for external_R_pij_rank_pv: nt=%i, nt2=%i, ns=%i",*ng,*nt,*ns)
	size_t	i,j;
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
	*ret=pij_rank_pv(mt,mt2,mp,(size_t)-1);
	//Copy data back
	if(!*ret)
		for(i=0;i<ngv;i++)
			for(j=0;j<ntv;j++)
				p[j*ngv+i]=(double)MATRIXFF(get)(mp,i,j);
	CLEANUP
#undef CLEANUP
}

void external_R_pij_rank(const int *ng,const int *nt,const int *ns,const double* t,const double* t2,double* p,const int* nodiag,int *ret)
{
#define	CLEANUP	CLEANMATF(mt)CLEANMATF(mt2)CLEANMATF(mp)
	LOG(12,"R interface for external_R_pij_rank: nt=%i, nt2=%i, ns=%i, nodiag=%i",*ng,*nt,*ns,*nodiag)
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

void external_R_netr_one_greedy(const int *nt,const double* p,const int* namax0,const int* nimax0,const int* nomax0,int* net,int *ret)
{
#define	CLEANUP	CLEANMATF(mp)CLEANMATUC(mnet)
	LOG(12,"R interface for external_R_netr_one_greedy: nt=%i, namax=%i, nimax=%i, nomax=%i",*nt,*namax0,*nimax0,*nomax0)
	size_t	i,j;
	size_t	ntv=(size_t)*nt,ret2;
	size_t	namax,nimax,nomax;
	MATRIXF		*mp;
	MATRIXUC	*mnet;
	
	mp=MATRIXFF(alloc)(ntv,ntv);
	mnet=MATRIXUCF(alloc)(ntv,ntv);
	if(!(mp&&mnet))
	{
		LOG(1,"Not enough memory.")
		CLEANUP
		*ret=1;
		return;
	}
	
	namax=(size_t)(*namax0<=0?-1:*namax0);
	nimax=(size_t)(*nimax0<=0?-1:*nimax0);
	nomax=(size_t)(*nomax0<=0?-1:*nomax0);
	
	//Copy data, R uses column major
	for(i=0;i<ntv;i++)
		for(j=0;j<ntv;j++)
			MATRIXFF(set)(mp,i,j,(FTYPE)(p[j*ntv+i]));
	
	//Calculation
	ret2=netr_one_greedy(mp,mnet,namax,nimax,nomax);
	*ret=(ret2==0);
	//Copy data back
	if(!*ret)
		for(i=0;i<ntv;i++)
			for(j=0;j<ntv;j++)
				net[j*ntv+i]=(int)MATRIXUCF(get)(mnet,i,j);
	CLEANUP
#undef CLEANUP
}


























