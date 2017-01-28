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
#include <assert.h>
#include <string.h>
#include "../base/general_alg.h"
#include "../base/macros.h"
#include "vg.h"

int cycle_vg_init(struct cycle_vg_system* restrict vg,size_t dim,size_t amax)
{
	int	ret;
	assert(vg);
	ret=0;
	vg->n=dim;
	vg->nam=amax;
	vg->nim=vg->nom=(size_t)-1;

	MALLOCSIZE(vg->lvf,dim);
	MALLOCSIZE(vg->lvb,dim);
	MALLOCSIZE(vg->go,dim);
	MALLOCSIZE(vg->goi,dim);
	ret=ret||data_ll_init(&vg->gao,amax)||data_ll_init(&vg->gai,amax);
	MALLOCSIZE(vg->gaof,dim);
	MALLOCSIZE(vg->gaif,dim);
	MALLOCSIZE(vg->gni,dim);
	MALLOCSIZE(vg->gno,dim);
	MALLOCSIZE(vg->lao,dim);
	MALLOCSIZE(vg->lai,dim);
	MALLOCSIZE(vg->buff,dim);
	MALLOCSIZE(vg->buff2,dim);
	ret=ret||data_heap_init(&vg->lvfl,dim)||data_heapdec_init(&vg->lvbl,dim);
	if(ret||!(vg->lvf&&vg->lvb&&vg->go&&vg->goi&&vg->gaof&&vg->gaif&&vg->gni&&vg->gno&&vg->lao&&vg->lai&&vg->buff&&vg->buff2))
	{
		cycle_vg_free(vg);
		LOG(1,"Not enough memory.")
		return 1;
	}
	//Initialize for empty graph.
	if(cycle_vg_empty(vg))
	{
		cycle_vg_free(vg);
		return 1;
	}
	return 0;
}

struct cycle_vg_system* cycle_vg_new(size_t dim,size_t amax)
{
	struct cycle_vg_system* vg;
	MALLOCSIZE(vg,1);
	if(!vg)
	{
		LOG(1,"Not enough memory.")
		return 0;
	}
	if(cycle_vg_init(vg,dim,amax))
	{
		free(vg);
		return 0;
	}
	return vg;
}

int cycle_vg_free(struct cycle_vg_system* restrict vg)
{
#define FREEMEM(X)	if(X){free(X);X=0;}
	FREEMEM(vg->lvf)
	FREEMEM(vg->lvb)
	FREEMEM(vg->go)
	FREEMEM(vg->goi)
	FREEMEM(vg->gaof)
	FREEMEM(vg->gaif)
	FREEMEM(vg->gni)
	FREEMEM(vg->gno)
	FREEMEM(vg->lao)
	FREEMEM(vg->lai)
	FREEMEM(vg->buff)
	FREEMEM(vg->buff2)
	data_ll_free(&vg->gao);
	data_ll_free(&vg->gai);
	data_heap_free(&vg->lvbl);
	data_heapdec_free(&vg->lvfl);
	return 0;
#undef FREEMEM
}

int cycle_vg_empty(struct cycle_vg_system* restrict vg)
{
	size_t i;
	vg->na=0;
	memset(vg->gaof,-1,vg->n*sizeof(*vg->gaof));
	memset(vg->gaif,-1,vg->n*sizeof(*vg->gaif));
	memset(vg->gni,0,vg->n*sizeof(*vg->gni));
	memset(vg->gno,0,vg->n*sizeof(*vg->gno));
	data_ll_empty(&vg->gao);
	data_ll_empty(&vg->gai);
	for(i=0;i<vg->n;i++)
	{
		vg->go[i]=i;
		vg->goi[i]=i;
	}
	return 0;
}

void cycle_vg_restore_order(struct cycle_vg_system* restrict vg,size_t vv)
{
	size_t	t;
	char	cond;
	size_t*	p[2];
	size_t* ps;
	
	t=vg->go[vv];
	cond=!!vg->lvfl.n;
	if(cond)
	{
		size_t t1;
		t1=data_heap_top(&vg->lvfl);
		if(t1<t)
			t=t1;
		else
			cond=0;
	}
	if(!cond)
	{
		p[0]=vg->buff;
		p[1]=vg->buff2;
		general_alg_categorize_embed(vg->goi,vg->lvf,p,t);
		*(p[0]++)=vg->goi[t];
		memcpy(p[0],vg->buff2,(size_t)(p[1]-vg->buff2)*sizeof(*p[0]));
		memcpy(vg->buff+t+1,vg->goi+t+1,(vg->n-t-1)*sizeof(*vg->buff));
		
		ps=vg->buff;
		vg->buff=vg->goi;
		vg->goi=ps;
		cycle_vg_fix_go(vg);
		return;
	}
	
	p[0]=vg->buff;
	p[1]=vg->buff2;
	general_alg_categorize_embed(vg->goi,vg->lvf,p,t);
	*(p[1]++)=vg->goi[t];
	ps=p[0];
	p[0]=p[1];
	p[1]=ps;
	general_alg_categorize_embed(vg->goi+t+1,vg->lvb,p,vg->n-t-1);
	memcpy(p[1],vg->buff2,(size_t)(p[0]-vg->buff2)*sizeof(*p[1]));

	ps=vg->buff;
	vg->buff=vg->goi;
	vg->goi=ps;
	cycle_vg_fix_go(vg);
	return;
}

int cycle_vg_add(struct cycle_vg_system* restrict vg,size_t v1,size_t v2)
{
	
	//Validity check
	assert(v1!=v2);
	if(vg->na>=vg->nam)
		return 1;
	if(vg->go[v1]<vg->go[v2])
		return cycle_vg_add_arc(vg,v1,v2);

	//Initialize
	data_heap_empty(&vg->lvfl);
	data_heapdec_empty(&vg->lvbl);
	memset(vg->lvf,0,vg->n*sizeof(*vg->lvf));  
	memset(vg->lvb,0,vg->n*sizeof(*vg->lvb));
	memset(vg->lao,-1,vg->n*sizeof(*vg->lao));
	memset(vg->lai,-1,vg->n*sizeof(*vg->lai));
	
	//Test loop
	//Enter function, line 1
	vg->lvf[v2]=1;
	vg->lvb[v1]=1;
	vg->lao[v2]=vg->gaof[v2];
	vg->lai[v1]=vg->gaif[v1];
	//line 2
	if(vg->gaof[v2]!=(size_t)-1)
		data_heap_push(&vg->lvfl,vg->go[v2]);
	//line 3
	if(vg->gaif[v1]!=(size_t)-1)
		data_heapdec_push(&vg->lvbl,vg->go[v1]);
	//line 4&5 (while)
	while((vg->lvfl.n>0)&&(vg->lvbl.n>0))
	{
		size_t vu,vx,vy,vz;
		
		vu=data_heap_top(&vg->lvfl);
		vz=data_heapdec_top(&vg->lvbl);
		if(vu>=vz)
			break;
		vu=vg->goi[vu];
		vz=vg->goi[vz];
		//Enter macro, line 1
		vx=data_ll_val(&vg->gao,vg->lao[vu]);
		vy=data_ll_val(&vg->gai,vg->lai[vz]);
		//line 2
		vg->lao[vu]=data_ll_child(&vg->gao,vg->lao[vu]);
		vg->lai[vz]=data_ll_child(&vg->gai,vg->lai[vz]);
		//line 3
		if(vg->lao[vu]==(size_t)-1)
			data_heap_pop(&vg->lvfl);
		if(vg->lai[vz]==(size_t)-1)
			data_heapdec_pop(&vg->lvbl);
		//line 4, first half
		if(vg->lvb[vx])
			return 1;
		//line 5-8 (if)
		if(!vg->lvf[vx])
		{
			//line 6,7
			vg->lvf[vx]=1;
			if(vg->gaof[vx]!=(size_t)-1)
			{
				vg->lao[vx]=vg->gaof[vx];
				data_heap_push(&vg->lvfl,vg->go[vx]);
			}
		}
		//line 4, second half
		if(vg->lvf[vy])
			return 1;
		//line 9-12 (if)
		if(!vg->lvb[vy])
		{
			//line 10,11
			vg->lvb[vy]=1;
			if(vg->gaif[vy]!=(size_t)-1)
			{
				vg->lai[vy]=vg->gaif[vy];
				data_heapdec_push(&vg->lvbl,vg->go[vy]);
			}
		}
	}
	
	//Add arc
	if(cycle_vg_add_arc(vg,v1,v2))
		return 1;
	
	//Recover ordering
	cycle_vg_restore_order(vg,v1);
	return 0;
}

int cycle_vg_test(const struct cycle_vg_system* restrict vg,size_t v1,size_t v2)
{
	LOG(0,"Not implemented.")
	return cycle_vg_add((struct cycle_vg_system*)vg,v1,v2);
}


void cycle_vg_extract_graph(const struct cycle_vg_system* restrict vg,MATRIXUC* g)
{
	size_t	i;
	size_t	t1;
	
	assert((vg->n==g->size1)&&(vg->n==g->size2));
	MATRIXUCF(set_zero)(g);
	for(i=0;i<vg->n;i++)
	{
		t1=vg->gaof[i];
		while(t1!=(size_t)-1)
		{
			MATRIXUCF(set)(g,i,data_ll_val(&vg->gao,t1),1);
			t1=data_ll_child(&vg->gao,t1);
		}
	}
}
