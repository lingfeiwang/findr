/* Copyright 2016-2018 Lingfei Wang
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
/* This lib contains the definitions of Vertex Guided Search and topological order maintenance.
 */
 
#ifndef _HEADER_LIB_CYCLE_VG_H_
#define _HEADER_LIB_CYCLE_VG_H_
#include "../base/config.h"
#include <stdlib.h>
#include "../base/data_struct.h"
#include "../base/logger.h"
#include "../base/types.h"
// #include "cycle_general.h"
#ifdef __cplusplus
extern "C"
{
#endif


struct cycle_vg_system
{
	//Constant parameters:
	//Number of vertices of the system
	size_t	n;
	//Maximum number of arcs of the system
	size_t 	nam;
	
	//Constant Parameters requiring manual enabling after initialization:
	//Maximum number of incoming arcs for each vertex
	size_t	nim;
	//Maximum number of outgoing arcs for each vertex
	size_t	nom;
	
	//Graph construction variables:
	//Current number of arcs
	size_t	na;
	//Order of vertices
	size_t* restrict go;
	//Inverse of go
	size_t* restrict goi;
	/* Graph representation for arcs out with linked list.
	 * Each value j in linked list i corresponds to arc (i,j).
	 * First item of each linked list i is specified in gaof.
	 */
	struct data_ll	gao;
	//First item of each linked list i of gao, or -1 of not exist.
	size_t* restrict	gaof;
	//Below for arcs in
	struct data_ll	gai;
	size_t* restrict	gaif;
	//Number of incoming arcs for each vertex
	size_t* restrict gni;
	//Number of outgoing arcs for each vertex
	size_t* restrict gno;
	
	
	
	//Loop detection temporary variables:
	//Vertices visitedness forward/backward, i.e. membership of F,B
	unsigned char* restrict	lvf;
	unsigned char* restrict	lvb;
	//Vertices to be visited forward/backward, i.e. membership of FL,BL
	struct data_heap	lvfl;
	struct data_heapdec	lvbl;
	/* Current arc id of those from/to a specific vertex.
	 * (i,lao[i]) is the current out arc from i during the search, indexed by gao.
	 * (lai[i],i) is the current in arc to i during the search, index by gai.
	 */
	size_t* restrict	lao;
	size_t* restrict	lai;
	//Buffer for calculation during loop detection and order maintenance.
 	size_t*	buff;
 	size_t*	buff2;
};

/* Initialize cycle detection system with vertex count and max number of arc count
 * vg:		Cycle detection system.
 * dim:		Number of vertices.
 * amax:	Max number of arcs.
 * Return:	0 on success.
 */
int cycle_vg_init(struct cycle_vg_system* restrict vg,size_t dim,size_t amax);
struct cycle_vg_system* cycle_vg_new(size_t dim,size_t amax);
int cycle_vg_free(struct cycle_vg_system* restrict vg);

/* Re-initialize existing cycle detection system to the same size.
 * vg:		Cycle detection system.
 * Return:	0 on success.
 */
int cycle_vg_empty(struct cycle_vg_system* restrict vg);

/* Obtain the number of vertices of the system
 * vg:		Cycle detection system.
 * Return:	Number of vertices on success.
 */
static inline size_t cycle_vg_dim(const struct cycle_vg_system* restrict vg);

/* Add arc v1->v2 to current graph in vg without loop checks.
 * vg:		Cycle detection system.
 * v1:		Source of arc
 * v2:		Destination of arc
 * Return:	1 if arc full, or otherwise 0 for success.
 */
static inline int cycle_vg_add_arc(struct cycle_vg_system* restrict vg,size_t v1,size_t v2);

// Fix vertex order array base on its inverse, and vice versa
static inline void cycle_vg_fix_go(struct cycle_vg_system* restrict vg);
static inline void cycle_vg_fix_goi(struct cycle_vg_system* restrict vg);

/* Restore order of vertices after adding a backward arc.
 * vg:		Cycle detection system.
 * vv:		Source of newly added arc
 */
void cycle_vg_restore_order(struct cycle_vg_system* restrict vg,size_t vv);

/* Try to add arc v1->v2 to current graph in vg.
 * vg:		Cycle detection system.
 * v1:		Source of arc
 * v2:		Destination of arc
 * Return:	0 if success, or 1 if failed because of loop or full arc.
 */
int cycle_vg_add(struct cycle_vg_system* restrict vg,size_t v1,size_t v2);

/* Extracts graph representation into matrix form.
 * vg:		Cycle detection system.
 * g:		(n,n) destination matrix. g[i,j]=1 if arc (i,j) exists, and 0 if not.
 */
void cycle_vg_extract_graph(const struct cycle_vg_system* restrict vg,MATRIXUC* g);





static inline size_t cycle_vg_dim(const struct cycle_vg_system* restrict vg)
{
	return vg->n;
}

static inline int cycle_vg_add_arc(struct cycle_vg_system* restrict vg,size_t v1,size_t v2)
{
	size_t	ret;
	
	if((vg->gno[v1]>=vg->nom)||(vg->gni[v2]>=vg->nim))
		return 1;
	ret=data_ll_insert_before(&vg->gao,vg->gaof[v1],v2);
	if(ret==(size_t)-1)
		return 1;
	vg->gaof[v1]=ret;
	ret=data_ll_insert_before(&vg->gai,vg->gaif[v2],v1);
	vg->gaif[v2]=ret;
	vg->na++;
	vg->gno[v1]++;
	vg->gni[v2]++;
	return 0;
}

static inline void cycle_vg_fix_go(struct cycle_vg_system* restrict vg)
{
	size_t 	i;
	for(i=0;i<vg->n;i++)
		vg->go[vg->goi[i]]=i;
}

static inline void cycle_vg_fix_goi(struct cycle_vg_system* restrict vg)
{
	size_t 	i;
	for(i=0;i<vg->n;i++)
		vg->goi[vg->go[i]]=i;
}


#ifdef __cplusplus
}
#endif
#endif
