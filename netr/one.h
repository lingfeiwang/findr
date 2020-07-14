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
/* This lib contains the implementation of network reconstruction
 * algorithms for a single network.
 */

#ifndef _HEADER_LIB_NETR_ONE_H_
#define _HEADER_LIB_NETR_ONE_H_
#include "../base/config.h"
#include <stdio.h>
#include "../base/types.h"

#ifdef __cplusplus
extern "C"
{
#endif

/* Construct a deterministic single best Direct Acyclic Graph from prior pij information,
 * and stop when the number of edges reaches threshold or no edge can be added.
 * This method sorts pij values and attempt to add edges from the most likely one,
 * therefore named 'greedy'.
 * p:		(n,n) for pij matrix
 * net:		(n,n) for constructed network. net[i,j]=1 if edge (i,j) exists, 0 if not.
 * nam:		Maximum number of edges. Set to (size_t)-1 for unlimited.
 * nimax:	Maximum number of incoming edges for each node. Set to (size_t)-1 for unlimited.
 * nomax:	Maximum number of outgoing edges for each node. Set to (size_t)-1 for unlimited.
 * Return:	Number of edges, or 0 if failed.
 */
size_t netr_one_greedy(const MATRIXF* p,MATRIXUC* net,size_t nam,size_t nimax,size_t nomax);

/* Construct a deterministic single best Direct Acyclic Graph from prior pij information,
 * and stop when the number of edges reaches threshold or no edge can be added.
 * This method sorts pij values and attempt to add edges from the most likely one.
 * Additional information is obtained in the output network variable.
 * p:	(n,n) for pij matrix
 * net:	(n,n) for constructed network. net[i,j]=0 indiates the edge is never tried.
 		net[i,j]!=0 indicates the edge has been tried. Its absolute values(=x) indicates
 		the edge is tried at the x-th edge addition attempt. net[i,j]>0 indicates successful
 		edge addition and <0 indicates failure.
 * time:(n,n) for CPU time passed from starting to add edges to finish trying
 * 		to add this edge in CPU seconds.
 * nam:	Maximum number of edges.
 * nimax:	Maximum number of incoming edges for each node. Set to (size_t)-1 for unlimited.
 * nomax:	Maximum number of outgoing edges for each node. Set to (size_t)-1 for unlimited.
 * Return:	Number of edges, or 0 if failed.
 */
size_t netr_one_greedy_info(const MATRIXF* p,MATRIXL* net,MATRIXD* time,size_t nam,size_t nimax,size_t nomax);

#ifdef __cplusplus
}
#endif
#endif
