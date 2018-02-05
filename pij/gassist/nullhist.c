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
#include "../../base/config.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "../../base/logger.h"
#include "../../base/gsl/histogram.h"
#include "../nullhist.h"
#include "nullhist.h"

int pij_gassist_nullhists(gsl_histogram** h[4],size_t nt,size_t ns,size_t nv,const FTYPE dmax[4])
{
	//Construct null density histograms
	h[0]=pij_nullhist((double)dmax[0],nv,nt,1,1,1,ns-2);
	h[1]=pij_nullhist((double)dmax[1],nv,nt,1,1,1,ns-3);
	h[2]=pij_nullhist((double)dmax[2],nv,nt,1,2,1,ns-3);
	h[3]=pij_nullhist((double)dmax[3],nv,nt,0,1,1,ns-3);
	if(h[0]&&h[1]&&h[2]&&h[3])
		return 0;

	LOG(1,"pij_nullhist failed.")
	return 1;
}





















