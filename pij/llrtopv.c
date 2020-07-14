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
#include "../base/config.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../base/macros.h"
#include "../base/threading.h"
#include "llrtopv.h"

void pij_llrtopvm(MATRIXF* p,size_t n1,size_t n2)
{
	#pragma omp parallel
	{
		size_t	m1,m2;
	
		threading_get_startend(p->size1,&m1,&m2);
		if(m2>m1)
		{
			MATRIXFF(view)	mvp=MATRIXFF(submatrix)(p,m1,0,m2-m1,p->size2);
			pij_llrtopvm_block(&mvp.matrix,n1,n2);
		}
	}
}














