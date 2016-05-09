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
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <tgmath.h>
#include "types.h"
#include "logger.h"
#include "math.h"


void math_cdf_quantile_calc(double start,double step,size_t n,double left,double right,double (*func)(double,const void*),const void* param,double eps,double* ans)
{
	size_t i,nleft,nright;
	double mid,midv,t1;
	
	while(n)
	{
		if(right-left<eps)
		{
			double width=(right-left)/(double)n;
			ans[0]=left+width/2;
			for(i=1;i<n;i++)
				ans[i]=ans[i-1]+width;
			return;
		}
		mid=(left+right)/2;
		midv=func(mid,param);
		t1=midv-start;
		nleft=t1>0?(size_t)ceil(t1/step):0;
		if(nleft>n)
			nleft=n;
		nright=n-nleft;
		//Small side first
		if(nleft<=nright)
		{
			if(nleft)
			{
				math_cdf_quantile_calc(start,step,nleft,left,mid,func,param,eps,ans);
				start+=step*(double)nleft;
				n=nright;
				ans+=nleft;
			}
			left=mid;
		}
		else
		{
			if(nright)
			{
				math_cdf_quantile_calc(start+(double)nleft*step,step,nright,mid,right,func,param,eps,ans+nleft);
				n=nleft;
			}
			right=mid;			
		}
	}
}

