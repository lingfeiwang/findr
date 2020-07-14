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
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <tgmath.h>
#include "gsl/sf.h"
#include "gsl/math.h"
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


/* This function is modified from  GNU Scientific Library (GSL) version 1.16.
 * See https://www.gnu.org/software/gsl/.
 */
int math_sf_2F1_m1(const double a, const double b, const double c,const double x, gsl_sf_result * result)
{
  double sum_pos = 0.0;
  double sum_neg = 0.0;
  double del_pos = 0.0;
  double del_neg = 0.0;
  double del = 0.0;
  double k = 0.0;
  int i = 0;

  if(fabs(c) < GSL_DBL_EPSILON) {
    result->val = 0.0; /* FIXME: ?? */
    result->err = 1.0;
    return 1;
  }

  do {
    if(++i > 30000) {
      result->val  = sum_pos - sum_neg;
      result->err  = del_pos + del_neg;
      result->err += 2.0 * GSL_DBL_EPSILON * (sum_pos + sum_neg);
      result->err += 2.0 * GSL_DBL_EPSILON * (2.0*sqrt(k)+1.0) * fabs(result->val);
      return 1;
    }
    del *= (a+k)*(b+k) * x / ((c+k) * (k+1.0));  /* Gauss series */

    if(del > 0.0) {
      del_pos  =  del;
      sum_pos +=  del;
    }
    else if(del == 0.0) {
      /* Exact termination (a or b was a negative integer).
       */
      del_pos = 0.0;
      del_neg = 0.0;
      break;
    }
    else {
      del_neg  = -del;
      sum_neg -=  del;
    }

    k += 1.0;
  } while(fabs((del_pos + del_neg)/(sum_pos-sum_neg)) > GSL_DBL_EPSILON);

  result->val  = sum_pos - sum_neg;
  result->err  = del_pos + del_neg;
  result->err += 2.0 * GSL_DBL_EPSILON * (sum_pos + sum_neg);
  result->err += 2.0 * GSL_DBL_EPSILON * (2.0*sqrt(k) + 1.0) * fabs(result->val);

  return 0;
}
