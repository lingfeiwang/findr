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
/* This lib contains mathematical functions:
 * 1:	Special functions
 * 2:	Cumulative density function related
 */

#ifndef _HEADER_LIB_MATH_H_
#define _HEADER_LIB_MATH_H_
#include "config.h"
#include <string.h>
#include <assert.h>
#include "gsl/math.h"
#include "gsl/sf.h"

#ifdef __cplusplus
extern "C"
{
#endif

/**************************************************
 * Special functions
 **************************************************/
// Calculates ln(Gamma(n/2))
static inline double math_sf_lngammahalf(size_t n);

// Calculates exp(x)-1, where x can be close to 0.
static inline double math_sf_expminusone(double x);

// Calculates log(x+1), where x can be close to 0.
static inline double math_sf_logplusone(double x);

// Calculates hypergeometric function minus 1, i.e. 2F1(a,b,c;x)-1
int math_sf_2F1_m1(const double a, const double b, const double c,const double x, gsl_sf_result * result);

/**************************************************
 * CDF related functions
 **************************************************/

/* Locate quantiles of CDF with binary search.
 * start:	Start quantile location to be calculated
 * step:	Step of quantile location
 * n:		Number of quantiles to calculate
 * left:	All quantiles are known >left.
 * right:	All quantiles are known <left.
 * func:	CDF function
 * param:	Parameter accepted by CDF function.
 * eps:		Calculation stops when right-left<eps.
 * ans:		(n) output location
 */
void math_cdf_quantile_calc(double start,double step,size_t n,double left,double right,double (*func)(double,const void*),const void* param,double eps,double* ans);

/* Locate quantiles of CDF with binary search.
 * n:		Number of quantiles to calculate
 * left:	All quantiles are known >left.
 * right:	All quantiles are known <left.
 * func:	CDF function
 * param:	Parameter accepted by CDF function.
 * eps:		Calculation stops when right-left<eps.
 * ans:		(n-1) output location
 */
static inline void math_cdf_quantile(size_t n,double left,double right,double (*func)(double,const void*),const void* param,double eps,double* ans);


/**************************************************
 * Static functions
 **************************************************/
 
static inline double math_sf_lngammahalf(size_t n)
{
	size_t n2=n/2;
	if(n%2)
	{
		if(!n2)
			return M_LNPI/2;
		else if(n2<GSL_SF_FACT_NMAX)
			return M_LNPI/2+gsl_sf_lndoublefact((unsigned int)(n-2))-(double)(n-1)*M_LN2/2;
		else
			return gsl_sf_lngamma(((double)n)/2);
	}
	else
	{
		if(!n2)
			return INFINITY;
		else if (n2<=GSL_SF_DOUBLEFACT_NMAX)
			return gsl_sf_lnfact((unsigned int)(n2-1));
		else
			return gsl_sf_lngamma((double)n2);
	}
}

static inline double math_sf_expminusone(double x)
{
	return fabs(x)>1E-4?exp(x)-1:x*(1+(x/2)*(1+(x/3)*(1+x/4)));
}


static inline double math_sf_logplusone(double x)
{
	return fabs(x)>1E-4?log(x+1):x*(1-(x/2)*(1+((x*2)/3)*(1-(x*3)/4)));
}

static inline void math_cdf_quantile(size_t n,double left,double right,double (*func)(double,const void*),const void* param,double eps,double* ans)
{
	double step=1./(double)n;
	assert(n>1);
	assert((func(left,param)<step)&&(func(right,param)>1-step));
	math_cdf_quantile_calc(step,step,n-1,left,right,func,param,eps,ans);
}

#ifdef __cplusplus
}
#endif
#endif
