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
/* This lib contains histogram function
 */

#ifndef _HEADER_LIB_HISTOGRAM_H_
#define _HEADER_LIB_HISTOGRAM_H_
#include "config.h"
#include "gsl/histogram.h"
#include "types.h"
#include "general_alg.h"
#include "math.h"
#ifdef __cplusplus
extern "C"
{
#endif


/*********************************************************
 * Histogram bin range constructions
 *********************************************************/

// Uniform bins
/* Construct uniform bin widths for a specific range.
 * dmin:	Lower bound of the range.
 * dmax:	Upper bound of the range.
 * n:		Number of bins to construct.
 * Return:	Array of bin range
 */
double* histogram_uniformbins(double dmin,double dmax,size_t n);

// Bins of (near) equal count
/* Construct histogram bin ranges of (approximately) equal count according to distribution with PDF given.
 * Approximate CDFs are obtain with interpolated PDF functions (series of order 0).
 * The approximation produces errors which are supposed to be fine for histogram range purposes.
 * n:		Number of bins to construct
 * binrange:(n+1) Return of constructed bin ranges. Must be prefilled
 			at [0] for minimum bin range and at [n] maximum bin range.
 * pdfs:	Function to calculate PDFs for multiple locations together
 * param:	Parameters of func
 * Return:	0 on success.
 */
int histogram_equalbins_fromnullpdfs(size_t n,double* binrange,void (*pdfs)(const VECTORD*,VECTORD*,const void*),const void* param);

/* Construct histogram bin ranges of (approximately) equal count according to distribution with CDF given.
 * Approximate CDFs are obtain with interpolated PDF functions (series of order 0).
 * The approximation produces errors which are supposed to be fine for histogram range purposes.
 * For parameters, see math_cdf_quantile.
 */
static inline void histogram_equalbins_fromnullcdf(size_t n,double left,double right,double (*func)(double,const void*),const void* param,double eps,double* ans);

/* Bins of with smooth transition from bins of equal count on left side
 * to bins of equal width on right side.
 */
/* Automatically determine the count of histogram bins.
 * n:	Number of samples.
 * Return:	Number of histogram bins determined.
 */
size_t histogram_unequalbins_param_count(size_t n);

/* Construct unequal bin ranges from known null distribution CDF.
 * Uses smooth transition from bins of equal count at left (-) side to
 * bins of equal width at right(+) side
 * n:		Number of bins
 * binrange:(n+1) Return of constructed bin ranges. Must be prefilled
 			at [0] for minimum bin range and at [n] maximum bin range.
 * cdf:	null distribution CDF
 * param:	Parameters of func
 */
void histogram_unequalbins_fromnullcdf(size_t n,double* binrange,double (*cdf)(double,const void*),const void* param);

/* Construct histogram unequal bin ranges which smoothly transits from
 * bins of equal count to bins of equal width, according to distribution with PDF given.
 * Uses histogram_equalbins_fromnullpdfs for generation of bin of equal counts.
 * Uses histogram_unequalbins_fromequalbins to convert to equal bins.
 * n:		Number of bins to construct
 * binrange:(n+1) Return of constructed bin ranges. Must be prefilled
 			at [0] for minimum bin range and at [n] maximum bin range.
 * pdfs:	Function to calculate PDFs for multiple locations together
 * param:	Parameters of func
 * Return:	0 on success.
 */
int histogram_unequalbins_fromnullpdfs(size_t n,double* binrange,void (*pdfs)(const VECTORD*,VECTORD*,const void*),const void* param);

/* Bins of with sharp transition from bins of equal count on left side
 * to bins of equal width on right side.
 */
//Struct used by histogram_unequalbins_exact
struct histogram_unequalbins_exact_param
{
	double ebwrate;
};

/* Applies histogram_unequalbins_param_sizing and histogram_unequalbins_exact_vf
 * to construct histogram ranges for data.
 * d:		data
 * param:	parameters to be passed on.
 * Return:	0 on failure.
 */
gsl_histogram* histogram_unequalbins_exact(const VECTORF* d,const void* param);

//Struct used by histogram_unequalbins_exp
struct histogram_unequalbins_exp_param
{
	size_t nsmnv;
	double ebwrate;
};

/* Applies histogram_unequalbins_param_sizing and histogram_unequalbins_exp_vf
 * to construct histogram ranges for data.
 * d:		data
 * param:	parameters to be passed on.
 * Return:	0 on failure.
 */
gsl_histogram* histogram_unequalbins_exp(const VECTORF* d,const void* param);

/* Produce unequal histogram bin parameters (for histogram_unequalbins_exp_vf) using the following strategy:
 * n=d.size
 * nbinsplit=1-ebwrate*3/sqrt(n)
 * nbin1=min(sqrt(n)-ebwrate*3,100)
 * d:			Data set
 * ebwrate:		See above
 * nbinsplit:	Return of nbinsplit (see histogram_unequalbins_exp_vf)
 * nbin1:		Return of nbin1 (see histogram_unequalbins_exp_vf)
 * Return:		0 on success
 */
int histogram_unequalbins_param_sizing(const VECTORF* d,double ebwrate,double* restrict nbinsplit, size_t* restrict nbin1);

/* Construct unequal histogram bin boundaries for data set for VECTORF format. For data within the <nbinsplit quantile,
 * apply bin widths for approximately equal count for a total of nbin1 bins. For data outside <nbinsplit quantile,
 * use equal bin width which is the last bin width of the previous case (which has unequal bin width).
 * Approximation avoids sorting and accelerates boundary decision.
 * The function first constructs a histogram >100 times larger than the desired one.
 * Within each bin, uniform density is assumed to obtain the approximate bin ranges as output.
 * d:			Data set
 * nsmnv:		Sample count minus genotype value count. Used to estimate null distribution as exponential:
 				pdf_null=exp(-(ns-nv)*x).
 * nbinsplit:	The split location of bin width strategies. >0 and <1.
 * nbin1:		The number of bins in Strategy 1 (<nbinsplit quantile).
 * nbin:		Location to return the count of bins
 * Return:		double* on success, or 0 on failure.
 */
double* histogram_unequalbins_exp_vf(const VECTORF* d,size_t nsmnv,double nbinsplit,size_t nbin1,size_t* restrict nbin);

/* Construct unequal histogram bin boundaries for data set for VECTORF format. For data within the <nbinsplit quantile,
 * apply bin widths for exact equal count for a total of nbin1 bins. For data outside <nbinsplit quantile,
 * use equal bin width which is the last bin width of the previous case (which has unequal bin width).
 * d:			Data set
 * nbinsplit:	The split location of bin width strategies. >0 and <1.
 * nbin1:		The number of bins in Strategy 1 (<nbinsplit quantile).
 * 				Also acts as return for actual number after construction.
 * nbin:		Location to return the count of bins
 * Return:		double* on success, or 0 on failure.
 */
double* histogram_unequalbins_exact_vf(const VECTORF* d,double nbinsplit,size_t* restrict nbin1,size_t* restrict nbin);


/*********************************************************
 * Histogram bin range operations
 *********************************************************/

/* Fluctuate inner bin edges and expand outer bin edges
 * range:	bin range to be fluctuated
 * n:		number of bins (=size of range-1)
 * diff:	amount of fluctuation/expansion relative to bin width
 */
void histogram_fluc_binrange(double* restrict range,size_t n,double diff);

/* Construct finer histogram ranges from existing sparse one.
 * Each bin is evenly split into n bins.
 * hnew[i*n+j]=(hold[i]*(n-j)+hold[i+1]*j)/n
 * 
 * hold:	(nbin+1) Input old (sparse) histogram ranges in VECTORDF format.
 * hnew:	(nbin*n+1) Output new (fine) histogram ranges.
 * n:		Multiplication factor for bin count.
 */
void histogram_finer_range(const VECTORD* hold,VECTORD* hnew,size_t n);

/* Construct finer histogram centrals from existing sparse bin ranges.
 * Each bin is evenly split into n centrals.
 * hnew[i*n+j]=hold[i]+(hold[i+1]-hold[i])*(j+1/2)/n
 * 
 * hold:	(nbin+1) Input old (sparse) histogram ranges in VECTORDF format.
 * hnew:	(nbin*n) Output new (fine) histogram ranges.
 * n:		Multiplication factor for bin count.
 */
void histogram_finer_central(const VECTORD* hold,VECTORD* hnew,size_t nsp);


/* Reshapes bin ranges of equal count into bin ranges of unequal width
 * according to following method:
 * 1.	Calculate mean width =(binrange[n]-binrange[0])/n
 * 2. 	Rescale each bin width according to their index as
 *			binwidth[i]=exp(log(binwidth[i])*(n-i)/n+log(meanwidth)*i/n).
 *			This achieves a smooth transition from uniform count to
 *			uniform width bin ranges.
 * 3.	All bin widths are rescaled identically so total width remain unchanged.
 * Parameters:
 * n:		Bin count
 * binrange:(n+1) original equal count bin range.
 */
void histogram_unequalbins_fromequalbins(size_t n,double* binrange);


/*********************************************************
 * Static functions
 *********************************************************/

static inline void histogram_equalbins_fromnullcdf(size_t n,double left,double right,double (*func)(double,const void*),const void* param,double eps,double* ans)
{
	math_cdf_quantile(n,left,right,func,param,eps,ans);
}









#ifdef __cplusplus
}
#endif
#endif

