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
/* This file contains the conversion from log likelihood ratio to probabilities
 *
 */

#ifndef _HEADER_LIB_PIJ_LLRTOPIJ_H_
#define _HEADER_LIB_PIJ_LLRTOPIJ_H_
#include "../base/config.h"
#include "../base/gsl/histogram.h"
#include "../base/types.h"
#ifdef __cplusplus
extern "C"
{
#endif

/* Smoothen true ratio histogram to make it monotonically increasing
 * Method:	1.	Construct increasing upper bound histogram (hup)
 * 				so every bin is max(self,left neighbor of self)
 * 			2.	Construct increasing lower bound histogram (hlow)
 * 				so every bin is min(self,right neighbor of self)
 * 			3.	Return (hup+hlow)/2
 * h:		[n] double histogram to be smoothened
 * 			Input as histogram of true ratio and,
 * 			output as smoothened histogram of true ratio
 * n:		size of histogram
 * Return:	0 on success.
 */
int pij_llrtopij_histogram_force_incremental(gsl_histogram *h);

/* Smoothen 1D histogram with buffer provided as the following:
 * 1. Calculate bin differences.
 * 2. Extrapolates with fixed boundary values to ensure same size of histogram
 * 3. Convolution with Gaussian filter.
 * 4. Calculate cumulative sum
 * 5. Scale and shift to original head-tail positions
 * WARNING:	Gaussian filter applied using bin index as distance measure,
 * 			not the actual histogram range.
 * h:		[n] double histogram to be smoothened
 * n:		size of histogram
 * sigma:	sigma of Gaussian filter
 * ncut:	size of Gaussian convolution vector is 2*ncut+1
 * vlarge:	[n+2*ncut-1] Buffer for data before convolution
 * vconv:	[2*ncut+1] Buffer for convolution mask.
 * Return:	0 on success
 */
void pij_llrtopij_histogram_smoothen_buffed(gsl_histogram *h,double sigma,size_t ncut,VECTORD *vlarge,VECTORD *vconv);
// Smoothen 1D histogram with self allocated buffer
int pij_llrtopij_histogram_smoothen(gsl_histogram *h,double sigma,size_t ncut);

/* Construct central value histogram from bounded histogram for interpolation.
 * In central value histogram, bin[i] is the value at range[i].
 * h:		(nbin) Input bounded histogram
 * hc:		(nbin+2) Output central value histogram
 */
void pij_llrtopij_histogram_to_central(const gsl_histogram *h,gsl_histogram* hc);

/* Construct central value histogram from bounded histogram for interpolation.
 * Uses pij_gassist_llrtopij_histogram_to_central.
 * h:		(nbin) Input bounded histogram
 * Return:	(nbin+2) Output central value histogram
 */
gsl_histogram* pij_llrtopij_histogram_central_from(const gsl_histogram *h);

/* Use central histogram to estimate distribution probabilities of any point
 * within the histogram range. Linear intepolation is used.
 * Points outside histogram range gives boundary output
 * hc:	central histogram for estimation
 * d:	data (x coordinates of histogram) to be estimated their probabilities
 * ans:	output of estimated probabilities
 */
void pij_llrtopij_histogram_interpolate_linear(const gsl_histogram *hc,const VECTORF* d,VECTORF* ans);
 
/* Allocate buffer for histogram conversion in pij_llrtopij_convert_histograms_buffed.
 * n:	Number of histogram bins. This must match pij_llrtopij_convert_histograms_buffed.
 * vb1,
 * vb2:	Output locations of allocated buffers.
 * Return:	0 on success.
 */
int pij_llrtopij_convert_histograms_make_buffs(size_t n,VECTORD** vb1,VECTORD** vb2);

/* Convert density histograms of null and real distribution into probability central
 * histogram with buffer provided. Both histograms must be distributions
 * (sum to unity and nonnegative).
 * hreal:	(n) Real density histogram to convert from. Also changed in calculation.
 * vnull:	(n) Null density histogram in vector format. Also changed in calculation.
 * hc:		(n+2) Central probability histogram as output.
 * vb1,
 * vb2:		Buffers needed for conversion. To allocate buffers, use
 			pij_llrtopij_convert_histograms_make_buffs.
 */
void pij_llrtopij_convert_histograms_buffed(gsl_histogram* hreal,VECTORD* vnull,gsl_histogram* hc,VECTORD* vb1,VECTORD* vb2);

/* Convert density histograms of null and real distribution into probability central histogram. Both histograms must be distributions (sum to unity and nonnegative).
 * hreal:	(n) Real density histogram to convert from. Also changed in calculation.
 * vnull:	(n) Null density histogram in vector format. Also changed in calculation.
 * hc:		(n+2) Central probability histogram as output.
 * Return:	0 if success.
 */
int pij_llrtopij_convert_histograms(gsl_histogram* hreal,VECTORD* vnull,gsl_histogram* hc);









#ifdef __cplusplus
}
#endif
#endif
