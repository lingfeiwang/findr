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
/* This file contains the conversion from log likelihood ratio to probabilities,
 * when all gene pairs are considered together.
 *
 */

#ifndef _HEADER_LIB_PIJ_LLRTOPIJ_TOT_H_
#define _HEADER_LIB_PIJ_LLRTOPIJ_TOT_H_
#include "../base/config.h"
#include "../base/gsl/histogram.h"
#include "../base/types.h"
#ifdef __cplusplus
extern "C"
{
#endif

/* Converts log likelihood ratio data into probabilities based on the reference null log likelihood ratios.
 * The approach assumes the real data is the outcome of a mixture of alternative and null hypotheses.
 * The approach first simulates null log likelihood ratios from permuted the real data, and then
 * draw histograms of log likelihood ratios for both real and only null situations.
 * Assuming on one tail of the histogram, the real data should contain mostly (if not only) null hypothesis.
 * The null fraction can be known at the tail, and therefore at every bin based on the simulated
 * null histogram shape. The final probability is defined to be the ratio of alternative hypothesis versus
 * total, only on that bin.
 * NOTE:		Real and simulated data should both be non-negative. Real data should reach true probability->0
 * 				at the limit of log likelihood ratio ->0+.
 * d:			Real log likelihood ratios, to be compared against permuted (null) data. Real data is assumed to contain
 *				both null and alternative hypotheses.
 * dconv:		Real log likelihood ratios to be converted to probabilities, based on the comparison between d and null distribution.
 * 				Same with d in most cases.
 * ans:			Output vector for the probabilities. Have same size with dconv, and each provides the estimated probability
 * 				of the corresponding entry of dconv.
 * fh:			Function to obtain real histogram from data d. It should automatically determine bin widths. Takes parameters:
 * 		Param 1:	Data points. (i.e. d).
 * 		Param 2:	Parameter of the function to be specified in a later parameter, and passed on by pij_llrtopij_convert.
 * 		Return:		Constructed gsl_histogram if success, or 0 if fail.
 * fnh:			Function to obtain null histogram. Can be simulation based or analytical. They are defined in pij_gassist/nullhist.h.
 * 		Param 1:	Parameter to pass on to the function and defined later in the parameters.
 * 		Param 2:	Existing empty histogram to fill for the function.
 * 		Param 3:	Number of parallel threads.
 * 		Return:		0 if success.
 * ph:			Parameter (2) to be passed onto fh.
 * pnh:			Parameter (1) to be passed onto fnh.
 * nullratio:	if not NULL, return the ratio of null data in real data into nullratio.
 * Return:		0 on success.
 */
int pij_llrtopij_tot_convert(const VECTORF* d,const VECTORF* dconv,VECTORF* ans,gsl_histogram* (*fh)(const VECTORF*,const void*),int (*fnh)(const void*,gsl_histogram*),const void* ph,const void* pnh,double* nullratio);


























#ifdef __cplusplus
}
#endif
#endif
